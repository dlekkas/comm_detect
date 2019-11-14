#include "../include/plm.h"
#include "../include/graph.h"
#include "../include/modularity.h"
#include <omp.h>

#include <map>

#define NUM_SPLIT 100

 void print(std::vector<int> const &input)
{
	for (int i = 0; i <(int)(input.size()); i++) {
		std::cout << input.at(i) << ' ';
	}
	std::cout << '\n';
}


void print_map(std::map<int,std::vector<int>> myMap) {
for(map<int, std::vector<int>>::const_iterator it = myMap.begin(); it != myMap.end(); ++it)
{
    std::cout << "key:" << it->first << endl;
    print(it-> second);
}
}

GraphComm PLM::coarsen(GraphComm* g_initial) {

	GraphComm g;
	std::vector<int> comm = (*g_initial).communities;

	std::vector<int> g_vertexes;
	int i, j;

	for (i=0; i<(*g_initial).n; i++) {
		int c = comm[i];
		
		if (std::find(g_vertexes.begin(), g_vertexes.end(), c) == g_vertexes.end()) {
			g_vertexes.push_back(c);
		}
	}
	sort(g_vertexes.begin(), g_vertexes.end());
	g.n = g_vertexes.size();
	//print(g_vertexes);

	 
	for (i=0; i<g.n; i++)	
		(*g_initial).com_map.insert(std::pair<int,int>(g_vertexes[i], i)); 

	network new_net;
	int new_net_array[g.n][g.n];

	for (i=0; i<g.n; i++)
		for (j=0; j<g.n; j++)
			new_net_array[i][j] = 0;
	
	for (i=0; i<(*g_initial).n; i++) {
		int c_i = (*g_initial).com_map[comm[i]];
		vector<pair<node_id, weight>> neighbors = g_initial->net[i];
		for (auto it=neighbors.begin(); it<neighbors.end(); ++it) {
			int c_j = (*g_initial).com_map[comm[it->first]];
			new_net_array[c_i][c_j] += it->second;
		}
		
	}


	for (i=0; i<g.n; i++) {
		std::vector<pair<node_id, weight>> v;
		weight i_volume = new_net_array[i][i];
		for (j=0; j<g.n; j++) {
			if (new_net_array[i][j] > 0) {
				v.push_back(make_pair(j, new_net_array[i][j]));
				i_volume += new_net_array[i][j];
			}
		}
		new_net.push_back(v);
		g.volumes.push_back(i_volume);
			
	}

	g.net = new_net;
	g.weight_net = (*g_initial).weight_net; //the sum of all edges remains the same

	return g;

}


std::vector<int> PLM::prolong(GraphComm g_initial, std::vector<int> coarsened_comm) {
	
	std::vector<int> init_comm = g_initial.communities;
	std::vector<int> new_comm(g_initial.n, 0);

	for (int i=0; i<g_initial.n; i++) {
		int i_comm = init_comm[i];
		new_comm[i] = coarsened_comm[g_initial.com_map[i_comm]];
	}
	return new_comm;
}

std::pair <int, float> max_pair_arg (std::pair <int, float> r, std::pair <int, float> n) {
        return (n.second >= r.second) ? n : r;
}

int PLM::ReturnCommunity(int i, GraphComm g) {

	std::vector<pair<node_id, weight>> n_i = g.net[i];
        std::unordered_map<int,float> mod_map;
	for (auto neighbor_it = n_i.begin(); neighbor_it < n_i.end(); ++neighbor_it) {
                // for each v, neighbor, compute mod_diff
		if (g.communities[i] != g.communities[neighbor_it->first]) {
                	float m = compute_modularity_difference(i, neighbor_it->first, g);
                	mod_map[neighbor_it->first] = m;
			//cout << neighbor_it->first <<  endl;
		}
        }

	std::pair<int, float> max_pair = std::make_pair(i, 0.0);

	for (size_t b = 0; b < mod_map.bucket_count(); b++) {
                for (auto bi = mod_map.begin(b); bi != mod_map.end(b); bi++)
                        max_pair = max_pair_arg(max_pair, *bi);
        }
	
	modularity max_diff = max_pair.second;
	if (max_diff > 0.0) {
		return g.communities[max_pair.first];
	}
	else 
		return g.communities[i]; 
}


void  PLM::Local_move(GraphComm* graph) {
	int unstable = 1;
	while (unstable) {
		//print((*graph).communities);
		//cout << "----------------------------" << endl;
		unstable = 0;
		for (int i=0; i<(*graph).n; i++) {
			int i_comm = (*graph).communities[i];
			int z = ReturnCommunity(i, *graph);
			//cout << z << endl;
			if (z != i_comm) { 
				(*graph).communities[i] = z;
				unstable=1;
			}
	
		}
	}
}

std::vector<int> PLM::Recursive_comm_detect(GraphComm g) {
	

	std::vector<int> c_singleton;
	
	for (int i=0; i<g.n; i++) {		
			c_singleton.push_back(i); 
	}
	g.communities = c_singleton;
	//print(c_singleton);
	Local_move(&g);
	//print(g.communities);
	
	if (g.communities != c_singleton) {
		GraphComm g_new = coarsen(&g);
		std::vector<int> c_coarsened = Recursive_comm_detect(g_new);
		g.communities = prolong(g, c_coarsened); 
	}
	return g.communities;

}

void get_weight_and_volumes(GraphComm* g) {
	weight sum = 0;
	std::vector<int> volumes((*g).n, 0);	
	int u=0;
	for (auto it = ((*g).net).begin(); it < ((*g).net).end(); ++it) {
		vector<pair<node_id, weight>> neighbors = *it;
		weight vol = 0;
		//add the weight of the neighbors to the total sum
		for (vector<pair<node_id, weight>>::iterator j = neighbors.begin(); j != neighbors.end(); ++j) {
			sum += j->second;
			vol += j->second;
			if (j->first == u)
				vol += j->second;
    		}	
		volumes[u++]=vol;
  	}
    	(*g).weight_net = sum;
	(*g).volumes = volumes;
}


void PLM::DetectCommunities() {

	cout << "detect comm" << endl;
	get_weight_and_volumes(&graph);
	cout << "weight: " << graph.weight_net << endl;
	graph.communities = Recursive_comm_detect(graph);
	cout << "final communities: ";
	print(graph.communities);

}

// the same for LP and Louvain. TODO: defined once 
void PLM::PrintCommunities(const std::string &file_name) {
	std::ofstream ofs;
	ofs.open(file_name, std::ios_base::out | std::ios_base::trunc);
	if (ofs.fail()) {
		std::cerr << "Error opening file " << file_name << std::endl;
		std::exit(2);
	}

	for (int i = 0; i < graph.n; i++) {
		ofs << graph.communities[i] << std::endl;
	}
}
