#include "../include/plm.h"
#include "../include/graph.h"
#include "../include/modularity.h"

#include <map>
 
void print(std::vector<int> const &input)
{
	for (int i = 0; i <(int)(input.size()); i++) {
		std::cout << input.at(i) << ' ';
	}
	std::cout << '\n';
}


void print_map(std::map<int,int> myMap) {
for(map<int, int>::const_iterator it = myMap.begin();
    it != myMap.end(); ++it)
{
    std::cout << it->first << " " << it->second << "\n";
}
}



int PLM::connected(int comm_1, int comm_2, std::vector<int> communities, std::vector<std::vector<int>> adj_list) {

	int n = communities.size();
	int weight=0;
	for (int i=0; i<n; i++) {
		if (communities[i] == comm_1) {
 			std::vector<int> neighbors_i = adj_list[i];
			for (int j=0; j<(int)(neighbors_i.size()); j++) {
				if (communities[neighbors_i[j]] == comm_2) {
 					weight+=1;
				}
			}
		}
	}
	return weight;
}


GraphComm PLM::coarsen(GraphComm* g_initial, std::vector<int> comm) {

	GraphComm g;

	std::vector<int> g_vertexes;
	int i;

	for (i=0; i<(*g_initial).n; i++) {
		int c = comm[i];
		if (std::find(g_vertexes.begin(), g_vertexes.end(), c) == g_vertexes.end()) {
			g_vertexes.push_back(c);
		}
	}
	sort(g_vertexes.begin(), g_vertexes.end());
	g.n = g_vertexes.size();
	cout << "coarsen!!!!!!!\n" ;
	print(g_vertexes);

	for (i=0; i<g.n; i++)	
		(*g_initial).com_map.insert(std::pair<int,int>(g_vertexes[i], i)); 

	network new_net;
	// TODO: make it quicker
	//g.PrintGraph();
	for (i=0; i<g.n; i++) {
		std::vector<int> neighbors_i;
		vector<pair<node_id, weight>> v;
		for (int j=0; j<g.n; j++) {
			int w = connected(g_vertexes[i], g_vertexes[j], comm, (*g_initial).adj_list); 
			if (w > 0) {
				neighbors_i.push_back(j);
				v.push_back(make_pair(j, w));
				//cout << "i: " << g_vertexes[i] << " j: " << g_vertexes[j] << " w: " << w << endl;
			}
		}	
		g.adj_list.push_back(neighbors_i);
		//print(v);
		new_net.push_back(v);		
	}

	g.net = new_net;
	return g;

}


std::vector<int> PLM::prolong(GraphComm g_initial, std::vector<int> coarsened_comm) {
	
	std::vector<int> init_comm = g_initial.communities;
	std::vector<int> new_comm;

	for (int i=0; i<g_initial.n; i++) {
		int i_comm = init_comm[i];
		new_comm.push_back(coarsened_comm[g_initial.com_map[i_comm]]);
	}
	return new_comm;
}


community PLM::get_community_vector(std::vector<int> communities, int comm) {
	community comm_vector;
	for (int i=0; i<(int)(communities.size()); i++) {
		if (communities[i] == comm)
			comm_vector.push_back(i);
	}
	return comm_vector;
}


std::vector<int> PLM::Local_move(GraphComm graph, std::vector<int> communities) {
	int unstable = 1;
	network net = graph.net;
	cout << "**************************************" << endl;
	while (unstable) {
		unstable = 0;
		print(communities);
		for (int i=0; i<graph.n; i++) {
			int i_comm = communities[i];
			community i_comm_vector = get_community_vector(communities, i_comm);
			int n_id = 0;
			modularity max_diff = -1.0;
			std::map<int,float> mod_map;
			for (int neighbor: graph.adj_list[i]) {
				int n_comm = communities[neighbor];
                       		community n_comm_vector = get_community_vector(communities, n_comm);
				// for each v, neighbor, compute mod_diff
				modularity mod_diff = compute_modularity_difference(i, i_comm_vector, n_comm_vector, net);
				mod_map.insert(std::pair<int,float>(neighbor, mod_diff));
			// choose the <id, d> = <neighbor_id, mod_diff> with the largest mod_diff
			}
			for(std::map<int, float>::iterator it=mod_map.begin(); it!=mod_map.end(); ++it) {
				if (max_diff <= it -> second) {
					max_diff = it -> second;
					n_id = it -> first;
				
				}

			}
			int z = communities[n_id];
			if (max_diff>0 and z != communities[i]) {
				communities[i] = z;
				unstable=1;
			}
	
		}
	}

	return communities;
}

std::vector<int> PLM::Recursive_comm_detect(GraphComm g) {
	

	std::vector<int> c_singleton;
	for (int i=0; i<g.n; i++)
		c_singleton.push_back(i);

	std::vector<int> c_new = Local_move(g, c_singleton);

	if (c_new != c_singleton) {
		g.communities = c_new;
		GraphComm g_new = coarsen(&g, c_new);
		//print_map(g.com_map);
		std::vector<int> c_coarsened = Recursive_comm_detect(g_new);
		c_new = prolong(g, c_coarsened); 
	}
	return c_new;

}

void PLM::DetectCommunities() {

 	graph.net = graph.CreateNetwork(); 
	graph.communities = Recursive_comm_detect(graph);
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
