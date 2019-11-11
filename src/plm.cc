#include "../include/plm.h"
#include "../include/graph.h"
#include "../include/modularity.h"
#include <omp.h>

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

	// TODO: parallel
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

	for (i=0; i<g.n; i++)	
		(*g_initial).com_map.insert(std::pair<int,int>(g_vertexes[i], i)); 

	network new_net;
	for (i=0; i<g.n; i++) {
		std::vector<int> neighbors_i;
		vector<pair<node_id, weight>> v;
		for (int j=0; j<g.n; j++) {
			int w = connected(g_vertexes[i], g_vertexes[j], comm, (*g_initial).adj_list); 
			if (w > 0) {
				neighbors_i.push_back(j);
				v.push_back(make_pair(j, w));
			}
		}	
		g.adj_list.push_back(neighbors_i);
		new_net.push_back(v);		
	}

	g.net = new_net;
	g.weight_net = weight_of_network(g.net);
	g.volumes = compute_node_volumes(g.n, g.net);
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
	// contains all elements of the community comm
	return comm_vector;
}

std::pair <int, float> max_pair_arg (std::pair <int, float> r, std::pair <int, float> n) {
        return (n.second >= r.second) ? n : r;
}

int PLM::ReturnCommunity(int i, int i_comm, community i_comm_vector, std::vector<int> adj_list_i, std::vector<int> communities, network net, weight w, std::vector<int> v) {
        std::unordered_map<int,float> mod_map;
	#pragma omp parallel for
        for (auto neighbor_it = adj_list_i.begin(); neighbor_it < adj_list_i.end(); ++neighbor_it) {
                community n_comm_vector = get_community_vector(communities, communities[*neighbor_it]);
                // for each v, neighbor, compute mod_diff
                mod_map[*neighbor_it] = compute_modularity_difference(i, i_comm_vector, n_comm_vector, net, w, v);
        }

	std::pair<int, float> max_pair = std::make_pair(i, -1.0);

	#pragma omp declare reduction \
        (maxpair : std::pair<int, float> : omp_out=max_pair_arg(omp_out,omp_in)) \
        initializer(omp_priv = omp_orig)
        #pragma omp parallel for shared(mod_map) reduction(maxpair:max_pair) schedule(static, 1)
        for (size_t b = 0; b < mod_map.bucket_count(); b++) {
                for (auto bi = mod_map.begin(b); bi != mod_map.end(b); bi++)
                        max_pair = max_pair_arg(max_pair, *bi);
        }
	
	modularity max_diff = max_pair.second;
	if (max_diff > 0)
		return communities[max_pair.first];
	else 
		return i_comm; 
}


std::vector<int> PLM::Local_move(GraphComm graph, std::vector<int> communities) {
	int unstable = 1;
	network net = graph.net;
	while (unstable) {
		//print(communities);
		unstable = 0;
		#pragma omp parallel for shared(unstable, communities, graph, net) schedule(static, 1)
		for (int i=0; i<graph.n; i++) {
			int i_comm = communities[i];
			community i_comm_vector = get_community_vector(communities, i_comm);
			int z = ReturnCommunity(i, i_comm, i_comm_vector, graph.adj_list[i], communities, net, graph.weight_net, graph.volumes);
			if (z != i_comm) {
				communities[i] = z;
				unstable=1;
			}
	
		}
	}

	return communities;
}

std::vector<int> PLM::Recursive_comm_detect(GraphComm g) {
	

	std::vector<int> c_singleton;
	#pragma omp parallel for ordered schedule(static, 1) 
	for (int i=0; i<g.n; i++) {
		#pragma omp ordered 
		{ 
			c_singleton.push_back(i); 
		}
	}
	std::vector<int> c_new = Local_move(g, c_singleton);

	if (c_new != c_singleton) {
		g.communities = c_new;
		GraphComm g_new = coarsen(&g, c_new);
		std::vector<int> c_coarsened = Recursive_comm_detect(g_new);
		c_new = prolong(g, c_coarsened); 
	}
	return c_new;

}

void PLM::DetectCommunities() {

 	graph.net = graph.CreateNetwork(); 
	graph.weight_net = weight_of_network(graph.net);
	graph.volumes = compute_node_volumes(graph.n, graph.net);
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
