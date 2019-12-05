#include "../include/plm.h"
#include "../include/graph.h"
#include "../include/modularity.h"
#include <omp.h>

#include <map>

#define NUM_SPLIT 100

 void print(std::vector<int> const &input)
{
	for (int i = 0; i < ((int) input.size()); i++) {
		std::cout << input.at(i) << ' ';
	}
	std::cout << '\n';
}

 void print_f(std::vector<float> const &input)
{
        for (int i = 0; i <(float)(input.size()); i++) {
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

	for (i=0; i<g.n; i++)
		(*g_initial).com_map.insert(std::pair<int,int>(g_vertexes[i], i));

	network new_net;

	std::vector<std::vector<int>> new_net_array(g.n, std::vector<int>(g.n, 0));
	std::vector<int> new_volumes(g.n, 0);

	#pragma omp parallel for shared(new_net_array) schedule(static, NUM_SPLIT)
	for (i=0; i<(*g_initial).n; i++) {
		int c_i = (*g_initial).com_map[comm[i]];
		vector<pair<node_id, weight>> neighbors = g_initial->net[i];
		for (auto it=neighbors.begin(); it<neighbors.end(); ++it) {
			int c_j = (*g_initial).com_map[comm[it->first]];
			#pragma omp atomic update
			new_net_array[c_i][c_j] += it->second;
		}

	}

	#pragma omp parallel for ordered shared(new_net) schedule(static, NUM_SPLIT)
	for (i=0; i<g.n; i++) {
		std::vector<pair<node_id, weight>> v;
		weight i_volume = new_net_array[i][i];
		for (j=0; j<g.n; j++) {
			if (new_net_array[i][j] > 0) {
				v.push_back(make_pair(j, new_net_array[i][j]));
				i_volume += new_net_array[i][j];
			}
		}
		new_volumes[i] = i_volume;
		#pragma omp ordered
		new_net.push_back(v);


	}

	g.volumes = new_volumes;
	g.net = new_net;
	g.weight_net = (*g_initial).weight_net; //the sum of all edges remains the same

	return g;

}


std::vector<int> PLM::prolong(GraphComm g_initial, std::vector<int> coarsened_comm) {

	std::vector<int> init_comm = g_initial.communities;
	std::vector<int> new_comm(g_initial.n, 0);

	#pragma omp parallel for schedule(static, NUM_SPLIT)
	for (int i=0; i<g_initial.n; i++) {
		int i_comm = init_comm[i];
		new_comm[i] = coarsened_comm[g_initial.com_map[i_comm]];
	}
	return new_comm;
}

std::pair <int, float> max_pair_arg (std::pair <int, float> r, std::pair <int, float> n) {
        return (n.second >= r.second) ? n : r;
}

std::pair<int, float> PLM::ReturnCommunity(int i, GraphComm g, int comm_size) {

	std::vector<pair<node_id, weight>> n_i = g.net[i];
	int j, c = g.communities[i];

	/* per thread */
        int threads = std::min(comm_size, omp_get_max_threads());
	//int threads=1;
	std::vector<std::vector<int>> weights_per_thread(threads, std::vector<int>(comm_size, 0));
	
	/* shared */
	std::vector<int> volumes(comm_size, 0);	
	std::vector<int> seen_comm(comm_size, 0);
	std::vector<pair<int, float>> results(threads, std::make_pair(c, 0.0)); /* each thread will write the best result it will find*/

	/* calculate volumes for all communties */
	//TODO: can we do it better? - compute it only once
	for (j=0; j<g.n; j++) 
		volumes[g.communities[j]] += g.volumes[j];

	/* iterate once over all neighbors and compute weights from i to all communities.
	   Update the weights for all threads */
	for (auto neighbor_it = n_i.begin(); neighbor_it < n_i.end(); ++neighbor_it) {
		int c_n = g.communities[neighbor_it->first];
		if ((int) neighbor_it->first != i) {
			if (seen_comm[c_n] == 0) {
				for (j=0; j<threads; j++)
					weights_per_thread[j][c_n] = neighbor_it->second;
				seen_comm[c_n] = 1;
			}
			else {
				for (j=0; j<threads; j++)
					weights_per_thread[j][c_n] += neighbor_it->second;
			}
		}
        }
	#pragma omp parallel num_threads(threads)
	{
	    /* Obtain thread number */
	    int tid = omp_get_thread_num();
	    std::vector<int> t_weights = weights_per_thread[tid];
	    std::vector<float> t_mod(comm_size, 0.0);
	    weight weight_c = t_weights[c];
	    weight volume_c = volumes[c] - g.volumes[i];
	    weight i_vol = g.volumes[i];
            weight n_w = g.weight_net;    
	
            /* find the id of communities that this thread will check */
	    std::vector<int> comm_to_check;
	    int c_number = tid;
	    while (c_number < comm_size) {
                comm_to_check.push_back(c_number);
		c_number += threads;
	    }    	
	    
	    for (std::vector<int>::iterator it = comm_to_check.begin() ; it != comm_to_check.end(); ++it) {
	       // compute difference in modularity for this community
	        
                float a =  ((1.0 * (t_weights[*it] - weight_c)) / n_w);
    	        float b = (1.0 * (volume_c - volumes[*it]) * i_vol) / (2 * n_w * n_w);
		t_mod[*it] = a + b;
		
	    }	   
	    t_mod[c] = 0.0;
	    std::pair<int, float> max_pair = std::make_pair(c, 0.0);
	    for (int k=0; k<comm_size; k++)
		if (t_mod[k] > max_pair.second) {
			max_pair.first = k;
			max_pair.second = t_mod[k];
		}
	    results[tid] = max_pair;

	}

	std::pair<int, float> max_pair = std::make_pair(c, 0.0);
        for (j=0; j<threads; j++) {
		if (results[j].second > max_pair.second)
			max_pair = results[j];
	}
	return max_pair;
}

std::map<int, int> PLM::Map_communities(GraphComm g) {

	std::vector<int> comms;
	std::map<int,int> com_map;
	for (int i=0; i<g.n; i++) {
		int c = g.communities[i];

		if (std::find(comms.begin(), comms.end(), c) == comms.end()) {
			comms.push_back(c);
		}
	}
	sort(comms.begin(), comms.end());

	for (int i=0; i < (int) comms.size(); i++)
		com_map.insert(std::pair<int,int>(comms[i], i));

	return com_map;

}


void  PLM::Local_move(GraphComm* graph) {


	int unstable = 1;
	while (unstable) {
		//print((*graph).communities);
		//cout << "----------------------------" << endl;
		unstable = 0;
		#pragma omp parallel for reduction(||: unstable)
		for (int i=0; i<(*graph).n; i++) {
			int i_comm = (*graph).communities[i];
			std::pair<int, float> res = ReturnCommunity(i, *graph, (*graph).n);
			if (res.first != i_comm) { // TODO: change iff the total modularity is optimized!
				(*graph).communities[i] = res.first;
				unstable=1;
			}

		}
	}

        std::map<int, int> com_map = Map_communities(*graph);
	for (int i=0; i<(*graph).n; i++)
                (*graph).communities[i] = com_map[(*graph).communities[i]];

}

std::vector<int> PLM::Recursive_comm_detect(GraphComm g) {


	std::vector<int> c_singleton(g.n, 0);
	# pragma omp parallel for schedule(static, NUM_SPLIT)
	for (int i=0; i<g.n; i++) {
		c_singleton[i] = i;
	}
	g.communities = c_singleton;
	Local_move(&g);
	
	if (g.communities != c_singleton) {
		GraphComm g_new = coarsen(&g);
		std::vector<int> c_coarsened = Recursive_comm_detect(g_new);
		g.communities = prolong(g, c_coarsened);
	}
	return g.communities;

}

void get_weight_and_volumes(GraphComm* g) {
	weight sum = 0;
	std::vector<int> volumes(g->n, 0);
	node_id u = 0;
	# pragma omp parallel for schedule(static, NUM_SPLIT)
	for (u = 0; u < (node_id) g->n; u++) {
		vector<pair<node_id, weight>> neighbors = g->net[u];
		weight vol = 0;
		weight my_weight=0;
		//add the weight of the neighbors to the total sum
		for (vector<pair<node_id, weight>>::iterator j = neighbors.begin(); j != neighbors.end(); ++j) {
			vol += j->second;
			if (j->first == u)
				my_weight = j->second;
		}
		volumes[u]=vol + my_weight;
		#pragma omp atomic update
		sum += vol;
    	}
    	g->weight_net = sum;
	g->volumes = volumes;
}


void PLM::DetectCommunities() {

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
