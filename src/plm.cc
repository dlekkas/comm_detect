#include "../include/plm.h"
#include "../include/graph.h"
#include "../include/modularity.h"
#include <omp.h>
#include <chrono>

#include <map>
#include <unordered_set>
#include <math.h>

#define NUM_SPLIT 100


//TODO: 
// 1) experiment with number of threads for local move at 2 levels + search for nested parallelism
// 2) check modularity
// 3) speed up modularity computations
std::vector<std::unordered_map<int,int>> new_net_array;	 
vector<vector<std::unordered_map<int,int>>> array_thread; 

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


void PLM::coarsen(GraphComm* g_initial, GraphComm* g, vector<vector<std::unordered_map<int,int>>>* array_thread, std::vector<std::unordered_map<int,int>>* new_net_array) {


	auto start_total = std::chrono::system_clock::now();

	auto start = std::chrono::system_clock::now();

	//GraphComm *g = new GraphComm;
	std::vector<int> comm = (*g_initial).communities;
	int n = *max_element(std::begin(comm), std::end(comm)) + 1;
	g->n = n;
		
	// create the new graph

	int threads=omp_get_max_threads();
    std::vector<int> new_volumes(n, 0);

	//cout << "Inside Coarsen with " << threads << " threads, create a new graph with " << n << " nodes " << endl;

	network new_net(n, std::vector<pair<node_id, weight>>());	
	// std::unordered_map<int,int> new_net_array[n];	 
	// std::vector<std::unordered_map<int,int>> new_net_array(n);	 
	// std::unordered_map<int,int> array_thread[threads][n];	 

	//vector<vector<std::unordered_map<int,int>>> array_thread(threads, vector<std::unordered_map<int,int>>(n)); 

	//std::unordered_map<int,int> array_thread[threads][n];

        auto end = std::chrono::system_clock::now();

        auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //std::cout << "-------------init  took time (in us) : " << total_time << std::endl;


	start = std::chrono::system_clock::now();
	
	#pragma omp parallel for num_threads(threads)
    	for (int i=0; i<(*g_initial).n; i++) {
		int tid = omp_get_thread_num();
		int c_i = comm[i];
		//std::unordered_map<int,int> *n_i = &(new_net_array[c_i]);
		vector<pair<node_id, weight>> neighbors = g_initial->net[i];
		for (auto it=neighbors.begin(); it<neighbors.end(); ++it) {
			int c_j = comm[it->first];
			//#pragma omp critical 
			//{
				//(*n_i)[c_j] += it->second;
			//	new_net_array[c_i][c_j] += it->second;
			// }
			
			(*array_thread)[tid][c_i][c_j] += it->second;
		}
	}

	end = std::chrono::system_clock::now();

	total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	//std::cout << "-------------per thread took time (in us) : " << total_time << std::endl;

	
	start = std::chrono::system_clock::now();
	#pragma omp parallel for num_threads(threads)
	for (int j=0; j<n; j++) {
		for (int i=0; i < threads; i++) {
    			for (auto r : (*array_thread)[i][j]) {
					(*new_net_array)[j][r.first] += r.second;
			}
    	}
	}

	end = std::chrono::system_clock::now();

	total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	//std::cout << "-------------reduce took time (in us) : " << total_time << std::endl;

	#pragma omp parallel for
	//int threads_c = std::min(n, omp_get_max_threads());
	for (int i=0; i<n; i++) {
		std::vector<pair<node_id, weight>> v;
		weight i_volume = 0; 
		std::unordered_map<int,int> net_i =  (*new_net_array)[i];
		for (auto it: net_i) {
			v.push_back(make_pair(it.first, it.second));
			i_volume += it.second;
			if (it.first == i)
				i_volume += it.second;
		}

		new_volumes[i] = i_volume;
		new_net[i]=v;
	}
	

	end = std::chrono::system_clock::now();

        total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //std::cout << "-------------network creation took time (in us) : " << total_time << std::endl;

	
	start = std::chrono::system_clock::now();

	g->volumes = new_volumes;
	g->net = new_net;
	g->weight_net = (*g_initial).weight_net; //the sum of all edges remains the same
	g->weight_sq = (*g_initial).weight_sq;


	end = std::chrono::system_clock::now();

        total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //std::cout << "-------------finalize took time (in us) : " << total_time << std::endl;

	//cout << "exit coarsen" << endl;
	//start = std::chrono::system_clock::now();
	//cout << std::chrono::duration_cast<std::chrono::microseconds>(start).count() << endl;

	auto end_total = std::chrono::system_clock::now();

        total_time = std::chrono::duration_cast<std::chrono::microseconds>(end_total - start_total).count();
        //std::cout << "-------------total coarsen took time (in us) : " << total_time << std::endl;
	
	//return g;

}


std::vector<int> PLM::prolong(GraphComm *g_initial, std::vector<int> coarsened_comm) {

	std::vector<int> init_comm = g_initial->communities;
	std::vector<int> new_comm(g_initial->n, 0);

	#pragma omp parallel for
	for (int i=0; i<g_initial->n; i++) {
		int i_comm = init_comm[i];
		new_comm[i] = coarsened_comm[i_comm];
	}
	return new_comm;
}
std::pair <int, float> max_pair_arg (std::pair <int, float> r, std::pair <int, float> n) {
        return (n.second >= r.second) ? n : r;
}



int PLM::ReturnCommunity(int i, GraphComm *g) {

	std::vector<pair<node_id, weight>> n_i = g->net[i];
	std::unordered_map<int,int> weights;
	std::unordered_map<int,int> volumes;
	int j, c=g->communities[i], c_n;	
	int comm_size = g->n;

	/* iterate once over all neighbors and compute weights from i to all communities.
	   Update the weights for all threads */
	for (auto neighbor_it = n_i.begin(); neighbor_it < n_i.end(); ++neighbor_it) {
		if ((int) neighbor_it->first != i) {
				c_n = g->communities[neighbor_it->first]; 
				volumes[c_n] = g->comm_volumes[c_n];
				weights[c_n] += neighbor_it->second;
		}
	}


	    volumes[c] = g->comm_volumes[c];
	    weight weight_c = weights[c];
	    weight volume_c = volumes[c] - g->volumes[i];
	    weight i_vol = g->volumes[i];
            float n_w_float = g->weight_net;
	    //long long double_sqr_n_w = g->weight_sq;
	    //float i_vol_divided = (float) i_vol / double_sqr_n_w;
	    float i_vol_divided = (float) i_vol/g->weight_net;
	    i_vol_divided = (float) i_vol_divided/(g->weight_sq);

	    float weight_c_divided = (float) weight_c / g->weight_net;

	    //if (i<5)
	    //	cout << "weight_float, i, i_vol, i_vol_divided, weight_c, weight_c_divided: " << n_w_float << ',' << i << ',' << i_vol << ',' <<  i_vol_divided << ',' << weight_c << ',' << weight_c_divided << endl;
		

	    std::pair<int, float> max_pair = std::make_pair(c, 0.0);
	    float a, b, dmod;
	
		for (auto c: weights) {
			a =  (c.second / n_w_float) - weight_c_divided;
			b = (volume_c - volumes[c.first]) * i_vol_divided;
			dmod = a + b;
			if (dmod > max_pair.second) {
				max_pair.first=c.first;
				max_pair.second=dmod;
			}

		}

	return max_pair.first;
}


std::unordered_map<int, int> PLM::Map_communities(GraphComm *g) {

	std::unordered_set<int> comms;
	std::unordered_map<int,int> com_map;
	for (int i=0; i<g->n; i++) {
		comms.insert(g->communities[i]);
	}
	// cout << "comms.size(): " << comms.size() << endl;

	int i = 0;
	for (auto it:comms) {
		com_map.insert(std::pair<int,int>(it, i));
		i++;
	}

	return com_map;

}


void  PLM::Local_move(GraphComm* graph) {


	int unstable = 1;
	int iterations=0;
	int updated = graph->n;
	int threads = std::min(graph->n, omp_get_max_threads());
	//cout << "Inside LM with " << threads << " threads, and " << graph->n << " nodes " << endl;


	/* */
	int threshold = graph->n * EPS;
	//cout << "threshold: " << threshold << endl;

	// when Local Move is called, each node alone is a community
	std::vector<int> volumes(graph->n, 0);	
	#pragma omp parallel for
	for (int j=0; j<graph->n; j++) 
		volumes[j] = graph->volumes[j];

	graph->comm_volumes = volumes;


	auto start = std::chrono::system_clock::now();	
	while (unstable && updated>threshold) {
		iterations++;
		//cout << "iteration:" << iterations << endl;
		//print(graph->communities);
		//cout << "----------------------------" << endl;
		unstable = 0;
		updated = 0;
		#pragma omp parallel for num_threads(threads)
		for (int i=0; i<graph->n; i++) {
			int i_comm = graph->communities[i];
			int new_comm = ReturnCommunity(i, graph);
			if (new_comm != i_comm) { 
				graph->communities[i] = new_comm;
				#pragma omp atomic write
				unstable=1;
				#pragma omp atomic update
				updated++;
				#pragma omp critical 
				{
					graph->comm_volumes[new_comm] += graph->volumes[i];
					graph->comm_volumes[i_comm] -= graph->volumes[i]; 
				
				}
			}

		}
		//cout << "updated: " << updated << endl;
		

	}
	auto end = std::chrono::system_clock::now();
	auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	//std::cout << "comm detection took time (in us) : " << total_time << std::endl;

	if (iterations > 1) {
		start = std::chrono::system_clock::now();
		std::unordered_map<int, int> com_map = Map_communities(graph);
		#pragma omp parallel for
		for (int i=0; i<(*graph).n; i++)
			(*graph).communities[i] = com_map[(*graph).communities[i]];
		end = std::chrono::system_clock::now();
		total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
		//std::cout << "map comm took time (in us) : " << total_time << std::endl;
	}
}


std::vector<int> PLM::Recursive_comm_detect(GraphComm *g) {


	std::vector<int> c_singleton(g->n, 0);
	# pragma omp parallel for schedule(static, NUM_SPLIT)
	for (int i=0; i<g->n; i++) {
		c_singleton[i] = i;
	}
	g->communities = c_singleton;
	Local_move(g);


	/*auto start = std::chrono::system_clock::now();
	GraphComm *g_new = coarsen(g);
	auto end = std::chrono::system_clock::now();

	auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "coarsen took time (in us) : " << total_time << std::endl;*/

    
	if (g->communities != c_singleton) {
		auto start = std::chrono::system_clock::now();
		GraphComm *g_new = new GraphComm;
		int n = *max_element(std::begin(g->communities), std::end(g->communities)) + 1;
		int threads=omp_get_max_threads();
		
		//vector<vector<std::unordered_map<int,int>>> array_thread(threads, vector<std::unordered_map<int,int>>(n)); 
		array_thread.resize(threads);
		for( auto &it : array_thread )
    		{	
        		it.resize(n);
    		}
		new_net_array.resize(n);	 

		coarsen(g, g_new, &array_thread, &new_net_array);
		auto end = std::chrono::system_clock::now();

        	auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        	//std::cout << "-------------coarsen took time (in us) : " << total_time << std::endl;

		start = std::chrono::system_clock::now();
	
		#pragma omp parallel for
		for (int j=0; j<n; j++) {
			for (int i=0; i < threads; i++)  
                		array_thread[i][j].clear();
		}

		#pragma omp parallel for
                for (int j=0; j<threads; j++) {
                 	array_thread[j].clear();
                }


		array_thread.clear();
		
		#pragma omp parallel for
		for (int i=0; i<n; i++) {
			new_net_array[i].clear();
		}
		new_net_array.clear();
		//new_net_array.resize(0);
		//new_net_array.shrink_to_fit();
		//vector<std::unordered_map<int,int>>().swap(new_net_array);


		end = std::chrono::system_clock::now();

        	total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        	//std::cout << "-------------cleaning took time (in us) : " << total_time << std::endl;
		
		std::vector<int> c_coarsened = Recursive_comm_detect(g_new);
		g->communities = prolong(g, c_coarsened);
	}
	return g->communities;

}

void get_weight_and_volumes(GraphComm* g) {
	long long sum = 0;
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

	get_weight_and_volumes(graph);
	//cout << "weight: " << graph->weight_net << endl;
	//graph->weight_sq = 2 * pow(graph->weight_net, 2); 
	graph->weight_sq=graph->weight_net*2;
	//cout << "weight_sq: " << graph->weight_sq << endl;
	graph->communities = Recursive_comm_detect(graph);
	//cout << "final communities: ";
	//print(graph->communities);
	//int found_comm = *std::max_element(std::begin(graph->communities), std::end(graph->communities)) + 1;
	//cout << "found " << found_comm << " communities" << endl;
}

// the same for LP and Louvain. TODO: defined once
void PLM::PrintCommunities(const std::string &file_name) {
	std::ofstream ofs;
	ofs.open(file_name, std::ios_base::out | std::ios_base::trunc);
	if (ofs.fail()) {
		std::cerr << "Error opening file " << file_name << std::endl;
		std::exit(2);
	}

	for (int i = 0; i < graph->n; i++) {
		ofs << graph->communities[i] << std::endl;
	}
}
