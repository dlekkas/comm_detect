#include "../include/graph.h"
#include "../include/plp.h"
#include "../include/modularity.h"

#include <omp.h>
#include <vector>
#include <iostream>
#include <chrono>
#include <unordered_set>

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <input-file>" << std::endl;
		std::exit(1);
	}
	std::string file_name = argv[1];

	int threads = omp_get_max_threads();

	#pragma omp parallel num_threads(threads)
	{
		/* Obtain thread number */
		int tid = omp_get_thread_num();

		if (tid == 0) {
			std::cout << "algo: PLP, threads: " << omp_get_num_threads() << ", ";
		}
	}

	/* initialize graph from file */
	GraphComm test_g;
	auto start = std::chrono::system_clock::now();
	test_g.Init(file_name);

	auto end = std::chrono::system_clock::now();
	auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	//cout << "Read file succesfully after " << total_time / 1000.0 << "s" << std::endl;
	// test_g.PrintNetwork();

	/* detect communities of graph */
	PLP test_plp(&test_g);

	/* benchmark time of community detection */
	start = std::chrono::system_clock::now();
	test_plp.DetectCommunities();
	end = std::chrono::system_clock::now();

	total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "time (in sec): " << total_time / 1000.0;

	/* print result to file */
	//test_plp.PrintCommunities("comm.out");
	// for (int i = 0; i < test_plp.graph.n; i++) {
	// 	cout << test_plp.graph.communities[i] << endl;
	// }

	cout << endl << "map communities . . ." <<  endl;

	// std::map<int, int> com_map = test_plp.Map_communities(test_plp.graph);
  	// std::unordered_set<int> cs;
	// int n = test_plp.graph->n;

	// for (int i = 0; i < test_plp.graph->n; i++) {
	// 	// cout << test_plp.graph.communities[i] << "->" << com_map[test_plp.graph.communities[i]] << endl;
	// 	cs.insert(com_map[test_plp.graph->communities[i]]);
	// }
	// cout << "set size: " << cs.size() << endl;
	// for (int i = 0; i < n; i++) {
	// 	cout << cs[i] << endl;
	// }

	// modularity mod;
  	// mod = compute_modularity_from_node_comm(cs, n,
	// 										test_g.net);
    // cout << ", modularity: " << mod << endl;

	return 0;
}
