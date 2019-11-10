#include "../include/graph.h"
#include "../include/plp.h"

#include <omp.h>
#include <vector>
#include <iostream>
#include <chrono>

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
			std::cout << "Parallel PLP with " << omp_get_num_threads() << " threads." << std::endl;
		}
	}

	/* initialize graph from file */
	GraphComm test_g;
	test_g.Init(file_name);

	/* detect communities of graph */
	PLP test_plp { test_g };

	/* benchmark time of community detection */
	auto start = std::chrono::system_clock::now();
	test_plp.DetectCommunities();
	auto end = std::chrono::system_clock::now();

	auto total_time = std::chrono::duration_cast<
			std::chrono::microseconds>(end - start).count();
	cout << "Detect Communities time (in ms): " << total_time / 1000.0 << endl;	

	/* print result to file */
	test_plp.PrintCommunities("comm.out");

	return 0;
}
