#include "../include/graph.h"
#include "../include/plp.h"

#include <omp.h>
#include <vector>
#include <iostream>


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

	/* initialize graph based on file and confirm correct parsing */
	GraphComm test_g;
	test_g.Init(file_name);
	test_g.PrintGraph();

	/* detect communities of graph */
	PLP test_plp { test_g };
	test_plp.DetectCommunities();

	/* print the result to file */
	test_plp.PrintCommunities("comm.out");

	return 0;
}
