#include "../include/graph.h"
#include "../include/plm.h"

#include <vector>
#include <sys/time.h>
#include <omp.h>
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
			std::cout << "algo: PLM, threads: " << omp_get_num_threads() << ", ";
		}
	}

	/* initialize graph based on file and confirm correct parsing */
	GraphComm test_g;
	test_g.Net_init(file_name);

	/* detect communities of graph */
	PLM test_plm { test_g };

        /* benchmark time of community detection */
	auto start = std::chrono::system_clock::now();
	test_plm.DetectCommunities();
	auto end = std::chrono::system_clock::now();

	auto total_time = std::chrono::duration_cast<
			std::chrono::milliseconds>(end - start).count();
	std::cout << "time (in sec) : " << total_time / 1000.0 << std::endl;

	return 0;
}
