#include "../include/graph.h"
#include "../include/plm.h"

#include <mpi.h>
#include <vector>
#include <sys/time.h>
#include <chrono>

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
    int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {

		if (argc != 2) {
			std::cerr << "Usage: " << argv[0] << " <input-file>" << std::endl;
			std::exit(1);
		}

		// std::cout << "algo: PLM_MPI, threads: " << size << ", ";

	}

	/* initialize graph based on file and confirm correct parsing */
	GraphComm test_g;
	if (rank == 0) {
		std::string file_name = argv[1];
		test_g.Net_init(file_name);
	}
	PLM_MPI test_plm { test_g };

	/* benchmark time of community detection */
	auto start = std::chrono::system_clock::now();

	test_plm.DetectCommunities(rank, size);

	auto end = std::chrono::system_clock::now();
	auto total_time = std::chrono::duration_cast<
			std::chrono::milliseconds>(end - start).count();

	if (rank == 0) {
		std::cout << "time (in sec) : " << total_time / 1000.0 << std::endl;
	}

	MPI_Finalize();
	return 0;
}
