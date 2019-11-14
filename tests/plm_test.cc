#include "../include/graph.h"
#include "../include/plm.h"

#include <vector>
#include <sys/time.h>
#include <omp.h>

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <input-file>" << std::endl;
		std::exit(1);
	}
	std::string file_name = argv[1];

	/* initialize graph based on file and confirm correct parsing */
	GraphComm test_g;
	test_g.Net_init(file_name);
	//test_g.PrintGraph();

	omp_set_num_threads(1);
	/* detect communities of graph */
	PLM test_plm { test_g };
	struct timeval start, end;

        gettimeofday(&start, NULL);

        test_plm.DetectCommunities();

        gettimeofday(&end, NULL);
        float duration = (1.0 * end.tv_sec * 1000 + (1.0 * end.tv_usec) / 1000) - (1.0 * start.tv_sec * 1000 + (1.0 * start.tv_usec) / 1000);

        cout << "Detect Communities time (in ms): " << duration << endl;


	/* print the result to file */
	test_plm.PrintCommunities("comm_plm.out");

	return 0;
}
