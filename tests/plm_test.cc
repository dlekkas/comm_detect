#include "../include/graph.h"
#include "../include/plm.h"

#include <vector>


int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <input-file>" << std::endl;
		std::exit(1);
	}
	std::string file_name = argv[1];

	/* initialize graph based on file and confirm correct parsing */
	GraphComm test_g;
	test_g.Init(file_name);
	//test_g.PrintGraph();

	/* detect communities of graph */
	PLM test_plm { test_g };
	test_plm.DetectCommunities();

	/* print the result to file */
	test_plm.PrintCommunities("comm_plm.out");

	return 0;
}
