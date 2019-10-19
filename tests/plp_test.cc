#include "../include/graph.h"
#include "../include/plp.h"

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
	test_g.PrintGraph();

	/* detect communities of graph */
	PLP test_plp { test_g };
	test_plp.DetectCommunities();

	/* print the result to stdout */
	test_g.ProduceResult("dummy.out");

	return 0;
}
