#include "../include/graph.h"

#include <vector>


int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <input-file>" << std::endl;
		std::exit(1);
	}

	std::string file_name = argv[1];

	GraphComm test_g;
	test_g.Init(file_name);
	test_g.PrintGraph();

	return 0;
}
