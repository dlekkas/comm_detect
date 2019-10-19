#include "../include/graph.h"

void GraphComm::Init(const std::string &file_name) {
	std::ifstream ifs;
	ifs.open(file_name, std::ios_base::in);
	if (ifs.fail()) {
		std::cerr << "Error opening file " << file_name << std::endl;
		std::exit(2);
	}

	int weighted;
	// read number of vertices and edges from 1st line of file
	ifs >> n >> m >> weighted;
	std::cout << "weighted: " << weighted << std::endl;

	// read neighours from the rest lines and populate adjacency list
	std::string tmp;
	std::getline(ifs, tmp);
	while (std::getline(ifs, tmp)) {
		std::istringstream buf(tmp);
		std::vector<int> line { std::istream_iterator<int>(buf),
				std::istream_iterator<int>()};
		adj_list.push_back(line);
	}
}


void GraphComm::PrintGraph() {
	for (int i = 0; i < n; i++) {
		// for each node print its neighbors
		std::cout << "Node[" << i << "] = ";
		for (auto j = adj_list[i].begin(); j != adj_list[i].end(); j++) {
			std::cout << *j << " ";
		}
		std::cout << std::endl;
	}
}


// TODO(dimlek): implement function to write communities to file
// in order to evaluate the output file on modularity
void GraphComm::ProduceResult(const std::string &file_name) {
}
