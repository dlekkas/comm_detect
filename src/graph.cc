#include "../include/graph.h"
#include "../include/network.h"


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

	// read neighours from the rest lines and populate adjacency list
	std::string tmp;
	std::getline(ifs, tmp);

	// unweighted graphs: adj_list is ok 
	/*
	while (std::getline(ifs, tmp)) {
		std::istringstream buf(tmp);
		std::vector<int> line { std::istream_iterator<int>(buf),
				std::istream_iterator<int>()};
		adj_list.push_back(line);
	}*/

	// weighted graphs: you need to keep weights
	while (std::getline(ifs, tmp)) {
        	std::istringstream buf(tmp);
        	std::vector<int> line { std::istream_iterator<int>(buf),
                              std::istream_iterator<int>()};
		std::vector<pair<node_id, weight>> v;
		if (weighted==1) {
			for (vector<int>::iterator it = line.begin(); it != line.end(); it+=2) {
		                int id = *it;
				int weight = *(it+1);
		                v.push_back(make_pair((node_id) (id-1), weight));
		        }
		}
		else {
			for (vector<int>::iterator it = line.begin(); it != line.end(); ++it) {
		                int id = *it;
		                v.push_back(make_pair((node_id) (id-1), 1));
		        }
		}

		net.push_back(v);
        }

}


void GraphComm::Net_init(const std::string &file_name) {
        std::ifstream ifs;
        ifs.open(file_name, std::ios_base::in);
        if (ifs.fail()) {
                std::cerr << "Error opening file " << file_name << std::endl;
                std::exit(2);
        }

        int weighted;
        //read number of vertices and edges from 1st line of file
        ifs >> n >> m >> weighted;
        //read neighours from the rest lines and populate adjacency list
        std::string tmp;
        std::getline(ifs, tmp);
        while (std::getline(ifs, tmp)) {
        	std::istringstream buf(tmp);
        	std::vector<int> line { std::istream_iterator<int>(buf),
                              std::istream_iterator<int>()};
		std::vector<pair<node_id, weight>> v;
		if (weighted==1) {
			for (vector<int>::iterator it = line.begin(); it != line.end(); it+=2) {
		                int id = *it;
				int weight = *(it+1);
		                v.push_back(make_pair((node_id) (id-1), weight));
		        }
		}
		else {
			for (vector<int>::iterator it = line.begin(); it != line.end(); ++it) {
		                int id = *it;
		                v.push_back(make_pair((node_id) (id-1), 1));
		        }
		}

		net.push_back(v);
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

network GraphComm::CreateNetwork() {
	network net;
	for (vector<vector<int>>::iterator it = adj_list.begin(); it != adj_list.end(); ++it) {
		vector<int> v_temp = *it;
		vector<pair<node_id, weight>> v;
		for (vector<int>::iterator it2 = v_temp.begin(); it2 != v_temp.end(); ++it2) {
				int id = *it2;
			v.push_back(make_pair((node_id) id, 1));
		}
        net.push_back(v);
	}
	return net;
}


void GraphComm::PrintNetwork() {
	for (int i = 0; i < n; i++) {
		std::cout << std::endl;
		std::cout << "Node[" << i << "] = ";
		vector<pair<node_id, weight>> my_neighbors = net[i]; 
		for (auto j = my_neighbors.begin(); j != my_neighbors.end(); j++) {
			std::cout << "node: " << j->first << ",  " << "weight: " << j->second << " | ";
		}
		std::cout << std::endl;
	}
}

// TODO(dimlek): implement function to write communities to file
// in order to evaluate the output file on modularity
void GraphComm::ProduceResult(const std::string &file_name) {
}
