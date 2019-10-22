#include "../include/plp.h"


int PLP::random_id(int n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0,n);
	return dist(gen);
}


int PLP::dominant_label(int node) {
	std::unordered_map<int,int> label_freq;
	for (auto const& neighbor: graph.adj_list[node]) {
		label_freq[graph.communities[neighbor]]++;
	}

	int max_val = 0, dom_label = -1;
	for(auto const& label: label_freq) {
		if (max_val <= label.second) {
			max_val = label.second;
			dom_label = label.first;
		}
	}

	return dom_label;
}


void PLP::DetectCommunities() {
	for (int i = 0; i < graph.n; i++) {
		graph.communities.push_back(i);
	}

	int updated = graph.n;
	while (updated > threshold) {
		updated = 0;
		for (int i = 0; i < graph.n; i++) {
			int new_label = dominant_label(i);
			if (new_label != graph.communities[i]) {
				graph.communities[i] = new_label;
				updated++;
			}

		}
	}

}

void PLP::PrintCommunities(const std::string &file_name) {
	std::ofstream ofs;
	ofs.open(file_name, std::ios_base::out | std::ios_base::trunc);
	if (ofs.fail()) {
		std::cerr << "Error opening file " << file_name << std::endl;
		std::exit(2);
	}

	for (int i = 0; i < graph.n; i++) {
		ofs << graph.communities[i] << std::endl;
	}
}


