#include <assert.h>

#include <vector>
#include <unordered_map>
#include <iostream>
#include <random>


#define EPS 0.0001


int random_id(int n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0,n);
	return dist(gen);
}


int dominant_label(const std::vector<int> &neighbors,
			const std::vector<int> &communities) {

	std::unordered_map<int,int> label_freq;
	for (auto const& node: neighbors) {
		label_freq[communities[node]]++;
	}

	int max_val = 0, dom_label = -1;
	for(auto const& label: label_freq) {
		if (max_val < label.second) {
			max_val = label.second;
			dom_label = label.first;
		}
	}

	return dom_label;
}




int main() {
	int n = 200;
	std::vector<int> communities;
	std::vector<std::vector<int>> adj_list;

	for (int i = 0; i < n; i++) {
		communities[i] = random_id(n);
	}

	/*
	int active_nodes[n];
	for (int i = 0; i < n; i++) {
		active_nodes[i] = 1;
	}
	*/

	int updated = n;
	int threshold = n * EPS;
	while (updated > threshold) {
		updated = 0;

		for (int i = 0; i < n; i++) {
			int new_label = dominant_label(adj_list[i], communities);
			if (new_label != communities[i]) {
				communities[i] = new_label;
				updated++;
			}

		}
	}

	return 0;
}


