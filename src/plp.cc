#include "../include/plp.h"

#include <omp.h>
#include <algorithm>

#define NUM_SPLIT 100

int PLP::random_id(int n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0,n);
	return dist(gen);
}

std::pair <int, int> max_with_arg (std::pair <int, int> r, std::pair <int, int> n) {
	// r is the already reduced value
	// n is the new value
	// fisrt of pair is the label and second the value
	return (n.second > r.second) ? n : r;
}

int PLP::dominant_label(int node) {
	//TODO: Replace the frequency array with an array of sum
	std::unordered_map<int,int> label_weights;

	std::vector<pair<node_id, weight>> neighbors = graph->net[node];
	// https://stackoverflow.com/a/17853547
	// parallel loops should be in the canonical form
	#pragma omp parallel for
	for (auto neighbor_it = neighbors.begin(); neighbor_it < neighbors.end(); ++neighbor_it)
	{
		#pragma omp atomic update
		label_weights[graph->communities[neighbor_it->first]]+=neighbor_it->second;
	}

	int max_val = 0, dom_label = graph->communities[node];
	std::pair<int, int> max_pair = std::make_pair(dom_label, max_val);
	for (auto bi: label_weights) {
		max_pair = max_with_arg(max_pair, bi);
	}

	/*
	#pragma omp parallel for shared(label_freq) reduction(max : max_label)
  	for (size_t b = 0; b < label_freq.bucket_count(); b++) {
		for (auto bi = label_freq.begin(b); bi != label_freq.end(b);bi++) {
			max_pair = max_with_arg(max_pair, *bi);
			max_label = max_pair.first;
		}
	*/

	// std::cout << "max: " << max_pair.second << " in: " << max_pair.first << std::endl;
	// std::unordered_map<int,int> label_freq2;

	// for (auto const& neighbor: graph->adj_list[node]) {
	// 	label_freq2[graph->communities[neighbor]]++;
	// }

	// int max_val2 = 0, dom_label2 = -1;
	// for(auto const& label2: label_freq2) {
	// 	if (max_val2 <= label2.second) {
	// 		max_val2 = label2.second;
	// 		dom_label2 = label2.first;
	// 	}
	// }
	// std::cout << "max: " << max_val2 << " in: " << dom_label2 << ", #2" << std::endl << std::endl;

	return max_pair.first;
	// return max_label;
}


void PLP::DetectCommunities() {
	graph->communities.reserve(graph->n);

	// TODO see https://stackoverflow.com/questions/13224155/how-does-the-omp-ordered-clause-work
	// play with the number of static scheduling
	# pragma omp parallel for schedule(static, NUM_SPLIT)
	for (int i = 0; i < graph->n; i++) {
		graph->communities[i] = i;
	}

	// // print to see if communities have been set in order
	// std::cout << "myvector contains:";
  	// std::vector<int> myvector = graph->communities;
	// for (auto it = myvector.begin() ; it != myvector.end(); ++it)
    // 	std::cout << ' ' << *it;m
  	// std::cout << '\n';

	int updated = graph->n;
	int updated_previous = 0;
	// while (updated > threshold && updated != updated_previous) {
	while ((updated > threshold) && (abs(updated_previous - updated) > threshold)) {
		updated_previous = updated;
		updated = 0;
		#pragma omp parallel for
		for (int i = 0; i < graph->n; i++) {
			int new_label = dominant_label(i);
			if (new_label != graph->communities[i]) {
				graph->communities[i] = new_label;
				updated++;
			}

		}
		cout << "updated_previous: " << updated_previous << endl;
		cout << "updated: " << updated << endl;
	}
	// cout << "final communities:" << endl;
	// for (int i = 0; i < graph->n; i++) {
	// 	std::cout << graph->communities[i] << ' ';
	// }
	// cout << endl;

}

std::map<int, int> PLP::Map_communities(GraphComm *g) {

	std::vector<int> comms;
	std::map<int,int> com_map;
	for (int i=0; i<g->n; i++) {
		int c = g->communities[i];

		if (std::find(comms.begin(), comms.end(), c) == comms.end()) {
			comms.push_back(c);
		}
	}
	sort(comms.begin(), comms.end());

	for (int i=0; i < (int) comms.size(); i++)
		com_map.insert(std::pair<int,int>(comms[i], i));

	return com_map;

}

void PLP::PrintCommunities(const std::string &file_name) {
	std::ofstream ofs;
	ofs.open(file_name, std::ios_base::out | std::ios_base::trunc);
	if (ofs.fail()) {
		std::cerr << "Error opening file " << file_name << std::endl;
		std::exit(2);
	}

	for (int i = 0; i < graph->n; i++) {
		ofs << graph->communities[i] << std::endl;
	}
}
