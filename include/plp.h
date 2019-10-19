#ifndef PLP_H
#define PLP_H

#include "../include/graph.h"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <random>

#define EPS 0.0001

class PLP {

	public:
		/* graph on which PLP will be applied */
		GraphComm graph;

		PLP(const GraphComm &init_graph): graph(init_graph),
				threshold(init_graph.n * EPS) {};

		~PLP() {};

		/* detect communities of graph */
		void DetectCommunities();

	private:
		/* number of node labels updated at a single iteration until we stop */
		int threshold;

		/* find the community which is dominating among this node's neighbours */
		int dominant_label(int node);

		/* get a truly random number from 0 to n-1 */
		int random_id(int n);

};


#endif 			// PLP_H
