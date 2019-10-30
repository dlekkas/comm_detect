#ifndef GRAPH_H
#define GRAPH_H

#include <assert.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include "../include/network.h"

class GraphComm {

	public:
		GraphComm(): n(0), m(0) {};

		~GraphComm() {}

		/* initialize graph based on a given input file */
		void Init(const std::string &in_file);

		/* write community assigned to each node in a file */
		void ProduceResult(const std::string &out_file);

		/* print the adjacency list (debugging purposes) */
		void PrintGraph();

		/* number of vertices */
		int n;
		/* number of edges */
		int m;

		/* adjacency list of graph */
		std::vector<std::vector<int>> adj_list;

		/* community of each node */
		std::vector<int> communities;

		/* turn graph into a network */
		network CreateNetwork();

		network net;


};



#endif 			// GRAPH_H
