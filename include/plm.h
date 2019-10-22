#ifndef PLM_H
#define PLM_H

#include "../include/graph.h"

#include <vector>
#include <unordered_map>
#include <iostream>

#define EPS 0.0001

class PLM {

	public:
		/* graph on which PLM will be applied */
		GraphComm graph;

		PLM(GraphComm &init_graph): graph(init_graph) {};

		~PLM() {};

		/* detect communities of graph */
		void DetectCommunities();
		
		/* do local moves for modularity gain */
		std::vector<int> Local_move(GraphComm graph, std::vector<int> communities); 

		/* Resulting communities are printed in the format specified by the
		 * DIMACS 10th challenge (i.e each line has the community number that
		 * corresponds to this node - line i has the community in which node i-1
		 * belongs) */
		void PrintCommunities(const std::string &file_name);

};


#endif 			// PLM_H
