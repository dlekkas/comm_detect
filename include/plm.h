#ifndef PLM_H
#define PLM_H

#include "../include/graph.h"
#include "../include/network.h"

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

		/* used by the 'coarsen' step of PLM: checks if 2 communities should be connected at the coarsened graph */		
		int connected(int comm_1, int comm_2, std::vector<int> communities, std::vector<std::vector<int>> adj_list);

		/* refresh the communities of the nodes of the input graph according to the new communities found by the coarsened step */
		std::vector<int> prolong(GraphComm g_initial, std::vector<int> coarsened_comm);

		/* turn each community of the graph taken as argument into a supernode and create a new graph */
		GraphComm coarsen(GraphComm* g_initial, std::vector<int> comm);


		/* do local moves for modularity gain */
		std::vector<int> Local_move(GraphComm graph, std::vector<int> communities); 

		/* */
		std::vector<int> Recursive_comm_detect(GraphComm g);

		/* */
		int ReturnCommunity(int node, int i_comm, community i_comm_vector, std::vector<int> adj_list_i, std::vector<int> communities, network net, weight w, std::vector<int> v);
		

		/* Resulting communities are printed in the format specified by the
		 * DIMACS 10th challenge (i.e each line has the community number that
		 * corresponds to this node - line i has the community in which node i-1
		 * belongs) */
		void PrintCommunities(const std::string &file_name);


		community get_community_vector(std::vector<int> communities, int comm);


};


#endif 			// PLM_H
