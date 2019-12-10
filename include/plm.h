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
		GraphComm coarsen(GraphComm* g_initial);

		/* do local moves for modularity gain */
		void Local_move(GraphComm* graph);

		/* */
		std::vector<int> Recursive_comm_detect(GraphComm g);

		/* */
		std::pair<int,float> ReturnCommunity(int i, GraphComm* g);

		/* */
		std::map<int, int> Map_communities(GraphComm g);

		/* Resulting communities are printed in the format specified by the
		 * DIMACS 10th challenge (i.e each line has the community number that
		 * corresponds to this node - line i has the community in which node i-1
		 * belongs) */
		void PrintCommunities(const std::string &file_name);

		community get_community_vector(std::vector<int> communities, int comm);

};


class PLM_MPI{

    public:
        /* graph on which PLM_MPI will be applied */
        GraphComm graph;

        PLM_MPI(GraphComm &init_graph): graph(init_graph) {};

        ~PLM_MPI() {};

        void DetectCommunities(int world_rank, int world_size);

		/* turn each community of the graph taken as argument into a supernode and create a new graph */
		GraphComm coarsen(GraphComm *g_initial, int world_rank, int world_size);

        int *Recursive_comm_detect(GraphComm g, int world_rank, int world_size, int recursions);

        /* do local moves for modularity gain */
		void Local_move(GraphComm *graph, int world_rank, int world_size);

		/* refresh the communities of the nodes of the input graph according to the new communities found by the coarsened step */
		int *prolong(GraphComm g_initial, int *coarsened_comm);

		//TODO: move it to the GraphComm class
		int *GetAdjacencyMatrix(GraphComm* g);

		int ReturnCommunity(int i, GraphComm g, int comm_size, int world_rank, int world_size);

		/* */
		std::map<int, int> Map_communities(GraphComm g);

		/* Resulting communities are printed in the format specified by the
		 * DIMACS 10th challenge (i.e each line has the community number that
		 * corresponds to this node - line i has the community in which node i-1
		 * belongs) */
		void PrintCommunities(const std::string &file_name);

};


#endif 			// PLM_H
