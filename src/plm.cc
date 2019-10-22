#include "../include/plm.h"
#include "../include/graph.h"

#include <map>
 

std::vector<int> PLM::Local_move(GraphComm graph, std::vector<int> communities) {
	int unstable = 1;

	while (unstable) {
		unstable = 0;
		for (int i; i<graph.n; i++) {
			// create a map <neighbor_id, mod_diff>
			int i_comm = communities[i];
			int n_id = 0, max_diff = -1;
			std::map<int,float> mod_map;
			for (int neighbor: graph.adj_list[i]) {
				int n_z = communities[neighbor];
				// for each v, neighbor, compute mod_diff
				float mod_diff = 0.0; //TODO: call function to compute difference on modularity (paper 3.1)
									 //when i is moved from i_comm to n_z
				mod_map[i] = mod_diff;
			// choose the <id, d> = <neighbor_id, mod_diff> with the largest mod_diff
			}
			for(std::map<int, float>::iterator it=mod_map.begin(); it!=mod_map.end(); ++it) {
				if (max_diff <= it -> second) {
					max_diff = it -> second;
					n_id = it -> first;
				}
			}
			int z = communities[n_id];
			if (max_diff>0) {
				communities[i] = z;
				unstable=1;
			}
	
		}
	}

	return communities;
}

