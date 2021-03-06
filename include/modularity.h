#ifndef MODULARITY_H
#define MODULARITY_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <unordered_map>
#include "./network.h"

#include <omp.h>

#define NUM_SPLIT 100

using namespace std;

inline bool belongs_to_community(node_id u, community c) { return find(c.begin(), c.end(), u) != c.end(); }

inline weight weight_from_node_to_community(node_id u, community c, network net) {
  // 2 possible ways:
  // a) iterate on neighbors of u and check if they belong to community and
  // then add the weight of their edge
  // b) iterate on community members and check if they are connected with u and
  // then add the weight of their edge

  // TODO Discuss a or b. (Trying a)
  weight sum = 0;

  // TODO check if u < n
  vector<pair<node_id, weight>> neighbors = net[u];
  // search if the neighbors of u belong to community; if so add the weight of
  // the edge to the total sum
  for (auto it = neighbors.begin(); it < neighbors.end(); ++it) {
    if (belongs_to_community(it->first, c)) {
      sum += it->second;
    }
  }

  return sum;
}

inline weight weight_of_community(community c, network net) {
  weight sum = 0;

  for (community::iterator it = c.begin(); it != c.end(); ++it) {
    sum += weight_from_node_to_community(*it, c, net);
  }

  return sum;
}

inline weight weight_of_network(network net) {
  weight sum = 0;

  for (auto it = net.begin(); it < net.end(); ++it) {
    vector<pair<node_id, weight>> neighbors = *it;
    // add the weight of the neighbors to the total sum
    for (vector<pair<node_id, weight>>::iterator j = neighbors.begin(); j != neighbors.end(); ++j) {
      sum += j->second;
    }
  }

  return sum;
}

inline weight volume_of_node(node_id u, network net) {
  weight sum = 0;

  // TODO check if u < n
  vector<pair<node_id, weight>> neighbors = net[u];
  for (auto it = neighbors.begin(); it < neighbors.end(); ++it) {
    sum += it->second;
    // edge to itself must be counted twice
    if (it->first == u)
    sum += it->second;
  }

  return sum;
}

inline std::vector<int> compute_node_volumes(int n, network net) {
    std::vector<int> volumes;
    for (int i=0; i<n; i++) {
        int v = volume_of_node(i, net);
        {
            volumes.push_back(v);
        }
    }
    return volumes;

}

inline weight volume_of_community(community c, network net) {
  weight sum = 0;

  for (community::iterator it = c.begin(); it != c.end(); ++it) {
    sum += volume_of_node(*it, net);
  }

  return sum;
}

inline weight volume_of_community_plm(community c, std::vector<int> volumes) {
  weight sum = 0;

  for (auto it = c.begin(); it < c.end(); ++it) {
    sum += volumes[*it];
  }

  return sum;
}

inline communities communities_from_file(const std::string &file_name) {
    std::ifstream ifs;
    ifs.open(file_name, std::ios_base::in);
    if (ifs.fail()) {
        std::cerr << "Error opening file " << file_name << std::endl;
        std::exit(2);
    }

  std::unordered_map<node_id, vector<node_id>> comms_map;
    std::string tmp;

  // fill a map that contains the participants of each community
  // the key of the map is the label of the community
  node_id i = 0;
  while (std::getline(ifs, tmp)) {
    node_id c = stoi(tmp);
        comms_map[c].push_back(i);
    i++;
    }

    communities comm_list;

  for (auto it: comms_map) {
    vector<node_id> participants_of_community;
    for (auto i: it.second) {
      participants_of_community.push_back(i);
    }
    comm_list.push_back(participants_of_community);
  }

  return comm_list;
}

inline void print_communities(communities cs) {

    for (auto i = cs.begin(); i != cs.end(); ++i) {
        // for each node print its neighbors
        std::cout << "Community[" << i - cs.begin() << "] = ";
        for (auto j = i->begin(); j != i->end(); j++) {
            std::cout << *j << " ";
        }
        std::cout << std::endl;
    }
}

inline void print_network(network net) {
  for (auto i = net.begin(); i != net.end(); ++i) {
    // for each node print its neighbors
    auto pos = i - net.begin() ;
    std::cout << "Node[" << pos << "] = ";
    for (auto j = i->begin(); j != i->end(); j++) {
            std::cout << "(" << j->first << ", " << j->second << ") ";
    }
    std::cout << endl;
  }

}

// TODO add neighbors_of_node function

inline modularity compute_modularity(communities z, network net) {
  modularity mod = 0;

  weight weight_net = weight_of_network(net); // TODO: this doesn't change, maybe compute it only once
  // cout << "weight_net: " << weight_net << endl;


  for (communities::iterator c = z.begin(); c != z.end(); ++c) {
    weight vol_c = volume_of_community(*c, net);
    // TODO move weight_net in denominator out of the sum
    // TODO compute separately the two components of the sum
    weight weight_comm = weight_of_community(*c, net);
    mod += (1.0 * weight_comm / weight_net) - ((1.0 * vol_c * vol_c) / (4 * weight_net * weight_net));
  }

  return mod;
}


inline modularity compute_modularity_from_node_comm(vector<int> communities, int n, network net) {
  	
  std::vector<int> comms;

  for (int i = 0; i < n; i++) {
		int c = communities[i];

		if (std::find(comms.begin(), comms.end(), c) == comms.end()) {
			comms.push_back(c);
		}
	}
	sort(comms.begin(), comms.end());

	std::vector<std::vector<int>> cs;
	for (int c = 0; c < ((int) comms.size()); c++) {
		std::vector<int> v;
		for (int i = 0; i < n; i++)
			if (c == communities[i])
				v.push_back(i);
		cs.push_back(v);
	}
	
	// cout << "comms.size(): " << comms.size() << endl;
	
	// cout << "communities for each node:" << endl;
	// for (int i = 0; i < ((int) communities.size()); i++) {
	// 	std::cout << communities.at(i) << " ";
	// }
	// cout << endl;
	
	// for (auto it = cs.begin(); it < cs.end(); it++) {
	// 	for (auto j = it->begin(); j < it->end(); j++)
	// 		cout << *j << " ";
	// 	cout << endl;
	// }

  return compute_modularity(cs, net);
}

inline std::vector<weight> compute_weights(node_id i, int c, int d, GraphComm g) {

    std::vector<weight> results(2, 0);
    network net = g.net;
    vector<pair<node_id, weight>> neighbors = net[i];

    for (auto it = neighbors.begin(); it < neighbors.end(); ++it) {
        if (it->first != i) {
            if (g.communities[it->first] == c)  {
                /* .... */
                results[0] += it->second;
            }
            if (g.communities[it->first] == d) {
                /* .... */
                results[1] += it->second;
            }
        }
    }

    return results;
}

inline modularity compute_modularity_difference(node_id i, node_id n, GraphComm g) {

    weight weight_net = g.weight_net;
    weight i_vol = g.volumes[i];
    int c = g.communities[i];
    int d = g.communities[n];

    weight volume_c = 0;
    weight volume_d = 0;

    for (int j=0; j<g.n; j++)
        if (j != (int) i) {
            if (g.communities[j] == c)
                volume_c += g.volumes[j];
            else if (g.communities[j] == d)
                volume_d += g.volumes[j];
        }

    std::vector<weight> weight_results = compute_weights(i, c, d, g);

    weight weight_c = weight_results[0];
    weight weight_d = weight_results[1];

    float a =  ((1.0 * (weight_d - weight_c)) / weight_net);
    float b = (1.0 * (volume_c - volume_d) * i_vol) / (2 * weight_net * weight_net);

    modularity mod_diff = a+b;

    return mod_diff;
}

inline std::vector<weight> compute_weights_array(node_id i, int c, int d, GraphComm g) {

    std::vector<weight> results(2, 0);
    int *net = g.adj_matrix;
    node_id n = (node_id) g.n;

    int *i_net = &(net[i * n]);

    for (node_id u = 0; u < n; u++) {
        // check if they are neighbors
        int weight_i_u = i_net[u];
        if (weight_i_u != 0) {
            if (u != i) {
                if (g.communities_array[u] == c)
                    results[0] += weight_i_u;
                if (g.communities_array[u] == d)
                    results[1] += weight_i_u;
            }
        }
    }

    return results;
}

inline modularity compute_modularity_difference_array(node_id i, node_id n, GraphComm g) {

    weight weight_net = g.weight_net;
    weight i_vol = g.volumes_array[i];
    int c = g.communities_array[i];
    int d = g.communities_array[n];

    weight volume_c = 0;
    weight volume_d = 0;

    for (int j=0; j<g.n; j++)
        if (j != (int) i) {
            if (g.communities_array[j] == c)
                volume_c += g.volumes_array[j];
            else if (g.communities_array[j] == d)
                volume_d += g.volumes_array[j];
        }

    std::vector<weight> weight_results = compute_weights_array(i, c, d, g);

    weight weight_c = weight_results[0];
    weight weight_d = weight_results[1];

    float a =  ((1.0 * (weight_d - weight_c)) / weight_net);
    float b = (1.0 * (volume_c - volume_d) * i_vol) / (2 * weight_net * weight_net);

    modularity mod_diff = a+b;

    return mod_diff;
}

#endif 			// MODULARITY_H
