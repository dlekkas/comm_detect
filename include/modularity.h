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

bool belongs_to_community(node_id u, community c) { return find(c.begin(), c.end(), u) != c.end(); }

weight weight_from_node_to_community(node_id u, community c, network net) {
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

weight weight_of_community(community c, network net) {
  weight sum = 0;

  for (community::iterator it = c.begin(); it != c.end(); ++it) {
    sum += weight_from_node_to_community(*it, c, net);
  }

  return sum;
}

weight weight_of_network(network net) {
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

weight volume_of_node(node_id u, network net) {
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

std::vector<int> compute_node_volumes(int n, network net) {
	std::vector<int> volumes;
	for (int i=0; i<n; i++) {
		int v = volume_of_node(i, net);
		{
			volumes.push_back(v);
		}
	}
	return volumes;

}

weight volume_of_community(community c, network net) {
  weight sum = 0;

  for (community::iterator it = c.begin(); it != c.end(); ++it) {
    sum += volume_of_node(*it, net);
  }

  return sum;
}

weight volume_of_community_plm(community c, std::vector<int> volumes) {
  weight sum = 0;

  for (auto it = c.begin(); it < c.end(); ++it) {  
    sum += volumes[*it];
  }

  return sum;
}

communities communities_from_file(const std::string &file_name) {
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

void print_communities(communities cs) {

  for (auto i = cs.begin(); i != cs.end(); ++i) {
		// for each node print its neighbors
		std::cout << "Community[" << i - cs.begin() << "] = ";
		for (auto j = i->begin(); j != i->end(); j++) {
			std::cout << *j << " ";
		}
		std::cout << std::endl;
	}
}

void print_network(network net) {
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

modularity compute_modularity(communities z, network net) {
  modularity mod = 0;

  weight weight_net = weight_of_network(net); // TODO: this doesn't change, maybe compute it only once
  cout << "weight_net: " << weight_net << endl;

  
  for (communities::iterator c = z.begin(); c != z.end(); ++c) {
    weight vol_c = volume_of_community(*c, net);
    // TODO move weight_net in denominator out of the sum
    // TODO compute separately the two components of the sum
    weight weight_comm = weight_of_community(*c, net);
    mod += (1.0 * weight_comm / weight_net) - ((1.0 * vol_c * vol_c) / (4 * weight_net * weight_net));
  }

  return mod;
}


std::vector<weight> compute_all_weights(node_id i, int c, int d, GraphComm g) {

   std::vector<weight> results(4, 0);
   network net = g.net;
   vector<pair<node_id, weight>> neighbors = net[i];

   for (auto it = neighbors.begin(); it < neighbors.end(); ++it) {
	if (it->first != i) {
		if (g.communities[it->first] == c)  {
			/* .... */
			results[0] += it->second;
			results[1] += g.volumes[it->first];
		}
		if (g.communities[it->first] == d) {
			/* .... */
			results[2] += it->second;
			results[3] += g.volumes[it->first];
		}
	}
   }	
   return results;

}

modularity compute_modularity_difference(node_id i, node_id n, GraphComm g) {

    weight weight_net = g.weight_net; 
    weight n_vol = g.volumes[n];

    int c = g.communities[i];
    int d = g.communities[n];

    std::vector<weight> weight_results = compute_all_weights(i, c, d, g);

    weight weight_c = weight_results[0];
    weight weight_d = weight_results[2];
    weight volume_c = weight_results[1];
    weight volume_d = weight_results[3];

    float a =  ((1.0 * (weight_d - weight_c)) / weight_net);
    float b = (1.0 * (volume_c - volume_d) * n_vol) / (2 * weight_net * weight_net);

    modularity mod_diff = a+b; 

    return mod_diff;
}

#endif 			// MODULARITY_H
