#ifndef MODULARITY_H
#define MODULARITY_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include "./network.h"

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
  for (vector<pair<node_id, weight>>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
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

  for (network::iterator it = net.begin(); it != net.end(); ++it) {
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
  for (vector<pair<node_id, weight>>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
    sum += it->second;
    // edge to itself must be counted twice
    if (it->first == u) sum += it->second;
  }

  return sum;
}

weight volume_of_community(community c, network net) {
  weight sum = 0;

  for (community::iterator it = c.begin(); it != c.end(); ++it) {
    sum += volume_of_node(*it, net);
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

	// read number of communities
	int n;
  ifs >> n;

  communities comm_list;
	std::string tmp;
  std::getline(ifs, tmp);
	while (std::getline(ifs, tmp)) {
		std::istringstream buf(tmp);
		std::vector<node_id> line { std::istream_iterator<int>(buf),
				std::istream_iterator<int>()};
		comm_list.push_back(line);
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

  weight weight_net = weight_of_network(net);
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

#endif 			// MODULARITY_H
