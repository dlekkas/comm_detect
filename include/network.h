#ifndef NETWORK_H
#define NETWORK_H

#include <iostream>
#include <vector>

using namespace std;

// type aliases
using modularity = double;
using node_id = int;
using weight = int;
using network = vector<vector<pair<node_id, weight>>>;
using community = vector<node_id>;
using communities = vector<community>;

#endif 			// NETWORK_H
