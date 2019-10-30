#include "../include/graph.h"
#include "../include/modularity.h"
#include "../include/network.h"

#include <vector>


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <graph-file> <communities-file>" << std::endl;
        std::exit(1);
    }

    std::string graph_file = argv[1], communities_file = argv[2];

    GraphComm g;
    g.Init(graph_file);
    vector<vector<int>> neighbors = g.adj_list;

    // transform graph to network - TODO: add this to graph file - change here
    network net;
    for (vector<vector<int>>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
        vector<int> v_temp = *it;
        vector<pair<node_id, weight>> v;
        for (vector<int>::iterator it2 = v_temp.begin(); it2 != v_temp.end(); ++it2) {
            int id = *it2;
            v.push_back(make_pair((node_id) id, 1));
        }
        net.push_back(v);
    }

    print_network(net);

    communities cs;
    cs = communities_from_file(communities_file);
    print_communities(cs);

    modularity mod;
    mod = compute_modularity(cs, net);
    cout << "Modularity: " << mod << endl;

    return 0;
}
