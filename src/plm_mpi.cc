#include "../include/plm.h"
#include "../include/graph.h"
#include "../include/modularity.h"

#include <mpi.h>
#include <map>

#define NUM_SPLIT 100

void print(std::vector<int> const &input)
{
    for (int i = 0; i <(int)(input.size()); i++) {
        std::cout << input.at(i) << ' ';
    }
    std::cout << '\n';
}


void print_map(std::map<int,std::vector<int>> myMap) {
for(map<int, std::vector<int>>::const_iterator it = myMap.begin(); it != myMap.end(); ++it)
{
    std::cout << "key:" << it->first << endl;
    print(it-> second);
}
}

GraphComm PLM::coarsen(GraphComm* g_initial) {


    GraphComm g;
    std::vector<int> comm = (*g_initial).communities;

    std::vector<int> g_vertexes;
    int i, j;

    // find how many different communities exist
    for (i=0; i<(*g_initial).n; i++) {
        int c = comm[i];

        if (std::find(g_vertexes.begin(), g_vertexes.end(), c) == g_vertexes.end()) {
            g_vertexes.push_back(c);
        }
    }

    // store the members of community for the prolong
    sort(g_vertexes.begin(), g_vertexes.end());
    g.n = g_vertexes.size();
    for (i=0; i<g.n; i++)
        (*g_initial).com_map.insert(std::pair<int,int>(g_vertexes[i], i));

    network new_net;

    std::vector<std::vector<int>> new_net_array(g.n, std::vector<int>(g.n, 0));
    std::vector<int> new_volumes(g.n, 0);

    #pragma omp parallel for shared(new_net_array) schedule(static, NUM_SPLIT)
    for (i=0; i<(*g_initial).n; i++) {
        int c_i = (*g_initial).com_map[comm[i]];
        vector<pair<node_id, weight>> neighbors = g_initial->net[i];
        for (auto it=neighbors.begin(); it<neighbors.end(); ++it) {
            int c_j = (*g_initial).com_map[comm[it->first]];
            #pragma omp atomic update
            new_net_array[c_i][c_j] += it->second;
        }

    }

    #pragma omp parallel for ordered shared(new_net) schedule(static, NUM_SPLIT)
    for (i=0; i<g.n; i++) {
        std::vector<pair<node_id, weight>> v;
        weight i_volume = new_net_array[i][i];
        for (j=0; j<g.n; j++) {
            if (new_net_array[i][j] > 0) {
                v.push_back(make_pair(j, new_net_array[i][j]));
                i_volume += new_net_array[i][j];
            }
        }
        new_volumes[i] = i_volume;
        #pragma omp ordered
        new_net.push_back(v);


    }

    g.volumes = new_volumes;
    g.net = new_net;
    g.weight_net = (*g_initial).weight_net; //the sum of all edges remains the same

    return g;

}


std::vector<int> PLM::prolong(GraphComm g_initial, std::vector<int> coarsened_comm) {

    std::vector<int> init_comm = g_initial.communities;
    std::vector<int> new_comm(g_initial.n, 0);

    #pragma omp parallel for schedule(static, NUM_SPLIT)
    for (int i=0; i<g_initial.n; i++) {
        int i_comm = init_comm[i];
        new_comm[i] = coarsened_comm[g_initial.com_map[i_comm]];
    }
    return new_comm;
}

std::pair <int, float> max_pair_arg (std::pair <int, float> r, std::pair <int, float> n) {
        return (n.second >= r.second) ? n : r;
}

int PLM_MPI::ReturnCommunity(int i, GraphComm g) {

    std::vector<pair<node_id, weight>> n_i = g.net[i];
        std::unordered_map<int,float> mod_map;
    std::vector<int> seen_communities(g.communities.size(), 0);
    seen_communities[g.communities[i]] = 1;

    #pragma omp parallel for shared(seen_communities) schedule(static, NUM_SPLIT)
    for (auto neighbor_it = n_i.begin(); neighbor_it < n_i.end(); ++neighbor_it) {
        // TODO: compute weights for all communities in a single iteration
        int c_n = g.communities[neighbor_it->first];
        if (seen_communities[c_n] == 0) {
                    mod_map[neighbor_it->first] = compute_modularity_difference(i, neighbor_it->first, g);
            #pragma omp atomic write
            seen_communities[c_n] = 1;
        }
        }

    std::pair<int, float> max_pair = std::make_pair(i, 0.0);

    #pragma omp declare reduction (maxpair : std::pair<int, float> : omp_out=max_pair_arg(omp_out,omp_in)) initializer(omp_priv = omp_orig)
    #pragma omp parallel for shared(mod_map) reduction(maxpair:max_pair) schedule(static, NUM_SPLIT)
    for (size_t b = 0; b < mod_map.bucket_count(); b++) {
                for (auto bi = mod_map.begin(b); bi != mod_map.end(b); bi++)
                        max_pair = max_pair_arg(max_pair, *bi);
        }

    modularity max_diff = max_pair.second;
    if (max_diff > 0.0) {
        return g.communities[max_pair.first];
    }
    else
        return g.communities[i];
}


void PLM_MPI::Local_move(GraphComm* graph) {
    int unstable = 1;
    while (unstable) {
        //print((*graph).communities);
        //cout << "----------------------------" << endl;
        unstable = 0;
        #pragma omp parallel for schedule(static, NUM_SPLIT)
        for (int i=0; i<(*graph).n; i++) {
            int i_comm = (*graph).communities[i];
            int z = ReturnCommunity(i, *graph);
            if (z != i_comm) { // TODO: change iff the total modularity is optimized!
                (*graph).communities[i] = z;
                unstable=1;
            }

        }
    }
}

std::vector<int> PLM_MPI::Recursive_comm_detect(GraphComm g, int n) {


    std::vector<int> c_singleton(g.n, 0);

    # pragma omp parallel for schedule(static, NUM_SPLIT)
    for (int i=0; i<g.n; i++) {
        c_singleton[i] = i;
    }
    g.communities = c_singleton;
    Local_move(&g);

    if (g.communities != c_singleton) {
        // GraphComm g_new = coarsen(&g);
        // std::vector<int> c_coarsened = Recursive_comm_detect(g_new);
        // g.communities = prolong(g, c_coarsened);
    }
    return g.communities;

}

std::pair<weight, weight*> get_weight_and_volumes(GraphComm* g, int world_rank, int world_size, MPI_Comm comm) {

    int *network = NULL;
    int n;

    if (world_rank == 0) {
        //TODO: examine potential use of https://www.open-mpi.org/doc/v1.4/man3/MPI_Type_vector.3.php
        // According to https://spcl.inf.ethz.ch/Teaching/2019-dphpc/recitations/mpi_part2.pdf pages 8-10
        // using the Vector will handle automatically the sending of it.
        // in order to avoid the creation of the adjacency matrix

        // create the adjacency matrix
        n=g->n;

        // calloc in order to zero-initialize
        network = (int *) calloc(n*n, sizeof(weight));
        cout << "network map" << endl;
        for (int u = 0; u < n; u++) {

            cout << u << ": ";
            for (auto it = (g->net[u]).begin() ; it != (g->net[u]).end(); ++it) {
                cout << "(" << it->first << ", " << it->second << ")";
                network[u * n + it->first] = it->second;
            }
            cout << endl;
        }

        cout << endl << "adjacency matrix:" << endl;
        for (int u = 0; u < n; u++) {
            cout << u << ": ";
            for (int i = 0; i < n; i++)
                cout << network[u * n + i] << " ";
            cout << endl;
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // scatter the nodes of the graph
    int rem = n % world_size; // nodes remaining after division among processes
    int displacement_index = 0;
    int chunk_size = (n / world_size) * n;
    int *sendcounts = (int *) malloc(sizeof(int) * world_size);
    int *displs = (int *) malloc(sizeof(int) * world_size);
    for (int i=0; i<world_size; ++i) {
        sendcounts[i] = chunk_size;
        if (rem > 0) {
            sendcounts[i] += n;
            rem--;
        }

        displs[i] = displacement_index;
        displacement_index += sendcounts[i];
    }

    // // print calculated send counts and displacements for each process
    // if (world_rank == 0 ) {
    //     for (int i = 0; i < size; i++) {
    //         printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
    //     }
    // }

    int elements_per_proc = sendcounts[world_rank];
    int nodes_per_proc = elements_per_proc / n;
    int *sub_network = (int *) malloc(sizeof(int) * elements_per_proc);
    MPI_Scatterv(network, sendcounts, displs, MPI_INT, sub_network, elements_per_proc, MPI_INT, 0, comm);

    int *partial_volumes = (int *) malloc(sizeof(int) * elements_per_proc);

    # pragma omp parallel for schedule(static, NUM_SPLIT)
    weight partial_weights_sum = 0;
    // in order to find the correct weights we need to know which node is exactly at the
    // original network
    int offset = displs[world_rank] / n;
    cout << "offset: " << offset << endl;
    for (int u = 0; u < nodes_per_proc; u++) {

        weight vol = 0;
        //add the weight of the neighbors to the total partial_weights_sum
        for (int i = 0; i < n; i++) {
            vol += sub_network[u*n + i];
        }
        partial_weights_sum += vol;
        // volume = weights of neighbors + 2 * weight connecting with itself
        partial_volumes[u] = vol + sub_network[u*n + (u + offset)];
    }

    // if (world_rank == 3)
    //     for (int u = 0; u < nodes_per_proc; u++) {
    //         cout << u << ": ";
    //         for (int i = 0; i < n; i++)
    //             cout << sub_network[u*n + i] << " ";
    //         cout << endl;
    //     }

    for (int u = 0; u < nodes_per_proc; u++) {

        printf("partial_volumes[%d]=%d\n", u+offset, partial_volumes[u]);
    }

    weight weight_net;
    MPI_Allreduce(&partial_weights_sum, &weight_net, 1, MPI_INT, MPI_SUM, comm);

    rem = n % world_size; // nodes remaining after division among processes
    displacement_index = 0;
    chunk_size = (n / world_size);
    int *recv_counts = (int *) malloc(sizeof(int) * world_size);
    int *recv_displs = (int *) malloc(sizeof(int) * world_size);
    for (int i = 0; i < world_size; ++i) {
        recv_counts[i] = chunk_size;
        if (rem > 0) {
            recv_counts[i]++;
            rem--;
        }

        recv_displs[i] = displacement_index;
        displacement_index += recv_counts[i];
    }

    // collect the volumes of graph nodes to each MPI node
    int *volumes = (int *) malloc(sizeof(int) * n);
    MPI_Allgatherv(partial_volumes, recv_counts[world_rank], MPI_INT, volumes, recv_counts, recv_displs, MPI_INT, comm);

    if (world_rank == 0) {
        // for (int i = 0; i < n; i++) {
        //     printf("volumes[%d] = %d\n", i, volumes[i]);
        // }
        g->weight_net = weight_net;
        //TODO: Maybe we do not need the map of the volumes
        // g->volumes = volumes;
    }

    //TODO: where we need the volumes
    return make_pair(weight_net, volumes);
}

void PLM_MPI::DetectCommunities(int world_rank, int world_size) {

    MPI_Comm comm = MPI_COMM_WORLD;
    cout << "world_rank: " << world_rank << endl;
    std::pair<weight, weight*> weight_vol_pair = get_weight_and_volumes(&graph, world_rank, world_size, comm);
    if (world_rank == 0) {
        cout << "weight: " << graph.weight_net << endl;
        // graph.communities = Recursive_comm_detect(graph, graph.n);
        //cout << "final communities: ";
        //print(graph.communities);
    }

}

// the same for LP and Louvain. TODO: defined once
void PLM::PrintCommunities(const std::string &file_name) {
    std::ofstream ofs;
    ofs.open(file_name, std::ios_base::out | std::ios_base::trunc);
    if (ofs.fail()) {
        std::cerr << "Error opening file " << file_name << std::endl;
        std::exit(2);
    }

    for (int i = 0; i < graph.n; i++) {
        ofs << graph.communities[i] << std::endl;
    }
}
