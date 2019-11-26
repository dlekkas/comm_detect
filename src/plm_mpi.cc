#include "../include/plm.h"
#include "../include/graph.h"
#include "../include/modularity.h"

#include <mpi.h>
#include <map>
#include <string.h>

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

GraphComm PLM_MPI::coarsen(GraphComm *g_initial, int world_rank, int world_size) {

    int new_n;
    int n = g_initial->n;

    int *network_initial = g_initial->adj_matrix;
    int *communities     = g_initial->communities_array;

    // store the mapping of previous nodes to new ones (useful for prolong)
    int *communities_map_array = (int *) malloc(n * sizeof(int));

        GraphComm g;

    if (world_rank == 0) {

        std::vector<int> g_vertexes;

        //TODO: maybe remove the use of g_vertexes
        // find how many different communities exist
        for (int i = 0; i < g_initial->n; i++) {
            int c = communities[i];

            if (std::find(g_vertexes.begin(), g_vertexes.end(), c) == g_vertexes.end()) {
                g_vertexes.push_back(c);
            }
        }

        // store the members of community for the prolong
        sort(g_vertexes.begin(), g_vertexes.end());

        new_n = g_vertexes.size();

        // store the mapping of previous nodes to new ones (useful for prolong)
        // e.g. we have from the initial graph communities 2, 5 survived; so
        // now we translate these to 1, 2 respectively. So 2 of initial graph -> 1 of new
        // and 5 of initial graph -> 2 of new
        for (int i = 0; i < n; i++)
            communities_map_array[g_vertexes[i]] = i;
    }

    MPI_Bcast(&new_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    g.n = new_n;

    MPI_Bcast(communities_map_array, n, MPI_INT, 0, MPI_COMM_WORLD);
    g_initial->com_map_array = communities_map_array;

    //TODO: check if new_n, n are correct
    int rem = new_n % world_size; // nodes remaining after division among processes
    int displacement_index = 0;
    int chunk_size = (new_n / world_size) * new_n;
    int *send_counts = (int *) malloc(world_size * sizeof(int));
    int *send_displs = (int *) malloc(world_size * sizeof(int));
    for (int i = 0; i < world_size; i++) {
        send_counts[i] = chunk_size;
        if (rem > 0) {
            send_counts[i] += new_n;
            rem--;
        }

        send_displs[i] = displacement_index;
        displacement_index += send_counts[i];
    }

    int elements_per_proc = send_counts[world_rank];
    int nodes_per_proc = elements_per_proc / new_n;

    int *new_partial_network = (int *) calloc(new_n * new_n, sizeof(int));

    // we need to know which node is exactly at the original network
    int node_offset = send_displs[world_rank] / new_n;
    cout << "node_offset: " << node_offset << endl;

    for (int i = node_offset; i < nodes_per_proc + node_offset; i++) {
        int c_i = communities_map_array[communities[i]];
        int* neighbors = &(network_initial[i * n]);
        for (int j = 0; j < g_initial->n; j++) {
            int weight_i_j = neighbors[j];
            if (weight_i_j != 0) {
                int c_j = communities_map_array[communities[j]];
                new_partial_network[c_i * new_n + c_j] += weight_i_j;
            }
        }
    }

    int *new_network = (int *) malloc(new_n * new_n * sizeof(int));
    //TODO: check if MPI_Allreduce could go after the computation of volumes or completely
    //alleviate it
    MPI_Allreduce(new_partial_network, new_network, new_n * new_n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // BCast essential data (com_map)
    int *new_partial_volumes = (int *)calloc(new_n, sizeof(int));
    for (int u = node_offset; u < nodes_per_proc + node_offset; u++) {

        weight u_volume = new_network[u * new_n + u];
        for (int j = 0; j < new_n; j++) {
            if (new_network[u * new_n + j] > 0) {
                u_volume += new_network[u * new_n + j];
            }
        }
        new_partial_volumes[u] = u_volume;
    }

    // for (int u = 0; u < nodes_per_proc; u++)
    //     printf("partial_volumes[%d]=%d\n", u+offset, partial_volumes[u]);

    rem = new_n % world_size; // nodes remaining after division among processes
    displacement_index = 0;
    chunk_size = new_n / world_size;
    for (int i = 0; i < world_size; ++i) {
        send_counts[i] = chunk_size;
        if (rem > 0) {
            send_counts[i]++;
            rem--;
        }

        send_displs[i] = displacement_index;
        displacement_index += send_counts[i];
    }

    // collect the volumes of graph nodes to each MPI node
    int *new_volumes = (int *) malloc(new_n * sizeof(int));
    MPI_Allgatherv(new_partial_volumes, send_counts[world_rank], MPI_INT, new_volumes, send_counts, send_displs,
                   MPI_INT, MPI_COMM_WORLD);

    // if (world_rank == 0)
    //     for (int u = 0; u < new_n; u++)
    //         cout << "new_volumes[" << u << "]=" << new_volumes[u] << endl;

    g.volumes_array = new_volumes;
    g.adj_matrix = new_network;
    g.weight_net = g_initial->weight_net; //the sum of all edges remains the same

    return g;
}


int *PLM_MPI::prolong(GraphComm g_initial, int *coarsened_comm) {

    int n_initial  = g_initial.n;
    int *init_comm = g_initial.communities_array;
    int *new_comm = (int *) malloc(n_initial * sizeof(int));
    int *communities_map = g_initial.com_map_array;

    #pragma omp parallel for schedule(static, NUM_SPLIT)
    for (int i = 0; i < n_initial; i++) {
        new_comm[i] = coarsened_comm[communities_map[init_comm[i]]];
    }
    return new_comm;
}

std::pair <int, float> max_pair_arg (std::pair <int, float> r, std::pair <int, float> n) {
        return (n.second >= r.second) ? n : r;
}

int PLM_MPI::ReturnCommunity(int i, GraphComm g) {

    int n = g.n;
    int *network = g.adj_matrix;

    int *n_i_arr = &(network[i*n]);
    std::vector<pair<node_id, weight>> n_i;
    for (int j = 0; j < n; j++)
    {
        if (n_i_arr[j] != 0)
            n_i.push_back(make_pair(j, n_i_arr[j]));
    }

    std::unordered_map<int,float> mod_map;
    std::vector<int> seen_communities(g.n, 0);
    seen_communities[g.communities_array[i]] = 1;

    for (auto neighbor_it = n_i.begin(); neighbor_it < n_i.end(); ++neighbor_it) {
        // TODO: compute weights for all communities in a single iteration
        int c_n = g.communities_array[neighbor_it->first];
        if (seen_communities[c_n] == 0) {
            mod_map[neighbor_it->first] = compute_modularity_difference(i, neighbor_it->first, g);
            seen_communities[c_n] = 1;
        }
    }

    std::pair<int, float> max_pair = std::make_pair(i, 0.0);

    for (size_t b = 0; b < mod_map.bucket_count(); b++) {
        for (auto bi = mod_map.begin(b); bi != mod_map.end(b); bi++)
                max_pair = max_pair_arg(max_pair, *bi);
    }

    modularity max_diff = max_pair.second;
    if (max_diff > 0.0) {
        return g.communities_array[max_pair.first];
    }
    else
        return g.communities_array[i];
}

int *PLM_MPI::GetAdjacencyMatrix(GraphComm* g) {

    //this function is called only by one MPI node
    int n = g->n;

    int *network = g->adj_matrix;
    // zero-initialize
    memset(network, 0, n * n * sizeof(*network));

    for (int u = 0; u < n; u++)
        for (auto it = (g->net[u]).begin() ; it != (g->net[u]).end(); ++it)
            network[u * n + it->first] = it->second;

    return network;

}

void PLM_MPI::Local_move(GraphComm *graph, int world_rank, int world_size) {

    int n            = graph->n;
    int *communities = graph->communities_array;

    int rem = n % world_size; // nodes remaining after division among processes
    int displacement_index = 0;
    int chunk_size = n / world_size;
    int *send_counts = (int *) malloc(world_size * sizeof(int));
    int *displs      = (int *) malloc(world_size * sizeof(int));
    for (int i = 0; i < world_size; ++i) {
        send_counts[i] = chunk_size;
        if (rem > 0) {
            send_counts[i]++;
            rem--;
        }

        displs[i] = displacement_index;
        displacement_index += send_counts[i];
    }

    int node_offset = displs[world_rank];
    int nodes_per_proc = send_counts[world_rank];

    int *partial_communities = (int *) malloc(nodes_per_proc * sizeof(int));
    for (int i = node_offset; i < nodes_per_proc + node_offset; i++)
        partial_communities[i - node_offset] = communities[i];

    int partial_unstable;
    int unstable = 1;
    while (unstable) {
        //print((*graph).communities);
        //cout << "----------------------------" << endl;
        partial_unstable = 0;
        for (int i = node_offset; i < nodes_per_proc + node_offset; i++) {
            int i_comm = communities[i];
            int z = ReturnCommunity(i, *graph);
            if (z != i_comm) { // TODO: change iff the total modularity is optimized!
                partial_communities[i - node_offset] = z;
                partial_unstable = 1;
            }
        }
        MPI_Allreduce(&partial_unstable, &unstable, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        if (unstable)
            MPI_Allgatherv(partial_communities, send_counts[world_rank], MPI_INT, communities, send_counts, displs,
                           MPI_INT, MPI_COMM_WORLD);

    }
}

int *PLM_MPI::Recursive_comm_detect(GraphComm g, int world_rank, int world_size) {

    int n = g.n;
    // TODO only world_rank == 0 should do this
    int *c_singleton = (int *) calloc(n, sizeof(int));

    for (int i = 0; i < n; i++) {
        c_singleton[i] = i;
    }

    g.communities_array = c_singleton;
    Local_move(&g, world_rank, world_size);

    int *communities = g.communities_array;
    bool equal = true;
    //TODO: maybe parallelize it if n > some number
    for (int i = 0; i < n; i++) {
        if (communities[i] != c_singleton[i]) {
            equal = false;
            break;
        }
    }

    if (!equal) {
        GraphComm g_new = coarsen(&g, world_rank, world_size);
        int *c_coarsened = Recursive_comm_detect(g_new, world_rank, world_size);
        g.communities_array = prolong(g, c_coarsened);
    }

    return g.communities_array;
}

void get_weight_and_volumes(GraphComm* g, int n, int world_rank, int world_size, MPI_Comm comm) {

    int *network = g->adj_matrix;
    // scatter the nodes of the graph
    int rem = n % world_size; // nodes remaining after division among processes
    int displacement_index = 0;
    int chunk_size = (n / world_size) * n;
    int *send_counts = (int *) malloc(world_size * sizeof(int));
    int *displs      = (int *) malloc(world_size * sizeof(int));
    for (int i = 0; i < world_size; i++) {
        send_counts[i] = chunk_size;
        if (rem > 0) {
            send_counts[i] += n;
            rem--;
        }

        displs[i] = displacement_index;
        displacement_index += send_counts[i];
    }

    // // print calculated send counts and displacements for each process
    // if (world_rank == 0 ) {
    //     for (int i = 0; i < size; i++) {
    //         printf("send_counts[%d] = %d\tdispls[%d] = %d\n", i, send_counts[i], i, displs[i]);
    //     }
    // }

    int elements_per_proc = send_counts[world_rank];
    int nodes_per_proc = elements_per_proc / n;
    int *sub_network = (int *) malloc(elements_per_proc * sizeof(int));

    //TODO: examine potential use of https://www.open-mpi.org/doc/v1.4/man3/MPI_Type_vector.3.php
    // According to https://spcl.inf.ethz.ch/Teaching/2019-dphpc/recitations/mpi_part2.pdf pages 8-10
    // using the Vector will handle automatically the sending of it.
    // in order to avoid the creation of the adjacency matrix
    MPI_Scatterv(network, send_counts, displs, MPI_INT, sub_network, elements_per_proc, MPI_INT, 0, comm);

    int *partial_volumes = (int *) malloc(elements_per_proc * sizeof(int));

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

    // for (int u = 0; u < nodes_per_proc; u++) {

    //     printf("partial_volumes[%d]=%d\n", u+offset, partial_volumes[u]);
    // }

    weight weight_net;
    MPI_Allreduce(&partial_weights_sum, &weight_net, 1, MPI_INT, MPI_SUM, comm);

    rem = n % world_size; // nodes remaining after division among processes
    displacement_index = 0;
    chunk_size = n / world_size;
    int *recv_counts = (int *) malloc(world_size * sizeof(int));
    int *recv_displs = (int *) malloc(world_size * sizeof(int));
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
    int *volumes = (int *) malloc(n * sizeof(int));
    MPI_Allgatherv(partial_volumes, recv_counts[world_rank], MPI_INT, volumes, recv_counts, recv_displs, MPI_INT,
                   comm);

    g->weight_net = weight_net;
    //TODO: Maybe we do not need the map of the volumes
    //TODO: where we need the volumes
    g->volumes_array = volumes;
}

void PLM_MPI::DetectCommunities(int world_rank, int world_size) {

    MPI_Comm comm = MPI_COMM_WORLD;

    int n;
    if (world_rank == 0) {
        n = graph.n;
        cout << "world_size: " << world_size << endl;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    graph.n = n;

    int *network = (int *) malloc (n * n * sizeof(int));
    graph.adj_matrix = network;
    if (world_rank == 0) {
        network = GetAdjacencyMatrix(&graph);

        cout << endl << "adjacency matrix:" << endl;
        int n = graph.n;
        for (int u = 0; u < n; u++) {
            cout << u << ": ";
            for (int i = 0; i < n; i++)
                cout << network[u * n + i] << " ";
            cout << endl;
        }
    }

    MPI_Bcast(network, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    get_weight_and_volumes(&graph, n, world_rank, world_size, comm);
    if (world_rank == 0)
        cout << "weight: " << graph.weight_net << endl;

    graph.communities_array = Recursive_comm_detect(graph, world_rank, world_size);
    // if (world_rank == 0) {
    //     cout << "final communities: ";
    //     print(graph.communities);
    // }

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
