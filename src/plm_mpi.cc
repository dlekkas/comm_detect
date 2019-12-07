#include "../include/plm.h"
#include "../include/graph.h"
#include "../include/modularity.h"

#include <mpi.h>
#include <map>
#include <string.h>
#include <chrono>

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

    GraphComm g;

    if (world_rank == 0) {

        std::vector<int> g_vertexes;

        //TODO: maybe remove the use of g_vertexes
        // find how many different communities exist
        for (int i = 0; i < g_initial->n; i++) {
            int c = communities[i];

            // cout << "communities[" << i << "]=" << communities[i] << endl;
            if (std::find(g_vertexes.begin(), g_vertexes.end(), c) == g_vertexes.end()) {
                g_vertexes.push_back(c);
            }
        }

        // store the members of community for the prolong
        sort(g_vertexes.begin(), g_vertexes.end());

        new_n = g_vertexes.size();

        cout << "new_n: " << new_n << endl;
    }

    MPI_Bcast(&new_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    g.n = new_n;


    int rem_init = n % world_size; // nodes remaining after division among processes
    int displacement_index_init = 0;
    int chunk_size_init = (n / world_size) * n;
    int *send_counts_init = (int *) malloc(world_size * sizeof(int));
    int *send_displs_init = (int *) malloc(world_size * sizeof(int));
    for (int i = 0; i < world_size; i++) {
        send_counts_init[i] = chunk_size_init;
        if (rem_init > 0) {
            send_counts_init[i] += n;
            rem_init--;
        }

        send_displs_init[i] = displacement_index_init;
        displacement_index_init += send_counts_init[i];
    }

    int elements_per_proc_init = send_counts_init[world_rank];
    int nodes_per_proc_init = elements_per_proc_init / n;

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
    int node_offset      = send_displs[world_rank] / new_n;
    int node_offset_init = send_displs_init[world_rank] / n;
    cout << "node_offset: " << node_offset << endl;
    cout << "nodes_per_proc: " << nodes_per_proc << endl;
    cout << "node_offset_init: " << node_offset_init << endl;
    cout << "nodes_per_proc_init: " << nodes_per_proc_init << endl;
    //FIXME please nodes_per_proc & node_offset should be in the initial n not new_n
    for (int i = node_offset_init; i < nodes_per_proc_init + node_offset_init; i++) {
        int c_i = communities[i];
        int* neighbors = &(network_initial[i * n]);
        for (int j = 0; j < g_initial->n; j++) {
            int weight_i_j = neighbors[j];
            // if (weight_i_j != 0) {
            int c_j = communities[j];
            new_partial_network[c_i * new_n + c_j] += weight_i_j;
            // }
        }
    }

    int *new_network = (int *) malloc(new_n * new_n * sizeof(int));
    //TODO: check if MPI_Allreduce could go after the computation of volumes or completely
    //alleviate it
    MPI_Allreduce(new_partial_network, new_network, new_n * new_n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    cout << "new network: ....................." << endl;
    for (int u = 0; u < new_n; u++) {
        cout << endl << "u = " << u << ":";
        for (int i = 0; i < new_n; i++) {
            cout << " " << new_network[u*new_n+i];
        }
    }
    cout << endl;


    // BCast essential data (com_map)
    int *new_partial_volumes = (int *)calloc(new_n, sizeof(int));
    for (int u = node_offset; u < nodes_per_proc + node_offset; u++) {

        weight u_volume = new_network[u * new_n + u];
        for (int j = 0; j < new_n; j++) {
            if (new_network[u * new_n + j] > 0) {
                u_volume += new_network[u * new_n + j];
            }
        }
        new_partial_volumes[u - node_offset] = u_volume;
    }

    // for (int u = 0; u < nodes_per_proc; u++)
    //     cout << "new_partial_volumes[" << u + node_offset << "]=" << new_partial_volumes[u] << endl;

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

    if (world_rank == 0) {
    cout << "----------------------" << endl;
        for (int u = 0; u < new_n; u++)
            cout << "new_volumes[" << u << "]=" << new_volumes[u] << endl;
    cout << "^^^^^^^^^^^^^^^^^^^^^^" << endl;
    }

    g.volumes_array = new_volumes;
    g.adj_matrix = new_network;
    g.weight_net = g_initial->weight_net; //the sum of all edges remains the same

    return g;
}


int *PLM_MPI::prolong(GraphComm g_initial, int *coarsened_comm) {

    int n_initial  = g_initial.n;
    int *init_comm = g_initial.communities_array;
    int *new_comm = (int *) malloc(n_initial * sizeof(int));

    #pragma omp parallel for schedule(static, NUM_SPLIT)
    for (int i = 0; i < n_initial; i++) {
        new_comm[i] = coarsened_comm[init_comm[i]];
    }
    return new_comm;
}

std::pair <int, float> max_pair_arg (std::pair <int, float> r, std::pair <int, float> n) {
        return (n.second >= r.second) ? n : r;
}

// https://www.imsc.res.in/~kabru/parapp/lecture3.pdf
// https://stackoverflow.com/questions/21726783/scatter-a-c-vector-of-pairs-with-mpi

// void Build_mpi_weight_pair(pair<int, float> w, MPI_Datatype* mytype)
// {
//     int array_of_blocklengths[2]={1,1};
//     MPI_Datatype array_of_types[2]={MPI_INT, MPI_FLOAT};
//     MPI_Aint i_addr, j_addr;
//     MPI_Aint array_of_displacements[2]={0};
//     MPI_Get_address(&(w.first), &i_addr);
//     MPI_Get_address(&(w.second), &j_addr);
//     array_of_displacements[1]=j_addr-i_addr;
//     MPI_Type_create_struct(2,array_of_blocklengths,array_of_displacements,array_of_types,mytype);
//     MPI_Type_commit(mytype);
// }

int PLM_MPI::ReturnCommunity(int i, GraphComm g, int comm_size, int world_rank, int world_size) {

    int n = g.n;
    int *network = g.adj_matrix;
    int c = g.communities_array[i];

    // each MPI thread will hold the same list of the weights per community
	std::vector<int> weights_per_thread(comm_size, 0);
	std::vector<int> volumes(comm_size, 0);

    /* calculate volumes for all communties */
	//TODO: can we do it better? - compute it only once
	for (int j = 0; j < n; j++)
		volumes[g.communities_array[j]] += g.volumes_array[j];

    // cout << "neighbors of " << i << " :";
	// for (int j=0; j<g.n; j++)
	// 	cout << " " << network[i*n + j];
    // cout << endl;

    // for (int j=0; j<g.n; j++)
	// 	cout << "volumes[" << j << "] :" << volumes[j] << " ";

    int *n_i_arr = &(network[i * n]);
    std::vector<pair<node_id, weight>> n_i;
    for (int j = 0; j < n; j++)
    {
        if (n_i_arr[j] != 0)
            n_i.push_back(make_pair(j, n_i_arr[j]));
    }

    std::unordered_map<int,float> mod_map;

    for (auto neighbor_it = n_i.begin(); neighbor_it < n_i.end(); ++neighbor_it) {
        // TODO: compute weights for all communities in a single iteration
        int c_n = g.communities_array[neighbor_it->first];
		if ((int) neighbor_it->first != i) {
            weights_per_thread[c_n] += neighbor_it->second;
        }
    }

    /* Obtain thread number */
    std::vector<float> t_mod(comm_size, 0.0);
    weight weight_c = weights_per_thread[c];
    weight i_vol = g.volumes_array[i];
    weight volume_c = volumes[c] - i_vol;
    weight n_w = g.weight_net;

    /* find the id of communities that this thread will check */
    std::vector<int> comm_to_check;
    int c_number = world_rank;
    while (c_number < comm_size) {
        comm_to_check.push_back(c_number);
        c_number += world_size;
    }

    for (std::vector<int>::iterator it = comm_to_check.begin() ; it != comm_to_check.end(); ++it) {
        // compute difference in modularity for this community
        float a = ((1.0 * (weights_per_thread[*it] - weight_c)) / n_w);
        float b = (1.0 * (volume_c - volumes[*it]) * i_vol) / (2 * n_w * n_w);
		t_mod[*it] = a + b;
    }
    t_mod[c] = 0.0;

    // pair<int, float> result = std::make_pair(c, 0.0);
    int max_comm = c;
	float max_mod_change = 0.0;

    for (int k = 0; k < comm_size; k++)
        if (t_mod[k] > max_mod_change) {
            max_comm = k;
			max_mod_change = t_mod[k];
		}

    // int *results_comm = NULL;
    // float *results_mod_change = NULL;
	// std::vector<pair<int, float>> results(world_size, std::make_pair(c, 0.0)); /* each thread will write the best result it will find*/

    // if (world_rank == 0) {
    //     results_comm       = (int *)   malloc(world_size * sizeof(int));
    //     results_mod_change = (float *) malloc(world_size * sizeof(float));
    // }

    //TODO use MPI_MAXLOC
    // MPI_Gather(&max_comm, 1, MPI_INT, results_comm, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Gather(&max_mod_change, 1, MPI_FLOAT, results_mod_change, world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    int new_community;
    new_community = max_comm;

    // if (world_rank == 0) {
        // float max_mod = -0.1; // At least remaining in the same community will give greater modularity difference
        // int arg_max = -1;
        // for (int j = 0; j < world_size; j++) {
        //     if (results_mod_change[j] > max_mod) {
        //         max_mod = results_mod_change[j];
        //         arg_max = j;
        //     }
        // }

        // new_community = results_comm[arg_max];
	// }

    // MPI_Bcast(&new_community, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return new_community;

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

std::map<int, int> PLM_MPI::Map_communities(GraphComm g) {

	std::vector<int> comms;
	std::map<int,int> com_map;
	for (int i = 0; i < g.n; i++) {
		int c = g.communities_array[i];

		if (std::find(comms.begin(), comms.end(), c) == comms.end()) {
			comms.push_back(c);
		}
	}
	sort(comms.begin(), comms.end());

	for (int i = 0; i < (int) comms.size(); i++)
		com_map.insert(std::pair<int,int>(comms[i], i));

	return com_map;
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
    int k = 0;
    while (unstable) {
        ++k;
        //print((*graph).communities);
        partial_unstable = 0;
        for (int i = node_offset; i < nodes_per_proc + node_offset; i++) {
            int i_comm = communities[i];
            int z = ReturnCommunity(i, *graph, graph->n, world_rank, world_size);
            // cout << "----------------------" << endl;
            // cout << "i: " << i << endl;
            // cout << "z: " << z << endl;
            // cout << "i_comm: " << i_comm << endl;
            if (z != i_comm) { // TODO: change iff the total modularity is optimized!
                partial_communities[i - node_offset] = z;
                partial_unstable = 1;
            }
            MPI_Allgatherv(partial_communities, send_counts[world_rank], MPI_INT, communities, send_counts, displs,
                           MPI_INT, MPI_COMM_WORLD);
            // cout << "partial communities: ................." << endl;
            // for (int i = node_offset; i < nodes_per_proc + node_offset; i++)
            //     cout << partial_communities[i - node_offset] << " ";
            // cout << "......................................" << endl;
            // if (world_rank == 0) {
            //     cout << "communities: ................." << endl;
            //     for (int i = 0; i < n; i++)
            //     cout << communities[i] << " ";
            // cout << "......................................" << endl;
            // }
        }
        // does not work properly
        MPI_Allreduce(&partial_unstable, &unstable, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        // Since it is already assigned there is no need for graph.communities = communities_array
        // graph->communities_array = communities;
    }
    cout << "iters for Local_move to converge: " << k << endl;

    std::map<int, int> com_map = Map_communities(*graph);
	for (int i = 0; i < n; i++)
        graph->communities_array[i] = com_map[graph->communities_array[i]];

    // reduce

}

int *PLM_MPI::Recursive_comm_detect(GraphComm g, int world_rank, int world_size, int recursions) {

    cout << "Recursive_comm_detect iteration: " << recursions << endl;
    int n = g.n;
    // TODO only world_rank == 0 should do this
    int *c_singleton = (int *) calloc(n, sizeof(int));

    for (int i = 0; i < n; i++) {
        c_singleton[i] = i;
    }

    g.communities_array = c_singleton;
    auto start = std::chrono::system_clock::now();
    Local_move(&g, world_rank, world_size);
	auto end = std::chrono::system_clock::now();
	auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Local_move took: " << total_time << "ms" << endl;

    int *communities = g.communities_array;
    bool equal = true;
    //TODO: maybe parallelize it if n > some number
    for (int i = 0; i < n; i++) {
        if (communities[i] != i) { // c_singleton[i]) {
            equal = false;
            break;
        }
    }

    if (!equal) {
        cout << "coarsen" << endl;
        GraphComm g_new = coarsen(&g, world_rank, world_size);
        cout << "Recursive_comm_detect" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        int *c_coarsened = Recursive_comm_detect(g_new, world_rank, world_size, ++recursions);
        cout << "prolong" << endl;
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
        cout << "n: " << n << endl;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    graph.n = n;

    int *network = (int *) malloc (n * n * sizeof(int));
    graph.adj_matrix = network;
    if (world_rank == 0) {
        network = GetAdjacencyMatrix(&graph);

        cout << endl << "adjacency matrix:" << endl;
        // int n = graph.n;
        // for (int u = 0; u < n; u++) {
        //     cout << u << ": ";
        //     for (int i = 0; i < n; i++)
        //         cout << network[u * n + i] << " ";
        //     cout << endl;
        // }
    }

    MPI_Bcast(network, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    get_weight_and_volumes(&graph, n, world_rank, world_size, comm);
    if (world_rank == 0)
        cout << "weight: " << graph.weight_net << endl;

    graph.communities_array = Recursive_comm_detect(graph, world_rank, world_size, 0);
    if (world_rank == 0) {
        cout << "final communities:";
        for (int u = 0; u < n; u++)
            cout << " " << graph.communities_array[u];
        cout << endl;
    }

}

// the same for LP and Louvain. TODO: defined once
void PLM_MPI::PrintCommunities(const std::string &file_name) {
    std::ofstream ofs;
    ofs.open(file_name, std::ios_base::out | std::ios_base::trunc);
    if (ofs.fail()) {
        std::cerr << "Error opening file " << file_name << std::endl;
        std::exit(2);
    }

    for (int i = 0; i < graph.n; i++) {
        ofs << graph.communities_array[i] << std::endl;
    }
}
