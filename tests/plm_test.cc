#include "../include/graph.h"
#include "../include/plm.h"
#include "../include/modularity.h"

#include <vector>
#include <sys/time.h>
#include <omp.h>
#include <chrono>
#include <math.h>

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <input-file>" << std::endl;
		std::exit(1);
	}
	std::string file_name = argv[1];

	int threads = omp_get_max_threads();

	#pragma omp parallel num_threads(threads)
	{
		/* Obtain thread number */
		int tid = omp_get_thread_num();

		if (tid == 0) {
			std::cout << "algo: PLM, threads: " << omp_get_num_threads() << ", ";

		}
	}

	/* initialize graph based on file and confirm correct parsing */
	GraphComm test_g;
	//test_g.PrintNetwork();
	
	auto start = std::chrono::system_clock::now();
        test_g.Net_init(file_name);

        auto end = std::chrono::system_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        cout << "Read file succesfully after " << total_time / 1000.0 << "s" << std::endl;


	/* detect communities of graph */
	PLM test_plm {&test_g };

        /* benchmark time of community detection */
	start = std::chrono::system_clock::now();
	test_plm.DetectCommunities();
	end = std::chrono::system_clock::now();

	total_time = std::chrono::duration_cast<
			std::chrono::milliseconds>(end - start).count();
	std::cout << "time (in sec): " << total_time / 1000.0 << endl;


	// Modularity!
	//
	
	/*int comms = *max_element(std::begin(test_plm.graph->communities), std::end(test_plm.graph->communities)) + 1;
	std::vector<int> w(comms, 0);
        std::vector<int> v(comms, 0);

	for (int i=0; i<test_plm.graph->n; i++) {
                int c_u = test_plm.graph->communities[i];
                v[c_u] += test_plm.graph->volumes[i];
                vector<pair<node_id, weight>> neighbors = test_plm.graph->net[i];
                for (auto it=neighbors.begin(); it<neighbors.end(); ++it) {
                        int c_j = test_plm.graph->communities[it->first];
                        if (c_u == c_j)
                                w[c_u] += it->second;
                }
        }

        cout << "comm size: " << comms << endl;
        cout << "network weight: " << test_plm.graph->weight_net << endl;
        test_plm.graph->weight_sq = pow(test_plm.graph->weight_net, 2);

        float mod=0.0;
        for (int i=0; i<(int) comms; i++) {
                float a = (float) w[i] / test_plm.graph->weight_net;
                float b =  (float) (v[i] * v[i]) / test_plm.graph->weight_sq;
                float mod_comm = a - b;
                mod += mod_comm;
        }

        //mod = compute_modularity(cs, &test_g);
        cout << "modularity: " << mod << endl;*/
        

	return 0;
}
