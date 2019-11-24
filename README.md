# Parallel Implementations of Community Detection Algorithms

## Introduction
This repository refers to a group project for the "Design of Parallel and High Performance Computing"
course at ETH Zurich. The parallelized algorithms presented here are highly influenced by the paper
"Engineering Parallel Algorithms for Community Detection in Massive Networks".

Our goal is to implement the serial version of state-of-the-art algorithms used for community detection
and then implement their parallelized versions to extensively evaluate and compare their performance.
We are aiming to evaluate multiple frameworks used for parallel programming (i.e pthreads, MPI, openMP).

The algorithms implemented appear below:
* Parallel Label Propagation (PLP)


## Experiments
### Datasets:
We will evaluate our algorithms on datasets widely used for community detection, which can be found
[here](https://www.cc.gatech.edu/dimacs10/downloads.shtml). Those datasets are comprised of various
graphs with different properties and structure.

Each input file, on its first line, has two obligatory entries representing the number
of vertices and the number of edges in the graph. The third entry is 0 for unweighted graphs.
In case of an unweighted graph, the remaining n lines list the adjacent nodes of each node
(i.e line k lists the neighbours of node k-1).


## MPI

1. `module load /cluster/apps/modules/new/gcc/6.3.0`
1. `module load open_mpi/1.6.5`
1. `bsub -I -n 4 mpirun ./bin/plm_test_mpi ./input/kron_g500-logn20.graph`
