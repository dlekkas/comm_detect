#!/bin/sh
set -x
algo=plm
input_file=netscience.graph 

job=$HOME/comm_detect/bin/$algo\_test
input_dir=$HOME/comm_detect/input

export LD_LIBRARY_PATH=/cluster/apps/gcc/6.3.0/lib64/:$LD_LIBRARY_PATH

for i in 1 2 4; do
	for j in $(seq 1 3); do
		echo "$i comma $j"
		export OMP_NUM_THREADS=$i
  		export OMP_NESTED=TRUE
		bsub -I -W 1:00 -n $((2 * $i)) "$job $input_dir/$input_file" >> tmp.res
		echo $(tail -n 1 tmp.res) >> $algo\_results_2.txt
		rm tmp.res
	done
done
