#!/bin/sh

algo=plp
input_file=kron_g500-logn20.graph

job=$HOME/comm_detect/bin/$algo\_test
input_dir=$HOME/comm_detect/input

export LD_LIBRARY_PATH=/cluster/apps/gcc/6.3.0/lib64/:$LD_LIBRARY_PATH

for i in 1 2 4 6 8 10 12 14 16 18 20; do
	for j in $(seq 1 5); do
		echo "$i comma $j"
		export OMP_NUM_THREADS=$i
  		export OMP_NESTED=TRUE
		bsub -I -W 1:00 -n $i "$job $input_dir/$input_file" >> tmp.res
		echo $(tail -n 1 tmp.res) >> $algo\_results.txt
		rm tmp.res
	done
done
