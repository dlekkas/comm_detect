import sys
import matplotlib.pyplot as plt
import collections
import numpy as np
import math

def parse_threads_and_time(file_name):
	res_d = {}
	with open(file_name, 'r') as f:
		for line in f:
			line = line.replace(' ', '')
			fields = line.split(',')

			threads_no = int(fields[1].split(':')[1])
			time = float(fields[2].split(':')[1])

			if threads_no in res_d:
				res_d[threads_no].append(time)
			else:
				res_d[threads_no] = [time]

	return res_d



def calculate_median(results_dict):
	
    res_sort = []
    median = []          
    res_median = {}
    n_cycles=0
    
    ### find median for each experiment #####
    for thread_no, times in results_dict.items():
        res_sort = sorted(results_dict[thread_no])
        
        num = len(times)
        if (num % 2) == 0: #even
            n_median = int((num/3)-1)
            
        else:
            n_median = int(((num + 1)/2)-1)
           
        median = (res_sort[n_median])
        res_median[thread_no] = median
        print(res_median,res_sort)
        n_cycles=n_cycles+1
    
    
    
    return res_median


# plot a bar plot showing the median execution time
# for different number of cores used


def plot_median_scale_diagram(results_dict, res_median, algo):
    fig_name = algo + '-time-median.png'
          
    D = collections.OrderedDict(sorted((res_median.items())))
    log_input = [20, 21, 22, 23, 24, 25]
   
    plt.xticks(range(len(D)), log_input)
    plt.xlabel("log[input size]", fontsize=15)
    plt.ylabel("Execution of time (s)", fontsize=15)
    plt.title("Weak Scaling PLP", fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.bar(range(len(D)), list(D.values()), color='darkred', zorder=2)   
    y_hor = range(10, 50, 10)
    for y in y_hor:
        plt.axhline(y=y, linewidth=0.8, color='black', linestyle='--', zorder=1)
    plt.savefig(fig_name)



file_name = sys.argv[1]
algo = sys.argv[2]
results_dict = parse_threads_and_time(file_name)
print(results_dict)
#max_num_threads = max(results_dict, key=int)
res_median = calculate_median(results_dict)
print(res_median)
plot_median_scale_diagram(results_dict, res_median, algo)
