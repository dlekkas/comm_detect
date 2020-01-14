import sys
import matplotlib.pyplot as plt
import collections
import numpy as np
import math


# parses a results file and returns a dictionary where each key
# is the number of threads and each value is a list of all the
# times that were recorded for this number of threads

def parse_threads_and_time(file_name):
    res_d = {}
    with open(file_name, 'r') as f:
         for line in f:
           line = line.replace(' ', '')
           fields = line.split(',')
                      
           threads_no = int(fields[1].split(':')[1])
           time = float(fields[2].split(':')[1])
           print(time)
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
            n_median = int((num/2)-1)
            median = (res_sort[n_median]+res_sort[n_median + 1])/2
            
        else:
            n_median = int(((num + 1)/2)-1)
            median = res_sort[n_median]
           
        #median = (res_sort[n_median])
        res_median[thread_no] = median
        #print(res_median,res_sort)
        n_cycles=n_cycles+1
    
    
    
    return res_median



def calculate_ci(results_dict,res_median):
	
    ### find CI lower and upper bounds ######
    
    z= 1.96 # for 95 % CI (ref. Le Boudec)
    #z= 2.58 # for 99 % CI (ref. Le Boudec)
    
    
    n_cycles=0
    for thread_no, times in results_dict.items():
         n_cycles = n_cycles + 1
    
    n = len(times) # number of repetitions of each experiment
    upper_lim = round(1+(n+z*math.sqrt(n))/2) - 1 # CI upper limit
    lower_lim = round((n-z*math.sqrt(n))/2) - 1
    print("upper_lim", upper_lim)
    print("lower_lim", lower_lim)
    res_ci = {}
    res_lolims = np.zeros(n_cycles)
    res_uplims = np.zeros(n_cycles)
    res_sort = []
    #lower_lim = 0 # if n > 5 comment these 2 lines
    #upper_lim = 3 # if n > 5 comment these 2 lines
    j=0
    
    for thread_no, times in results_dict.items():
		
        res_sort = sorted(results_dict[thread_no])
        uplims = res_sort[upper_lim]-res_median[thread_no] # upper CI limit difference form median
        print("res_median[thread_no]", res_median[thread_no])
        print("res_sort", res_sort)
        lolims = res_median[thread_no]-res_sort[lower_lim] # lower CI limit difference from median

        res_lolims[j] = lolims
        res_uplims[j] = uplims
        j=j+1
    print("upper_lim", res_uplims)
    print("lower_lim", res_lolims)   
    
    return res_lolims, res_uplims




def filter_dict(old_dict, limit):
	new_dict = {}
	for key, value in old_dict.items():
		if key < limit:
			new_dict[key] = value

	return new_dict



#### Calculate number of measurements 


##### Best that number of threads follows powers of two #####


# plot a bar plot showing the median execution time
# for different number of cores used


def plot_median_scale_diagram(results_dict, res_median, algo):
    fig_name = algo + '-time-mean.png'
    
    res_lolims, res_uplims = calculate_ci(results_dict,res_median)
    
    res_ideal = {}
    
    f=0
    res_ideal[1]=res_median[1]
    for k in list(res_median.keys()):
          
          res_ideal[k] = (1-f)*res_ideal[1]/k + f*res_ideal[1] # Amdahl's law
    
    
    D = collections.OrderedDict(sorted((res_median.items())))
    D1 = collections.OrderedDict(sorted((res_ideal.items())))
    
    plt.xticks(range(len(D)), list(D.keys()))
    plt.xlabel("Number of cores",fontsize=15)
    plt.ylabel("Execution Time (s)",fontsize=15)
    plt.title("[Cores-Time Diagram Median]",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    asymmetric_error = [res_lolims, res_uplims]  

    #print("upper_lim2", res_uplims)
    #print("lower_lim2", res_lolims)
   
    plt.plot(range(len(D)), list(D1.values()), color='g',label='Ideal Linear Bound', marker='o',linestyle= '--')  
    plt.errorbar(range(len(D)), list(D.values()),color='r',ecolor='black',label='Measurement Results',xerr=0.0,yerr=asymmetric_error,marker='s',linewidth=2)
    
    plt.legend(fontsize=12)

    plt.show()
    plt.savefig(fig_name)



	
def plot_speedup(results_dict, algo):
    fig_name = algo + '-speedup.png'
    
    res_median = calculate_median(results_dict)

    speedup_dict={}
    f=0.01 # How to calculate this?
    speedup_ideal = {}
    speedup_amd={}
    for k in list(res_median.keys()):
        
        speedup_dict[k] = res_median[1]/res_median[k]
        speedup_amd[k] = 1/(f+(1-f)/k) # Amdahl's law
        speedup_ideal[k] = k # Amdahl's law for f = 0. In the ideal case we assume 100 % paralelization
        
    
    
    #speedup_dict = filter_dict(speedup_dict, 12)
    D = collections.OrderedDict(sorted((speedup_dict.items())))
    D1 = collections.OrderedDict(sorted((speedup_ideal.items())))
    D2 = collections.OrderedDict(sorted((speedup_amd.items())))
    
    plt.plot(range(len(D)), list(D.values()),  color='r',marker='s', linewidth=2, label="Measurement Results")
    plt.plot(range(len(D)), list(D1.values()),  color='g',marker='o', linewidth=2, label="Ideal Linear Bound",linestyle= '--')
    #plt.plot(range(len(D)), list(D2.values()), marker='o', linewidth=2, label="results")
    
    plt.xticks(range(len(D)), list(D.keys()))
    plt.xlabel("Number of cores", fontsize=15)
    plt.ylabel("Speedup", fontsize=15)
    plt.title("[Cores-Speedup Diagram Median]",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)

    plt.show()
    plt.savefig(fig_name)	
	
	
algos = ['plp']
#algos = ['plm']
number_uk = '-2002'
for algo in algos:
	#file_name = algo + '_results_uk'+number_uk+'.txt'
	file_name = algo + '_results_uk'+'.txt'
	results_dict = parse_threads_and_time(file_name)
	
	
	
	max_num_threads = max(results_dict, key=int)
	#plot_scale_diagram(results_dict, algo)
	res_median = calculate_median(results_dict)
	plot_median_scale_diagram(results_dict, res_median, algo)
	#plot_speedup(results_dict, algo)
	
	
	
	#calculate_ci(results_dict)
	calculate_median(results_dict)
