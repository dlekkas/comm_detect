import sys
import matplotlib.pyplot as plt
import collections

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

			if threads_no in res_d:
				res_d[threads_no].append(time)
			else:
				res_d[threads_no] = [time]

	return res_d


# TODO: calculate whatever statistics we like
def calculate_statistics(results_dict):
	pass

def filter_dict(old_dict, limit):
	new_dict = {}
	for key, value in old_dict.items():
		if key < limit:
			new_dict[key] = value

	return new_dict

# plot a bar plot showing the average execution time
# for different number of cores used
def plot_scale_diagram(results_dict, algo):
	fig_name = algo + '-time.png'
	avg_dict = {}
	for thread_no, times in results_dict.items():
		if len(times) != 0:
			avg_dict[thread_no] = sum(times) / len(times)

	# keep only keys that are lower than 12
	avg_dict = filter_dict(avg_dict, 12)

	D = collections.OrderedDict(sorted((avg_dict.items())))
	#plt.bar(range(len(D)), D.values(), align='center')
	plt.plot(range(len(D)), list(D.values()), marker='o', linewidth=2)
	plt.xticks(range(len(D)), list(D.keys()))
	plt.xlabel("Number of cores")
	plt.ylabel("Execution of time (s)")
	plt.title("[Cores-Time Diagram]")
	plt.savefig(fig_name)


# plot a bar plot showing the average execution time
# for different number of cores used
def plot_speedup(results_dict, algo):
	fig_name = algo + '-speedup.png'
	avg_dict = {}
	for thread_no, times in results_dict.items():
		if len(times) != 0:
			avg_dict[thread_no] = sum(times) / len(times)

	speedup_dict={}
	for k in list(avg_dict.keys()):
		speedup_dict[k] = avg_dict[1]/avg_dict[k]
	
	speedup_dict = filter_dict(speedup_dict, 12)
	D = collections.OrderedDict(sorted((speedup_dict.items())))
	
	plt.plot(range(len(D)), list(D.values()), marker='o', linewidth=2)
	plt.xticks(range(len(D)), list(D.keys()))
	plt.xlabel("Number of cores")
	plt.ylabel("Speedup")
	plt.title("[Cores-Speedup Diagram]")
	plt.savefig(fig_name)
	


algos = ['plp']

for algo in algos:
	file_name = algo + '_results.txt'
	results_dict = parse_threads_and_time(file_name)
	max_num_threads = max(results_dict, key=int)
	plot_scale_diagram(results_dict, algo)
	#plot_speedup(results_dict, algo)
	calculate_statistics(results_dict)




