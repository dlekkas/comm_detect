import sys
import matplotlib.pyplot as plt


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


# plot a bar plot showing the average execution time
# for different number of cores used
def plot_scale_diagram(results_dict):
	avg_dict = {}
	for thread_no, times in results_dict.items():
		if len(times) != 0:
			avg_dict[thread_no] = sum(times) / len(times)

	D = collections.OrderedDict(sorted((avg_dict.items())))
	plt.bar(range(len(D)), D.values(), align='center')
	plt.xticks(range(len(D)), list(D.keys()))
	plt.xlabel("Number of cores")
	plt.ylabel("Execution of time (s)")
	plt.title("[Cores-Time Diagram]")
	plt.show()



algos = ['plp']

for algo in algos:
	file_name = algo + '_results.txt'
	results_dict = parse_threads_and_time(file_name)
	max_num_threads = max(results_dict, key=int)
	plot_scale_diagram(results_dict)
	calculate_statistics(results_dict)




