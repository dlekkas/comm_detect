import sys
import numpy as np
import matplotlib.pyplot as plt


def get_value():
	for i in range(3):
		line= fp.readline()
	tokens=line.split(':')
	a = float(tokens[1])	
	return a

max_num_threads=8
stride=1
x=range(1,max_num_threads+1,stride)
x_l = [str(i) for i in x]
y=[]
fp = open("results_time")

for i in range(max_num_threads):
	y.append(get_value())

y_s = [y[0]/i for i in y]
print(y_s)

fig, ax = plt.subplots()
index = np.arange(len(x))
bar_width = 0.08
opacity = 0.8

r1=plt.bar(x, y_s, align='center', color='darkred')

plt.xlabel('Number of Threads')
plt.ylabel('Speedup')

ax.set_xticks(x)
ax.set_xticklabels(x_l)

save="threads_speedup.png"
plt.savefig(save)
