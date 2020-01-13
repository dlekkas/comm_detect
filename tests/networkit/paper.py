from networkit import *

#networkit.community.detectCommunities(G, algo=None, inspect=True)

f = open('jazz.graph', 'r')
lines = f.readlines()

info = lines[0].split()
n_nodes = int(info[0])
n_edges = int(info[1])

G = Graph(n_nodes, weighted=True)

cnt = 0
for line in lines[1:]:
    for edge in line.split():
        G.addEdge(cnt, int(edge)-1, w=1.0)

print(overview(G))

print("#################################")
print("#################################")
print("#################################")

# recurse (bool, optional) – use recursive coarsening
# refine (bool, optional) – Add a second move phase to refine the communities.
# gamma (double) – Multi-resolution modularity parameter: 1.0 -> standard mod
plm = community.PLM(G, gamma=1.0, refine=True, recurse=True)

communities = community.detectCommunities(G, algo=plm, inspect=True)

print("#################################")
print("#################################")
print("#################################")

# updateThreshold (int) – number of nodes that have to be changed in each iteration so that a new iteration starts.
plp = community.PLP(G, updateThreshold=3)

communities = community.detectCommunities(G, algo=plp, inspect=True)
