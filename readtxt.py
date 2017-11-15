import numpy as np
import copy
from scipy import linalg as la
from collections import OrderedDict, defaultdict
from operator import itemgetter
import csv
import networkx as nx
import matplotlib.pyplot as plt

# read both predator-prey relation file and label file
# construct a original graph as well as a modified graph
def readFile():
	graph = defaultdict(list)
	label = defaultdict(list)
	with open("Predator-Prey Relations.txt") as f1:
		for line in f1:
			(predator, prey) = line.split(';')
			graph[int(prey)].append(int(predator))
	with open("Label.txt") as f2:
		for line in f2:
			(species, num) = line.split(';')
			label[int(num)] = species

	# get the original graph by adding unconnected vertices
	for key in label:
		graph[key]
	# deep copy the graph as original graph
	originalGraph = copy.deepcopy(graph)

	# get the modified graph by
	# adding edges from every node to decaying material
	for key in label:
		if key != 1:
			graph[key].append(1)
	return originalGraph, graph, label

# plot given graph with labels
def plot(graph, label):
	G=nx.DiGraph(graph)
	newLabel = dict()
	# truncate the label
	for num in label:
		data = label[num]
		if len(data) > 15:
			newLabel[num] = (data[:15] + '...')
		else:
			newLabel[num] = data
	# replace labels of numbers with labels of names
	H=nx.relabel_nodes(G,label)
	# set size of canvas
	plt.figure(figsize=(14,10))
	# pos = nx.random_layout(H) # random layout
	pos = nx.spring_layout(H,k=1,iterations=30) # Force-directed graph drawing
	nx.draw(H,pos,node_size=300,font_size=9,node_color='#A0CBE2',edge_color='#BB0000',width=0.5,edge_cmap=plt.cm.Blues,arrows=True,with_labels=True)
	plt.show()

# get the limiting distribution of given transition matrix
# return the rank of node by ranking the limiting distribution
def rank(G, label):
	# The transition matrix of a Markov Chain always has an eigenvalue 1
	# find the eigenvector corresponding to the eigenvalue 1
	e_vals, e_vecs = la.eig(G.T)
	i = 0
	while np.abs(e_vals[i] - 1.) > 1e-8:
		i = i + 1
	limiting_distribution = e_vecs[:, i]
	# normalize the eigenvector we found
	limiting_distribution = limiting_distribution / np.sum(limiting_distribution)
	d = dict()
	for i in label:
		d[label[i]] = limiting_distribution[i - 1]
	rank = OrderedDict(sorted(d.items(), key=itemgetter(1), reverse=True))
	return rank

# create a markov chain based on the given graph
# return a transition matrix of the markov chain
def get_markov_chain(graph, N):
	transition_matrix = np.zeros(shape = (N, N))
	for i in range(N):
		for j in range(N):
			edges = graph[i+1]
			if j+1 in edges:
				degree = len(edges)
				transition_matrix[i][j] = 1. / degree
	return transition_matrix

# construct a new matrix with damping factor
# by using pagerank algorithm
def get_google_matrix(T, N):
	d = 0.9 # damping factor
	G = d * T + (1-d)/N * np.ones((N, N))
	return G

# return a rank by degree of vertices
def rankByDegree(graph, label):
	d = dict()
	for i in label:
		d[label[i]] = len(graph[i])
	rank = OrderedDict(sorted(d.items(), key=itemgetter(1), reverse=True))
	return rank

# output the rank as a csv file
def output(dict1):
	with open('output.csv', 'wb') as output:
		writer = csv.writer(output)
		for key, value in dict1.iteritems():
		    writer.writerow([key, value])

if __name__ == '__main__':
	originalGraph, graph, label = readFile()
	N = len(graph)
	T = get_markov_chain(graph, N)
	G = get_google_matrix(T, N)
	rank1 = rankByDegree(originalGraph, label)
	rank2 = rankByDegree(graph, label)
	rank3 = rank(T, label)
	rank4 = rank(G, label)
	# plot(graph, label)
	# output(rank1)

