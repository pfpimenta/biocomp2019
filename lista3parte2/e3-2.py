#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 3 parte 2 de Biologia Computacional
# Pedro Foletto Pimenta, outubro de 2019
###


import sys

# tree class
class Tree:

	# initialization
	def __init__(self):
        	self.left = None
        	self.right = None
        	self.left_dist = None
        	self.right_dist = None

	# distance from this point to the leaves
	def get_leaves_dist(self):

		# if it is ultrametric, it does not matter if we check left or right
		if(isinstance(self.left, Tree)):
			return self.left_dist + self.left.get_leaves_dist()
		else:
			return self.left_dist
	
	# printing function ( CALL THIS )
	def print_tree(self): 
		self.print_subtree(self.right, 1, self.right_dist)
		self.print_subtree(self.left, 1, self.left_dist)


	# auxiliary printing function
	def print_subtree(self, subtree, level, dist):
		if(isinstance(subtree, Tree)):
			print(level*'   '  + '-')
			self.print_subtree(subtree.left, level+1, subtree.left_dist)
			print(level*'   ' + level*'-' + str(dist) + level*'-')
			self.print_subtree(subtree.right, level+1, subtree.right_dist)
			print(level*'   '  + '-')
		elif(isinstance(subtree, str)):
			print(2*level*'   ' + level*' - '+str(dist)+ level*' - ' +'\t'+ subtree)


# returns the q_matrix generated from dist_matrix
# according to the neighbour joining algorithm
def get_q_matrix(dist_matrix, otu_list):
	# dist_matrix : dicionario com as distancias entre as OTUs
	# otu_list : lista das "OTUs" no passo atual ("clusters")

	# q_matrix is a dict
	q_matrix = {}
	
	# create a q_matrix element for each dist_matrix element
	for otu_pair in dist_matrix.keys():
		
		# get otu pair
		otu_a, otu_b = otu_pair
		
		# sum of the distances for otu_a and otu_b
		sum_a, sum_b = 0, 0
		for otu in otu_list:
			if(otu != otu_a):
				sum_a = sum_a + dist_matrix[(otu_a, otu)]
			if(otu != otu_b):
				sum_b = sum_b + dist_matrix[(otu_b, otu)]

		# q_matrix formula
		q_matrix[otu_pair] = (len(otu_list) - 2) * dist_matrix[otu_pair] - sum_a - sum_b

	return q_matrix

# returns dist_matrix with otu_a and otu_b fused into a new otu
def update_dist_matrix(dist_matrix, otu_list, otu_a, otu_b):
	# dist_matrix : dicionario com as distancias entre as OTUs
	# otu_list : lista das "OTUs" no passo atual ("clusters")
	# otu_a, otu_b : OTUs a serem fusionadas em uma nova OTU

	# add new otu
	new_otu = otu_a + '-' + otu_b
	for otu in otu_list:
		if otu != otu_a and otu != otu_b:
			new_dist = (dist_matrix[(otu, otu_a)] + dist_matrix[(otu, otu_b)] - dist_matrix[(otu_a, otu_b)])/2
			new_key = (new_otu, otu)
			dist_matrix[new_key] = new_dist
			new_key = (otu, new_otu)
			dist_matrix[new_key] = new_dist

	# remove otu_a and otu_b
	for otu in otu_list:
		# remove otu_a distances
		key = (otu, otu_a)
		dist_matrix.pop(key, None)
		key = (otu_a, otu)
		dist_matrix.pop(key, None)
		# remove otu_b distances
		key = (otu, otu_b)
		dist_matrix.pop(key, None)
		key = (otu_b, otu)
		dist_matrix.pop(key, None)

	return dist_matrix

# Agglomerative methods for ultrametric trees (Neighbour Joining)
def neighbour_joining(dist_matrix):
	# dist_matrix : dicionario com as distancias entre as OTUs

	# inicializacao da lista de OTUs
	otu_list = []
	for key in dist_matrix:
		if(key[0] not in otu_list):
			otu_list.append(key[0])
		if(key[1] not in otu_list):
			otu_list.append(key[1])
	
	# initialize tree clusters
	tree_clusters = {}
	for otu in otu_list:
		tree_clusters[otu] = otu
	
	# enquanto a arvore nao tiver completa
	while(len(otu_list)>1):

		# calcular matriz Q
		q_matrix = get_q_matrix(dist_matrix, otu_list)

		# find smallest distance for clustering
		otu_a, otu_b = min(q_matrix, key=q_matrix.get)

		# branch lenght estimation
		sum_a, sum_b = 0, 0
		for otu in otu_list:
			if(otu != otu_a):
				sum_a = sum_a + dist_matrix[(otu_a, otu)]
			if(otu != otu_b):
				sum_b = sum_b + dist_matrix[(otu_b, otu)]
		branch_lenght_a = (dist_matrix[(otu_a, otu_b)])/2 + (1.0/(2*len(otu_list) - 2))*(sum_b - sum_a)
		branch_lenght_b = (dist_matrix[(otu_a, otu_b)]) - branch_lenght_a
		# versao com distancias arredondadas (mais paredida com a do wikipedia mas nao funciona para o nosso caso)
		#branch_lenght_a = int(round( (dist_matrix[(otu_a, otu_b)])/2 + (1.0/(2*len(otu_list) - 2))*(sum_b - sum_a) ))
		#branch_lenght_b = (dist_matrix[(otu_a, otu_b)]) - branch_lenght_a

		# update distance matrix
		dist_matrix = update_dist_matrix(dist_matrix, otu_list, otu_a, otu_b)
		
		# update OTU list
		new_otu = otu_a + '-' + otu_b
		otu_list.append(new_otu)
		otu_list.remove(otu_a)
		otu_list.remove(otu_b)

		# update tree : new tree node
		new_tree_node = Tree()
		new_tree_node.right = tree_clusters[otu_a]
		new_tree_node.right_dist = branch_lenght_a
		new_tree_node.left = tree_clusters[otu_b]
		new_tree_node.left_dist = branch_lenght_b
		
		# update tree clusters
		tree_clusters.pop(otu_a)
		tree_clusters.pop(otu_b)
		tree_clusters[new_otu] = new_tree_node 
	
		# DEBUG
		#print("DEBUG otu list: " + str(otu_list))

	return tree_clusters[otu_list[0]]	

# objetivos:
# - construcao de arvores filogeneticas
# - implementacao do metodo Agglomerative methods for ultrametric trees (UPGMA)


# matriz de distancias dos 5 primatas (Gorila, Orangotango, Humano, Chipanzé, Gibão):
# é um dict na real

dist_matrix = {}

#dist_matrix[('gor', 'gor')] = 0.0
dist_matrix[('gor', 'ora')] = 0.189
dist_matrix[('gor', 'hum')] = 0.11
dist_matrix[('gor', 'chi')] = 0.113
dist_matrix[('gor', 'gib')] = 0.215

dist_matrix[('ora', 'gor')] = 0.189
#dist_matrix[('ora', 'ora')] = 0.0
dist_matrix[('ora', 'hum')] = 0.179
dist_matrix[('ora', 'chi')] = 0.192
dist_matrix[('ora', 'gib')] = 0.211

dist_matrix[('hum', 'gor')] = 0.11
dist_matrix[('hum', 'ora')] = 0.179
#dist_matrix[('hum', 'hum')] = 0.0
dist_matrix[('hum', 'chi')] = 0.09405
dist_matrix[('hum', 'gib')] = 0.205

dist_matrix[('chi', 'gor')] = 0.113
dist_matrix[('chi', 'ora')] = 0.192
dist_matrix[('chi', 'hum')] = 0.09405
#dist_matrix[('chi', 'chi')] = 0.0
dist_matrix[('chi', 'gib')] = 0.214

dist_matrix[('gib', 'gor')] = 0.215
dist_matrix[('gib', 'ora')] = 0.211
dist_matrix[('gib', 'hum')] = 0.205
dist_matrix[('gib', 'chi')] = 0.214
#dist_matrix[('gib', 'gib')] = 0.0


# construcao de uma arvore ultrametrica filogenetica a partir da matriz de distancias usando UPGMA
tree = neighbour_joining(dist_matrix)

# printar resultados
print("printing resulting tree:")
tree.print_tree()
