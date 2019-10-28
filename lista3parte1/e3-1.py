#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 3 parte 2 de Biologia Computacional
# Pedro Foletto Pimenta, outubro de 2019
###


import sys

# ultrametric tree class
class U_Tree:

	# initialization
	def __init__(self):
        	self.left = None
        	self.right = None
        	self.left_dist = None
        	self.right_dist = None

	# verify if tree is ultrametric
	def is_ultrametric(self):
		# returns True if self is ultrametric
		# returns False if self is not ultrametric
		
		# sub trees must also be ultrametric
		if(isinstance(self.left, U_Tree)):
			if(not self.left.is_ultrametric()):
				return False
			left_dist_leaves = self.left_dist + self.left.get_leaves_dist()
		else:
			left_dist_leaves = self.left_dist

		# sub trees must also be ultrametric
		if(isinstance(self.right, U_Tree)):
			if(not self.right.is_ultrametric()):
				return False
			right_dist_leaves = self.right_dist + self.right.get_leaves_dist()
		else:
			right_dist_leaves = self.right_dist

		# left distance to the leaves == right distance to the leaves
		if(left_dist_leaves != right_dist_leaves):
			return False

		return True

	# distance from this point to the leaves
	def get_leaves_dist(self):

		# if it is ultrametric, it does not matter if we check left or right
		if(isinstance(self.left, U_Tree)):
			return self.left_dist + self.left.get_leaves_dist()
		else:
			return self.left_dist
	
	# printing function ( CALL THIS )
	def print_tree(self): 
		self.print_subtree(self.left, 1, self.left_dist)
		self.print_subtree(self.right, 1, self.right_dist)

	# auxiliary printing function
	def print_subtree(self, subtree, level, dist):
		if(isinstance(subtree, U_Tree)):
			self.print_subtree(subtree.left, level+1, subtree.left_dist)
			print(level*'   '+'-' + str(dist) + '-')
			self.print_subtree(subtree.right, level+1, subtree.right_dist)
		elif(isinstance(subtree, str)):
			print(level*'   ' + level*'--'+ '-'+str(dist)+ level*'--' +'\t'+ subtree)


# find smallest distance in the distance matrix
def get_smallest_dist(dist_matrix):
	pass # TODO
	# se nao tiver os 0s na matriz, nao vou precisar dessa funcao

# TODO terminar :
# returns dist_matrix with otu_a and otu_b fused into a new otu
def merge_matrix_otus(dist_matrix, otu_list, otu_a, otu_b):

	# remove otu_a and otu_b
	for otu in otu_list:
		if otu != otu_a:
			key = (otu, otu_a)
			dist_matrix.pop(key, None)
			key = (otu_a, otu)
			dist_matrix.pop(key, None)
		if otu != otu_b:
			key = (otu, otu_b)
			dist_matrix.pop(key, None)
			key = (otu_b, otu)
			dist_matrix.pop(key, None)

	# add new otu
	new_otu = otu_a + '-' + otu_b
	for otu in otu_list:
		if otu != otu_a and otu != otu_b:
			new_dist = 1  # TODO
			new_key = (new_otu, otu)
			dist_matrix[new_key] = new_dist
			new_key = (otu, new_otu)
			dist_matrix[new_key] = new_dist

	return dist_matrix


# agglomerative method for ultrametric trees (UPGMA)
def upgma(otu_list, dist_matrix):
	# otus : lista das OTUs no passo atual
	# dist_matrix : dicionario com as distancias entre as OTUs
	
	tree_clusters = {}

	# enquanto a arvore nao tiver completa
	while(len(otu_list)>1):

		# find smallest distance for clustering
		otu_a, otu_b = min(dist_matrix, key=dist_matrix.get)

		# branch length estimation
		branch_lenght = dist_matrix[(otu_a, otu_b)]/2

		# update distance matrix
		#print("DEBUG 1 : " + str(dist_matrix)) # DEBUG
		dist_matrix = merge_matrix_otus(dist_matrix, otu_list, otu_a, otu_b)
		#print("DEBUG 2 : " + str(dist_matrix)) # DEBUG
		
		# update OTU list
		new_otu = otu_a + '-' + otu_b
		otu_list.append(new_otu)
		otu_list.remove(otu_a)
		otu_list.remove(otu_b)

		# update tree
		tree = U_Tree()
		tree.right = otu_a
		tree.right_dist = branch_lenght
		tree.left = otu_b
		tree.left_dist = branch_lenght

		# DEBUG
		print("DEBUG otu list: " + str(otu_list))

	return tree
	
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

# DEBUG : test tree
test_tree = U_Tree()
test_tree.left = "folhaLOKA"
test_tree.left_dist = 5
test_tree.right = U_Tree()
test_tree.right_dist = 1
test_tree.right.right = "oloko"
test_tree.right.right_dist = 4
test_tree.right.left = "bixo"
test_tree.right.left_dist = 4

# DEBUG print
print("printing test tree:")
test_tree.print_tree()

# DEBUG distance to the leaves
leaves_dist = test_tree.get_leaves_dist()
print("leaves dist : " + str(leaves_dist))

if(test_tree.is_ultrametric()):
	print("tree is ultrametric")
else:
	print("tree is not ultrametric")


# inicializacao da lista de OTUs
otu_list = []
for key in dist_matrix:
		if(key[0] not in otu_list):
			otu_list.append(key[0])
		if(key[1] not in otu_list):
			otu_list.append(key[1])

# TODO: metodo UPGMA, e dai usar ele pra construir a arvore
tree = upgma(otu_list, dist_matrix)

# printar resultados
print("printing resulting tree:")
tree.print_tree()

