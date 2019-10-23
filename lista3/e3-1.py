#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 3 de Biologia Computacional
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
	# talvez seja melhor refazer essa funcao
		if(isinstance(subtree, U_Tree)):
			#print(level*'-' + subtree.name)
			self.print_subtree(subtree.left, level+1, subtree.left_dist)
			print(level*'   '+'-' + str(dist) + '-')
			self.print_subtree(subtree.right, level+1, subtree.right_dist)
		elif(isinstance(subtree, str)):
			#print(level*'----' + subtree)
			print(level*'   ' + level*'--'+ '-'+str(dist)+ level*'--' +'\t'+ subtree) #versao melhorada, mas tem q conseguir a dist


# Agglomerative methods for ultrametric trees (UPGMA)
def upgma(dist_matrix, tree):

	# se a matriz soh tem um elemento, acabou
	if(len(dist_matrix)==1):
		return tree

	# TODO

	
# objetivos:
# - construcao de arvores filogeneticas
# - implementacao do metodo Agglomerative methods for ultrametric trees (UPGMA)

# matriz de distancias dos 5 primatas (Gorila, Orangotango, Humano, Chipanzé, Gibão):
# é um dict na real

dist_matrix = {}

dist_matrix[('gor', 'gor')] = 0.0
dist_matrix[('gor', 'ora')] = 0.189
dist_matrix[('gor', 'hum')] = 0.11
dist_matrix[('gor', 'chi')] = 0.113
dist_matrix[('gor', 'gib')] = 0.215

dist_matrix[('ora', 'gor')] = 0.189
dist_matrix[('ora', 'ora')] = 0.0
dist_matrix[('ora', 'hum')] = 0.179
dist_matrix[('ora', 'chi')] = 0.192
dist_matrix[('ora', 'gib')] = 0.211

dist_matrix[('hum', 'gor')] = 0.11
dist_matrix[('hum', 'ora')] = 0.179
dist_matrix[('hum', 'hum')] = 0.0
dist_matrix[('hum', 'chi')] = 0.09405
dist_matrix[('hum', 'gib')] = 0.205

dist_matrix[('chi', 'gor')] = 0.113
dist_matrix[('chi', 'ora')] = 0.192
dist_matrix[('chi', 'hum')] = 0.09405
dist_matrix[('chi', 'chi')] = 0.0
dist_matrix[('chi', 'gib')] = 0.214

dist_matrix[('gib', 'gor')] = 0.215
dist_matrix[('gib', 'ora')] = 0.211
dist_matrix[('gib', 'hum')] = 0.205
dist_matrix[('gib', 'chi')] = 0.214
dist_matrix[('gib', 'gib')] = 0.0

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
test_tree.print_tree()

# DEBUG distance to the leaves
leaves_dist = test_tree.get_leaves_dist()
print("leaves dist : " + str(leaves_dist))

if(test_tree.is_ultrametric()):
	print("tree is ultrametric")
else:
	print("tree is not ultrametric")


# criar arvore vazia

# TODO: metodo UPGMA, e dai usar ele pra construir a arvore

# find smallest distance

