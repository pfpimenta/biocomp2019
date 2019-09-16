#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio 2a da lista 2 de Biologia Computacional
# Pedro Foletto Pimenta, setembro de 2019
###


import sys
import numpy as np
from alpha_sequences import * # load data

# prints an alignment in an organized way
def print_alignment( seq1, seq2):
	# TODO
	pass

# returns the match score according to the BLOSUM62 matrix
def get_blosum62_score(seq1, seq2):
	return blosum62[seq1][seq2]

# needleman_wunsch algorithm: global alignment of two sequences
def needleman_wunsch(seq1, seq2):

	# size of the grid
	size1 = len(seq1) + 1 # +1 for the extra column/row in the points table
	size2 = len(seq2) + 1


	# construct the grid
	point_matrix = np.zeros((size1, size2))
	for i in range(size1):
		for j in range(size2):
			if(i==0): 	# first row
				point_matrix[i][j] = -j
			elif(j==0): 	# first column
				point_matrix[i][j] = -i
			else: 		# rest of the table
				if(seq1[i-1] == seq2[j-1]): # match
					match_score = 5
				else: # mismatch
					match_score = -3
				expr1 = point_matrix[i-1][j-1] + match_score
				expr2 = point_matrix[i][j-1] -4 # -4 == +gap(seq1) ????? 
				expr3 = point_matrix[i-1][j] -4 # -4 == +gap(seq2) ????? 
				point_matrix[i][j] = max(expr1, expr2, expr3)

	# DEBUG
	print("point_matrix:")
	print(point_matrix.shape)

	# find the best alignment
	# TODO
	# start at the bottom right
	i, j = size1-1, size2-1
	while (i != 0 or j != 0):
		up = point_matrix[i-1][j]
		down = point_matrix[i][j-1]
		diag = point_matrix[i-1][j-1]
		pass
			
	# print alignment table
	print(point_matrix)
	
	# print aligned sequences
	#print_alignment( ???

	# print "identidade do alinhamento" ?
	#print_alignment( ???
	
	# print score
	#print(score)

	



# main

# DEBUG
seq_a = "COCOCOCO"
seq_b = "ACACOCOO"


#needleman_wunsch(seq_human, seq_horse) # DEBUG
needleman_wunsch(seq_a, seq_b) # DEBUG

# qual das especies apresenta a maior semelhanca, em termos de sequencia, com a espécie humana (homo sapiens)
for seq in seq_list.remove(seq_human):
	#needleman_wunsch(seq, seq_human)
	pass
