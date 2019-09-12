#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio 1a da lista 2 de Biologia Computacional
# Pedro Foletto Pimenta, setembro de 2019
###


import sys
import numpy as np

# prints an alignment in an organized way
def print_alignment( seq1, seq2):
	# TODO
	pass

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

	print(point_matrix)

	# find the best alignment
	# TODO
			
	# print score, alignment
	# TODO



# main

# homo sapiens 
seq_homo = "VLSPADKTNVKAA"
# horse
seq_horse = "VLSAADKTNVKAA"

needleman_wunsch(seq_homo, seq_horse) # DEBUG
