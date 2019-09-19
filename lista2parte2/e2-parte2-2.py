#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio 2 da lista 2 parte 2 de Biologia Computacional
# Pedro Foletto Pimenta, setembro de 2019
###

import numpy as np
from sequences import * # load data


# "regras" do alinhamento
GAP_SCORE = -4
MATCH_SCORE = 5
MISMATCH_SCORE = -3

# prints an alignment in an organized way
def print_alignment(alignment0, alignment1, alignment2):
	print("Alinhamento final:")
	print(alignment0)
	print(alignment1)
	print(alignment2)

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
				point_matrix[i][j] = j * GAP_SCORE
			elif(j==0): 	# first column
				point_matrix[i][j] = i * GAP_SCORE
			else: 		# rest of the table
				if(seq1[i-1] == seq2[j-1]): # match
					match_score = MATCH_SCORE
				else: # mismatch
					match_score = MISMATCH_SCORE
				expr1 = point_matrix[i-1][j-1] + match_score
				expr2 = point_matrix[i][j-1] + GAP_SCORE # -4 == +gap(seq1)
				expr3 = point_matrix[i-1][j] + GAP_SCORE # -4 == +gap(seq2) 
				point_matrix[i][j] = max(expr1, expr2, expr3)

	# find the best alignment
	alignment0 = "" # string q indica se ouve match, mismatch ou gap para cada posicao
	alignment1 = ""
	alignment2 = ""

	match_count = 0

	# start at the bottom right
	i, j = size1-1, size2-1
	alignment_score = point_matrix[i][j] # init total score

	# trace back to the upper left corner
	while (i > 0 and j > 0):
		up = point_matrix[i-1][j]
		left = point_matrix[i][j-1]
		diag = point_matrix[i-1][j-1]
		# find out best direction		
		if(diag >= up and diag >= left):
			i = i - 1
			j = j - 1
			alignment1 = seq1[i] + alignment1
			alignment2 = seq2[j] + alignment2
			if( seq1[i] ==  seq2[j]):
				alignment0 = "O" + alignment0 # match
			else:
				alignment0 = "X" + alignment0 # mismatch
		elif(up >= diag and up >= left):
			i = i - 1
			alignment0 = "-" + alignment0 # gap
			alignment1 = seq1[i] + alignment1
			alignment2 = "-" + alignment2
		else: # (left >= diag and left >= up)
			j = j - 1
			alignment0 = "-" + alignment0 # gap
			alignment1 = "-" + alignment1
			alignment2 = seq2[j] + alignment2
		# count score
		alignment_score = alignment_score + point_matrix[i][j]

	# get last remaining letters
	while(i > 0):
		i = i - 1
		alignment1 = seq1[i] + alignment1
		alignment2 = "-" + alignment2
	while(j > 0):
		j = j - 1
		alignment1 = "-" + alignment1
		alignment2 = seq2[j] + alignment2

	# identidade
	match_count = 0
	for i in range(len(alignment1)):
		if(alignment1[i] == alignment2[i]):
			match_count += 1
			
	# print point matrix
	print("...printing point matrix of dimensions " + str(point_matrix.shape))
	print(point_matrix)
	
	# print aligned sequences
	print_alignment(alignment0, alignment1, alignment2)

	# print "identidade do alinhamento" ?
	print("Identidade (?): " + str(match_count))
	
	# print score
	print("Alignment score: " + str(alignment_score))

	return alignment_score

	


#################################
# main


score = needleman_wunsch(seq_musmusculus, seq_musmusculus_mut)

