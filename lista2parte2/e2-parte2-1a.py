#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio 1 da lista 2 parte 2 de Biologia Computacional
# Pedro Foletto Pimenta, setembro de 2019
###

import numpy as np
from sequences import * # load data


# "regras" do alinhamento
GAP_SCORE = -2
MATCH_SCORE = 1
MISMATCH_SCORE = -1

# prints an alignment in an organized way
def print_alignment(alignment1, alignment2):
	print("Final alignment:")
	print(alignment1)
	print(alignment2)


# retorna os indices i e j correspondentes a posicao com o maior valor da matriz
def index_highest_value(point_matrix):
	size_i, size_j = np.shape(point_matrix) # tamanho da matriz
	highest_value = 0 # maior valor da matriz achado ate agora
	highest_value_i = 0
	highest_value_j = 0
	
	# percorre a matriz
	for i in range(size_i):
		for j in range(size_j):
			# se for maior q o maior valor ate agora, atualiza
			if(point_matrix[i][j] > highest_value):
				highest_value = point_matrix[i][j]
				highest_value_i = i
				highest_value_j = j
	
	return highest_value_i, highest_value_j
				

# Smith–Waterman algorithm: global alignment of two sequences
def smith_waterman(seq1, seq2):

	# size of the grid
	size1 = len(seq1) + 1 # +1 for the extra column/row in the points table
	size2 = len(seq2) + 1

	# construct the grid
	point_matrix = np.zeros((size1, size2))
	for i in range(size1):
		for j in range(size2):
			if(i==0 or j == 0): 	# first row or first column
				point_matrix[i][j] = 0
			else: 			# rest of the table
				if(seq1[i-1] == seq2[j-1]): # match
					match_score = MATCH_SCORE
				else: 		# mismatch
					match_score = MISMATCH_SCORE
				expr1 = point_matrix[i-1][j-1] + match_score
				expr2 = point_matrix[i][j-1] + GAP_SCORE # -4 == +gap(seq1)
				expr3 = point_matrix[i-1][j] + GAP_SCORE # -4 == +gap(seq2) 
				point_matrix[i][j] = max(0, expr1, expr2, expr3)

	# find the best alignment
	alignment1 = ""
	alignment2 = ""

	match_count = 0

	# start at the bottom right
	i, j = index_highest_value(point_matrix)
	end_i, end_j = i, j
	alignment_score = point_matrix[i][j] # inicializar total score

	# traceback starting at the element with the highest score until 0 is encountered
	while (point_matrix[i][j] != 0):
		up = point_matrix[i-1][j]
		left = point_matrix[i][j-1]
		diag = point_matrix[i-1][j-1]
		# find out best direction		
		if(diag >= up and diag >= left):
			i = i - 1
			j = j - 1
			alignment1 = seq1[i] + alignment1
			alignment2 = seq2[j] + alignment2
		elif(up >= diag and up >= left):
			i = i - 1
			alignment1 = seq1[i] + alignment1
			alignment2 = "-" + alignment2
		else: # (left >= diag and left >= up)
			j = j - 1
			alignment1 = "-" + alignment1
			alignment2 = seq2[j] + alignment2
		# count score
		alignment_score = alignment_score + point_matrix[i][j]


	# guardar os indices de começo do alinhamento local
	begin_i, begin_j = i, j

	# identidade
	match_count = 0
	for i in range(len(alignment1)):
		if(alignment1[i] == alignment2[i]):
			match_count += 1
			
	# print point matrix
	print("...printing point matrix of dimensions " + str(point_matrix.shape))
	print(point_matrix)
	
	# print aligned sequences
	print("Alinhamento final de tamanho "+str(end_i-begin_i)+" :")
	print(alignment1 + " (indices de " + str(begin_i) + " ate " + str(end_i) + ")")
	print(alignment2 + " (indices de " + str(begin_j) + " ate " + str(end_j) + ")")

	# print "identidade do alinhamento" (?)
	print("Identidade: " + str(match_count))
	
	# print score
	print("Alignment score: " + str(alignment_score))

	return alignment_score

	


#################################
# main


score = smith_waterman(seq_haemoglobin, seq_cytoglobin)

		
