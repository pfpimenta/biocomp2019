#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 1 de Biologia Computacional
# Pedro Foletto Pimenta, setembro de 2019
###


import sys


#print 'Number of arguments:', len(sys.argv), 'arguments.' #DEBUG
if(len(sys.argv) != 2):
	print("argumento especificando o arquivo FASTA faltando")
	exit()

filename = str(sys.argv[1]) #+ ".fasta"

# open file containing Cromossomo 7
fasta = open(filename)
fastaSeq = fasta.read()

# remove "\n"
fastaSeq = fastaSeq.replace("\n", "") 

# quando nao mutada esta subsequencia e unica no Cromossomo 7
mutated_subseq = "CAATTGAATAATTG"


subseq_occurences = fastaSeq.count(mutated_subseq)
#print("subseq_occurences " + str(subseq_occurences)) # DEBUG

possible_chars = ['A', 'C', 'G', 'T'] # characters that can appear in a DNA sequence
# test mutations to find original subsequence
for i in range(len(mutated_subseq)): # len(mutated_subseq) == 14
# for each position in the subsequence
	for c in possible_chars:
		# mutation:
		s = list(mutated_subseq)
		s[i] = c
		subseq = ''.join(s)
		subseq_occurences = fastaSeq.count(subseq)
		print(subseq + " ... " + str(subseq_occurences) + " occurences")
