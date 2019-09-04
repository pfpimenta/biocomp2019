#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio c da lista 1 de Biologia Computacional
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

# remove first line
fastaSeq = fastaSeq.splitlines(True)[1:]
fastaSeq = ''.join(fastaSeq)
# remove "\n"
fastaSeq = fastaSeq.replace("\n", "")


ocorrencias_subseqs37 = {}

# iterate through chromossome 7
for i in range(len(fastaSeq)-37):
	# get subsequence:
	subseq = fastaSeq[i:i+37]
	# count occurence
	if(subseq in ocorrencias_subseqs37.keys()):
		ocorrencias_subseqs37[subseq] += 1
	else:
		ocorrencias_subseqs37[subseq] = 1

# numero de subsequencias de tamanho 37
print(len(ocorrencias_subseq37.keys()))
# numero de ocorrencias para cada subsequencia de tamanho 37
#print(ocorrencias_subseq37)
