#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio b da lista 1 de Biologia Computacional
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

#numPalindromos9 = 0
#numPalindromos11 = 0
#listaPalindromos9 = []
#listaPalindromos11 = []
ocorrenciasPalindromos9 = {}
ocorrenciasPalindromos11 = {}

for i in range(len(fastaSeq)-9):
	subseq = fastaSeq[i:i+9]
	if(subseq == subseq[::-1]): # check if it is a palindrome
		if(subseq in ocorrenciasPalindromos9.keys()):
			ocorrenciasPalindromos9[subseq] += 1
		else:
			ocorrenciasPalindromos9[subseq] = 1

#print(ocorrenciasPalindromos9) # DEBUG


for i in range(len(fastaSeq)-11):
	subseq = fastaSeq[i:i+11]
	if(subseq == subseq[::-1]): # check if it is a palindrome
		if(subseq in ocorrenciasPalindromos11.keys()):
			ocorrenciasPalindromos11[subseq] += 1
		else:
			ocorrenciasPalindromos11[subseq] = 1

#print(ocorrenciasPalindromos11) # DEBUG

