#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio d da lista 1 de Biologia Computacional
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

occurencesDict = {}

# iterate through chromossome 7
for i in range(len(fastaSeq)):
	# get char:
	c = fastaSeq[i]
	# count occurence
	if(c in occurencesDict.keys()):
		occurencesDict[c] += 1
	else:
		occurencesDict[c] = 1

print(occurencesDict)
