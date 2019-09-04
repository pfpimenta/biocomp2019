#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio e da lista 1 de Biologia Computacional
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


# construct complementar DNA
complementaryDNA = ''

for i in range(len(fastaSeq)):
	if(fastaSeq[i] == 'G'):
		complementaryDNA += 'C'
	elif(fastaSeq[i] == 'C'):
		complementaryDNA += 'G'
	elif(fastaSeq[i] == 'T'):
		complementaryDNA += 'A'
	elif(fastaSeq[i] == 'A'):
		complementaryDNA += 'T'
	elif(fastaSeq[i] == 'N'):
		complementaryDNA += 'N'

print(complementaryDNA)

