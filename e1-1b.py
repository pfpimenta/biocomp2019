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

# remove first line
fastaSeq = fastaSeq.splitlines(True)[1:]
fastaSeq = ''.join(fastaSeq)
# remove "\n"
fastaSeq = fastaSeq.replace("\n", "")


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


# quantidade de palindromos de tamanho 9 e 11
numPalindromos9 = len(ocorrenciasPalindromos9.keys())
numPalindromos11 = len(ocorrenciasPalindromos11.keys())

print("num de palindromos de tamanho 9: " + str(numPalindromos9))
print("num de palindromos de tamanho 11: " + str(numPalindromos11))

# numero de ocorencias de cada palindromo:
# informacao contida nos dicionarios ocorrenciasPalindromos9 e ocorrenciasPalindromos11




