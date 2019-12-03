#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 5 de Biologia Computacional
# Pedro Foletto Pimenta, novembro de 2019

# objetivo: TODO

#######################################################
## imports
import random
import numpy as np
import pandas
import time

#######################################################
## funcoes

def consensus(sequences):
    # TODO
    return 

#######################################################
## main


# data
sequences_1 =
['cccctgatagacgctatctggctatccacgtacataggtcctctgtgcgaatctatgcgtttccaaccat',
'agtttactggtgtacatttgatacgtacgtacaccggcaacctgaaacaaacgctcagaaccagaagtgc',
'aaaggagtccgtgcaccctctttcttcgtggctctggccaacgagggctgatgtataagacgaaaatttt',
'agcccctccgatgtaagtcatagctgtaactattacctgccacccctattacatcttacgtacgtataca',
'ctgggttatacaacgcgtcatggcggggtatgcgttttggtcgtcgtacgctcgatcgttaacgtaggtc']

sequences_2 =
['gtcacgcttctgcataccatcctgactactcgtggcgaatacggttcgtctcagaacattgacgagtaggacctccatgtacacgtgagttcgccagtagagggcagaactagaggcccgagctcgttacccagtatatgtactcggcacacactgggatataatactacacgggatactaatagtggcatatcacgccg',
'atccctctaacaagttgttttgacggaccgtatttccaaatgtgctcggcttcagaaacaacctttctgccctctactggcgacgtcacaacgacgacaacagaccatatggagtggaccctactcatgtaattgagaccgtcgcatgtagttgatttatgtaaacatatggctctagtttcaggcccctgtaaaggtaa',
'ttacataggttccttcacgtcactccttgtccgcgatatctcctcttacccttactaccaagcgtttcctgaaaggcaatgaaaagttgccatgcgctgtcgccagtagagggcagaataccaaggcgcttcagacaactgtcgctgttcgtgggtgggagggattgtatctataatataggatagttcgtatcgaaaaa',
'ttatcgaccgccactttctcgccagtagagggcagaaccacaaagtgactccccgagcaatggctgacctactagttatccggcatcacatcggcacatatacgggcgagaccgagccctctccgtaaccaccagtcccactacttcacaggcatatcctgtatcaatgaaatcacaaacgttcgcatgaagataatcgt',
'gacggcacattttaacggcccaggttggcagacaggaaatctacgatggtgctactgctttcccgagctctgccacgatgccacacagcacaattctgccctctactggcgatcacctcgaataaaaccgaatgcaagaccgagtaacagcggctggtaacatgcgggaggacgcgctttccgcaagtatattaataggt',
'tgcatcataggttagtaagagttataaatcttcgatccctaagtgtggtgcactactcggtcgacctcgcattgacacaacgcgagagtcgccagtagagggcagaagccggcacttttgacctcttctatagaaggtagaccgtgagatcgcgcccgaaggggcccgacggtctccaaggtggaacgtattaggtaatc',
'gccgtgtataggcctccgatcgtgcggttctgccctctactggcgaaaggggcatttgctattccaatcgcatagattaccaaataaaaaacgaaagaaggccgtccttgcaaagcttagtccttaaactgagatgcttggcgaccggccataagctccactcgcttgagcacatcaccaagaatcaaagtagcaaaccc',
'tgtccgctctgcagacgtccgggatctacgttggtgttctctctagtaacagtacggcagttctttttcggtccgaagcgaatcccacccgccaaggttacataagcattatctgaaggcaaccatacgaactctcattggctcgccagtagagggcagaaggacatcgtgtcatagcacatgcccacagaggagattcg',
'ccggtctcaatagccgaacgaggatcgactggtaggcgtgtcgggtgtgtgtggaccggctttggaagaaccacacttctctggccctcaattggccaaaagtcatcttaaggatcggttggccggcagcccctcgtgactacgataaccccaggttctgccctctactggcgaccttgacagagcacttacccactgta',
'aaaagagtagtgatgagttagaaagaatttaagggacatcctcttgatttggacggctatccccaggaaatcgtaggcgggggtgcacatggatatctttaggtattaattcccccattccctcttcgttctgccctctactggcgaatgttgctcgcaatactacagcctcctcaatacaggtagggtattttacatat']

print("rodando algoritmo de consensus...")
startTime = time.time() # medir o tempo de execucao a partir daqui

# get melhor combinaçao de 3572 genes
#best_solution, best_solution_score = consensus(points, labels)
print("\n\n...Melhor selecao de genes encontrada:\n" + str(best_solution))
print("...score: " + str(best_solution_score))
endTime = time.time()
totalTime = endTime - startTime
print("...tempo de execucao: %.3f segundos"%(totalTime))