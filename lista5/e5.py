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

# dadas as sequencias (array de strings) (sequences),
# o tamanho do padrao (motif_lenght),
# e o numero maximo de mutacoes aceitas por padrao achado(mutations_accepted)
# retorna
# a string consenso (consensus_string)
# e o score do consenso (score)
def find_motif(sequences, motif_lenght):
    num_sequences = len(sequences)

    # TODO entender como fazer isso
    
    consensus_string = "coco"
    score = 10
    
    return consensus_string, score

# dadas
# as sequencias (sequences),
# posicoes iniciais do motif (start_positions)
# e seu tamanho (motif_lenght)
# retorna
# a string consenso (consensus_string)
# e o score do consenso (score)
def consensus(sequences, start_positions, motif_lenght):
    num_sequences = len(sequences)
    assert(num_sequences == len(start_positions))

    # corta as sequencias de acordo com start_positions e motif_lenght
    for i in range(num_sequences):
        start_index = int(start_positions[i])
        end_index = start_index + motif_lenght
        seq = sequences[i]
        sequences[i] = seq[start_index:end_index]
    
    # fazer tabela de frequencias
    pfm = get_pfm(sequences) # position frequency matrix
    print(pfm) # DEBUG

    consensus_string, score = get_consensus_string(pfm)
   

    return consensus_string, score

# given an array of sequences (array of strings)
# returns the pfm (pattern frequency matrix)
def get_pfm(sequences):
    num_sequences = len(sequences)
    motif_lenght = len(sequences[0])
    
    pfm = np.zeros((4, motif_lenght)) # init table

    # fill table
    for seq in sequences:
        for i in range(motif_lenght):
            # TODO : refazer esse if feio de um jeito elegante kk
            if(seq[i] == 'a'):
                pfm[0, i] += 1
            elif(seq[i] == 'c'):
                pfm[1, i] += 1
            elif(seq[i] == 'g'):
                pfm[2, i] += 1
            elif(seq[i] == 't'):
                pfm[3, i] += 1
            else:
                print("ihhh deu merda")
    return pfm

# given a pattern frequency matrix
# returns the consensus string
# and its score
def get_consensus_string(pfm):
    _, motif_lenght = np.shape(pfm)
    
    # init 
    consensus_string = ''
    score = 0
    
    for i in range(len(pfm)):
            letra = np.argmax(pfm[:, i]) # acha letra do consenso na posicao i
            score += pfm[letra, i]
            # TODO : mmudar esse if else feio pra algo bunitu
            if(letra == 0): # a
                consensus_string += 'a'
            elif(letra == 1): # c
                consensus_string += 'c'
            elif(letra == 2): # g
                consensus_string += 'g'
            elif(letra == 3): # t
                consensus_string += 't'
            else:
                print("ihhh deu merda") # nunca chega aqui

    return consensus_string, score

#######################################################
## main


# data
sequences_1 = ['cccctgatagacgctatctggctatccacgtacataggtcctctgtgcgaatctatgcgtttccaaccat',
'agtttactggtgtacatttgatacgtacgtacaccggcaacctgaaacaaacgctcagaaccagaagtgc',
'aaaggagtccgtgcaccctctttcttcgtggctctggccaacgagggctgatgtataagacgaaaatttt',
'agcccctccgatgtaagtcatagctgtaactattacctgccacccctattacatcttacgtacgtataca',
'ctgggttatacaacgcgtcatggcggggtatgcgttttggtcgtcgtacgctcgatcgttaacgtaggtc']

sequences_2 = ['gtcacgcttctgcataccatcctgactactcgtggcgaatacggttcgtctcagaacattgacgagtaggacctccatgtacacgtgagttcgccagtagagggcagaactagaggcccgagctcgttacccagtatatgtactcggcacacactgggatataatactacacgggatactaatagtggcatatcacgccg',
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


# DEBUG TEST
# testando o algoritmo consensus "basico"
# tu da sequences, starting_positions, e  
motif_lenght = 4
num_sequences = len(sequences_1)
start_positions = np.zeros(num_sequences) # DEBUG: pegar um vetor s "aleatorio"
consensus_string, score = consensus(sequences_1, start_positions, motif_lenght)
# DEBUG:
print("\n\n...resultados:")
print("consensus: " + str(consensus_string))
print("score: " + str(score))


## exercicio a:

# Encontrar o motivo de tamanho 8 aceitando 2 mutações;
motif_lenght = 8
mutations_accepted = 2
#consensus_string, score = motif_finder(sequences, motif_lenght, mutations_accepted)

# Encontrar o motivo de tamanho 5 aceitando 3 mutações;
motif_lenght = 5
mutations_accepted = 3
# consensus_string, score = consensus(sequences_1, )

# Encontrar o motivo de tamanho 3 aceitando 1 mutação;
motif_lenght = 3
mutations_accepted = 1
# consensus_string, score = consensus(sequences_1, )

## exercicio b:
# Encontrar o motivo de tamanho 3 aceitando 1 mutações;
motif_lenght = 3
mutations_accepted = 1
# Encontrar o motivo de tamanho 5 aceitando 2 mutações;
motif_lenght = 5
mutations_accepted = 2

endTime = time.time()
totalTime = endTime - startTime
print("...tempo de execucao: %.3f segundos"%(totalTime))