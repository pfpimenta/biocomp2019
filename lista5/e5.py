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

# dadas as sequencias (array de strings) (sequences)
# e o tamanho do padrao (motif_lenght)
# retorna
# a string consenso (consensus_string)
# e o score do consenso (score)
def find_motif(sequences, motif_lenght):
    num_sequences = len(sequences)

    # get start positions que maximizam o score
    # start_positions = start_positions(sequences_1) # TODO : a parte mais dificil eh aqui eu acho
    # DEBUG: pegar um vetor s "aleatorio"
    start_positions = np.zeros(num_sequences)
    # print(start_positions) # DEBUG

    consensus_string, score = consensus(sequences, start_positions, motif_lenght)
    
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
    consensus_string = get_consensus_string(pfm)
    
    # get score
    score = 10 # TODO

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
def get_consensus_string(pfm):
    _, motif_lenght = np.shape(pfm)
    consensus_string = ''  # init 
    for i in range(len(pfm)):
            letra = np.argmax(pfm[:, i]) # acha letra do consenso na posicao i
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

    return consensus_string

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
consensus_string, score = find_motif(sequences_1, motif_lenght)
# DEBUG:
print("\n\n...resultados:")
print("consensus: " + str(consensus_string))
print("score: " + str(score))


## exercicio a:

# Encontrar o motivo de tamanho 8 aceitando 2 mutações;
motif_lenght = 8
mutations_accepted = 2
#consensus_string, score = motif_finder()

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