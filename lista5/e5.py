#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 5 de Biologia Computacional
# Pedro Foletto Pimenta, novembro de 2019

# objetivo: TODO
# consensus algorithm by Hertz, Stromo (1989)

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
def consensus(sequences, motif_lenght, mutations_accepted):
    # algorithm by Hertz, Stromo (1989)
    num_sequences = len(sequences)
    sequence_lenght = len(sequences[0])
    
    # TODO entender como fazer isso
    
    # The  algorithm  starts  by  forming  a  matrix  for  each  of  the
    # L-mers  in the first sequence 
    num_matrix = sequence_lenght - motif_lenght + 1
    matrix = np.zeros((4, motif_lenght)) # just to initiate the matrix_list
    matrix_list = np.array([matrix for i in range(num_matrix)])
    score_list = np.zeros(num_matrix) # na verdade uma lista dos valores de "information content" (Schneider et al., 1986)
    for i in range(num_matrix):
        cropped_seq = sequences[0][i:i+motif_lenght]
        print(cropped_seq)
        m = get_pfm_seq(cropped_seq)
        score = get_pfm_information_content(m)
        matrix_list[i] = m
        score_list[i] = score
    
    # DEBUG
    print("debuggy")
    print(matrix_list)
    print(np.shape(matrix_list))
    
    # Each  of these matrices is  then  combined  with  each  L-mer in  the
    # second  sequence  to form    new   matrices   containing   two   L-mers   
    # TODO
    for i in range(num_matrix):
        cropped_seq = sequences[1][i:i+motif_lenght]
        m = get_pfm_seq(cropped_seq)
        for i in range(num_matrix):
                new_matrix = m + matrix_list[i]
                
                #best_i_matrix
    
    # only  save  the  matrix  with  the  lowest  probability  of 
    # occurring  by chance,  i.e.  the 'best'  matrix
    # TODO
    
    consensus_string = "coco"
    score = 10
    
    return consensus_string, score

# dadas
# as sequencias (sequences),
# posicoes iniciais do motif (start_positions)
# e seu tamanho (motif_lenght)
# retorna
# a padrao (motif)
# e o score do consenso (score)
def get_motif(sequences, start_positions, motif_lenght):
    num_sequences = len(sequences)
    assert(num_sequences == len(start_positions))

    # corta as sequencias de acordo com start_positions e motif_lenght
    cropped_sequences = []
    for i in range(num_sequences):
        start_index = int(start_positions[i])
        end_index = start_index + motif_lenght
        seq = sequences[i]
        cropped_sequences.append(seq[start_index:end_index])
    
    # fazer tabela de frequencias
    pfm = get_pfm_seqs(cropped_sequences) # position frequency matrix
    
    motif, score = get_motif_string(pfm)
    return motif, score

# given an array of sequences (array of strings) already cropped
# returns the pfm (pattern frequency matrix)
def get_pfm_seqs(cropped_sequences):
    motif_lenght = len(cropped_sequences[0])
    
    pfm = np.zeros((4, motif_lenght)) # init table

    # fill table
    for seq in cropped_sequences:
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

# given a sequence already cropped (string)
# returns the pfm (pattern frequency matrix)
def get_pfm_seq(cropped_seq):
    motif_lenght = len(cropped_seq)
    
    pfm = np.zeros((4, motif_lenght)) # init table

    # fill table
    for i in range(motif_lenght):
        # TODO : refazer esse if feio de um jeito elegante kk
        if(cropped_seq[i] == 'a'):
            pfm[0, i] += 1
        elif(cropped_seq[i] == 'c'):
            pfm[1, i] += 1
        elif(cropped_seq[i] == 'g'):
            pfm[2, i] += 1
        elif(cropped_seq[i] == 't'):
            pfm[3, i] += 1
        else:
            print("ihhh deu merda")
    return pfm

def get_pfm_information_content(pfm):
    
    _, motif_lenght = np.shape(pfm)
    print(pfm)
    print(motif_lenght)
    
    for i in range(motif_lenght):
        pass
    return 0

# given a pattern frequency matrix
# returns the motif string ("consensus")
# and its score
def get_motif_string(pfm):
    _, motif_lenght = np.shape(pfm)
    
    # init 
    motif = ''
    score = 0
    
    for i in range(len(pfm)):
            letra = np.argmax(pfm[:, i]) # acha letra do consenso na posicao i
            score += pfm[letra, i]
            # TODO : mmudar esse if else feio pra algo bunitu
            if(letra == 0): # a
                motif += 'a'
            elif(letra == 1): # c
                motif += 'c'
            elif(letra == 2): # g
                motif += 'g'
            elif(letra == 3): # t
                motif += 't'
            else:
                print("ihhh deu merda") # nunca chega aqui

    return motif, score

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
# testando o algoritmo find_motif (consensus "basico")
# tu da sequences, starting_positions, e  
motif_lenght = 4
num_sequences = len(sequences_1)
start_positions = np.zeros(num_sequences) # DEBUG: pegar um vetor s "aleatorio"
motif, score = get_motif(sequences_1, start_positions, motif_lenght)
# DEBUG:
print("\n\n...resultados:")
print("motif: " + str(motif))
print("score: " + str(score))

## DEBUG test
#for i in range(num_sequences):
#    num_seqs = i+1
#    print(num_seqs)
#    seqs = sequences_1[:num_seqs]
#    start_positions = np.zeros(num_seqs)
#    motif, score = get_motif(seqs, start_positions, motif_lenght)


## exercicio a:

# Encontrar o motivo de tamanho 8 aceitando 2 mutações;
motif_lenght = 8
mutations_accepted = 2
consensus_string, score = consensus(sequences_1, motif_lenght, mutations_accepted)

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