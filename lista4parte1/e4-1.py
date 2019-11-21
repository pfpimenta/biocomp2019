#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 4 parte 1 de Biologia Computacional
# Pedro Foletto Pimenta, novembro de 2019

#######################################################
## imports
import random
import csv
import numpy as np
import pandas

#######################################################
## funcoes

# retorna o score de um agrupamento de pontos/vetores
# OBS: so funciona com clustering em dois grupos (k=2 : ALL e AML)
def get_clustering_score(classes, labels):

    assert(np.shape(classes )== np.shape(labels))
    vecloko = np.ones(np.shape(classes)) # vetor de 1s para calculo do score

    # caso A -> ALL: 0, AML: 1
    matches_A = np.sum(vecloko[np.logical_and(labels=='ALL',classes==0)])
    matches_A = matches_A + np.sum(vecloko[np.logical_and(labels=='AML',classes==1)])
    
    # caso B -> ALL: 1, AML: 0
    matches_B = np.sum(vecloko[np.logical_and(labels=='ALL',classes==1)])
    matches_B = matches_B + np.sum(vecloko[np.logical_and(labels=='AML',classes==0)])

    # score = matches / num_points
    score = float(max(matches_A, matches_B)) / np.size(classes)

    return score

# retorna um numpy array com vetores normalizados de 0 a 1
def normalize_points(points):

    num_points, num_dim = np.shape(points)
    normalized_points = np.zeros((num_points, num_dim))

    # normalizar valores em cada dimensao
    for i in range(num_dim):
        # x = x - min(x)
        normalized_points[:, i] = points[:, i] -np.min(points[:, i])
        # x = x/max(x)
        normalized_points[:, i] = normalized_points[:, i] / np.max(normalized_points[:, i])

    return normalized_points

# retorna um numpy array contendo a posicao dos k centroides
def update_centroids(points, classes, k):

    num_points, num_dim = np.shape(points)
    centroids = np.zeros((k, num_dim))

    for i in range(k):
        if(points[classes==i].size == 0):
            # classe vazia -> reinicia cetroide aleatoriamente perto da media dos pontos
            centroids[i] = np.mean(points, axis=0) + 0.1 * np.random.rand(k, num_dim)
        else:
            # ajusta centroide para a media dos pontos contidos na sua classe
            centroids[i] = np.mean(points[classes==i], axis=0)
    return centroids

# retorna um numpy array onde cada posiçao i contem a classe atribuida ao ponto i
def classify_points(points, centroids, num_points, k):
    classes = np.zeros(num_points)
    for i in range(num_points):
        p = points[i]
        min_dist = 99999
        for j in range(k):
            c = centroids[j]
            dist = np.linalg.norm(p-c) # distancia entre ponto p e centroide c
            if(dist < min_dist): # atualiza classe do centroide mais perto do ponto p
                min_dist = dist
                classes[i] = j
    return classes

# retorna as classes de cada ponto e os centroides de cada classe 
# (dois numpy arrays)
def k_means(k, points):
    # points[num_points, num_dim] : pontos a classificar em k grupos
    # k : numero de grupos para classificar os pontos
    # classes[num_points] : classificacao dada a cada ponto de points
    # centroids[k, num_dim] : centroides de cada classe
    # old_centroids[k, num_dim] : centroides de cada classe na ultima iteracao (para criterio de parada)

    num_points, num_dim = np.shape(points)
    
    # gerar centroides iniciais aleatoriamente mas perto da media dos pontos
    centroids = [np.mean(points, axis=0) + 0.1 * np.random.rand(num_dim) for i in range(k)]
    old_centroids = np.zeros((k, num_dim))

    # classificacao inicial dos pontos
    classes = classify_points(points, centroids, num_points, k)
    
    # iterar ate convergencia
    while(not np.array_equal(centroids,old_centroids)): # enquanto os centroides nao convergirem
        old_centroids = centroids
        # atualiza centroides
        centroids = update_centroids(points, classes, k)
        # atualiza classificacao de cada ponto
        classes = classify_points(points, centroids, num_points, k)
    
    return classes, centroids

#######################################################
## main

# carregar dados do arquivo csv
df = pandas.read_csv('leukemia_big.csv', header=None)

# get data from dataframe
labels = np.array(df.iloc[0].values) # get labels (first row)
df = df.drop(0) # remove labels from dataframe
df = df.T # transpose data
points = (df.values).astype(np.float) # convert strings to floats and put it in a numpy array
points = normalize_points(points) # normalize to [0,1] interval

k = 2
classes, centroids = k_means(k, points)
score = get_clustering_score(classes, labels)
print("\nResultado para K=2")
print("...classes: "+str(classes))
print("...centroides: "+str(centroids))
print("...score: " + str(score))

k = 3
classes, centroids = k_means(k, points)
print("\nResultado para K=3")
print("...classes: "+str(classes))
print("...centroides: "+str(centroids))

# teste para ver score medio
k = 2
score_list = []
for i in range(100):
    classes, centroids = k_means(k, points)
    score = get_clustering_score(classes, labels)
    score_list.append(score)
print("\n\n\nscore medio: "+ str(np.mean(score_list)))