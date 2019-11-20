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

# retorna um numpy array de tamanho k
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
    classes = classes + k # sei la kkkk só pra evitar erro
    for i in range(num_points):
        p = points[i]
        min_dist = 99999
        for j in range(k):
            c = centroids[j]
            dist = np.linalg.norm(p-c)
            if(dist < min_dist):
                min_dist = dist
                classes[i] = j
    return classes

# retorna os centroides de cada classe 
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

# # dummy data
# points = [[0,1], [-2,-5], [-1,42], [10,1], [2,-3], [1,-10]]
# labels = ['ALL','AML','ALL', 'ALL','AML','AML']

# carregar dados do arquivo csv
df = pandas.read_csv('leukemia_big.csv', header=None)

# get data from dataframe
labels = df.iloc[0] # get labels (first row)
df = df.drop(0) # remove labels from dataframe
df = df.T # transpose data
points = (df.values).astype(np.float) # convert strings to floats and put it in a numpy array
# TODO : normalize points

# # debug data print
# print(df)
# print(df.head())
# print(df.at([10, 10]))
# print("DEBUG points: ")
# print(points)
# print("DEBUG labels: ")
# print(labels)
# print("DEBUG points[0]: ")
# print(points[0])

k = 2
classes, centroids = k_means(k, points)
print("\nResultado para K=2")
print(classes)
print(centroids)

k = 3
classes, centroids = k_means(k, points)
print("\nResultado para K=3")
print(classes)
print(centroids)