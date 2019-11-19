#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 4 parte 1 de Biologia Computacional
# Pedro Foletto Pimenta, novembro de 2019

#######################################################
## imports
import pandas
import csv
import numpy as np

#######################################################
## funcoes
def k_means(k, points):
    # separa os pontos em k grupos/clusters
    pass # TODO


def cluster_points(X, mu):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
                    for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters
 
def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu
 
def has_converged(mu, oldmu):
    return  (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]))
 
def find_centers(X, K):
    # Initialize to K random centers
    oldmu = random.sample(X, K)
    mu = random.sample(X, K)
    while not has_converged(mu, oldmu):
        oldmu = mu
        # Assign all points in X to clusters
        clusters = cluster_points(X, mu)
        # Reevaluate centers
        mu = reevaluate_centers(oldmu, clusters)
    return(mu, clusters)


#######################################################
## main

# carregar dados do arquivo csv

#df = pandas.read_csv('leukemia_big.csv')#, index_col='Name')

#

points = np.array([])
with open('leukemia_big.csv') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    
    # get first line ('ALL' or 'AML)
    labels = next(csvreader)
    labels = labels[0].split(',')

    for row in csvreader:
        coco = row[0].split(',')
        coco = list(map(float, coco))
        print(coco)
        input()



# debug data print
#print(df)
#print(df.head())
#print(df.at([10, 10]))

points = df
points = [[0,1], [-2,-5], [-1,42]]

k = 2
k_means(k, points)
k = 3
k_means(k, points)