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
def k_means(k, points):
    # separa os pontos em k grupos/clusters
    pass # TODO

##
# tiradas de https://datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/
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
##

#######################################################
## main

# carregar dados do arquivo csv

# # versao sem pandas
# points = np.array([])
# with open('leukemia_big.csv') as csvfile:
#     csvreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    
#     # get first line ('ALL' or 'AML)
#     labels = next(csvreader)
#     labels = labels[0].split(',')

#     # get num and dim of points
#     num_points = len(labels)
#     dim = len([0 for row in csvreader])
#     #print("DEBUG num_points e dim: " +str(num_points)+" e "+str(dim)) #DEBUG

#     #points = np.empty((num_points, dim)) # init empty array
#     points = [ [] for i in range(dim)] # init empty array
#     print(points)

#     # fill points array
#     for row in csvreader:
#         coco = row[0].split(',')
#         coco = list(map(float, coco))
#         print(coco)
#         input()

# versao com pandas
df = pandas.read_csv('leukemia_big.csv', header=None)

# # dummy data
# points = [[0,1], [-2,-5], [-1,42]]
# labels = ['ALL','AML','ALL']

# get data from dataframe
labels = df.iloc[0]
df = df.drop(0) # remove labels from dataframe
points = df.values

# debug data print
#print(df)
#print(df.head())
#print(df.at([10, 10]))
print("DEBUG pointss: ")# + labels)
print(points)
print("DEBUG labels: ")# + labels)
print(labels)


k = 2
k_means(k, points)
k = 3
k_means(k, points)