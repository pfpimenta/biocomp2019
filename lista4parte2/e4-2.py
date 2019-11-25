#!/usr/bin/python
# -*- coding: utf-8 -*-

# exercicio a da lista 4 parte 2 de Biologia Computacional
# Pedro Foletto Pimenta, novembro de 2019

# objetivo: encontrar um conjunto reduzido de genes
# que separe o conjunto de amostras e que aumente
# a predominância de rótulos ALL e AML em um dos grupos.

#######################################################
## imports
import random
import csv
import numpy as np
import pandas
import time

#######################################################
## funcoes

# retorna um conjunto reduzido de genes
# que separe o melhor possivel o conjunto de amostras
# de acordo com os rotulos (ALL e AML)
def alg_genetico(points, labels):
    
    # parametros fixos do alg genetico
    new_num_dim = 3572 # num de genes (features) escolhidos do total
    population_size = 50
    num_generations = 100
    # parametros 'variaveis' do alg genetico
    prob_mutacao = 0.25 # chance de ocorrer uma mutacao em um novo individuo

    num_points, num_dim = np.shape(points)

    # populacao aleatoria inicial
    population = np.array([random.sample(range(num_dim), new_num_dim) for i in range(population_size)])
    # print(np.shape(population)) # (50, 3572) # DEBUG

    for generation in range(num_generations):
        # avalia populacao
        scores = evaluate_population(population, points, labels)
        print("DEBUG scores: "+str(np.shape(scores))+" mean score: % .4f max score: % .4f min score: % .4f"%(np.mean(scores),np.max(scores), np.min(scores)))
        # gera nova populacao
        population = generate_new_population(population, scores, prob_mutacao)
        
    # sort population to get the best solution
    population = sort_population(population, scores)
    best_solution = population[0]
    classes, _ = k_means(2, points[:, best_solution])
    best_solution_score = get_clustering_score(classes, labels)

    # TODO: O algoritmo deve ser executado em sextuplicata sendo calculado a média e o desvio padrão da aptidão do indivı́duo encontrado na última geração.

    return best_solution#, best_solution_score

# gera nova populacao de solucoes com base na performance da ultima
def evaluate_population(population, points, labels):
    
    population_size = len(population)
    scores = np.zeros(population_size)
    
    for i in range(population_size):
        solution = population[i]
        # get only the features we want
        new_points = points[:, solution.astype(int)]
        # evaluate solution
        classes, _ = k_means(2, new_points)
        scores[i] = get_clustering_score(classes, labels)

    return scores

# gera nova populacao de solucoes com base na performance da ultima
def generate_new_population(population, scores, prob_mutacao):
    
    population, scores = sort_population(population, scores)

    population_size = np.shape(population)[0]
    #print("DEBUG population_size: "+str(population_size))

    new_population = np.zeros(np.shape(population))

    # primeiros 10% sao copiados dos 10% melhores da geracao passad
    new_population[0:int(population_size/10)] = population[0:int(population_size/10)]

    # proximos 30% sao cruzas de 
    # boas solucoes da populacao passada 
    # + chance de mutacao
    for i in range(int(population_size/10), int(population_size*4/10)):
        index_a = random.randint(0, int(population_size*5/10)) # um dos melhores 50%
        index_b = random.randint(0, int(population_size*2/10)) # um dos melhores 20%
        new_solution = crossover(population[index_a], population[index_b])
        # chance de mutação
        if(choseWithProb(prob_mutacao)):
            new_solution = mutacao(new_solution)
        new_population[i] = new_solution

    # proximos 30% sao cruzas de 
    # uma boa solucao da populacao passada
    # e uma solucao aleatoria da populacao passada
    # + chance de mutacao
    for i in range(int(population_size*4/10), int(population_size*7/10)):
        index_a = random.randint(0, int(population_size*3/10)) # um dos melhores 30%
        index_b = random.randint(0, population_size-1) # qualquer um
        new_solution = crossover(population[index_a], population[index_b])
        # chance de mutação
        if(choseWithProb(prob_mutacao)):
            new_solution = mutacao(new_solution)
        new_population[i] = new_solution

    # proximos 30% sao cruzas de duas solucoes aleatorias da populacao passada
    # + 100% de chance de mutacao
    for i in range(int(population_size*7/10), population_size):
        #print("DEBUG i "+str(i))
        index_a = random.randint(0, population_size-1) # qualquer um
        index_b = random.randint(0, population_size-1) # qualquer um
        #print("DEBUG index a e b : "+str(index_a)+" e "+str(index_b))
        new_solution = crossover(population[index_a], population[index_b])
        # 100% de chance de mutação
        new_solution = mutacao(new_solution)
        new_population[i] = new_solution

    return new_population


# retorna uma nova solucao gerada pela combinacao
# das solucoes a e b
def crossover(solucao_a, solucao_b):
    assert(len(solucao_a) == len(solucao_b))

    return solucao_a # TODO

# aplica uma mutacao em uma solucao
def mutacao(solucao):
    # TODO
    return solucao

# ordena as solucoes de uma populacao com base nos seus scores
def sort_population(population, scores):
    
    sorted_indexes = np.argsort(scores)
    sorted_population = population[sorted_indexes, :]
    sorted_scores = scores[sorted_indexes]

    return sorted_population, sorted_scores

# retorna o score de um agrupamento de pontos/vetores
# OBS: so funciona com clustering em dois grupos (k=2 : ALL e AML)
def get_clustering_score(classes, labels):

    assert(np.shape(classes)== np.shape(labels))
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

# retorna True com probabilidade oneProb
# e False com probabilidade (1 - oneProb)
def choseWithProb( oneProb ):
    zeroProb = 1 - oneProb
    result = np.random.choice([False,True], 1, p= [zeroProb, oneProb ])[0]
    return result

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

print("carregando dados...")

# carregar dados do arquivo csv
df = pandas.read_csv('leukemia_big.csv', header=None)

# get data from dataframe
labels = np.array(df.iloc[0].values) # get labels (first row)
df = df.drop(0) # remove labels from dataframe
df = df.T # transpose data
points = (df.values).astype(np.float) # convert strings to floats and put it in a numpy array
points = normalize_points(points) # normalize to [0,1] interval

k=2

print("rodando algoritmo genetico...")
startTime = time.time() # medir o tempo de execucao a partir daqui

# get melhor combinaçao de 3572 genes
best_solution = alg_genetico(points, labels)
#best_solution, best_solution_score = alg_genetico(points, labels) # TODO
# TODO: O algoritmo deve ser executado em sextuplicata sendo calculado a média e o desvio padrão da aptidão do indivı́duo encontrado na última geração.
# TODO: Faça um gráfico demonstrando a convergência ao longo das 100 gerações (considere a repetição que obteve melhor aptidão ao final). 
print("\n\n...Melhor selecao de genes encontrada: " + str(best_solution))
#print("...score: " + str(best_solution_score))
endTime = time.time()
totalTime = endTime - startTime
print("...tempo de execucao: " + str(totalTime))