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
import numpy as np
import pandas
import time

#######################################################
## parametros

# parametros fixos do alg genetico
POPULATION_SIZE = 50 # tem q ser 50
NUM_GENERATIONS = 100 # tem q ser 100
NEW_NUM_DIM = 23 # num de genes (features) escolhidos do total
# 5 -> 0.95 # 10 -> 0.97 # 15 -> 1.0 # 50 -> 1.0 # 150 # -> 1.0 # 500 -> 0.97 # 1500 -> 0.95 # 3572 -> 0.93

# parametros 'variaveis' do alg genetico
PROB_MUTACAO = 0.2 # chance de ocorrer uma mutacao em um novo individuo
NUM_MUTACOES = 10 # numero de genes mutados no evento de uma mutacao
NUM_SCORE_MEAN = 50 # numero de vezes a validar uma solucao pra fazer o score

#######################################################
## funcoes

# printa os parametros do algoritmo genetico
def print_params():
    print("Parametros do algoritmo genetico:")
    print("POPULATION_SIZE %i" % POPULATION_SIZE)
    print("NUM_GENERATIONS %i" % NUM_GENERATIONS)
    print("NEW_NUM_DIM %i" % NEW_NUM_DIM)
    print("PROB_MUTACAO %f" % PROB_MUTACAO)
    print("NUM_MUTACOES %i" % NUM_MUTACOES)
    print("NUM_SCORE_MEAN %i" % NUM_SCORE_MEAN)

# retorna um conjunto reduzido de genes
# que separe o melhor possivel o conjunto de amostras
# de acordo com os rotulos (ALL e AML)
def alg_genetico(points, labels):

    print_params()

    num_points, num_dim = np.shape(points)

    # populacao aleatoria inicial
    population = np.array([random.sample(range(num_dim), NEW_NUM_DIM) for i in range(POPULATION_SIZE)])
    melhores_scores = []
    scores_medias = []

    for generation in range(NUM_GENERATIONS):
        # avalia populacao
        scores = evaluate_population(population, points, labels)
        print("Generation "+str(generation)+" scores: "+str(np.shape(scores))+"  mean score: % .4f  max score: % .4f  min score: % .4f"%(np.mean(scores),np.max(scores), np.min(scores)))
        # gera nova populacao
        population = generate_new_population(population, scores, num_dim)

        # pra gerar o grafico q o dorn pediu:
        melhores_scores.append(np.max(scores))
        scores_medias.append(np.mean(scores))
        
    # sort population to get the best solution
    population, scores = sort_population(population, scores)
    best_solution = population[0].astype(int)
    classes, _ = k_means(2, points[:, best_solution])
    best_solution_score = scores[0]

    print("melhores_scores")
    print(melhores_scores)
    print("scores_medias")
    print(scores_medias)

    return best_solution, best_solution_score

# gera nova populacao de solucoes com base na performance da ultima
def evaluate_population(population, points, labels):
    
    #population_size = len(population)
    scores = np.zeros(POPULATION_SIZE)
    
    for i in range(POPULATION_SIZE):
        solution = population[i]
        # get only the features we want
        new_points = points[:, solution.astype(int)]
        # evaluate solution
        classes, _ = k_means(2, new_points)
        # avalia cada solucao em sextuplicata
        scores[i] = np.mean([get_clustering_score(classes, labels) for j in range(30)])

    return scores

# gera nova populacao de solucoes com base na performance da ultima
def generate_new_population(population, scores, num_dim):
    
    population, scores = sort_population(population, scores)

    population_size, new_num_dim = np.shape(population)
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
        if(choseWithProb(PROB_MUTACAO)):
            new_solution = mutacao(new_solution, num_dim)
        # ordena e remove repeticoes
        new_solution = ajusta_solucao(new_solution, num_dim)
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
        if(choseWithProb(PROB_MUTACAO)):
            new_solution = mutacao(new_solution, num_dim)
        # ordena e remove repeticoes
        new_solution = ajusta_solucao(new_solution, num_dim)
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
        new_solution = mutacao(new_solution, num_dim)
        # ordena e remove repeticoes
        new_solution = ajusta_solucao(new_solution, num_dim)
        new_population[i] = new_solution

    return new_population

# retorna a solucao sem repeticoes
def ajusta_solucao(solucao, num_dim):
    solucao = np.unique(solucao) # remove repeticoes
    # arruma tamanho do array se necessario 
    while(len(solucao) != NEW_NUM_DIM): # enquanto nao tiver no tamanho certo
        # completa array
        falta = NEW_NUM_DIM - len(solucao)
        extra_values = np.array(random.sample(range(num_dim), falta))
        solucao = np.append(solucao, extra_values)
        # ordena e remove repeticoes
        solucao = np.unique(solucao)
    return solucao

# retorna uma nova solucao gerada pela combinacao
# das solucoes a e b
def crossover(solucao_a, solucao_b):
    assert(len(solucao_a) == len(solucao_b))
    size = len(solucao_a)
    resultado = np.empty((size)).astype(int)

    #separacao = int(size/2) # crossover simples
    separacao = random.randint(1, size-2) # crossover 2.0 ihaaaa

    resultado[:separacao] = solucao_a[:separacao]
    resultado[separacao:] = solucao_a[separacao:]

    return resultado

# aplica uma mutacao em uma solucao
def mutacao(solucao, num_dim):
    # muda NUM_MUTACOES valores aleatoriamente
    for i in range(NUM_MUTACOES):
        random_index = random.randint(0, len(solucao)-1)
        solucao[random_index] = random.randint(0, num_dim-1)
    return solucao

# ordena as solucoes de uma populacao com base nos seus scores
def sort_population(population, scores):
    
    sorted_indexes = np.argsort(scores)
    sorted_population = np.array(population[sorted_indexes, :])
    sorted_scores = np.array(scores[sorted_indexes])

    return sorted_population[::-1], sorted_scores[::-1]

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
best_solution, best_solution_score = alg_genetico(points, labels)
print("\n\n...Melhor selecao de genes encontrada:\n" + str(best_solution))
print("...score: " + str(best_solution_score))
endTime = time.time()
totalTime = endTime - startTime
print("...tempo de execucao: %.3f segundos"%(totalTime))