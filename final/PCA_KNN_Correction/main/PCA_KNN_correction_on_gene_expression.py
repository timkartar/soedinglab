import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Trans-eQTL joint analysis from all SNPs (TEJAS)')

    parser.add_argument('--expression',
                        type=str,
                        dest='expression_filename',
                        metavar='FILE',
                        help='input expression file')

    parser.add_argument('--output',
                        type=str,
                        dest='output_fileprefix',
                        metavar='FILE',
                        help='output file prefix')

    parser.add_argument('--k',
                        type=int,
                        dest='neighbours',
                        help='No of nearest neighbours')

    opts = parser.parse_args()
    return opts

def gene_distance(a, b):
    return np.linalg.norm(a - b)

opts=parse_args()
print(opts)

#read the gene expression data
genes = list()
expr = list()
sample_ID = list()

with open(opts.expression_filename,"r") as expression:
    header = expression.readline()
    sample_ID = header.strip().split("\t")[1:]
    for line in expression:
        linesplit = line.strip().split("\t")
        genes.append(linesplit[0])
        expr.append(np.array(linesplit[1:],dtype=float))
expr = np.array(expr) 

######### normalisation #########
expr = ((expr.T - np.mean(expr.T,axis = 0))/(np.std(expr.T,axis = 0))).T
print(np.mean(np.mean(expr,axis = 1)), np.mean(np.std(expr,axis = 1)))

#reduce the dimension to find the distance matrix between neighbours
pca = PCA(n_components=40)
pca.fit(np.transpose(expr))
expr_pca = pca.transform(np.transpose(expr))

nsample = expr_pca.shape[0]

# Create a distance matrix
distance_matrix = np.zeros((nsample, nsample))
for i in range(nsample):
    for j in range(i+1, nsample):
        d = gene_distance(expr_pca[i,:], expr_pca[j,:])
        #d = cosine_similarity(expr[:, i], expr[:, j])
        distance_matrix[i, j] = d
        distance_matrix[j, i] = d

#number of neighbours to select
neighbors_K = opts.neighbours
corrected_expression = np.zeros_like(expr)

for i in range(nsample):
    distance = distance_matrix[i, :]
    #neighbors = distance.argsort()[::-1][:neighbors_K]
    neighbors = np.argsort(distance)[: neighbors_K+1]
    corrected_expression[:, i] = expr[:, i] - np.mean(expr[:, neighbors[1:]], axis = 1)


distance_matrix = pd.DataFrame(distance_matrix)
distance_matrix.index = sample_ID
distance_matrix.columns = sample_ID
(distance_matrix).to_csv("sample_distance_matrix.csv",sep="\t")
corrected_expression = pd.DataFrame(corrected_expression)
corrected_expression.index = genes
corrected_expression.columns = sample_ID
corrected_expression.to_csv(opts.output_fileprefix,"\t")

