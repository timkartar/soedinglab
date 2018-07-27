#!/usr/bin/env python
import numpy as np
import ctypes
from scipy.stats import chi2
from statsmodels.sandbox.stats.multicomp import fdrcorrection_twostage as fdr
import os
import argparse
import time
import matplotlib.pyplot as plt
from iotools.dosageparser import DosageParser
from iotools.readVCF import ReadVCF 
from scipy.optimize import minimize
from scipy.special import erfinv
from inference.regulizer_optimizer import OptimizeRegularizer
from select_genes_per_snp import linear_reg
import hist_qq
import regnull
import regstat
from tqdm import tqdm

def parse_args():

    parser = argparse.ArgumentParser(description='Trans-Eqtls from Joint Association AnalysiS (TEJAAS)')

    parser.add_argument('--genotype',
                        type=str,
                        dest='genotype_filename',
                        metavar='FILE',
                        help='input genotype file')

    parser.add_argument('--transgeno',
                        type=str,
                        dest='trans_genofile',
                        metavar='FILE',
                        default = " ",
                        help='input potential trans-eQTLs genotype filename for optimization.')

    parser.add_argument('--optimize',
                        type=int,
                        dest='optimize_sigbeta',
                        default = 0,
                        help='1/0 option to optimize sigmabeta')

    parser.add_argument('--sigbeta',
                        type=float,
                        dest='sigbeta',
                        default = 0.006,
                        help='User given sigbeta')
    
    parser.add_argument('--sample',
                        type=str,
                        dest='sample_filename',
                        metavar='FILE',
                        help='input fam file')

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

    parser.add_argument('--start',
                        type=int,
                        dest='startsnp',
                        help='starting SNP index')

    parser.add_argument('--end',
                        type=int,
                        dest='endsnp',
                        help='ending SNP index')

    opts = parser.parse_args()
    return opts



def read_expression(filename):
    gene_names = list()
    expression = list()
    with open(filename, 'r') as genefile:
        header = genefile.readline()
        donorids = header.strip().split("\t")[1:]
        for line in genefile:
            linesplit = line.strip().split("\t")
            expression.append(np.array(linesplit[1:], dtype=float))
            gene_names.append(linesplit[0])
    expression = np.array(expression)
    return donorids, expression, gene_names


def norm_binom(genotype, f):
    genotype = (genotype - (2 * f)) / np.sqrt(2 * f * (1 - f))
    return genotype

def parse_geno (genotype_filename, sample_filename, startsnp, endsnp):                              #Read genotype here
    ds        = DosageParser(genotype_filename, sample_filename, startsnp, endsnp)
    dosage    = ds.dosage
    snpinfo   = ds.snpinfo
    donorids  = ds.sample_id
    nsnps     = ds.nsnps
    nsample   = ds.nsample
    return dosage, snpinfo, donorids, nsnps, nsample

                                                                                                    #quality check, matching
def quality_check ( sampleids, donorids):
    choose_ids     = [x for x in sampleids if x in donorids]
    dosage_indices = [i for i, x in enumerate(donorids)  if x in choose_ids]
    exprsn_indices = [i for i, x in enumerate(sampleids) if x in choose_ids]
    return dosage_indices, exprsn_indices

def read_distance(filename):
    samplelist = list()
    distances = list()
    with open(filename, 'r') as mfile:
        mfile.readline()
        for line in mfile:
            linestrip = line.strip().split("\t")
            samplelist.append(linestrip[0])
            distances.append(np.array(linestrip[1:],dtype=float))
    return np.array(distances), samplelist


opts = parse_args()                                                                                
out_fileprefix      = opts.output_fileprefix
genotype_filename   = opts.genotype_filename
sample_filename     = opts.sample_filename
expression_filename = opts.expression_filename
startsnp            = opts.startsnp
endsnp              = opts.endsnp
optim               = opts.optimize_sigbeta
user_sigbeta        = opts.sigbeta
transgeno           = opts.trans_genofile
                                                                                                    #Read gene expression here
sampleids, expression, gene_names = read_expression(expression_filename)


'''
       Read dosage maf > 0.1 and < 0.9, monoalleleic SNPs. \  
       Check the Dosageparser file before doing this as there \ 
       is dependencies on MAF.  
'''


                                                                                                    #Check if optimization is demanded.
if optim:
    print ("\nSigma beta optimization started. Reading from the provided trans-genotype file.")
    tic = time.time()
    try:
        trans_dosage, trans_snpinfo, trans_donorids, trans_nsnps, trans_nsample = parse_geno (transgeno, sample_filename, 0, 50000)
        dosage_indices, exprsn_indices = quality_check (sampleids , trans_donorids)
        conv_dosage = np.round(trans_dosage)                                                        #Best performance if dosages are rounded off
        geno        = conv_dosage[:, dosage_indices]
                                                   
                                                                                                    #Check the normalization method. This is done using data mean now
                                                                                                    #freq        = np.array([x.freq for x in trans_snpinfo])
        expr        = expression[:, exprsn_indices]
                                                                                                    #geno        = norm_binom(geno, freq)
        geno = ((geno.T - np.mean(geno.T,axis = 0)) / np.std(geno.T,axis = 0)).T
        
        optimize_sigbeta   = OptimizeRegularizer(geno, expr)
        optimize_sigbeta.update()
        optimize_sigbeta.niter
        sigbeta            = optimize_sigbeta.sigmareg
        toc = time.time()
        print ("Sigma beta optimization completed in :", toc - tic , "seconds\n")
        print ("Optimized sigma beta value is:" , sigbeta,"\n")
        del geno, conv_dosage, trans_dosage
    except OSError as e:
        print("\n",e, ". Trans-eQTLs file not provided for Optimization\n")
        raise SystemExit


else:
    print("\n=============================================================")
    print("\nsigma beta optimization not requested")
    if user_sigbeta == 0.006:
        print ("\nsigma beta not provided. Using default value of 0.006")
    else:
        print("\nsigma beta value provided as : " , user_sigbeta)
    print("=============================================================")
    sigbeta         = user_sigbeta


tic = time.time()


print ("\nReading Genotype")


dosage, snpinfo, donorids, nsnps, nsample = parse_geno ( genotype_filename, sample_filename, startsnp, endsnp)
dosage_indices, exprsn_indices = quality_check (sampleids , donorids)
conv_dosage = np.round(dosage)                                                                      #Best performance if dosages are rounded off
geno        = conv_dosage[:, dosage_indices]                                                    
freq        = np.array([x.freq for x in snpinfo])
#geno = ((geno.T - np.mean(geno.T,axis = 0)) / np.std(geno.T,axis = 0)).T

#Scaling and normalization
geno  = norm_binom(geno.T, freq).T  #commented out to later check it 
expr  = expression[:, exprsn_indices]
                                                                                                    # SXN
print ("Completed data loading and processing\n")

print(geno.shape, expr.shape)
##################################################################################    
'''choose_ids = [x for x in sampleids if x in donorids]
distance, sample = read_distance("PCA_KNN_Correction/main/sample_distance_matrix.csv")
distance_matrix = distance[exprsn_indices,:][:,exprsn_indices]

k=12
 

corrected_geno = np.zeros_like(geno)
for i,j in enumerate(choose_ids):
    distances = distance_matrix[i,:]
    neighbors =np.argsort(distances)[: k+1]
    corrected_geno[:, i] = geno[:, i] - np.mean(geno[:, neighbors[1:]], axis = 1)

geno = corrected_geno'''
#####################################################################################

#geno = ((geno.T - np.mean(geno.T,axis = 0)) / np.std(geno.T,axis = 0)).T

#####################################################################################

ori_shuffled_genom = geno.copy()

for i in range(ori_shuffled_genom.shape[0]):
    idx = np.random.permutation(np.arange(0,ori_shuffled_genom.shape[1]))
    np.random.shuffle(idx)
    ori_shuffled_genom[i,:] = ori_shuffled_genom[i,idx]
                                                                                                    #randoms = np.random.choice(geno.shape[0], 30000)
                                                                                                    #genos = geno[randoms,:]

genotype = ori_shuffled_genom
Rsvd = np.zeros(genotype.shape[0])
scaled_Rsvd = np.zeros(genotype.shape[0])
pval_R = np.zeros(genotype.shape[0])
pval_R_old = np.zeros(geno.shape[0])
#print("genotype",genotype.shape)
start = time.time()

for snp in tqdm(range(genotype.shape[0])):
#    print(genotype[snp:snp+1,:].shape)
    l_reg = linear_reg(genotype[snp:snp+1,:],expr,1000)
    selected_genes = l_reg.selected_variables
    #Calculate Qscores here
    sigma_geno = 1
    sigmax2    = sigma_geno**2
    sigmabeta2 = sigbeta**2
    
    Rscore, S, Wsvd = regstat.rscore_maf(geno[snp:snp+1,:], expr[selected_genes], sigmabeta2, sigmax2)
    #U, S, V_t = np.linalg.svd(np.transpose(expr[selected_genes]),full_matrices=False)                   #, Check this
    fact = sigma_geno**2 / sigbeta**2
    diagw = np.square(S) / (np.square(S) + fact)
    #Wsvd = np.dot(U, np.dot(np.diag(diagw), U.T))                                                       # Reverse because in the methods
    #Rsvd[snp]  = np.diag(np.dot(genotype[snp:snp+1,:], np.dot(Wsvd,  genotype[snp:snp+1,:].T)))         # geno should be SNP X Sample
    Rsvd[snp] = Rscore[0]
    #print(Rscore.shape)                                                                                                    #Inference parameters

    ############ after update on 22 May 2018 ##############
    pval_R_old[snp] = regnull.pvals(S, Rscore, sigmabeta2, sigmax2, corr=False)[0]
    pval_R[snp] = regnull.pvals(S, Rscore, sigmabeta2, sigmax2, corr=True, maf=np.array([freq[snp]]),W = Wsvd)
    #Rmean    = np.sum(diagw)
    #Rvar     = (2 * np.sum(np.square(diagw))) + (np.mean(genotype[snp:snp+1,:]**4, axis = 1) - 3) * np.sum(np.square(np.diag(Wsvd)))
    #Rvar = Rvar[0]
    #scaled_Rsvd[snp] = (Rsvd[snp] - Rmean)/np.sqrt(Rvar)
    #alpha    = Rvar / (2 * Rmean)
    #f_degree = Rmean / alpha

    #pval_R[snp] = 1 - chi2.cdf(Rsvd[snp]/alpha,f_degree)
end = time.time()
print((end - start)/60)
################## Plotting ###########################

fig = plt.figure(figsize = (10, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
nbins = 100
sizefactor = 1
hist_qq.plot(ax1, ax2, pval_R_old, nbins, 'uniform', loc = 0, scale = 1, size = sizefactor)
hist_qq.plot(ax3, ax4, pval_R, nbins, 'uniform', loc = 0, scale = 1, size = sizefactor)

plt.tight_layout()
plt.show()
