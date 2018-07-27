import numpy as np
import regnull
'''
  gene expression: Matrix of shape G x N
  genotype: Matrix of shape I x N
'''

## Qscore calculation
## Qscore calculation

def rscore_maf(GT, GX, sigmabeta2, sigmax2):
    nsnps = GT.shape[0]
    nsamples = GT.shape[1]
    ngenes = GX.shape[0]

    Yt = GX.T # shape N x G
    U, S, Vt = np.linalg.svd(Yt, full_matrices=False)
    S2 = np.square(S)
    S2mod = S2 + sigmax2 / sigmabeta2
    W = np.dot(U, np.dot(np.diag(S2 / S2mod), U.T))
    ## replace this with a vector array (pythonic version)
    Rscore = [None for i in range(nsnps)]
    for i in range(nsnps):
        Rscore[i] = np.sum(np.square(np.dot(U.T, GT[i,:])) * S2 / S2mod)
    Rscore = np.array(Rscore)

    return Rscore, S, W


def rscore(GT, GX, sigmabeta2, sigmax2, svd_corr = False):
    nsnps = GT.shape[0]
    nsamples = GT.shape[1]
    ngenes = GX.shape[0]
    Rscore = [None for i in range(nsnps)]
    pvals_ultra = [None for i in range(nsnps)]
    Yt = GX.T # shape N x G
    U, S, Vt = np.linalg.svd(Yt, full_matrices=False)
    #if svd_corr:
    #    GT = np.dot(U,GT.T).T
    #    GX = np.dot(U,GX.T).T
    S2 = np.square(S)
    for i in range(nsnps):
        S2mod = S2 + sigmax2[i] / sigmabeta2
        W = np.dot(U, np.dot(np.diag(S2 / S2mod), U.T)) #/ sigmax2[i]
        ## replace this with a vector array (pythonic version)
        Rscore[i] = np.sum(np.square(np.dot(U.T, GT[i,:])) * S2 / S2mod)
        pvals_ultra[i] = regnull.ultra_new_pvals(GT[i:i+1,:],np.array([Rscore[i]]),W,sigmabeta2)[0]
    Rscore = np.array(Rscore)
    pvals_ultra = np.array(pvals_ultra)

    return Rscore, S, pvals_ultra

def rscore_wfull(GT, GX, sigmabeta2, sigmax2):
    nsnps = GT.shape[0]
    nsamples = GT.shape[1]
    ngenes = GX.shape[0]

    Yt = GX.T # shape N x G
    U, S, Vt = np.linalg.svd(Yt, full_matrices=False)
    S2 = np.square(S)
    S2mod = S2 + sigmax2 / sigmabeta2

    Rscore = [None for i in range(nsnps)]
    Lk = np.diag(S2 / S2mod)
    W = np.dot(U, np.dot(Lk, U.T))
    Rscore  = np.diag(np.dot(GT, np.dot(W,  GT.T)))

    return Rscore, S, W

    # == DO NOT REMOVE ==
    #A = np.dot(Yt.T, Yt)
    #A[np.diag_indices(ngenes)] += sigmax2 / sigmabeta2
    #logdetA = np.linalg.slogdet(A)
