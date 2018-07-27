# Author: Raktim Mitra (email: timkartar7879@gmail.com)

import numpy as np

def genotype(GX, nsamples = 400, ntrans = 1000, nsnps = 40000, sigmabeta = 0.05, ngenes = 1000, fmin = 0.1, fmax = 0.9, cfrac = 0.1):

    maf = np.random.uniform(fmin, fmax, nsnps)
    nc  = np.random.gamma(ngenes*cfrac, scale = 1.0, size = ntrans).astype(int)
    dosage = np.zeros((nsnps, nsamples))
    betaTy = np.zeros((ntrans, nsamples))
    for i in range(nsnps - ntrans):
        dosage[i, :] = np.random.binomial(2, maf[i], nsamples)
    
    for i in range(ntrans):
        itrans =  nsnps - ntrans + i
        ncausal = min(ngenes, nc[i])
        betas  =  np.random.normal(0, sigmabeta, ncausal)
        choose =  np.random.choice(ngenes, ncausal, replace=False)
        beta0  = -np.log ((1 / maf[itrans]) - 1)
        betaTy[i, :] = np.dot(GX[choose, :].T, betas) + beta0

        prob = (1 / (1 + np.exp(-betaTy[i, :])))
        for j, f in enumerate(prob):
            bins = [(1 - f)**2, (1 - f**2)]
            x = np.random.rand()
            dosage[itrans, j] = np.digitize(x, bins) # if you don't understand, ask the author.
   
    return dosage, maf, betaTy
