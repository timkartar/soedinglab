import numpy as np
from scipy.stats import chi2

def chisquare(S, R, sigmabeta2, sigmax2, corr=False, maf=None, W=None):
    S2 = np.square(S)
    S2mod = S2 + sigmax2 / sigmabeta2
    mu = np.sum(S2 / S2mod)
    sigma = 2 * np.sum(np.square(S2 / S2mod))
    if corr:
        x4 = gt_fourth_moment(maf)
        corr_maf = (x4 - 3) * np.sum(np.diag(np.square(W)))
        #corr_maf = (np.mean(GT ** 4, axis = 1) - 3) * np.sum(np.diag(np.square(W)))
        sigma = np.mean(sigma + corr_maf)
    mscale = sigma / mu / 2.0
    df = mu / mscale
    Rscaled = R / mscale
    return Rscaled, df

def pvals(S, R, sigmabeta2, sigmax2, corr=False, maf=None, W=None):
    S2 = np.square(S)
    S2mod = S2 + sigmax2 / sigmabeta2
    mu = np.sum(S2 / S2mod)
    sigma = 2 * np.sum(np.square(S2 / S2mod))
    if corr:
        x4 = gt_fourth_moment(maf)
        corr_maf = (x4 - 3) * np.sum(np.diag(np.square(W)))
        sigma += corr_maf
    mscale = sigma / mu / 2.0
    df = mu / mscale
    Rscaled = R / mscale
    p = 1 - chi2.cdf(Rscaled, df)
    return p

def gt_fourth_moment(maf):
    f0 = np.square(1 - maf)
    f1 = 2.0 * maf * (1 - maf)
    f2 = np.square(maf)
    mu = 2 * maf
    sig = np.sqrt(f1)
    x4 = f0 *  (-mu/sig) ** 4 + f1 * ((1 - mu)/sig)**4 + f2 * ((2 - mu)/sig)**4
    return x4

def ultra_new_pvals(GT, R, W, sigmabeta2):
    mu2, mu4 = moment_data(GT)
    N = GT.shape[1]
    q11 = np.sum(W)
    q2  = np.sum(np.diag(W))
    muQ = (mu2/(N-1))*(N*q2 - q11)

    v31 = -1*(mu4/(N-1))
    v22 = (N*(mu2**2) - mu4)/(N-1)
    v211 = -1*(v31 + v22)/(N-2)
    v1111 = -3*v211/(N-3)

    q31 = np.dot(np.diag(W),np.sum(W,axis = 1))
    q4 = np.sum(np.square(np.diag(W)))
    q22 = np.sum(np.square(W))
    q211 = np.sum(np.square(np.sum(W,axis = 1)))

    sigma2 = v1111*(q11**2 - 2*q2*q11 - 4*q211 + 8*q31 + 2*q22 + q2**2 - 6*q4) + 2*v211*(q2*q11 + 2*q211 - 6*q31 - 2*q22 - q2**2 + 6*q4) + v22*(q2**2 + 2*q22 - 3*q4) + 4*v31*(q31 - q4) + mu4*q4
    
    sigma2 = sigma2 - muQ**2
    mscale = sigma2 / muQ / 2.0                                                                             
    df = muQ / mscale                                                                                      
    Rscaled = R / mscale                                                                                  
    p = 1 - chi2.cdf(Rscaled, df)                                                                         
    return p

def moment_data(GT):   #GT ixN
    GT2 = np.square(GT)
    GT4 = np.square(GT2)
    mu2 = np.mean(GT2, axis = 1)
    mu4 = np.mean(GT4, axis = 1)
    return mu2, mu4
    
