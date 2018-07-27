import numpy as np
from scipy.optimize import minimize, fmin_l_bfgs_b
import regstat
from tqdm import tqdm
import qstats.crrstat as crrstat



class OptimizeRegularizer:

    def __init__(self, genotype, gene_expr, sigmax = 1, tol = 0.001, sigmabeta = 2):
        self._genotype = genotype
        self._expr = gene_expr
        self._sigmax = sigmax
        self._sigmareg = sigmabeta
        self._tolerance = tol
        self._niter = 0
        self._U, self._S, self._vt = self._svd(self._expr)
        
        #self._A =  np.dot(self._expr,self._expr.T) 
        #self._detA = np.linalg.det(self._A)
        #self._trA  = np.trace(self._A)
        #self._trA_2 = np.trace(self._A)**2
        #self._tr_A2 = np.trace(np.dot(self._A,self._A))

    def _svd(self, expr):
        U, S, Vt = np.linalg.svd(np.transpose(expr),full_matrices=False)
        return (U, S, Vt)

    @property
    def sigmareg(self):
        return self._sigmareg


    @property
    def niter(self):
        return self._niter
    
    def _logml(self, logsigmab, *args):
        sigmabeta = np.e ** logsigmab
        X, Y, S, U, sigmax = args
        nsnps      = 1
        nsamples   = X.shape[1]
        ngenes     = Y.shape[0]
        sigmabeta2 = sigmabeta * sigmabeta
        sigmax2    = sigmax**2
        lml = 0
        #eps = sigmabeta2/sigmax2
        #log_identity_term = - np.log(eps) * ngenes 
        #logdetA = log_identity_term + np.log(1 + self._detA + eps * self._trA + 0.5 * eps * eps * self._trA_2 - 0.5 * eps * eps *self._tr_A2)
        #A[np.diag_indices(ngenes)] += sigmax2 / sigmabeta2
        #logdetA = np.linalg.slogdet(A)
        logdetA = np.sum(np.log(np.square(S) + sigmax2 / sigmabeta2) ) + (ngenes - nsamples)*np.log(sigmax2 / sigmabeta2)
    
        Smod = np.diag(np.square(S) / (np.square(S) + sigmax2 / sigmabeta2))
        
        W = np.dot(U, np.dot(Smod, U.T))
            
        const_term = - 0.5 * ngenes * np.log(2 * np.pi * sigmabeta2) - 0.5 * logdetA            #[1]

        snp_term = 0.5 * np.dot(X, np.dot(W, X.T)) / sigmax2

        lml = lml + (const_term) + snp_term

        return -lml 
    
    def _grad_logml(self, logsigmab, *args):
        sigmabeta = np.e ** logsigmab          #x is logsigmabeta
        X, Y, S, U, sigmax = args
        nsnps = 1
        nsamples = X.shape[1]
        ngenes = Y.shape[0]
        sigmabeta2 = sigmabeta * sigmabeta
        sigmax2    = sigmax**2
    
        term1 = -ngenes

    
    
        der = 0
        Smod = np.square(S) * sigmabeta2 / sigmax2[i]
        term2 = (ngenes - nsamples) + np.sum(1 / (Smod + 1))
        term3 = 0
        for k in range(nsamples):
            uk = U[:, k]
            sk = S[k]
            smod = sk * sk * sigmabeta2 / sigmax2
            term3 += smod * np.square(np.dot(uk, X)) / sigmax2 / np.square(smod + 1)
        der += term1 + term2 + term3
        return der

    def update (self):
        sigmareg_old = np.e ** self._sigmareg
        #N = 20                                                                                         # expected number of true trans-eQTLs. Hard-coded
        iterate = True
        sigmax2 = self._sigmax ** 2
        optimised_betas = []
        for snp in tqdm(range(self._genotype.shape[0])):
            
            arguments = (self._genotype[snp:snp+1,:], self._expr, self._S, self._U, self._sigmax[snp])
        
            #print("optimizing...")
            #print(self._logml(self._sigmareg, *arguments))
            #res = fmin_l_bfgs_b(self._logml, [2], fprime = self._grad_logml, args = arguments)# jac=self._grad_logml)
            res = minimize(self._logml, [2], method="L-BFGS-B", args = arguments, options={'disp': True}) 
            #res = minimize(self._logml, [0], method="nelder-mead", args = arguments) 
            sigbeta_new = np.e**res.x[0]
            #print("optimised sigmabeta for snp ", snp, sigbeta_new)
            optimised_betas.append(sigbeta_new)
            #checksigma = self.check_convergence(sigbeta_new, sigmareg_old)

            #if checksigma:
            #    iterate = False
            #    sigma_optimized = sigbeta_new
                
            #self._niter += 1
            #sigmareg_old = sigbeta_new 
            
        self._sigmareg = optimised_betas
        
    def check_convergence(self, x, xold):
        check = False
        tol = self._tolerance
        diff = x - xold
        if diff == 0:
            check = True
        if not check and xold != 0.0:
            diff_percent = abs(diff) / xold
            if diff_percent < tol:
                check = True
            
        return check

