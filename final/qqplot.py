import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from scipy.stats import norm
from scipy.stats import uniform
from scipy.stats.mstats import mquantiles

def normal_quantile_plot(data, mu = 0, sigma = 1):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = mu + sigma * np.sqrt(2) * special.erfinv(2 * mq - 1) # qunatile function for normal distribution (see wiki)
    return x, y

def uniform_quantile_plot(data, a = 0, b = 1):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = a + mq * (b - a) # quantile function for uniform distribution: a + p(b - a)
    return x, y

def plot_hist_qq(ax1, ax2, data, nbins, dist, loc = 0, scale = 1):
    n = data.shape[0]
    xmin = np.min(data)
    xmax = np.max(data)
    bins = np.linspace(xmin, xmax, nbins)
    xvals = [(bins[i] + bins[i+1]) / 2 for i in range(nbins - 1)]
    h, _ = np.histogram(data, bins=bins, density=True)
    ax1.fill_between(xvals, h, 0, color='blue', alpha = 0.2)

    xvals = np.linspace(xmin, xmax, n)
    if dist == 'normal':
        yvals = norm.pdf(xvals, loc = loc, scale = scale)
        x, y = normal_quantile_plot(data)
        xmin = -10
        xmax = 10
    elif dist == 'uniform':
        yvals = uniform.pdf(xvals, loc = loc, scale = scale)
        x, y = uniform_quantile_plot(data)
        x = - np.log10(x)
        y = - np.log10(y)
        xmin = -np.log10(xmin)
        xmax = -np.log10(xmax)
    else:
        print ("No distribution found with name {:s}".format(dist))
    
    ax1.plot(xvals, yvals, lw = 4, color='blue')

    ax2.plot([xmin, xmax], [xmin, xmax], lw = 2, ls = 'dashed', color='gray')
    ax2.scatter(x, y, s = 10, color='red', alpha = 0.5)
    
    ax1.set_xlabel('x', {'size': 15}, labelpad = 10)
    ax1.set_ylabel('Probability density function', {'size': 15}, labelpad = 10)
    
    ax2.set_xlabel('Expected', {'size': 15}, labelpad = 10)
    ax2.set_ylabel('Observed', {'size': 15}, labelpad = 10)
    
    return None
    
'''    
## ============ MAIN ==============================

fig = plt.figure(figsize = (10, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

n = 100000

## ========== Normal distribution ================
mu = 0
sigma = 1
Y = np.random.normal(mu, sigma, n)
nbins = 50
plot_hist_qq(ax1, ax2, Y, nbins, 'normal', loc = mu, scale = sigma)

## ========== Uniform distribution ================
a = 0
b = 1
Y = np.random.uniform(a, b, n)
nbins = 50
plot_hist_qq(ax3, ax4, Y, nbins, 'uniform', loc = a, scale = b)


plt.tight_layout()
plt.show()
'''
