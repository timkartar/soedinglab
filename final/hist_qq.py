import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import chi2
from scipy.stats import uniform
from scipy import special

def normal_quantile_plot(data, mu = 0, sigma = 1):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = mu + sigma * np.sqrt(2) * special.erfinv(2 * mq - 1) # qunatile function for normal distribution (see wiki)
    return x, y

def uniform_quantile_plot_old(data, a = 0, b = 1):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = a + mq * (b - a) # quantile function for uniform distribution: a + p(b - a)
    return x, y

def uniform_quantile_plot(data, a = 0, b = 1):
    y = data[np.argsort(data)]
    xrand = np.random.uniform(0, 1, size=y.shape[0])
    x = xrand[np.argsort(xrand)]
    return x, y


def chisquare_quantile_plot(data, df = 1):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = chi2.ppf(mq, df) # quantile function for chisquare distribution from python
    return x, y

def get_pdf(dist, x, df = 1, loc = 0, scale = 1):
    if dist == 'normal':
        y = norm.pdf(x, loc = loc, scale = scale)
    elif dist == 'uniform':
        y = uniform.pdf(x, loc = loc, scale = scale)
    elif dist == 'chi2':
        y = chi2.pdf(x, df)
    else:
        print ("No distribution found with name {:s}".format(dist))
        y = np.zeros_like(x)
    return y

def get_quantile(dist, data, df = 1, loc = 0, scale = 1):
    if dist == 'normal':
        x, y = normal_quantile_plot(data)
    elif dist == 'uniform':
        x, y = uniform_quantile_plot(data)
        x = - np.log10(x)
        y = - np.log10(y)
    elif dist == 'chi2':
        x, y = chisquare_quantile_plot(data, df = df)
    else:
        print ("No distribution found with name {:s}".format(dist))
        y = np.zeros_like(x)
    return x, y

def plot(ax1, ax2, data, nbins, dist, df = 0, loc = 0, scale = 1, size = 1):

    axisfontsize = 20 * size
    labelfontsize = 15 * size
    padwidth = 10 * size
    msize = 10 * size
    wsize = 4 * size
    ticklen = 5 * size
    borderwidth = 2 * size
    bordercolor = 'black'

    banskt_colors_hex = [
        '#2D69C4', # blue 
        '#FFB300', # Vivid Yellow
        '#93AA00', # Vivid Yellowish Green
        '#CC2529', # red
        '#535154', # gray
        '#6B4C9A', # purple
        '#922428', # dark brown
        '#948B3D', # olive
        ]
    colors = banskt_colors_hex
    

    xmin = min(0, np.min(data))
    xmax = max(1, np.max(data))
    bins = np.linspace(xmin, xmax, nbins)
    xvals = [(bins[i] + bins[i+1]) / 2 for i in range(nbins - 1)]
    h, _ = np.histogram(data, bins=bins, density=True)
    ax1.fill_between(xvals, h, 0, color=colors[0], alpha = 0.2)

    xvals = np.linspace(xmin, xmax, data.shape[0])
    yvals = get_pdf(dist, xvals, df = df, loc = loc, scale = scale)
    ax1.plot(xvals, yvals, lw = wsize, color=colors[0])

    x, y = get_quantile(dist, data, df = df, loc = loc, scale = scale)
    xmin = min(np.min(x), np.min(y))
    xmax = max(np.max(x), np.max(y))
    ax2.scatter(x, y, s = msize, color=colors[3], alpha = 0.5)
    ax2.plot([xmin, xmax], [xmin, xmax], lw = wsize / 4, ls = 'dashed', color=colors[4])

    ax1.set_xlabel('x', {'size': axisfontsize}, labelpad = padwidth)
    ax1.set_ylabel('PDF', {'size': axisfontsize}, labelpad = padwidth)
    
    ax2.set_xlabel('Expected', {'size': axisfontsize}, labelpad = padwidth)
    ax2.set_ylabel('Observed', {'size': axisfontsize}, labelpad = padwidth)

    xticks = ax2.get_xticks()
    ax2.set_yticks(xticks)
    ax2.set_xticks(xticks)
    for ax in [ax1, ax2]:
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        ax.tick_params(axis='both', which = 'major',
                       length = ticklen, width = borderwidth, pad=padwidth,
                       labelsize = labelfontsize,
                       color = bordercolor,
                       labelcolor = bordercolor,
                       bottom = True, top = False, left = True, right = False,
                      )
        for side, border in ax.spines.items():
            border.set_linewidth(borderwidth)
            border.set_color(bordercolor)
            
    font_properties = {'family':'sans-serif', 'weight': 'bold', 'size': labelfontsize}
    
    return None
