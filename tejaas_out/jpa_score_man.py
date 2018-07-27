from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

kelly_colors_hex = [
    '#FFB300', # Vivid Yellow
    '#803E75', # Strong Purple
    '#FF6800', # Vivid Orange
    '#A6BDD7', # Very Light Blue
    '#C10020', # Vivid Red
    '#CEA262', # Grayish Yellow
    '#817066', # Medium Gray

    # The following don't work well for people with defective color vision
    '#007D34', # Vivid Green
    '#F6768E', # Strong Purplish Pink
    '#00538A', # Strong Blue
    '#FF7A5C', # Strong Yellowish Pink
    '#53377A', # Strong Violet
    '#FF8E00', # Vivid Orange Yellow
    '#B32851', # Strong Purplish Red
    '#F4C800', # Vivid Greenish Yellow
    '#7F180D', # Strong Reddish Brown
    '#93AA00', # Vivid Yellowish Green
    '#593315', # Deep Yellowish Brown
    '#F13A13', # Vivid Reddish Orange
    '#232C16', # Dark Olive Green
    ]
kch = kelly_colors_hex #['red' for i in range(200)]
figsize = (20, 5)
axis_font_size = 15
label_font_size = 20
legend_font_size = 15
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)
#ax.axhline(y=-np.log10(0.00024), ls='dashed', color = 'darkgrey', lw=2)
#ax.axhline(y=-np.log10(0.00000005), ls='dashed', color = 'darkgrey', lw=2)

colors = kch#[kch[1], kch[2], kch[19]]
bgcolors = ['#FFB300','#00538A','#007D34','#B32851','#232C16' ]
lastbp = 0
xlabels = list()
xlabels_pos = list()
for i in range(22):#tqdm(range(22)):
    chrm = i + 1
    dirc = Path.cwd()/('cardio_output/chr' + str(chrm))
    pths = [pth for pth in dirc.iterdir() if "_rr.txt" in pth.name]
    #jpapths = [pth for pth in dirc.iterdir() if "_jpa.txt" in pth.name]
    xvals = [] #[lastbp + x.bp_location for x in unused[i]]
    yvals = [] # [-np.log10(x.p) for x in unused[i]]
    l = []
    bpmax = 0
    for pth in pths:
        l = open(str(pth),"r").readlines()
        ljpa = open(str(pth).replace("rr","jpa"),"r").readlines()
        pos = [int(line.split("\t")[1]) for line in l[1:]]
        tempmax = max(pos)
        xvals = xvals + [lastbp + x for x in pos] #int(line.split("\t")[1]) for line in l[1:]]
        bpmax = max([tempmax,bpmax])
        lpvals = [float(line.split("\t")[1]) for line in ljpa[1:]]
        yvals = yvals + [x for x in lpvals]
        #for i,lpval in enumerate(lpvals):
        #    if(lpval > 10):
        #        print(chrm, l[1:][i].split("\t")[0], lpval)
    ax.scatter(xvals, yvals, color=bgcolors[i % len(bgcolors)], s = 2, alpha=0.8)
    
    #bpmax = int(l[-1].split("\t")[1])
    xlabels.append('%i' % (i + 1))
    xlabels_pos.append(lastbp + bpmax / 2 )
    
    #loci = [x for x in newloci if x[0].chrm == chrm]
    #for j, locus in enumerate(loci):
    #    bp = [lastbp + x.bp_location for x in locus]
    #    p = [min(-np.log10(x.p), 9.8) for x in locus]
    #    ax.scatter(bp, p, color = colors[j % len(colors)], s = 10, alpha = 0.6)
        
    lastbp += bpmax
    #print(lastbp)

ax.set_xticks(xlabels_pos)
ax.set_xticklabels(xlabels)
ax.set_xlim([0, lastbp])
ax.set_ylim([0, max(yvals)])
ax.set_xlabel('Chromosome', {'size': axis_font_size}, labelpad = 15)
ax.set_ylabel('JPA score', {'size': axis_font_size}, labelpad = 15)
#ax.set_ylabel(r'$-\log_{10}(p)$', {'size': axis_font_size}, labelpad = 15)

plt.tight_layout()
plt.savefig("cardio_jpa_manhattan.png", bbox_inches='tight')
plt.show()
