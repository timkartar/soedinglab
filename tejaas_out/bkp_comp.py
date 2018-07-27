from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import random

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
#ax.axhline(y=-np.log10(0.00000005), ls='dashed', color = 'darkgrey', lw=2)

colors = kch#[kch[1], kch[2], kch[19]]
bgcolors = ['#FFB300','#00538A','#007D34','#B32851','#232C16' ]
lastbp = 0
xlabels = list()
xlabels_pos = list()
gtex_dict = dict()
cardio_dict = dict()
for i in tqdm(range(22)):
    chrm = i + 1
    dirc = Path.cwd()/('output/chr' + str(chrm))
    pths = [pth for pth in dirc.iterdir() if "_rr.txt" in pth.name]
    for pth in pths:
        l = open(str(pth),"r").readlines()
        rsids = [line.split("\t")[0] for line in l[1:]]
        pvals = [-np.log10(float(line.split("\t")[5])) for line in l[1:]]
        for i in range(len(rsids)):
            gtex_dict[rsids[i]] = pvals[i]

    dirc = Path.cwd()/('cardio_output/chr' + str(chrm))
    pths = [pth for pth in dirc.iterdir() if "_rr.txt" in pth.name]
    for pth in pths:
        l = open(str(pth),"r").readlines()
        rsids = [line.split("\t")[0] for line in l[1:]]
        pvals = [-np.log10(float(line.split("\t")[5])) for line in l[1:]]
        for i in range(len(rsids)):
            cardio_dict[rsids[i]] = pvals[i]

print(len(gtex_dict.keys()))
toplot = []
sig_thres = -np.log10(0.05)
thress = [(15 - 0.2*i) for i in range(38)] + [-np.log10(0.00000005)]
#x = thress
x      = [-x for x in thress]
random_toplot = []
checkset_sizes = []
for thres in thress:
    check_snps = [key for key, val in cardio_dict.items() if val > thres]
    checkset_sizes.append(len(check_snps))
    random_check_snps = random.sample(list(cardio_dict.keys()), len(check_snps))
    positives = []
    random_positives = []
    for i in check_snps:
        try:
            if(gtex_dict[i] > sig_thres):
                positives.append(i)
        except:
            pass

    for i in random_check_snps:
        try:
            if(gtex_dict[i] > sig_thres):
                random_positives.append(i)
        except:
            pass

    toplot.append(len(positives))
    random_toplot.append(len(random_positives))

print(checkset_sizes)
print(toplot)
print(random_toplot)
plt.plot(x,toplot, label='Top Comparison')
plt.plot(x,random_toplot, label='Random Comparison')
plt.tight_layout()
plt.legend()
plt.xlabel("cardio_neg_log_pval_thres")
plt.ylabel("gtex_validation_count")
plt.savefig("cardio_gtex_comparison.png", bbox_inches='tight')
plt.show()
