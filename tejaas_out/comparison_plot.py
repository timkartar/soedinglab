from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import random
import operator

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
figsize = (5, 5)
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
shuf_cardio_dict = dict()

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
    
    dirc = Path.cwd()/('shuf_cardio_output/chr' + str(chrm))
    pths = [pth for pth in dirc.iterdir() if "_rr.txt" in pth.name]
    for pth in pths:
        l = open(str(pth),"r").readlines()
        rsids = [line.split("\t")[0] for line in l[1:]]
        pvals = [-np.log10(float(line.split("\t")[5])) for line in l[1:]]
        for i in range(len(rsids)):
            shuf_cardio_dict[rsids[i]] = pvals[i]

toplot = []
random_toplot = []
shuf_toplot = []

sig_thres = -np.log10(0.05)

check_x = []
check_random_x = []
check_shuf_x = []

cardio_snps = len(cardio_dict.keys())
cardio_tuples = sorted(cardio_dict.items(), key=operator.itemgetter(1),reverse=True)

shuf_cardio_snps = len(shuf_cardio_dict.keys())
shuf_cardio_tuples = sorted(shuf_cardio_dict.items(), key=operator.itemgetter(1),reverse=True)

random_check_snps = list(cardio_dict.keys())

random.shuffle(random_check_snps)

positives = []
i = 0 

while(i < cardio_snps):
    try:
        if(gtex_dict[cardio_tuples[i][0]] > sig_thres):
            positives.append(i)
    except:
        pass
    i = i + 1
    while(i < cardio_snps and cardio_tuples[i][1] == cardio_tuples[i-1][1]):
        try:
            if(gtex_dict[cardio_tuples[i][0]] > sig_thres):
                positives.append(i)
        except:
            pass
        i = i + 1
    check_x.append(i)
    
    toplot.append(len(positives))


shuf_cardio_positives = []
i = 0

while(i < shuf_cardio_snps):
    try:
        if(gtex_dict[shuf_cardio_tuples[i][0]] > sig_thres):
            shuf_cardio_positives.append(i)
    except:
        pass
    i = i + 1
    while(i < shuf_cardio_snps and shuf_cardio_tuples[i][1] == shuf_cardio_tuples[i-1][1]):
        try:
            if(gtex_dict[shuf_cardio_tuples[i][0]] > sig_thres):
                shuf_cardio_positives.append(i)
        except:
            pass
        i = i + 1
    check_shuf_x.append(i)
    
    shuf_toplot.append(len(shuf_cardio_positives))

random_positives = []
i =0
while(i < cardio_snps):
    try:
        if(gtex_dict[random_check_snps[i]] > sig_thres):
            random_positives.append(i)
    except:
        pass
    i = i + 1 
    while(i < cardio_snps and cardio_dict[random_check_snps[i]] == cardio_dict[random_check_snps[i-1]]):
        try:
            if(gtex_dict[random_check_snps[i]] > sig_thres):
                random_positives.append(i)
        except:
            pass
        i = i + 1
    check_random_x.append(i)
    
    random_toplot.append(len(random_positives))

plt.xlim([1500000, 3000000])
plt.ylim([75000,175000])
plt.plot(check_x,toplot, label='Top Comparison')
plt.plot(check_random_x,random_toplot, label='Random Comparison')
plt.plot(check_shuf_x, shuf_toplot, label='shuffled cardio comparison')
plt.tight_layout()
plt.legend()
plt.xlabel("cardio_neg_log_pval_thres")
plt.ylabel("gtex_validation_count")

plt.savefig("cardio_gtex_comparison.png", bbox_inches='tight')
plt.show()
