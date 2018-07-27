from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

fig = plt.figure()

yvals = []
for i in range(22):
    chrm = i + 1
    dirc = Path.cwd()/('gtex_maf_01_output/chr' + str(chrm))
    pths = [pth for pth in dirc.iterdir() if "_rr.txt" in pth.name]
    for pth in pths:
        l = open(str(pth),"r").readlines()
        try:
            yvals = yvals + [float(line.split("\t")[5]) for line in l[1:]]
        except:
            print(pth)
    #print(chrm)
#print("\n".join([str(x) for x in yvals]))
#np.save("maf.npy",yvals)
yvals = np.array(yvals)
plt.hist(yvals[~np.isnan(yvals)], bins=400)
plt.show()
