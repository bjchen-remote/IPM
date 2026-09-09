"""Standalone scientific figure from saved curves; no numerical field recomputation."""
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root=Path(sys.argv[1]).resolve()
report=json.loads((root/"report.json").read_text())
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
fig,axes=plt.subplots(2,2,figsize=(10.5,7),layout="constrained")
colors=["#2166ac","#b2182b","#238b45"]
curves={}
for k,c in enumerate(report["cases"]):
    for direction in (1,2):
        z=np.loadtxt(root/f'{c["label"]}_{direction}_curve.csv',delimiter=',');curves[k,direction]=z
        label=f'{c["label"][0]}: H{[8,16,8][k]}, target {[32,32,42][k]}, tau={c["canonicalTime"]:.4f}'
        axes[0,direction-1].plot(z[:,0],z[:,1],color=colors[k],lw=1.8,label=label)
        axes[0,direction-1].plot(z[:,0],z[:,2],color=colors[k],lw=.8,ls='--',alpha=.7)
for j,(a,b) in enumerate([(0,1),(0,2),(1,2)]):
    for direction in (1,2):
        aa,bb=curves[a,direction],curves[b,direction]
        assert np.array_equal(aa[:,0],bb[:,0])
        ax=axes[1,direction-1]
        ax.plot(aa[:,0],bb[:,1]-aa[:,1],color=colors[j],lw=1.5,label=f'{"ABC"[b]} - {"ABC"[a]}')
        ax.plot(aa[:,0],bb[:,2]-aa[:,2],color=colors[j],lw=.7,ls='--',alpha=.65)
for direction,title in enumerate(["Wall: x = a + wall-width * xi, y = 0", "Vertical: x = a, y = vertical-width * eta"]):
    axes[0,direction].set_title(title)
    axes[0,direction].set_ylabel("rho_x / continuous peak")
    axes[1,direction].set_ylabel("Aligned shape difference")
    axes[1,direction].axhline(0,color=".5",lw=.6)
    for row in (0,1):
        axes[row,direction].set_xlabel("xi" if direction==0 else "eta")
        axes[row,direction].grid(alpha=.2)
    axes[0,direction].legend(fontsize=8,loc="best")
    axes[1,direction].legend(fontsize=9)
fig.suptitle(f'Actual common physical time t = {report["samePhysicalEndpoint"]:.12f}\nOwn C1 peak + connected 90% widths; solid Hermite / dashed linear',fontsize=13)
fig.savefig(root/"ABC_aligned_inner_shapes.png",dpi=180)
fig.savefig(root/"ABC_aligned_inner_shapes.pdf")
print(root/"ABC_aligned_inner_shapes.png")
