#!/usr/bin/env python3
"""Summarize + visualize a rebuilt Omega_mid cell_def (address-quadtree)."""
import csv, math, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon

path = sys.argv[1] if len(sys.argv)>1 else '/tmp/cell_def_omega_mid_rebuilt.csv'
cells=[]
with open(path) as f:
    for r in csv.DictReader(f):
        cells.append((r['address'], float(r['x_inf']), float(r['x_sup']),
                      float(r['theta_inf']), float(r['theta_sup'])))
levels=[len(a) for a,*_ in cells]
print(f"cells: {len(cells)}   levels min/median/max: "
      f"{min(levels)}/{int(np.median(levels))}/{max(levels)}")
from collections import Counter
for lv,c in sorted(Counter(levels).items()):
    print(f"  level {lv:2d}: {c:8d} cells")

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(15,7))
# (x,theta) plane
rects=[Rectangle((xi,ti),xs-xi,ts-ti) for _,xi,xs,ti,ts in cells]
pc=PatchCollection(rects,facecolor='none',edgecolor='b',linewidth=0.05)
ax1.add_collection(pc); ax1.set_xlim(0.5,1.0); ax1.set_ylim(0,math.pi/3)
ax1.set_xlabel('x'); ax1.set_ylabel('theta'); ax1.set_title(f'Omega_mid cells in (x,theta)  [{len(cells)} cells]')
# (x,y) plane: each cell -> quad with y=x tan(theta)
polys=[]
for _,xi,xs,ti,ts in cells:
    pts=[(xi,xi*math.tan(ti)),(xs,xs*math.tan(ti)),(xs,xs*math.tan(ts)),(xi,xi*math.tan(ts))]
    polys.append(Polygon(pts,closed=True))
pc2=PatchCollection(polys,facecolor='none',edgecolor='g',linewidth=0.05)
ax2.add_collection(pc2)
th=np.linspace(0,math.pi/2,200); ax2.plot(np.cos(th),np.sin(th),'r-',lw=0.8,label='x^2+y^2=1')
ax2.plot([0.5,0.5],[0,math.sqrt(3)/2],'k--',lw=0.5)
ax2.scatter([0.5],[math.sqrt(3)/2],c='r',s=20,zorder=5,label='p0 (equilateral)')
ax2.axhline(0.04,color='orange',lw=0.6,ls=':',label='y=eps_down')
ax2.set_xlim(0.45,1.0); ax2.set_ylim(0,0.9); ax2.set_xlabel('x'); ax2.set_ylabel('y')
ax2.set_title('mapped triangles apex (x,y)'); ax2.legend(fontsize=8)
plt.tight_layout(); plt.savefig('/tmp/omega_mid_cells.png',dpi=130)
print("wrote /tmp/omega_mid_cells.png")
