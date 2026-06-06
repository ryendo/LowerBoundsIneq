import csv, math, sys
P0=(0.5,math.sqrt(3)/2); EPS_UP=0.122
def yv(x,t): return x*math.tan(t)
def mind(c):
    xi,xs,ti,ts=c; return min(math.hypot(x-0.5,yv(x,t)-P0[1]) for x,t in ((xi,ti),(xs,ti),(xi,ts),(xs,ts)))
rows=list(csv.DictReader(open('/tmp/cell_def_omega_mid_rebuilt.csv')))
print(f"cells: {len(rows)}")
for DLG in (0.135,0.14,0.15,0.18,0.22,0.26):
    nLG=sum(1 for r in rows if mind((float(r['x_inf']),float(r['x_sup']),float(r['theta_inf']),float(r['theta_sup'])))<DLG)
    print(f"  DIST_LG={DLG:.3f}: LG starts {nLG:7d} ({100*nLG/len(rows):5.1f}%), CR starts {len(rows)-nLG:7d} ({100*(len(rows)-nLG)/len(rows):5.1f}%)")
# write with chosen DIST_LG (arg) and CR mesh gradation
DLG=float(sys.argv[1]) if len(sys.argv)>1 else 0.14
MESH_CR=0.02; MESH_CR_FINE=0.01; MESH_LG=0.10; ORD_LG=2
hdr=['address','x_inf','x_sup','theta_inf','theta_sup','mesh_size_lower_cr','isLG','mesh_size_lower_LG','fem_order_lower_LG']
nLG=0
with open('/tmp/cell_def_omega_mid_rebuilt.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(hdr)
    for r in rows:
        c=(float(r['x_inf']),float(r['x_sup']),float(r['theta_inf']),float(r['theta_sup'])); d=mind(c)
        if d<DLG: isLG=1; mcr=MESH_CR_FINE; nLG+=1
        else:     isLG=0; mcr=(MESH_CR_FINE if d<DLG+0.12 else MESH_CR)
        w.writerow([r['address'],r['x_inf'],r['x_sup'],r['theta_inf'],r['theta_sup'],f'{mcr:g}',isLG,f'{MESH_LG:g}',ORD_LG])
print(f"\nwrote cell_def with DIST_LG={DLG}: LG starts {nLG} ({100*nLG/len(rows):.1f}%), CR-only starts {len(rows)-nLG} ({100*(len(rows)-nLG)/len(rows):.1f}%)")
print("(rigorous pass: CR-start cells escalate to LG only if CR fails; LG-start = immediate disk neighborhood)")
