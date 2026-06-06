#!/usr/bin/env python3
"""
Integrated adaptive Omega_mid cell builder (address quadtree), CORRECT vertex.

Iteration 1 : reuse rigorous lambda_1_lb from ryendo/DatabaseTriangle for the
              database leaf cells (no FEM).
Iteration 2+: any cell that fails J at the correct vertex is subdivided into its
              4 address-children; the children's lambda_1 is estimated by the fast
              float FEM table (preliminary). Recurse until both J1,J2 pass or
              max_depth. The rigorous INTLAB pass later recomputes every leaf.

Lower bound used (Lemma 4.9):  J_k(T^p) >= B_k( p_{i,j} ; lambda_1_lb_cell ),
geometry at the inner corner p_{i,j}=(x_inf,theta_inf), eigenvalue lower bound
valid over the whole cell (so <= lambda_1 at the outer corner p_{i+1,j+1}).
"""
import csv, math, sys, time
from multiprocessing import Pool
from omega_mid_lambda1 import _lam1_n

PI=math.pi
C1=PI**2/16.0; C2=7.0*math.sqrt(3.0)*PI**2/12.0
CSTAR=4.0*PI**2/(3.0+math.sqrt(PI*math.sqrt(3.0)))**2
EPS_UP=0.122; EPS_DOWN=0.04; P0=(0.5,math.sqrt(3.0)/2.0)
THETA_MAX=PI/3.0
HAIRCUT=0.0            # accurate where it matters (see lam1_float); DB lb used as-is
MAX_DEPTH=18
Y_FLOOR=EPS_DOWN                       # degenerate-side floor = paper's eps_down
Y_TOP=math.sqrt(3.0)/2.0 - EPS_UP      # equilateral-side ceiling = Omega_up strip floor
STRADDLE_MAXLEN=13                     # cap boundary (straddle) refinement at this address length
DB_CSV='/tmp/DatabaseTriangle/results/database_best_cover.csv'

# ---- aspect-aware float lambda1 for iteration>=2 refined cells ----
# An interpolation table was tried first but its ~1% noise exceeds the small J
# near the disk boundary -> spurious non-convergence. Direct FEM instead, but:
#   * near-equilateral cells (aspect<3, small J): Richardson (2 solves) for accuracy;
#   * thin/near-degenerate cells (aspect>=3, large J): one coarse solve is plenty,
#     since B is hugely positive there and high-n FEM would be the bottleneck.
def lam1_float(x, theta):
    y = x*math.tan(theta); aspect = max(x,1.0)/max(y,1e-9)
    if aspect < 3.0:
        l1 = _lam1_n(x, y, 44); l2 = _lam1_n(x, y, 88)
        return (4.0*l2 - l1)/3.0
    n = min(int(26*math.sqrt(aspect)/2.0), 64); n = max(n, 24)
    return _lam1_n(x, y, n)

def yv(x,t): return x*math.tan(t)
def s_omega(x,t): return x*x/math.cos(t)**2
def d2(x,t): y=yv(x,t); return (x-P0[0])**2+(y-P0[1])**2

def classify(c):
    # The builder covers the SLAB  Y_FLOOR <= y <= Y_TOP  inside Omega:
    #   Y_FLOOR = eps_down            (below -> Omega_down, handled analytically),
    #   Y_TOP   = sqrt3/2 - eps_up    (above -> Omega_up; the upper boundary of
    #             Omega_mid is the CIRCLE dist(.,p0)=eps_up, which lies entirely at
    #             y >= Y_TOP, hence inside Omega_up. So Omega_up covers the whole
    #             circular upper part of Omega_mid AND the small-J disk interior;
    #             our cells never need to approach p0 -> no blow-up).
    # A cell is DISCARDed only if it has NO Omega_mid portion (fully below eps_down,
    # fully above Y_TOP, or fully outside Omega). Every other cell -- INCLUDING those
    # straddling y=eps_down, y=Y_TOP or the Omega arc -- has an Omega_mid portion and
    # must be VERIFIED (eval_cell evaluates B over that portion, clipped to y>=eps_down),
    # so there are no unverified boundary slivers (rounding-level coverage).
    xi,xs,ti,ts=c
    y_min=yv(xi,ti); y_max=yv(xs,ts)
    if y_max < Y_FLOOR:      return 'discard'    # fully below eps_down -> Omega_down
    if y_min > Y_TOP:        return 'discard'    # fully above Y_TOP    -> Omega_up
    if s_omega(xi,ti) > 1.0: return 'discard'    # fully outside Omega (inner corner already out)
    return 'verify'

def Bk(x,t,lam):
    y=yv(x,t); A=y/2.0
    P=1.0+math.hypot(x,y)+math.hypot(1.0-x,y)
    return (lam*A-C1*P*P/A-C2, lam*A-CSTAR*(P+math.sqrt(4.0*PI*A))**2/(4.0*A))

def children(cell,addr):
    xi,xs,ti,ts=cell; xm=0.5*(xi+xs); tm=0.5*(ti+ts)
    return [(addr+'1',(xm,xs,tm,ts)),(addr+'2',(xi,xm,tm,ts)),
            (addr+'3',(xi,xm,ti,tm)),(addr+'4',(xm,xs,ti,tm))]

def portion_eval_points(cell):
    # Lower-left boundary of the cell's Omega_mid portion (y >= eps_down), where
    # B (increasing in x and theta, Lemma 4.9) attains its infimum over the portion.
    # - interior / top / arc cells: the inner corner (x_inf,theta_inf) is the min;
    # - bottom-straddle cells (y_inner < eps_down): along y = eps_down the area is
    #   fixed and the perimeter grows with x, so B decreases in x -> sample the
    #   y=eps_down curve at several x and take the worst.
    xi,xs,ti,ts=cell; y_inner=yv(xi,ti); pts=[]
    if y_inner >= Y_FLOOR:
        pts.append((xi,ti))
    if y_inner < Y_FLOOR < yv(xs,ts):
        for xx in (xi,0.5*(xi+xs),xs):
            th=math.atan(Y_FLOOR/xx)
            if ti-1e-15 <= th <= ts+1e-15:
                pts.append((xx,min(max(th,ti),ts)))
    if not pts:
        pts.append((xi,min(max(math.atan(Y_FLOOR/xi),ti),ts)))
    return pts

def eval_cell(item):
    addr,cell,db_lam=item
    if classify(cell)=='discard':
        return {'addr':addr,'cell':cell,'status':'discard'}
    xs,ts=cell[1],cell[3]
    if db_lam is not None:               # iteration 1: rigorous DB lower bound
        lam_use=db_lam; src='db'
    else:                                # iteration >=2: float estimate at outer corner
        lam_use=lam1_float(xs,ts)*(1.0-HAIRCUT); src='fem'
    B1=B2=math.inf
    for (x,t) in portion_eval_points(cell):
        if s_omega(x,t) > 1.0+1e-12: continue     # eval point not in Omega -> skip
        b1,b2=Bk(x,t,lam_use); B1=min(B1,b1); B2=min(B2,b2)
    ok=(B1>0.0)and(B2>0.0)
    return {'addr':addr,'cell':cell,'status':'leaf' if ok else 'split',
            'B1':B1,'B2':B2,'src':src}

def load_db_seed():
    seed=[]
    with open(DB_CSV) as f:
        for r in csv.DictReader(f):
            try:
                cell=(float(r['x_inf']),float(r['x_sup']),float(r['theta_inf']),float(r['theta_sup']))
            except: continue
            lam=r.get('lambda_1_lb','')
            db_lam=None if lam in ('','NaN','nan') else float(lam)
            seed.append((r['address'],cell,db_lam))
    return seed

def build(nproc=None,max_depth=MAX_DEPTH):
    frontier=load_db_seed()
    print(f"[build] seeded {len(frontier)} database leaf cells (iteration 1)",flush=True)
    leaves=[]; boundary=[]; unverified=[]; discarded=0
    pool=Pool(nproc); depth=0
    while frontier and depth<=max_depth:
        t0=time.time()
        res=pool.map(eval_cell,frontier,chunksize=max(1,len(frontier)//(8*(nproc or 8))+1))
        nxt=[]; nl=ns=nd=0; ndb=0
        for r in res:
            if r['status']=='discard': nd+=1; discarded+=1
            elif r['status']=='leaf':
                leaves.append(r); nl+=1; ndb+= (r.get('src')=='db')
            else:                                   # split (B not yet positive)
                ns+=1
                if depth<max_depth:
                    nxt.extend([(a,c,None) for a,c in children(r['cell'],r['addr'])])
                else:
                    unverified.append(r)            # still failing at max depth
        tag='iter1(DB)' if depth==0 else f'iter{depth+1}'
        print(f"  {tag:10s} in={len(frontier):8d} leaf={nl:7d}(db={ndb:6d}) "
              f"split={ns:7d} discard={nd:7d} ({time.time()-t0:.1f}s)",flush=True)
        frontier=nxt; depth+=1
    pool.close(); pool.join()
    return leaves,boundary,unverified,discarded

def write_cell_def(rows,path):
    hdr=['address','x_inf','x_sup','theta_inf','theta_sup',
         'mesh_size_lower_cr','isLG','mesh_size_lower_LG','fem_order_lower_LG']
    with open(path,'w',newline='') as f:
        w=csv.writer(f); w.writerow(hdr)
        for r in rows:
            xi,xs,ti,ts=r['cell']
            w.writerow([r['addr'],f'{xi:.17g}',f'{xs:.17g}',f'{ti:.17g}',f'{ts:.17g}',
                        '0.04','1','0.1249','2'])

if __name__=='__main__':
    nproc=int(sys.argv[1]) if len(sys.argv)>1 else None
    md=int(sys.argv[2]) if len(sys.argv)>2 else MAX_DEPTH
    t0=time.time()
    print(f"[build] integrated Omega_mid build  eps_up={EPS_UP} eps_down={EPS_DOWN} "
          f"haircut={HAIRCUT} max_depth={md} nproc={nproc}",flush=True)
    leaves,boundary,unverified,discarded=build(nproc,md)
    ndb=sum(1 for r in leaves if r.get('src')=='db')
    print(f"[build] DONE {time.time()-t0:.1f}s",flush=True)
    print(f"[build] verified leaf cells           : {len(leaves)}  (DB-reused {ndb}, FEM-refined {len(leaves)-ndb})")
    print(f"[build] boundary slivers (max depth)  : {len(boundary)}")
    print(f"[build] UNVERIFIED at max depth       : {len(unverified)}")
    print(f"[build] discarded outside Omega_mid   : {discarded}")
    if leaves:
        d=sorted(len(r['addr']) for r in leaves)
        print(f"[build] leaf level min/median/max     : {d[0]}/{d[len(d)//2]}/{d[-1]}")
    write_cell_def(leaves,'/tmp/cell_def_omega_mid_rebuilt.csv')
    print(f"[build] wrote /tmp/cell_def_omega_mid_rebuilt.csv ({len(leaves)} cells)",flush=True)
    if unverified: write_cell_def(unverified,'/tmp/cell_def_unverified.csv')
    if boundary:   write_cell_def(boundary,'/tmp/cell_def_boundary.csv')
