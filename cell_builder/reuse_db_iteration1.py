#!/usr/bin/env python3
"""
ITERATION 1 of the Omega_mid rebuild — REUSE the rigorous lambda_1 lower bounds
already in ryendo/DatabaseTriangle (no FEM).

For every leaf cell of the spectral database (address quadtree) that lies in
Omega_mid, evaluate the CORRECT-vertex lower bound
      J_k(T^p) >= B_k( p_{i,j} ; lambda_1_lb )     (Lemma 4.9),
where p_{i,j}=(x_inf,theta_inf) is the inner corner and lambda_1_lb is the
database's rigorous lower bound for lambda_1 over the cell (hence <= lambda_1 at
the outer corner p_{i+1,j+1}). A cell PASSES iff B1>0 and B2>0; failing cells are
handed to iteration 2 (adaptive subdivision + fresh FEM).
"""
import csv, math, sys
PI=math.pi
C1=PI**2/16.0; C2=7.0*math.sqrt(3.0)*PI**2/12.0
CSTAR=4.0*PI**2/(3.0+math.sqrt(PI*math.sqrt(3.0)))**2
EPS_UP=0.122; EPS_DOWN=0.04; P0=(0.5,math.sqrt(3.0)/2.0)

def yv(x,t): return x*math.tan(t)
def s_omega(x,t): return x*x/math.cos(t)**2
def d2(x,t): y=yv(x,t); return (x-P0[0])**2+(y-P0[1])**2

def classify(c):
    xi,xs,ti,ts=c
    if s_omega(xi,ti)>1.0: return 'discard'
    if yv(xs,ts)<EPS_DOWN: return 'discard'
    cs=[(xi,ti),(xs,ti),(xi,ts),(xs,ts)]; dd=[d2(x,t) for x,t in cs]
    if max(dd)<EPS_UP*EPS_UP: return 'discard'
    if s_omega(xs,ts)<=1.0 and yv(xi,ti)>=EPS_DOWN and min(dd)>=EPS_UP*EPS_UP:
        return 'mid_full'
    return 'straddle'

def Bk(xi,ti,lam):
    y=yv(xi,ti); A=y/2.0
    P=1.0+math.hypot(xi,y)+math.hypot(1.0-xi,y)
    B1=lam*A - C1*P*P/A - C2
    B2=lam*A - CSTAR*(P+math.sqrt(4.0*PI*A))**2/(4.0*A)
    return B1,B2

def main(dbcsv):
    tot=mid_full=straddle=discard=skipnan=0
    npass=nfail=0
    passed=[]; failed=[]; boundary=[]
    worstB=math.inf
    with open(dbcsv) as f:
        for r in csv.DictReader(f):
            tot+=1
            try:
                xi=float(r['x_inf']); xs=float(r['x_sup'])
                ti=float(r['theta_inf']); ts=float(r['theta_sup'])
            except: continue
            cell=(xi,xs,ti,ts)
            cls=classify(cell)
            if cls=='discard': discard+=1; continue
            if cls=='straddle': straddle+=1; boundary.append((r['address'],cell)); continue
            mid_full+=1
            lam_s=r.get('lambda_1_lb','')
            if lam_s in ('','NaN','nan'): skipnan+=1; continue
            lam=float(lam_s)
            B1,B2=Bk(xi,ti,lam)
            if B1>0 and B2>0:
                npass+=1; passed.append((r['address'],cell,B1,B2,lam))
                worstB=min(worstB,B1,B2)
            else:
                nfail+=1; failed.append((r['address'],cell,B1,B2,lam))
    print(f"[iter1] database leaves total           : {tot}")
    print(f"[iter1] in Omega_mid (fully inside)      : {mid_full}")
    print(f"[iter1]   PASS (B1>0 and B2>0)           : {npass}")
    print(f"[iter1]   FAIL (need iteration 2)        : {nfail}")
    print(f"[iter1]   skipped NaN lambda (theta0)    : {skipnan}")
    print(f"[iter1] straddle (boundary, -> refine)   : {straddle}")
    print(f"[iter1] discarded (outside Omega_mid)    : {discard}")
    if passed:
        print(f"[iter1] min B over PASSED cells          : {worstB:.4e}")
    if failed:
        fb=[min(b1,b2) for _,_,b1,b2,_ in failed]
        print(f"[iter1] FAILED cells min B range         : [{min(fb):.3e}, {max(fb):.3e}]")
        xs_f=[c[0] for _,c,*_ in failed]; ys_f=[yv(c[0],c[2]) for _,c,*_ in failed]
        print(f"[iter1] FAILED band x in [{min(xs_f):.3f},{max(xs_f):.3f}], "
              f"y in [{min(ys_f):.3f},{max(ys_f):.3f}]")
    # write outputs
    def wr(path,rows):
        with open(path,'w',newline='') as f:
            w=csv.writer(f); w.writerow(['address','x_inf','x_sup','theta_inf','theta_sup','B1','B2','lambda_1_lb'])
            for a,c,b1,b2,l in rows:
                w.writerow([a,f'{c[0]:.17g}',f'{c[1]:.17g}',f'{c[2]:.17g}',f'{c[3]:.17g}',f'{b1:.6e}',f'{b2:.6e}',f'{l:.10g}'])
    wr('/tmp/iter1_verified.csv',passed)
    wr('/tmp/iter1_failed.csv',failed)
    print("[iter1] wrote /tmp/iter1_verified.csv, /tmp/iter1_failed.csv")

if __name__=='__main__':
    main(sys.argv[1] if len(sys.argv)>1 else '/tmp/DatabaseTriangle/results/database_best_cover.csv')
