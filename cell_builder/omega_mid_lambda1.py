# Fast floating-point lambda1 estimate for triangle (0,0),(1,0),(x, x*tan(theta)).
# P1 FEM on a reference-triangle mesh mapped to the target; Richardson-extrapolated
# in 1/n^2 to approximate the true lambda1.  This is the PRELIMINARY (non-rigorous)
# estimate used only to decide the adaptive cell structure; the rigorous pass
# recomputes everything with INTLAB.
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from functools import lru_cache
import math

def _assemble(ax, ay, n):
    # nodes on unit ref triangle (0,0),(1,0),(0,1), uniform barycentric grid
    idx = {}; pts = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            idx[(i, j)] = len(pts); pts.append((i / n, j / n))
    pts = np.asarray(pts)
    tris = []
    for i in range(n):
        for j in range(n - i):
            tris.append((idx[(i, j)], idx[(i + 1, j)], idx[(i, j + 1)]))
            if j < n - i - 1:
                tris.append((idx[(i + 1, j)], idx[(i + 1, j + 1)], idx[(i, j + 1)]))
    tris = np.asarray(tris)
    # affine map ref->target: (r,s) -> (r*1 + s*ax, s*ay)
    X = pts[:, 0] + pts[:, 1] * ax
    Y = pts[:, 1] * ay
    nodes = np.column_stack([X, Y]); N = len(nodes)

    # ---- vectorized P1 stiffness/mass assembly (COO) ----
    p = nodes[tris]                                   # T x 3 x 2
    v1 = p[:, 1] - p[:, 0]; v2 = p[:, 2] - p[:, 0]    # T x 2
    det = v1[:, 0]*v2[:, 1] - v1[:, 1]*v2[:, 0]       # T
    area = np.abs(det)/2.0
    Binv = np.empty((tris.shape[0], 2, 2))
    Binv[:, 0, 0] = v2[:, 1];  Binv[:, 0, 1] = -v2[:, 0]
    Binv[:, 1, 0] = -v1[:, 1]; Binv[:, 1, 1] = v1[:, 0]
    Binv /= det[:, None, None]
    g0 = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])     # 3 x 2
    g = np.einsum('ij,tjk->tik', g0, Binv)            # T x 3 x 2  (grad of each P1 basis)
    Ke = area[:, None, None] * np.einsum('tik,tjk->tij', g, g)        # T x 3 x 3
    Mloc = np.array([[2., 1., 1.], [1., 2., 1.], [1., 1., 2.]])/12.0
    Me = area[:, None, None] * Mloc[None]             # T x 3 x 3

    rows = np.repeat(tris, 3, axis=1).ravel()         # T*9
    cols = np.tile(tris, (1, 3)).ravel()              # T*9
    K = coo_matrix((Ke.reshape(-1), (rows, cols)), shape=(N, N)).tocsr()
    M = coo_matrix((Me.reshape(-1), (rows, cols)), shape=(N, N)).tocsr()

    onb = [k for k, (r, s) in enumerate(pts)
           if r < 1e-12 or s < 1e-12 or abs(r + s - 1) < 1e-12]
    free = np.array(sorted(set(range(N)) - set(onb)))
    return K[free][:, free], M[free][:, free]

def _lam1_n(ax, ay, n):
    K, M = _assemble(ax, ay, n)
    w = eigsh(K, k=1, M=M, sigma=0, which='LM', return_eigenvectors=False)
    return float(w[0])

def lam1(x, theta, base_n=48, richardson=True):
    """Approximate true lambda1 of triangle (0,0),(1,0),(x, x tan theta).
    P1 gives an upper bound ~ lam_true + C/n^2; Richardson with (n,2n) removes
    the leading 1/n^2 term. Mesh is refined for thin (small-theta) triangles."""
    ax = x; ay = x * math.tan(theta)
    aspect = max(ax, 1.0) / max(ay, 1e-9)           # base/height
    n = int(base_n * max(1.0, math.sqrt(aspect) / 3.0))
    n = min(max(n, 24), 220)
    l1 = _lam1_n(ax, ay, n)
    if not richardson:
        return l1
    l2 = _lam1_n(ax, ay, 2 * n)
    # P1 error ~ C/n^2 : extrapolate  lam ~ (4*l2 - l1)/3
    return (4.0 * l2 - l1) / 3.0

if __name__ == '__main__':
    # validate against the rigorous lower bounds found earlier
    tests = [
        # (x, theta, expected approx lambda1, note)
        (0.5, math.atan(math.sqrt(3)), 16*math.pi**2/3, 'equilateral -> 16pi^2/3=52.638'),
        (0.5640625, 0.9225504643977442, 61.7, 'cell 73917 outer corner (lam_low=61.644)'),
        (0.80, math.atan(0.04007/0.80), 7155.0, 'thin cell 69026 inner corner (~7155)'),
    ]
    for x, th, exp, note in tests:
        est = lam1(x, th)
        print(f"x={x:.5f} th={th:.5f}  lam1_est={est:10.3f}   expected~{exp:9.3f}   [{note}]")
