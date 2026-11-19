"""
Independent check of Impedora's slot solver (slotInductance in src/physics/fieldSolver2d.ts).

A zero-thickness trace of width w at height h over a plane with a slot of half-width s centred under
it, and a second plane `Hc` below the first, all in air. Second-order finite elements (scikit-fem) on
an unstructured Gmsh mesh refined at the trace and slot edges; returns L' = mu0*eps0 / C'air.

Usage:  python3 slot_fem.py <w mm> <h mm> <Hc mm> <s mm> [<s mm> ...]
Example (case ms2):  python3 slot_fem.py 0.19 0.1 0.535 0.05 0.1 0.2 0.3 0.4
Needs gmsh and scikit-fem (pip install gmsh scikit-fem).
"""
import sys, json, numpy as np, gmsh
from skfem import MeshTri, Basis, ElementTriP2, BilinearForm, condense, solve, asm
from skfem.helpers import dot, grad
EPS0 = 8.8541878128e-12; MU0 = 4e-7 * np.pi
def lprime(w, h, s, Hc, top=None, lc_scale=1.0):
    # half domain x >= 0; units: micrometres
    U = 1e6; w, h, s, Hc = w * U, h * U, s * U, Hc * U
    X = 60 * max(h, w, s) + 20 * Hc
    Yt = (top * U) if top else h + 60 * max(h, w) + 20 * Hc
    gmsh.initialize(); gmsh.option.setNumber('General.Terminal', 0); gmsh.model.add('s')
    geo = gmsh.model.geo
    lcf = min(w, h, s if s > 0 else 1e9) / 6 * lc_scale
    lcF = X / 8
    pts = {}
    def P(x, y, lc): 
        k = (round(x, 6), round(y, 6))
        if k not in pts: pts[k] = geo.addPoint(x, y, 0, lc)
        return pts[k]
    # outer boundary polygon (counter-clockwise): (0,-Hc) -> (X,-Hc) -> (X,Yt) -> (0,Yt) -> back down the axis
    a0 = P(0, -Hc, lcF); a1 = P(X, -Hc, lcF); a2 = P(X, Yt, lcF); a3 = P(0, Yt, lcF)
    # axis points: trace at (0,h), slot/plane at (0,0) if s>0
    ah = P(0, h, lcf / 4)
    a00 = P(0, 0, lcf)
    L = []
    L.append(geo.addLine(a0, a1)); L.append(geo.addLine(a1, a2)); L.append(geo.addLine(a2, a3))
    axis_pts = [a3, ah, a00, a0]
    for i in range(len(axis_pts) - 1): L.append(geo.addLine(axis_pts[i], axis_pts[i + 1]))
    loop = geo.addCurveLoop(L); surf = geo.addPlaneSurface([loop])
    # trace: from (0,h) to (w/2,h); plane: from (s,0) to (X,0) (if s==0 the plane touches the axis)
    te = P(w / 2, h, lcf / 8)
    tr = geo.addLine(ah, te)
    ps = P(s, 0, lcf / 8) if s > 0 else a00
    pe = P(X * 0.999, 0, lcF)
    pl = geo.addLine(ps, pe)
    geo.synchronize()
    gmsh.model.mesh.embed(1, [tr, pl], 2, surf)
    # refinement near the trace edge and slot edge
    f = gmsh.model.mesh.field
    d = f.add('Distance'); f.setNumbers(d, 'PointsList', [te, ps, ah]); f.setNumber(d, 'Sampling', 100)
    t = f.add('Threshold'); f.setNumber(t, 'InField', d); f.setNumber(t, 'SizeMin', lcf / 20 * lc_scale); f.setNumber(t, 'SizeMax', lcF)
    f.setNumber(t, 'DistMin', lcf / 4); f.setNumber(t, 'DistMax', 30 * max(h, w, s, 1))
    d2 = f.add('Distance'); f.setNumbers(d2, 'CurvesList', [tr, pl]); f.setNumber(d2, 'Sampling', 400)
    t2 = f.add('Threshold'); f.setNumber(t2, 'InField', d2); f.setNumber(t2, 'SizeMin', lcf * lc_scale); f.setNumber(t2, 'SizeMax', lcF); f.setNumber(t2, 'DistMin', lcf); f.setNumber(t2, 'DistMax', 20 * max(h, w, s, 1))
    mn = f.add('Min'); f.setNumbers(mn, 'FieldsList', [t, t2]); f.setAsBackgroundMesh(mn)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0); gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0); gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.model.mesh.generate(2)
    nt, xyz, _ = gmsh.model.mesh.getNodes(); xy = xyz.reshape(-1, 3)[:, :2]
    et, etags, enodes = gmsh.model.mesh.getElements(2)
    tri = np.array(enodes[0], dtype=np.int64).reshape(-1, 3)
    idx = {int(k): i for i, k in enumerate(nt)}
    tri = np.vectorize(idx.get)(tri)
    gmsh.finalize()
    m = MeshTri(xy.T.copy(), tri.T.copy())
    basis = Basis(m, ElementTriP2())
    @BilinearForm
    def lap(u, v, _): return dot(grad(u), grad(v))
    K = lap.assemble(basis)
    xd = basis.doflocs
    tolr = 1e-6 * max(h, w)
    on_tr = (np.abs(xd[1] - h) < tolr) & (xd[0] <= w / 2 + tolr)
    on_pl = (np.abs(xd[1]) < tolr) & (xd[0] >= s - tolr)
    on_out = (np.abs(xd[1] + Hc) < tolr) | (np.abs(xd[0] - X) < tolr) | (np.abs(xd[1] - Yt) < tolr)
    V = np.zeros(basis.N); V[on_tr] = 1.0
    D = np.where(on_tr | on_pl | on_out)[0]
    V = solve(*condense(K, x=V, D=D))
    Wn = V @ (K @ V)  # = C/eps0 for the half domain (energy 1/2 C V^2, with integral |grad V|^2)
    C = 2 * EPS0 * Wn  # both halves, F/m (units cancel: 2D)
    return MU0 * EPS0 / C, int(basis.N)
if __name__ == '__main__':
    w, h, Hc = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    top = None
    for s in [0.0] + [float(x) for x in sys.argv[4:]]:
        Lp, n = lprime(w * 1e-3, h * 1e-3, s * 1e-3, Hc * 1e-3)
        print(f's {s:.3f} mm  L\' {Lp*1e6:.4f} pH/um  dofs {n}')
