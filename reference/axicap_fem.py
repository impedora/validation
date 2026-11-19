"""
Independent check of Impedora's axisymmetric barrel capacitance (barrelStackCapacitanceEx in
src/physics/axisym.ts) for a blind via that ends on a stripline pad inside a plane pair.

Second-order finite elements (scikit-fem) on a Gmsh mesh in (r, z), refined at every conductor edge.
The barrel and its pad are at 1 V, the planes at 0 V, the plane beyond the barrel's end closes the
bottom of the domain. The charge on each conductor is the nodal reaction of the assembled system.

Geometry (mm, z up, top of the stack at 0): a plane 0 to -tP, a first plane pair of height h1, a
second plane, then the barrel ends with a pad of radius rp at depth d below the second plane, and
the solid plane closes the pair at depth h2 below the second plane.

Usage:  python3 axicap_fem.py <rb> <rp> <ra> <tP> <h1> <h2> <d> <er> [pad: 1|0]
Example (case bl1):  python3 axicap_fem.py 0.125 0.25 0.35 0.035 0.3 0.6 0.2 4.2
Needs gmsh and scikit-fem (pip install gmsh scikit-fem).
"""
import sys, json
import numpy as np
import gmsh
from skfem import MeshTri, Basis, ElementTriP2, BilinearForm, asm, condense, solve
from skfem.helpers import dot, grad

EPS0 = 8.8541878128e-12


def capacitance(rb, rp, ra, tP, h1, h2, d, er, pad=True, R=6.0, lc_scale=1.0):
    gmsh.initialize(); gmsh.option.setNumber('General.Terminal', 0); gmsh.model.add('axi')
    occ = gmsh.model.occ
    zP1a, zP1b = -(tP + h1), -(2 * tP + h1)          # second plane: top and bottom face
    zPad = zP1b - d                                   # pad level (barrel end)
    zBot = zP1b - h2                                  # top face of the closing plane
    dom = occ.addRectangle(0, zBot, 0, R, -zBot)
    barrel = occ.addRectangle(0, zPad, 0, rb, -zPad)
    p0 = occ.addRectangle(ra, -tP, 0, R - ra, tP)
    p1 = occ.addRectangle(ra, zP1b, 0, R - ra, tP)
    cut, _ = occ.cut([(2, dom)], [(2, barrel), (2, p0), (2, p1)])
    tools = []
    if pad:
        a = occ.addPoint(rb, zPad, 0); b = occ.addPoint(rp, zPad, 0)
        tools.append((1, occ.addLine(a, b)))
    if tools:
        cut, _ = occ.fragment(cut, tools)
    occ.synchronize()
    edges = [(rb, zPad), (rp, zPad), (ra, 0.0), (ra, -tP), (ra, zP1a), (ra, zP1b)]
    pts = [occ.addPoint(x, y, 0) for x, y in edges]
    occ.synchronize()
    lmin = min(rb, rp - rb if pad else rb, ra - rp, tP, d, h2 - d) / 12 * lc_scale
    f = gmsh.model.mesh.field
    dist = f.add('Distance'); f.setNumbers(dist, 'PointsList', pts)
    th = f.add('Threshold'); f.setNumber(th, 'InField', dist); f.setNumber(th, 'SizeMin', lmin / 8)
    f.setNumber(th, 'SizeMax', R / 20); f.setNumber(th, 'DistMin', lmin / 4); f.setNumber(th, 'DistMax', 4 * (h1 + h2))
    th2 = f.add('MathEval'); f.setString(th2, 'F', f'{lmin} + 0.25*sqrt(x*x)')
    mn = f.add('Min'); f.setNumbers(mn, 'FieldsList', [th, th2]); f.setAsBackgroundMesh(mn)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0); gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.model.mesh.generate(2)
    _, xyz, _ = gmsh.model.mesh.getNodes(); xy = xyz.reshape(-1, 3)[:, :2]
    _, _, enodes = gmsh.model.mesh.getElements(2)
    tri = np.array(enodes[0], dtype=np.int64).reshape(-1, 3) - 1
    gmsh.finalize()
    used = np.unique(tri); remap = -np.ones(len(xy), dtype=np.int64); remap[used] = np.arange(len(used))
    m = MeshTri(xy[used].T.copy(), remap[tri].T.copy())
    basis = Basis(m, ElementTriP2())

    @BilinearForm
    def lap(u, v, w):
        return dot(grad(u), grad(v)) * w.x[0]
    K = asm(lap, basis)
    X = basis.doflocs
    tol = 1e-9
    on = lambda c: np.abs(c) < tol
    barrel_d = (X[0] <= rb + tol) & (X[1] >= zPad - tol)
    pad_d = on(X[1] - zPad) & (X[0] <= rp + tol) if pad else np.zeros(X.shape[1], bool)
    p0_d = (X[0] >= ra - tol) & (X[1] >= -tP - tol)
    p1_d = (X[0] >= ra - tol) & (X[1] <= zP1a + tol) & (X[1] >= zP1b - tol)
    bot_d = on(X[1] - zBot)
    hot = barrel_d | pad_d
    groups = {'P0': p0_d & ~hot, 'P1': p1_d & ~hot, 'end': bot_d & ~hot}
    D = hot | groups['P0'] | groups['P1'] | groups['end']
    u = np.zeros(X.shape[1]); u[hot] = 1.0
    u = solve(*condense(K, x=u, D=np.where(D)[0]))
    q = K @ u
    scale = 2 * np.pi * EPS0 * er * 1e-3            # lengths in mm
    out = {k: float(-q[v].sum() * scale) for k, v in groups.items()}
    out['total'] = float(u @ (K @ u) * scale)
    # P1 split into its upper and lower face, as Impedora books it (by the mid-thickness).
    zmid = (zP1a + zP1b) / 2
    out['P1 up'] = float(-q[groups['P1'] & (X[1] >= zmid)].sum() * scale)
    out['P1 dn'] = float(-q[groups['P1'] & (X[1] < zmid)].sum() * scale)
    out['dofs'] = int(X.shape[1])
    return out


if __name__ == '__main__':
    if len(sys.argv) < 9:
        sys.exit(__doc__)
    a = [float(x) for x in sys.argv[1:9]]
    pad = len(sys.argv) < 10 or sys.argv[9] != '0'
    for s in (1.0, 0.6):
        r = capacitance(*a, pad=pad, lc_scale=s)
        print(json.dumps({k: (round(v * 1e15, 3) if isinstance(v, float) else v) for k, v in r.items()}), '(fF)')
