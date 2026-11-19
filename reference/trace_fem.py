#!/usr/bin/env python3
"""
Independent check of the 2D trace solver (src/physics/fieldSolver2d.ts).

Solves the same electrostatic problem with a different method and different code: second-order
(P2) finite elements on an unstructured triangle mesh from Gmsh, refined towards the trace corners,
full cross-section (no symmetry assumed), assembled with scikit-fem. The numbers it prints are
stored in src/physics/reference.ts as the expected values for the thick-trace cases.

  pip install scikit-fem scipy numpy      (plus Gmsh's Python module)
  python3 scripts/reference/trace_fem.py

Each case is run on two meshes; the difference between them is printed as the uncertainty.
"""
import json, sys
import numpy as np
import gmsh
from skfem import MeshTri, Basis, ElementTriP2, ElementTriP0, BilinearForm, condense, solve
from skfem.helpers import dot, grad

C0 = 299792458.0
EPS0 = 8.8541878128e-12

CASES = {
    # name: slabs [(h, er)] from the lower plane up, topPlane, y0, w, t, s (None = single), mask (t, er[, tBase]) or None
    "ms":        dict(slabs=[(0.15e-3, 4.2)], top=False, y0=0.15e-3, w=0.3e-3, t=35e-6, s=None, mask=None),
    "ms_mask":   dict(slabs=[(0.15e-3, 4.2)], top=False, y0=0.15e-3, w=0.3e-3, t=35e-6, s=None, mask=(25e-6, 3.5)),
    "ms_diff":   dict(slabs=[(0.1e-3, 4.2)], top=False, y0=0.1e-3, w=0.15e-3, t=35e-6, s=0.15e-3, mask=None),
    "ms_diff_mask": dict(slabs=[(0.1e-3, 4.2)], top=False, y0=0.1e-3, w=0.15e-3, t=35e-6, s=0.15e-3, mask=(20e-6, 3.5)),
    "sl_offset": dict(slabs=[(0.1e-3, 4.4), (0.218e-3, 3.9)], top=True, y0=0.1e-3, w=0.1e-3, t=18e-6, s=None, mask=None),
    "sl_diff":   dict(slabs=[(0.1e-3, 4.4), (0.218e-3, 3.9)], top=True, y0=0.1e-3, w=0.1e-3, t=18e-6, s=0.12e-3, mask=None),
    "ms_embedded": dict(slabs=[(0.2e-3, 4.3), (0.035e-3, 4.3), (0.2e-3, 4.3)], top=False, y0=0.2e-3, w=0.3e-3, t=35e-6, s=None, mask=None),
    # Etched (trapezoidal) traces: `narrow` is the total width lost on the narrow face.
    "ms_trap":      dict(slabs=[(0.15e-3, 4.2)], top=False, y0=0.15e-3, w=0.3e-3, t=35e-6, s=None, mask=None, narrow=30e-6),
    "ms_trap_mask": dict(slabs=[(0.15e-3, 4.2)], top=False, y0=0.15e-3, w=0.3e-3, t=35e-6, s=None, mask=(25e-6, 3.5), narrow=30e-6),
    "ms_diff_trap": dict(slabs=[(0.1e-3, 4.2)], top=False, y0=0.1e-3, w=0.15e-3, t=35e-6, s=0.15e-3, mask=(20e-6, 3.5), narrow=35e-6),
    "sl_trap_down": dict(slabs=[(0.1e-3, 4.4), (0.218e-3, 3.9)], top=True, y0=0.1e-3, w=0.1e-3, t=18e-6, s=0.12e-3, mask=None, narrow=18e-6, narrowAt="bottom"),
    # Mask given as (on copper, er, on laminate): a fabricator's usual thin-on-copper, thick-on-laminate coat.
    "ms_diff_fab":  dict(slabs=[(0.2e-3, 4.3)], top=False, y0=0.2e-3, w=0.2e-3, t=40.64e-6, s=0.2e-3, mask=(15.24e-6, 3.8, 30.48e-6), narrow=17.78e-6),
    # Coplanar ground on the trace's layer, `gap` from the trace edge.
    "cpwg":         dict(slabs=[(0.2e-3, 4.3)], top=False, y0=0.2e-3, w=0.3e-3, t=35e-6, s=None, mask=None, gap=0.15e-3),
    "cpwg_mask_trap": dict(slabs=[(0.2e-3, 4.3)], top=False, y0=0.2e-3, w=0.3e-3, t=35e-6, s=None, mask=(25e-6, 3.5), gap=0.15e-3, narrow=30e-6),
    "cpwg_diff":    dict(slabs=[(0.1e-3, 4.2)], top=False, y0=0.1e-3, w=0.15e-3, t=35e-6, s=0.15e-3, mask=None, gap=0.2e-3),
    "cpw_inner":    dict(slabs=[(0.1e-3, 4.4), (0.218e-3, 3.9)], top=True, y0=0.1e-3, w=0.1e-3, t=18e-6, s=None, mask=None, gap=0.1e-3),
}


def polygons(case, U=1e6):
    """Conductor and mask outlines in micrometres. Returns traces, grounds, mask pieces, box."""
    slabs = [(h * U, er) for h, er in case["slabs"]]
    top, y0, w, t = case["top"], case["y0"] * U, case["w"] * U, case["t"] * U
    s = case["s"] * U if case["s"] else None
    mt = case["mask"][0] * U if case["mask"] else 0.0            # mask on the copper
    mb = case["mask"][2] * U if case["mask"] and len(case["mask"]) > 2 else mt  # mask on the laminate
    d = case.get("narrow", 0.0) * U / 2            # lean of each side wall
    wide_low = case.get("narrowAt", "top") == "top"
    gap = case.get("gap", 0.0) * U
    H = sum(h for h, _ in slabs)
    spans = [(-w / 2, w / 2)] if s is None else [(-s / 2 - w, -s / 2), (s / 2, s / 2 + w)]
    span = max(abs(x) for tr in spans for x in tr)
    X = span + (4 * H if top else 200 * (y0 + t + span))
    Y = H if top else y0 + t + 200 * (y0 + t + span)
    yt = y0 + t
    lo, hi = (0.0, d) if wide_low else (d, 0.0)      # inset of the lower and upper face
    traces = [[(a + lo, y0), (b - lo, y0), (b - hi, yt), (a + hi, yt)] for a, b in spans]
    grounds, masks = [], []
    if gap:
        xg = span + gap
        grounds = [[(xg + lo, y0), (X, y0), (X, yt), (xg + hi, yt)], [(-X, y0), (-xg - lo, y0), (-xg - hi, yt), (-X, yt)]]
    if mt:
        masks.append([(-X, H), (X, H), (X, H + mb), (-X, H + mb)])
        for a, b in spans:  # wall parallel to the copper up to its top, then straight up
            masks.append([(a + lo - mt, y0), (b - lo + mt, y0), (b - hi + mt, yt), (b - hi + mt, yt + mt), (a + hi - mt, yt + mt), (a + hi - mt, yt)])
        if gap:
            masks.append([(xg + lo - mt, y0), (X, y0), (X, yt + mt), (xg + hi - mt, yt + mt), (xg + hi - mt, yt)])
            masks.append([(-X, y0), (-xg - lo + mt, y0), (-xg - hi + mt, yt), (-xg - hi + mt, yt + mt), (-X, yt + mt)])
    return slabs, top, traces, grounds, masks, (X, Y, H), (y0, t, max(mt, H + mb - y0 - t))


def inside(pt, poly, tol):
    """Point strictly inside a convex or simple polygon (ray casting), with a bounding-box margin."""
    x, y = pt
    n, c = len(poly), False
    for i in range(n):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            c = not c
    return c


def on_outline(mid, poly, tol):
    """Boolean array: which points (2 x N) lie on the outline of the polygon."""
    out = np.zeros(mid.shape[1], dtype=bool)
    for i in range(len(poly)):
        p1, p2 = np.array(poly[i]), np.array(poly[(i + 1) % len(poly)])
        e = p2 - p1
        L = np.hypot(*e)
        rel = mid - p1[:, None]
        along = (rel[0] * e[0] + rel[1] * e[1]) / L
        across = np.abs(rel[0] * e[1] - rel[1] * e[0]) / L
        out |= (across < tol) & (along > -tol) & (along < L + tol)
    return out


def build(case, fine):
    """Mesh the full cross-section. Returns points, triangles, er per triangle and the outlines."""
    slabs, top, traces, grounds, masks, (X, Y, H), (y0, t, mt) = polygons(case)
    tol = 1e-3
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ

    def surface(poly):
        pts = [occ.addPoint(x, y, 0) for x, y in poly]
        lines = [occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
        return occ.addPlaneSurface([occ.addCurveLoop(lines)])

    tools, y = [], 0.0
    for h, _ in slabs:
        tools.append((2, occ.addRectangle(-X, y, 0, 2 * X, h)))
        y += h
    if not top:
        tools.append((2, occ.addRectangle(-X, H, 0, 2 * X, Y - H)))
    tools += [(2, surface(m)) for m in masks]
    cond = [(2, surface(c)) for c in traces + grounds]
    occ.synchronize()
    occ.fragment(tools, cond)
    occ.synchronize()

    def bbox_in(tag, poly):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(2, tag)
        px, py = [q[0] for q in poly], [q[1] for q in poly]
        return xmin >= min(px) - tol and xmax <= max(px) + tol and ymin >= min(py) - tol and ymax <= max(py) + tol

    surf = []  # (tag, er)
    y_tops = np.cumsum([h for h, _ in slabs])
    for _, tag in gmsh.model.getEntities(2):
        cx, cy, _ = occ.getCenterOfMass(2, tag)
        if any(bbox_in(tag, c) and inside((cx, cy), c, tol) for c in traces + grounds):
            continue  # copper: not part of the field region
        er = 1.0
        for (h, e), ytop in zip(slabs, y_tops):
            if cy < ytop:
                er = e
                break
        else:
            # Above the laminate, the mask is the only piece that stays below its own top surface.
            if mt and gmsh.model.getBoundingBox(2, tag)[4] <= y0 + t + mt + tol:
                er = case["mask"][1]
        surf.append((tag, er))

    corners = [q for c in traces + grounds for q in c if abs(abs(q[0]) - X) > tol]
    pts = []
    for _, tag in gmsh.model.getEntities(0):
        px, py, _ = gmsh.model.getValue(0, tag, [])
        if any(abs(px - qx) < tol and abs(py - qy) < tol for qx, qy in corners):
            pts.append(tag)
    f = gmsh.model.mesh.field
    dist = f.add("Distance"); f.setNumbers(dist, "PointsList", pts)
    feat = min(abs(c[i][j] - c[(i + 1) % len(c)][j]) for c in traces for i in range(len(c)) for j in (0, 1) if abs(c[i][j] - c[(i + 1) % len(c)][j]) > tol)
    feat = min(feat, y0, case.get("gap", 1) * 1e6 if grounds else feat)
    hmin = feat / (400 if fine else 200)
    m = f.add("MathEval"); f.setString(m, "F", f"{hmin} + {0.12 if fine else 0.2}*F{dist}")
    f.setAsBackgroundMesh(m)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    idx = np.zeros(int(node_tags.max()) + 1, dtype=np.int64)
    idx[node_tags.astype(np.int64)] = np.arange(len(node_tags))
    p = coords.reshape(-1, 3)[:, :2].T.copy()
    tris, ers = [], []
    for tag, er in surf:
        types, _, nodes = gmsh.model.mesh.getElements(2, tag)
        for et, nn in zip(types, nodes):
            if et == 2:
                tri = idx[nn.astype(np.int64)].reshape(-1, 3)
                tris.append(tri)
                ers.append(np.full(len(tri), er))
    gmsh.finalize()
    tt = np.vstack(tris).T.copy()
    used = np.unique(tt)
    remap = -np.ones(p.shape[1], dtype=np.int64); remap[used] = np.arange(len(used))
    return p[:, used], remap[tt], np.concatenate(ers), traces, grounds, (X, Y, H)


def capacitance(case, fine, volts, use_er=True):
    p, tt, er, traces, grounds, (X, Y, H) = build(case, fine)
    mesh = MeshTri(p, tt)
    basis = Basis(mesh, ElementTriP2())
    eps = basis.with_element(ElementTriP0()).interpolate(er if use_er else np.ones_like(er))

    @BilinearForm
    def a(u, v, w):
        return w.eps * dot(grad(u), grad(v))

    A = a.assemble(basis, eps=eps)
    tol = 1e-3
    x = np.zeros(basis.N)
    bf = mesh.boundary_facets()
    mid = mesh.p[:, mesh.facets[:, bf]].mean(axis=1)
    dofs = lambda sel: basis.get_dofs(bf[sel]).all()
    ground = np.abs(mid[1]) < tol
    if case["top"]:
        ground |= np.abs(mid[1] - Y) < tol
    for gnd in grounds:
        ground |= on_outline(mid, gnd, tol)
    D = [dofs(ground)]
    for tr, v in zip(traces, volts):
        d = dofs(on_outline(mid, tr, tol)); x[d] = v; D.append(d)
    D = np.unique(np.concatenate(D))
    x = solve(*condense(A, x=x, D=D))
    return EPS0 * float(x @ (A @ x)), mesh.t.shape[1]


def run(name):
    case = CASES[name]
    pair = case["s"] is not None
    modes = {"single": [1.0]} if not pair else {"odd": [1.0, -1.0], "even": [1.0, 1.0]}
    out = {}
    for mode, volts in modes.items():
        res = []
        for fine in (False, True):
            div = 2.0 if pair else 1.0  # energy of two lines
            c, n = capacitance(case, fine, volts, True)
            ca, _ = capacitance(case, fine, volts, False)
            c /= div; ca /= div
            res.append((1 / (C0 * np.sqrt(c * ca)), c / ca, n))
        (z1, e1, n1), (z2, e2, n2) = res
        out[mode] = dict(z0=z2, eeff=e2, dz=abs(z2 - z1), deeff=abs(e2 - e1), triangles=n2)
        print(f"{name:14s} {mode:6s} Z0 {z2:.4f} Ω (±{abs(z2 - z1):.4f})  εeff {e2:.5f} (±{abs(e2 - e1):.5f})  {n1}/{n2} triangles", flush=True)
    return out


if __name__ == "__main__":
    names = sys.argv[1:] or list(CASES)
    results = {n: run(n) for n in names}
    print(json.dumps(results, indent=1))
