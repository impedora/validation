"""
Independent full-wave check of round (axisymmetric) via structures: a body-of-revolution finite-element
solve, shares no code or method with Impedora's via model.

The field of a round via between round planes has no azimuthal variation: only E_r, E_z and H_phi
(TM0). The unknown is H_phi(r, z) on second-order triangles (Gmsh), with the r-weighted weak form

    ∫ (1/εr) [ (∂H/∂r + H/r)(∂v/∂r + v/r) + ∂H/∂z ∂v/∂z ] r dA − k0² ∫ H v r dA = −jωε0 ∮ (n × E)_φ v r dl

Metal is perfect (natural boundary), the axis has H = 0, a coaxial TEM port in an antipad annulus has
H = I/(2πr), and the outer wall r = R is the same first-order absorbing wall as in the 3D cases
(E_z = −η H_phi). Each port is driven with the others open (H = 0 on their aperture); the port voltage
comes from the weak-form residual on the port with the test function 1/r. Output: the Z matrix.

What it can check: anything round, so a via without stitching vias, its pads, a barrel that ends
inside a plane pair. What it cannot: stitching vias, traces, pairs.

    python3 via_bor_fem.py check     the two documented checks (docs/VALIDATION.md), about a minute
    python3 via_bor_fem.py <geometry.json> <out.json>

Geometry (mm): {"R", "er", "planes": [{"z0", "z1", "ra" or null}], "barrel": {"a", "z0", "z1"},
"pads": [{"z", "r", "t"}], "ports": ["top" | "bottom"], "near", "far", "fmax", "fstep", "refine": [[r, z]]}
with z up and the ports on the outer faces of the outermost planes.
Needs gmsh, numpy and scipy (pip install gmsh numpy scipy).
"""
import json
import math
import sys

import gmsh
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

EPS0, MU0, C0 = 8.8541878128e-12, 1.25663706212e-6, 299792458.0
MM = 1e-3


def solve(g):
    R, er = g["R"], g["er"]
    zmin = min(p["z0"] for p in g["planes"])
    zmax = max(p["z1"] for p in g["planes"])
    ba = g["barrel"]
    a = ba["a"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("bor")
    occ = gmsh.model.occ
    dom = occ.addRectangle(0, zmin, 0, R, zmax - zmin)
    cuts = [(2, occ.addRectangle(p["ra"] or 0.0, p["z0"], 0, R - (p["ra"] or 0.0), p["z1"] - p["z0"])) for p in g["planes"]]
    cuts.append((2, occ.addRectangle(0, ba["z0"], 0, a, ba["z1"] - ba["z0"])))
    for pd in g.get("pads", []):  # thin disk, overlapping the barrel a little
        cuts.append((2, occ.addRectangle(a - 1e-4, pd["z"] - pd["t"] / 2, 0, pd["r"] - a + 1e-4, pd["t"])))
    occ.cut([(2, dom)], cuts)
    occ.synchronize()
    # Mesh: fine near the given corners (sheet and antipad edges, where the field is singular).
    pids = [tag for (_, tag) in gmsh.model.getEntities(0)
            if any(abs(gmsh.model.getValue(0, tag, [])[0] - pr) < 0.004 and abs(gmsh.model.getValue(0, tag, [])[1] - pz) < 0.004 for pr, pz in g.get("refine", []))]
    f1 = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f1, "PointsList", pids)
    f2 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f2, "InField", f1)
    gmsh.model.mesh.field.setNumber(f2, "SizeMin", g["near"])
    gmsh.model.mesh.field.setNumber(f2, "SizeMax", g["far"])
    gmsh.model.mesh.field.setNumber(f2, "DistMin", 2 * g["near"])
    gmsh.model.mesh.field.setNumber(f2, "DistMax", 1.5)
    gmsh.model.mesh.field.setAsBackgroundMesh(f2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.model.mesh.generate(2)
    nt, xyz, _ = gmsh.model.mesh.getNodes()
    X = xyz.reshape(-1, 3)[:, :2].copy()
    idx = {int(n): i for i, n in enumerate(nt)}
    et, _, en = gmsh.model.mesh.getElements(2)
    tri = np.array([idx[int(n)] for n in en[list(et).index(9)]]).reshape(-1, 6)
    bt, _, bn = gmsh.model.mesh.getElements(1)
    lines = np.array([idx[int(n)] for n in bn[list(bt).index(8)]]).reshape(-1, 3)
    gmsh.finalize()
    used = np.unique(tri.ravel())
    remap = -np.ones(len(X), dtype=int)
    remap[used] = np.arange(len(used))
    X, tri, lines = X[used], remap[tri], remap[lines]
    lines = lines[np.all(lines >= 0, axis=1)]
    N = len(X)

    # Element matrices, vectorised over the elements (6-point quadrature, exact for P2 here up to r).
    qp = np.array([[0.445948490915965, 0.445948490915965], [0.445948490915965, 0.108103018168070], [0.108103018168070, 0.445948490915965],
                   [0.091576213509771, 0.091576213509771], [0.091576213509771, 0.816847572980459], [0.816847572980459, 0.091576213509771]])
    qw = np.array([0.223381589678011] * 3 + [0.109951743655322] * 3) / 2
    P = X[tri]
    J00, J01 = P[:, 1, 0] - P[:, 0, 0], P[:, 2, 0] - P[:, 0, 0]
    J10, J11 = P[:, 1, 1] - P[:, 0, 1], P[:, 2, 1] - P[:, 0, 1]
    det = J00 * J11 - J01 * J10
    iJ = np.stack([np.stack([J11, -J01], -1), np.stack([-J10, J00], -1)], -2) / det[:, None, None]
    Ke = np.zeros((len(tri), 6, 6))
    Me = np.zeros((len(tri), 6, 6))
    for (x, y), w in zip(qp, qw):
        l1, l2, l3 = 1 - x - y, x, y
        Nn = np.array([l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), l3 * (2 * l3 - 1), 4 * l1 * l2, 4 * l2 * l3, 4 * l3 * l1])
        G = np.stack([[-(4 * l1 - 1), 4 * l2 - 1, 0, 4 * (l1 - l2), 4 * l3, -4 * l3],
                      [-(4 * l1 - 1), 0, 4 * l3 - 1, -4 * l2, 4 * l2, 4 * (l1 - l3)]], -1)
        gr = np.einsum("kd,edc->ekc", G, iJ)
        r = (Nn[None, :] * P[:, :, 0]).sum(1)
        A = gr[:, :, 0] + Nn[None, :] / r[:, None]
        wt = w * np.abs(det) * r
        Ke += (A[:, :, None] * A[:, None, :] + gr[:, :, 1][:, :, None] * gr[:, :, 1][:, None, :]) * (wt / er)[:, None, None]
        Me += (Nn[:, None] * Nn[None, :])[None] * wt[:, None, None]
    rows, cols = np.repeat(tri, 6, axis=1).ravel(), np.tile(tri, (1, 6)).ravel()
    K = sp.csr_matrix((Ke.ravel(), (rows, cols)), shape=(N, N)) * MM          # lengths in mm → SI
    M = sp.csr_matrix((Me.ravel(), (rows, cols)), shape=(N, N)) * MM ** 3

    tol = 1e-7
    gl = [(-0.774596669241483, 5 / 9), (0.0, 8 / 9), (0.774596669241483, 5 / 9)]
    brow, bcol, bv = [], [], []
    port_nodes = {k: set() for k in range(len(g["ports"]))}
    axis = set()
    ptop = next(p for p in g["planes"] if p["z1"] == zmax)
    pbot = next(p for p in g["planes"] if p["z0"] == zmin)
    for Ln in lines:
        Q = X[Ln]
        if np.all(abs(Q[:, 0] - R) < tol):  # absorbing wall: ∫ φ φ R dz
            Le = abs(Q[1, 1] - Q[0, 1])
            Be = sum(np.outer(n, n) * R * w * Le / 2 for s, w in gl for n in [np.array([s * (s - 1) / 2, s * (s + 1) / 2, 1 - s * s])])
            brow += list(np.repeat(Ln, 3)); bcol += list(np.tile(Ln, 3)); bv += list(Be.ravel())
        for k, side in enumerate(g["ports"]):
            zf, pl = (zmax, ptop) if side == "top" else (zmin, pbot)
            if pl["ra"] is not None and np.all(abs(Q[:, 1] - zf) < tol) and np.all(Q[:, 0] >= a - tol) and np.all(Q[:, 0] <= pl["ra"] + tol):
                port_nodes[k] |= set(Ln.tolist())
        if np.all(abs(Q[:, 0]) < tol):
            axis |= set(Ln.tolist())
    B = sp.csr_matrix((bv, (brow, bcol)), shape=(N, N)) * MM ** 2
    ports = [np.array(sorted(port_nodes[k]), dtype=int) for k in range(len(g["ports"]))]
    sgn = [1 if s == "top" else -1 for s in g["ports"]]  # port normal ±z
    axis = np.array(sorted(axis), dtype=int)
    eta = math.sqrt(MU0 / (EPS0 * er))
    dirich = np.unique(np.concatenate([*ports, axis])).astype(int)
    free = np.setdiff1d(np.arange(N), dirich)
    out = {"f": [], "Z": [], "nodes": int(N)}
    for fg in np.arange(g["fstep"], g["fmax"] + 1e-9, g["fstep"]):
        w = 2 * math.pi * fg * 1e9
        k0 = w / C0
        A = ((K - k0 * k0 * M).astype(complex) + 1j * w * EPS0 * eta * B).tocsr()
        lu = spl.splu(A[free][:, free].tocsc())
        Z = [[0j] * len(ports) for _ in ports]
        for k in range(len(ports)):
            U = np.zeros(N, dtype=complex)
            U[ports[k]] = -sgn[k] / (X[ports[k], 0] * MM)  # 2π A into the structure
            sol = U.copy()
            sol[free] = lu.solve(-(A[free][:, dirich] @ U[dirich]))
            res = A @ sol
            for j in range(len(ports)):
                V = (res[ports[j]] @ (1.0 / (X[ports[j], 0] * MM))) / (-1j * w * EPS0 * sgn[j])
                Z[j][k] = V / (2 * math.pi)
        out["f"].append(float(fg))
        out["Z"].append([[[z.real, z.imag] for z in row] for row in Z])
    return out


def blind(h, d, a, ra, rp, R=3.0, er=4.2, t=0.035):
    """A barrel from the top face of the upper plane down to a pad d below it, in a plane pair of height h; solid plane below."""
    zp = h - d
    return {"R": R, "er": er, "planes": [{"z0": -t, "z1": 0.0, "ra": None}, {"z0": h, "z1": h + t, "ra": ra}],
            "barrel": {"a": a, "z0": zp, "z1": h + t}, "pads": [{"z": zp, "r": rp, "t": 0.002}], "ports": ["top"],
            "near": 0.003, "far": 0.08, "fmax": 50, "fstep": 1,
            "refine": [[a, zp], [rp, zp], [ra, h], [ra, h + t], [a, h + t], [0.0, zp]]}


def check():
    # 1. Through via in one plane pair, no stitching vias (3D case p2g0): S against Palace.
    p2g0 = {"R": 3.0, "er": 4.2, "planes": [{"z0": 0.0, "z1": 0.035, "ra": 0.35}, {"z0": 0.835, "z1": 0.87, "ra": 0.35}],
            "barrel": {"a": 0.125, "z0": 0.0, "z1": 0.87}, "pads": [], "ports": ["top", "bottom"],
            "near": 0.004, "far": 0.08, "fmax": 50, "fstep": 1,
            "refine": [[0.35, 0.035], [0.35, 0.835], [0.35, 0.0], [0.35, 0.87], [0.125, 0.0], [0.125, 0.87]]}
    r = solve(p2g0)
    print("p2g0, through via, 50 ohm ports: f GHz, S11, S21")
    for i, f in enumerate(r["f"]):
        if int(f) in (1, 10, 20, 30, 40, 50):
            Z = np.array([[complex(*r["Z"][i][j][k]) for k in range(2)] for j in range(2)])
            S = (Z - 50 * np.eye(2)) @ np.linalg.inv(Z + 50 * np.eye(2))
            print(f"{f:4.0f}  {S[0, 0].real:+.4f}{S[0, 0].imag:+.4f}j  {S[1, 0].real:+.4f}{S[1, 0].imag:+.4f}j")
    # 2. Blind via end (bl2's barrel and pad, no trace, no stitching vias): input impedance at the top.
    r = solve(blind(0.4, 0.15, 0.1, 0.3, 0.225))
    print("blind end, h 0.4 mm, pad 0.15 mm deep: f GHz, Z ohm, C seen = -1/(ωX) fF")
    for i, f in enumerate(r["f"]):
        if int(f) in (1, 5, 10, 15, 20, 25, 30, 40, 50):
            z = complex(*r["Z"][i][0][0])
            print(f"{f:4.0f}  {z.real:8.4f} {z.imag:+10.3f}j  {-1 / (2 * math.pi * f * 1e9 * z.imag) * 1e15:8.2f}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "check":
        check()
    else:
        json.dump(solve(json.load(open(sys.argv[1]))), open(sys.argv[2], "w"))
