#!/usr/bin/env python3
"""
A plain trace between two trace ports, no via: the "thru" that calibrates the trace ports of the
via cases with a launch or an exit (cases.json, "launch" and "exit").

A lumped port sheet at the end of a trace is not a perfect transmission-line port: the current has
to turn from the trace into the sheet, which adds a little series inductance and shunt capacitance
of its own. The via cases with a trace have one such port. Solving the same trace with a sheet port
at both ends, built exactly the same way, and comparing it with the ideal line from the 2D trace
solver measures that parasitic, so it can be taken out of the via comparison (src/physics/fullwave.ts).

The geometry copies via_case.py: same dielectric, same trace width and height, same port sheets
(microstrip: one sheet from the trace down to the plane; stripline: one sheet up and one down, in
parallel), same mesh sizes near the trace. The domain is a box with absorbing side walls; the trace
runs along x from −len to +len.

For a differential pair the thru is the two traces side by side at the via pitch (centre to
centre), with ports 1 and 2 at one end and 3 and 4 at the other; the comparison uses its odd mode.

A second thru at half the length ("half") separates what belongs to the ports from what belongs to
the trace: the 3D model's trace is not exactly the 2D solver's (mesh error, a few percent in
impedance on the coarse mesh), and with two lengths the calibration measures it instead of
assuming it.

Usage:  python3 trace_thru.py <case-id> <level> <out-dir> [half]     case = a cases.json entry with launch or exit
"""
import json
import os
import sys

import gmsh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from via_case import LEVELS, load_case, port_fence, sheet_edges, size_field, strip_split  # noqa: E402


def build(c, level, path, half=False):
    lv = LEVELS[level]
    ms = c.get("launch")
    ex = c.get("exit")
    if bool(ms) == bool(ex):
        sys.exit("trace_thru: the case needs exactly one of launch or exit.")
    er = c["cavities"][0]["er"]
    w = (ms or ex)["w"]
    ln = (ms or ex)["len"] * (0.5 if half else 1.0)
    if ms:
        h = ms["h"]
        z_tr, z_top = h, h + ms["air"]
        H = h
    else:
        H = c["cavities"][ex.get("pair", len(c["cavities"]) - 1)]["h"]  # the plane pair the stripline is in
        z_tr = H - ex["depth"]
        z_top = H
    X = ln + 1.0  # the box ends 1 mm beyond each port
    centers = [-c["diff"]["pitch"] / 2, c["diff"]["pitch"] / 2] if "diff" in c else [0.0]  # trace centre lines (y)
    Y = 1.5 + 3 * w + max(centers)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("thru")
    occ = gmsh.model.occ
    diel = occ.addBox(-X, -Y, 0, 2 * X, 2 * Y, H)
    vols = [(3, diel)]
    if ms:
        vols.append((3, occ.addBox(-X, -Y, H, 2 * X, 2 * Y, z_top - H)))
    traces = [occ.addRectangle(-ln, yc - w / 2, z_tr, 2 * ln, w) for yc in centers]
    if ex:
        # The same fence vias as at the via case's stripline port, at both ends (PEC, plane to plane).
        rb = c["drill"] / 2
        posts = [(3, occ.addCylinder(x, y, -0.01, 0, 0, H + 0.02, rb)) for xp in (-ln, ln) for (x, y) in port_fence(c, xp)]
        vols, _ = occ.cut(vols, posts)

    def wall(x, za, zb, ya=-w / 2, yb=w / 2):
        pp = [occ.addPoint(x, ya, za), occ.addPoint(x, yb, za), occ.addPoint(x, yb, zb), occ.addPoint(x, ya, zb)]
        return occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(pp[i], pp[(i + 1) % 4]) for i in range(4)])])

    # Port sheets at both ends. Microstrip: one sheet down to the plane. Stripline: strips up and
    # down in the ratio of the return currents, exactly as in via_case.py (strip_split).
    k_up, k_dn = strip_split(ex["depth"], H) if ex else (0, 1)
    cuts = lambda k, yc: [(yc - w / 2 + i * w / k, yc - w / 2 + (i + 1) * w / k) for i in range(k)]
    # One port per trace end, in the order: near end of each trace, then far end of each trace.
    sheets = [[wall(x, z_tr, H, ya, yb) for ya, yb in cuts(k_up, yc)] + [wall(x, 0.0, z_tr, ya, yb) for ya, yb in cuts(k_dn, yc)] for x in (-ln, ln) for yc in centers]
    vol_out, parts = occ.fragment(vols, [(2, t_) for t_ in traces] + [(2, s) for end in sheets for s in end])
    occ.synchronize()
    n_tr = len(traces)
    tr_parts = sorted({abs(t) for d, t in [e for q in parts[len(vols):len(vols) + n_tr] for e in q] if d == 2})
    flat = [e for q in parts[len(vols) + n_tr:] for e in q]
    sheet_parts, k = [], 0
    for end in sheets:
        row = []
        for _ in end:
            row.append(sorted({abs(t) for d, t in parts[len(vols) + n_tr + k] if d == 2}))
            k += 1
        sheet_parts.append(row)
    vtags = [t for d, t in vol_out if d == 3]
    if ms:
        air = [v for v in vtags if gmsh.model.getBoundingBox(3, v)[2] > H - 1e-6]
        dvol = [v for v in vtags if v not in air]
    else:
        air, dvol = [], vtags

    eps = 1e-6
    pec, absorbing = list(tr_parts), []
    for _, s in gmsh.model.getBoundary([(3, v) for v in vtags], combined=True, oriented=False):
        s = abs(s)
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, s)
        if abs(z1) < eps or (not ms and abs(z0 - H) < eps):
            pec.append(s)  # the plane(s)
        elif ex and (x1 - x0) < 2.2 * c["drill"] / 2 and (y1 - y0) < 2.2 * c["drill"] / 2:
            pec.append(s)  # a fence via
        else:
            absorbing.append(s)  # side walls, and the top of the air
    gmsh.model.addPhysicalGroup(3, dvol, 1, "dielectric")
    if air:
        gmsh.model.addPhysicalGroup(3, air, 2, "air")
    gmsh.model.addPhysicalGroup(2, pec, 2, "pec")
    gmsh.model.addPhysicalGroup(2, absorbing, 3, "absorbing")
    # One attribute per element: port 1 from 10, port 2 from 30, … (upper strips first).
    for end in range(len(sheet_parts)):
        for i, tags in enumerate(sheet_parts[end]):
            gmsh.model.addPhysicalGroup(2, tags, 10 + 20 * end + i, f"port{end + 1}-{i + 1}")

    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", tr_parts + [t for end in sheet_parts for row in end for t in row])
    f_th = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_th, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_th, "SizeMin", lv["near"])
    gmsh.model.mesh.field.setNumber(f_th, "SizeMax", lv["far"])
    gmsh.model.mesh.field.setNumber(f_th, "DistMin", 0.5 * c["antipad"] / 2)
    gmsh.model.mesh.field.setNumber(f_th, "DistMax", 0.5 * c.get("domainR", 3.0))
    size_field(f_th, sheet_edges(tr_parts), lv)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
    gmsh.model.mesh.generate(3)
    n_tet = len(gmsh.model.mesh.getElementsByType(11)[0])
    gmsh.write(path)
    gmsh.finalize()
    return n_tet, (k_up, k_dn), er


def write_config(c, level, mesh_name, split, er, path):
    k_up, k_dn = split
    strip = k_up > 0

    def port(i):
        base = 10 + 20 * (i - 1)
        # Trace positive: strips to the plane above point up, strips to the plane below point down.
        els = [{"Attributes": [base + j], "Direction": "+Z" if j < k_up else "-Z"} for j in range(k_up + k_dn)]
        return {"Index": i, "R": 50.0, "Elements": els} if len(els) > 1 else {"Index": i, "R": 50.0, **els[0]}
    p1 = port(1)
    p1["Excitation"] = 1
    cfg = {
        "Problem": {"Type": "Driven", "Verbose": 1, "Output": "postpro"},
        "Model": {"Mesh": mesh_name, "L0": 1.0e-3},
        "Domains": {"Materials": [{"Attributes": [1], "Permeability": 1.0, "Permittivity": er, "LossTan": c["tanD"]}]
                    + ([] if strip else [{"Attributes": [2], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0}])},
        "Boundaries": {"PEC": {"Attributes": [2]}, "Absorbing": {"Attributes": [3], "Order": 1},
                       "LumpedPort": [p1] + [port(i) for i in range(2, (4 if "diff" in c else 2) + 1)]},
        "Solver": {
            "Order": LEVELS[level]["order"], "Device": "CPU",
            "Driven": {"Samples": [{"Type": "Linear", "MinFreq": c["f"]["step"], "MaxFreq": c["f"]["max"], "FreqStep": c["f"]["step"]}],
                       "AdaptiveTol": 1.0e-4, "Save": []},
            "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-9, "MaxIts": 400},
        },
    }
    with open(path, "w") as f:
        json.dump(cfg, f, indent=2)


if __name__ == "__main__":
    if len(sys.argv) not in (4, 5) or (len(sys.argv) == 5 and sys.argv[4] != "half"):
        sys.exit(__doc__)
    cid, level, out = sys.argv[1:4]
    half = len(sys.argv) == 5
    c = load_case(cid)
    os.makedirs(out, exist_ok=True)
    n_tet, split, er = build(c, level, os.path.join(out, "mesh.msh"), half)
    write_config(c, level, "mesh.msh", split, er, os.path.join(out, "config.json"))
    print(json.dumps({"case": cid, "thru": "half" if half else True, "level": level, "tetrahedra": n_tet, "gmsh": gmsh.__version__}))
