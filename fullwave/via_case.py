#!/usr/bin/env python3
"""
Build the Palace model (mesh + config) for one via case in fullwave/cases.json.

The full-wave model uses the same idealisation as Impedora's physics-based via model, so the
comparison tests the physics (plate wave, radial modes, barrel-to-plane capacitance) and not
differences in how the problem is set up:

  * a stack of solid copper planes (PEC, real thickness) with dielectric between them,
  * one signal via barrel (PEC cylinder) through the whole stack, with an antipad in every plane,
  * stitching vias (PEC cylinders) that short every plane together,
  * a coaxial lumped port in the antipad annulus on the top face of the top plane and on the
    bottom face of the bottom plane (barrel to plane, 50 ohm), which is where Impedora's ports sit;
    a differential pair has four (via 1 top, via 2 top, via 1 bottom, via 2 bottom). With
    "port2Plane", port 2 sits on the top face of an inner plane instead and the rest of the barrel
    is an open stub,
  * with "launch": a layer of the same dielectric on top of the top plane with a microstrip on it.
    The barrel comes up to that layer and ends in a pad; the trace runs from the pad towards −x.
    Port 1 is then a 50 ohm sheet under the far end of the trace (trace to plane), not a coaxial
    port, and there is air above with an absorbing top. Trace and pad are zero-thickness sheets.
    A differential pair (vias on the x axis) gets one pad and trace per via, the traces running
    side by side towards −y, and ports 1 and 2 under their ends,
  * with "exit": the barrel ends inside the lowest plane pair (a blind via, or a backdrilled one with
    a short stub left), in a pad with a stripline running towards −x at that depth. Port 2 is then
    two 50 ohm sheets at the far end of the trace, one up to the plane above and one down to the
    plane below, in parallel. The plane below is solid. With "pair" (counted from the top) the
    stripline sits in a higher plane pair and the barrel runs on through the planes below it, as a
    stub, to an open end at the bottom face of the bottom plane (PMC across its antipad),
  * a first-order absorbing boundary on the outer cylindrical wall. The via model is given the
    same wall for the comparison (`edge` in PhysViaInput), so the domain size drops out.

Usage:  python3 via_case.py <case-id> <level> <out-dir>      level = coarse | fine
"""
import json
import math
import os
import sys

import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))

# Mesh levels. Sizes are in mm. "near" applies at via surfaces, "far" at the domain edge.
# The fine level has about 2.5 times as many elements; comparing the two gives the numerical
# uncertainty of the reference itself.
LEVELS = {
    "coarse": {"near": 0.06, "far": 0.40, "order": 2},
    "fine": {"near": 0.045, "far": 0.30, "order": 2},
    # Fine, plus a finer mesh along the edges of the trace and pad sheets: the field is singular at a
    # sheet's edge, and the two levels above leave the capacitance where a trace leaves its pad several
    # fF too high (docs/VALIDATION.md). "edge" is the size there (mm).
    "edge": {"near": 0.045, "far": 0.30, "order": 2, "edge": 0.008},
}


def size_field(base, curves, lv):
    """Use the size field `base`, refined along `curves` (the sheet edges) when the level has an edge size."""
    e = lv.get("edge")
    if not e or not curves:
        gmsh.model.mesh.field.setAsBackgroundMesh(base)
        return
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", sorted(set(curves)))
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 300)
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", e)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", lv["near"])
    gmsh.model.mesh.field.setNumber(ft, "DistMin", 0.5 * e)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", 0.06)
    gmsh.model.mesh.field.setNumber(ft, "StopAtDistMax", 1)  # beyond it, the base field alone
    fm = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fm, "FieldsList", [base, ft])
    gmsh.model.mesh.field.setAsBackgroundMesh(fm)


def sheet_edges(surfaces):
    """Boundary curves of a set of sheet surfaces."""
    return [abs(t) for _, t in gmsh.model.getBoundary([(2, s) for s in surfaces], combined=False, oriented=False)]


def strip_split(d, h):
    """Strips (up, down) for a stripline port at depth d in a plane pair of height h.

    A stripline returns (h − d)/h of its current through the plane above and d/h through the plane
    below. Palace gives every element of a multi-element lumped port the same resistance (n·R), so
    the current splits by element count. The port's two sheets are therefore cut into side-by-side
    strips in that ratio; with one sheet each, the port would force a 50/50 split and launch a
    parallel-plate wave between the planes (it showed up as radiation loss in a lossless stripline).
    """
    for n in range(2, 9):
        for k_up in range(1, n):
            k_dn = n - k_up
            if abs(k_up / n - (h - d) / h) < 0.005:
                return k_up, k_dn
    sys.exit(f"exit: a stripline at depth {d} in {h} needs a return split that 8 strips cannot represent; pick a depth with a simple ratio.")


def load_case(case_id):
    with open(os.path.join(HERE, "cases.json")) as f:
        cases = json.load(f)["cases"]
    for c in cases:
        if c["id"] == case_id:
            return c
    sys.exit(f"Unknown case '{case_id}'. Known: {', '.join(c['id'] for c in cases)}")


FENCE_GAP = 0.25  # mm from the trace edge to the centre of a fence via (see port_fence)


def port_fence(c, x_port):
    """Ground vias around a stripline port at x = x_port: one on each side of the trace and one behind
    its end. They short the plane pair right at the port, so the little parallel-plate field the
    lumped port stirs up stays local. A lone lumped port between two planes leaks power into the
    plane pair, and the leak depends on what else is nearby, which would make the port's
    calibration (trace_thru.py) not carry over to the via case."""
    w = c["exit"]["w"]
    yf = w / 2 + FENCE_GAP
    back = x_port - yf if x_port < 0 else x_port + yf
    return [(x_port, yf), (x_port, -yf), (back, 0.0)]


def signal_positions(c):
    """One via at the origin, or a differential pair on the x axis."""
    if "diff" in c:
        return [(-c["diff"]["pitch"] / 2, 0.0), (c["diff"]["pitch"] / 2, 0.0)]
    return [(0.0, 0.0)]


def ground_positions(c):
    """Same layout rules as viaPositions() in src/physics/viaModel.ts (fullwave.ts passes the angle and
    the fence vias to it)."""
    n, d = c["ground"]["n"], c["ground"]["dist"]
    if "diff" not in c:
        # With a stripline exit the trace runs along −x inside the plane pair, where a stitching via
        # would sit. The ring is turned by half a step, and the port gets its fence.
        a0 = math.pi / n if (c.get("exit") and n) else 0.0
        ring = [(d * math.cos(a0 + 2 * math.pi * i / n), d * math.sin(a0 + 2 * math.pi * i / n)) for i in range(n)]
        return ring + (port_fence(c, -c["exit"]["len"]) if c.get("exit") else [])
    half = c["diff"]["pitch"] / 2
    g = []
    if n == 1:
        g.append((0.0, d))
    if n >= 2:
        g += [(-(half + d), 0.0), (half + d, 0.0)]
    if n == 4:
        g += [(0.0, d), (0.0, -d)]
    return g


def mesh_level(c, level):
    """Mesh sizes for a case: the level's sizes times the case's optional "meshScale"."""
    k = c.get("meshScale", 1.0)
    lv = dict(LEVELS[level])
    lv["near"] *= k
    lv["far"] *= k
    return lv


def build_mesh(c, level, path):
    lv = mesh_level(c, level)
    t = c["planeT"]
    hs = [cv["h"] for cv in c["cavities"]]  # listed from the top, as in cases.json and the via model
    n_planes = len(hs) + 1
    R = c["domainR"]
    rb, ra = c["drill"] / 2, c["antipad"] / 2

    # z of the bottom face of each plane, building upwards from z = 0 at the bottom of the stack.
    # z_planes[0] is the bottom plane; the cavity above it is the last one in the list.
    z_planes, z = [], 0.0
    for i in range(n_planes):
        z_planes.append(z)
        z += t + (hs[n_planes - 2 - i] if i < n_planes - 1 else 0.0)
    T = z  # total height = top face of the top plane
    # Plane spacing seen from the top must match the case, cavity by cavity.
    gaps = [z_planes[n_planes - 1 - q] - (z_planes[n_planes - 2 - q] + t) for q in range(n_planes - 1)]
    assert all(abs(g - h) < 1e-9 for g, h in zip(gaps, hs)), (gaps, hs)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(c["id"])
    occ = gmsh.model.occ

    la = c.get("launch")
    pair_launch = bool(la) and len(signal_positions(c)) == 2
    zs = T + (la["h"] if la else 0.0)  # top of the dielectric: the trace layer when there is a launch
    ex = c.get("exit")
    if ex and (len(signal_positions(c)) != 1 or la or c.get("port2Plane") is not None):
        sys.exit("An exit is supported for a single via with a coaxial port on top.")
    pair = ex.get("pair", len(hs) - 1) if ex else len(hs) - 1
    thru = bool(ex) and pair < len(hs) - 1      # the barrel runs on below the stripline's plane pair
    if thru and ex.get("stub", 0) > 0:
        sys.exit("exit: a stub length is for an exit in the lowest pair; above it the barrel runs through.")
    z_lo, z_hi = z_planes[n_planes - 2 - pair] + t, z_planes[n_planes - 1 - pair]  # the exit's plane pair: faces of the dielectric
    z_ex = z_hi - ex["depth"] if ex else None    # trace level inside it
    z_end = z_ex - ex["stub"] if (ex and not thru) else -0.01   # lower end of the signal barrel
    if ex and not (z_lo + 0.02 < z_ex < z_hi - 0.02 and (thru or z_lo + 0.02 < z_end <= z_ex)):
        sys.exit("exit: the trace and the barrel end must lie inside the plane pair.")
    domain = occ.addCylinder(0, 0, 0, 0, 0, zs, R)
    cut = []
    sigs = signal_positions(c)
    for q_, zp in enumerate(z_planes):
        slab = occ.addCylinder(0, 0, zp, 0, 0, t, R * 1.01)
        # The plane below an exit in the lowest pair is solid: the barrel does not reach it.
        holes = [] if (ex and not thru and q_ == 0) else [(3, occ.addCylinder(x, y, zp, 0, 0, t, ra)) for (x, y) in sigs]
        out = [(3, slab)] if not holes else occ.cut([(3, slab)], holes)[0]
        cut += out
    # Signal barrels end at the trace layer, stitching vias at the top plane.
    for (x, y) in sigs:
        cut.append((3, occ.addCylinder(x, y, z_end, 0, 0, zs + (0.0 if la else 0.01) - z_end, rb)))
    for (x, y) in ground_positions(c):
        cut.append((3, occ.addCylinder(x, y, -0.01, 0, 0, T + (0.01 if la else 0.02), rb)))
    fluid, _ = occ.cut([(3, domain)], cut)
    air_vol = None
    if la:
        # Air above, the trace and pad as sheets on the interface, the port sheet under the trace end.
        air = occ.addCylinder(0, 0, zs, 0, 0, la["air"], R)
        rp, w, ln = la["pad"] / 2, la["w"], la["len"]
        if not pair_launch:
            sheet, _ = occ.fuse([(2, occ.addDisk(0, 0, zs, rp, rp))], [(2, occ.addRectangle(-ln, -w / 2, zs, ln, w))])
            pp = [occ.addPoint(-ln, -w / 2, T), occ.addPoint(-ln, w / 2, T), occ.addPoint(-ln, w / 2, zs), occ.addPoint(-ln, -w / 2, zs)]
            port_sheets = [occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(pp[i], pp[(i + 1) % 4]) for i in range(4)])])]
        else:
            # One pad and trace per via; the traces run side by side towards −y, port sheets under their ends.
            parts_ = [(2, occ.addDisk(x, 0, zs, rp, rp)) for (x, _) in sigs] + [(2, occ.addRectangle(x - w / 2, -ln, zs, w, ln)) for (x, _) in sigs]
            sheet, _ = occ.fuse(parts_[:1], parts_[1:])
            port_sheets = []
            for (x, _) in sigs:
                pp = [occ.addPoint(x - w / 2, -ln, T), occ.addPoint(x + w / 2, -ln, T), occ.addPoint(x + w / 2, -ln, zs), occ.addPoint(x - w / 2, -ln, zs)]
                port_sheets.append(occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(pp[i], pp[(i + 1) % 4]) for i in range(4)])]))
        fluid, parts = occ.fragment(fluid + [(3, air)], sheet + [(2, q_) for q_ in port_sheets])
        # The sheets are inside the dielectric, not on its boundary.
        port_sheet_parts = [e for q_ in parts[len(parts) - len(port_sheets):] for e in q_ if e[0] == 2]

    # Port 2 inside the stack ("port2Plane": q, counted from the top plane): an annulus across the
    # antipad on the top face of that plane, embedded in the volume. The barrel below it is a stub,
    # and its end at the bottom face of the bottom plane is left open (PMC).
    q2 = c.get("port2Plane")
    z_bot_port = 0.0
    if q2 is not None:
        z_bot_port = z_planes[n_planes - 1 - q2] + t
        rings = []
        for (x, y) in sigs:
            ring, _ = occ.cut([(2, occ.addDisk(x, y, z_bot_port, ra, ra))], [(2, occ.addDisk(x, y, z_bot_port, rb, rb))])
            rings += ring
        fluid, _ = occ.fragment(fluid, rings)
    occ.synchronize()

    ex_parts = {"cu": [], "strips": []}  # strips: one list of faces per port element, upper strips first
    if ex:
        # Trace and pad as a sheet at the exit depth, and the two port sheets at the trace end.
        rp, w, ln = ex["pad"] / 2, ex["w"], ex["len"]
        sheet, _ = occ.fuse([(2, occ.addDisk(0, 0, z_ex, rp, rp))], [(2, occ.addRectangle(-ln, -w / 2, z_ex, ln, w))])
        if thru or ex["stub"] > 0:  # the barrel passes through: the pad is a ring around it
            sheet, _ = occ.cut(sheet, [(2, occ.addDisk(0, 0, z_ex, rb, rb))])

        def wall(za, zb, ya, yb):
            pp = [occ.addPoint(-ln, ya, za), occ.addPoint(-ln, yb, za), occ.addPoint(-ln, yb, zb), occ.addPoint(-ln, ya, zb)]
            return occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(pp[i], pp[(i + 1) % 4]) for i in range(4)])])

        k_up, k_dn = strip_split(z_hi - z_ex, z_hi - z_lo)
        cuts = lambda k: [(-w / 2 + i * w / k, -w / 2 + (i + 1) * w / k) for i in range(k)]
        strips = [wall(z_ex, z_hi, ya, yb) for ya, yb in cuts(k_up)] + [wall(z_lo, z_ex, ya, yb) for ya, yb in cuts(k_dn)]
        fluid, parts = occ.fragment(fluid, sheet + [(2, q_) for q_ in strips])
        n_sheet, n_str = len(sheet), len(strips)
        ex_parts["cu"] = [e for q_ in parts[len(parts) - n_str - n_sheet:len(parts) - n_str] for e in q_ if e[0] == 2]
        ex_parts["strips"] = [[e for e in q_ if e[0] == 2] for q_ in parts[len(parts) - n_str:]]
        ex_parts["k_up"] = k_up
        occ.synchronize()
    vols = [tag for (dim, tag) in fluid if dim == 3]
    if la:
        air_vol = [v for v in vols if gmsh.model.getBoundingBox(3, v)[2] > zs - 1e-6]
        vols = [v for v in vols if v not in air_vol]
        assert len(air_vol) == 1 and vols, "air / dielectric volumes not found"
    # Boundaries volume by volume, so a face shared by two volumes (an embedded port) is kept.
    surfs = gmsh.model.getBoundary([(3, v) for v in vols + (air_vol or [])], combined=False, oriented=False)
    if la:
        surfs = surfs + port_sheet_parts
    if ex:
        # Sheets inside the dielectric are not on its boundary, so they are added by hand.
        ex_parts["cu"] = sorted({abs(e[1]) for e in ex_parts["cu"]})
        ex_parts["strips"] = [sorted({abs(e[1]) for e in q_}) for q_ in ex_parts["strips"]]
        surfs = surfs + [(2, tg) for tg in ex_parts["cu"]] + [(2, tg) for q_ in ex_parts["strips"] for tg in q_]
    surfs = sorted({abs(tag) for (_, tag) in surfs})

    # Ports: top apertures first (one per signal via), then bottom apertures, as in the via model.
    ports = [[] for _ in range(2 * len(sigs))]
    absorbing, pec, open_end, la_cu = [], [], [], []
    port2_more = [[] for _ in ex_parts["strips"][1:]] if ex else []  # port 2 elements after the first
    eps = 1e-6
    for s in surfs:
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, s)
        flat = abs(z1 - z0) < eps
        near = None
        if flat:
            for i, (sx, sy) in enumerate(sigs):
                if max(math.hypot(x - sx, y - sy) for x in (x0, x1) for y in (y0, y1)) <= ra * 1.5:
                    near = i
        if ex and s in ex_parts["strips"][0]:
            ports[1].append(s)  # first strip, trace to the plane above
        elif ex and any(s in q_ for q_ in ex_parts["strips"][1:]):
            port2_more[next(i for i, q_ in enumerate(ex_parts["strips"][1:]) if s in q_)].append(s)
        elif ex and s in ex_parts["cu"]:
            pec.append(s)
        elif pair_launch and abs(y0 + la["len"]) < eps and abs(y1 + la["len"]) < eps and z0 > T - eps and z1 < zs + eps:
            ports[min(range(len(sigs)), key=lambda i: abs((x0 + x1) / 2 - sigs[i][0]))].append(s)  # under a trace end
        elif la and not pair_launch and abs(x0 + la["len"]) < eps and abs(x1 + la["len"]) < eps and z0 > T - eps and z1 < zs + eps:
            ports[0].append(s)  # the sheet under the trace end
        elif la and flat and abs(z0 - zs) < eps:
            # The dielectric-air interface: copper where the trace and pad are, nothing elsewhere.
            if pair_launch:
                xs = [x for (x, _) in sigs]
                if x0 > min(xs) - la["pad"] / 2 - eps and x1 < max(xs) + la["pad"] / 2 + eps and y0 > -la["len"] - eps and y1 < la["pad"] / 2 + eps:
                    pec.append(s)
                    la_cu.append(s)
            elif x0 > -la["len"] - eps and x1 < la["pad"] / 2 + eps and max(abs(y0), abs(y1)) < max(la["pad"], la["w"]) / 2 + eps:
                pec.append(s)
                la_cu.append(s)
        elif la and flat and abs(z0 - (zs + la["air"])) < eps:
            absorbing.append(s)  # top of the air region
        elif near is not None and abs(z0 - T) < eps and not la:
            ports[near].append(s)
        elif thru and near is not None and abs(z0) < eps:
            open_end.append(s)  # the stub's end at the bottom of the stack
        elif near is not None and abs(z0 - z_bot_port) < eps:
            ports[near + len(sigs)].append(s)
        elif near is not None and abs(z0) < eps:
            open_end.append(s)
        elif not flat and min(abs(x0), abs(x1), abs(y0), abs(y1)) > R * 0.9:
            absorbing.append(s)
        else:
            pec.append(s)
    if os.environ.get("VIA_CASE_DEBUG"):
        for s_ in surfs:
            print([round(v, 3) for v in gmsh.model.getBoundingBox(2, s_)], "port" if any(s_ in pt for pt in ports) else "abs" if s_ in absorbing else "pec" if s_ in pec else "open" if s_ in open_end else "-", file=sys.stderr)
    assert all(ports) and absorbing and all(port2_more), "surface classification failed"
    assert bool(open_end) == (q2 is not None or thru), "open stub end not found"

    gmsh.model.addPhysicalGroup(3, vols, 1, "dielectric")
    if air_vol:
        gmsh.model.addPhysicalGroup(3, air_vol, 2, "air")
    gmsh.model.addPhysicalGroup(2, pec, 2, "pec")
    gmsh.model.addPhysicalGroup(2, absorbing, 3, "absorbing")
    for i, pt in enumerate(ports):
        gmsh.model.addPhysicalGroup(2, pt, 4 + i, f"port{i + 1}")
    if open_end:
        gmsh.model.addPhysicalGroup(2, open_end, 4 + len(ports), "open")
    for i, q_ in enumerate(port2_more):  # attributes 20, 21, … (see write_config)
        gmsh.model.addPhysicalGroup(2, q_, 20 + i, f"port2-strip{i + 2}")

    # Size field: fine near every via barrel and inside the antipads, growing outwards.
    via_surfs = []
    for s in pec:
        x0, y0, z0, x1, y1, _ = gmsh.model.getBoundingBox(2, s)
        if ((x1 - x0) < 2.2 * ra and (y1 - y0) < 2.2 * ra) or (la and abs(z0 - zs) < eps) or (ex and s in ex_parts["cu"]):
            via_surfs.append(s)  # barrels, antipad walls, and the trace and pad
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", via_surfs + [q for pt in ports for q in pt] + [q for st in port2_more for q in st])
    f_th = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_th, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_th, "SizeMin", lv["near"])
    gmsh.model.mesh.field.setNumber(f_th, "SizeMax", lv["far"])
    gmsh.model.mesh.field.setNumber(f_th, "DistMin", 0.5 * ra)
    gmsh.model.mesh.field.setNumber(f_th, "DistMax", 0.5 * R)
    size_field(f_th, sheet_edges(ex_parts["cu"] if ex else la_cu), lv)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    # Second-order (curved) elements: barrels and antipads are true circles, not polygons. With straight
    # elements the coaxial port's radii are off by about 1%, which shows up directly in S.
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
    gmsh.model.mesh.generate(3)
    n_tet = len(gmsh.model.mesh.getElementsByType(11)[0])  # 10-node tetrahedra
    gmsh.write(path)
    gmsh.finalize()
    return n_tet


def exit_thru(c):
    """True when the exit sits above the lowest plane pair, so the barrel runs on below it."""
    ex = c.get("exit")
    return bool(ex) and ex.get("pair", len(c["cavities"]) - 1) < len(c["cavities"]) - 1


def exit_elements(c):
    """Port 2 of a stripline exit: strips up to the plane above (+Z), then down to the plane below (−Z)."""
    h = c["cavities"][c["exit"].get("pair", len(c["cavities"]) - 1)]["h"]
    k_up, k_dn = strip_split(c["exit"]["depth"], h)
    first = [{"Attributes": [5], "Direction": "+Z"}]
    more = [{"Attributes": [20 + i], "Direction": "+Z" if i + 1 < k_up else "-Z"} for i in range(k_up + k_dn - 1)]
    return first + more


def write_config(c, level, mesh_name, path):
    lv = LEVELS[level]
    cfg = {
        "Problem": {"Type": "Driven", "Verbose": 1, "Output": "postpro"},
        "Model": {"Mesh": mesh_name, "L0": 1.0e-3},
        "Domains": {
            "Materials": [
                {
                    "Attributes": [1],
                    "Permeability": 1.0,
                    "Permittivity": c["cavities"][0]["er"],
                    "LossTan": c["tanD"],
                }
            ] + ([{"Attributes": [2], "Permeability": 1.0, "Permittivity": 1.0, "LossTan": 0.0}] if c.get("launch") else [])
        },
        "Boundaries": {
            "PEC": {"Attributes": [2]},
            "Absorbing": {"Attributes": [3], "Order": 1},
            **({"PMC": {"Attributes": [4 + 2 * len(signal_positions(c))]}} if (c.get("port2Plane") is not None or exit_thru(c)) else {}),
            # Only port 1 is excited (one column of S). The cases are mirror-symmetric, top to bottom and
            # between the two vias of a pair, so the other columns hold nothing new. A case that is not
            # symmetric can set "exciteAll": true.
            "LumpedPort": [
                dict({"Index": i + 1, "R": 50.0},
                     # Signal conductor positive everywhere: at a coaxial port the field points outwards
                     # from the barrel, at a trace port away from the trace towards its plane.
                     **({"Elements": exit_elements(c)} if (c.get("exit") and i == 1) else
                        {"Attributes": [4 + i], "Direction": "-Z" if (c.get("launch") and i < len(signal_positions(c))) else "+R"}),
                     **({"Excitation": i + 1} if (c.get("exciteAll") or i == 0) else {}))
                for i in range(2 * len(signal_positions(c)))
            ],
        },
        "Solver": {
            "Order": lv["order"],
            "Device": "CPU",
            "Driven": {
                "Samples": [
                    {
                        "Type": "Linear",
                        "MinFreq": c["f"]["step"],
                        "MaxFreq": c["f"]["max"],
                        "FreqStep": c["f"]["step"],
                    }
                ],
                "AdaptiveTol": 1.0e-4,
                "Save": [],
            },
            "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-9, "MaxIts": 400},
        },
    }
    with open(path, "w") as f:
        json.dump(cfg, f, indent=2)


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    case_id, level, out = sys.argv[1:]
    if level not in LEVELS:
        sys.exit(f"level must be one of {', '.join(LEVELS)}")
    c = load_case(case_id)
    if len(set(cv["er"] for cv in c["cavities"])) != 1:
        sys.exit("This generator supports one dielectric constant per case.")
    os.makedirs(out, exist_ok=True)
    n_tet = build_mesh(c, level, os.path.join(out, "via.msh"))
    write_config(c, level, "via.msh", os.path.join(out, "config.json"))
    print(json.dumps({"case": case_id, "level": level, "tetrahedra": n_tet, "gmsh": gmsh.__version__}))


if __name__ == "__main__":
    main()
