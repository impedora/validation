#!/usr/bin/env python3
"""
Run the full-wave via cases inside the container and store the S-parameters in
fullwave/results.json (the file Impedora's tests read).

Usage (normally through `npm run fullwave`, which starts the container):
  python3 fullwave/run.py                    every case on the coarse mesh that has no result yet
  python3 fullwave/run.py p2g4 p4g4          selected cases
  python3 fullwave/run.py --level fine       the fine mesh (slow); results.json holds both meshes
  python3 fullwave/run.py --level both       both meshes
  python3 fullwave/run.py --np 4             MPI processes (default: all CPUs)
  python3 fullwave/run.py p4g4 --force       solve again although a result is stored

A case with a trace ("launch" or "exit") has two companions, "<id>-thru" and "<id>-thru-half": the
same trace between two of the same trace ports, without the via (trace_thru.py), at the case's
trace length and at half of it. Together they calibrate the ports and measure the 3D model's own
trace (src/physics/fullwave.ts, portBox). They are listed and run like any other case, and naming
the case runs its thrus too.

A case and mesh level that already has a stored result is skipped, so adding a case and running
the whole set solves only the new one. Every run is saved as soon as it finishes, so an
interrupted batch keeps what it has done.
Working files go to fullwave/runs/ (not committed).
"""
import argparse
import cmath
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys

from via_case import LEVELS  # noqa: E402 (same folder)
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results.json")
RUNS = os.path.join(HERE, "runs")


def palace_version():
    exe = shutil.which("palace")
    if not exe:
        sys.exit("palace not found. Run this inside the Impedora full-wave container (npm run fullwave).")
    out = subprocess.run([exe, "--version"], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("Palace version:"):
            return "Palace " + line.split(":", 1)[1].strip().replace("-dirty", "")
    return f"Palace {os.environ.get('PALACE_VERSION', 'unknown version')}"


def read_s(path):
    """port-S.csv -> (freqs in GHz, excited ports, S[f][row][col] as [re, im]).

    Palace writes one column of S per excited port. Only those columns are stored: col runs over
    `excited` (1-based port numbers), row over all ports.
    """
    with open(path) as f:
        rows = [r for r in csv.reader(f) if r]
    head = [h.strip() for h in rows[0]]
    cols = {}
    for k, h in enumerate(head):
        mt = re.match(r"\|S\[(\d+)\]\[(\d+)\]\|", h)
        if mt:
            cols[(int(mt.group(1)), int(mt.group(2)))] = (k, k + 1)  # |S| dB, then arg(S) deg
    n = max(i for i, _ in cols)
    excited = sorted({j for _, j in cols})
    freqs, S = [], []
    for r in rows[1:]:
        v = [float(x) for x in r]
        freqs.append(v[0])
        m = [[[0.0, 0.0] for _ in excited] for _ in range(n)]
        for (i, j), (km, kp) in cols.items():
            z = cmath.rect(10 ** (v[km] / 20), math.radians(v[kp]))
            m[i - 1][excited.index(j)] = [z.real, z.imag]
        S.append(m)
    return freqs, excited, S


def load_results():
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            return json.load(f)
    return {"results": []}


def save_run(case_id, level, solver, run):
    data = load_results()
    data["_comment"] = "Written by npm run fullwave. Do not edit by hand."
    entry = next((r for r in data["results"] if r["id"] == case_id), None)
    if entry is None:
        entry = {"id": case_id, "solver": solver, "runs": {}}
        data["results"].append(entry)
    entry["solver"] = solver
    entry["runs"][level] = run
    data["results"].sort(key=lambda r: r["id"])
    tmp = RESULTS + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=1)
        f.write("\n")
    os.replace(tmp, RESULTS)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases", nargs="*")
    ap.add_argument("--level", choices=[*LEVELS, "both"], default="coarse", help="mesh level (via_case.py LEVELS); both = coarse and fine")
    ap.add_argument("--np", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--force", action="store_true", help="solve again even if a result is already stored")
    a = ap.parse_args()

    with open(os.path.join(HERE, "cases.json")) as f:
        cases = json.load(f)["cases"]
    all_ids = []
    for c in cases:
        all_ids.append(c["id"])
        if (c.get("launch") or c.get("exit")) and not c.get("thru"):  # "thru": share another case's thrus
            all_ids += [c["id"] + "-thru", c["id"] + "-thru-half"]
    ids = a.cases or all_ids
    for i in ids:
        if i not in all_ids:
            sys.exit(f"Unknown case '{i}'. Known: {', '.join(all_ids)}")
    # A case with its own thrus brings them along when it is named (they calibrate its trace ports).
    for c in cases:
        if a.cases and c["id"] in ids and (c.get("launch") or c.get("exit")) and not c.get("thru"):
            ids += [x for x in (c["id"] + "-thru", c["id"] + "-thru-half") if x not in ids]
    levels = ["coarse", "fine"] if a.level == "both" else [a.level]
    solver = palace_version()
    print(f"{solver}, {a.np} MPI processes")

    stored = {r["id"]: r["runs"] for r in load_results()["results"]}
    for cid in ids:
        for lv in levels:
            if not a.force and lv in stored.get(cid, {}):
                print(f"→ {cid} [{lv}] already stored, skipped (--force to solve again)")
                continue
            d = os.path.join(RUNS, cid, lv)
            shutil.rmtree(d, ignore_errors=True)
            m = re.fullmatch(r"(.+)-thru(-half)?", cid)
            script, args = ("trace_thru.py", [m[1], lv, d] + (["half"] if m[2] else [])) if m else ("via_case.py", [cid, lv, d])
            gen = subprocess.run([sys.executable, os.path.join(HERE, script), *args],
                                 capture_output=True, text=True)
            if gen.returncode:
                sys.exit(f"Mesh generation failed for {cid}/{lv}:\n{gen.stderr}")
            info = json.loads(gen.stdout.strip().splitlines()[-1])
            tets = info["tetrahedra"]
            tool = f"{solver}, Gmsh {info['gmsh']}"
            print(f"→ {cid} [{lv}] {tets} tetrahedra … ", end="", flush=True)
            t0 = time.time()
            log = open(os.path.join(d, "palace.log"), "w")
            p = subprocess.run(["palace", "-np", str(a.np), "config.json"], cwd=d, stdout=log, stderr=subprocess.STDOUT)
            log.close()
            dt = time.time() - t0
            s_csv = os.path.join(d, "postpro", "port-S.csv")
            if p.returncode or not os.path.exists(s_csv):
                print("failed")
                sys.exit(f"Palace failed for {cid}/{lv}. Log: fullwave/runs/{cid}/{lv}/palace.log")
            freqs, excited, S = read_s(s_csv)
            save_run(cid, lv, tool, {"tetrahedra": tets, "seconds": round(dt, 1), "excited": excited, "f": freqs, "S": S})
            print(f"done in {dt / 60:.1f} min")
    print("Saved to fullwave/results.json. Check the comparison with: npm run diagnose")


if __name__ == "__main__":
    main()
