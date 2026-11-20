# Impedora validation data

The evidence behind Impedora's via model: every 3D
case it is compared with, the pipeline that solved them, the stored 3D results, and Impedora's own
output for the same cases, so anyone can rerun the comparison and check it.

Impedora is a set of free PCB signal-integrity calculators that run in the browser. Its via model is
not a 3D solver: it is a fast physics-based model (the radial-mode cavity method, with the pads and
traces attached) that is compared openly with one. The calculators' code is closed; this data is open.

## What is here

| Path | Contents |
|---|---|
| `fullwave/cases.json` | Every 3D case: stackup, via, pads, antipads, stitching vias, traces, frequency sweep. Lengths in mm. |
| `fullwave/results.json` | The stored 3D solver results (Palace, finite elements), per case and mesh level (`coarse`, `fine`, `edge`): frequencies in GHz and the S-parameters of the excited ports, real and imaginary. |
| `fullwave/*.py`, `fullwave/Dockerfile` | The pipeline that made them: the 3D model builder (Gmsh), the solver run and the container recipe (Palace). `fullwave/README.md` says how to rerun a case. |
| `impedora/` | Impedora's S-parameters for each case on the same frequencies (Touchstone), as compared on Impedora's Validation page. |
| `impedora/comparison.csv` | Frequency by frequency: the largest difference \|S_Impedora − S_3D\| over the compared entries, and \|S11\|, \|S21\| of both in dB. |
| `reference/` | Independent second-order finite-element checks (Gmsh with scikit-fem or SciPy): trace cross-sections, a trace over a slotted plane, a blind via's static capacitance, and the full-wave response of a round via (body of revolution, 1 to 50 GHz). |
| `docs/` | How the comparison is made and what it shows (`VALIDATION.md`), and the via model itself (`PHYSICS.md`). |

## How to read the comparison

- Cases with their ports in the antipads (`p2g4`, `p2g2`, …) compare the via alone, S11 and S21
  (Sdd for the pair `p2d2`), from 1 to 50 GHz.
- Cases with pads and traces (`ms…`, `bl…`, `bd…`, `th…`, `md1`) have their ports at the far end of
  each trace. Each has two trace-only runs (`<case>-thru`, `<case>-thru-half`, or those of the case in
  its `thru` field) that calibrate the 3D ports and the 3D trace out; Impedora's files for these cases
  include that calibration, exactly as on the Validation page. `docs/VALIDATION.md` explains it.
- Where a case is solved on two meshes, the difference between them is the 3D result's own
  uncertainty.

## Contributing reference data

Measured S-parameters or 3D simulations of traces or vias, with full geometry, are the most useful
thing anyone can add. See `CONTRIBUTING.md`.

## Licence

- **Data and documents** (`impedora/`, `docs/`, `fullwave/cases.json`, `fullwave/results.json` and
  every Markdown file): Creative Commons Attribution 4.0 International, in `LICENSE`. Reuse them for
  any purpose, commercial included, with credit to "Impedora validation data, impedora.com".
- **Code** (every `.py` file and `fullwave/Dockerfile`): MIT, in `LICENSE-MIT`.

The 3D solver itself, [Palace](https://github.com/awslabs/palace), is Apache-2.0 and is installed by
the container recipe; it is not part of this repository.
