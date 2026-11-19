> Copied from Impedora's repository, where it is docs/VALIDATION.md. File paths under `src/` and `scripts/` refer to that repository; the 3D cases and results are in `fullwave/` here.

# Validation

All reference cases live in `src/physics/reference.ts`. Each has a source, an expected value and a tolerance. They run:

- in `npm test` (one test per case),
- in `npm run diagnose` (got / expected / error % table, exit code 1 on failure),
- on the app's **Validation** page, live in the visitor's browser.

Extra unit tests (Bessel functions vs SciPy, Touchstone round-trip, KiCad parser, stackup geometry) live next to the code as `*.test.ts`.

## What is checked

| Model | Checked against |
|---|---|
| Elliptic integrals | Abramowitz & Stegun |
| Trace field solver | Exact zero-thickness stripline and coupled stripline (within 0.01%); exact zero-thickness coplanar waveguide; Hammerstad–Jensen microstrip in air (0.01%); independent finite-element solver for thick, etched, masked and coplanar traces (next section, within 0.007%); grid convergence; mirror symmetry |
| Trace loss | Exact TEM dielectric attenuation (0.005%); shape derivative vs really shrinking the conductors (0.006%); wide-stripline limit Rs/w; Kramers–Kronig slope of the Debye model; roughness formula |
| Trace S-parameters | Reciprocity, passivity, Sdd21 equals the odd-mode line, |S21| equals the loss when matched, Touchstone round trip |
| Microstrip formula | Pozar, *Microwave Engineering* 4e, Example 3.7 (Z0 within 0.1%, εeff within 0.6%); air-substrate limit |
| Stripline formula | Exact conformal-mapping solution (Wheeler within 0.3%); Pozar's wide-strip formula |
| Coupled stripline | Uncoupled limit equals single stripline exactly; Zodd < Zeven |
| Via, quick estimate | Johnson & Graham hand calculation; energy conservation; stub notch position |
| Bessel J, Y, H⁽²⁾, I0, K0 | SciPy 1.18 and Abramowitz & Stegun, real and complex arguments, error < 1e-9 |
| Capacitance solver | Exact coaxial capacitance per length (0.6%); grid refinement (0.1%) |
| Plate-pair impedance | Thin-wire loop inductance (exact); exact two-wire inductance with proximity effect, (μh/π)·acosh(s/2a), within 0.02% |
| Graf's addition theorem | Both sides evaluated numerically, error < 1e-9 |
| Radial modes | Static capacitance from the mode sum vs the field solver, two independent methods, within 1%; low-frequency limit equals the textbook plate impedance |
| Physics-based via | 3D field solver (next section); S21 → 0 dB at low frequency; reciprocity (S21 = S12 to 1e-15); passivity; more stitching vias never worse; mirror symmetry of a symmetric stackup; odd-mode symmetry for differential pairs; a stripline launch at the top end equals the flipped stack's bottom end (1e-6); a pair with its launch converts no mode |
| Current | IPC-2221 hand calculation; R = ρL/A |

## Trace field solver vs an independent solver

Exact solutions exist only for zero-thickness traces in one dielectric. For real cross-sections the field solver is compared with a second solver that shares no code or method with it: second-order finite elements on an unstructured triangle mesh refined towards the trace corners, the full cross-section with no symmetry assumed, built with Gmsh and scikit-fem (`scripts/reference/trace_fem.py`). Each reference value was computed on two meshes that differ by less than 0.0025% in Z0.

| Case | Impedora | Independent FEM | Difference |
|---|---|---|---|
| Microstrip, 300 µm wide, 35 µm copper, 150 µm of εr 4.2 | 47.0525 Ω | 47.0550 Ω | 0.005% |
| Same, under 25 µm solder mask (εr 3.5) | 45.2394 Ω | 45.2425 Ω | 0.007% |
| Microstrip pair, 150 µm traces and gap, odd / even | 47.0129 / 59.8268 Ω | 47.0147 / 59.8291 Ω | 0.004% |
| Same pair under 20 µm solder mask, odd / even | 43.8565 / 57.3554 Ω | 43.8594 / 57.3586 Ω | 0.007% |
| Offset stripline, two dielectrics, 18 µm copper | 50.5016 Ω | 50.5045 Ω | 0.006% |
| Offset stripline pair, odd / even | 44.0161 / 56.4365 Ω | 44.0175 / 56.4386 Ω | 0.004% |
| Embedded microstrip (dielectric over the trace) | 48.4297 Ω | 48.4330 Ω | 0.007% |
| Etched microstrip, 300 µm base, 270 µm top | 47.6810 Ω | 47.6838 Ω | 0.006% |
| Same, under 25 µm solder mask | 46.0189 Ω | 46.0219 Ω | 0.007% |
| Etched microstrip pair under solder mask, odd / even | 45.9763 / 58.7165 Ω | 45.9793 / 58.7197 Ω | 0.007% |
| Etched stripline pair, narrow face down, odd / even | 45.9863 / 58.2536 Ω | 45.9879 / 58.2558 Ω | 0.004% |
| Etched pair, mask 15 µm on copper and 30 µm on laminate, odd / even | 51.4662 / 75.0424 Ω | 51.4693 / 75.0466 Ω | 0.006% |
| Coplanar waveguide with ground, 150 µm gaps | 49.1852 Ω | 49.1868 Ω | 0.003% |
| Same, etched and under solder mask | 47.8998 Ω | 47.9025 Ω | 0.006% |
| Microstrip pair with coplanar ground, odd / even | 46.6778 / 58.7401 Ω | 46.6790 / 58.7413 Ω | 0.003% |
| Offset stripline with coplanar ground | 47.3538 Ω | 47.3546 Ω | 0.002% |

Effective εr agrees within 0.004% in every case. With the grid refined three times the field solver reproduces the FEM values to the last printed digit, so the remaining difference is the field solver's own grid error, and it is on the low side by 0.003 to 0.007%.

To rerun: `pip install scikit-fem scipy numpy`, install Gmsh's Python module, then `python3 scripts/reference/trace_fem.py`. It takes about a minute and a half.

What this does and does not show: both solvers solve the same mathematical problem (ideal rectangular copper, uniform dielectrics, static fields), and they agree to five digits. It says nothing about how well that problem describes a manufactured board. That needs measurements.

## Physics-based via model vs a 3D field solver

The via model is compared with [Palace](https://github.com/awslabs/palace), an open-source 3D finite-element solver. Geometries are in `fullwave/cases.json`, stored solver results in `fullwave/results.json`, and `npm run fullwave` reruns them (see `fullwave/README.md`). Both sides solve the same structure: solid planes, the via or pair, stitching vias, 50 Ω coaxial ports between barrel and plane, and an absorbing wall where the 3D domain ends.

The figure of merit is the largest complex difference |S_model − S_3D| over every compared S-parameter and every frequency from 1 to 50 GHz. S lies between 0 and 1, so 0.01 means 1% of full scale, or an error floor of −40 dB. For the differential pair all four entries of the first column of S are compared (by symmetry they determine every mixed-mode term).

| Case | Compared | Largest difference | Single-node model | 3D mesh uncertainty |
|---|---|---|---|---|
| `p2g4` One plane pair, four stitching vias | S11, S21, 1–50 GHz | 0.009 | 0.59 | 0.001 |
| `p2g2` One thick plane pair, two stitching vias | S11, S21, 1–50 GHz | 0.003 | 0.53 | 0.001 |
| `p2g1` One plane pair, a single stitching via on one side | S11, S21, 1–50 GHz | 0.008 | 0.37 | 0.001 |
| `p2g0` One plane pair, no stitching vias (the return current spreads into the planes) | S11, S21, 1–50 GHz | 0.004 | 0.32 | 0.001 |
| `p4g4` Four planes (three cavities), lossy dielectric, four stitching vias | S11, S21, 1–50 GHz | 0.030 | 1.08 | 0.003 |
| `p2d2` Differential pair, one plane pair, two stitching vias, lossy dielectric | Sdd: all four ports, 1–50 GHz | 0.010 | 0.47 | 0.001 |
| `p4s` Via with a stub: port at the third plane, open stub through a 1.2 mm plane pair | S11, S21, 1–50 GHz | 0.029 | 0.77 | 0.005 |

The single-node column is the classic form of this method (one node for the barrel, vias as line currents, static capacitances), which is what v0.2 used. `npm run diagnose` prints how much each part of the model contributes.

What this shows:

- For a single plane pair the model and the 3D solver agree to better than 0.01 up to 50 GHz, with four, two, one or no stitching vias.
- The differential pair agrees to 0.010, crosstalk terms included. Sdd21 at 50 GHz: −1.6 dB from the model, −1.8 dB from the 3D solver, −13 dB from the single-node model. The coupling of the two barrels' fringing fields is what brings it there (0.043 without it).
- The stub notch is at 23 to 24 GHz in both, −19 to −20 dB deep. The single-node model puts it near 20 GHz. Around the notch the difference peaks at 0.018 (a sharp resonance turns a small shift into a large difference in S), and stays below 0.01 from 29 to 40 GHz.
- For four planes 0.3 mm apart the difference is 0.006 up to 30 GHz, 0.013 up to 40 GHz and 0.030 at 50 GHz. Refining the 3D mesh moved it slightly away from the via model (0.027 on the coarse mesh), so this is a model error, not 3D noise. The likely cause is listed under assumptions in `docs/PHYSICS.md`: with thin planes close together, the fringing fields of neighbouring plane pairs reach each other through the antipad hole, and the model includes only the static part of that. The stub case, which has two 0.3 mm plane pairs, shows the same growth near 50 GHz (0.029), and 0.022 around its notch at 21 to 30 GHz.

Accuracy of the 3D results themselves: second-order curved elements, adaptive frequency sweep to a tolerance of 1e-4. Every case is stored on two meshes, coarse (16,000 to 70,000 tetrahedra) and fine (37,000 to 155,000). The comparison uses the fine mesh; the last column is the largest difference between the two meshes, at most 0.005, so the reference is good to well under the model differences it is used to measure. The coarse runs were computed on a 2-core Linux machine, the fine runs on an M4 MacBook Air in the container from `fullwave/Dockerfile` (Palace v0.19.0, Gmsh 4.12.1). Rerunning reproduces them to within the mesh uncertainty, not bit for bit.

The pass mark in the reference cases is 0.035. It was set after these comparisons were made, just above the worst case, as a regression guard: a change that makes any case worse than today's worst fails the tests. It was 0.03 until the fine-mesh runs put the four-plane case at 0.030; it was raised rather than the band cut. The accuracy claim is the table above, not the pass mark.

## The via with its pads and traces vs a 3D field solver

The cases above have their ports in the antipads. On a board a via is fed by a trace through a pad, and often ends inside the stack or runs on past its trace as a stub. Eighteen more 3D cases add that (`fullwave/cases.json`); the via model's launch is on for them, with its ports on the traces as in the 3D model (docs/PHYSICS.md).

| Case | What | Trace over the opening beyond the pad |
|---|---|---|
| `ms1` | Microstrip launch, 0.2 mm layer over the top plane, 0.38 mm trace, 0.5 mm pad, 0.7 mm antipad, one 0.8 mm plane pair | 0.1 mm |
| `ms1a`, `ms1b` | As `ms1` with a 0.5 mm and a 1.0 mm antipad | none, 0.25 mm |
| `ms2` | Microstrip launch on a thin layer: 0.1 mm, 0.19 mm trace, 0.45 mm pad, 0.6 mm antipad, 0.2 mm drill, 0.5 mm pair | 0.075 mm |
| `ms2a`, `ms2b` | As `ms2` with a 0.45 mm and a 0.85 mm antipad | none, 0.2 mm |
| `bl1` | Blind via: the barrel ends on a stripline 0.2 mm into the second plane pair (0.6 mm) | |
| `bl2` | Blind via into the first plane pair: stripline 0.15 mm below the top plane, 0.4 mm pair | |
| `bd1` | Backdrilled: as `bl1` with 0.2 mm of barrel left below the stripline | |
| `bl3` | As `bl2` with a 0.26 mm pad on the 0.2 mm drill: the trace runs 0.17 mm under the antipad opening | |
| `bl4` | As `bl2` with the four stitching vias twice as far out, at 1.6 mm | |
| `bl5` | As `bl2` with the stripline 0.25 mm deep instead of 0.15 (0.15 mm above the plane below) | |
| `bl6` | Thicker low-loss board, εr 3.4: the barrel crosses a 0.4 mm pair and ends on a 0.28 mm stripline 0.2 mm into the next (0.5 mm), 0.5 mm pad, 0.75 mm antipad | |
| `bd2` | Backdrilled in the same board: as `bl6` with 0.25 mm of barrel left below the stripline | |
| `th1` | Through via: the stripline of `bl1`, with the barrel running on through a third plane pair (0.3 mm) as a 0.8 mm stub | |
| `th2` | Through via with a long stub: the stripline of `bl2`, with the barrel running on through a 1.0 mm pair | |
| `th3` | Through via in the low-loss board: the stripline of `bl6`, with the barrel running on through a 0.8 mm pair | |
| `md1` | Differential pair, 0.8 mm pitch: a microstrip launch on each via (0.2 mm layer, 0.5 mm pads, 0.38 mm traces side by side), 0.6 mm antipads, one 0.8 mm plane pair | 0.05 mm |

**Ports and calibration.** The 3D ports are 50 Ω sheets at the far end of each trace. A stripline port is cut into strips up and down in the ratio of the stripline's return currents (a 50/50 split leaks 1.4% of the power into the plane pair at 10 GHz), with three fence vias around it. Each port still adds a little inductance and capacitance, and the 3D trace is not quite the trace solver's: a zero-thickness trace comes out 1 to 6% low in impedance and up to 3% high in εeff on these meshes, which is mesh error of the reference. So each case has two trace-only runs, at its trace length and at half of it. A fit of both gives each port as a small π box and the 3D trace's own impedance, εeff and leakage; the via model is then compared as port box · 3D trace · via. The box is a fit, not a circuit, so its two shunt capacitances may take either sign: the stripline ports want about −20 fF at the port and +40 fF on the trace side, and with the sign free the fit reproduces their thrus to 0.002 below 10 GHz instead of 0.011. The fit's residual up to 30 GHz is the last column below: the comparison's own uncertainty. For the pair the thrus are the two traces side by side, fitted and compared on their odd mode (Sdd).

**The 3D mesh where a trace leaves its pad.** Trace and pad are zero-thickness sheets, and the field at a sheet's edge is singular. On the coarse and fine meshes (second-order elements, 45 to 60 µm at the via) the 3D model gets too little of that edge field, and where the trace leaves the pad that shows as too much capacitance at the pad: a static (electrostatic) Palace solve of `bl1`'s barrel, pad and trace, taken apart into the pad, the line and the junction, puts the junction at about 7 fF on the coarse mesh and 11 to 13 fF as the mesh along the sheet edges is refined, still rising slowly. So there is a third mesh level, `edge`: the fine mesh with 8 µm elements along every trace and pad edge, which costs about 15% more. Every case with a trace has edge-mesh results and is compared on them. Up to each case's limit the fine and edge results differ by 0.011 to 0.030 (the column below), most of it this capacitance; between coarse and fine, which share the problem, the difference was 0.006 to 0.014.

Largest |ΔS| per band (`npm run diagnose` prints the same table, the Validation page shows it live):

| Case | 1–10 GHz | 11–20 | 21–30 | 31–40 | 41–50 | Mesh | Between meshes | Calibration |
|---|---|---|---|---|---|---|---|---|
| `ms1` | 0.002 | 0.004 | 0.013 | 0.062 | 0.139 | edge | 0.012 | 0.007 |
| `ms1a` | 0.007 | 0.017 | 0.027 | 0.041 | 0.044 | edge | 0.011 | 0.007 |
| `ms1b` | 0.014 | 0.031 | 0.063 | 0.141 | 0.352 | edge | 0.012 | 0.007 |
| `ms2` | 0.008 | 0.019 | 0.038 | 0.076 | 0.161 | edge | 0.021 | 0.002 |
| `ms2a` | 0.012 | 0.026 | 0.047 | 0.089 | 0.168 | edge | 0.020 | 0.002 |
| `ms2b` | 0.010 | 0.028 | 0.062 | 0.119 | 0.223 | edge | 0.018 | 0.002 |
| `bl1` | 0.004 | 0.021 | 0.087 | 0.187 | 0.234 | edge | 0.022 | 0.015 |
| `bl2` | 0.006 | 0.023 | 0.062 | 0.120 | 0.215 | edge | 0.032 | 0.005 |
| `bl3` | 0.019 | 0.032 | 0.034 | 0.038 | 0.120 | edge | 0.029 | 0.005 |
| `bl4` | 0.007 | 0.025 | 0.085 | 0.093 | 0.134 | edge | 0.031 | 0.005 |
| `bl5` | 0.010 | 0.032 | 0.084 | 0.159 | 0.380 | edge | 0.030 | 0.005 |
| `bd1` | 0.008 | 0.042 | 0.130 | 0.266 | 0.524 | edge | 0.021 | 0.015 |
| `bl6` | 0.009 | 0.024 | 0.054 | 0.122 | 0.229 | edge | 0.017 | 0.007 |
| `bd2` | 0.009 | 0.048 | 0.177 | 0.592 | 0.752 | edge | 0.014 | 0.007 |
| `th1` | 0.004 | 0.008 | 0.034 | 0.034 | 0.060 | edge | 0.017 | 0.015 |
| `th2` | 0.006 | 0.008 | 0.016 | 0.047 | 0.058 | edge | 0.020 | 0.005 |
| `th3` | 0.009 | 0.020 | 0.026 | 0.030 | 0.053 | edge | 0.011 | 0.007 |
| `md1` (Sdd) | 0.010 | 0.023 | 0.042 | 0.071 | 0.099 | edge | 0.011 | 0.007 |

What this shows:

- **The antipad opening was the largest part of the microstrip gap.** Without the opening term the microstrip launch differed from 3D by 0.036 at 1–10 GHz and 0.076 at 11–20 GHz on `ms2b` (0.030 and 0.063 on `ms1b`, fine mesh, v0.6); with it, about 0.01 and 0.03. The a/b variants were run to find out: shrinking the antipad to the pad took away most of the inductance the model was missing (18 pH on `ms2`, 4 pH on `ms2a`), and the 0.2 mm opening of `ms2b` raised it to 50 pH. One number in the opening term is fitted, the share of the slot inductance that counts across the pad, on `ms1`, `ms2`, `ms2a` and `ms2b`: 0.5 on the fine mesh (0.3 to 0.8 fitted about as well), 0.8 on the edge mesh (0.8 and 0.9 are within 0.001 of each other). Of the cases left out of the fit, two came closer to 3D with it, up to 20 GHz `ms1b` from 0.036 to 0.031 and `md1` from 0.026 to 0.023; `ms1a` went from 0.015 to 0.017 there and stays within 0.035 up to 36 GHz instead of 35.
- **Microstrip launch: checked up to 20 GHz.** All six cases stay within 0.035 up to 21 GHz (worst `ms1b`, 0.031 at 20 GHz), `ms1` up to 35 GHz. Above 30 GHz the open microstrip radiates, and the 3D loses more power than the model.
- **The pad–trace junction is in the model (0.8.1).** Where a stripline leaves its pad, the pad's disk and the trace's line share some field, and the model used to count it twice: 4 to 16 fF too much at the pad, the main error of every stripline launch below 15 GHz. `junction.ts` now solves that region in 3D in the browser (barrel, pad, planes with their antipads, the trace) and takes it off: 9.8 fF on `bl2`, 11.4 on `th2`, 13.5 on `bl1`, 15.3 on `bd1`, 16.4 on `th1`. Against the static 3D solve of `bl1` refined along the sheet edges (above) it is within about 1 fF; its own grid is within 0.2 fF of one twice as fine. The 1–10 GHz difference of the five stripline cases drops from up to 0.025 to 0.004–0.009.
- **Stripline on a through via: checked up to 30 GHz (13 before).** With the junction in, `th1` stays within 0.035 up to 45 GHz and `th2` up to 35. The stub's notch in S21 lands where 3D puts it: 31 against 32 GHz on `th1` (both about −27 dB), 23 GHz in both on `th2`. The held-out `th3`, the low-loss board's stripline with a 0.8 mm pair of stub below it, holds to 42 GHz (30 without the opening term).
- **Blind via ending on its stripline: checked up to 20 GHz. Backdrilled, with a stub left in the pair: 17 GHz.** One limit of 18 GHz covered both in 0.8. 0.9 adds five cases: three variants of `bl2`, `bl3` (a 0.26 mm pad, barely larger than the 0.2 mm drill), `bl4` (the stitching vias twice as far out) and `bl5` (the trace 0.25 mm deep instead of 0.15), and a thicker low-loss board, εr 3.4, as a blind (`bl6`) and a backdrilled (`bd2`) via. On the edge mesh the 0.8.2 model held them to 15, 21, 20, 22 and 17 GHz: `bl3` and `bd2` broke the 18 GHz claim. With the stripline's first stretch under the antipad opening in the model (next point), the blind cases hold `bl1` 22, `bl2` 23, `bl3` 39, `bl4` 22, `bl5` 20 and `bl6` 24 GHz; the backdrilled ones `bd1` 18 and `bd2` 17. `bl6` and `bd2` were run after the model was settled, as a held-out check; the opening term moved `bl6` from 22 to 24 and left `bd2` where it was. The stub is where the backdrilled cases part from 3D first; a series inductance for the stub's own current makes it slightly worse.
- **The stripline's first stretch under the antipad opening.** On `bl3` the trace runs 170 µm under the opening of the plane above before the stripline has both its planes, and its return there detours round the opening. `stripOpeningL` integrates the 2D solve of a stripline under a slot along the opening's chord, as the microstrip opening term does: 10.8 pH on `bl3`, 3 to 6 pH on the other cases, and it takes `bl3` from 15 to 39 GHz while the others move by a GHz either way. Across the pad itself nothing is added; the microstrip's share of 0.8 there cost the other cases 1 to 2 GHz.
- **What a blind end misses above 20 GHz.** Seen from its top port, the capacitance at a blind end rises with frequency, and the model misses part of the rise. A body-of-revolution full-wave solve (`scripts/reference/via_bor_fem.py`, second-order elements in r and z, which reproduces the 3D case `p2g0` within 0.008 up to 50 GHz) of `bl2`'s barrel and pad on their own, round, without trace or stitching vias: static 105.05 fF in both, but 116.5 fF at 20 GHz against the model's 110.0. A Palace run of the same round geometry (a research variant of `via_case.py` without the trace and its port) rises the same way, 1.6 fF higher throughout (the pad's sheet edge on the edge mesh). The charge that lands on the plane below has to come round through the plane pair, so all of it drives the pair's parallel-plate wave: 45% of the via's charge (47 fF), where the depth rule (`PHYSICS.md`) drives it with 27 fF. Stamping each part on the plane it lands on matches the round solve within 0.4 fF and the round 3D run within 0.8 fF up to 30 GHz, and a 3D run with `bl2`'s four stitching vias, still without the trace, within 0.6 fF up to 25 GHz.
  With a trace on the pad the same change does not hold. Below 20 GHz it makes the cases worse (`bl2` 19 GHz instead of 23), and a run of `bl2` without its stitching vias, trace kept, does the same (17 instead of 19, fine mesh): the trace, not the stitching vias, sets that. Attaching the stripline's own return to the plane below at the pad edge instead of at the barrel brings the cases back (`bd1` 20, the others 20 to 27 GHz), but above 30 GHz the end then couples far too strongly to the plane pair's resonances between the stitching vias: band errors of 0.3 to 1.3 above 30 GHz where the depth rule gives 0.1 to 0.5. The round 3D run with stitching vias shows that resonance too, much weaker than this split gives it. So the depth rule stays until the end's coupling is understood with a trace on it; the round references and these runs are the starting point. `launch.test.ts` holds the model to the round solve's static value (0.3%) and to its rise within 7% at 20 GHz.
- **Differential pair, microstrip launch: checked up to 20 GHz.** Compared on Sdd on the edge mesh, `md1` stays within 0.023 up to 20 GHz and 0.035 up to 27 GHz. One case so far; a pair's stripline launch has no case yet, so there the single-via stripline cases are the check.
- These limits are `F_LAUNCH_CHECKED` in the via model: above them the Via page shades the sweep as indicative, and `launch.test.ts` fails if any case leaves 0.035 below its limit.

These cases are measured, not pass/fail: they are reported band by band instead of against a single mark.

**Past the limits.** A limit says where every case of a launch stays within 0.035; it does not say how far off the model is above it. So the Via page quotes that too: above a launch's limit, the largest |ΔS| over that launch's cases up to each of 20, 25, 30, 40 and 50 GHz above the limit, for the via's own launches (the worse one at each frequency when it has two), in the notes, the spec table and the Touchstone header. The numbers (bold: above the limit, what is quoted; rounded up):

| Launch | Cases | Checked to | ≤ 20 GHz | ≤ 25 GHz | ≤ 30 GHz | ≤ 40 GHz | ≤ 50 GHz |
|---|---|---|---:|---:|---:|---:|---:|
| Microstrip | ms1, ms2, ms1a, ms1b, ms2a, ms2b | 20 GHz | 0.031 | **0.045** | **0.064** | **0.15** | **0.36** |
| Blind, ends on its stripline | bl1, bl2, bl3, bl4, bl5, bl6 | 20 GHz | 0.033 | **0.054** | **0.087** | **0.19** | **0.38** |
| Backdrilled, stub left in the pair | bd1, bd2 | 17 GHz | **0.048** | **0.096** | **0.18** | **0.60** | **0.76** |
| Stripline on a through via | th1, th2, th3 | 30 GHz | 0.021 | 0.022 | 0.034 | **0.047** | **0.061** |
| Pair, microstrip launch | md1 | 20 GHz | 0.023 | **0.031** | **0.042** | **0.072** | **0.10** |

`npm run launch-error` writes them to `src/physics/launchError.ts` from the stored 3D results, and `launch.test.ts` fails if they no longer match. Each case has a launch on one end only; on a via with a launch on both ends the two ends' differences can add, and the Via page says so.

## What is not validated yet

- The launch is checked up to 20 GHz (microstrip), 20 GHz (a blind via that ends on its stripline), 17 GHz (a backdrilled one with a stub left in the pair), 30 GHz (a stripline on a through via) and 20 GHz (a pair's microstrip launch, one case). Above those a blind or backdrilled end parts from 3D (above). A pair's stripline launch has no 3D case yet; shared (oval) antipads of a pair are not modelled.
- Trace loss has no independent check of the complete number: its parts are checked against exact limits, but the conductor-loss rule (Wheeler) is itself an approximation that needs a skin depth well below the copper thickness, and roughness models are fits. Expect a few percent on smooth copper and more on rough copper.
- Coplanar ground is ideal copper tied to the plane. Real ground pours with sparse stitching vias are not modelled.
- Nothing has been compared with a measured board: no trace impedance, no trace loss, no via.
- Above 50 GHz there is no comparison; the app shades that part of every via sweep.

If you can share measured or simulated S-parameters of traces or vias with full geometry, please open a **Reference data** issue. Each data set becomes a reference case.

## Adding a reference case

```ts
{
  id: 'ms.myCase', group: 'Microstrip', what: 'short description of inputs → expected',
  source: 'where the expected value comes from',
  run: () => ({ got: microstrip(...).z0, expected: 50.0, unit: 'Ω', tolPct: 1 }),
},
```

Use `tolPct` for relative tolerance or `tolAbs` for absolute. Run `npm run diagnose` to see it.
