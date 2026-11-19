> Copied from Impedora's repository, where it is docs/PHYSICS.md. File paths under `src/` and `scripts/` refer to that repository; the 3D cases and results are in `fullwave/` here.

# Physics models

All code is in `src/physics/`. Units are SI (metres, ohms, hertz) everywhere inside the physics code; conversion to mm/mil happens only in the UI.

Accuracy is enforced by the reference cases in `src/physics/reference.ts`. Run `npm run diagnose` to see current errors.

## Trace impedance: 2D field solver
Files: `fieldSolver2d.ts` → `solveTrace(geometry)`, `traceModel.ts` → `solveTraceStack(input)`

This is what the Trace page shows. The closed-form formulas further down are kept for comparison and to speed up the width search.

**Problem solved.** The cross-section of the trace is solved as an electrostatic problem, ∇·(ε∇V) = 0, with the trace at 1 V and the planes at 0 V: copper of real thickness, rectangular or etched to a trapezoid, every dielectric layer between the trace and its planes with its own εr, one or two reference planes, optional solder mask (one thickness on the copper and its side walls, another over the laminate), optional ground copper on the trace's own layer (coplanar waveguide with ground), a single trace or an edge-coupled pair. Two solves per mode, one with the dielectrics and one with air everywhere, give the capacitance per length C and C_air:

    Z0 = 1 / (c·√(C·C_air)),   εeff = C / C_air,   L = μ0·ε0 / C_air

For a pair, the odd mode has an electric wall between the traces and the even mode a magnetic wall. Zdiff = 2·Zodd, Zcommon = Zeven/2.

**Method.**
- Finite volumes on a rectangular grid with a grid line on every edge. On this grid the scheme is identical to linear finite elements on right triangles, so the capacitance (from the field energy) is always slightly too high and falls towards the true value as the grid is refined.
- The field is singular at the trace corners. The first cell there is 1/4096 of the smallest dimension (1/65536 for a zero-thickness trace, whose edge is a stronger singularity), and cells grow by a factor 1.5 with distance. Because the growth is geometric, a tiny first cell and a very large box cost only a few more grid lines.
- Open structures (microstrip) are closed by a box 200 times the size of the structure. 100 and 1600 give the same result to 1 part in 100,000.
- Each mode is solved on two grids, the second with every cell halved. The energy error falls four-fold, so (4·fine − coarse)/3 removes the leading error (Richardson extrapolation). Alone, the fine grid is about 0.3% off; extrapolated, the result is within 0.01% of the exact and independent references.
- Etched traces: the grid columns lean with the copper walls. Column i is displaced sideways by lean(x)·τ(y), where τ goes from 0 at the wide face to 1 at the narrow face, so rows stay horizontal and only the cells beside the copper become slanted. The elements there are general triangles (a 7-point stencil instead of 5). The preconditioner is built on the unsheared grid, which keeps it an M-matrix.
- Coplanar ground: copper at 0 V on the trace's layer from `gap` outward on both sides, with the same thickness and etched walls, and solder mask over it when the mask is on.
- Linear system: conjugate gradients with a modified incomplete Cholesky preconditioner, about 40 iterations per solve. A single trace takes roughly 0.1 s, a pair 0.3 s. In the app the solve runs in a Web Worker (`traceWorker.ts`, jobs in `traceModel.ts` → `runTraceJob`), so the page stays responsive; a new edit replaces the job in flight.

**Where the stackup goes** (`store.svelte.ts` → `traceStack`). From the signal layer, walk to the nearest reference plane on each side. Each dielectric layer is kept with its own thickness, εr and loss tangent. Two planes: stripline. One plane: microstrip, with any dielectric on the other side of the trace kept (embedded microstrip). Non-plane copper crossed on the way is treated as dielectric of the next layer and flagged. On an inner layer the resin beside the trace takes the average εr of the two adjacent dielectrics.

**Assumptions and limits.**
- Quasi-static (TEM). Valid while the cross-section is small against the wavelength; for PCB traces that holds far beyond 50 GHz.
- Etched shape is a straight-walled trapezoid given by the width lost on the narrow face. Real etch profiles are slightly curved; the difference is small against the uncertainty of the etch itself.
- Coplanar ground is ideal: continuous copper at the potential of the reference plane. That needs stitching vias close to the gap, spaced well under a tenth of a wavelength. A pour with sparse vias behaves differently and can resonate.
- Planes are wide (no slots or edges near the trace) and other traces are absent, apart from the partner of a pair.
- The result is only as good as the inputs. Finished trace width, plated copper thickness, pressed dielectric thickness and the εr of that resin content usually matter far more than the numerical error.

**Checked** (`reference.ts`, group "Trace solver"): exact zero-thickness stripline and edge-coupled stripline (within 0.01%), Hammerstad–Jensen microstrip in air (0.01%), sixteen thick-trace geometries against an independent second-order finite-element solver (within 0.007%; they include etched traces with either face narrow, solder mask with one and with two thicknesses, coplanar ground, pairs), exact zero-thickness coplanar waveguide (unit test), grid convergence, mirror symmetry. See VALIDATION.md.

## Trace loss
File: `traceLoss.ts` → `lineAt(mode, setup, f)`

R, L, G and C per length at each frequency, then γ = √((R + jωL)(G + jωC)) and Zc = √((R + jωL)/(G + jωC)). Insertion loss of a matched line is 8.686·Re(γ)·length dB. For a pair the odd mode is used, which gives the differential insertion loss.

**Dielectric (C, G).** The solver returns the share g of the field energy in each material, so C = ε0·Σ εr·g. With complex, frequency-dependent permittivities the same sum gives C(f) and G(f) without another field solve; because the energy is stationary this is exact to first order in the change of ε. Each dielectric follows the wideband Debye model (Djordjevic, Biljic, Likar-Smiljanic, Sarkar, IEEE Trans. EMC 43, 2001):

    ε(f) = ε∞ + a·ln((f2 + jf)/(f1 + jf)),   f1 = 1 kHz, f2 = 1 THz

with ε∞ and a chosen so that ε equals the entered εr·(1 − j·tanδ) at the frequency you enter. Loss tangent is then nearly flat and εr falls by 1.47·εr·tanδ per decade, which is what causality requires.

**Conductor (R, L).** Wheeler's incremental-inductance rule: R_ac = (Rs/μ0)·∂L/∂n, Rs = √(π·f·μ0·ρ), where ∂L/∂n is the change of inductance when every conductor surface recedes. The solver gets it from the same field solution by re-evaluating the energy on a grid whose conductor edges are moved by 0.01% (a shape derivative, no extra solve). It includes current crowding at the trace edges, the proximity effect in a pair and the loss in the reference planes. The skin layer also adds internal inductance R_ac/ω, which raises the impedance slightly at low GHz.

- Low frequency: the rule needs a skin depth well below the copper thickness. Below that, R is blended into the DC resistance of the trace, R = √(R_dc² + R_ac²). The page warns when the skin depth exceeds a third of the thickness.
- Roughness multiplies R_ac. Cannonball–Huray (Huray 2010; Simonovich 2015): K = 1 + (7π/3) / (1 + δ/r + δ²/2r²), r = 0.06·Rz. Hammerstad's K = 1 + (2/π)·atan(1.4·(Rq/δ)²) is in the code as well. The same roughness is applied to every surface.
- Copper at 20 °C, ρ = 1.72·10⁻⁸ Ω·m.

**Not included:** plating finish (ENIG nickel on outer layers adds loss), different roughness on the two sides of the foil, glass-weave effects, radiation, temperature.

**Checked** (group "Trace loss"): dielectric attenuation of a TEM line against π·f·√εr·tanδ/c (0.005%), the shape derivative against really shrinking every conductor by 50 nm, for rectangular and for etched copper with coplanar ground (0.006%), the wide-stripline limit Rs/w (0.14%), the Kramers–Kronig slope of the Debye model, the roughness formula.

### Touchstone export for a trace
File: `traceSparams.ts` → `traceSweep(...)`

A length of line as a two-port: A = D = cosh(γℓ), B = Zc·sinh(γℓ), C = sinh(γℓ)/Zc, converted to S at 50 Ω. A pair is a four-port from its odd and even modes (each solved and given loss the same way): S11 = (S11e + S11o)/2, S12 = (S11e − S11o)/2, S13 = (S21e + S21o)/2, S14 = (S21e − S21o)/2, ports 1, 2 at the near end and 3, 4 at the far end. Mode conversion is zero because the pair is symmetric. Checked: reciprocity, passivity, Sdd21 equals the odd-mode line, Touchstone round trip.

## Trace impedance: closed-form formulas

Shown next to the field-solver result for comparison.

### Microstrip: Hammerstad & Jensen (1980)
File: `microstrip.ts` → `microstrip(w, h, t, er)`

Quasi-static Z₀ and εeff with the HJ finite-thickness correction.

- Source: E. Hammerstad, Ø. Jensen, "Accurate Models for Microstrip Computer-Aided Design", IEEE MTT-S Digest, 1980, pp. 407–409.
- Valid: 0.01 ≤ w/h ≤ 100, εr ≤ 128. The paper quotes ≤ 0.2% for εeff (and better for the air impedance) in that range, zero thickness.
- Not modelled: solder mask, trapezoidal etch, dispersion, Dk vs frequency.
- Checked: Pozar Example 3.7 (εr 2.2, w/h 3.08 → 50 Ω, εeff 1.87): matches within 0.1% (Z₀).

### Stripline, zero thickness: exact
File: `stripline.ts` → `striplineExact(w, b, er)`

Z₀ = η0/(4√εr) · K(k)/K(k′), k = sech(πw/2b). Exact by conformal mapping (Cohn 1954). η0/4 = 94.18 Ω; the textbook form 30π = 94.25 Ω rounds η0 to 120π and is 0.07% high.

### Stripline, finite thickness: Wheeler (1978)
File: `stripline.ts` → `stripline(w, b, t, er)`

- Source: H. A. Wheeler, "Transmission-Line Properties of a Strip Line Between Parallel Planes", IEEE Trans. MTT-26, 1978.
- Checked: as t → 0 it matches the exact solution within 0.3% for w/b 0.2–1.0.
- Assumes centred trace and uniform dielectric. For off-centre traces and mixed εr the field solver is the answer; this formula then only serves as a rough comparison (thickness-weighted εr).

### Edge-coupled stripline: Cohn (1955), exact for t = 0
File: `stripline.ts` → `coupledStripline(w, s, b, er)`

Even/odd-mode impedances via elliptic integrals; Zdiff = 2·Zodd, Zcm = Zeven/2.
Copper thickness not included (typically lowers Z 5–15% in thin cores).

### Edge-coupled microstrip: IPC-2141 coupling (approximate)
File: `microstrip.ts` → `coupledMicrostripIPC(...)`

Zdiff = 2·Z₀·(1 − 0.48·e^(−0.96·s/h)), Z₀ from Hammerstad-Jensen. Empirical; several % error, worse for tight coupling and with solder mask.

## Via: physics-based model
Files: `viaModel.ts`, `viaModes.ts`, `axisym.ts`, `bessel.ts`

The board is cut at every reference plane in the via's return path. Each pair of neighbouring planes is a cavity. Three effects are modelled in each cavity.

**1. Plate impedance and proximity** (`plateScattering.ts`). The via current drives the parallel-plate wave. For thin vias, each one is a line current and they couple through the impedance of an infinite parallel-plate waveguide:

  Z_ab = (ω·μ0·h / 4) · H0⁽²⁾(k·ρ_ab),  k = ω·√(μ0·ε0·εr·(1 − j·tanδ))

where ρ_ab is the distance between vias a and b, and ρ_aa is the via's own radius. Stitching (ground) vias short the cavity at their positions, which reduces the signal-via impedance to Z_ss − Z_sg·Z_gg⁻¹·Z_gs.

Real barrels are conducting cylinders, and current crowds on the side that faces a neighbour. Impedora therefore expands the plate wave around every via in cylindrical harmonics e^{jmφ}, |m| ≤ 2, moves them between vias with Graf's addition theorem, and sets every harmonic of the total field to zero on each barrel. Harmonic 0 alone is the line-current formula above. The error of stopping at harmonic m, for the loop inductance of two equal barrels against the exact (μ0·h/π)·acosh(s/2a):

| Spacing / radius | m = 0 (line currents) | m ≤ 1 | m ≤ 2 |
|---|---|---|---|
| 8 | 0.8% | 0.006% | < 0.001% |
| 6 | 1.6% | 0.03% | < 0.001% |
| 4 | 5.3% | 0.2% | 0.01% |
| 3 | 14% | 1.2% | 0.1% |
| 2.5 | 32% | 4.7% | 0.9% |

**2. Barrel-to-plane capacitance, static** (`axisym.ts`). One finite-volume electrostatic solve in (r, z) for the whole barrel: barrel and its pads at 1 V, every plane it passes (with its antipad) at 0 V, the stackup's dielectrics in between. Non-uniform grid that conforms to every edge, preconditioned conjugate gradients. The charge each plane collects is split at the plane's mid-thickness into the part that arrives from the cavity above and the part from the cavity below. Solving the whole stack at once matters when planes are close together: the fringing fields of two antipads then overlap and a plane-by-plane solve double counts.

**3. Radial modes** (`viaModes.ts`). The fringing field around an antipad is a sum of evanescent plate modes, cos(nπz/h) across the cavity, n ≥ 1. They are not static: each one decays radially as exp(−κ_n·ρ) with κ_n = √((nπ/h)² − k²), so the capacitance rises with frequency, and they link the antipad at the top of the cavity to the one at the bottom. The older single-node barrel model ignores both. Impedora solves the modes for a barrel of radius a through antipads of radius b:

- Each antipad is a coaxial aperture with the TEM field E_r = V/(r·ln(b/a)).
- Between the plates, H_φ(r, z) = Σ g_n(r)·cos(nπz/h). Each g_n obeys Bessel's equation of order 1 with k_n² = k² − (nπ/h)², driven by the aperture field, with E_z = 0 on the barrel and an outgoing (or decaying) wave for r > b.
- The magnetic field is projected back on the aperture's TEM mode (Galerkin), which keeps the result reciprocal.

The admittance seen at one aperture from mode n is

  y_n = −j·2π·ω·ε·ε_n / (h·ln(b/a)·k_n²) · [1 − T_n·H0⁽²⁾(k_n·b)/H0⁽²⁾(k_n·a)]
  T_n = π·[Y0(k_n·b)·J0(k_n·a) − J0(k_n·b)·Y0(k_n·a)] / (2·ln(b/a)),  ε_0 = 1, ε_n = 2

and (−1)ⁿ·y_n between the two apertures. For n ≥ 1 this is a capacitance, y_n = jω·C_n:

  C_n = 4π·ε / (h·ln(b/a)·κ_n²) · {1 − [I0(κ_n·b)·K0(κ_n·b) − I0(κ_n·a)·K0(κ_n·b)²/K0(κ_n·a)] / ln(b/a)}

As a Π network with the barrel as the common node:

| Branch | Value | Meaning |
|---|---|---|
| Barrel to each plate | 2·Σ C_n over odd n | barrel-to-plane capacitance, frequency dependent |
| Plate to plate | y_0 + jω·Σ (−1)ⁿ·C_n over n ≥ 1 | plate impedance in parallel with a negative capacitance |

At low frequency y_0 → 1/Z_pp, the plate impedance of (1). The negative capacitance is what makes a via in a thin cavity behave like a short transmission line: the effective inductance falls with frequency while the capacitance rises.

How the pieces are combined: the static capacitance comes from the field solver (2), because it handles plane thickness and pads. The mode sums add the frequency-dependent change, 2·Σ [C_n(f) − C_n(0)], and the plate-to-plate capacitance. The first 40 modes carry the frequency dependence; 400 are summed for the static limit. In the zero-thickness limit the two independent methods agree on the static capacitance within about 1% (reference case `modes.static`).

Coupling to other vias. Let R_n be the reaction matrix of (1) at wavenumber k_n: the wave of mode n that arrives at each signal via per unit wave it sends out, after scattering off every other barrel. Mode n's admittance matrix over the signal vias is then

  Y_n = y_n·I − c_n · R_n · (H0⁽²⁾(k_n·a)·I + J0(k_n·a)·R_n)⁻¹,  c_n = 4ω·ε·ε_n·T_n² / (k_n²·h·H0⁽²⁾(k_n·a))

For n = 0 this reduces to (Z_ss − Z_sg·Z_gg⁻¹·Z_gs)⁻¹ at low frequency. For n ≥ 1 the modes are evanescent and R_n falls off as exp(−κ_n·d); it is included for every mode that has not decayed by e⁻⁷ at the nearest other via. In thick plane pairs this is what couples the two barrels of a differential pair through their fringing fields.

**Plane thickness.** The antipad hole through a plane of thickness t is a short coaxial section. Its capacitance is part of (2); its inductance μ0·t·ln(b/a)/(2π) sits in series between the two faces of the plane.

**Loss.** Dielectric loss enters through the complex k and through a factor (1 − j·tanδ) on every capacitance. Conductor loss is not modelled.

Ports sit between the barrel and the plane nearest to each trace layer, on the face towards that layer. The circuit is solved by nodal analysis at each frequency; differential pairs give mixed-mode Sdd = ½(S_aa − S_ab − S_ba + S_bb).

Method: R. Rimolo-Donadio et al., "Physics-Based Via and Trace Models for Efficient Link Simulation on Multilayer Structures Up to 40 GHz", IEEE Trans. MTT 57(8), 2009, pp. 2072–2083. Y.-J. Zhang, J. Fan, "An Intrinsic Circuit Model for Multiple Vias in an Irregular Plate Pair Through Rigorous Electromagnetic Analysis", IEEE Trans. MTT 58(8), 2010, pp. 2251–2265. The mode solution above follows the second paper's idea and is derived independently in `viaModes.ts`.

**Bessel functions** (`bessel.ts`): power series for |z| ≤ 12, Hankel asymptotic expansion above; complex argument, so dielectric loss is exact. Integer orders by Miller's backward recurrence (J) and forward recurrence (H⁽²⁾). I0 and K0 come from J0 and H0⁽²⁾ at imaginary argument, with a scaled asymptotic form for large arguments.

Assumptions (also shown in the app):
- Planes are infinite: no board-edge cavity resonances. (An absorbing circular wall can be switched on for comparison with a 3D solver, `edge` in `PhysViaInput`.)
- All reference planes are tied together by the stitching vias.
- Validity limit: 50 GHz (the highest frequency compared with a 3D solver) or 90% of the cut-off of the first higher-order plate mode in the thickest plane pair, f = c/(2h√εr), whichever is lower. Above it the plot is shaded as indicative.
- Antipads must not overlap: vias closer than an antipad diameter are outside the model.
- The antipad carries only the coaxial TEM field. With thin planes close together, the fringing fields of neighbouring cavities also reach each other through the antipad hole. The static part of that is in (2); its frequency dependence is not modelled. This is the largest known error: 0.030 in S at 50 GHz for four planes 0.3 mm apart (0.013 up to 40 GHz) (see docs/VALIDATION.md).
- Pads on unconnected inner layers are assumed removed. Conductor loss and copper roughness of the via are not modelled. The pads and traces at the two ends are the launch, next section.
- Without the launch, a stripline port is referenced to its nearest plane only.
- When the barrel ends inside a cavity (backdrill, blind via), that cavity falls back to the single-node model with the barrel length inside it. This is approximate, and the validity limit there is a quarter wavelength of barrel. With the launch on, such a barrel ends on a stripline and the cavity is handled as described below.

Touchstone export (`touchstone.ts`): Touchstone 1.1, `# Hz S RI R <z0>`. 2-port order S11 S21 S12 S22. 4-port (differential) rows S11…S14, S21…; ports 1/2 = via A/B on the first layer, 3/4 = via A/B on the second layer. Single-ended data; convert to mixed mode in your tool if needed.

## Via: the launch (pads and traces)
Files: `viaModel.ts` (`launch`), `channel.ts`

Without the launch the via's ports sit in the antipads, between the barrel and a plane. A real via is fed by a trace through a pad. With **Include pads and the traces** on, each port moves onto its trace at the edge of the pad, the pad and the short sections around it become part of the via, and the traces (any length, from the trace solver, with their loss) are cascaded on: trace · via · trace as one two-port, exported as one Touchstone file.

**Microstrip end** (outer layer over a plane):

| Part | Model |
|---|---|
| Trace up to the pad edge | Lossy line from the 2D field solver and the loss model (`traceLoss.ts`): A = D = cosh γℓ, B = Zc·sinh γℓ, C = sinh γℓ / Zc |
| Path across the pad, from its edge to the barrel | The trace's own inductance per length L′ = Z0·√εeff / c over rp − a, split as a T around the pad node. Its capacitance is the pad's, counted below |
| Pad | Its share of the barrel-to-plane capacitance: the axisymmetric solve with the pad minus the solve without it, at the pad node. Counted once |
| Barrel from the pad down to the plane | Half the partial inductance of a post twice as long (the plane's image): ½·Lp(2h, a), Lp(l, r) = (μ0·l/2π)·[asinh(l/r) − √(1 + (r/l)²) + r/l] |
| Pad and trace over the antipad opening | See below |

From the plane on, the via is the model above.

**The antipad opening.** Between the barrel and the edge of the plane's clearance there is no plane under the pad, nor under the trace where the clearance is larger than the pad. There the return current cannot run under the conductor, and the field reaches down through the hole to the next plane. At a distance x from the via axis, the opening's chord under the trace is a slot of half-width s = √(ra² − x²) in the plane. A 2D solve of the trace over a plane with such a slot, and the next plane of the stack below (`slotInductance` in `fieldSolver2d.ts`, the same finite-volume method as the trace solver), gives its inductance per length L′(s), and the stretch adds

ΔL = ∫ [L′(s(x)) − L′(0)] dx

from the pad edge out to the plane's edge (in series with the trace), plus a share of the same integral from the barrel to the pad edge (in series with the crossing). Across the pad the current spreads over the pad's full width instead of the trace's, and the 3D cases put the share that counts at 0.8 on the edge-refined mesh (0.5 on the fine mesh, where 0.3 to 0.8 fitted about equally well). L′(s) is solved at three slot widths and interpolated (monotone cubic). For case `ms2b` (0.1 mm layer, 0.2 mm of opening beyond the pad) this is 56 pH, for `ms1b` (0.2 mm layer, 0.25 mm) 38 pH: the largest single correction of the launch, and the one that made the 3D cases agree (docs/VALIDATION.md). An independent second-order FEM solve of the same slot problem (`scripts/reference/slot_fem.py`) agrees with the solver within 0.05%.

**Stripline end** (an inner trace between two planes). The barrel either ends on it (a blind or backdrilled via) or runs on through the plane beyond it (a through via: the rest of the barrel is a stub). Either end of the via can be a stripline; the upper one is the mirror image of the lower one, and the tests check that the two agree.

- A stripline at depth d below its upper plane, in a plane pair of height h, returns (h − d)/h of its current through the upper plane and d/h through the lower one. The port's voltage is measured against the same weighted mean of the two planes, which is what keeps it orthogonal to the plate wave.
- At the via, the part of the return current that has to change planes crosses the plane pair as displacement current and drives the parallel-plate wave. Where the barrel ends on the trace, this cavity carries the full-height plate impedance between the two planes (with the stitching vias), not one scaled by the barrel length. Where the barrel runs on, the cavity is the radial-mode model of a barrel through the whole pair, and the current the barrel carries past the trace (into the stub) crosses the full height on its own; the network adds the two.
- Where the trace leaves the barrel, the barrel's current steps, and the step stirs the higher-order plate modes. Their inductance is Σ_{n≥1} (μ/πh)·S_n²·I0(κ_n·a)·K0(κ_n·a), S_n = (h/nπ)·sin(nπd/h), κ_n = nπ/h. Only the trace's current makes the step (a current through the whole height excites no higher-order mode), so it sits in the trace's branch, the same whether the barrel ends there or runs on.
- Where the barrel ends, its end and any stub below have an end capacitance, from the axisymmetric solve with the plane beyond closed under the barrel. The stripline's pad is solved on its own: the barrel's static solve without the pad, then with it, and the difference is the pad's. A barrel that runs on has no end there: its pad is part of the barrel's capacitance, and the stub is the rest of the via model.
- Depth rule (v0.8). Inside a plane pair every capacitance is fed by conduction current along the barrel or the trace, and its field lands on both planes: a charge at depth z below the upper plane, in a pair of height h, ends (h − z)/h on the upper plane and z/h on the lower one. So the end capacitance (at the end's depth) and the pad (at the trace's depth d) are stamped against that weighted mean of the two plate nodes, not against the plane beyond. Before v0.8 they were, which put the end resonance some 25% too low in frequency. The share is exact for a charge between two whole planes. At a blind end it is low: the antipad of the upper plane sits right above the pad, so more of the pad's and the end's field reaches the lower plane, and that charge is fed through the plane pair's own inductance. The static solve puts 45% of `bl2`'s 105 fF on the plane below, 47 fF, where the rule drives the plate wave with 27. A round full-wave solve (`scripts/reference/via_bor_fem.py`) and a Palace 3D run of `bl2`'s barrel and pad, round, no trace, both show the missing part as capacitance that rises with frequency (6.5 fF of 116 at 20 GHz), and stamping each part on the plane it lands on matches both within 0.4 and 0.8 fF up to 30 GHz. With a trace on the pad the same change makes the 3D cases worse below 20 GHz, and above 30 GHz it couples the end far too strongly to the plane pair's resonances between the stitching vias; it waits until that is understood (VALIDATION.md).
- The junction (0.8.1). Where the trace leaves the pad, the pad's capacitance (from the axisymmetric solve, a whole disk) and the trace's (the line model from the pad edge on, full C′ from the start) share part of the same field: the pad's fringe towards the trace, the trace's own field next to the pad, and the first stretch of trace under the antipad of the plane above, where one of its planes is missing. That shared part was counted twice. `junction.ts` solves the region in 3D: a finite-difference Laplace solve on a graded Cartesian grid over half the geometry (the trace's centre plane is a symmetry plane), barrel and pad at 1 V once without the trace and once with it, the trace running out of the grid as a line; the junction is C(pad + trace) − C(pad) − C′·(trace length in the grid), with C′ from a 2D solve on the same grid so the grid's own error cancels. It is negative, 10 to 16 fF on the 3D cases, and is added to the pad's capacitance. Planes are slabs with their thickness and antipads; the grid reaches one pair height past the near planes. A tenth of a second to about one second per geometry, cached; the same junction upside down is solved once, so the two ends of a flipped stack agree exactly. Checked against a static 3D solve refined along the sheet edges (VALIDATION.md) and against its own grid refined twice.
- The path across the pad is the same T as above, with the stripline's L′.
- The stripline's first stretch under the antipad opening (0.9). From the pad edge out to the antipad edge, the plane the barrel passed through has its opening over the trace: at distance x from the axis a slot of half-width √(ra² − x²). Its return there detours round the opening. The extra inductance is ∫ [L′(s(x)) − L′(0)] dx from the pad edge to the antipad edge, with L′(s) from the 2D solve of a stripline under a slot (`slotPlanesC`, both planes of the pair slotted where the barrel passes through both, closed by the next plane beyond or 20 pair heights away), in series with the trace. 10.8 pH on a 0.26 mm pad under a 0.6 mm antipad (`bl3`), 3 to 6 pH on the usual pads. Across the pad nothing is added, unlike the microstrip pad over its opening: the current spreads over the whole disc there and returns on the plane on the pad's other side, and the 3D cases agree best without it. Its capacitance is in the junction.

**Differential pair.** Each via of the pair gets its own launch, as above. The two barrels between pad and plane of a microstrip end couple through their images: mutual inductance ½·Mp(2h, p), with Mp the partial mutual inductance of two parallel wires at the pitch p (the formula for Lp with r → p). In the odd mode that takes it off each via's inductance. The traces' coupling across the pads, the pads' mutual capacitance and the higher-order plate modes of one via reaching the other are left out: at practical pitches they are small next to what is kept (the plate modes decay as K0(nπp/h)). On the Via page a pair's ports sit on the traces at the pad edge; attaching a length of coupled traces comes later.

**Not modelled yet** (the app says so where it applies):

- Above 20 GHz a barrel that ends on its stripline parts from 3D at its end (the depth rule's share, above); with a backdrill stub left in the pair above 17 GHz (`bd1`, `bd2`).
- Shared (oval) antipads of a pair.
- Radiation from the open microstrip launch above about 30 GHz.

The validity limit includes the launch: above the frequency up to which the 3D checks agree within about 0.035 (`F_LAUNCH_CHECKED`), the sweep is shaded as indicative, and the Via page quotes how far off the launch's 3D cases are there, up to 20, 25, 30, 40 and 50 GHz (`launchError.ts`, docs/VALIDATION.md).

## Via: Johnson & Graham lumped model (quick estimate)
File: `via.ts`

- L = 5.08·h·[ln(4h/d) + 1] nH, C = 1.41·εr·T·D₁/(D₂ − D₁) pF (inches; h via length, d drill, T board thickness, D₁ pad, D₂ antipad).
- Barrel impedance √(L/C) uses the full board thickness (where C is defined). The signal section gets L(its length) and C × (its length / T). Stubs are open-circuited lines of the barrel impedance, velocity c/√εr. This avoids double-counting stub capacitance (regression case `via.nodouble`).
- Stub resonance: f = c / (4·ℓ·√εr).
- Source: H. Johnson, M. Graham, *High-Speed Digital Design*, 1993, ch. 7.
- Limits: first-order. Valid while the via is electrically short. No plane-pair cavity modes, no stitching vias, single-ended only. The physics-based via model (v0.2) replaces this.

## Current: IPC-2221
File: `current.ts`

I = k·ΔT^0.44·A^0.725 (A in mil², k = 0.048 external, 0.024 internal).
Superseded by IPC-2152; internal-layer values are known to be conservative. Resistance: R = ρL/(w·t), ρ = 1.72×10⁻⁸ Ω·m at 20 °C, α = 0.00393 /°C.

## KiCad stackup import
File: `src/import/kicad.ts`

Reads `(layers …)` and `(setup (stackup …))` from a `.kicad_pcb` (KiCad 6+). Copper and dielectric layers are kept in order; thickness converted from mm. Dielectric sub-layers (`addsublayer`) are merged into one layer with thickness-weighted εr. Inner copper layers of KiCad type `power`, or whose name looks like GND/PWR/VCC/plane, are marked as reference planes; outer layers never are. The loss tangent of each dielectric is imported when the file has it. Solder mask is reported but not imported; switch it on in the Trace page. Fixture and tests: `src/import/fixtures/`, `kicad.test.ts`.

## Stackup geometry
File: `src/state/store.svelte.ts` → `geometryFor(layerId)`

From the signal layer, walk up and down to the nearest copper marked **Reference plane**. Signal copper crossed on the way adds height and is flagged (dual-stripline coupling not modelled). One plane found → microstrip; two → stripline. `traceStack(layerId)` keeps every dielectric layer separately for the field solver; `geometryFor(layerId)` lumps them into one height and a thickness-weighted εr for the closed-form formulas and the drawing.
