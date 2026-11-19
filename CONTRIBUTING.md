# Contributing reference data

Impedora's models are checked against exact solutions, independent solvers and a 3D field solver.
None of it has been compared with a measured board yet. Measured or simulated data from you can
change that.

## What helps most

- **Measured S-parameters** (Touchstone) of a via, a trace or a test coupon, with the board's
  stackup as built: copper and dielectric thicknesses, the laminate (and its εr and loss tangent at
  the frequency you have), trace widths, drill, pad and antipad sizes, stitching vias, backdrill depth.
  Say how the fixture was de-embedded (2x-thru, TRL, none).
- **3D simulations** (HFSS, CST, openEMS, Palace, …) with the same geometry details, the port setup
  and the mesh or convergence settings.
- A cross-section photo or the fabricator's stackup report makes measured data far more useful.

## How to send it

Open an issue with the files attached, or a pull request that adds a folder `contrib/<short-name>/`
with the Touchstone files and a `case.md` describing the geometry. Say how you would like to be
credited (name, company, or anonymous). Only send data you are allowed to share; it will be
published here under CC BY 4.0.

## What happens next

The case is rebuilt in Impedora's stackup format, compared, and added to the Validation page with the
difference, whatever it is. If Impedora is off, that is the most valuable result.
