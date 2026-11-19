# The 3D reference runs

Impedora's via model is compared with a 3D finite-element solver, [Palace](https://github.com/awslabs/palace)
(AWS, Apache-2.0). Everything needed to reproduce the stored results is in this folder.

| File | What it is |
|---|---|
| `cases.json` | The via geometries, the single source for both the 3D model and Impedora. |
| `via_case.py` | Builds the Palace mesh (Gmsh) and configuration for one case. |
| `trace_thru.py` | Builds the trace-only "thru" runs that calibrate the trace ports of the cases with a trace. |
| `run.py` | Runs cases inside the container and writes `results.json`. |
| `results.json` | The stored S-parameters. |
| `Dockerfile` | Ubuntu 24.04, Palace (pinned version) and Gmsh. Builds natively on arm64 and amd64. |

## Rerun a case

With Docker installed, from the repository's top folder:

```
docker build -t impedora-fullwave fullwave
docker run --rm -v "$PWD:/work" -w /work impedora-fullwave python3 fullwave/run.py p2g4 --force
```

The image build takes 30 to 40 minutes the first time (it compiles Palace). `run.py` without case names
solves every case that has no stored result; `--level fine` uses the fine mesh (slow), `--level edge` the fine mesh refined along every trace and pad edge (for the cases with a trace and their thrus), `--np 4` sets
the number of MPI processes, `--force` solves a stored run again. Each finished run is written to
`results.json` straight away; working files go to `fullwave/runs/`.

A case with a trace has two companions, `<id>-thru` and `<id>-thru-half` (or those of the case named
in its `thru` field): the same trace between the same ports without the via, at the case's trace length
and at half of it. They calibrate the 3D ports and the 3D trace out of the comparison
(`docs/VALIDATION.md`). Naming a case runs its thrus too.

Compare a new run with Impedora's output in `impedora/` (same frequencies), or post it in an issue.
