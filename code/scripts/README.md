How to initialize calculation using results from a previous run

- Prefer copying Delta.next.dat, param.eps, and param.mu.next. Delta.next.dat
  stores the authoritative Gamma=-ImDelta table.
- ReDelta.next.dat and ImDelta.next.dat are regenerated from Delta.next.dat.
- A legacy ImDelta.dat regular file is accepted only when Delta.dat is absent;
  it is converted to Gamma during initialization.
- Start DMFT script

Occupancy control

- `occupancy_control --accurate` is the default loop evaluator;
  `--fast` rigidly shifts the represented spectrum without extrapolation.
- `occupancy.log` keeps its historical four-column format.
  `OCCUPANCY_METRICS` records the filling residual, proposed or previously
  applied mu step, spectral weights, solve status, and evaluation count.
- Accurate standard updates stage their final `H_0` real/spectral pair for
  one-time reuse by `dmftDOS-stable`; its marker binds the pair to the mu,
  clipping value, frozen self-energy, and DOS.
