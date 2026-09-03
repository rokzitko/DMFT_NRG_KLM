# Running the DMFT loop

- Use NRG Ljubljana 2026.09 at commit `b3a800e0` or later and Bubble 1.14 or
  later.
- Reproduce the checked-in Bethe tables by running `mkDOS` followed by
  `mkPHI`.
- Modify model and solver settings in `param.loop` if needed.
- `occupancy_mode=accurate` recomputes the lattice spectrum at trial chemical
  potentials with the current self-energy frozen. Use `fast` only when the
  rigid-spectrum approximation is intentionally desired.
- Use `Delta.next.dat` as the restart hybridization; ReDelta and ImDelta are
  reconstructed from it and `param.eps`.
- Start the calculation with `START` from this directory.

## Wilson-chain length

Set either `Nmax` in the range 1 through 998 or a positive `Tmin` in the
`[param]` block. Explicit `Nmax` and `Tmin` are mutually exclusive. The
standalone initializer derives `Nmax` with the same rule as `nrginit`: starting
from zero, it increments while `SCALE[Nmax+1] >= Tmin`, where

```text
SCALE(n) = bandrescale * A * Lambda^(1-z-(n-1)/2)
A = (1+1/Lambda)/2              for discretization Y
A = (1-1/Lambda)/log(Lambda)    for discretizations C and Z
```

Resolution occurs after each `z` value is rendered, so different twists can
have different chain lengths. With the checked-in parameters, `Tmin=1e-7`
gives `Nmax=54,54,53,53` for `z=1/4,1/2,3/4,1`.

A positive `Tmin_ratio` together with a positive explicit `T` first sets
`Tmin=T*Tmin_ratio` and replaces an explicit `Nmax`. An explicit `Tmin` takes
precedence over that provisional value, matching `nrginit` behavior.
