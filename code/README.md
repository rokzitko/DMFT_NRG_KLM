# Running the DMFT loop

- Use NRG Ljubljana 2026.09 at commit `b3a800e0` or later and Bubble 1.14 or
  later.
- Reproduce the checked-in Bethe tables by running `mkDOS` followed by
  `mkPHI`.
- Modify model and solver settings in `param.loop` if needed.
- Use `Delta.next.dat` as the restart hybridization; ReDelta and ImDelta are
  reconstructed from it and `param.eps`.
- Start the calculation with `START` from this directory.
