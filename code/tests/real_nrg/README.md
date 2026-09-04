# Reduced real-NRG test

`run` copies `code/` to a fresh directory, seeds a particle-hole-symmetric
Bethe bath, and runs two DMFT cycles with the real NRG Ljubljana solver. `Nz=2`
exercises z-averaging; the deliberately insufficient two-row convergence window
and inclusive `maxiter=2` make `STOP` the deterministic terminal result.

The fixture is excluded from the fast TAP suite. With NRG Ljubljana 2026.09 and
its source scripts on `PATH`, run it explicitly as:

```sh
code/tests/real_nrg/run
```

Set `DMFT_REAL_NRG_WORKDIR` to an empty directory to retain the calculation at a
known path, or `DMFT_REAL_NRG_TIMEOUT` to change the default `75m` timeout.
