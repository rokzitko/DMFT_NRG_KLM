How to initialize calculation using results from a previous run

- Prefer copying Delta.next.dat, param.eps, and param.mu.next. Delta.next.dat
  stores the authoritative Gamma=-ImDelta table.
- ReDelta.next.dat and ImDelta.next.dat are regenerated from Delta.next.dat.
- A legacy ImDelta.dat regular file is accepted only when Delta.dat is absent;
  it is converted to Gamma during initialization.
- Start DMFT script
