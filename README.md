Example code for dynamical mean-field theory (DMFT) calculation using the [NRG
Ljubljana](https://github.com/rokzitko/nrgljubljana) code as
the impurity solver

Model: Kondo lattice model (KLM), S=1/2; model description in template/.

Features:
- use of templates (i.e., Mathematica is not required at run time)
- improved NRG discretization scheme (Zitko, Pruschke, 2008)
- improved estimator for the self-energy (Kugel, 2022)
- support for arbitrary density of states (tabulated in file DOS.dat)
- robust band occupancy control by shifting chemical potential mu
- simple mixing (implemented at the level of hybridisation function)
- adaptive grid for better capturing sharp spectral features
- transport calculation using external [bubble](https://github.com/rokzitko/bubble) code

Requirements:
- NRG Ljubljana with associated tools (hilb, kk, adapt, nrgchain, broaden, resample, matrix, diag, unitary)
- associated scripts (getparam, scaley, getiter, newiter, subtracty...), in github repo rokzitko/nrgljubljana under scripts/.
- perl
- Python (for occupancy_control)
- m4 macro processor
- bubble (optional)

Two modes of operation:
- local: a script named "mynrgrun" must exist; a minimal version just calls "nrg", but typically you will want to set up
  the environment (number of threads, working directory, piping of output to files); an example is provided
- slurm: a script named "subslurm" must exist to create a job script and submit it to the cluster for execution

Rok Zitko, 2026

# Conventions

This document defines the normalization and sign conventions used by this
repository.  The formulas describe the current unpolarized, single-channel
implementation.  In particular, transport outputs are dimensionless
code-normalized quantities.  They are not conductivities in SI units.

## Units and energy variables

- The half-bandwidth is denoted by `D` and is set to `D = 1`.
- The full bandwidth is therefore `2D = 2`.
- Literature that denotes the Bethe half-bandwidth by `W`, including
  Arsenault and Tremblay, has `W = D` in this document; `W` is not the full
  bandwidth here.
- `k_B = hbar = 1`.
- Temperature, fermionic frequency, optical frequency, chemical potential,
  self-energy, hybridization, interaction strengths, and kinetic energy are
  all expressed in units of `D`.
- Before setting `D = 1`, densities of states and spectral functions have
  units `1/D`.  Their tabulated values are numerical values in the `D = 1`
  system.
- `omega` is a fermionic frequency measured relative to the chemical
  potential, so the Fermi function is

  ```math
  f(\omega)=\frac{1}{1+\exp(\omega/T)}.
  ```

- `Omega` is the external bosonic frequency used for optical conductivity.
  It is distinct from the internal integration frequency `omega`.
- `epsilon` denotes a bare band energy relative to the center of `DOS.dat`.
- `epsilon_d`, stored in `param.eps`, is the common on-site shift of the
  conduction band.  The lattice and transport routines therefore pass
  `mu - epsilon_d` to `hilb` and `bubble`.

The NRG parameter `bandrescale` is an internal solver rescaling and is not the
physical bandwidth.  With `data_has_rescaled_energies=false`, output energies
are in the physical `D = 1` convention above.

Under the usual infinite-coordination Bethe-lattice notation, in which the
scaled hopping is `t_*` and `D = 2 t_*`, this repository has

```math
t_*=\frac12, \qquad t_*^2=\frac14.
```

This hopping convention is inferred from the semicircular DOS; the repository
does not define a finite-coordination bond hopping.

## Model and spin convention

In energy representation, the grand-canonical lattice Hamiltonian represented
by the impurity problem is

```math
K = \sum_{k\sigma}
  (\epsilon_d+\epsilon_k-\mu)c^\dagger_{k\sigma}c_{k\sigma}
  + U\sum_i n_{i\uparrow}n_{i\downarrow}
  + J_K\sum_i \mathbf S_i\mathbin{\cdot}\mathbf s_i
  + B\sum_i(S_i^z+s_i^z).
```

Here

```math
\mathbf S\mathbin{\cdot}\mathbf s
=S^zs^z+\frac12(S^+s^-+S^-s^+).
```

Important details are:

- There is one conduction channel with two spin components.
- The default localized spin is `S = 1/2`.
- `J_K > 0` is antiferromagnetic: the local singlet has exchange energy
  `-3 J_K/4` and the triplet has `J_K/4`.
- The optional `U` acts on the conduction orbital.  The default is `U = 0`.
- The field term has the sign `+B(S^z+s^z)`, rather than the frequently used
  `-B(S^z+s^z)` convention.
- `Himp` in the expectation-value output is grand canonical because it
  contains `(epsilon_d - mu)n`.
- `Hpot` contains only `U n_up n_down + J_K S.s`; it excludes the level,
  chemical-potential, and field terms.

The current DMFT and transport pipeline uses one scalar, spin-averaged
self-energy and `QS` symmetry.  Its supported convention is therefore the
unpolarized `B = 0` case.  The parameter `spin=1/2` refers to the localized
Kondo spin and is unrelated to factors of two from conduction-electron spin.

Unless explicitly stated otherwise:

- `DOS.dat`, one-particle Green functions, and spectral functions are per
  conduction-electron spin.
- The filling `n` is summed over the two conduction-electron spins and ranges
  from zero to two.
- `n_d` is the total two-spin occupancy.
- `n_d_ud` is the double occupancy.
- `ekin.dat` is spin-summed.
- The transport bubbles contain one spin copy and no explicit spin sum.

## Retarded Green functions

The retarded convention is

```math
G^R_{AB}(t)=-i\theta(t)\langle\{A(t),B(0)\}\rangle,
\qquad
G^R_{AB}(\omega)=\int_{-\infty}^{\infty}dt\,
e^{i\omega t}G^R_{AB}(t).
```

The band-energy-resolved lattice Green function is

```math
G^R_\epsilon(\omega)=
\frac{1}{
  \omega+\mu-\epsilon_d-\epsilon-\Sigma^R(\omega)
}.
```

The local Green function and spectral function are

```math
G^R_{\mathrm{loc}}(\omega)
=\int d\epsilon\,\rho_0(\epsilon)G^R_\epsilon(\omega),
\qquad
A(\omega)=-\frac{1}{\pi}\mathrm{Im}\,G^R(\omega).
```

Thus `Im Sigma^R <= 0`, `Im Delta^R <= 0`, and a normalized one-particle
spectral function obeys

```math
\int_{-\infty}^{\infty}d\omega\,A(\omega)=1
```

per spin.

The reconstructed impurity propagator is

```math
G^R_{\mathrm{imp}}(\omega)=
\frac{1}{
  \omega+\mu-\epsilon_d-\Delta^R(\omega)-\Sigma^R(\omega)
}.
```

The improved self-energy estimator uses the auxiliary retarded correlators
`F` and `I`:

```math
\Sigma^R(\omega)=\Sigma_H+I^R(\omega)
-\frac{[F^R(\omega)]^2}{G^R(\omega)}.
```

`Sigma_H` and the self-energy are averaged over up and down spins before they
enter the scalar DMFT loop.

## Bethe density of states

For a general half-bandwidth `D`, the normalized semicircular DOS is

```math
\rho_0(\epsilon)=
\frac{2}{\pi D^2}\sqrt{D^2-\epsilon^2}\,
\Theta(D-|\epsilon|).
```

With `D = 1`, as generated by `code/mkDOS`,

```math
\rho_0(\epsilon)=
\frac{2}{\pi}\sqrt{1-\epsilon^2}\,
\Theta(1-|\epsilon|).
```

Its normalization is per spin:

```math
\int d\epsilon\,\rho_0(\epsilon)=1,
\qquad
\int d\epsilon\,\epsilon^2\rho_0(\epsilon)=\frac14.
```

For this DOS, the Bethe DMFT self-consistency relation may be written

```math
\Delta^R(\omega)=t_*^2G^R_{\mathrm{loc}}(\omega)
=\frac14G^R_{\mathrm{loc}}(\omega).
```

The implementation normally evaluates the more general Hilbert transform
using the tabulated `DOS.dat`, then constructs

```math
\mathcal{G}_0^{-1}=(G^R_{\mathrm{loc}})^{-1}+\Sigma^R,
\qquad
\Delta^R=\omega+\mu-\epsilon_d-\mathcal{G}_0^{-1}.
```

The last equation is the intended convention.  The current
`code/scripts/dmftDOS` implementation omits `-epsilon_d` in this final step.
For nonzero `param.eps`, the freshly written `ReDelta.dat` is consequently
shifted by `+epsilon_d`.  `ImDelta.dat`, and hence the bath spectrum used by
NRG, is unchanged; after mixing, `ReDelta.dat` is regenerated from
`ImDelta.dat` by Kramers-Kronig transformation.  The discrepancy is absent for
the default `param.eps = 0`.

## Numerical file conventions

Several historical filenames use `im` and `re` for quantities that have
already been multiplied by `-1/pi`.  The distinction is important.

| File | Second column |
|---|---|
| `DOS.dat` | Bare `rho_0(epsilon)`, normalized per spin |
| `PHI.dat` | Code-normalized transport function `Phi(epsilon)` |
| `res/c-imG.dat`, `imaw.dat` | `-Im G^R/pi = A`, not `Im G^R` |
| `res/c-reG.dat`, `reaw.dat` | `-Re G^R/pi`, not `Re G^R` |
| `res/c-imF.dat`, `res/c-imI.dat` | `-Im F^R/pi`, `-Im I^R/pi` |
| `res/c-reF.dat`, `res/c-reI.dat` | `-Re F^R/pi`, `-Re I^R/pi` |
| `imsigma.dat`, `resigma.dat` | Actual `Im Sigma^R`, `Re Sigma^R` |
| `ImDelta.dat`, `ReDelta.dat` | Actual `Im Delta^R`, `Re Delta^R` |
| `Delta.dat` | `Gamma=-Im Delta^R >= 0`, used as NRG input |
| `self.dat`, `res/c-self.dat` | Reconstructed impurity spectral function, not the self-energy |
| `dos.dat` | Interacting local DOS written by `bubble`; distinct from uppercase `DOS.dat` |
| `cond.opt-PHI.dat` | `Omega`, `sigma_code(Omega)` |
| `ekin.dat` | Spin-summed kinetic-energy scalar |

`mesh.dat` has two columns for compatibility, but its second column is ignored.
It is zero on a freshly generated mesh but can contain copied hybridization
values after a restart.  Most other numerical tables are whitespace-separated
and have no header.

`DOS.dat` and `PHI.dat` are intended to be strictly increasing, real,
two-column tables.  The currently generated tables contain a
Mathematica-format complex roundoff residue in the row just beyond
`epsilon = 1`.  This row has zero intended physical weight.  The C++ stream
reader in `bubble` 1.5 stops at this non-C++ numeric token, silently ignores the
remaining zero-padding rows, and uses this row as the effective upper table
endpoint.  No nonzero support is lost, but regenerated tables should clamp such
endpoint roundoff to a real zero.

## Bethe transport function

Before setting `D = 1`, the dimensionless kernel corresponding to the
code-normalized longitudinal transport function is

```math
\Phi(\epsilon)
=\frac{D^2-\epsilon^2}{D}\rho_0(\epsilon).
```

For `D = 1`, `code/mkPHI` writes

```math
\Phi(\epsilon)=
\frac{2}{\pi}(1-\epsilon^2)^{3/2}\,
\Theta(1-|\epsilon|).
```

The factor `2/pi` is already part of `PHI.dat`.  The identities needed for the
sum rule are

```math
\Phi(\epsilon)=(1-\epsilon^2)\rho_0(\epsilon),
\qquad
\Phi'(\epsilon)=-3\epsilon\rho_0(\epsilon),
\qquad
\int d\epsilon\,\Phi(\epsilon)=\frac34.
```

This is not the directional normalization used by Arsenault and Tremblay.  In
their convention for a Bethe lattice with `D = 1`,

```math
\Phi_{xx}^{\mathrm{AT}}(\epsilon)
=\frac{1}{3d}(1-\epsilon^2)\rho_0(\epsilon),
```

and hence

```math
\Phi(\epsilon)=3d\,\Phi_{xx}^{\mathrm{AT}}(\epsilon).
```

Here `d` is the formal dimension/current-direction normalization in the
Arsenault-Tremblay Bethe prescription.  It is not a physical dimension or a
parameter used by this repository.

The code has absorbed the directional `1/d` and the additional `1/3` into its
conductivity normalization.  This is natural for a Bethe lattice, which has no
ordinary crystal momentum or uniquely normalized Cartesian current, but it
means that `Phi` must not be mistaken for a physical velocity-squared DOS.

The built-in `bubble` kernel `m=5` is only
`(1-epsilon^2)^(3/2)` and does not contain `2/pi`.  Therefore the two wrappers
use equivalent prefactors:

```math
\begin{aligned}
\texttt{bbl}:&\quad (2/\pi)\pi^2=2\pi, \\[0pt]
\texttt{bbl-PHI}:&\quad \pi^2,
\end{aligned}
```

because the latter obtains `2/pi` from `PHI.dat`.

## Bubble and conductivity normalization

For `D = 1`, define

```math
A_\epsilon(\omega)=-\frac{1}{\pi}
\mathrm{Im}\,G^R_\epsilon(\omega).
```

For `Omega > 0`, the raw finite-frequency result returned by `bubble` with a
tabulated kernel and `n=2`, `o=0` is

```math
B(\Omega)=
\int d\epsilon\,\Phi(\epsilon)
\int d\omega\,
\frac{f(\omega)-f(\omega+\Omega)}{\Omega}
A_\epsilon(\omega)A_\epsilon(\omega+\Omega).
```

At `Omega = 0`, `bubble` uses the corresponding `-df/domega` limit directly.

There is no explicit charge, spin sum, lattice spacing, volume, direction, or
velocity prefactor in `B`.  `code/scripts/cond.opt-PHI` and
`code/scripts/bbl-PHI` define

```math
\boxed{\sigma_{\mathrm{code}}(\Omega)=\pi^2B(\Omega)}.
```

The `pi^2` converts the two normalized spectral functions
`A=-Im G/pi` to the historical convention in which the integrand contains two
factors of `-Im G`.

For comparison, the regular absorptive part in the usual Kubo convention is

```math
\mathrm{Re}\,\sigma_{xx}^{\mathrm{Kubo,reg}}(\Omega)
=\pi e^2\sum_\sigma
\int d\epsilon\,\Phi_{xx}(\epsilon)
\int d\omega\,
\frac{f(\omega)-f(\omega+\Omega)}{\Omega}
A_{\epsilon\sigma}(\omega)A_{\epsilon\sigma}(\omega+\Omega).
```

This spectral-product expression is for `Omega > 0`; it is not the full
complex conductivity and does not display a possible zero-frequency Drude
distribution.

For two degenerate spins and the Arsenault-Tremblay `Phi_xx`, this gives the
purely algebraic relation

```math
\sigma_{\mathrm{code}}
=\frac{3\pi d}{2e^2}\,
\mathrm{Re}\,\sigma_{xx}^{\mathrm{Kubo,reg}},
```

before restoring lattice-spacing and volume factors.  This comparison explains
the normalization, but does not define an SI conversion for the Bethe lattice.

Files named `condMIR.dat` and `rhoMIR.dat` are therefore best read as
dimensionless code-normalized conductivity and its reciprocal.  This
repository does not define `sigma_MIR`, `e`, a lattice spacing, a unit-cell
volume, or an SI restoration factor.

`bubble` evaluates the one-particle bubble with a local scalar self-energy.  It
does not add vertex corrections, a separate diamagnetic term, or a separate
Drude delta function.  The finite-frequency result is the regular absorptive
bubble contribution.  In the infinite-dimensional one-band setting the
longitudinal vertex corrections normally vanish, but this remains part of the
model assumption.

## Kinetic energy

Define the band-energy occupation per spin by

```math
n_\epsilon=\int d\omega\,f(\omega)A_\epsilon(\omega).
```

`code/scripts/ekin` uses `bubble -f`, the kernel
`epsilon rho_0(epsilon)`, and one spectral function.  `bubble` therefore first
computes the one-spin quantity

```math
K_\uparrow=\int d\epsilon\,\epsilon\rho_0(\epsilon)n_\epsilon.
```

The script then explicitly multiplies it by two:

```math
\boxed{E_{\mathrm{kin}}
=2\int d\epsilon\,\epsilon\rho_0(\epsilon)n_\epsilon}.
```

`E_kin` is the conduction-band kinetic energy per DOS-normalized lattice site,
summed over the two spins.  It excludes the on-site `epsilon_d n` term, the
grand-canonical `-mu n` term, the Hubbard interaction, and the Kondo exchange.

## Optical f-sum rule

For one spin, a local causal self-energy, and complete frequency integrals, the
full-frequency bubble identity is

```math
\int_{-\infty}^{\infty}d\Omega\,B(\Omega)
=-\int_{\epsilon_-}^{\epsilon_+}d\epsilon\,
\Phi(\epsilon)\frac{\partial n_\epsilon}{\partial\epsilon}.
```

Integration by parts gives

```math
\int_{-\infty}^{\infty}d\Omega\,B(\Omega)
=\int d\epsilon\,\Phi'(\epsilon)n_\epsilon
-[\Phi(\epsilon)n_\epsilon]_{\epsilon_-}^{\epsilon_+}.
```

For the Bethe kernel, `Phi` vanishes at both band edges, so the boundary term
is zero.  Equivalently, `Phi` may be extended by zero and differentiated in the
distributional sense.

The real optical conductivity is even in `Omega`, so the positive-frequency
form is

```math
\int_0^\infty d\Omega\,B(\Omega)
=\frac12\int d\epsilon\,\Phi'(\epsilon)n_\epsilon.
```

Using the project definition `sigma_code=pi^2 B` and the Bethe identity
`Phi'=-3 epsilon rho_0` gives

```math
\begin{aligned}
\int_0^\infty d\Omega\,\sigma_{\mathrm{code}}(\Omega)
&=\frac{\pi^2}{2}
  \int d\epsilon\,\Phi'(\epsilon)n_\epsilon \\[0pt]
&=-\frac{3\pi^2}{2}
  \int d\epsilon\,\epsilon\rho_0(\epsilon)n_\epsilon \\[0pt]
&=-\frac{3\pi^2}{4}E_{\mathrm{kin}}.
\end{aligned}
```

Thus the check in `code/scripts/sumrule` is

```math
\boxed{
\frac{
  \int_0^\infty d\Omega\,\sigma_{\mathrm{code}}(\Omega)
}{-E_{\mathrm{kin}}}
=\frac{3\pi^2}{4}
}.
```

The factors have distinct origins:

- `3` comes from `Phi'=-3 epsilon rho_0`.
- `pi^2` comes from the project's two-`Im G` conductivity convention.
- One factor `1/2` comes from integrating only positive optical frequencies.
- The other factor `1/2` comes from comparing a one-spin optical bubble with
  the spin-summed `E_kin`.

Equivalently, with the conventional directional transport function and an
explicit two-spin sum, the familiar form is

```math
\int_0^\infty d\Omega\,
\mathrm{Re}\,\sigma_{xx}^{\mathrm{Kubo}}(\Omega)
=-\frac{\pi e^2}{2d}E_{\mathrm{kin}}.
```

Here a zero-frequency Drude distribution, if present, is included with its
positive-frequency half-weight.

Multiplying this expression by the code normalization
`3 pi d/(2 e^2)` gives `-3 pi^2 E_kin/4`.

### Numerical interpretation

`code/scripts/sumrule` is a numerical diagnostic, not an exact pass/fail test.
It uses the external `integrate` script, which applies the trapezoidal rule to
the positive-frequency table without extrapolating to zero or infinity.

For the checked-in reference result,

```math
E_{\mathrm{kin}}=-0.335592587787968,
```

and trapezoidal integration gives

```math
\int d\Omega\,\sigma_{\mathrm{code}}(\Omega)
=2.48796949266544,
\qquad
r=1.00154780145565.
```

The approximately `0.155%` excess is primarily the trapezoidal error on a
geometric grid whose spacing grows by ten percent.  Applying
`scipy.integrate.simpson` to the same irregular grid gives
`r = 1.00001502846936`.

Smaller residual errors can also arise because `DOS.dat` and `PHI.dat` are
interpolated independently, the kinetic and optical calls use different finite
internal-frequency ranges, and numerical clipping changes `Im Sigma` without
reconstructing `Re Sigma`.  The exact identity assumes mutually consistent
kernels, a causal spectral representation, and complete integrals.

The tabulated mesh starts at `Omega=1e-7`, ends at
`Omega=4.54618133178882`, and contains neither zero nor negative frequencies.
The missing ranges are negligible for the checked-in result but need not be
negligible for other parameters.  In a clean limit, any true Drude delta at
zero must also be counted with the appropriate positive-frequency half-weight;
the tabulated regular bubble does not represent a separate delta function.

## DC and thermoelectric outputs

For moment index `o=0,1,2`, let `a_o` denote the raw `bubble` moment after the
same `pi^2` normalization used for conductivity.  The wrappers define

```math
A_0=T a_0, \qquad A_1=a_1, \qquad A_2=\frac{a_2}{T}.
```

They then write

```math
\begin{aligned}
\texttt{condMIR.dat}:&\quad \sigma_{\mathrm{code}}=a_0, \\[0pt]
\texttt{rhoMIR.dat}:&\quad \rho_{\mathrm{code}}=1/a_0, \\[0pt]
\texttt{thermopowerS.dat}:&\quad S=-A_1/A_0, \\[0pt]
\texttt{kappa.dat}:&\quad \kappa=A_2-A_1^2/A_0, \\[0pt]
\texttt{LL.dat}:&\quad L_{\mathrm{raw}}=A_2/A_0, \\[0pt]
\texttt{ZT.dat}:&\quad ZT=S^2a_0T/\kappa.
\end{aligned}
```

`LL.dat` is not the conventional open-circuit Lorenz ratio.  With the code's
own definitions,

```math
\frac{\kappa}{\sigma_{\mathrm{code}}T}
=\frac{A_2}{A_0}-\left(\frac{A_1}{A_0}\right)^2
=L_{\mathrm{raw}}-S^2.
```

Charge, `k_B`, and physical conductivity/thermal-conductivity units are not
restored in these files.

## Numerical floors and integration domains

- `code/scripts/sigmatrick` enforces `Im Sigma^R <= -1e-6`.  This is a
  numerical causality floor, not a physical scattering rate.
- The `clip` parameter in `param.loop` is described there as an `ImSigma`
  floor, but it is actually used as the intended floor for
  `Gamma=-Im Delta` during bath construction and mixing.  The post-mixing
  `clipy` call enforces `Gamma >= 1e-5`.  The preliminary test in
  `code/scripts/dmftDOS` compares `Im Delta` with `+1e-5`, rather than
  `-1e-5`, so the same bound is not guaranteed for the raw or initial bath.
- `bubble` has its own default floor `-Im Sigma >= 1e-8`.  It normally has no
  effect here because the stored self-energy has already been clipped at
  `1e-6`.
- The optical wrappers use a `20T` internal-frequency cutoff and request
  absolute and relative quadrature tolerances `1e-8` and `1e-7`.
- In optical mode, the internal integration interval is approximately
  `[-20T-Omega, 20T]`, and the shifted propagator also requires self-energy
  data through `20T+Omega`.
- The kinetic-energy call uses `bubble -f` defaults.  Its upper occupied-energy
  cutoff is `15T`; its lower endpoint is the first self-energy frequency.
- `bubble` and the auxiliary NRG Ljubljana command-line tools are external and
  are not version-pinned by this repository.  Exact numerical quadrature and
  interpolation behavior can therefore depend on the installed versions.

For a non-Bethe or otherwise modified `DOS.dat`, the coefficient
`3 pi^2/4` is not automatic.  A corresponding `PHI.dat` must use the same
energy scale and satisfy the appropriate relation between `Phi'` and the
stress or kinetic-energy kernel.  It must also vanish at the integration
boundaries, or the boundary term in the general sum rule must be retained.

## References and implementation

- [`bubble` documentation](https://github.com/rokzitko/bubble)
- L.-F. Arsenault and A.-M. S. Tremblay, *Phys. Rev. B* **88**, 205109
  (2013), [doi:10.1103/PhysRevB.88.205109](https://doi.org/10.1103/PhysRevB.88.205109)
- [`code/mkDOS`](code/mkDOS) and [`code/mkPHI`](code/mkPHI)
- [`code/scripts/ekin`](code/scripts/ekin)
- [`code/scripts/cond.opt-PHI`](code/scripts/cond.opt-PHI)
- [`code/scripts/bbl-PHI`](code/scripts/bbl-PHI)
- [`code/scripts/sumrule`](code/scripts/sumrule)
