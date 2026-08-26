# Quasiparticle renormalization factor `Z`

`Z` estimates the quasiparticle renormalization factor from a two-column
real-part of self-energy file.

It computes

```text
Z = 1 / (1 - d ReSigma(omega) / d omega | omega=0).
```

## Usage

Run the tool from a calculation directory containing `resigma.dat`.

The default output is only the numeric value of `Z`, followed by a newline. This
makes it suitable for shell scripts. A different input file may be supplied as a
positional argument:

```sh
scripts/Z path/to/resigma.dat
```

Print details about the selected range and fit:

```sh
scripts/Z --report
```

Generate `Z.pdf` in the current directory:

```sh
scripts/Z --plot
```

Options may be combined:

```sh
scripts/Z --report --plot path/to/resigma.dat
```

If no stable low-frequency range can be identified, the tool normally prints an
error to standard error, produces no value on standard output, and exits with a
non-zero status. To obtain the best available estimate anyway, use:

```sh
scripts/Z --force
```

A forced estimate is marked in `--report` output, and the reason it was forced is
always written to standard error. Invalid or insufficient input cannot be
overridden with `--force`.

## Range selection

The positive and negative frequency branches are paired to form local symmetric
slope estimates,

```text
[ReSigma(+omega) - ReSigma(-omega)] / (2 omega).
```

For a smooth low-frequency expansion these estimates are linear in `omega^2` up
to cubic order in `ReSigma`. The tool searches logarithmically spaced candidate
windows for the first stable low-frequency plateau. Candidate fits are checked
against a higher-order fit and against fits to the two halves of the window.
Huber weighting and median-absolute-deviation clipping suppress isolated NRG
artifacts. The final value is obtained from a robust cubic fit to the original
`ReSigma` data in the selected range.

The input must contain finite, strictly increasing frequencies on both sides of
zero. Slightly different positive and negative meshes are handled by
interpolation over their common range.

## Report fields

`--report` prints:

- `Z`
- `dReSigma/domega`, the fitted derivative at zero
- the selected absolute-frequency range
- retained and rejected fit-point counts
- the spread of `Z` across accepted neighboring ranges
- fit status and any warnings

Values outside the conventional Fermi-liquid interval `0 < Z <= 1` are reported
but produce a warning on standard error.

## Plot

`--plot` writes a four-panel PDF showing the raw data, selected points, rejected
points, cubic fit, and tangent at zero. The panels use:

1. An overview with asymptotic high-frequency tails removed.
2. The full frequency range on a logarithmic `|omega|` axis, with separate lines
   for positive and negative frequencies.
3. The range containing the 100 input points nearest zero by `|omega|`.
4. Ten times the frequency range of the third panel.

For the overview, the high-frequency limit is estimated by fitting the outer 10%
of the data to a quadratic in `1/omega`. The plot retains frequencies where the
departure from that limit is at least 10% of the robust maximum departure. The
overview is never narrower than the fourth panel, and all ranges are clipped to
the available input data.

NumPy is required. Matplotlib is required only when `--plot` is used.
