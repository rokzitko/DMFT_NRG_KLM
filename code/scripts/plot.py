#!/usr/bin/env python3
"""Generate diagnostic plots for a single DMFT calculation directory."""

from __future__ import annotations

import argparse
import configparser
from dataclasses import dataclass, field
import logging
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, PowerNorm, TwoSlopeNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.integrate import cumulative_simpson


logging.getLogger("fontTools.ttLib.tables._h_e_a_d").setLevel(logging.ERROR)


PNG_DPI = 200
SPECTRAL_ZOOM_T = 10.0
SELF_ENERGY_ZOOM_T = 5.0
MESH_BIN_WIDTH_DECADES = 0.1
MAP_EPSILON_POINTS = 700
MAP_OMEGA_POINTS = 1100
MAP_COLOR_PERCENTILE = 99.5

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7A5195"
SKY = "#56B4E9"
VERMILLION = "#CC3311"
GRAY = "#666666"


class PlotError(RuntimeError):
    """Raised when a run cannot be plotted safely."""


@dataclass(frozen=True)
class RunParameters:
    temperature: float
    lambda_: float
    broaden_min: float
    broaden_max: float
    broaden_ratio: float
    bandrescale: float
    nz: int
    goal: float
    conveps: float
    clip_delta: float
    clip_sigma: float
    clip_g: float
    alpha: float
    opt_ymin: float
    jk: float
    u: float
    field: float
    spin: str


@dataclass
class RunContext:
    root: Path
    output: Path
    parameters: RunParameters
    tables: dict[str, np.ndarray] = field(default_factory=dict)
    scalars: dict[str, float] = field(default_factory=dict)

    def table(
        self,
        name: str,
        *,
        columns: int | None = None,
        minimum_rows: int = 1,
        increasing_x: bool = False,
    ) -> np.ndarray:
        if name not in self.tables:
            path = self.root / name
            if not path.is_file():
                raise PlotError(f"required root-level input is missing: {path}")
            try:
                data = np.loadtxt(path, dtype=float, ndmin=2)
            except (OSError, ValueError) as exc:
                raise PlotError(f"failed to read {path}: {exc}") from exc
            if data.ndim != 2 or not np.all(np.isfinite(data)):
                raise PlotError(f"expected a finite numeric table in {path}")
            self.tables[name] = data

        data = self.tables[name]
        if columns is not None and data.shape[1] != columns:
            raise PlotError(
                f"expected {columns} columns in {self.root / name}, "
                f"found {data.shape[1]}"
            )
        if data.shape[0] < minimum_rows:
            raise PlotError(
                f"expected at least {minimum_rows} rows in {self.root / name}, "
                f"found {data.shape[0]}"
            )
        if increasing_x and np.any(np.diff(data[:, 0]) <= 0.0):
            raise PlotError(f"first column is not strictly increasing in {self.root / name}")
        return data

    def scalar(self, name: str) -> float:
        if name not in self.scalars:
            path = self.root / name
            if not path.is_file():
                raise PlotError(f"required root-level input is missing: {path}")
            try:
                values = np.asarray(np.loadtxt(path, dtype=float)).reshape(-1)
            except (OSError, ValueError) as exc:
                raise PlotError(f"failed to read {path}: {exc}") from exc
            if values.size != 1 or not np.isfinite(values[0]):
                raise PlotError(f"expected one finite value in {path}")
            self.scalars[name] = float(values[0])
        return self.scalars[name]


def parse_parameters(path: Path) -> RunParameters:
    if not path.is_file():
        raise PlotError(f"parameter file is missing: {path}")

    parser = configparser.ConfigParser(interpolation=None)
    try:
        with path.open(encoding="ascii") as handle:
            parser.read_file(handle)
    except (OSError, configparser.Error, UnicodeError) as exc:
        raise PlotError(f"failed to parse {path}: {exc}") from exc

    def value(section: str, option: str) -> str:
        if not parser.has_option(section, option):
            raise PlotError(f"missing [{section}] {option} in {path}")
        return parser.get(section, option)

    def number(section: str, option: str) -> float:
        raw = value(section, option)
        try:
            result = float(raw)
        except ValueError as exc:
            raise PlotError(f"invalid [{section}] {option}={raw!r} in {path}") from exc
        if not np.isfinite(result):
            raise PlotError(f"non-finite [{section}] {option} in {path}")
        return result

    def clip_value(option: str) -> float:
        if parser.has_option("dmft", option):
            return number("dmft", option)
        return number("dmft", "clip")

    result = RunParameters(
        temperature=number("param", "T"),
        lambda_=number("param", "Lambda"),
        broaden_min=number("param", "broaden_min"),
        broaden_max=number("param", "broaden_max"),
        broaden_ratio=number("param", "broaden_ratio"),
        bandrescale=number("param", "bandrescale"),
        nz=int(number("dmft", "Nz")),
        goal=number("dmft", "goal"),
        conveps=number("dmft", "conveps"),
        clip_delta=clip_value("clipDelta"),
        clip_sigma=clip_value("clipSigma"),
        clip_g=clip_value("clipG"),
        alpha=number("dmft", "alpha"),
        opt_ymin=number("dmft", "opt_ymin"),
        jk=number("extra", "JK"),
        u=number("extra", "U1"),
        field=number("extra", "B"),
        spin=value("extra", "spin"),
    )

    positive = {
        "T": result.temperature,
        "Lambda": result.lambda_,
        "broaden_min": result.broaden_min,
        "broaden_max": result.broaden_max,
        "broaden_ratio": result.broaden_ratio,
        "bandrescale": result.bandrescale,
        "Nz": result.nz,
        "conveps": result.conveps,
        "clipDelta": result.clip_delta,
        "clipSigma": result.clip_sigma,
        "clipG": result.clip_g,
        "opt_ymin": result.opt_ymin,
    }
    for name, item in positive.items():
        if item <= 0.0:
            raise PlotError(f"expected positive parameter {name}, found {item}")
    if result.lambda_ <= 1.0:
        raise PlotError(f"expected Lambda > 1, found {result.lambda_}")
    if result.broaden_ratio <= 1.0:
        raise PlotError(f"expected broaden_ratio > 1, found {result.broaden_ratio}")
    return result


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.0,
            "axes.titlesize": 12.0,
            "axes.labelsize": 10.5,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
            "legend.fontsize": 8.8,
            "lines.linewidth": 1.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white",
            "axes.formatter.use_mathtext": True,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(axis: plt.Axes, *, grid: bool = True) -> None:
    if grid:
        axis.grid(True, color="#000000", alpha=0.12, linewidth=0.6)
    axis.tick_params(direction="out", length=3.5, width=0.7)


def format_math_number(value: float, significant_digits: int = 6) -> str:
    """Format a number for use inside a Matplotlib mathtext expression."""
    text = f"{value:.{significant_digits}g}"
    if "e" not in text.lower():
        return text
    mantissa, exponent_text = text.lower().split("e", maxsplit=1)
    exponent = int(exponent_text)
    if mantissa == "1":
        return rf"10^{{{exponent}}}"
    if mantissa == "-1":
        return rf"-10^{{{exponent}}}"
    return rf"{mantissa}\times 10^{{{exponent}}}"


def save_figure(context: RunContext, figure: plt.Figure, stem: str, title: str) -> None:
    metadata = {"Title": title, "Creator": "plot.py"}
    figure.savefig(context.output / f"{stem}.png", dpi=PNG_DPI, metadata=metadata)
    figure.savefig(context.output / f"{stem}.pdf", metadata=metadata)
    plt.close(figure)


def require_matching_x(first: np.ndarray, second: np.ndarray, description: str) -> None:
    if first.shape[0] != second.shape[0] or not np.allclose(
        first[:, 0], second[:, 0], rtol=1e-11, atol=0.0
    ):
        raise PlotError(f"incompatible frequency meshes for {description}")


def bare_band_support(context: RunContext) -> tuple[float, float]:
    dos = context.table("DOS.dat", columns=2, minimum_rows=3, increasing_x=True)
    epsilon, density = dos.T
    threshold = max(float(np.max(np.abs(density))) * 1e-13, np.finfo(float).tiny)
    indices = np.flatnonzero(density > threshold)
    if indices.size == 0:
        raise PlotError(f"DOS has no positive support in {context.root / 'DOS.dat'}")
    lower_index = max(int(indices[0]) - 1, 0)
    upper_index = min(int(indices[-1]) + 1, epsilon.size - 1)
    lower = float(epsilon[lower_index])
    upper = float(epsilon[upper_index])
    if lower >= upper:
        raise PlotError("could not infer a nonzero bare-band interval from DOS.dat")
    return lower, upper


def padded_limits(lower: float, upper: float, fraction: float) -> tuple[float, float]:
    span = upper - lower
    return lower - fraction * span, upper + fraction * span


def clip_limits(
    limits: tuple[float, float], available: tuple[float, float]
) -> tuple[float, float]:
    lower = max(limits[0], available[0])
    upper = min(limits[1], available[1])
    if lower >= upper:
        return available
    return lower, upper


def overview_limits(context: RunContext, x: np.ndarray) -> tuple[float, float]:
    band_lower, band_upper = bare_band_support(context)
    half_width = 0.5 * (band_upper - band_lower)
    requested = band_lower - half_width, band_upper + half_width
    return clip_limits(requested, (float(x[0]), float(x[-1])))


def low_frequency_limits(
    context: RunContext, x: np.ndarray, temperature_multiple: float
) -> tuple[float, float]:
    radius = temperature_multiple * context.parameters.temperature
    requested = -radius, radius
    limits = clip_limits(requested, (float(x[0]), float(x[-1])))
    if np.count_nonzero((x >= limits[0]) & (x <= limits[1])) < 10:
        nearest = np.argsort(np.abs(x))[: min(20, x.size)]
        radius = float(np.max(np.abs(x[nearest])))
        limits = -radius, radius
    return limits


def significant_support_limits(
    x: np.ndarray, y: np.ndarray, fraction: float
) -> tuple[float, float]:
    threshold = fraction * float(np.max(y))
    indices = np.flatnonzero(y > threshold)
    if indices.size == 0:
        raise PlotError(f"no values exceed {fraction:g} of the tabulated maximum")
    lower_index = max(int(indices[0]) - 1, 0)
    upper_index = min(int(indices[-1]) + 1, x.size - 1)
    return float(x[lower_index]), float(x[upper_index])


def cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    increments = np.diff(x) * (y[1:] + y[:-1]) * 0.5
    return np.concatenate(([0.0], np.cumsum(increments)))


def fermi_function(omega: np.ndarray, temperature: float) -> np.ndarray:
    result = np.empty_like(omega, dtype=float)
    positive = omega >= 0.0
    exp_negative = np.exp(-omega[positive] / temperature)
    result[positive] = exp_negative / (1.0 + exp_negative)
    exp_positive = np.exp(omega[~positive] / temperature)
    result[~positive] = 1.0 / (1.0 + exp_positive)
    return result


def load_custom_average(context: RunContext) -> dict[str, float]:
    path = context.root / "customfdm.avg"
    data = context.table("customfdm.avg", minimum_rows=1)
    if data.shape[0] != 1:
        raise PlotError(f"expected one averaged data row in {path}")

    try:
        comments = [
            line.lstrip()[1:].strip().split()
            for line in path.read_text(encoding="ascii").splitlines()
            if line.lstrip().startswith("#")
        ]
    except (OSError, UnicodeError) as exc:
        raise PlotError(f"failed to read headers from {path}: {exc}") from exc

    names: list[str] | None = None
    for tokens in reversed(comments):
        if len(tokens) == data.shape[1] and "T" in tokens:
            names = tokens
            break
    if names is None:
        raise PlotError(f"could not identify column names in {path}")
    result = dict(zip(names, data[0], strict=True))
    required = {"T", "Himp", "Hpot", "SdSk", "n_d", "n_d_ud"}
    missing = sorted(required - result.keys())
    if missing:
        raise PlotError(f"missing columns in {path}: {', '.join(missing)}")
    return {name: float(value) for name, value in result.items()}


def add_colored_line(
    axis: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    color_values: np.ndarray,
    *,
    cmap: str,
    norm: Normalize,
    linewidth: float = 2.0,
) -> LineCollection:
    points = np.column_stack((x, y)).reshape(-1, 1, 2)
    segments = np.concatenate((points[:-1], points[1:]), axis=1)
    collection = LineCollection(segments, cmap=cmap, norm=norm, linewidth=linewidth)
    collection.set_array(0.5 * (color_values[:-1] + color_values[1:]))
    axis.add_collection(collection)
    axis.update_datalim(np.column_stack((x, y)))
    axis.autoscale_view()
    return collection


def plot_convergence(context: RunContext) -> None:
    consecutive = context.table("DIFFS_C", columns=2, minimum_rows=1)
    closure = context.table("DIFFS_LatLoc", columns=2, minimum_rows=1)
    if np.any(consecutive[:, 1] <= 0.0) or np.any(closure[:, 1] <= 0.0):
        raise PlotError("convergence errors must be positive for a logarithmic axis")

    figure, axis = plt.subplots(figsize=(7.2, 4.5), layout="constrained")
    axis.plot(
        consecutive[:, 0],
        consecutive[:, 1],
        marker="o",
        markersize=4.0,
        color=BLUE,
        label="consecutive lattice GF spectra",
    )
    axis.plot(
        closure[:, 0],
        closure[:, 1],
        marker="s",
        markersize=3.7,
        color=ORANGE,
        label="lattice vs. local spectrum",
    )
    tolerance = context.parameters.conveps
    axis.axhline(tolerance, color=GRAY, linestyle="--", linewidth=1.1)
    axis.annotate(
        rf"convergence tolerance $={format_math_number(tolerance, 2)}$",
        xy=(0.985, tolerance),
        xycoords=("axes fraction", "data"),
        xytext=(-2, 5),
        textcoords="offset points",
        ha="right",
        va="bottom",
        color=GRAY,
        fontsize=8.5,
    )
    axis.set_yscale("log")
    axis.set_xlabel("DMFT step")
    axis.set_ylabel("integrated error")
    axis.set_title("DMFT convergence")
    axis.xaxis.set_major_locator(MaxNLocator(integer=True))
    axis.legend(frameon=False, loc="upper right")
    style_axis(axis)
    save_figure(context, figure, "01_convergence", "DMFT convergence")


def plot_hybridization(context: RunContext) -> None:
    delta = context.table("Delta.dat", columns=2, minimum_rows=3, increasing_x=True)
    omega, gamma = delta.T
    if np.any(gamma <= 0.0):
        raise PlotError("Delta.dat must contain positive Gamma=-Im Delta")
    support_limits = significant_support_limits(omega, gamma, 1e-3)

    figure, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), layout="constrained")
    axes[0].plot(omega, gamma, color=BLUE)
    axes[0].set_xlim(*support_limits)
    axes[0].set_xlabel(r"$\omega$")
    axes[0].set_ylabel(r"$\Gamma(\omega)=-\mathrm{Im}\,\Delta(\omega)$")
    axes[0].set_title(r"$\Gamma>10^{-3}\Gamma_{\max}$ window")
    style_axis(axes[0])

    axes[1].plot(omega, gamma, color=BLUE)
    axes[1].axhline(
        context.parameters.clip_delta,
        color=ORANGE,
        linestyle="--",
        linewidth=1.1,
        label=rf"$\mathrm{{clipDelta}}={format_math_number(context.parameters.clip_delta, 2)}$",
    )
    axes[1].set_xscale(
        "symlog",
        linthresh=max(context.parameters.temperature, context.parameters.broaden_min),
    )
    axes[1].set_yscale("log")
    axes[1].set_xlim(float(omega[0]), float(omega[-1]))
    axes[1].set_xlabel(r"$\omega$")
    axes[1].set_ylabel(r"$\Gamma(\omega)$")
    axes[1].set_title("full tabulated range")
    axes[1].legend(frameon=False, loc="center")
    style_axis(axes[1], grid=False)

    figure.suptitle("hybridization function")
    save_figure(context, figure, "02_hybridization", "hybridization function")


def plot_discretization(context: RunContext) -> None:
    positive = context.table("FSOL.dat", columns=2, minimum_rows=2, increasing_x=True)
    negative = context.table("FSOLNEG.dat", columns=2, minimum_rows=2, increasing_x=True)
    require_matching_x(positive, negative, "FSOL positive and negative branches")

    asymptote = (1.0 - 1.0 / context.parameters.lambda_) / math.log(
        context.parameters.lambda_
    )
    figure, axis = plt.subplots(figsize=(7.2, 4.5), layout="constrained")
    axis.plot(positive[:, 0], positive[:, 1], color=BLUE, label=r"positive $\omega$")
    axis.plot(negative[:, 0], negative[:, 1], color=ORANGE, label=r"negative $\omega$")
    axis.axhline(
        asymptote,
        color=GRAY,
        linestyle="--",
        linewidth=1.0,
        label=rf"$(1-\Lambda^{{-1}})/\ln\Lambda={format_math_number(asymptote, 4)}$",
    )
    axis.set_xlabel(r"discretization coordinate $x$")
    axis.set_ylabel(r"$f_\pm(x)$")
    axis.set_title(r"discretization functions $f(x)$")
    axis.legend(frameon=False, loc="lower right")
    style_axis(axis)
    save_figure(
        context,
        figure,
        "03_discretization_functions",
        "discretization functions f(x)",
    )


def plot_optical_conductivity(context: RunContext) -> None:
    optical = context.table(
        "cond.opt-PHI.dat", columns=2, minimum_rows=3, increasing_x=True
    )
    omega, conductivity = optical.T
    sigma_dc = context.scalar("condMIR.dat")
    resistivity = context.scalar("rhoMIR.dat")
    thermopower = context.scalar("thermopowerS.dat")
    kappa = context.scalar("kappa.dat")
    zt = context.scalar("ZT.dat")
    if sigma_dc <= 0.0 or np.any(omega <= 0.0) or np.any(conductivity <= 0.0):
        raise PlotError("optical and DC conductivities must be positive")

    normalized = conductivity / sigma_dc
    cumulative = cumulative_trapezoid(conductivity, omega)
    if cumulative[-1] > 0.0:
        index = int(np.searchsorted(cumulative, 0.995 * cumulative[-1]))
        index = min(index, omega.size - 1)
        linear_upper = min(float(omega[-1]), 1.1 * float(omega[index]))
    else:
        linear_upper = float(omega[-1])

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(13.8, 4.3),
        layout="constrained",
        gridspec_kw={"width_ratios": (1.25, 1.25, 0.9)},
    )
    axes[0].plot(np.concatenate(([0.0], omega)), np.concatenate(([1.0], normalized)), color=BLUE)
    axes[0].axhline(1.0, color=GRAY, linestyle="--", linewidth=1.0)
    axes[0].set_xlim(0.0, linear_upper)
    axes[0].set_ylim(bottom=0.0)
    axes[0].set_xlabel(r"$\Omega$")
    axes[0].set_ylabel(r"$\sigma(\Omega)/\sigma_{\mathrm{dc}}$")
    style_axis(axes[0])

    axes[1].loglog(omega, normalized, color=BLUE)
    axes[1].axhline(1.0, color=GRAY, linestyle="--", linewidth=1.0)
    axes[1].set_xlabel(r"$\Omega$")
    axes[1].set_ylabel(r"$\sigma(\Omega)/\sigma_{\mathrm{dc}}$")
    style_axis(axes[1], grid=False)

    axes[2].axis("off")
    values = [
        rf"$T = {format_math_number(context.parameters.temperature, 7)}$",
        rf"$\sigma_{{\mathrm{{dc}}}} = {format_math_number(sigma_dc, 7)}$",
        rf"$\rho_{{\mathrm{{dc}}}} = {format_math_number(resistivity, 7)}$",
        rf"$S = {format_math_number(thermopower, 7)}$",
        rf"$\kappa_\mathrm{{e}} = {format_math_number(kappa, 7)}$",
        rf"$ZT = {format_math_number(zt, 7)}$",
    ]
    axes[2].text(
        0.05,
        0.93,
        "transport",
        transform=axes[2].transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
    )
    axes[2].text(
        0.05,
        0.79,
        "\n".join(values),
        transform=axes[2].transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        linespacing=1.55,
    )
    figure.suptitle("optical conductivity")
    save_figure(
        context, figure, "04_optical_conductivity", "optical conductivity"
    )


def plot_local_spectral_function(context: RunContext) -> None:
    spectrum = context.table("imaw.dat", columns=2, minimum_rows=3, increasing_x=True)
    omega, spectral = spectrum.T
    custom = load_custom_average(context)
    occupancy = context.table("occupancy.log", columns=4, minimum_rows=1)
    lattice_n = float(occupancy[-1, 3])
    chemical_potential = context.scalar("param.mu")
    kinetic_energy = context.scalar("ekin.dat")

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(13.8, 4.4),
        layout="constrained",
        gridspec_kw={"width_ratios": (1.35, 1.15, 1.0)},
    )
    axes[0].plot(omega, spectral, color=BLUE)
    axes[0].set_xlim(*overview_limits(context, omega))
    axes[0].set_ylim(bottom=0.0)
    axes[0].axvline(0.0, color=GRAY, linestyle=":", linewidth=1.0)
    axes[0].set_xlabel(r"$\omega$")
    axes[0].set_ylabel(r"$A_{\mathrm{lat}}(\omega)$")
    axes[0].set_title("overview")
    style_axis(axes[0])

    axes[1].plot(omega, spectral, color=BLUE)
    axes[1].set_xlim(*low_frequency_limits(context, omega, SPECTRAL_ZOOM_T))
    axes[1].set_ylim(bottom=0.0)
    axes[1].axvline(0.0, color=GRAY, linestyle=":", linewidth=1.0)
    axes[1].set_xlabel(r"$\omega$")
    axes[1].set_ylabel(r"$A_{\mathrm{lat}}(\omega)$")
    axes[1].set_title(rf"Fermi-level view ($|\omega|\leq {SPECTRAL_ZOOM_T:g}T$)")
    style_axis(axes[1])

    axes[2].axis("off")
    values = [
        rf"$T = {format_math_number(context.parameters.temperature, 7)}$",
        rf"$J_K = {format_math_number(context.parameters.jk, 7)}$",
        rf"$\mu = {format_math_number(chemical_potential, 9)}$",
        rf"$\langle n\rangle_\mathrm{{lat}} = {format_math_number(lattice_n, 9)}$",
        rf"$\langle n_d\rangle_\mathrm{{FDM}} = {format_math_number(custom['n_d'], 9)}$",
        rf"$E_\mathrm{{kin}} = {format_math_number(kinetic_energy, 9)}$",
        rf"$\langle H_\mathrm{{imp}}\rangle = {format_math_number(custom['Himp'], 9)}$",
        rf"$\langle H_\mathrm{{pot}}\rangle = {format_math_number(custom['Hpot'], 9)}$",
        rf"$\langle \mathbf{{S}}_d\!\cdot\!\mathbf{{S}}_K\rangle = {format_math_number(custom['SdSk'], 9)}$",
        rf"$\langle n_{{d\uparrow}}n_{{d\downarrow}}\rangle = {format_math_number(custom['n_d_ud'], 9)}$",
    ]
    axes[2].text(
        0.02,
        0.97,
        "\n".join(values),
        transform=axes[2].transAxes,
        ha="left",
        va="top",
        fontsize=9.3,
        linespacing=1.45,
    )

    figure.suptitle("local spectral function")
    save_figure(
        context,
        figure,
        "05_local_spectral_function",
        "local spectral function",
    )


def load_self_energy(context: RunContext) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    real = context.table("resigma.dat", columns=2, minimum_rows=3, increasing_x=True)
    imaginary = context.table(
        "imsigma.dat", columns=2, minimum_rows=3, increasing_x=True
    )
    require_matching_x(real, imaginary, "real and imaginary self-energy")
    return real[:, 0], real[:, 1], imaginary[:, 1]


def self_energy_figure(
    context: RunContext,
    limits: tuple[float, float],
    stem: str,
    subtitle: str,
    *,
    show_component_titles: bool = True,
) -> None:
    omega, real, imaginary = load_self_energy(context)
    figure, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), layout="constrained")
    axes[0].plot(omega, real, color=BLUE)
    axes[0].axhline(0.0, color=GRAY, linestyle=":", linewidth=0.9)
    axes[0].axvline(0.0, color=GRAY, linestyle=":", linewidth=0.9)
    axes[0].set_xlim(*limits)
    axes[0].set_xlabel(r"$\omega$")
    axes[0].set_ylabel(r"$\mathrm{Re}\,\Sigma(\omega)$")
    if show_component_titles:
        axes[0].set_title("real part")
    style_axis(axes[0])

    axes[1].plot(omega, imaginary, color=ORANGE)
    axes[1].axhline(0.0, color=GRAY, linestyle=":", linewidth=0.9)
    axes[1].axvline(0.0, color=GRAY, linestyle=":", linewidth=0.9)
    axes[1].set_xlim(*limits)
    axes[1].set_xlabel(r"$\omega$")
    axes[1].set_ylabel(r"$\mathrm{Im}\,\Sigma(\omega)$")
    if show_component_titles:
        axes[1].set_title("imaginary part")
    style_axis(axes[1])

    title = f"self-energy: {subtitle}"
    figure.suptitle(title)
    save_figure(context, figure, stem, title)


def plot_self_energy(context: RunContext) -> None:
    omega, _, _ = load_self_energy(context)
    self_energy_figure(
        context,
        overview_limits(context, omega),
        "06_self_energy_overview",
        "overview",
        show_component_titles=False,
    )
    self_energy_figure(
        context,
        low_frequency_limits(context, omega, SELF_ENERGY_ZOOM_T),
        "07_self_energy_low_frequency",
        "low frequencies",
    )


def branch_density(log_frequency: np.ndarray, edges: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(log_frequency, bins=edges)
    return counts / np.diff(edges)


def local_linear_density(omega: np.ndarray) -> np.ndarray:
    if omega.size < 3:
        raise PlotError("at least three branch points are needed for local mesh density")
    spacing = np.gradient(omega)
    if np.any(spacing <= 0.0):
        raise PlotError("mesh branch is not strictly increasing")
    return 1.0 / spacing


def plot_mesh_density(context: RunContext) -> None:
    mesh = context.table("mesh.dat", columns=2, minimum_rows=3, increasing_x=True)
    delta = context.table("Delta.dat", columns=2, minimum_rows=3, increasing_x=True)
    omega = mesh[:, 0]
    positive = omega[omega > 0.0]
    signed_negative = omega[omega < 0.0]
    negative = -signed_negative[::-1]
    if positive.size < 3 or signed_negative.size < 3:
        raise PlotError("mesh.dat must contain positive and negative frequency branches")

    all_absolute = np.concatenate((positive, negative))
    lower_decade = math.floor(float(np.log10(np.min(all_absolute))) / MESH_BIN_WIDTH_DECADES) * MESH_BIN_WIDTH_DECADES
    upper_decade = math.ceil(float(np.log10(np.max(all_absolute))) / MESH_BIN_WIDTH_DECADES) * MESH_BIN_WIDTH_DECADES
    edges = np.arange(
        lower_decade,
        upper_decade + 1.5 * MESH_BIN_WIDTH_DECADES,
        MESH_BIN_WIDTH_DECADES,
    )
    centers = 10.0 ** (0.5 * (edges[:-1] + edges[1:]))
    positive_density = branch_density(np.log10(positive), edges)
    negative_density = branch_density(np.log10(negative), edges)
    baseline = math.log(10.0) / math.log(context.parameters.broaden_ratio)

    signed_negative_density = local_linear_density(signed_negative)
    positive_linear_density = local_linear_density(positive)
    base_factor = math.log(context.parameters.broaden_ratio)
    negative_base = 1.0 / (np.abs(signed_negative) * base_factor)
    positive_base = 1.0 / (positive * base_factor)
    linear_limits = significant_support_limits(delta[:, 0], delta[:, 1], 1e-3)

    figure, axes = plt.subplots(1, 2, figsize=(13.8, 4.7), layout="constrained")
    logarithmic_axis, linear_axis = axes
    if positive.shape == negative.shape and np.allclose(positive, negative, rtol=1e-11):
        logarithmic_axis.step(
            centers,
            positive_density,
            where="mid",
            color=BLUE,
            label="positive and negative branches",
        )
        logarithmic_axis.fill_between(
            centers, positive_density, step="mid", color=BLUE, alpha=0.12
        )
    else:
        logarithmic_axis.step(
            centers, positive_density, where="mid", color=BLUE, label=r"$\omega>0$"
        )
        logarithmic_axis.step(
            centers, negative_density, where="mid", color=ORANGE, label=r"$\omega<0$"
        )
    logarithmic_axis.axhline(
        baseline,
        color=GRAY,
        linestyle="--",
        linewidth=1.0,
        label=rf"base mesh: ${format_math_number(baseline, 4)}$ points/decade",
    )
    logarithmic_axis.axvline(
        context.parameters.temperature,
        color=GREEN,
        linestyle=":",
        linewidth=1.1,
        label=rf"$T={format_math_number(context.parameters.temperature, 3)}$",
    )
    logarithmic_axis.set_xscale("log")
    logarithmic_axis.set_xlim(float(np.min(all_absolute)), float(np.max(all_absolute)))
    logarithmic_axis.set_ylim(bottom=0.0)
    logarithmic_axis.set_xlabel(r"$|\omega|$")
    logarithmic_axis.set_ylabel("mesh points per decade")
    logarithmic_axis.set_title(r"logarithmic $|\omega|$ scale")
    logarithmic_axis.legend(frameon=False, loc="upper left")
    style_axis(logarithmic_axis, grid=False)

    linear_axis.plot(
        signed_negative,
        signed_negative_density,
        color=ORANGE,
        label=r"mesh density, $\omega<0$",
    )
    linear_axis.plot(
        positive, positive_linear_density, color=BLUE, label=r"mesh density, $\omega>0$"
    )
    linear_axis.plot(
        signed_negative,
        negative_base,
        color=GRAY,
        linestyle="--",
        linewidth=1.1,
        label="base geometric density",
    )
    linear_axis.plot(positive, positive_base, color=GRAY, linestyle="--", linewidth=1.1)
    linear_axis.axvline(0.0, color="#000000", linestyle=":", linewidth=0.9)
    linear_axis.set_yscale("log")
    linear_axis.set_xlim(*linear_limits)
    linear_axis.set_xlabel(r"$\omega$")
    linear_axis.set_ylabel(r"mesh points per unit $\omega$")
    linear_axis.set_title(r"linear $\omega$ scale")
    linear_axis.legend(frameon=False, loc="upper right")
    style_axis(linear_axis, grid=False)

    figure.suptitle("frequency-mesh density")
    save_figure(context, figure, "08_mesh_density", "frequency-mesh density")


def plot_bare_inputs(context: RunContext) -> None:
    dos = context.table("DOS.dat", columns=2, minimum_rows=3, increasing_x=True)
    phi = context.table("PHI.dat", columns=2, minimum_rows=3, increasing_x=True)
    require_matching_x(dos, phi, "DOS and transport function")
    epsilon = dos[:, 0]
    dos_integral = float(np.trapezoid(dos[:, 1], epsilon))
    phi_integral = float(np.trapezoid(phi[:, 1], epsilon))

    figure, axis = plt.subplots(figsize=(7.2, 4.5), layout="constrained")
    axis.plot(epsilon, dos[:, 1], color=BLUE, label=r"$\rho_0(\epsilon)$")
    axis.plot(epsilon, phi[:, 1], color=ORANGE, label=r"$\Phi(\epsilon)$")
    axis.set_xlim(float(epsilon[0]), float(epsilon[-1]))
    axis.set_ylim(bottom=0.0)
    axis.set_xlabel(r"bare-band energy $\epsilon$")
    axis.set_ylabel("tabulated function")
    axis.set_title("bare density of states and transport function")
    axis.legend(frameon=False, loc="upper left")
    axis.text(
        0.98,
        0.95,
        rf"trapezoid: $\int\rho_0\,d\epsilon={format_math_number(dos_integral, 8)}$"
        + "\n"
        + rf"$\int\Phi\,d\epsilon={format_math_number(phi_integral, 8)}$",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.8,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "#cccccc"},
    )
    style_axis(axis)
    save_figure(
        context,
        figure,
        "09_bare_dos_transport",
        "bare density of states and transport function",
    )


def epsilon_resolved_data(
    epsilon: np.ndarray,
    omega: np.ndarray,
    real_sigma: np.ndarray,
    imaginary_sigma: np.ndarray,
    chemical_potential: float,
) -> np.ndarray:
    gamma = -imaginary_sigma[:, None]
    if np.any(gamma <= 0.0):
        raise PlotError("Im Sigma must be negative when constructing A(epsilon, omega)")
    denominator = (
        omega[:, None]
        + chemical_potential
        - real_sigma[:, None]
        - epsilon[None, :]
    ) ** 2 + gamma**2
    return gamma / (np.pi * denominator)


def plot_epsilon_resolved_spectrum(context: RunContext) -> None:
    omega, real_sigma, imaginary_sigma = load_self_energy(context)
    chemical_potential = context.scalar("param.mu")
    epsilon_lower, epsilon_upper = bare_band_support(context)
    epsilon = np.linspace(epsilon_lower, epsilon_upper, MAP_EPSILON_POINTS)
    broad_limits = overview_limits(context, omega)
    zoom_limits = low_frequency_limits(context, omega, SPECTRAL_ZOOM_T)
    broad_omega = np.linspace(broad_limits[0], broad_limits[1], MAP_OMEGA_POINTS)
    zoom_omega = np.linspace(zoom_limits[0], zoom_limits[1], MAP_OMEGA_POINTS)

    broad_real = np.interp(broad_omega, omega, real_sigma)
    broad_imaginary = np.interp(broad_omega, omega, imaginary_sigma)
    zoom_real = np.interp(zoom_omega, omega, real_sigma)
    zoom_imaginary = np.interp(zoom_omega, omega, imaginary_sigma)
    broad_spectral = epsilon_resolved_data(
        epsilon,
        broad_omega,
        broad_real,
        broad_imaginary,
        chemical_potential,
    )
    zoom_spectral = epsilon_resolved_data(
        epsilon,
        zoom_omega,
        zoom_real,
        zoom_imaginary,
        chemical_potential,
    )
    color_max = float(np.percentile(broad_spectral, MAP_COLOR_PERCENTILE))
    if color_max <= 0.0:
        raise PlotError("epsilon-resolved spectrum has no positive intensity")
    color_norm = PowerNorm(gamma=0.5, vmin=0.0, vmax=color_max, clip=True)

    figure, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), layout="constrained")
    meshes = []
    for axis, plot_omega, plot_real, spectral, subtitle in (
        (axes[0], broad_omega, broad_real, broad_spectral, "overview"),
        (
            axes[1],
            zoom_omega,
            zoom_real,
            zoom_spectral,
            rf"Fermi-level view ($|\omega|\leq {SPECTRAL_ZOOM_T:g}T$)",
        ),
    ):
        mesh = axis.pcolormesh(
            epsilon,
            plot_omega,
            spectral,
            shading="auto",
            cmap="magma",
            norm=color_norm,
            rasterized=True,
        )
        meshes.append(mesh)
        ridge = plot_omega + chemical_potential - plot_real
        inside = (ridge >= epsilon_lower) & (ridge <= epsilon_upper)
        axis.plot(
            ridge[inside],
            plot_omega[inside],
            color="white",
            linewidth=1.0,
            alpha=0.9,
            label=r"$\epsilon=\omega+\mu-\mathrm{Re}\,\Sigma$",
        )
        bare_omega = epsilon - chemical_potential
        visible = (bare_omega >= plot_omega[0]) & (bare_omega <= plot_omega[-1])
        axis.plot(
            epsilon[visible],
            bare_omega[visible],
            color=SKY,
            linestyle="--",
            linewidth=1.0,
            alpha=0.9,
            label="bare dispersion",
        )
        sigma_zero = float(np.interp(0.0, omega, real_sigma))
        epsilon_fermi = chemical_potential - sigma_zero
        if epsilon_lower <= epsilon_fermi <= epsilon_upper:
            axis.axvline(epsilon_fermi, color="white", linestyle=":", linewidth=0.9, alpha=0.8)
        axis.axhline(0.0, color="white", linestyle=":", linewidth=0.9, alpha=0.8)
        axis.set_xlim(epsilon_lower, epsilon_upper)
        axis.set_ylim(float(plot_omega[0]), float(plot_omega[-1]))
        axis.set_xlabel(r"bare-band energy $\epsilon$")
        axis.set_ylabel(r"frequency $\omega$")
        axis.set_title(subtitle)
        axis.tick_params(direction="out", length=3.5, width=0.7)

    axes[0].legend(
        frameon=True,
        facecolor="white",
        framealpha=0.78,
        fontsize=8.0,
        loc="upper left",
    )
    colorbar = figure.colorbar(meshes[0], ax=axes, shrink=0.96, pad=0.02)
    colorbar.set_ticks(np.linspace(0.0, color_max, 5))
    colorbar.set_label(r"$A(\epsilon,\omega)$")
    figure.suptitle(r"$\epsilon$-resolved spectral function")
    save_figure(
        context,
        figure,
        "10_epsilon_resolved_spectrum",
        "epsilon-resolved spectral function",
    )


def plot_effective_medium(context: RunContext) -> None:
    omega, real_sigma, imaginary_sigma = load_self_energy(context)
    spectrum = context.table("imaw.dat", columns=2, minimum_rows=3, increasing_x=True)
    require_matching_x(
        np.column_stack((omega, real_sigma)), spectrum, "self-energy and lattice spectrum"
    )
    chemical_potential = context.scalar("param.mu")
    limits = overview_limits(context, omega)
    selected = (omega >= limits[0]) & (omega <= limits[1])
    w = omega[selected]
    zeta_real = w + chemical_potential - real_sigma[selected]
    zeta_imaginary = -imaginary_sigma[selected]
    spectral = spectrum[selected, 1]
    band_lower, band_upper = bare_band_support(context)

    if w[0] < 0.0 < w[-1]:
        color_norm: Normalize = TwoSlopeNorm(vmin=float(w[0]), vcenter=0.0, vmax=float(w[-1]))
    else:
        color_norm = Normalize(vmin=float(w[0]), vmax=float(w[-1]))

    figure, axes = plt.subplots(1, 2, figsize=(11.2, 4.5), layout="constrained")
    colored = add_colored_line(
        axes[0],
        zeta_real,
        zeta_imaginary,
        w,
        cmap="coolwarm",
        norm=color_norm,
        linewidth=2.2,
    )
    axes[0].plot(
        [band_lower, band_upper],
        [0.0, 0.0],
        color="black",
        linewidth=3.0,
        solid_capstyle="butt",
        label="bare-band support",
    )
    zeta_zero = (
        chemical_potential - float(np.interp(0.0, omega, real_sigma)),
        -float(np.interp(0.0, omega, imaginary_sigma)),
    )
    axes[0].scatter(*zeta_zero, s=34, color="white", edgecolor="black", zorder=4, label=r"$\omega=0$")
    axes[0].set_xlabel(r"$\mathrm{Re}\,\zeta(\omega)$")
    axes[0].set_ylabel(r"$\mathrm{Im}\,\zeta(\omega)=-\mathrm{Im}\,\Sigma(\omega)$")
    axes[0].set_title(r"$\zeta(\omega)=\omega+\mu-\Sigma(\omega)$")
    axes[0].legend(frameon=False, loc="upper right")
    style_axis(axes[0])

    add_colored_line(
        axes[1],
        w,
        spectral,
        w,
        cmap="coolwarm",
        norm=color_norm,
        linewidth=2.2,
    )
    axes[1].axvline(0.0, color=GRAY, linestyle=":", linewidth=1.0)
    axes[1].set_xlim(*limits)
    axes[1].set_ylim(bottom=0.0)
    axes[1].set_xlabel(r"$\omega$")
    axes[1].set_ylabel(r"$A_{\mathrm{lat}}(\omega)$")
    axes[1].set_title("corresponding lattice spectrum")
    style_axis(axes[1])

    colorbar = figure.colorbar(colored, ax=axes, shrink=0.93, pad=0.02)
    colorbar.set_label(r"frequency $\omega$")
    figure.suptitle("complex effective medium")
    save_figure(
        context, figure, "11_effective_medium", "complex effective medium"
    )


def plot_optical_sum_rule(context: RunContext) -> None:
    optical = context.table(
        "cond.opt-PHI.dat", columns=2, minimum_rows=3, increasing_x=True
    )
    omega, conductivity = optical.T
    kinetic_energy = context.scalar("ekin.dat")
    expected_weight = 3.0 * np.pi**2 * (-kinetic_energy) / 4.0
    if expected_weight <= 0.0:
        raise PlotError("optical sum rule requires a negative kinetic energy")
    cumulative = cumulative_trapezoid(conductivity, omega)
    ratio = cumulative / expected_weight
    crossing = np.flatnonzero(ratio >= 0.99)

    figure, axis = plt.subplots(figsize=(7.3, 4.5), layout="constrained")
    axis.plot(omega, ratio, color=BLUE)
    axis.axhline(1.0, color=GRAY, linestyle="--", linewidth=1.1, label="sum-rule value")
    if crossing.size:
        index = int(crossing[0])
        axis.scatter(omega[index], ratio[index], color=ORANGE, s=30, zorder=3)
        axis.annotate(
            rf"99% at $\Omega={format_math_number(omega[index], 4)}$",
            xy=(omega[index], ratio[index]),
            xytext=(8, -18),
            textcoords="offset points",
            fontsize=8.5,
            color=ORANGE,
        )
    axis.text(
        0.98,
        0.08,
        rf"trapezoidal final ratio $={format_math_number(ratio[-1], 8)}$",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.0,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "#cccccc"},
    )
    axis.set_xscale("log")
    axis.set_xlabel(r"optical frequency $\Omega$")
    axis.set_ylabel(r"$\int_0^\Omega \sigma(\Omega')d\Omega' / [(3\pi^2/4)(-E_{\mathrm{kin}})]$")
    axis.set_title("trapezoidal cumulative optical weight")
    axis.legend(frameon=False, loc="center left")
    style_axis(axis, grid=False)
    save_figure(
        context, figure, "12_optical_sum_rule", "cumulative optical sum rule"
    )


def plot_spectral_closure(context: RunContext) -> None:
    lattice = context.table("imaw.dat", columns=2, minimum_rows=3, increasing_x=True)
    local = context.table("self.dat", columns=2, minimum_rows=3, increasing_x=True)
    independent = context.table("ldos.dat", columns=2, minimum_rows=3, increasing_x=True)
    require_matching_x(lattice, local, "lattice and reconstructed local spectra")
    omega = lattice[:, 0]
    limits = overview_limits(context, omega)
    independent_selected = (independent[:, 0] >= omega[0]) & (independent[:, 0] <= omega[-1])
    independent_omega = independent[independent_selected, 0]
    independent_spectral = independent[independent_selected, 1]
    lattice_on_independent = np.interp(independent_omega, omega, lattice[:, 1])

    lattice_weight = cumulative_simpson(lattice[:, 1], x=omega, initial=0.0)
    local_weight = cumulative_simpson(local[:, 1], x=omega, initial=0.0)
    fermi = fermi_function(omega, context.parameters.temperature)
    lattice_occupancy = 2.0 * cumulative_simpson(
        lattice[:, 1] * fermi, x=omega, initial=0.0
    )
    local_occupancy = 2.0 * cumulative_simpson(
        local[:, 1] * fermi, x=omega, initial=0.0
    )

    figure, axes = plt.subplots(2, 2, figsize=(11.4, 8.0), layout="constrained")
    axes[0, 0].plot(omega, lattice[:, 1], color=BLUE, label="lattice, adaptive mesh")
    axes[0, 0].plot(omega, local[:, 1], color=ORANGE, linestyle="--", label="reconstructed local")
    axes[0, 0].plot(
        independent_omega,
        independent_spectral,
        color=GREEN,
        linewidth=1.2,
        alpha=0.9,
        label="lattice, independent linear grid",
    )
    axes[0, 0].set_xlim(*limits)
    axes[0, 0].set_ylim(bottom=0.0)
    axes[0, 0].set_xlabel(r"$\omega$")
    axes[0, 0].set_ylabel(r"$A(\omega)$")
    axes[0, 0].set_title("spectral functions")
    axes[0, 0].legend(frameon=False, loc="center left")
    style_axis(axes[0, 0])

    axes[0, 1].plot(omega, local[:, 1] - lattice[:, 1], color=PURPLE)
    axes[0, 1].axhline(0.0, color=GRAY, linestyle=":", linewidth=0.9)
    axes[0, 1].set_xlim(*limits)
    axes[0, 1].set_xlabel(r"$\omega$")
    axes[0, 1].set_ylabel(r"$A_{\mathrm{local}}-A_{\mathrm{lat}}$")
    axes[0, 1].set_title("DMFT closure residual")
    axes[0, 1].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    style_axis(axes[0, 1])

    axes[1, 0].plot(
        independent_omega,
        independent_spectral - lattice_on_independent,
        color=GREEN,
    )
    axes[1, 0].axhline(0.0, color=GRAY, linestyle=":", linewidth=0.9)
    axes[1, 0].set_xlim(*limits)
    axes[1, 0].set_xlabel(r"$\omega$")
    axes[1, 0].set_ylabel(r"$A_{\mathrm{linear}}-A_{\mathrm{adaptive}}$")
    axes[1, 0].set_title("independent lattice-method residual")
    axes[1, 0].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    style_axis(axes[1, 0])

    axes[1, 1].plot(
        omega,
        local_weight - lattice_weight,
        color=BLUE,
        label=r"$\Delta\int^\omega A(\omega')d\omega'$",
    )
    axes[1, 1].plot(
        omega,
        local_occupancy - lattice_occupancy,
        color=ORANGE,
        label=r"$\Delta[2\int^\omega fA\,d\omega']$",
    )
    axes[1, 1].axhline(0.0, color=GRAY, linestyle=":", linewidth=0.9)
    axes[1, 1].set_xlim(*limits)
    axes[1, 1].set_xlabel(r"$\omega$")
    axes[1, 1].set_ylabel("local minus lattice cumulative value")
    axes[1, 1].set_title("cumulative sum-rule residuals")
    axes[1, 1].legend(frameon=False, loc="upper left")
    axes[1, 1].text(
        0.98,
        0.05,
        (
            rf"$W_\mathrm{{lat}}={format_math_number(lattice_weight[-1], 9)}$, "
            rf"$W_\mathrm{{local}}={format_math_number(local_weight[-1], 9)}$"
            "\n"
            rf"$n_\mathrm{{lat}}={format_math_number(lattice_occupancy[-1], 9)}$, "
            rf"$n_\mathrm{{local}}={format_math_number(local_occupancy[-1], 9)}$"
        ),
        transform=axes[1, 1].transAxes,
        ha="right",
        va="bottom",
        fontsize=8.0,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.85, "edgecolor": "#cccccc"},
    )
    style_axis(axes[1, 1])

    figure.suptitle("spectral consistency and DMFT closure")
    save_figure(
        context,
        figure,
        "13_spectral_closure",
        "spectral consistency and DMFT closure",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate PNG and PDF diagnostics for one DMFT run"
    )
    parser.add_argument(
        "run_directory",
        nargs="?",
        default=Path.cwd(),
        type=Path,
        help="calculation directory containing param.loop (default: current directory)",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = args.run_directory.expanduser().resolve()
    if not root.is_dir():
        raise PlotError(f"calculation directory does not exist: {root}")
    parameters = parse_parameters(root / "param.loop")
    output = root / "plots"
    output.mkdir(parents=True, exist_ok=True)
    context = RunContext(root=root, output=output, parameters=parameters)

    configure_matplotlib()
    plot_convergence(context)
    plot_hybridization(context)
    plot_discretization(context)
    plot_optical_conductivity(context)
    plot_local_spectral_function(context)
    plot_self_energy(context)
    plot_mesh_density(context)
    plot_bare_inputs(context)
    plot_epsilon_resolved_spectrum(context)
    plot_effective_medium(context)
    plot_optical_sum_rule(context)
    plot_spectral_closure(context)

    print(f"Wrote 13 PNG and 13 PDF figures to {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PlotError as exc:
        print(f"plot.py: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
