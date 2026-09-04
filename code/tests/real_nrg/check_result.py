#!/usr/bin/env python3
"""Check structural, causal, and broad physical invariants of the nightly run."""

import argparse
import math
from pathlib import Path
import re
import sys


NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")


class CheckError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise CheckError(message)


def finite_number(text, context):
    require(NUMBER.fullmatch(text) is not None, f"invalid number in {context}: {text!r}")
    value = float(text)
    require(math.isfinite(value), f"non-finite number in {context}: {text!r}")
    return value


def data_lines(path):
    require(path.is_file(), f"missing file: {path.name}")
    lines = [
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    require(lines, f"no data rows in {path.name}")
    return lines


def scalar(path):
    fields = " ".join(data_lines(path)).split()
    require(len(fields) == 1, f"expected one value in {path.name}")
    return finite_number(fields[0], path.name)


def table(path, columns=2):
    rows = []
    for line_number, line in enumerate(data_lines(path), 1):
        fields = line.split()
        require(
            len(fields) == columns,
            f"expected {columns} columns in {path.name} data row {line_number}",
        )
        rows.append(
            [finite_number(field, f"{path.name} data row {line_number}") for field in fields]
        )
    return rows


def same_mesh(reference, candidate, context):
    require(len(reference) == len(candidate), f"mesh length mismatch in {context}")
    for index, (first, second) in enumerate(zip(reference, candidate), 1):
        scale = max(abs(first[0]), abs(second[0]), 1.0)
        require(
            abs(first[0] - second[0]) <= 1e-11 * scale,
            f"frequency mismatch in {context} at row {index}",
        )


def check_hybridization(root, suffix):
    gamma = table(root / f"Delta.{suffix}.dat")
    real = table(root / f"ReDelta.{suffix}.dat")
    imaginary = table(root / f"ImDelta.{suffix}.dat")
    mesh = table(root / f"mesh.{suffix}.dat")
    same_mesh(gamma, real, f"ReDelta.{suffix}.dat")
    same_mesh(gamma, imaginary, f"ImDelta.{suffix}.dat")
    same_mesh(gamma, mesh, f"mesh.{suffix}.dat")
    require(
        all(first[0] < second[0] for first, second in zip(gamma, gamma[1:])),
        f"Delta.{suffix}.dat mesh is not strictly increasing",
    )
    require(gamma[0][1] == 0.0, f"lower Delta.{suffix}.dat endpoint is not zero")
    require(gamma[-1][1] == 0.0, f"upper Delta.{suffix}.dat endpoint is not zero")
    for index, ((_, gamma_value), (_, imaginary_value)) in enumerate(
        zip(gamma, imaginary), 1
    ):
        require(gamma_value >= 0.0, f"negative Gamma at {suffix} row {index}")
        if index not in (1, len(gamma)):
            require(
                gamma_value >= 1e-6 * (1.0 - 1e-11),
                f"Gamma below clipDelta at {suffix} row {index}",
            )
        scale = max(abs(gamma_value), abs(imaginary_value), 1.0)
        require(
            abs(gamma_value + imaginary_value) <= 1e-11 * scale,
            f"ImDelta != -Gamma at {suffix} row {index}",
        )
    return gamma


def key_values(path):
    values = {}
    for line in data_lines(path):
        fields = line.split("=", 1)
        require(len(fields) == 2 and fields[0] and fields[1], f"malformed {path.name}")
        require(fields[0] not in values, f"duplicate {fields[0]} in {path.name}")
        values[fields[0]] = fields[1]
    return values


def custom_values(path):
    lines = path.read_text(encoding="ascii").splitlines()
    headers = [line.lstrip()[1:].split() for line in lines if line.lstrip().startswith("#")]
    names = next((fields for fields in headers if "T" in fields and "n_d" in fields), None)
    require(names is not None, f"missing named header in {path.name}")
    rows = table(path, len(names))
    require(len(rows) == 1, f"expected one averaged row in {path.name}")
    return dict(zip(names, rows[0]))


def check(root):
    require(data_lines(root / "ITER") == ["2"], "nightly run did not publish ITER=2")
    require((root / "STOP").exists(), "STOP marker is missing")
    for name in (
        "CONVERGED",
        "ERROR",
        "POSTPROCESSING_FAILED",
        "NEXT_PREPARATION_FAILED",
    ):
        require(not (root / name).exists(), f"unexpected failure marker: {name}")
    require(not (root / "res").exists(), "staging directory was not removed")
    numeric_directories = [
        path.name for path in root.iterdir() if path.is_dir() and path.name.isdigit()
    ]
    require(not numeric_directories, f"NRG directories remain: {numeric_directories}")

    calls = data_lines(root / "nrg.calls")
    require(calls == ["1", "2", "1", "2"], f"unexpected NRG call order: {calls}")

    occupancy = table(root / "occupancy.log", 4)
    require(len(occupancy) == 2, "occupancy.log does not contain two cycles")
    require(abs(occupancy[-1][3] - 1.0) <= 0.03, "final lattice filling is not near one")
    used_mu = scalar(root / "param.mu.used")
    next_mu = scalar(root / "param.mu.next")
    require(abs(used_mu) <= 0.02, "used mu is not near zero")
    require(abs(next_mu) <= 0.02, "next mu is not near zero")
    require(abs(occupancy[-1][0] - used_mu) <= 1e-10, "occupancy log used mu is stale")
    require(abs(occupancy[-1][1] - next_mu) <= 1e-10, "occupancy log next mu is stale")
    require(
        abs(occupancy[-1][1] - occupancy[-1][0] - occupancy[-1][2]) <= 1e-10,
        "occupancy log step is inconsistent",
    )

    for name in ("DIFFS_C", "DIFFS_LatLoc"):
        rows = table(root / name)
        require(len(rows) == 1, f"{name} must contain one difference row")
        require(rows[0][0] == 2.0, f"{name} does not describe iteration 2")
        require(rows[0][1] >= 0.0, f"{name} contains a negative difference")
    for name in ("logNEW.1", "logUNION.1", "logNEW.2", "logUNION.2"):
        require((root / name).is_file(), f"missing mesh diagnostic: {name}")

    metrics = key_values(root / "OCCUPANCY_METRICS")
    require(metrics.get("version") == "1", "unexpected occupancy metrics version")
    require(metrics.get("iteration") == "2", "occupancy metrics are not from iteration 2")
    require(metrics.get("mode") == "accurate", "accurate occupancy mode was not exercised")
    require(metrics.get("status") == "within_tolerance", "unexpected occupancy update status")
    for key in ("mu_old", "mu_new", "n_old", "n_new", "error_new", "goal", "weight_new"):
        require(key in metrics, f"missing {key} in occupancy metrics")
    metric_mu_old = finite_number(metrics["mu_old"], "occupancy metrics mu_old")
    metric_mu_new = finite_number(metrics["mu_new"], "occupancy metrics mu_new")
    metric_n_old = finite_number(metrics["n_old"], "occupancy metrics n_old")
    metric_n_new = finite_number(metrics["n_new"], "occupancy metrics n_new")
    metric_error_new = finite_number(metrics["error_new"], "occupancy metrics error_new")
    metric_goal = finite_number(metrics["goal"], "occupancy metrics goal")
    metric_weight_new = finite_number(metrics["weight_new"], "occupancy metrics weight_new")
    require(abs(metric_mu_old - used_mu) <= 1e-10, "metrics used mu is stale")
    require(abs(metric_mu_new - next_mu) <= 1e-10, "metrics next mu is stale")
    require(abs(metric_n_old - occupancy[-1][3]) <= 1e-10, "metrics filling is stale")
    require(abs(metric_n_new - 1.0) <= 0.03, "next-cycle filling is not near one")
    require(metric_goal == 1.0, "occupancy goal is not one")
    require(
        abs(metric_error_new - (metric_n_new - metric_goal)) <= 1e-10,
        "occupancy metrics residual is inconsistent",
    )
    require(abs(metric_weight_new - 1.0) <= 0.01, "next-cycle spectral weight is invalid")

    aliases = {
        "ReDelta.dat": "ReDelta.next.dat",
        "ImDelta.dat": "ImDelta.next.dat",
        "Delta.dat": "Delta.next.dat",
        "param.mu": "param.mu.next",
        "mesh.dat": "mesh.next.dat",
    }
    for alias, target in aliases.items():
        path = root / alias
        require(path.is_symlink(), f"{alias} is not a symlink")
        require(path.readlink() == Path(target), f"{alias} does not point to {target}")

    used_gamma = check_hybridization(root, "used")
    check_hybridization(root, "next")

    self_energy = table(root / "imsigma.dat")
    same_mesh(used_gamma, self_energy, "imsigma.dat")
    require(
        all(first[0] < second[0] for first, second in zip(self_energy, self_energy[1:])),
        "imsigma.dat mesh is not strictly increasing",
    )
    require(
        all(value <= -1e-12 * (1.0 - 1e-11) for _, value in self_energy),
        "ImSigma violates the configured causal floor",
    )

    custom = custom_values(root / "customfdm.avg")
    for name in ("Himp", "Hpot", "SdSk", "SigmaHd", "n_d"):
        require(name in custom and math.isfinite(custom[name]), f"missing finite {name}")
    require(abs(custom["n_d"] - 1.0) <= 0.03, "impurity filling is not near one")
    require(abs(custom["SigmaHd"]) <= 1e-8, "U=0 Hartree shift is not zero")
    require(-0.75 <= custom["SdSk"] < -0.05, "Kondo correlation is outside smoke bounds")
    require(
        abs(custom["Hpot"] - 0.4 * custom["SdSk"]) <= 5e-4,
        "Hpot is inconsistent with JK*SdSk",
    )
    require(
        abs(custom["Himp"] - (custom["Hpot"] - used_mu * custom["n_d"])) <= 5e-3,
        "Himp is inconsistent with the impurity Hamiltonian",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", nargs="?", type=Path, default=Path("."))
    arguments = parser.parse_args()
    try:
        check(arguments.workdir.resolve())
    except (CheckError, OSError, UnicodeError) as exc:
        print(f"check_result: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print("Real-NRG nightly result is structurally, causally, and physically valid.")


if __name__ == "__main__":
    main()
