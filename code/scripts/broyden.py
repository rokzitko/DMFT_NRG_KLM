#!/usr/bin/env python3
"""Modified Broyden mixing following D. D. Johnson, PRB 38, 12807 (1988)."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


W0 = 0.01
PARAM_FILE = "param.loop"


class BroydenError(RuntimeError):
    pass


def command_output(command):
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise BroydenError(f"Failed to run {command[0]}: {exc}") from exc

    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        raise BroydenError(
            f"{' '.join(command)} failed with status {result.returncode}: {detail}"
        )
    return result.stdout.strip()


def get_parameter(name):
    result = subprocess.run(
        ["getparam", name, PARAM_FILE],
        capture_output=True,
        text=True,
        check=False,
    )
    value = result.stdout.strip()
    return value if result.returncode == 0 and value else None


def get_float_parameter(name, default=None):
    value = get_parameter(name)
    if value is None:
        if default is not None:
            return default
        raise BroydenError(f"Required parameter {name} is not defined in {PARAM_FILE}")
    try:
        number = float(value)
    except ValueError as exc:
        raise BroydenError(f"Invalid {name}={value!r} in {PARAM_FILE}") from exc
    if not np.isfinite(number):
        raise BroydenError(f"Invalid non-finite {name}={value!r} in {PARAM_FILE}")
    return number


def get_max_history():
    value = get_parameter("broydenM")
    if value is None:
        return None
    try:
        maximum = int(value)
    except ValueError as exc:
        raise BroydenError(f"Invalid M={value!r} in {PARAM_FILE}") from exc
    if maximum < 1 or str(maximum) != value.strip():
        raise BroydenError(f"M must be a positive integer, got {value!r}")
    return maximum


def load_table(path):
    try:
        table = np.loadtxt(path, dtype=float, ndmin=2)
    except (OSError, ValueError) as exc:
        raise BroydenError(f"Failed to read {path}: {exc}") from exc

    if table.ndim != 2 or table.shape[0] < 2 or table.shape[1] != 2:
        raise BroydenError(f"Expected a two-column table with at least two rows in {path}")
    if not np.all(np.isfinite(table)):
        raise BroydenError(f"Non-finite value found in {path}")
    if not np.all(np.diff(table[:, 0]) > 0.0):
        raise BroydenError(f"Frequency mesh is not strictly increasing in {path}")
    return table[:, 0], table[:, 1]


def same_mesh(first, second):
    if first.shape != second.shape:
        return False
    scale = max(float(np.max(np.abs(first))), float(np.max(np.abs(second))), 1.0)
    return np.allclose(first, second, rtol=1e-12, atol=1e-13 * scale)


def interpolate(path, mesh):
    source_mesh, values = load_table(path)
    scale = max(float(np.max(np.abs(source_mesh))), float(np.max(np.abs(mesh))), 1.0)
    tolerance = 1e-12 * scale
    if mesh[0] < source_mesh[0] - tolerance or mesh[-1] > source_mesh[-1] + tolerance:
        raise BroydenError(f"Mesh in {path} does not cover the current frequency range")
    return np.interp(mesh, source_mesh, values)


def load_scalar(path):
    try:
        fields = path.read_text(encoding="ascii").split()
    except OSError as exc:
        raise BroydenError(f"Failed to read {path}: {exc}") from exc
    if len(fields) != 1:
        raise BroydenError(f"Expected one numeric value in {path}")
    try:
        value = float(fields[0])
    except ValueError as exc:
        raise BroydenError(f"Invalid numeric value in {path}: {fields[0]!r}") from exc
    if not np.isfinite(value):
        raise BroydenError(f"Non-finite value found in {path}")
    return value


def load_mu_occupancy(path):
    try:
        fields = path.read_text(encoding="ascii").split()
    except OSError as exc:
        raise BroydenError(f"Failed to read {path}: {exc}") from exc
    if len(fields) != 2:
        raise BroydenError(f"Expected 'mu occupancy' in {path}")
    try:
        mu, occupancy = (float(field) for field in fields)
    except ValueError as exc:
        raise BroydenError(f"Invalid mu/occupancy data in {path}") from exc
    if not np.isfinite(mu) or not np.isfinite(occupancy):
        raise BroydenError(f"Non-finite mu/occupancy data in {path}")
    return mu, occupancy


def atomic_copy(source, destination):
    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, destination)
    except OSError as exc:
        temporary.unlink(missing_ok=True)
        raise BroydenError(f"Failed to archive {source} as {destination}: {exc}") from exc


def write_table(path, mesh, values):
    temporary = path.with_name(f".{path.name}.tmp")
    try:
        np.savetxt(temporary, np.column_stack((mesh, values)), fmt="%.17g")
        os.replace(temporary, path)
    except OSError as exc:
        temporary.unlink(missing_ok=True)
        raise BroydenError(f"Failed to write {path}: {exc}") from exc


def write_scalar(path, value):
    if path.is_symlink():
        target = os.readlink(path)
        path = path.parent / target
    temporary = path.with_name(f".{path.name}.tmp")
    try:
        temporary.write_text(f"{value:.17g}\n", encoding="ascii")
        os.replace(temporary, path)
    except OSError as exc:
        temporary.unlink(missing_ok=True)
        raise BroydenError(f"Failed to write {path}: {exc}") from exc


def input_re_path(history_dir, iteration):
    return history_dir / f"{iteration}-ReDelta.dat"


def input_im_path(history_dir, iteration):
    return history_dir / f"{iteration}-ImDelta.dat"


def input_mu_path(history_dir, iteration):
    return history_dir / f"{iteration}-param.mu"


def raw_re_path(history_dir, iteration):
    return history_dir / f"{iteration}-ReDelta.raw.dat"


def raw_im_path(history_dir, iteration):
    return history_dir / f"{iteration}-ImDelta.raw.dat"


def occupancy_path(history_dir, iteration):
    return history_dir / f"{iteration}-mu-occup.dat"


def required_history_paths(history_dir, iteration, mix_delta, control_mu):
    paths = []
    if mix_delta:
        paths.extend(
            (
                input_re_path(history_dir, iteration),
                input_im_path(history_dir, iteration),
                raw_re_path(history_dir, iteration),
                raw_im_path(history_dir, iteration),
            )
        )
    if control_mu:
        paths.extend(
            (
                input_mu_path(history_dir, iteration),
                occupancy_path(history_dir, iteration),
            )
        )
    return paths


def history_indices(history_dir, iteration, mix_delta, control_mu, maximum):
    indices = []
    for candidate in range(iteration, 0, -1):
        missing = [
            path
            for path in required_history_paths(
                history_dir, candidate, mix_delta, control_mu
            )
            if not path.is_file() or path.stat().st_size == 0
        ]
        if missing:
            if candidate == iteration:
                names = ", ".join(str(path) for path in missing)
                raise BroydenError(f"Current Broyden history is incomplete: {names}")
            break
        indices.append(candidate)

    indices.reverse()
    if maximum is not None:
        indices = indices[-maximum:]
    return indices


def build_history(history_dir, indices, mesh, mix_delta, control_mu, goal):
    inputs = []
    residuals = []
    current_occupancy = None

    for iteration in indices:
        input_parts = []
        residual_parts = []

        if mix_delta:
            input_re = interpolate(input_re_path(history_dir, iteration), mesh)
            input_im = interpolate(input_im_path(history_dir, iteration), mesh)
            output_re = interpolate(raw_re_path(history_dir, iteration), mesh)
            output_im = interpolate(raw_im_path(history_dir, iteration), mesh)
            input_parts.extend((input_re, input_im))
            residual_parts.extend((output_re - input_re, output_im - input_im))

        if control_mu:
            input_mu = load_scalar(input_mu_path(history_dir, iteration))
            recorded_mu, occupancy = load_mu_occupancy(
                occupancy_path(history_dir, iteration)
            )
            if not np.isclose(input_mu, recorded_mu, rtol=1e-11, atol=1e-13):
                raise BroydenError(
                    f"Chemical potential mismatch in iteration {iteration}: "
                    f"{input_mu} != {recorded_mu}"
                )
            input_parts.append(np.array([input_mu]))
            residual_parts.append(np.array([goal - occupancy]))
            if iteration == indices[-1]:
                current_occupancy = occupancy

        inputs.append(np.concatenate(input_parts))
        residuals.append(np.concatenate(residual_parts))

    return np.stack(inputs), np.stack(residuals), current_occupancy


def broyden_step(inputs, residuals, initial_factors):
    current_input = inputs[-1]
    current_residual = residuals[-1]
    mixed = current_input + initial_factors * current_residual

    delta_residuals = []
    correction_vectors = []
    for previous in range(len(inputs) - 1):
        delta_residual = residuals[previous + 1] - residuals[previous]
        norm = np.linalg.norm(delta_residual)
        scale = max(
            np.linalg.norm(residuals[previous]),
            np.linalg.norm(residuals[previous + 1]),
            1.0,
        )
        if norm <= 10.0 * np.finfo(float).eps * scale:
            print(
                f"Skipping degenerate Broyden pair "
                f"{previous + 1}->{previous + 2}"
            )
            continue

        normalized_residual = delta_residual / norm
        normalized_input = (inputs[previous + 1] - inputs[previous]) / norm
        delta_residuals.append(normalized_residual)
        correction_vectors.append(
            initial_factors * normalized_residual + normalized_input
        )

    if not delta_residuals:
        return mixed, 0, 1.0

    # This is the matrix form of the double sum in the old broyden.m code.
    delta_matrix = np.column_stack(delta_residuals)
    correction_matrix = np.column_stack(correction_vectors)
    system = delta_matrix.T @ delta_matrix
    system += W0**2 * np.eye(system.shape[0])
    projection = delta_matrix.T @ current_residual
    coefficients = np.linalg.solve(system, projection)
    mixed -= correction_matrix @ coefficients
    return mixed, len(delta_residuals), float(np.linalg.cond(system))


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Apply modified Johnson-Broyden mixing to the current DMFT iteration"
    )
    parser.add_argument("--delta", action="store_true", help="mix ReDelta and ImDelta")
    parser.add_argument("--mu", action="store_true", help="control param.mu")
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("."),
        help="directory containing current Delta files (default: .)",
    )
    parser.add_argument(
        "--iteration",
        type=int,
        help="DMFT iteration number (default: getiter output)",
    )
    parser.add_argument(
        "--mu-input",
        type=Path,
        default=Path("param.mu"),
        help="current chemical-potential file (default: param.mu)",
    )
    parser.add_argument(
        "--mu-output",
        type=Path,
        default=Path("param.mu"),
        help="next chemical-potential file (default: param.mu)",
    )
    parser.add_argument(
        "--history-dir",
        type=Path,
        default=Path("dmft"),
        help="Broyden history directory (default: dmft)",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=Path("occupancy.log"),
        help="occupancy update log (default: occupancy.log)",
    )
    arguments = parser.parse_args()
    if not arguments.delta and not arguments.mu:
        parser.error("at least one of --delta or --mu is required")
    return arguments


def main():
    arguments = parse_arguments()
    if not arguments.history_dir.is_dir():
        raise BroydenError(
            f"{arguments.history_dir} does not exist or is not a directory; "
            "Broyden mixing requires track=true"
        )

    iteration = arguments.iteration
    if iteration is None:
        try:
            iteration = int(command_output(["getiter"]))
        except ValueError as exc:
            raise BroydenError("getiter did not return an integer") from exc
    if iteration < 1:
        raise BroydenError(f"Invalid DMFT iteration {iteration}")

    alpha = None
    if arguments.delta:
        alpha = get_float_parameter("alpha")
        if not 0.0 < alpha <= 1.0:
            raise BroydenError(f"alpha must satisfy 0 < alpha <= 1, got {alpha}")

    under = None
    goal = None
    maxdx = None
    if arguments.mu:
        under = get_float_parameter("under", 0.5)
        goal = get_float_parameter("goal")
        maxdx = get_float_parameter("maxdx")
        if under <= 0.0:
            raise BroydenError(f"under must be positive, got {under}")
        if maxdx <= 0.0:
            raise BroydenError(f"maxdx must be positive, got {maxdx}")

    mesh = None
    if arguments.delta:
        raw_re_source = arguments.workdir / "ReDelta.dat.NEW-common"
        raw_im_source = arguments.workdir / "ImDelta.dat.NEW-common"
        raw_re_history = raw_re_path(arguments.history_dir, iteration)
        raw_im_history = raw_im_path(arguments.history_dir, iteration)
        atomic_copy(raw_re_source, raw_re_history)
        atomic_copy(raw_im_source, raw_im_history)
        mesh, _ = load_table(raw_re_history)
        im_mesh, _ = load_table(raw_im_history)
        old_mesh, _ = load_table(
            arguments.workdir / "ReDelta.dat.OLD-common"
        )
        old_im_mesh, _ = load_table(
            arguments.workdir / "ImDelta.dat.OLD-common"
        )
        if not same_mesh(mesh, im_mesh) or not same_mesh(mesh, old_mesh):
            raise BroydenError("ReDelta/ImDelta OLD/NEW common meshes do not agree")
        if not same_mesh(mesh, old_im_mesh):
            raise BroydenError("ReDelta/ImDelta OLD/NEW common meshes do not agree")

    maximum = get_max_history()
    indices = history_indices(
        arguments.history_dir,
        iteration,
        arguments.delta,
        arguments.mu,
        maximum,
    )
    inputs, residuals, current_occupancy = build_history(
        arguments.history_dir,
        indices,
        mesh,
        arguments.delta,
        arguments.mu,
        goal,
    )

    delta_size = 2 * len(mesh) if arguments.delta else 0
    initial_factors = np.empty(inputs.shape[1])
    if arguments.delta:
        initial_factors[:delta_size] = alpha
    if arguments.mu:
        initial_factors[delta_size] = under

    mixed, pair_count, condition = broyden_step(
        inputs, residuals, initial_factors
    )
    if not np.all(np.isfinite(mixed)):
        raise BroydenError("Broyden update produced a non-finite value")

    if arguments.delta:
        write_table(
            arguments.workdir / "ReDelta.dat.TEMP",
            mesh,
            mixed[: len(mesh)],
        )
        write_table(
            arguments.workdir / "ImDelta.dat.TEMP",
            mesh,
            mixed[len(mesh) : delta_size],
        )

    if arguments.mu:
        old_mu = inputs[-1, delta_size]
        current_mu = load_scalar(arguments.mu_input)
        if not np.isclose(old_mu, current_mu, rtol=1e-11, atol=1e-13):
            raise BroydenError(
                f"{arguments.mu_input} changed after occupancy measurement: "
                f"{old_mu} != {current_mu}"
            )
        requested_dx = mixed[delta_size] - old_mu
        dx = float(np.clip(requested_dx, -maxdx, maxdx))
        if dx != requested_dx:
            print(f"Clipping mu step from {requested_dx:.17g} to {dx:.17g}")
        new_mu = old_mu + dx
        write_scalar(arguments.mu_output, new_mu)
        if not arguments.log.parent.is_dir():
            raise BroydenError(
                f"Output directory for {arguments.log} does not exist: "
                f"{arguments.log.parent}"
            )
        with arguments.log.open("a", encoding="ascii") as output:
            output.write(
                f"{old_mu:.17g} {new_mu:.17g} {dx:.17g} "
                f"{current_occupancy:.17g}\n"
            )
        print(
            f"mu: old={old_mu:.17g} new={new_mu:.17g} "
            f"occupancy={current_occupancy:.17g} goal={goal:.17g}"
        )

    print(
        f"Broyden iteration={iteration} history={indices[0]}..{indices[-1]} "
        f"pairs={pair_count} residual={np.linalg.norm(residuals[-1]):.8g} "
        f"condition={condition:.8g}"
    )


if __name__ == "__main__":
    try:
        main()
    except (BroydenError, OSError, ValueError, np.linalg.LinAlgError) as exc:
        print(f"broyden.py: {exc}", file=sys.stderr)
        sys.exit(1)
