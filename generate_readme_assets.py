"""Generate deterministic README figures from the quick model and FEMM outputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import sys
import types
from typing import Any, Iterable

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Wedge

import design_pmsm
from motor_config import (
    ConfigError,
    MotorConfig,
    canonical_config_sha256,
    load_motor_config,
    motor_config_from_dict,
    physical_model_fingerprint,
)


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_ACCEPTANCE_DIR = REPO_ROOT / "output_femm_fix_smoke"
DEFAULT_FEMM_IMAGE_DIR = REPO_ROOT / "output_readme_femm"
DEFAULT_ASSET_DIR = REPO_ROOT / "docs" / "assets"
PNG_METADATA = {"Software": "motor-design generate_readme_assets.py"}

PHASE_COLORS = {"A": "#C44E52", "B": "#4C72B0", "C": "#55A868"}
MAGNET_COLORS = {"N": "#DD8452", "S": "#8172B3"}


def _configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "figure.titlesize": 13,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return resolved.as_posix()


def _require_file(path: Path, description: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {description}: {path}")
    return path


def _load_resolved_config(path: Path, description: str) -> tuple[str, str]:
    _require_file(path, description)
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid resolved configuration JSON {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"Resolved configuration must be a JSON object: {path}")
    declared_hash = payload.get("shared_config_sha256")
    if not isinstance(declared_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", declared_hash):
        raise ValueError(
            f"{path} does not contain a valid shared_config_sha256 value."
        )
    raw_config = payload.get("shared_config")
    if not isinstance(raw_config, dict):
        raise ValueError(f"{path} does not contain a shared_config object.")
    try:
        resolved_config = motor_config_from_dict(raw_config)
    except (ConfigError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid shared_config in {path}: {exc}") from exc
    computed_hash = canonical_config_sha256(resolved_config)
    if declared_hash != computed_hash:
        raise ValueError(
            f"Resolved configuration hash mismatch in {path}: "
            f"declared {declared_hash}, computed {computed_hash}."
        )
    physical_hash = physical_model_fingerprint(
        resolved_config.machine,
        resolved_config.materials,
        resolved_config.conventions,
    )
    declared_physical = payload.get("physical_model_sha256")
    if declared_physical is not None and declared_physical != physical_hash:
        raise ValueError(
            f"Resolved physical-model hash mismatch in {path}: "
            f"declared {declared_physical}, computed {physical_hash}."
        )
    return computed_hash, physical_hash


def _require_matching_resolved_config(
    path: Path,
    *,
    expected_physical_hash: str,
    description: str,
) -> Path:
    _, actual_physical_hash = _load_resolved_config(path, description)
    if actual_physical_hash != expected_physical_hash:
        raise ValueError(
            f"{description} was generated from a different FEMM physical model: "
            f"expected {expected_physical_hash}, found {actual_physical_hash} "
            f"in {path}."
        )
    return path


def _manifest_value_matches(actual: Any, expected: Any) -> bool:
    if isinstance(expected, float):
        try:
            return math.isclose(float(actual), expected, rel_tol=0.0, abs_tol=1e-12)
        except (TypeError, ValueError):
            return False
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(
            _manifest_value_matches(actual_item, expected_item)
            for actual_item, expected_item in zip(actual, expected)
        )
    return actual == expected


def _require_provenanced_file(
    path: Path,
    expected_physical_hash: str,
    *,
    required_analyses: Iterable[str] = (),
    expected_arguments: dict[str, Any] | None = None,
) -> Path:
    """Require an exact file hash in a successful FEMM run manifest."""
    _require_file(path, "FEMM result")
    manifest_path = path.parent / "femm_run_manifest.json"
    _require_file(manifest_path, f"run manifest for {path.name}")
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid FEMM run manifest {manifest_path}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError(f"Unsupported FEMM run manifest: {manifest_path}")
    runs = payload.get("runs")
    if not isinstance(runs, list):
        raise ValueError(f"Invalid run list in FEMM manifest: {manifest_path}")

    expected_file_hash = _sha256(path)
    relative = path.relative_to(path.parent).as_posix()
    required_analysis_set = set(required_analyses)
    expected_arguments = {} if expected_arguments is None else expected_arguments
    for run in runs:
        if not isinstance(run, dict):
            continue
        if run.get("physical_model_sha256") != expected_physical_hash:
            continue
        analyses = run.get("analyses")
        if not isinstance(analyses, list) or not required_analysis_set <= set(analyses):
            continue
        arguments = run.get("arguments")
        if not isinstance(arguments, dict) or any(
            not _manifest_value_matches(arguments.get(name), expected)
            for name, expected in expected_arguments.items()
        ):
            continue
        declared_out = arguments.get("out")
        if not isinstance(declared_out, str):
            continue
        declared_out_path = Path(declared_out)
        if not declared_out_path.is_absolute():
            declared_out_path = REPO_ROOT / declared_out_path
        if declared_out_path.resolve() != path.parent.resolve():
            continue
        files = run.get("files")
        if not isinstance(files, list):
            continue
        for item in files:
            if not isinstance(item, dict):
                continue
            if item.get("path") == relative and item.get("sha256") == expected_file_hash:
                return manifest_path
    raise ValueError(
        f"No successful FEMM run in {manifest_path} binds {path.name} "
        f"(sha256 {expected_file_hash}) to physical model {expected_physical_hash} "
        f"and the required analysis contract."
    )


def _read_csv_rows(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    _require_file(path, "CSV input")
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or ())
        missing = sorted(set(required) - fields)
        if missing:
            raise ValueError(f"{path} is missing column(s): {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"CSV input has no data rows: {path}")
    return rows


def _float_column(rows: list[dict[str, str]], name: str) -> list[float]:
    try:
        values = [float(row[name]) for row in rows]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Column {name!r} contains a non-numeric value.") from exc
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"Column {name!r} contains a non-finite value.")
    return values


def _save_figure(fig: matplotlib.figure.Figure, destination: Path, dpi: int) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.stem}.tmp.png")
    try:
        fig.savefig(
            temporary,
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
            metadata=PNG_METADATA,
        )
        temporary.replace(destination)
    finally:
        plt.close(fig)
        if temporary.exists():
            temporary.unlink()
    return destination


def _png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Expected a PNG image: {path}")
    return struct.unpack(">II", header[16:24])


def _asset_record(
    path: Path,
    *,
    kind: str,
    sources: Iterable[Path],
) -> dict[str, Any]:
    width, height = _png_dimensions(path)
    return {
        "name": path.name,
        "path": _display_path(path),
        "kind": kind,
        "sha256": _sha256(path),
        "bytes": path.stat().st_size,
        "width_px": width,
        "height_px": height,
        "sources": [_display_path(source) for source in sources],
    }


def _source_record(path: Path, role: str) -> dict[str, Any]:
    return {
        "path": _display_path(path),
        "role": role,
        "sha256": _sha256(path),
        "bytes": path.stat().st_size,
    }


def generate_quick_sweep(
    config: MotorConfig,
    asset_dir: Path,
) -> Path:
    result = design_pmsm.run_quick_analysis(config)
    destination = design_pmsm.plot_results(
        result,
        asset_dir,
        save_png=True,
        show=False,
    )
    if destination is None:
        raise RuntimeError("The quick model did not return a PNG path.")
    return _require_file(destination, "quick-model plot")


def _install_geometry_import_stubs() -> None:
    """Import FEMM geometry helpers without opening or requiring FEMM COM."""
    if importlib.util.find_spec("femm") is None and "femm" not in sys.modules:
        sys.modules["femm"] = types.ModuleType("femm")
    if importlib.util.find_spec("tqdm") is None and "tqdm" not in sys.modules:
        module = types.ModuleType("tqdm")
        module.tqdm = lambda iterable=None, *args, **kwargs: iterable
        sys.modules["tqdm"] = module


def generate_winding_layout(
    config: MotorConfig,
    destination: Path,
    dpi: int,
) -> Path:
    _install_geometry_import_stubs()
    import femm_spm_template as femm_model

    machine = femm_model.machine_from_motor_config(config)
    layout = femm_model.winding_layout(machine)
    if len(layout) != machine.slots:
        raise ValueError("FEMM winding layout length does not match the slot count.")

    fig, ax = plt.subplots(figsize=(8.2, 8.2), constrained_layout=True)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Stator steel, rotor steel, and nonmagnetic bore establish the radial proportions.
    ax.add_patch(
        Wedge(
            (0.0, 0.0),
            machine.r_stator_outer_mm,
            0.0,
            360.0,
            width=machine.r_stator_outer_mm - machine.r_stator_inner_mm,
            facecolor="#D9D9D9",
            edgecolor="#555555",
            linewidth=1.0,
        )
    )
    ax.add_patch(
        Wedge(
            (0.0, 0.0),
            machine.r_rotor_mm,
            0.0,
            360.0,
            width=machine.r_rotor_mm - machine.r_shaft_mm,
            facecolor="#A7A7A7",
            edgecolor="#4A4A4A",
            linewidth=1.0,
        )
    )
    ax.add_patch(
        Wedge(
            (0.0, 0.0),
            machine.r_shaft_mm,
            0.0,
            360.0,
            facecolor="white",
            edgecolor="#4A4A4A",
            linewidth=1.0,
        )
    )

    pole_count = 2 * machine.pole_pairs
    pole_pitch_deg = 360.0 / pole_count
    magnet_span_deg = machine.magnet_arc_ratio * pole_pitch_deg
    for pole in range(pole_count):
        center = pole * pole_pitch_deg
        polarity = "N" if pole % 2 == 0 else "S"
        ax.add_patch(
            Wedge(
                (0.0, 0.0),
                machine.r_mag_outer_mm,
                center - 0.5 * magnet_span_deg,
                center + 0.5 * magnet_span_deg,
                width=machine.mag_thickness_mm,
                facecolor=MAGNET_COLORS[polarity],
                edgecolor="white",
                linewidth=0.45,
            )
        )

    slot_pitch_deg = 360.0 / machine.slots
    slot_span_deg = 0.55 * slot_pitch_deg
    slot_inner = machine.r_stator_inner_mm + 0.2
    slot_outer = machine.r_slot_outer_mm
    label_radius = 0.5 * (slot_inner + slot_outer)
    for slot, (phase, sign) in enumerate(layout):
        center = slot * slot_pitch_deg
        hatch = "///" if sign > 0 else "\\\\\\"
        ax.add_patch(
            Wedge(
                (0.0, 0.0),
                slot_outer,
                center - 0.5 * slot_span_deg,
                center + 0.5 * slot_span_deg,
                width=slot_outer - slot_inner,
                facecolor=PHASE_COLORS[phase],
                edgecolor="white",
                linewidth=0.65,
                hatch=hatch,
                alpha=0.92,
            )
        )
        angle = math.radians(center)
        ax.text(
            label_radius * math.cos(angle),
            label_radius * math.sin(angle),
            f"{slot + 1}\n{phase}{'+' if sign > 0 else '-'}",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            fontweight="bold",
        )

    limit = machine.r_stator_outer_mm * 1.07
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_title(
        "18-slot / 20-pole single-layer winding\n"
        f"{machine.turns_per_slot} turns/slot | "
        f"{config.machine.turns_per_phase} series turns/phase",
        pad=12,
    )
    legend = [
        Patch(facecolor=PHASE_COLORS[phase], label=f"Phase {phase}")
        for phase in ("A", "B", "C")
    ]
    legend.extend(
        Patch(facecolor=MAGNET_COLORS[polarity], label=f"Magnet {polarity}")
        for polarity in ("N", "S")
    )
    legend.extend(
        (
            Patch(facecolor="#D9D9D9", edgecolor="#555555", label="Stator steel"),
            Patch(facecolor="#A7A7A7", edgecolor="#4A4A4A", label="Rotor steel"),
            Patch(facecolor="white", edgecolor="#4A4A4A", label="Nonmagnetic bore"),
        )
    )
    ax.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.04),
        ncol=4,
        frameon=False,
    )
    return _save_figure(fig, destination, dpi)


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in {"true", "false"}:
        return value.strip().lower() == "true"
    raise ValueError(f"Expected a boolean acceptance value, got {value!r}.")


def _normalize_metrics(raw: Any) -> dict[str, dict[str, Any]]:
    if isinstance(raw, dict) and "metrics" in raw:
        raw = raw["metrics"]
    if isinstance(raw, list):
        normalized: dict[str, dict[str, Any]] = {}
        for item in raw:
            if not isinstance(item, dict) or "metric" not in item:
                raise ValueError("Acceptance metric lists require a 'metric' field.")
            normalized[str(item["metric"])] = dict(item)
        return normalized
    if isinstance(raw, dict) and all(isinstance(value, dict) for value in raw.values()):
        return {str(name): dict(value) for name, value in raw.items()}
    raise ValueError("Unsupported electromagnetic acceptance metrics JSON structure.")


def load_acceptance_metrics(path: Path) -> dict[str, dict[str, Any]]:
    _require_file(path, "electromagnetic acceptance metrics")
    if path.suffix.lower() == ".json":
        try:
            raw = json.loads(
                path.read_text(encoding="utf-8-sig"),
                parse_constant=lambda value: (_ for _ in ()).throw(
                    ValueError(f"Non-finite JSON constant: {value}")
                ),
            )
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid acceptance JSON {path}: {exc}") from exc
        metrics = _normalize_metrics(raw)
    else:
        rows = _read_csv_rows(
            path,
            ("metric", "value", "criterion", "limit", "passed"),
        )
        metrics = {row["metric"]: dict(row) for row in rows}

    for name, metric in metrics.items():
        if not {"value", "criterion", "limit", "passed"} <= set(metric):
            raise ValueError(f"Acceptance metric {name!r} is incomplete.")
        metric["value"] = float(metric["value"])
        limit = metric["limit"]
        metric["limit"] = math.nan if limit is None else float(limit)
        metric["criterion"] = str(metric["criterion"])
        metric["passed"] = _parse_bool(metric["passed"])
    return metrics


def generate_acceptance_plot(
    validation_csv: Path,
    metrics_path: Path,
    destination: Path,
    dpi: int,
) -> Path:
    columns = (
        "rotor_angle_deg",
        "psi_d0_wb",
        "psi_q0_wb",
        "cogging_torque_nm",
        "torque_q_pos_nm",
        "torque_q_neg_nm",
        "Ld_H",
        "Lq_H",
    )
    rows = _read_csv_rows(validation_csv, columns)
    metrics = load_acceptance_metrics(metrics_path)
    angle = _float_column(rows, "rotor_angle_deg")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    ax = axes[0, 0]
    psi_d_mwb = [
        1000.0 * value for value in _float_column(rows, "psi_d0_wb")
    ]
    psi_q_mwb = [
        1000.0 * value for value in _float_column(rows, "psi_q0_wb")
    ]
    line_d = ax.plot(angle, psi_d_mwb, "o-", label="Psi d", color="#4C72B0")
    ax.set_title("No-load dq flux linkage")
    ax.set_xlabel("Rotor mechanical angle [deg]")
    ax.set_ylabel("Psi d [mWb]", color="#4C72B0")
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax_q = ax.twinx()
    line_q = ax_q.plot(
        angle,
        psi_q_mwb,
        "s--",
        label="Psi q",
        color="#DD8452",
    )
    ax_q.set_ylabel("Psi q [mWb]", color="#DD8452")
    ax_q.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.legend(line_d + line_q, ["Psi d", "Psi q"], loc="best")

    ax = axes[0, 1]
    cogging = _float_column(rows, "cogging_torque_nm")
    positive_torque = [
        total - cogging_value
        for total, cogging_value in zip(
            _float_column(rows, "torque_q_pos_nm"), cogging
        )
    ]
    negative_torque_magnitude = [
        -(total - cogging_value)
        for total, cogging_value in zip(
            _float_column(rows, "torque_q_neg_nm"), cogging
        )
    ]
    ax.plot(angle, positive_torque, "o-", label="+Iq electromagnetic torque")
    ax.plot(
        angle,
        negative_torque_magnitude,
        "s--",
        label="Magnitude at -Iq",
    )
    ax.set_title("Cogging-subtracted torque symmetry")
    ax.set_xlabel("Rotor mechanical angle [deg]")
    ax.set_ylabel("Torque [N m]")
    ax.legend()

    ax = axes[1, 0]
    ax.plot(angle, [1000.0 * value for value in _float_column(rows, "Ld_H")], "o-", label="Ld")
    ax.plot(angle, [1000.0 * value for value in _float_column(rows, "Lq_H")], "s-", label="Lq")
    ax.set_title("Centered incremental inductance")
    ax.set_xlabel("Rotor mechanical angle [deg]")
    ax.set_ylabel("Inductance [mH]")
    ax.legend()

    gated = [
        (name, metric)
        for name, metric in metrics.items()
        if metric["criterion"] != "report"
    ]
    if not gated:
        raise ValueError("Acceptance metrics contain no gated checks.")
    table_rows = []
    row_colors = []
    for name, metric in gated:
        table_rows.append(
            [
                name.replace("_", " "),
                f"{metric['value']:.4g}",
                f"{metric['criterion']} {metric['limit']:.4g}",
                "PASS" if metric["passed"] else "FAIL",
            ]
        )
        status_color = "#DCEFD8" if metric["passed"] else "#F3D2D2"
        row_colors.append(["white", "white", "white", status_color])
    ax = axes[1, 1]
    ax.set_axis_off()
    ax.set_title("Electromagnetic acceptance gates")
    table = ax.table(
        cellText=table_rows,
        colLabels=["Metric", "Value", "Criterion", "Status"],
        cellColours=row_colors,
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.49, 0.15, 0.19, 0.12],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.3)
    table.scale(1.0, 1.35)

    passed = sum(bool(metric["passed"]) for _, metric in gated)
    fig.suptitle(
        f"Reduced FEMM electromagnetic regression ({len(rows)} positions) | "
        f"{passed}/{len(gated)} gates passed"
    )
    return _save_figure(fig, destination, dpi)


def generate_convergence_plot(
    convergence_csv: Path,
    integral_csv: Path,
    destination: Path,
    dpi: int,
) -> Path:
    summary_columns = (
        "study",
        "level",
        "loaded_mean_nm",
        "cogging_pp_nm",
        "loaded_mean_error",
        "cogging_pp_error",
    )
    summary = _read_csv_rows(convergence_csv, summary_columns)
    integral_columns = (
        "airgap_fraction",
        "points",
        "loaded_increment_maxwell_nm",
        "loaded_increment_wst_nm",
    )
    integral = _read_csv_rows(integral_csv, integral_columns)

    labels = [f"{row['study']}\n{row['level']}" for row in summary]
    positions = list(range(len(summary)))
    loaded = _float_column(summary, "loaded_mean_nm")
    cogging = _float_column(summary, "cogging_pp_nm")
    loaded_error = [100.0 * value for value in _float_column(summary, "loaded_mean_error")]
    cogging_error = [100.0 * value for value in _float_column(summary, "cogging_pp_error")]

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(positions, loaded, "o-", color="#4C72B0")
    ax.set_xticks(positions, labels)
    ax.set_title("Loaded mean torque")
    ax.set_ylabel("Torque [N m]")

    ax = axes[0, 1]
    ax.plot(positions, cogging, "o-", color="#C44E52")
    ax.set_xticks(positions, labels)
    ax.set_title("Cogging peak-to-peak torque")
    ax.set_ylabel("Torque [N m]")

    ax = axes[1, 0]
    ax.plot(positions, loaded_error, "o-", label="Loaded mean error")
    ax.plot(positions, cogging_error, "s-", label="Cogging p-p error")
    ax.axhline(2.0, color="#4C72B0", linestyle=":", label="Loaded criterion (2%)")
    ax.axhline(10.0, color="#C44E52", linestyle=":", label="Cogging criterion (10%)")
    ax.set_xticks(positions, labels)
    ax.set_title("Change relative to study reference")
    ax.set_ylabel("Relative change [%]")
    ax.legend(fontsize=7)

    ax = axes[1, 1]
    fractions = sorted({float(row["airgap_fraction"]) for row in integral})
    for fraction in fractions:
        selected = sorted(
            (row for row in integral if float(row["airgap_fraction"]) == fraction),
            key=lambda row: int(float(row["points"])),
        )
        ax.plot(
            [int(float(row["points"])) for row in selected],
            [float(row["loaded_increment_maxwell_nm"]) for row in selected],
            "o-",
            label=f"Maxwell, gap fraction {fraction:g}",
        )
    wst_values = _float_column(integral, "loaded_increment_wst_nm")
    ax.axhline(
        sum(wst_values) / len(wst_values),
        color="black",
        linestyle="--",
        linewidth=1.0,
        label="Weighted stress tensor",
    )
    ax.set_title("Air-gap torque integral convergence")
    ax.set_xlabel("Circumferential sample points")
    ax.set_ylabel("Loaded torque increment [N m]")
    ax.legend(fontsize=7)

    mesh_points = sorted(
        {int(float(row["points"])) for row in summary if row["study"] == "mesh"}
    )
    angle_points = sorted(
        {int(float(row["points"])) for row in summary if row["study"] == "angle"}
    )
    mesh_label = "/".join(str(value) for value in mesh_points)
    angle_label = "/".join(str(value) for value in angle_points)
    fig.suptitle(
        "Reduced FEMM convergence regression | "
        f"mesh: {mesh_label} positions, angle: {angle_label} positions"
    )
    return _save_figure(fig, destination, dpi)


FEMM_IMAGE_SPECS = (
    ("field_density", "field_density_*.png", "femm_field_density.png"),
    ("field_lines", "field_lines_*.png", "femm_flux_lines.png"),
    ("airgap_B", "airgap_B_*.png", "femm_airgap_flux_density.png"),
)


def _select_femm_image(
    image_dir: Path,
    prefix: str,
    pattern: str,
    case_tag: str | None,
) -> Path:
    if case_tag:
        source = image_dir / f"{prefix}_{case_tag}.png"
        return _require_file(source, f"FEMM {prefix} image")
    matches = sorted(image_dir.glob(pattern), key=lambda path: path.name)
    if not matches:
        raise FileNotFoundError(
            f"No FEMM image matching {pattern!r} in {image_dir}; "
            "generate the snapshot or use --skip-femm."
        )
    if len(matches) > 1:
        names = ", ".join(path.name for path in matches)
        raise ValueError(
            f"Multiple FEMM images match {pattern!r}: {names}. "
            "Select one reproducibly with --femm-case-tag."
        )
    return matches[0]


def copy_femm_images(
    image_dir: Path,
    asset_dir: Path,
    case_tag: str | None,
) -> list[tuple[Path, Path]]:
    if not image_dir.is_dir():
        raise FileNotFoundError(
            f"FEMM image directory does not exist: {image_dir}; use --skip-femm."
        )
    selected = [
        (
            _select_femm_image(image_dir, prefix, pattern, case_tag),
            asset_dir / stable_name,
        )
        for prefix, pattern, stable_name in FEMM_IMAGE_SPECS
    ]
    for source, _ in selected:
        _png_dimensions(source)
    copied: list[tuple[Path, Path]] = []
    for source, destination in selected:
        shutil.copyfile(source, destination)
        copied.append((destination, source))
    return copied


def _resolve_input(override: Path | None, directory: Path, filename: str) -> Path:
    return override if override is not None else directory / filename


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate stable README assets without opening FEMM or a GUI."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=REPO_ROOT / "motor_config.json",
        help="Shared motor configuration JSON.",
    )
    parser.add_argument(
        "--acceptance-dir",
        type=Path,
        default=DEFAULT_ACCEPTANCE_DIR,
        help="Directory containing FEMM validation and convergence tables.",
    )
    parser.add_argument("--validation-csv", type=Path)
    parser.add_argument("--acceptance-metrics", type=Path)
    parser.add_argument("--convergence-csv", type=Path)
    parser.add_argument("--integral-csv", type=Path)
    parser.add_argument(
        "--mesh-check-json",
        type=Path,
        default=(
            REPO_ROOT
            / "output_femm_meshcheck"
            / "electromagnetic_mesh_comparison_medium_vs_fine.json"
        ),
        help="Medium/fine global-mesh comparison used by the README.",
    )
    parser.add_argument(
        "--femm-image-dir",
        type=Path,
        default=DEFAULT_FEMM_IMAGE_DIR,
        help="Directory containing field_density, field_lines, and airgap_B PNGs.",
    )
    parser.add_argument(
        "--femm-case-tag",
        help=(
            "Case suffix such as rpm600_iq5_a0 when the FEMM directory has "
            "multiple cases."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_ASSET_DIR,
        help="Stable README asset directory.",
    )
    parser.add_argument(
        "--skip-femm",
        action="store_true",
        help="Generate analytical/acceptance assets without field and air-gap PNGs.",
    )
    parser.add_argument("--dpi", type=int, default=180, help="DPI for generated figures.")
    return parser


def generate_assets(args: argparse.Namespace) -> Path:
    if args.dpi < 72:
        raise ValueError("--dpi must be at least 72.")
    config_path = _require_file(args.config, "motor configuration")
    acceptance_dir = args.acceptance_dir
    validation_csv = _resolve_input(
        args.validation_csv,
        acceptance_dir,
        "electromagnetic_validation_samples.csv",
    )
    convergence_csv = _resolve_input(
        args.convergence_csv,
        acceptance_dir,
        "torque_convergence_summary.csv",
    )
    integral_csv = _resolve_input(
        args.integral_csv,
        acceptance_dir,
        "torque_integral_convergence.csv",
    )
    mesh_check_json = args.mesh_check_json

    warnings: list[str] = []
    if args.acceptance_metrics is not None:
        metrics_path = args.acceptance_metrics
    else:
        metrics_path = acceptance_dir / "electromagnetic_acceptance_metrics.json"
        if not metrics_path.is_file():
            legacy = acceptance_dir / "electromagnetic_acceptance.csv"
            _require_file(legacy, "acceptance metrics JSON or legacy CSV")
            metrics_path = legacy
            warnings.append(
                "electromagnetic_acceptance_metrics.json was absent; used legacy "
                "electromagnetic_acceptance.csv."
            )

    for path, label in (
        (validation_csv, "validation samples"),
        (metrics_path, "acceptance metrics"),
        (convergence_csv, "torque convergence summary"),
        (integral_csv, "torque integral convergence"),
        (mesh_check_json, "medium/fine global-mesh comparison"),
    ):
        _require_file(path, label)

    asset_dir = args.out_dir
    asset_dir.mkdir(parents=True, exist_ok=True)
    config = load_motor_config(config_path)
    expected_physical_hash = physical_model_fingerprint(
        config.machine,
        config.materials,
        config.conventions,
    )
    validation_contract = {
        "mesh_level": "medium",
        "validation_iq": 1.0,
        "validation_delta_current": 1.0,
        "validation_steps": 3,
    }
    convergence_contract = {
        "convergence_iq": 5.0,
        "convergence_mesh_points": 3,
        "convergence_angle_points": ["3,6"],
    }
    mesh_check_contract = {
        "mesh_check_reference_level": "medium",
        "mesh_check_level": "fine",
        "validation_iq": 1.0,
        "validation_delta_current": 1.0,
    }
    result_specs = (
        (validation_csv, {"validate"}, validation_contract),
        (metrics_path, {"validate"}, validation_contract),
        (convergence_csv, {"convergence"}, convergence_contract),
        (integral_csv, {"convergence"}, convergence_contract),
        (mesh_check_json, {"meshcheck"}, mesh_check_contract),
    )
    resolved_inputs: set[Path] = set()
    run_manifests: set[Path] = set()
    for result_path, required_analyses, expected_arguments in result_specs:
        resolved_path = result_path.parent / "resolved_femm_config.json"
        resolved_inputs.add(
            _require_matching_resolved_config(
                resolved_path,
                expected_physical_hash=expected_physical_hash,
                description=f"resolved configuration for {result_path.name}",
            )
        )
        run_manifests.add(
            _require_provenanced_file(
                result_path,
                expected_physical_hash,
                required_analyses=required_analyses,
                expected_arguments=expected_arguments,
            )
        )

    image_resolved: Path | None = None
    selected_femm_inputs: list[Path] = []
    snapshot_metadata: Path | None = None
    if not args.skip_femm:
        selected_femm_inputs = [
            _select_femm_image(
                args.femm_image_dir,
                prefix,
                pattern,
                args.femm_case_tag,
            )
            for prefix, pattern, _ in FEMM_IMAGE_SPECS
        ]
        if args.femm_case_tag:
            snapshot_metadata = _require_file(
                args.femm_image_dir
                / f"field_snapshot_{args.femm_case_tag}.csv",
                "FEMM snapshot metadata",
            )
            selected_femm_inputs.append(snapshot_metadata)
        image_resolved = _require_matching_resolved_config(
            args.femm_image_dir / "resolved_femm_config.json",
            expected_physical_hash=expected_physical_hash,
            description="FEMM field-snapshot resolved configuration",
        )
        snapshot_contract = {
            "mesh_level": "medium",
            "rpm_list": ["600"],
            "iq_list": ["5"],
            "snapshot_angle_deg": 0.0,
            "field_radial_points": 24,
            "field_angular_points": 120,
            "airgap_points": 360,
        }
        for index, result_path in enumerate(selected_femm_inputs):
            required_analyses = (
                {"field"}
                if index < 2
                else {"airgap"}
                if index == 2
                else {"field", "airgap"}
            )
            run_manifests.add(
                _require_provenanced_file(
                    result_path,
                    expected_physical_hash,
                    required_analyses=required_analyses,
                    expected_arguments=snapshot_contract,
                )
            )
    _configure_plot_style()

    assets: list[dict[str, Any]] = []
    quick_path = generate_quick_sweep(config, asset_dir)
    assets.append(
        _asset_record(
            quick_path,
            kind="generated:quick-model",
            sources=(config_path, REPO_ROOT / "design_pmsm.py"),
        )
    )

    winding_path = generate_winding_layout(
        config,
        asset_dir / "winding_18s20p.png",
        args.dpi,
    )
    assets.append(
        _asset_record(
            winding_path,
            kind="generated:femm-winding-layout",
            sources=(config_path, REPO_ROOT / "femm_spm_template.py"),
        )
    )

    acceptance_path = generate_acceptance_plot(
        validation_csv,
        metrics_path,
        asset_dir / "electromagnetic_acceptance.png",
        args.dpi,
    )
    assets.append(
        _asset_record(
            acceptance_path,
            kind="generated:femm-acceptance",
            sources=(validation_csv, metrics_path),
        )
    )

    convergence_path = generate_convergence_plot(
        convergence_csv,
        integral_csv,
        asset_dir / "torque_convergence.png",
        args.dpi,
    )
    assets.append(
        _asset_record(
            convergence_path,
            kind="generated:femm-convergence",
            sources=(convergence_csv, integral_csv),
        )
    )

    femm_sources: list[Path] = []
    if args.skip_femm:
        for _, _, stable_name in FEMM_IMAGE_SPECS:
            stale_snapshot = asset_dir / stable_name
            if stale_snapshot.is_file():
                stale_snapshot.unlink()
        warnings.append("FEMM field and air-gap PNG copying was skipped by --skip-femm.")
    else:
        for destination, source in copy_femm_images(
            args.femm_image_dir,
            asset_dir,
            args.femm_case_tag,
        ):
            femm_sources.append(source)
            assets.append(
                _asset_record(
                    destination,
                    kind="copied:femm-snapshot",
                    sources=(source,),
                )
            )
    source_roles = {
        config_path: "shared motor configuration",
        REPO_ROOT / "motor_config.py": "configuration schema",
        REPO_ROOT / "design_pmsm.py": "quick operating-point model",
        REPO_ROOT / "femm_spm_template.py": "FEMM geometry and winding definition",
        validation_csv: "electromagnetic validation samples",
        metrics_path: "electromagnetic acceptance metrics",
        convergence_csv: "torque convergence summary",
        integral_csv: "torque integral convergence",
        mesh_check_json: "medium/fine global-mesh sensitivity comparison",
        Path(__file__).resolve(): "asset generator",
    }
    for resolved_path in resolved_inputs:
        source_roles[resolved_path] = "FEMM result resolved configuration"
    for manifest_path in run_manifests:
        source_roles[manifest_path] = "per-file FEMM run provenance"
    if image_resolved is not None:
        source_roles[image_resolved] = "field-snapshot resolved configuration"
    source_roles.update({source: "FEMM rendered snapshot" for source in femm_sources})
    if snapshot_metadata is not None:
        source_roles[snapshot_metadata] = "FEMM snapshot operating point"
    sources = [
        _source_record(path, role)
        for path, role in sorted(
            source_roles.items(), key=lambda item: _display_path(item[0])
        )
    ]

    manifest = {
        "schema_version": 1,
        "generator": _display_path(Path(__file__)),
        "backend": "Agg",
        "custom_figure_dpi": args.dpi,
        "matplotlib_version": matplotlib.__version__,
        "femm_snapshots": {
            "skipped": bool(args.skip_femm),
            "source_directory": _display_path(args.femm_image_dir),
            "case_tag": args.femm_case_tag,
        },
        "warnings": warnings,
        "sources": sources,
        "assets": sorted(assets, key=lambda item: item["name"]),
    }
    manifest_path = asset_dir / "manifest.json"
    temporary_manifest = asset_dir / ".manifest.tmp.json"
    temporary_manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temporary_manifest.replace(manifest_path)
    return manifest_path


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        manifest = generate_assets(args)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(f"Generated README assets: {manifest.parent.resolve()}")
    print(f"Manifest: {manifest.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
