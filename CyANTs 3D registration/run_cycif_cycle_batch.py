#!/usr/bin/env python3
# Developed by Alex Wong
# Cite: Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in
# Preprint: https://doi.org/10.64898/2026.05.17.725158
# Registration method: ANTsX/ANTsPy - https://github.com/ANTsX/ANTsPy

"""Run a small CyCIF ROI registration batch across two cycles."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile

from ants_roi_quicksyn import default_registered_output, format_duration, progress_heartbeat
from extract_ims_roi_channels import parse_channels


def parse_optional_channels(value: str) -> list[int]:
    value = value.strip()
    if not value:
        return []
    return parse_channels(value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "One-command ROI workflow: refresh/apply Cycle1 DAPI transform, apply Cycle1 "
            "other channels, register Cycle2 DAPI, apply Cycle2 other channels, then write "
            "one multichannel OME-TIFF per cycle."
        )
    )
    parser.add_argument(
        "--project-root",
        "--proj",
        default="",
        help=(
            "Registration project root. If supplied, standard paths are inferred from "
            "reference/roi and cycles/cycle_###."
        ),
    )
    parser.add_argument(
        "--acquisition-root",
        "--acq",
        default="",
        help="Raw acquisition root containing Cycle1/Cycle2 folders. Default: parent of --project-root.",
    )
    parser.add_argument("--run-id", "--run", default="run_0001", help="Registration run id under each cycle. Default: run_0001")
    parser.add_argument("--cycle1-id", default="cycle_001", help="Cycle1 folder id under cycles/. Default: cycle_001")
    parser.add_argument("--cycle2-id", default="cycle_002", help="Cycle2 folder id under cycles/. Default: cycle_002")
    parser.add_argument("--fixed-crop", "--fixed", default="", help="Reference Cycle0 DAPI ROI crop NRRD")
    parser.add_argument("--channels", "--ch", required=True, type=parse_channels, help="Other channels to apply, for example 1-4")
    parser.add_argument("--spacing", default="0.711", help="Spacing override passed to .ims ROI extraction. Default: 0.711")
    parser.add_argument("--threads", "--t", type=int, default=32, help="CPU threads for ANTs/ITK. Default: 32")
    parser.add_argument(
        "--progress-interval",
        "--pi",
        type=float,
        default=30.0,
        help="Elapsed-time heartbeat interval for long operations. Use 0 to disable. Default: 30.",
    )
    parser.add_argument("--downsample-factor", "--ds", type=float, default=4.0, help="Cycle2 DAPI registration factor. Default: 4")
    parser.add_argument(
        "--intra-cycle-align-channels",
        "--intra-ch",
        type=parse_optional_channels,
        default=[],
        help=(
            "Channels to align to each cycle's DAPI crop before applying the cycle transform. "
            "Use for channel-specific lateral/chromatic offsets, for example 1,3."
        ),
    )
    parser.add_argument(
        "--pre-crop-intra-cycle-align-channels",
        "--precrop-intra-ch",
        "--precrop-ch",
        type=parse_optional_channels,
        default=[],
        help=(
            "Channels to align to each cycle's DAPI from padded .ims crops before recropping "
            "to the ROI, for example 2,3."
        ),
    )
    parser.add_argument(
        "--intra-cycle-transform",
        "--intra-tx",
        default="TRSAA",
        help="ANTsPy transform for intra-cycle channel-to-DAPI alignment. Default: TRSAA.",
    )
    parser.add_argument(
        "--intra-cycle-downsample-factor",
        "--intra-ds",
        type=float,
        default=8.0,
        help="Downsample factor for intra-cycle channel-to-DAPI alignment. Default: 8.",
    )
    parser.add_argument(
        "--intra-padding-xy",
        "--pad",
        type=int,
        default=256,
        help="XY padding in voxels for pre-crop intra-cycle alignment. Default: 256.",
    )
    parser.add_argument(
        "--save-pre-crop-aligned-roi",
        "--save-precrop",
        action="store_true",
        help="Save recropped pre-crop intra-cycle aligned ROI intermediates in the cycle ROI folder.",
    )
    parser.add_argument(
        "--keep-pre-crop-temp",
        "--keep-precrop-temp",
        action="store_true",
        help="Keep temporary padded DAPI/channel NRRDs used for pre-crop intra-cycle alignment.",
    )
    parser.add_argument("--open-qc", action="store_true", help="Open DAPI QC overlay PNGs after writing them")
    parser.add_argument(
        "--output-pixeltype",
        choices=["uint16", "float32"],
        default="uint16",
        help="Registered output pixel type. Default: uint16",
    )
    parser.add_argument(
        "--registered-format",
        "--reg-format",
        choices=["nrrd", "tif", "both"],
        default="nrrd",
        help="Registered channel output format for non-DAPI channels. Default: nrrd.",
    )
    parser.add_argument(
        "--uint16-scaling",
        "--u16",
        choices=["minmax", "robust", "clip"],
        default="minmax",
        help="How to convert transformed float images to uint16. Default: minmax",
    )
    parser.add_argument(
        "--tiff-compression",
        default="",
        help="Optional TIFF compression, for example zlib. Default: uncompressed for speed.",
    )
    parser.add_argument("--tmp-dir", default="", help="Temporary directory for multichannel TIFF staging")
    parser.add_argument(
        "--skip-ome",
        action="store_true",
        help="Do not build aggregate OME-TIFFs; useful when converting per-channel TIFF stacks to .ims.",
    )
    parser.add_argument(
        "--registered-tiff-manifest",
        default="",
        help="JSON manifest for registered per-channel TIFF stacks. Default: <project-root>/registered_tiff_manifest.json",
    )
    parser.add_argument(
        "--imaris-wrapper-cmd",
        default="",
        help=(
            "Optional command template to run after TIFF stacks are written. Placeholders: "
            "{manifest}, {output}, {project}, {cycle1_output}, {cycle2_output}. "
            "Example: matlab -batch \"XR_imaris_conversion_data_wrapper('{manifest}','{output}')\""
        ),
    )
    parser.add_argument(
        "--imaris-output",
        default="",
        help="Output .ims path placeholder for --imaris-wrapper-cmd. Default: <project-root>/registered_all_cycles.ims",
    )
    parser.add_argument(
        "--xr-imaris",
        action="store_true",
        help="Generate and run a MATLAB XR_imaris_conversion_data_wrapper call from registered TIFF outputs.",
    )
    parser.add_argument(
        "--xr-output-path",
        default="",
        help="Output path passed to XR_imaris_conversion_data_wrapper. Default: <project-root>/xr_imaris",
    )
    parser.add_argument(
        "--xr-pixel-sizes",
        default="",
        help="Pixel sizes passed to XR_imaris_conversion_data_wrapper as x,y,z. Default: --spacing.",
    )
    parser.add_argument(
        "--matlab",
        default="matlab",
        help="MATLAB executable for --xr-imaris. Default: matlab.",
    )
    parser.add_argument(
        "--matlab-start",
        action="store_true",
        help="Launch MATLAB using Windows start / matlab -r and do not wait. Default uses matlab -batch and waits.",
    )

    parser.add_argument(
        "--cycle1-source-ims",
        "--c1-src",
        default="",
        help=(
            "Original Cycle1 source .ims path, usually ch0. The runner searches from the "
            "containing Cycle1 folder to find other channels."
        ),
    )
    parser.add_argument(
        "--cycle1-ims-template",
        "--c1-template",
        default="",
        help="Cycle1 source path template with {channel}. Paths can be .ims or .tif/.tiff.",
    )
    parser.add_argument(
        "--cycle1-ims-map",
        "--c1-map",
        "--cycle1-source-map",
        "--c1-source-map",
        default="",
        help="Cycle1 explicit channel=path source map separated by semicolons. Paths can be .ims or .tif/.tiff.",
    )
    parser.add_argument(
        "--cycle1-ims-override",
        "--c1-ims",
        "--c1-source",
        action="append",
        default=[],
        help=(
            "Cycle1 single channel source override as channel=path. Can be repeated. "
            "The path can be .ims or a corrected full-volume .tif/.tiff stack."
        ),
    )
    parser.add_argument("--cycle1-roi", "--c1-roi", default="", help="Cycle1 ImageJ ROI file")
    parser.add_argument("--cycle1-roi-csv", "--c1-csv", default="", help="Cycle1 XY coordinate CSV ROI file")
    parser.add_argument(
        "--cycle1-intra-cycle-align-channels",
        "--c1-intra-ch",
        type=parse_optional_channels,
        default=None,
        help="Cycle1-specific intra-cycle alignment channels. Defaults to --intra-ch.",
    )
    parser.add_argument(
        "--cycle1-pre-crop-intra-cycle-align-channels",
        "--c1-precrop-intra-ch",
        "--c1-pre",
        type=parse_optional_channels,
        default=None,
        help="Cycle1-specific pre-crop intra-cycle alignment channels. Defaults to --precrop-intra-ch.",
    )
    parser.add_argument("--cycle1-moving-crop", "--c1-crop", default="", help="Cycle1 DAPI ROI crop NRRD")
    parser.add_argument("--cycle1-transform-dir", "--c1-tx", default="", help="Cycle1 existing transform directory")
    parser.add_argument("--cycle1-raw-output-dir", "--c1-raw", default="", help="Cycle1 raw ROI NRRD output directory")
    parser.add_argument(
        "--cycle1-registered-output-dir",
        default="",
        help="Cycle1 registered output directory. Default: <cycle1 run>/outputs",
    )
    parser.add_argument(
        "--cycle1-tiff-output",
        default="",
        help="Cycle1 multichannel OME-TIFF output. Default: <cycle1 run>/outputs/cycle_001_registered_multichannel.ome.tif",
    )
    parser.add_argument(
        "--skip-cycle1-ch0-apply",
        action="store_true",
        help="Do not refresh Cycle1 registered ch0 with the existing transform.",
    )

    parser.add_argument(
        "--cycle2-source-ims",
        "--c2-src",
        default="",
        help=(
            "Original Cycle2 source .ims path, usually ch0. The runner searches from the "
            "containing Cycle2 folder to find other channels."
        ),
    )
    parser.add_argument(
        "--cycle2-ims-template",
        "--c2-template",
        default="",
        help="Cycle2 source path template with {channel}. Paths can be .ims or .tif/.tiff.",
    )
    parser.add_argument(
        "--cycle2-ims-map",
        "--c2-map",
        "--cycle2-source-map",
        "--c2-source-map",
        default="",
        help="Cycle2 explicit channel=path source map separated by semicolons. Paths can be .ims or .tif/.tiff.",
    )
    parser.add_argument(
        "--cycle2-ims-override",
        "--c2-ims",
        "--c2-source",
        action="append",
        default=[],
        help=(
            "Cycle2 single channel source override as channel=path. Can be repeated. "
            "The path can be .ims or a corrected full-volume .tif/.tiff stack."
        ),
    )
    parser.add_argument("--cycle2-roi", "--c2-roi", default="", help="Cycle2 ImageJ ROI file")
    parser.add_argument("--cycle2-roi-csv", "--c2-csv", default="", help="Cycle2 XY coordinate CSV ROI file")
    parser.add_argument(
        "--cycle2-intra-cycle-align-channels",
        "--c2-intra-ch",
        type=parse_optional_channels,
        default=None,
        help="Cycle2-specific intra-cycle alignment channels. Defaults to --intra-ch.",
    )
    parser.add_argument(
        "--cycle2-pre-crop-intra-cycle-align-channels",
        "--c2-precrop-intra-ch",
        "--c2-pre",
        type=parse_optional_channels,
        default=None,
        help="Cycle2-specific pre-crop intra-cycle alignment channels. Defaults to --precrop-intra-ch.",
    )
    parser.add_argument("--cycle2-moving-crop", "--c2-crop", default="", help="Cycle2 DAPI ROI crop NRRD")
    parser.add_argument("--cycle2-transform-dir", "--c2-tx", default="", help="Cycle2 transform directory")
    parser.add_argument("--cycle2-raw-output-dir", "--c2-raw", default="", help="Cycle2 raw ROI NRRD output directory")
    parser.add_argument(
        "--cycle2-registered-output-dir",
        default="",
        help="Cycle2 registered output directory. Default: <cycle2 run>/outputs",
    )
    parser.add_argument(
        "--cycle2-tiff-output",
        default="",
        help="Cycle2 multichannel OME-TIFF output. Default: <cycle2 run>/outputs/cycle_002_registered_multichannel.ome.tif",
    )
    parser.add_argument(
        "--skip-cycle2-dapi-registration",
        action="store_true",
        help="Skip Cycle2 DAPI registration and only use the existing transform files.",
    )
    return parser


def registration_run_dir(transform_dir: Path) -> Path:
    transform_name = transform_dir.name.lower()
    if transform_name == "transforms" or transform_name.endswith("_transforms"):
        return transform_dir.parent
    return transform_dir


def record_timing(timings: list[dict], stage: str, start: float, **extra) -> float:
    elapsed = time.time() - start
    entry = {"stage": stage, "seconds": elapsed, "duration": format_duration(elapsed)}
    entry.update(extra)
    timings.append(entry)
    print(f"[timing] {stage}: {format_duration(elapsed)}", flush=True)
    return elapsed


def append_timing(timings: list[dict], stage: str, elapsed: float, **extra) -> None:
    entry = {"stage": stage, "seconds": elapsed, "duration": format_duration(elapsed)}
    entry.update(extra)
    timings.append(entry)
    print(f"[timing] {stage}: {format_duration(elapsed)}", flush=True)


def run_step(label: str, command: list[str]) -> float:
    print("", flush=True)
    print(f"===== {label} =====", flush=True)
    print(" ".join(f'"{part}"' if " " in part else part for part in command), flush=True)
    start = time.time()
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed with exit code {completed.returncode}")
    elapsed = time.time() - start
    print(f"===== finished {label} in {format_duration(elapsed)} =====", flush=True)
    return elapsed


def run_shell_step(label: str, command: str) -> float:
    print("", flush=True)
    print(f"===== {label} =====", flush=True)
    print(command, flush=True)
    start = time.time()
    completed = subprocess.run(command, shell=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed with exit code {completed.returncode}")
    elapsed = time.time() - start
    print(f"===== finished {label} in {format_duration(elapsed)} =====", flush=True)
    return elapsed


def output_dir_or_default(value: str, transform_dir: Path) -> Path:
    return Path(value).expanduser().resolve() if value else (registration_run_dir(transform_dir) / "outputs").resolve()


def tiff_output_or_default(value: str, transform_dir: Path, cycle_name: str) -> Path:
    if value:
        return Path(value).expanduser().resolve()
    return output_dir_or_default("", transform_dir) / f"{cycle_name}_registered_multichannel.ome.tif"


def imaris_output_or_default(value: str, project_root: Path) -> Path:
    if value:
        return Path(value).expanduser().resolve()
    return project_root / "registered_all_cycles.ims"


def script_path(name: str) -> str:
    return str((Path(__file__).resolve().parent / name).resolve())


def require_value(value: str, label: str) -> str:
    if not value:
        raise ValueError(f"{label} is required unless it can be inferred from --project-root")
    return value


def one_match(paths: list[Path], label: str) -> Path:
    existing = [p for p in paths if p.exists()]
    if len(existing) == 1:
        return existing[0]
    if not existing:
        raise FileNotFoundError(f"Could not infer {label}; no matching files found")
    raise ValueError(f"Could not infer {label}; multiple matches found: {existing}")


def first_glob(root: Path, patterns: list[str], label: str, exclude_registered: bool = True) -> Path:
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(sorted(root.glob(pattern)))
    if exclude_registered:
        matches = [p for p in matches if "registered" not in p.name.lower()]
    return one_match(matches, label)


def first_dapi_crop(root: Path, label: str) -> Path:
    matches: list[Path] = []
    for pattern in ("*level0-1.nrrd", "*level0*.nrrd", "*ch0*.nrrd"):
        matches.extend(sorted(root.glob(pattern)))
    filtered = []
    seen = set()
    for path in matches:
        key = str(path.resolve()).lower()
        if key in seen:
            continue
        seen.add(key)
        name = path.name.lower()
        if "registered" in name:
            continue
        if name.endswith("_roi.nrrd") or "_roi_" in name:
            continue
        filtered.append(path)
    return one_match(filtered, label)


def first_existing_glob(root: Path, patterns: list[str], exclude_registered: bool = True) -> Path | None:
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(sorted(root.glob(pattern)))
    if exclude_registered:
        matches = [p for p in matches if "registered" not in p.name.lower()]
    existing = [p for p in matches if p.exists()]
    return existing[0] if existing else None


def parse_ims_overrides(entries: list[str]) -> dict[int, Path]:
    overrides: dict[int, Path] = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError(f"Invalid source override '{entry}'. Use channel=path.")
        channel_s, path_s = entry.split("=", 1)
        channel = int(channel_s.strip())
        path = Path(path_s.strip().strip('"')).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Source override for channel {channel} does not exist: {path}")
        if path.suffix.lower() not in {".ims", ".tif", ".tiff"}:
            raise ValueError(f"Source override for channel {channel} must be .ims, .tif, or .tiff: {path}")
        overrides[channel] = path
    return overrides


def source_ims_candidates(root: Path) -> list[Path]:
    candidates: list[Path] = []
    for pattern in ("**/*ch0*.ims", "**/*.ims"):
        candidates.extend(sorted(root.glob(pattern)))
    out: list[Path] = []
    seen = set()
    for path in candidates:
        key = str(path.resolve()).lower()
        if key in seen or not path.is_file():
            continue
        seen.add(key)
        out.append(path)
    return out


def resolve_source_ims_path(source_ims: str, cycle_label: str) -> Path:
    source_path = Path(source_ims).expanduser().resolve()
    if not source_path.exists():
        raise FileNotFoundError(f"{cycle_label} source .ims does not exist: {source_path}")
    if source_path.is_dir():
        resolved = one_match(source_ims_candidates(source_path), f"{cycle_label} source .ims under {source_path}")
        print(f"[batch] resolved {cycle_label} source folder to .ims: {resolved}", flush=True)
        return resolved
    if source_path.suffix.lower() != ".ims":
        raise ValueError(f"{cycle_label} source must be a .ims file or a folder containing one: {source_path}")
    return source_path


def source_cycle_root(source_ims: str, cycle_label: str) -> Path:
    source_path = resolve_source_ims_path(source_ims, cycle_label)
    target = cycle_label.lower()
    for parent in source_path.parents:
        if parent.name.lower() == target:
            print(f"[batch] using {cycle_label} source root from --{cycle_label.lower()}-src: {parent}", flush=True)
            return parent
    # Fallback for layouts like Cycle1/session/channel_folder/file.ims.
    fallback = source_path.parents[2] if len(source_path.parents) > 2 else source_path.parent
    print(f"[batch] using {cycle_label} fallback source root from source .ims: {fallback}", flush=True)
    return fallback


def source_session_root(source_ims: str, cycle_label: str) -> Path | None:
    if not source_ims:
        return None
    source_path = resolve_source_ims_path(source_ims, cycle_label)
    if len(source_path.parents) > 1:
        session_root = source_path.parents[1]
        print(f"[batch] using {cycle_label} session root from source .ims: {session_root}", flush=True)
        return session_root
    return source_path.parent


def infer_standard_paths(args) -> None:
    if not args.project_root:
        args.fixed_crop = require_value(args.fixed_crop, "--fixed-crop")
        if not args.cycle1_roi and not args.cycle1_roi_csv:
            args.cycle1_roi = require_value(args.cycle1_roi, "--cycle1-roi or --cycle1-roi-csv")
        args.cycle1_moving_crop = require_value(args.cycle1_moving_crop, "--cycle1-moving-crop")
        args.cycle1_transform_dir = require_value(args.cycle1_transform_dir, "--cycle1-transform-dir")
        args.cycle1_raw_output_dir = require_value(args.cycle1_raw_output_dir, "--cycle1-raw-output-dir")
        if not args.cycle2_roi and not args.cycle2_roi_csv:
            args.cycle2_roi = require_value(args.cycle2_roi, "--cycle2-roi or --cycle2-roi-csv")
        args.cycle2_moving_crop = require_value(args.cycle2_moving_crop, "--cycle2-moving-crop")
        args.cycle2_transform_dir = require_value(args.cycle2_transform_dir, "--cycle2-transform-dir")
        args.cycle2_raw_output_dir = require_value(args.cycle2_raw_output_dir, "--cycle2-raw-output-dir")
        return

    project_root = Path(args.project_root).expanduser().resolve()
    cycles_root = project_root / "cycles"
    cycle1_root = cycles_root / args.cycle1_id
    cycle2_root = cycles_root / args.cycle2_id
    if not args.fixed_crop:
        args.fixed_crop = str(first_dapi_crop(project_root / "reference" / "roi", "fixed reference ROI crop"))
    if not args.cycle1_roi_csv:
        csv_path = first_existing_glob(cycle1_root / "roi", ["XY_Coordinates.csv", "*.csv"], exclude_registered=False)
        if csv_path is not None:
            args.cycle1_roi_csv = str(csv_path)
    if not args.cycle1_roi and not args.cycle1_roi_csv:
        args.cycle1_roi = str(first_glob(cycle1_root / "roi", ["*.roi", "*.zip"], "Cycle1 ImageJ ROI", exclude_registered=False))
    if not args.cycle1_moving_crop:
        args.cycle1_moving_crop = str(first_dapi_crop(cycle1_root / "roi", "Cycle1 DAPI ROI crop"))
    if not args.cycle1_transform_dir:
        args.cycle1_transform_dir = str(cycle1_root / "registration" / args.run_id / "quicksyn_transforms")
    if not args.cycle1_raw_output_dir:
        args.cycle1_raw_output_dir = str(cycle1_root / "roi")
    if not args.cycle2_roi_csv:
        csv_path = first_existing_glob(cycle2_root / "roi", ["XY_Coordinates.csv", "*.csv"], exclude_registered=False)
        if csv_path is not None:
            args.cycle2_roi_csv = str(csv_path)
    if not args.cycle2_roi and not args.cycle2_roi_csv:
        args.cycle2_roi = str(first_glob(cycle2_root / "roi", ["*.roi", "*.zip"], "Cycle2 ImageJ ROI", exclude_registered=False))
    if not args.cycle2_moving_crop:
        args.cycle2_moving_crop = str(first_dapi_crop(cycle2_root / "roi", "Cycle2 DAPI ROI crop"))
    if not args.cycle2_transform_dir:
        args.cycle2_transform_dir = str(cycle2_root / "registration" / args.run_id / "quicksyn_transforms")
    if not args.cycle2_raw_output_dir:
        args.cycle2_raw_output_dir = str(cycle2_root / "roi")

    acquisition_root = Path(args.acquisition_root).expanduser().resolve() if args.acquisition_root else project_root.parent
    if args.cycle1_source_ims:
        args.cycle1_source_ims = str(resolve_source_ims_path(args.cycle1_source_ims, "Cycle1"))
    if args.cycle2_source_ims:
        args.cycle2_source_ims = str(resolve_source_ims_path(args.cycle2_source_ims, "Cycle2"))
    cycle1_source_root = source_cycle_root(args.cycle1_source_ims, "Cycle1") if args.cycle1_source_ims else acquisition_root / "Cycle1"
    cycle2_source_root = source_cycle_root(args.cycle2_source_ims, "Cycle2") if args.cycle2_source_ims else acquisition_root / "Cycle2"
    cycle1_overrides = parse_ims_overrides(args.cycle1_ims_override)
    cycle2_overrides = parse_ims_overrides(args.cycle2_ims_override)
    if args.cycle1_ims_map and cycle1_overrides:
        args.cycle1_ims_map = merge_ims_map(args.cycle1_ims_map, cycle1_overrides, "Cycle1")
    if args.cycle2_ims_map and cycle2_overrides:
        args.cycle2_ims_map = merge_ims_map(args.cycle2_ims_map, cycle2_overrides, "Cycle2")
    if args.cycle1_source_ims and not args.cycle1_ims_template and not args.cycle1_ims_map:
        args.cycle1_ims_map = merge_ims_map("", cycle1_overrides, "Cycle1") if cycle1_overrides else ""
    elif not args.cycle1_ims_template and not args.cycle1_ims_map:
        args.cycle1_ims_map = infer_channel_ims_map(
            cycle1_source_root,
            args.channels,
            "Cycle1",
            cycle1_overrides,
            source_ims=args.cycle1_source_ims,
            preferred_root=source_session_root(args.cycle1_source_ims, "Cycle1"),
        )
    if args.cycle2_source_ims and not args.cycle2_ims_template and not args.cycle2_ims_map:
        args.cycle2_ims_map = merge_ims_map("", cycle2_overrides, "Cycle2") if cycle2_overrides else ""
    elif not args.cycle2_ims_template and not args.cycle2_ims_map:
        args.cycle2_ims_map = infer_channel_ims_map(
            cycle2_source_root,
            args.channels,
            "Cycle2",
            cycle2_overrides,
            source_ims=args.cycle2_source_ims,
            preferred_root=source_session_root(args.cycle2_source_ims, "Cycle2"),
        )


def merge_ims_map(value: str, overrides: dict[int, Path], label: str) -> str:
    pairs: dict[int, str] = {}
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Invalid {label} source map entry: {item}")
        channel_s, path_s = item.split("=", 1)
        pairs[int(channel_s.strip())] = path_s.strip().strip('"')
    for channel, path in overrides.items():
        print(f"[batch] using {label} channel {channel} override: {path}", flush=True)
        pairs[channel] = str(path)
    return ";".join(f"{channel}={path}" for channel, path in sorted(pairs.items()))


def infer_channel_ims_map(
    cycle_source_dir: Path,
    channels: list[int],
    label: str,
    overrides: dict[int, Path] | None = None,
    source_ims: str = "",
    preferred_root: Path | None = None,
) -> str:
    if not cycle_source_dir.exists():
        raise FileNotFoundError(f"Could not infer {label} .ims files; folder not found: {cycle_source_dir}")
    overrides = overrides or {}
    search_root = preferred_root if preferred_root is not None and preferred_root.exists() else cycle_source_dir
    all_ims = sorted(search_root.rglob("*.ims"))
    if not all_ims and search_root != cycle_source_dir:
        all_ims = sorted(cycle_source_dir.rglob("*.ims"))
    source_path = Path(source_ims).expanduser().resolve() if source_ims else None
    fallback_by_channel = sequential_channel_candidates(all_ims, channels, source_path)
    entries = []
    for channel in channels:
        if channel in overrides:
            match = overrides[channel]
            print(f"[batch] using {label} channel {channel} override: {match}", flush=True)
        else:
            matches = [p for p in all_ims if ims_matches_channel(p, channel)]
            if matches:
                match = one_match(matches, f"{label} channel {channel} .ims")
                print(f"[batch] inferred {label} channel {channel} .ims: {match}", flush=True)
            elif channel in fallback_by_channel:
                match = fallback_by_channel[channel]
                print(
                    f"[batch] inferred {label} channel {channel} .ims by sibling order: {match}",
                    flush=True,
                )
            else:
                raise FileNotFoundError(
                    f"Could not infer {label} channel {channel} .ims under {search_root}. "
                    f"Pass --{label.lower()}-ims channel=path for this channel; the override may be .ims, .tif, or .tiff."
                )
        entries.append(f"{channel}={match}")
    return ";".join(entries)


def sequential_channel_candidates(
    all_ims: list[Path],
    channels: list[int],
    source_path: Path | None,
) -> dict[int, Path]:
    if source_path is not None:
        source_resolved = source_path.resolve()
        candidates = [p for p in all_ims if p.resolve() != source_resolved]
    else:
        candidates = all_ims[:]
    candidates = sorted(candidates, key=channel_sort_key)
    if len(candidates) < len(channels):
        return {}
    return {channel: candidates[idx] for idx, channel in enumerate(channels)}


def channel_sort_key(path: Path) -> tuple[int, str]:
    lowered = str(path).lower()
    match = re.search(r"ex[_-]?(\d+)", lowered)
    if match:
        return (int(match.group(1)), lowered)
    match = re.search(r"em[_-]?(\d+)", lowered)
    if match:
        return (int(match.group(1)), lowered)
    match = re.search(r"ch[_-]?(\d+)", lowered)
    if match:
        return (int(match.group(1)), lowered)
    return (10**9, lowered)


def ims_matches_channel(path: Path, channel: int) -> bool:
    lowered_stem = path.stem.lower()
    lowered_path = str(path).lower()
    patterns = (
        rf"(?<!\d)ch[_-]?0*{channel}(?!\d)",
        rf"(?<!\d)em[_-]?0*{channel}(?!\d)",
        rf"(?<!\d)channel[_ -]?0*{channel}(?!\d)",
    )
    return any(re.search(pattern, lowered_stem) or re.search(pattern, lowered_path) for pattern in patterns)


def read_uint16_zyx(path: Path) -> np.ndarray:
    img = sitk.ReadImage(str(path))
    arr = sitk.GetArrayFromImage(img)
    if arr.dtype != np.uint16:
        arr = np.rint(np.clip(arr, 0, 65535)).astype(np.uint16)
    return arr


def write_tiff_stack_from_image(input_path: Path, output_path: Path, progress_interval: float) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[tiff-stack] converting registered image to TIFF stack: {input_path} -> {output_path}", flush=True)
    with progress_heartbeat("tiff-stack read/write", progress_interval):
        arr = read_uint16_zyx(input_path)
        tifffile.imwrite(output_path, arr, imagej=True, metadata={"axes": "ZYX"})
    print(f"[tiff-stack] wrote TIFF stack: {output_path}", flush=True)
    return output_path


def registered_tiff_for_image(input_path: Path) -> Path:
    return input_path.with_suffix(".tif")


def channel_outputs_from_manifest(manifest_path: Path) -> dict[int, Path]:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    return {int(item["channel"]): Path(item["registered"]) for item in payload.get("outputs", [])}


def channel_tiff_outputs_from_manifest(manifest_path: Path) -> dict[int, Path]:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    out: dict[int, Path] = {}
    for item in payload.get("outputs", []):
        channel = int(item["channel"])
        path_s = item.get("registered_tif") or item.get("registered")
        path = Path(path_s)
        if path.suffix.lower() in {".tif", ".tiff"}:
            out[channel] = path
    return out


def write_multichannel_tiff(
    label: str,
    ch0_registered: Path,
    channel_manifest: Path,
    channels: list[int],
    output_path: Path,
    compression: str,
    tmp_dir: Path | None,
    progress_interval: float,
) -> None:
    print("", flush=True)
    print(f"===== {label}: write multichannel OME-TIFF =====", flush=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    channel_paths = {0: ch0_registered}
    channel_paths.update(channel_outputs_from_manifest(channel_manifest))
    ordered_channels = [0] + channels

    missing = [ch for ch in ordered_channels if ch not in channel_paths or not channel_paths[ch].exists()]
    if missing:
        raise FileNotFoundError(f"Missing registered channel output(s) for {label}: {missing}")

    first = read_uint16_zyx(channel_paths[0])
    z_size, y_size, x_size = first.shape
    temp_parent = tmp_dir or output_path.parent
    temp_parent.mkdir(parents=True, exist_ok=True)
    tmp_file = tempfile.NamedTemporaryFile(
        prefix=f"{output_path.stem}_",
        suffix=".npy",
        dir=str(temp_parent),
        delete=False,
    )
    tmp_file.close()
    tmp_path = Path(tmp_file.name)
    print(f"[tiff] staging uint16 stack as ZCYX memmap: {tmp_path}", flush=True)
    stack = np.lib.format.open_memmap(
        tmp_path,
        mode="w+",
        dtype=np.uint16,
        shape=(z_size, len(ordered_channels), y_size, x_size),
    )
    try:
        for c_index, channel in enumerate(ordered_channels):
            path = channel_paths[channel]
            print(f"[tiff] loading channel {channel}: {path}", flush=True)
            with progress_heartbeat(f"tiff load channel {channel}", progress_interval):
                arr = read_uint16_zyx(path)
            if arr.shape != (z_size, y_size, x_size):
                raise ValueError(
                    f"Channel {channel} shape mismatch: expected {(z_size, y_size, x_size)}, got {arr.shape}"
                )
            stack[:, c_index, :, :] = arr
            stack.flush()
            print(f"[tiff] staged channel {channel} ({c_index + 1}/{len(ordered_channels)})", flush=True)

        metadata = {
            "axes": "ZCYX",
            "Channel": {"Name": [f"ch{channel}" for channel in ordered_channels]},
        }
        kwargs = {
            "bigtiff": True,
            "ome": True,
            "metadata": metadata,
        }
        if compression:
            kwargs["compression"] = compression
        print(f"[tiff] writing OME-TIFF: {output_path}", flush=True)
        with progress_heartbeat("tiff write multichannel", progress_interval):
            tifffile.imwrite(output_path, stack, **kwargs)
        print(f"[tiff] wrote OME-TIFF: {output_path}", flush=True)
    finally:
        del stack
        try:
            tmp_path.unlink()
            print(f"[tiff] removed temp stack: {tmp_path}", flush=True)
        except FileNotFoundError:
            pass


def write_registered_tiff_manifest(
    project_root: Path,
    manifest_path: Path,
    cycle_entries: list[dict],
) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "project_root": str(project_root),
        "cycles": cycle_entries,
    }
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[imaris] wrote registered TIFF manifest: {manifest_path}", flush=True)


def collect_cycle_tiff_entry(
    cycle_name: str,
    ch0_tiff: Path,
    channel_manifest: Path,
    channels: list[int],
) -> dict:
    channel_paths = {0: ch0_tiff}
    channel_paths.update(channel_tiff_outputs_from_manifest(channel_manifest))
    ordered_channels = [0] + channels
    missing = [ch for ch in ordered_channels if ch not in channel_paths or not channel_paths[ch].exists()]
    if missing:
        raise FileNotFoundError(f"Missing registered TIFF stack(s) for {cycle_name}: {missing}")
    return {
        "cycle": cycle_name,
        "channels": [
            {"channel": int(channel), "path": str(channel_paths[channel])}
            for channel in ordered_channels
        ],
    }


def run_imaris_wrapper_if_requested(
    command_template: str,
    manifest_path: Path,
    output_path: Path,
    project_root: Path,
    cycle1_output_dir: Path,
    cycle2_output_dir: Path,
) -> None:
    if not command_template:
        print("[imaris] no --imaris-wrapper-cmd supplied; skipping .ims conversion", flush=True)
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    command = command_template.format(
        manifest=str(manifest_path),
        output=str(output_path),
        project=str(project_root),
        cycle1_output=str(cycle1_output_dir),
        cycle2_output=str(cycle2_output_dir),
    )
    run_shell_step("run Imaris/PetaKit conversion wrapper", command)


def parse_float_triplet(value: str, default_value: str) -> tuple[float, float, float]:
    raw = value.strip() if value else default_value.strip()
    parts = [p.strip() for p in raw.split(",") if p.strip()]
    if len(parts) == 1:
        v = float(parts[0])
        return (v, v, v)
    if len(parts) == 3:
        return (float(parts[0]), float(parts[1]), float(parts[2]))
    raise ValueError(f"Expected one scalar or x,y,z pixel sizes, got: {raw}")


def matlab_quote(value: str) -> str:
    return value.replace("'", "''")


def matlab_path(value: Path) -> str:
    return matlab_quote(str(value))


def matlab_pattern_path(path: Path, root: Path) -> str:
    try:
        rel = path.resolve().relative_to(root.resolve())
        value = str(rel)
    except ValueError:
        value = str(path.resolve())
    return matlab_quote(value)


def flatten_tiff_channel_patterns(cycle_entries: list[dict], root: Path) -> list[str]:
    patterns: list[str] = []
    for cycle in cycle_entries:
        for channel in cycle["channels"]:
            patterns.append(matlab_pattern_path(Path(channel["path"]), root))
    return patterns


def write_xr_imaris_script(
    project_root: Path,
    output_path: Path,
    pixel_sizes: tuple[float, float, float],
    cpus_per_task: int,
    channel_patterns: list[str],
    script_path: Path,
) -> Path:
    script_path.parent.mkdir(parents=True, exist_ok=True)
    pattern_lines = "\n".join(f"    '{pattern}';" for pattern in channel_patterns)
    pixel_text = " ".join(f"{v:g}" for v in pixel_sizes)
    body = f"""% Auto-generated by run_cycif_cycle_batch.py
inputPath = '{matlab_path(project_root)}';
outputPath = '{matlab_path(output_path)}';
channelPatterns = {{
{pattern_lines}
}};
XR_imaris_conversion_data_wrapper(inputPath, ...
    'outputPath', outputPath, ...
    'pixelSizes', [{pixel_text}], ...
    'cpusPerTask', {int(cpus_per_task)}, ...
    'channelPatterns', channelPatterns);
"""
    script_path.write_text(body, encoding="utf-8")
    print(f"[xr-imaris] wrote MATLAB wrapper script: {script_path}", flush=True)
    return script_path


def run_xr_imaris_conversion(
    enabled: bool,
    matlab_exe: str,
    matlab_start: bool,
    project_root: Path,
    output_path: Path,
    pixel_sizes: tuple[float, float, float],
    cpus_per_task: int,
    cycle_entries: list[dict],
) -> None:
    if not enabled:
        return
    output_path.mkdir(parents=True, exist_ok=True)
    script_path = project_root / "run_xr_imaris_conversion.m"
    patterns = flatten_tiff_channel_patterns(cycle_entries, project_root)
    write_xr_imaris_script(project_root, output_path, pixel_sizes, cpus_per_task, patterns, script_path)
    script_for_matlab = matlab_path(script_path)
    if matlab_start and os.name == "nt":
        cmd = f'start "" "{matlab_exe}" -r "run(\'{script_for_matlab}\')"'
        run_shell_step("launch XR Imaris conversion in MATLAB", cmd)
    else:
        cmd = f'"{matlab_exe}" -batch "run(\'{script_for_matlab}\')"'
        run_shell_step("run XR Imaris conversion in MATLAB", cmd)


def run_apply_channels(
    label: str,
    source_ims: str,
    ims_template: str,
    ims_map: str,
    roi: str,
    roi_csv: str,
    channels: list[int],
    fixed_crop: str,
    moving_crop: str,
    intra_cycle_align_channels: list[int],
    pre_crop_intra_cycle_align_channels: list[int],
    intra_cycle_transform: str,
    intra_cycle_downsample_factor: float,
    intra_padding_xy: int,
    transform_dir: Path,
    raw_output_dir: str,
    registered_output_dir: Path,
    spacing: str,
    threads: int,
    progress_interval: float,
    output_pixeltype: str,
    registered_format: str,
    uint16_scaling: str,
    open_qc: bool,
    save_pre_crop_aligned_roi: bool,
    keep_pre_crop_temp: bool,
) -> float:
    command = [
        sys.executable,
        script_path("apply_ims_roi_channels.py"),
        "--channels",
        ",".join(str(ch) for ch in channels),
        "--fixed-crop",
        fixed_crop,
        "--moving-crop-reference",
        moving_crop,
        "--transform-dir",
        str(transform_dir),
        "--raw-output-dir",
        raw_output_dir,
        "--registered-output-dir",
        str(registered_output_dir),
        "--spacing",
        spacing,
        "--threads",
        str(threads),
        "--progress-interval",
        str(progress_interval),
        "--output-pixeltype",
        output_pixeltype,
        "--registered-format",
        registered_format,
        "--uint16-scaling",
        uint16_scaling,
    ]
    if intra_cycle_align_channels:
        command.extend(["--intra-cycle-align-channels", ",".join(str(ch) for ch in intra_cycle_align_channels)])
    if pre_crop_intra_cycle_align_channels:
        command.extend(
            [
                "--pre-crop-intra-cycle-align-channels",
                ",".join(str(ch) for ch in pre_crop_intra_cycle_align_channels),
            ]
        )
        command.extend(["--intra-padding-xy", str(intra_padding_xy)])
        if save_pre_crop_aligned_roi:
            command.append("--save-pre-crop-aligned-roi")
        if keep_pre_crop_temp:
            command.append("--keep-pre-crop-temp")
    if intra_cycle_align_channels or pre_crop_intra_cycle_align_channels:
        command.extend(["--intra-cycle-transform", intra_cycle_transform])
        command.extend(["--intra-cycle-downsample-factor", str(intra_cycle_downsample_factor)])
        if open_qc:
            command.append("--open-qc")
    if roi_csv:
        command.extend(["--roi-csv", roi_csv])
    elif roi:
        command.extend(["--roi-imagej", roi])
    else:
        raise ValueError(f"{label} needs a CSV ROI or ImageJ ROI")
    if source_ims:
        command.extend(["--input-ims", source_ims])
    if ims_map:
        command.extend(["--input-ims-map", ims_map])
    elif ims_template:
        command.extend(["--input-ims-template", ims_template])
    elif not source_ims:
        raise ValueError(f"{label} needs an .ims source, template, or map")

    return run_step(
        f"{label}: extract ROI and apply transforms to channels",
        command,
    )


def run_dapi_apply(
    label: str,
    fixed_crop: str,
    moving_crop: str,
    transform_dir: Path,
    threads: int,
    progress_interval: float,
    output_pixeltype: str,
    uint16_scaling: str,
    open_qc: bool,
) -> float:
    expected_output = default_registered_output(Path(moving_crop).expanduser().resolve(), transform_dir)
    if expected_output.exists() and expected_output.stat().st_size > 0:
        print(
            f"[batch] skipping {label} ch0 apply; registered output already exists: {expected_output}",
            flush=True,
        )
        return 0.0
    cmd = [
        sys.executable,
        script_path("ants_roi_quicksyn.py"),
        "--mode",
        "apply",
        "--fixed-crop",
        fixed_crop,
        "--moving-crop",
        moving_crop,
        "--transform-dir",
        str(transform_dir),
        "--threads",
        str(threads),
        "--progress-interval",
        str(progress_interval),
        "--output-pixeltype",
        output_pixeltype,
        "--uint16-scaling",
        uint16_scaling,
    ]
    if open_qc:
        cmd.append("--open-qc")
    return run_step(f"{label}: apply existing DAPI transform", cmd)


def run_dapi_registration(
    label: str,
    fixed_crop: str,
    moving_crop: str,
    transform_dir: Path,
    downsample_factor: float,
    threads: int,
    progress_interval: float,
    output_pixeltype: str,
    uint16_scaling: str,
    open_qc: bool,
) -> float:
    cmd = [
        sys.executable,
        script_path("ants_roi_quicksyn.py"),
        "--fixed-crop",
        fixed_crop,
        "--moving-crop",
        moving_crop,
        "--transform-dir",
        str(transform_dir),
        "--downsample-factor",
        str(downsample_factor),
        "--threads",
        str(threads),
        "--progress-interval",
        str(progress_interval),
        "--output-pixeltype",
        output_pixeltype,
        "--uint16-scaling",
        uint16_scaling,
        "--verbose",
    ]
    if open_qc:
        cmd.append("--open-qc")
    return run_step(f"{label}: register DAPI ROI", cmd)


def main() -> int:
    args = build_parser().parse_args()
    infer_standard_paths(args)
    start = time.time()
    cycle1_transform_dir = Path(args.cycle1_transform_dir).expanduser().resolve()
    cycle2_transform_dir = Path(args.cycle2_transform_dir).expanduser().resolve()
    cycle1_output_dir = output_dir_or_default(args.cycle1_registered_output_dir, cycle1_transform_dir)
    cycle2_output_dir = output_dir_or_default(args.cycle2_registered_output_dir, cycle2_transform_dir)
    cycle1_tiff = tiff_output_or_default(args.cycle1_tiff_output, cycle1_transform_dir, "cycle_001")
    cycle2_tiff = tiff_output_or_default(args.cycle2_tiff_output, cycle2_transform_dir, "cycle_002")
    project_root = Path(args.project_root).expanduser().resolve() if args.project_root else Path.cwd().resolve()
    registered_tiff_manifest = (
        Path(args.registered_tiff_manifest).expanduser().resolve()
        if args.registered_tiff_manifest
        else (project_root / "registered_tiff_manifest.json").resolve()
    )
    imaris_output = imaris_output_or_default(args.imaris_output, project_root).resolve()
    xr_output_path = (
        Path(args.xr_output_path).expanduser().resolve()
        if args.xr_output_path
        else (project_root / "xr_imaris").resolve()
    )
    xr_pixel_sizes = parse_float_triplet(args.xr_pixel_sizes, args.spacing)
    tmp_dir = Path(args.tmp_dir).expanduser().resolve() if args.tmp_dir else None

    print("[batch] starting CyCIF cycle batch", flush=True)
    print(f"[batch] channels={args.channels} spacing={args.spacing} threads={args.threads}", flush=True)
    print(f"[batch] registered format={args.registered_format} skip_ome={args.skip_ome}", flush=True)
    if args.xr_imaris:
        print(
            f"[batch] XR Imaris conversion enabled outputPath={xr_output_path} "
            f"pixelSizes={xr_pixel_sizes}",
            flush=True,
        )
    print(f"[batch] fixed crop={args.fixed_crop}", flush=True)
    print(f"[batch] cycle1 roi={args.cycle1_roi}", flush=True)
    print(f"[batch] cycle1 roi csv={args.cycle1_roi_csv}", flush=True)
    print(f"[batch] cycle1 moving crop={args.cycle1_moving_crop}", flush=True)
    print(f"[batch] cycle1 transform dir={cycle1_transform_dir}", flush=True)
    print(f"[batch] cycle2 roi={args.cycle2_roi}", flush=True)
    print(f"[batch] cycle2 roi csv={args.cycle2_roi_csv}", flush=True)
    print(f"[batch] cycle2 moving crop={args.cycle2_moving_crop}", flush=True)
    print(f"[batch] cycle2 transform dir={cycle2_transform_dir}", flush=True)
    print(f"[batch] cycle1 outputs={cycle1_output_dir}", flush=True)
    print(f"[batch] cycle2 outputs={cycle2_output_dir}", flush=True)
    cycle1_intra_channels = (
        args.cycle1_intra_cycle_align_channels
        if args.cycle1_intra_cycle_align_channels is not None
        else args.intra_cycle_align_channels
    )
    cycle2_intra_channels = (
        args.cycle2_intra_cycle_align_channels
        if args.cycle2_intra_cycle_align_channels is not None
        else args.intra_cycle_align_channels
    )
    cycle1_pre_crop_intra_channels = (
        args.cycle1_pre_crop_intra_cycle_align_channels
        if args.cycle1_pre_crop_intra_cycle_align_channels is not None
        else args.pre_crop_intra_cycle_align_channels
    )
    cycle2_pre_crop_intra_channels = (
        args.cycle2_pre_crop_intra_cycle_align_channels
        if args.cycle2_pre_crop_intra_cycle_align_channels is not None
        else args.pre_crop_intra_cycle_align_channels
    )
    print(f"[batch] cycle1 intra-cycle alignment channels={cycle1_intra_channels}", flush=True)
    print(f"[batch] cycle2 intra-cycle alignment channels={cycle2_intra_channels}", flush=True)
    print(f"[batch] cycle1 pre-crop intra-cycle alignment channels={cycle1_pre_crop_intra_channels}", flush=True)
    print(f"[batch] cycle2 pre-crop intra-cycle alignment channels={cycle2_pre_crop_intra_channels}", flush=True)
    timing_log: list[dict] = []

    if not args.skip_cycle1_ch0_apply:
        elapsed = run_dapi_apply(
            "cycle1",
            args.fixed_crop,
            args.cycle1_moving_crop,
            cycle1_transform_dir,
            args.threads,
            args.progress_interval,
            args.output_pixeltype,
            args.uint16_scaling,
            args.open_qc,
        )
        append_timing(timing_log, "cycle1 apply existing DAPI transform", elapsed)
    else:
        print("[batch] skipping Cycle1 ch0 apply", flush=True)
        append_timing(timing_log, "cycle1 apply existing DAPI transform skipped", 0.0)

    elapsed = run_apply_channels(
        "cycle1",
        args.cycle1_source_ims,
        args.cycle1_ims_template,
        args.cycle1_ims_map,
        args.cycle1_roi,
        args.cycle1_roi_csv,
        args.channels,
        args.fixed_crop,
        args.cycle1_moving_crop,
        cycle1_intra_channels,
        cycle1_pre_crop_intra_channels,
        args.intra_cycle_transform,
        args.intra_cycle_downsample_factor,
        args.intra_padding_xy,
        cycle1_transform_dir,
        args.cycle1_raw_output_dir,
        cycle1_output_dir,
        args.spacing,
        args.threads,
        args.progress_interval,
        args.output_pixeltype,
        args.registered_format,
        args.uint16_scaling,
        args.open_qc,
        args.save_pre_crop_aligned_roi,
        args.keep_pre_crop_temp,
    )
    append_timing(timing_log, "cycle1 extract/apply channels", elapsed)
    cycle1_ch0_nrrd = default_registered_output(Path(args.cycle1_moving_crop).expanduser().resolve(), cycle1_transform_dir)
    cycle1_ch0_tiff = registered_tiff_for_image(cycle1_ch0_nrrd)
    if args.registered_format in {"tif", "both"}:
        stage_start = time.time()
        write_tiff_stack_from_image(cycle1_ch0_nrrd, cycle1_ch0_tiff, args.progress_interval)
        record_timing(timing_log, "cycle1 write registered ch0 TIFF", stage_start)
    if not args.skip_ome:
        stage_start = time.time()
        write_multichannel_tiff(
            "cycle1",
            cycle1_ch0_nrrd,
            cycle1_output_dir / "registered_channel_outputs.json",
            args.channels,
            cycle1_tiff,
            args.tiff_compression,
            tmp_dir,
            args.progress_interval,
        )
        record_timing(timing_log, "cycle1 write multichannel OME-TIFF", stage_start)

    if not args.skip_cycle2_dapi_registration:
        elapsed = run_dapi_registration(
            "cycle2",
            args.fixed_crop,
            args.cycle2_moving_crop,
            cycle2_transform_dir,
            args.downsample_factor,
            args.threads,
            args.progress_interval,
            args.output_pixeltype,
            args.uint16_scaling,
            args.open_qc,
        )
        append_timing(timing_log, "cycle2 register DAPI ROI", elapsed)
    else:
        print("[batch] skipping Cycle2 DAPI registration", flush=True)
        append_timing(timing_log, "cycle2 register DAPI ROI skipped", 0.0)

    elapsed = run_apply_channels(
        "cycle2",
        args.cycle2_source_ims,
        args.cycle2_ims_template,
        args.cycle2_ims_map,
        args.cycle2_roi,
        args.cycle2_roi_csv,
        args.channels,
        args.fixed_crop,
        args.cycle2_moving_crop,
        cycle2_intra_channels,
        cycle2_pre_crop_intra_channels,
        args.intra_cycle_transform,
        args.intra_cycle_downsample_factor,
        args.intra_padding_xy,
        cycle2_transform_dir,
        args.cycle2_raw_output_dir,
        cycle2_output_dir,
        args.spacing,
        args.threads,
        args.progress_interval,
        args.output_pixeltype,
        args.registered_format,
        args.uint16_scaling,
        args.open_qc,
        args.save_pre_crop_aligned_roi,
        args.keep_pre_crop_temp,
    )
    append_timing(timing_log, "cycle2 extract/apply channels", elapsed)
    cycle2_ch0_nrrd = default_registered_output(Path(args.cycle2_moving_crop).expanduser().resolve(), cycle2_transform_dir)
    cycle2_ch0_tiff = registered_tiff_for_image(cycle2_ch0_nrrd)
    if args.registered_format in {"tif", "both"}:
        stage_start = time.time()
        write_tiff_stack_from_image(cycle2_ch0_nrrd, cycle2_ch0_tiff, args.progress_interval)
        record_timing(timing_log, "cycle2 write registered ch0 TIFF", stage_start)
    if not args.skip_ome:
        stage_start = time.time()
        write_multichannel_tiff(
            "cycle2",
            cycle2_ch0_nrrd,
            cycle2_output_dir / "registered_channel_outputs.json",
            args.channels,
            cycle2_tiff,
            args.tiff_compression,
            tmp_dir,
            args.progress_interval,
        )
        record_timing(timing_log, "cycle2 write multichannel OME-TIFF", stage_start)

    if args.registered_format in {"tif", "both"}:
        stage_start = time.time()
        cycle_entries = [
            collect_cycle_tiff_entry(
                "cycle_001",
                cycle1_ch0_tiff,
                cycle1_output_dir / "registered_channel_outputs.json",
                args.channels,
            ),
            collect_cycle_tiff_entry(
                "cycle_002",
                cycle2_ch0_tiff,
                cycle2_output_dir / "registered_channel_outputs.json",
                args.channels,
            ),
        ]
        write_registered_tiff_manifest(project_root, registered_tiff_manifest, cycle_entries)
        record_timing(timing_log, "write registered TIFF manifest", stage_start)
        stage_start = time.time()
        run_imaris_wrapper_if_requested(
            args.imaris_wrapper_cmd,
            registered_tiff_manifest,
            imaris_output,
            project_root,
            cycle1_output_dir,
            cycle2_output_dir,
        )
        record_timing(timing_log, "run optional Imaris wrapper", stage_start)
        stage_start = time.time()
        run_xr_imaris_conversion(
            args.xr_imaris,
            args.matlab,
            args.matlab_start,
            project_root,
            xr_output_path,
            xr_pixel_sizes,
            args.threads,
            cycle_entries,
        )
        record_timing(timing_log, "run optional XR Imaris conversion", stage_start)
    elif args.xr_imaris:
        raise ValueError("--xr-imaris requires --reg-format tif or --reg-format both")

    total_elapsed = record_timing(timing_log, "batch total", start)
    timing_path = project_root / "batch_stage_timing_log.json"
    timing_path.write_text(json.dumps({"timings": timing_log, "total_seconds": total_elapsed}, indent=2), encoding="utf-8")
    print(f"[timing] wrote batch timing log: {timing_path}", flush=True)
    print(f"[batch] finished all cycles in {format_duration(total_elapsed)}", flush=True)
    return 0


if __name__ == "__main__":
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    raise SystemExit(main())
