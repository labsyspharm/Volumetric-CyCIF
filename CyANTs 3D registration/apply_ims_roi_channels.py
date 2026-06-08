#!/usr/bin/env python3
# Developed by Alex Wong
# Cite: Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in
# Preprint: https://doi.org/10.64898/2026.05.17.725158
# Registration method: ANTsX/ANTsPy - https://github.com/ANTsX/ANTsPy

"""Extract ROI channels from .ims or TIFF sources and apply an existing ANTs ROI transform."""

from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile

from ants_roi_quicksyn import (
    ants_to_sitk_uint16,
    configure_ants_runtime,
    describe_image,
    downsample_by_spacing,
    export_center_slice_qc,
    format_duration,
    load_transform_manifest,
    progress_heartbeat,
    registration_run_dir,
    write_ants_image,
)
from extract_ims_roi_channels import parse_channels, parse_imagej_roi, parse_z_range
from cyants_io import (
    Tile,
    build_ims_volume_spec,
    parse_spacing_override,
    read_ims_tile_from_spec,
    stem_for_output,
)

ants = None

TIFF_SUFFIXES = {".tif", ".tiff"}


@dataclass(frozen=True)
class ChannelSourceSpec:
    path: Path
    kind: str
    dataset_key: str
    axis_order: str
    shape_zyx: tuple[int, int, int]
    spacing_xyz: tuple[float, float, float]
    ims_spec: object | None = None


def parse_optional_channels(value: str) -> list[int]:
    value = value.strip()
    if not value:
        return []
    return parse_channels(value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Use a saved ImageJ ROI to extract full-resolution channel crops from .ims or TIFF sources, "
            "then apply an existing ROI registration transform to each crop."
        )
    )
    parser.add_argument(
        "--input-ims",
        "--input-source",
        default="",
        help="Source .ims or .tif/.tiff stack. Use this when all channels share one source .ims.",
    )
    parser.add_argument(
        "--input-ims-template",
        "--input-source-template",
        default="",
        help=(
            "Optional per-channel source path template with {channel}, for example "
            "G:\\...\\img_0002_ch{channel}.ims or G:\\...\\ch{channel}.tif."
        ),
    )
    parser.add_argument(
        "--input-ims-map",
        "--input-source-map",
        "--input-map",
        default="",
        help=(
            "Optional explicit per-channel source map as channel=path entries separated by semicolons. "
            "Each path can be .ims, .tif, or .tiff, for example "
            "1=G:\\...\\ch1_corrected.tif;2=G:\\...\\img_0002_ch0.ims."
        ),
    )
    parser.add_argument("--roi-imagej", default="", help="ImageJ/Fiji .roi or ROI .zip file")
    parser.add_argument(
        "--roi-csv",
        default="",
        help=(
            "CSV containing full-image XY coordinates. Expected columns: X,Y. "
            "The min/max X/Y bounds are used; Z defaults to full source depth unless --z-range is supplied."
        ),
    )
    parser.add_argument(
        "--z-range",
        type=parse_z_range,
        default=None,
        help="Optional ImageJ 1-based z slice range as first_slice,last_slice. Default: full source z depth.",
    )
    parser.add_argument("--channels", required=True, type=parse_channels, help="Channels, for example 1-4")
    parser.add_argument("--fixed-crop", required=True, help="Fixed/reference full-resolution ROI crop NRRD")
    parser.add_argument(
        "--moving-crop-reference",
        default="",
        help=(
            "Moving DAPI ROI crop NRRD used to create the registration transform. "
            "Extracted channel ROIs inherit this image geometry before transforms are applied."
        ),
    )
    parser.add_argument("--transform-dir", required=True, help="Directory containing fwdtransforms.json")
    parser.add_argument("--raw-output-dir", required=True, help="Directory for extracted raw channel ROI NRRDs")
    parser.add_argument(
        "--no-reuse-raw-roi",
        action="store_true",
        help="Always re-extract raw ROI NRRDs from source data even when matching raw ROI files already exist.",
    )
    parser.add_argument(
        "--intra-cycle-align-channels",
        default="",
        type=parse_optional_channels,
        help=(
            "Optional channel list to align to the moving DAPI crop before applying the cycle transform, "
            "for example 1,3. Intended for per-channel lateral/chromatic offsets."
        ),
    )
    parser.add_argument(
        "--pre-crop-intra-cycle-align-channels",
        "--precrop-ch",
        default="",
        type=parse_optional_channels,
        help=(
            "Optional channel list to align to same-cycle DAPI using a padded ROI before "
            "recropping to the requested ROI. Use when offset channels would be clipped by "
            "cropping first, for example 2,3."
        ),
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
        help="Also save the recropped pre-crop intra-cycle aligned ROI intermediate. Default: keep it in memory only.",
    )
    parser.add_argument(
        "--keep-pre-crop-temp",
        "--keep-precrop-temp",
        action="store_true",
        help="Keep temporary padded DAPI/channel NRRDs used for pre-crop intra-cycle alignment.",
    )
    parser.add_argument(
        "--intra-cycle-transform",
        default="TRSAA",
        help=(
            "ANTsPy transform for intra-cycle channel-to-DAPI alignment. "
            "Default: TRSAA (translation, rigid, similarity, affine)."
        ),
    )
    parser.add_argument(
        "--intra-cycle-downsample-factor",
        type=float,
        default=8.0,
        help="Downsample factor for intra-cycle channel-to-DAPI alignment. Default: 8.",
    )
    parser.add_argument(
        "--qc-output-dir",
        default="",
        help="Directory for intra-cycle before/after QC PNGs. Default: <registration-run>/qc/intra_cycle_channel_alignment",
    )
    parser.add_argument(
        "--no-qc",
        dest="qc",
        action="store_false",
        help="Disable intra-cycle center-slice QC PNG output.",
    )
    parser.set_defaults(qc=True)
    parser.add_argument("--open-qc", action="store_true", help="Open intra-cycle QC PNGs after writing them")
    parser.add_argument(
        "--qc-max-panel-side",
        type=int,
        default=900,
        help="Maximum displayed side length per QC panel. Use 0 for original downsampled size. Default: 900.",
    )
    parser.add_argument(
        "--registered-output-dir",
        default="",
        help="Directory for registered channel NRRDs. Default: <registration-run>/outputs",
    )
    parser.add_argument("--prefix", default="", help="Optional output filename prefix")
    parser.add_argument("--ims-resolution-level", type=int, default=0, help=".ims resolution level. Default: 0")
    parser.add_argument("--ims-timepoint", type=int, default=0, help=".ims timepoint. Default: 0")
    parser.add_argument(
        "--ims-dataset-template",
        default="",
        help=(
            "Optional HDF5 dataset template with {level}, {timepoint}, and {channel}; "
            "for example DataSet/ResolutionLevel {level}/TimePoint {timepoint}/Channel {channel}/Data"
        ),
    )
    parser.add_argument(
        "--ims-axis-order",
        choices=["auto", "zyx", "xyz"],
        default="auto",
        help="Axis order for .ims dataset interpretation. Default: auto",
    )
    parser.add_argument(
        "--spacing",
        type=parse_spacing_override,
        default=None,
        help="Override spacing in x,y,z or isotropic scalar, for example 0.711",
    )
    parser.add_argument(
        "--interpolator",
        default="linear",
        help="Interpolator for applying transforms. Use nearestNeighbor for labels. Default: linear",
    )
    parser.add_argument(
        "--output-pixeltype",
        choices=["uint16", "float32"],
        default="uint16",
        help="Pixel type for registered channel outputs. Default: uint16",
    )
    parser.add_argument(
        "--registered-format",
        "--reg-format",
        choices=["nrrd", "tif", "both"],
        default="nrrd",
        help="Registered channel output format. Default: nrrd.",
    )
    parser.add_argument(
        "--uint16-scaling",
        choices=["minmax", "robust", "clip"],
        default="minmax",
        help="How to convert transformed float images to uint16. Default: minmax",
    )
    parser.add_argument("--threads", type=int, default=0, help="Number of CPU threads for ANTs/ITK. Default: 0")
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=30.0,
        help="Print elapsed-time heartbeat every N seconds during long operations. Use 0 to disable. Default: 30.",
    )
    return parser


def parse_input_ims_map(value: str) -> dict[int, Path]:
    mapping: dict[int, Path] = {}
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Invalid --input-map entry: {item}")
        channel_s, path_s = item.split("=", 1)
        mapping[int(channel_s.strip())] = Path(path_s.strip().strip('"')).expanduser().resolve()
    return mapping


def channel_ims_path(args, channel: int, input_map: dict[int, Path] | None = None) -> Path:
    if input_map and channel in input_map:
        return input_map[channel]
    if args.input_ims_template:
        return Path(args.input_ims_template.format(channel=channel)).expanduser().resolve()
    if not args.input_ims:
        raise ValueError("Provide --input-source, --input-source-template, or --input-map")
    return Path(args.input_ims).expanduser().resolve()


def tiff_stack_shape_zyx(path: Path) -> tuple[int, int, int]:
    with tifffile.TiffFile(str(path)) as tif:
        if not tif.series:
            raise ValueError(f"TIFF has no image series: {path}")
        series = tif.series[0]
        axes = str(getattr(series, "axes", "")).upper()
        shape = tuple(int(v) for v in series.shape)
        if {"Z", "Y", "X"}.issubset(set(axes)):
            return (shape[axes.index("Z")], shape[axes.index("Y")], shape[axes.index("X")])
        if len(tif.pages) > 1 and len(tif.pages[0].shape) == 2:
            y_size, x_size = (int(v) for v in tif.pages[0].shape)
            return (len(tif.pages), y_size, x_size)
        if len(shape) == 3:
            return (shape[0], shape[1], shape[2])
        if len(shape) == 2:
            return (1, shape[0], shape[1])
    raise ValueError(f"Could not interpret TIFF as a 3D ZYX stack: {path}")


def build_tiff_source_spec(args, input_path: Path) -> ChannelSourceSpec:
    shape_zyx = tiff_stack_shape_zyx(input_path)
    spacing_xyz = args.spacing if args.spacing is not None else (1.0, 1.0, 1.0)
    with tifffile.TiffFile(str(input_path)) as tif:
        series = tif.series[0]
        axes = str(getattr(series, "axes", ""))
        dataset_key = f"TIFF series 0 axes={axes or 'unknown'} pages={len(tif.pages)}"
    return ChannelSourceSpec(
        path=input_path,
        kind="tiff",
        dataset_key=dataset_key,
        axis_order="zyx",
        shape_zyx=shape_zyx,
        spacing_xyz=spacing_xyz,
    )


def build_channel_spec(args, input_path: Path, channel: int) -> ChannelSourceSpec:
    if not input_path.exists():
        raise FileNotFoundError(f"Channel {channel} source does not exist: {input_path}")
    suffix = input_path.suffix.lower()
    if suffix in TIFF_SUFFIXES:
        return build_tiff_source_spec(args, input_path)
    if suffix != ".ims":
        raise ValueError(
            f"Unsupported channel source extension for channel {channel}: {input_path}. "
            "Use .ims, .tif, or .tiff."
        )
    dataset_path = (
        args.ims_dataset_template.format(
            level=args.ims_resolution_level,
            timepoint=args.ims_timepoint,
            channel=channel,
        )
        if args.ims_dataset_template
        else ""
    )
    spec = build_ims_volume_spec(
        path=input_path,
        res_level=args.ims_resolution_level,
        timepoint=args.ims_timepoint,
        channel=channel,
        dataset_path=dataset_path,
        axis_order=args.ims_axis_order,
    )
    if args.spacing is not None:
        spec = type(spec)(
            path=spec.path,
            dataset_key=spec.dataset_key,
            axis_order=spec.axis_order,
            shape_zyx=spec.shape_zyx,
            spacing_xyz=args.spacing,
        )
    return ChannelSourceSpec(
        path=input_path,
        kind="ims",
        dataset_key=spec.dataset_key,
        axis_order=spec.axis_order,
        shape_zyx=spec.shape_zyx,
        spacing_xyz=spec.spacing_xyz,
        ims_spec=spec,
    )


def tiff_series_array_to_zyx(arr: np.ndarray, axes: str) -> np.ndarray:
    axes = axes.upper()
    if not axes:
        if arr.ndim == 2:
            return arr[np.newaxis, :, :]
        if arr.ndim == 3:
            return arr
        raise ValueError(f"Cannot infer TIFF axes for shape {arr.shape}")

    arr_work = arr
    axes_work = axes
    for axis in list(axes_work):
        if axis in {"Z", "Y", "X"}:
            continue
        axis_index = axes_work.index(axis)
        if arr_work.shape[axis_index] != 1:
            raise ValueError(f"Cannot read multi-dimensional TIFF with non-singleton axis {axis}: shape={arr.shape} axes={axes}")
        arr_work = np.take(arr_work, 0, axis=axis_index)
        axes_work = axes_work[:axis_index] + axes_work[axis_index + 1 :]

    if "Z" not in axes_work and {"Y", "X"}.issubset(set(axes_work)):
        y_axis = axes_work.index("Y")
        x_axis = axes_work.index("X")
        arr_yx = np.moveaxis(arr_work, (y_axis, x_axis), (0, 1))
        return arr_yx[np.newaxis, :, :]
    if not {"Z", "Y", "X"}.issubset(set(axes_work)):
        raise ValueError(f"Cannot read TIFF as ZYX; shape={arr.shape} axes={axes}")
    return np.moveaxis(
        arr_work,
        (axes_work.index("Z"), axes_work.index("Y"), axes_work.index("X")),
        (0, 1, 2),
    )


def validate_tile_in_source(spec: ChannelSourceSpec, tile: Tile) -> None:
    z_size, y_size, x_size = spec.shape_zyx
    if tile.x0 < 0 or tile.y0 < 0 or tile.z0 < 0 or tile.x1 > x_size or tile.y1 > y_size or tile.z1 > z_size:
        raise ValueError(
            f"ROI tile {tile} extends outside {spec.kind} source bounds "
            f"xyz-size=({x_size}, {y_size}, {z_size}): {spec.path}"
        )


def read_tiff_tile_from_spec(spec: ChannelSourceSpec, tile: Tile) -> sitk.Image:
    validate_tile_in_source(spec, tile)
    with tifffile.TiffFile(str(spec.path)) as tif:
        if len(tif.pages) > 1 and len(tif.pages[0].shape) == 2:
            first_plane = np.asarray(tif.pages[tile.z0].asarray())
            arr = np.empty(
                (tile.z1 - tile.z0, tile.y1 - tile.y0, tile.x1 - tile.x0),
                dtype=first_plane.dtype,
            )
            arr[0] = first_plane[tile.y0 : tile.y1, tile.x0 : tile.x1]
            for out_z, page_index in enumerate(range(tile.z0 + 1, tile.z1), start=1):
                plane = np.asarray(tif.pages[page_index].asarray())
                arr[out_z] = plane[tile.y0 : tile.y1, tile.x0 : tile.x1]
        else:
            series = tif.series[0]
            arr_zyx = tiff_series_array_to_zyx(np.asarray(series.asarray()), str(getattr(series, "axes", "")))
            arr = arr_zyx[tile.z0 : tile.z1, tile.y0 : tile.y1, tile.x0 : tile.x1]

    img = sitk.GetImageFromArray(np.ascontiguousarray(arr))
    img.SetSpacing(spec.spacing_xyz)
    img.SetOrigin(
        (
            tile.x0 * spec.spacing_xyz[0],
            tile.y0 * spec.spacing_xyz[1],
            tile.z0 * spec.spacing_xyz[2],
        )
    )
    return img


def read_channel_tile_from_spec(spec: ChannelSourceSpec, tile: Tile) -> sitk.Image:
    if spec.kind == "ims":
        validate_tile_in_source(spec, tile)
        if spec.ims_spec is None:
            raise ValueError(f"Missing .ims metadata for source: {spec.path}")
        return read_ims_tile_from_spec(spec.ims_spec, tile)
    if spec.kind == "tiff":
        return read_tiff_tile_from_spec(spec, tile)
    raise ValueError(f"Unsupported source kind: {spec.kind}")


def parse_xy_coordinates_csv(
    path: Path,
    z_range: tuple[int, int] | None,
    default_z_range: tuple[int, int],
) -> tuple[int, int, int, int, int, int]:
    xs = []
    ys = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        normalized = {name.strip().lower(): name for name in reader.fieldnames}
        if "x" not in normalized or "y" not in normalized:
            raise ValueError(f"CSV must contain X,Y columns: {path}")
        for row in reader:
            xs.append(float(row[normalized["x"]]))
            ys.append(float(row[normalized["y"]]))
    if not xs or not ys:
        raise ValueError(f"CSV contains no coordinates: {path}")
    x0 = int(math.floor(min(xs)))
    x1 = int(math.ceil(max(xs)))
    y0 = int(math.floor(min(ys)))
    y1 = int(math.ceil(max(ys)))
    z0, z1 = z_range if z_range is not None else default_z_range
    if x1 <= x0 or y1 <= y0 or z1 <= z0:
        raise ValueError(f"Invalid CSV ROI bounds from {path}: {(x0, y0, z0, x1, y1, z1)}")
    return (x0, y0, z0, x1, y1, z1)


def print_roi_summary(
    roi_path: Path,
    roi_kind: str,
    roi: tuple[int, int, int, int, int, int],
    spec,
    z_defaulted: bool,
) -> None:
    x0, y0, z0, x1, y1, z1 = roi
    print(f"[apply-channels] ROI source ({roi_kind}): {roi_path}", flush=True)
    print(
        "[apply-channels] ROI parsed as full-image voxel bounds "
        f"xyz-exclusive=(x0={x0}, y0={y0}, z0={z0}, x1={x1}, y1={y1}, z1={z1})",
        flush=True,
    )
    print(
        f"[apply-channels] ROI size xyz=({x1 - x0}, {y1 - y0}, {z1 - z0})",
        flush=True,
    )
    print(
        f"[apply-channels] first source kind={spec.kind} dataset={spec.dataset_key} "
        f"shape_zyx={spec.shape_zyx} spacing_xyz={spec.spacing_xyz}",
        flush=True,
    )
    print(
        f"[apply-channels] ROI z range {'defaulted to full source depth' if z_defaulted else 'came from ROI/z-range'}",
        flush=True,
    )
    z_size, y_size, x_size = spec.shape_zyx
    if x1 > x_size or y1 > y_size or z1 > z_size:
        print(
            "[apply-channels] WARNING ROI extends beyond first source bounds "
            f"xyz-size=({x_size}, {y_size}, {z_size})",
            flush=True,
        )
    if x0 == 0 and y0 == 0:
        print(
            "[apply-channels] WARNING ROI starts at x=0/y=0; if this ROI was saved from a cropped image, "
            "its coordinates may be crop-local rather than full-source coordinates.",
            flush=True,
        )


def print_sitk_stats(label: str, img: "sitk.Image") -> None:
    arr = sitk.GetArrayViewFromImage(img)
    print_array_stats(label, arr)


def print_ants_stats(label: str, img) -> None:
    arr = img.numpy()
    print_array_stats(label, arr)


def print_array_stats(label: str, arr: np.ndarray) -> None:
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        print(f"{label}: no finite voxels", flush=True)
        return
    nonzero = int(np.count_nonzero(arr))
    total = int(arr.size)
    pct = np.percentile(finite, [0, 1, 50, 99, 99.9, 100])
    print(
        f"{label}: dtype={arr.dtype} shape={arr.shape} nonzero={nonzero}/{total} "
        f"p0={pct[0]:.6g} p1={pct[1]:.6g} p50={pct[2]:.6g} "
        f"p99={pct[3]:.6g} p99.9={pct[4]:.6g} max={pct[5]:.6g}",
        flush=True,
    )
    if nonzero == 0:
        print(f"{label}: WARNING all voxels are zero/blank", flush=True)


def record_timing(timings: list[dict], stage: str, start: float, **extra) -> float:
    elapsed = time.time() - start
    entry = {"stage": stage, "seconds": elapsed, "duration": format_duration(elapsed)}
    entry.update(extra)
    timings.append(entry)
    print(f"[timing] {stage}: {format_duration(elapsed)}", flush=True)
    return elapsed


def ants_to_array_for_tiff(img, pixeltype: str, uint16_scaling: str) -> np.ndarray:
    if pixeltype == "uint16":
        sitk_img = ants_to_sitk_uint16(img, uint16_scaling)
        return sitk.GetArrayFromImage(sitk_img)
    if pixeltype == "float32":
        arr = img.numpy()
        return np.asarray(arr.T, dtype=np.float32)
    raise ValueError(f"Unsupported TIFF pixel type: {pixeltype}")


def registered_output_paths(output_dir: Path, prefix: str, registered_format: str) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    if registered_format in {"nrrd", "both"}:
        paths["nrrd"] = output_dir / f"{prefix}_roi_registered.nrrd"
    if registered_format in {"tif", "both"}:
        paths["tif"] = output_dir / f"{prefix}_roi_registered.tif"
    return paths


def write_registered_outputs(
    img,
    paths: dict[str, Path],
    pixeltype: str,
    progress_interval: float,
    uint16_scaling: str,
) -> None:
    if "nrrd" in paths:
        write_ants_image(img, paths["nrrd"], pixeltype, progress_interval, uint16_scaling)
        print(f"[apply-channels] wrote registered NRRD: {paths['nrrd']}", flush=True)
    if "tif" in paths:
        paths["tif"].parent.mkdir(parents=True, exist_ok=True)
        print(f"[apply-channels] writing registered TIFF stack: {paths['tif']}", flush=True)
        with progress_heartbeat("write registered tiff", progress_interval):
            arr = ants_to_array_for_tiff(img, pixeltype, uint16_scaling)
            tifffile.imwrite(paths["tif"], arr, imagej=True, metadata={"axes": "ZYX"})
        print(f"[apply-channels] wrote registered TIFF stack: {paths['tif']}", flush=True)


def copy_geometry_from_reference(img: "sitk.Image", reference: "sitk.Image", label: str) -> "sitk.Image":
    if img.GetSize() != reference.GetSize():
        print(
            f"[apply-channels] WARNING not copying {label} geometry because size differs: "
            f"channel={img.GetSize()} reference={reference.GetSize()}",
            flush=True,
        )
        return img
    img.CopyInformation(reference)
    print(
        f"[apply-channels] copied {label} geometry: "
        f"origin={img.GetOrigin()} spacing={img.GetSpacing()} direction={img.GetDirection()}",
        flush=True,
    )
    return img


def tile_from_roi(roi: tuple[int, int, int, int, int, int]) -> Tile:
    x0, y0, z0, x1, y1, z1 = roi
    return Tile(z0=z0, z1=z1, y0=y0, y1=y1, x0=x0, x1=x1)


def padded_tile(base_tile: Tile, spec, padding_xy: int) -> Tile:
    padding_xy = max(0, int(padding_xy))
    z_size, y_size, x_size = spec.shape_zyx
    tile = Tile(
        z0=base_tile.z0,
        z1=base_tile.z1,
        y0=max(0, base_tile.y0 - padding_xy),
        y1=min(y_size, base_tile.y1 + padding_xy),
        x0=max(0, base_tile.x0 - padding_xy),
        x1=min(x_size, base_tile.x1 + padding_xy),
    )
    if tile.x1 <= tile.x0 or tile.y1 <= tile.y0 or tile.z1 <= tile.z0:
        raise ValueError(f"Invalid padded tile produced from {base_tile}: {tile}")
    return tile


def set_tile_geometry_from_reference(
    img: "sitk.Image",
    reference: "sitk.Image",
    tile: Tile,
    base_tile: Tile,
    label: str,
) -> "sitk.Image":
    spacing = tuple(float(v) for v in reference.GetSpacing())
    origin = np.asarray(reference.GetOrigin(), dtype=float)
    direction = np.asarray(reference.GetDirection(), dtype=float).reshape((3, 3))
    delta_index_xyz = np.asarray(
        (
            tile.x0 - base_tile.x0,
            tile.y0 - base_tile.y0,
            tile.z0 - base_tile.z0,
        ),
        dtype=float,
    )
    delta_physical = direction @ (delta_index_xyz * np.asarray(spacing, dtype=float))
    img.SetSpacing(spacing)
    img.SetDirection(reference.GetDirection())
    img.SetOrigin(tuple(float(v) for v in origin + delta_physical))
    print(
        f"[apply-channels] set {label} tile geometry: tile={tile} base={base_tile} "
        f"origin={img.GetOrigin()} spacing={img.GetSpacing()}",
        flush=True,
    )
    return img


def apply_geometry_to_ants_image(img, reference: "sitk.Image", label: str):
    reference_size = tuple(int(v) for v in reference.GetSize())
    image_size = tuple(int(v) for v in img.shape)
    if image_size != reference_size:
        print(
            f"[apply-channels] WARNING not setting {label} ANTs geometry because size differs: "
            f"channel={image_size} reference={reference_size}",
            flush=True,
        )
        return img

    origin = tuple(float(v) for v in reference.GetOrigin())
    spacing = tuple(float(v) for v in reference.GetSpacing())
    img.set_origin(origin)
    img.set_spacing(spacing)
    direction_flat = tuple(float(v) for v in reference.GetDirection())
    dim = len(reference_size)
    try:
        img.set_direction(np.asarray(direction_flat, dtype=float).reshape((dim, dim)))
        direction_text = f" direction={direction_flat}"
    except Exception as exc:  # pragma: no cover - depends on ANTsPy build details
        direction_text = f" direction=unchanged ({exc})"
    print(
        f"[apply-channels] set {label} ANTs geometry: origin={origin} spacing={spacing}{direction_text}",
        flush=True,
    )
    return img


def maybe_downsample_for_intra(img, factor: float, label: str, progress_interval: float):
    if factor <= 1.0:
        print(f"[intra-align] {label}: using full resolution for intra-cycle registration", flush=True)
        return img
    print(f"[intra-align] {label}: downsampling by spacing factor {factor:g}", flush=True)
    with progress_heartbeat(f"intra-align downsample {label}", progress_interval):
        ds = downsample_by_spacing(img, factor)
    describe_image(f"[intra-align] {label} downsampled", ds)
    return ds


def intra_cycle_align_to_dapi(
    channel: int,
    moving_crop,
    moving_reference_path: Path,
    transform_dir: Path,
    transform_type: str,
    downsample_factor: float,
    qc_dir: Path | None,
    open_qc: bool,
    qc_max_panel_side: int,
    progress_interval: float,
):
    start = time.time()
    intra_dir = transform_dir.parent / "intra_cycle_channel_alignment"
    intra_dir.mkdir(parents=True, exist_ok=True)
    outprefix = intra_dir / f"ch{channel}_to_ch0_"

    print(
        f"[intra-align] channel {channel}: loading moving DAPI reference for channel-to-DAPI alignment: "
        f"{moving_reference_path}",
        flush=True,
    )
    with progress_heartbeat(f"intra-align load ch0 reference for ch{channel}", progress_interval):
        moving_reference = ants.image_read(str(moving_reference_path))
    describe_image(f"[intra-align] channel {channel} moving DAPI reference", moving_reference)
    describe_image(f"[intra-align] channel {channel} input ROI", moving_crop)

    fixed_ds = maybe_downsample_for_intra(
        moving_reference,
        downsample_factor,
        f"channel {channel} DAPI reference",
        progress_interval,
    )
    moving_ds = maybe_downsample_for_intra(
        moving_crop,
        downsample_factor,
        f"channel {channel}",
        progress_interval,
    )

    print(
        f"[intra-align] channel {channel}: running ANTs {transform_type} to align channel -> moving DAPI",
        flush=True,
    )
    with progress_heartbeat(f"intra-align register ch{channel}", progress_interval):
        tx = ants.registration(
            fixed=fixed_ds,
            moving=moving_ds,
            type_of_transform=transform_type,
            outprefix=str(outprefix),
            aff_metric="mattes",
            singleprecision=True,
            verbose=True,
        )
    transform_json = intra_dir / f"ch{channel}_to_ch0_fwdtransforms.json"
    transform_json.write_text(json.dumps({"fwdtransforms": tx["fwdtransforms"]}, indent=2), encoding="utf-8")
    print(f"[intra-align] channel {channel}: wrote transform manifest: {transform_json}", flush=True)

    if qc_dir is not None:
        print(f"[intra-align] channel {channel}: building downsampled before/after QC", flush=True)
        with progress_heartbeat(f"intra-align qc apply ch{channel}", progress_interval):
            aligned_ds = ants.apply_transforms(
                fixed=fixed_ds,
                moving=moving_ds,
                transformlist=tx["fwdtransforms"],
                interpolator="linear",
                singleprecision=True,
            )
        export_center_slice_qc(
            fixed_ds,
            moving_ds,
            aligned_ds,
            qc_dir,
            f"intra_cycle_ch{channel}_to_ch0",
            open_qc,
            qc_max_panel_side,
        )
        del aligned_ds

    if fixed_ds is not moving_reference:
        del fixed_ds
    if moving_ds is not moving_crop:
        del moving_ds
    gc.collect()

    print(f"[intra-align] channel {channel}: applying intra-cycle channel-to-DAPI transform", flush=True)
    with progress_heartbeat(f"intra-align apply ch{channel}", progress_interval):
        aligned = ants.apply_transforms(
            fixed=moving_reference,
            moving=moving_crop,
            transformlist=tx["fwdtransforms"],
            interpolator="linear",
            singleprecision=True,
        )
    print_ants_stats(f"[intra-align] channel {channel} aligned-to-DAPI stats", aligned)
    del moving_reference
    gc.collect()
    print(
        f"[intra-align] channel {channel}: finished in {format_duration(time.time() - start)}",
        flush=True,
    )
    return aligned


def write_sitk_temp_for_ants(img: "sitk.Image", path: Path, label: str, progress_interval: float):
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[pre-crop-align] writing temporary {label}: {path}", flush=True)
    with progress_heartbeat(f"pre-crop write {label}", progress_interval):
        sitk.WriteImage(img, str(path))
    print(f"[pre-crop-align] loading temporary {label} with ANTs: {path}", flush=True)
    with progress_heartbeat(f"pre-crop load {label}", progress_interval):
        ants_img = ants.image_read(str(path))
    return ants_img, path


def pre_crop_intra_cycle_align_to_dapi(
    channel: int,
    channel_input_path: Path,
    channel_spec,
    dapi_input_path: Path,
    dapi_spec,
    base_tile: Tile,
    moving_reference_path: Path,
    moving_reference_sitk: "sitk.Image",
    transform_dir: Path,
    raw_output_dir: Path,
    prefix: str,
    transform_type: str,
    downsample_factor: float,
    padding_xy: int,
    save_aligned_roi: bool,
    keep_temp: bool,
    qc_dir: Path | None,
    open_qc: bool,
    qc_max_panel_side: int,
    progress_interval: float,
):
    start = time.time()
    intra_dir = transform_dir.parent / "intra_cycle_channel_alignment"
    tmp_dir = intra_dir / "pre_crop_tmp"
    intra_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    padded = padded_tile(base_tile, channel_spec, padding_xy)
    outprefix = intra_dir / f"ch{channel}_to_ch0_precrop_"

    print(
        f"[pre-crop-align] channel {channel}: correcting on padded source data before ROI recrop; "
        f"base_tile={base_tile} padded_tile={padded}",
        flush=True,
    )
    print(f"[pre-crop-align] channel {channel}: DAPI source={dapi_input_path}", flush=True)
    print(f"[pre-crop-align] channel {channel}: moving source={channel_input_path}", flush=True)

    with progress_heartbeat(f"pre-crop extract ch0 padded for ch{channel}", progress_interval):
        dapi_padded_sitk = read_channel_tile_from_spec(dapi_spec, padded)
    dapi_padded_sitk = set_tile_geometry_from_reference(
        dapi_padded_sitk,
        moving_reference_sitk,
        padded,
        base_tile,
        f"channel {channel} padded DAPI",
    )
    print_sitk_stats(f"[pre-crop-align] channel {channel} padded DAPI stats", dapi_padded_sitk)

    with progress_heartbeat(f"pre-crop extract ch{channel} padded", progress_interval):
        moving_padded_sitk = read_channel_tile_from_spec(channel_spec, padded)
    moving_padded_sitk = set_tile_geometry_from_reference(
        moving_padded_sitk,
        moving_reference_sitk,
        padded,
        base_tile,
        f"channel {channel} padded moving",
    )
    print_sitk_stats(f"[pre-crop-align] channel {channel} padded moving stats", moving_padded_sitk)

    dapi_padded, dapi_temp_path = write_sitk_temp_for_ants(
        dapi_padded_sitk,
        tmp_dir / f"ch{channel}_ch0_padded.nrrd",
        f"ch{channel} padded DAPI",
        progress_interval,
    )
    moving_padded, moving_temp_path = write_sitk_temp_for_ants(
        moving_padded_sitk,
        tmp_dir / f"ch{channel}_moving_padded.nrrd",
        f"ch{channel} padded moving",
        progress_interval,
    )
    del dapi_padded_sitk, moving_padded_sitk
    gc.collect()

    fixed_ds = maybe_downsample_for_intra(
        dapi_padded,
        downsample_factor,
        f"channel {channel} padded DAPI",
        progress_interval,
    )
    moving_ds = maybe_downsample_for_intra(
        moving_padded,
        downsample_factor,
        f"channel {channel} padded moving",
        progress_interval,
    )

    print(
        f"[pre-crop-align] channel {channel}: running ANTs {transform_type} on padded source crops",
        flush=True,
    )
    with progress_heartbeat(f"pre-crop register ch{channel}", progress_interval):
        tx = ants.registration(
            fixed=fixed_ds,
            moving=moving_ds,
            type_of_transform=transform_type,
            outprefix=str(outprefix),
            aff_metric="mattes",
            singleprecision=True,
            verbose=True,
        )
    transform_json = intra_dir / f"ch{channel}_to_ch0_precrop_fwdtransforms.json"
    transform_json.write_text(json.dumps({"fwdtransforms": tx["fwdtransforms"]}, indent=2), encoding="utf-8")
    print(f"[pre-crop-align] channel {channel}: wrote transform manifest: {transform_json}", flush=True)

    if qc_dir is not None:
        print(f"[pre-crop-align] channel {channel}: building padded before/after QC", flush=True)
        with progress_heartbeat(f"pre-crop qc apply ch{channel}", progress_interval):
            aligned_ds = ants.apply_transforms(
                fixed=fixed_ds,
                moving=moving_ds,
                transformlist=tx["fwdtransforms"],
                interpolator="linear",
                singleprecision=True,
            )
        export_center_slice_qc(
            fixed_ds,
            moving_ds,
            aligned_ds,
            qc_dir,
            f"pre_crop_intra_cycle_ch{channel}_to_ch0",
            open_qc,
            qc_max_panel_side,
        )
        del aligned_ds

    if fixed_ds is not dapi_padded:
        del fixed_ds
    if moving_ds is not moving_padded:
        del moving_ds
    gc.collect()

    print(
        f"[pre-crop-align] channel {channel}: recropping aligned channel into original ROI geometry: "
        f"{moving_reference_path}",
        flush=True,
    )
    with progress_heartbeat(f"pre-crop load recrop reference ch{channel}", progress_interval):
        recrop_reference = ants.image_read(str(moving_reference_path))
    with progress_heartbeat(f"pre-crop recrop/apply ch{channel}", progress_interval):
        aligned_roi = ants.apply_transforms(
            fixed=recrop_reference,
            moving=moving_padded,
            transformlist=tx["fwdtransforms"],
            interpolator="linear",
            singleprecision=True,
        )
    print_ants_stats(f"[pre-crop-align] channel {channel} recropped aligned ROI stats", aligned_roi)

    prealigned_path = None
    if save_aligned_roi:
        prealigned_path = raw_output_dir / f"{prefix}_roi_precrop_intracycle_aligned.nrrd"
        print(f"[pre-crop-align] writing recropped aligned ROI for channel {channel}: {prealigned_path}", flush=True)
        with progress_heartbeat(f"pre-crop write aligned roi ch{channel}", progress_interval):
            ants.image_write(aligned_roi, str(prealigned_path))
        print(f"[pre-crop-align] channel {channel}: wrote recropped aligned ROI: {prealigned_path}", flush=True)
    else:
        print(
            f"[pre-crop-align] channel {channel}: keeping recropped aligned ROI in memory only",
            flush=True,
        )

    del dapi_padded, moving_padded, recrop_reference
    if keep_temp:
        print(
            f"[pre-crop-align] keeping temporary padded NRRDs: {dapi_temp_path} ; {moving_temp_path}",
            flush=True,
        )
    else:
        for temp_path in (dapi_temp_path, moving_temp_path):
            try:
                temp_path.unlink()
                print(f"[pre-crop-align] removed temporary file: {temp_path}", flush=True)
            except FileNotFoundError:
                pass
    gc.collect()
    print(
        f"[pre-crop-align] channel {channel}: finished pre-crop correction in "
        f"{format_duration(time.time() - start)}",
        flush=True,
    )
    return aligned_roi, prealigned_path


def main() -> int:
    run_start = time.time()
    args = build_parser().parse_args()
    configure_ants_runtime(args.threads)
    global ants
    import ants as ants_module

    ants = ants_module

    transform_dir = Path(args.transform_dir).expanduser().resolve()
    manifest = load_transform_manifest(transform_dir)
    input_map = parse_input_ims_map(args.input_ims_map) if args.input_ims_map else {}
    raw_output_dir = Path(args.raw_output_dir).expanduser().resolve()
    registered_output_dir = (
        Path(args.registered_output_dir).expanduser().resolve()
        if args.registered_output_dir
        else (registration_run_dir(transform_dir) / "outputs").resolve()
    )
    raw_output_dir.mkdir(parents=True, exist_ok=True)
    registered_output_dir.mkdir(parents=True, exist_ok=True)
    timing_log: list[dict] = []
    qc_dir = None
    if args.qc:
        qc_dir = (
            Path(args.qc_output_dir).expanduser().resolve()
            if args.qc_output_dir
            else (registration_run_dir(transform_dir) / "qc" / "intra_cycle_channel_alignment").resolve()
        )
        print(f"[apply-channels] intra-cycle QC output dir: {qc_dir}", flush=True)

    first_channel = args.channels[0]
    first_path = channel_ims_path(args, first_channel, input_map)
    first_spec = build_channel_spec(args, first_path, first_channel)
    full_z_range = (0, first_spec.shape_zyx[0])
    if args.roi_csv:
        roi_path = Path(args.roi_csv).expanduser().resolve()
        roi_kind = "csv"
        roi = parse_xy_coordinates_csv(roi_path, args.z_range, full_z_range)
    elif args.roi_imagej:
        roi_path = Path(args.roi_imagej).expanduser().resolve()
        roi_kind = "imagej"
        roi = parse_imagej_roi(
            roi_path,
            args.z_range,
            default_z_range=full_z_range,
        )
    else:
        raise ValueError("Provide either --roi-csv or --roi-imagej")
    x0, y0, z0, x1, y1, z1 = roi
    tile = tile_from_roi(roi)

    print(f"[apply-channels] ROI xyz exclusive: {roi}", flush=True)
    print_roi_summary(
        roi_path,
        roi_kind,
        roi,
        first_spec,
        z_defaulted=bool(args.z_range is None),
    )
    print(f"[apply-channels] raw output dir: {raw_output_dir}", flush=True)
    print(f"[apply-channels] registered output dir: {registered_output_dir}", flush=True)
    print(f"[apply-channels] loading fixed crop: {args.fixed_crop}", flush=True)
    stage_start = time.time()
    with progress_heartbeat("apply-channels load fixed crop", args.progress_interval):
        fixed_crop = ants.image_read(str(Path(args.fixed_crop).expanduser().resolve()))
    record_timing(timing_log, "load fixed crop", stage_start)
    moving_reference_path = None
    moving_reference_sitk = None
    if args.moving_crop_reference:
        moving_reference_path = Path(args.moving_crop_reference).expanduser().resolve()
        print(f"[apply-channels] loading moving crop geometry reference: {moving_reference_path}", flush=True)
        stage_start = time.time()
        with progress_heartbeat("apply-channels load moving geometry reference", args.progress_interval):
            moving_reference_sitk = sitk.ReadImage(str(moving_reference_path))
        record_timing(timing_log, "load moving crop geometry reference", stage_start)
        print(
            f"[apply-channels] moving reference size={moving_reference_sitk.GetSize()} "
            f"origin={moving_reference_sitk.GetOrigin()} spacing={moving_reference_sitk.GetSpacing()}",
            flush=True,
        )

    outputs = []
    intra_channels = set(args.intra_cycle_align_channels)
    pre_crop_intra_channels = set(args.pre_crop_intra_cycle_align_channels)
    if intra_channels and moving_reference_path is None:
        raise ValueError("--intra-cycle-align-channels requires --moving-crop-reference")
    if pre_crop_intra_channels and moving_reference_path is None:
        raise ValueError("--pre-crop-intra-cycle-align-channels requires --moving-crop-reference")
    if intra_channels:
        print(
            f"[apply-channels] intra-cycle channel-to-DAPI alignment enabled for channels: "
            f"{sorted(intra_channels)}",
            flush=True,
        )
    dapi_input_path = None
    dapi_spec = None
    if pre_crop_intra_channels:
        dapi_input_path = channel_ims_path(args, 0, input_map)
        dapi_spec = build_channel_spec(args, dapi_input_path, 0)
        print(
            f"[apply-channels] pre-crop intra-cycle channel-to-DAPI alignment enabled for channels: "
            f"{sorted(pre_crop_intra_channels)} padding_xy={args.intra_padding_xy}",
            flush=True,
        )
        print(
            f"[apply-channels] pre-crop alignment DAPI source: {dapi_input_path} "
            f"kind={dapi_spec.kind} dataset={dapi_spec.dataset_key} shape_zyx={dapi_spec.shape_zyx}",
            flush=True,
        )
    for channel in args.channels:
        channel_start = time.time()
        input_path = channel_ims_path(args, channel, input_map)
        spec = first_spec if channel == first_channel and input_path == first_path else build_channel_spec(args, input_path, channel)
        prefix = args.prefix or stem_for_output(input_path)
        if not args.prefix and not prefix.lower().endswith(f"_ch{channel}"):
            prefix = f"{prefix}_ch{channel}"
        raw_path = raw_output_dir / f"{prefix}_roi.nrrd"
        registered_paths = registered_output_paths(registered_output_dir, prefix, args.registered_format)
        primary_registered_path = registered_paths.get("tif") or registered_paths.get("nrrd")
        manifest_raw_path = raw_path

        if channel in pre_crop_intra_channels:
            stage_start = time.time()
            moving_crop, prealigned_path = pre_crop_intra_cycle_align_to_dapi(
                channel,
                input_path,
                spec,
                dapi_input_path,
                dapi_spec,
                tile,
                moving_reference_path,
                moving_reference_sitk,
                transform_dir,
                raw_output_dir,
                prefix,
                args.intra_cycle_transform,
                args.intra_cycle_downsample_factor,
                args.intra_padding_xy,
                args.save_pre_crop_aligned_roi,
                args.keep_pre_crop_temp,
                qc_dir,
                args.open_qc,
                args.qc_max_panel_side,
                args.progress_interval,
            )
            record_timing(timing_log, f"channel {channel} pre-crop intra-cycle alignment", stage_start, channel=channel)
            manifest_raw_path = prealigned_path or raw_path
        elif raw_path.exists() and raw_path.stat().st_size > 0 and not args.no_reuse_raw_roi:
            print(f"[apply-channels] reusing existing raw channel {channel} ROI: {raw_path}", flush=True)
            print(f"[apply-channels] applying transform to channel {channel}", flush=True)
            stage_start = time.time()
            with progress_heartbeat(f"load channel {channel} roi", args.progress_interval):
                moving_crop = ants.image_read(str(raw_path))
            record_timing(timing_log, f"channel {channel} load existing ROI", stage_start, channel=channel)
        else:
            print(f"[apply-channels] extracting channel {channel} ROI from: {input_path}", flush=True)
            print(
                f"[apply-channels] channel {channel} source_kind={spec.kind} dataset={spec.dataset_key} "
                f"shape_zyx={spec.shape_zyx} spacing_xyz={spec.spacing_xyz}",
                flush=True,
            )
            stage_start = time.time()
            with progress_heartbeat(f"extract channel {channel} roi", args.progress_interval):
                raw_img = read_channel_tile_from_spec(spec, tile)
            record_timing(timing_log, f"channel {channel} extract ROI from source", stage_start, channel=channel)
            if moving_reference_sitk is not None:
                raw_img = copy_geometry_from_reference(raw_img, moving_reference_sitk, "moving crop reference")
            print_sitk_stats(f"[apply-channels] raw channel {channel} ROI stats", raw_img)
            stage_start = time.time()
            sitk.WriteImage(raw_img, str(raw_path))
            record_timing(timing_log, f"channel {channel} write raw ROI", stage_start, channel=channel)
            print(f"[apply-channels] wrote raw channel {channel} ROI: {raw_path}", flush=True)
            del raw_img
            gc.collect()
            print(f"[apply-channels] applying transform to channel {channel}", flush=True)
            stage_start = time.time()
            with progress_heartbeat(f"load channel {channel} roi", args.progress_interval):
                moving_crop = ants.image_read(str(raw_path))
            record_timing(timing_log, f"channel {channel} load extracted ROI", stage_start, channel=channel)
        if moving_reference_sitk is not None and channel not in pre_crop_intra_channels:
            moving_crop = apply_geometry_to_ants_image(moving_crop, moving_reference_sitk, "moving crop reference")
        print_ants_stats(f"[apply-channels] loaded channel {channel} ROI stats", moving_crop)
        if channel in intra_channels and channel not in pre_crop_intra_channels:
            stage_start = time.time()
            aligned_crop = intra_cycle_align_to_dapi(
                channel,
                moving_crop,
                moving_reference_path,
                transform_dir,
                args.intra_cycle_transform,
                args.intra_cycle_downsample_factor,
                qc_dir,
                args.open_qc,
                args.qc_max_panel_side,
                args.progress_interval,
            )
            del moving_crop
            gc.collect()
            moving_crop = aligned_crop
            record_timing(timing_log, f"channel {channel} intra-cycle alignment", stage_start, channel=channel)
        stage_start = time.time()
        with progress_heartbeat(f"apply channel {channel} transform", args.progress_interval):
            registered = ants.apply_transforms(
                fixed=fixed_crop,
                moving=moving_crop,
                transformlist=manifest["fwdtransforms"],
                interpolator=args.interpolator,
                singleprecision=True,
            )
        record_timing(timing_log, f"channel {channel} apply cycle transform", stage_start, channel=channel)
        print_ants_stats(f"[apply-channels] registered channel {channel} stats before write", registered)
        stage_start = time.time()
        write_registered_outputs(
            registered,
            registered_paths,
            args.output_pixeltype,
            args.progress_interval,
            args.uint16_scaling,
        )
        record_timing(timing_log, f"channel {channel} write registered outputs", stage_start, channel=channel)
        print(f"[apply-channels] wrote registered channel {channel}: {primary_registered_path}", flush=True)
        channel_elapsed = record_timing(timing_log, f"channel {channel} total", channel_start, channel=channel)
        outputs.append(
            {
                "channel": channel,
                "input_ims": str(input_path),
                "input_source": str(input_path),
                "input_source_kind": spec.kind,
                "raw_roi": str(manifest_raw_path),
                "pre_crop_intracycle_aligned": bool(channel in pre_crop_intra_channels),
                "pre_crop_aligned_roi_saved": str(prealigned_path) if channel in pre_crop_intra_channels and prealigned_path else "",
                "registered": str(primary_registered_path),
                "registered_nrrd": str(registered_paths.get("nrrd", "")),
                "registered_tif": str(registered_paths.get("tif", "")),
                "elapsed_seconds": channel_elapsed,
                "elapsed": format_duration(channel_elapsed),
            }
        )
        del moving_crop, registered
        gc.collect()

    total_elapsed = record_timing(timing_log, "apply channels total", run_start)
    timing_path = registered_output_dir / "stage_timing_apply_channels.json"
    timing_path.write_text(json.dumps({"timings": timing_log, "total_seconds": total_elapsed}, indent=2), encoding="utf-8")
    print(f"[timing] wrote apply-channel timing log: {timing_path}", flush=True)

    manifest_path = registered_output_dir / "registered_channel_outputs.json"
    manifest_path.write_text(
        json.dumps(
            {
                "roi_imagej": str(Path(args.roi_imagej).expanduser().resolve()) if args.roi_imagej else "",
                "roi_csv": str(Path(args.roi_csv).expanduser().resolve()) if args.roi_csv else "",
                "roi_source": roi_kind,
                "roi_xyz_exclusive": {"x0": x0, "y0": y0, "z0": z0, "x1": x1, "y1": y1, "z1": z1},
                "fixed_crop": str(Path(args.fixed_crop).expanduser().resolve()),
                "transform_dir": str(transform_dir),
                "fwdtransforms": manifest["fwdtransforms"],
                "timing_log": str(timing_path),
                "timings": timing_log,
                "outputs": outputs,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"[apply-channels] wrote manifest: {manifest_path}", flush=True)
    del fixed_crop
    gc.collect()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
