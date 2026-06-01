#!/usr/bin/env python3
"""Whole-.ims global registration + tiled QuickSyN refinement for very large channel volumes."""

from __future__ import annotations

import argparse
import gc
import html
import json
import math
import os
import sys
import tempfile
import time
import webbrowser
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile

from ants_roi_quicksyn import (
    compute_similarity_metrics_from_images,
    configure_ants_runtime,
    export_center_slice_qc,
    format_duration,
    progress_heartbeat,
    write_metrics_json,
)
from extract_ims_roi_channels import parse_channels
from cyants_io import (
    Tile,
    build_ims_volume_spec,
    generate_tiles,
    parse_spacing_override,
    raised_cosine_weights,
    read_ims_downsampled_from_spec,
    read_ims_tile_from_spec,
)

TIFF_SUFFIXES = {".tif", ".tiff"}


@dataclass(frozen=True)
class SourceSpec:
    path: Path
    kind: str
    channel: int
    dataset_key: str
    axis_order: str
    shape_zyx: tuple[int, int, int]
    spacing_xyz: tuple[float, float, float]
    ims_spec: object | None = None


def parse_triplet(value: str, label: str) -> tuple[int, int, int]:
    parts = [int(v.strip()) for v in value.split(",")]
    if len(parts) != 3 or any(v <= 0 for v in parts):
        raise ValueError(f"--{label} must be z,y,x positive integers")
    return tuple(parts)


def parse_overlap(value: str) -> tuple[int, int, int]:
    parts = [int(v.strip()) for v in value.split(",")]
    if len(parts) != 3 or any(v < 0 for v in parts):
        raise ValueError("--overlap must be z,y,x non-negative integers")
    return tuple(parts)


def choose_xy_grid(tile_count: int, shape_zyx: tuple[int, int, int]) -> tuple[int, int]:
    if tile_count < 1:
        raise ValueError("--tiles must be >= 1")
    _, y_size, x_size = shape_zyx
    aspect = x_size / max(1, y_size)
    candidates = []
    for rows in range(1, int(math.sqrt(tile_count)) + 1):
        if tile_count % rows:
            continue
        cols = tile_count // rows
        grid_aspect = cols / rows
        score = abs(math.log(max(grid_aspect, 1e-9) / max(aspect, 1e-9)))
        candidates.append((score, rows, cols))
    _, rows, cols = min(candidates)
    return rows, cols


def auto_axis_tile(length: int, count: int, overlap_fraction: float | None, overlap_pixels: int = 0) -> tuple[int, int]:
    if count <= 1:
        return length, 0
    if overlap_fraction is not None:
        if not (0.0 <= overlap_fraction < 1.0):
            raise ValueError("fractional --overlap with --tiles must be >= 0 and < 1")
        denominator = count - (count - 1) * overlap_fraction
        tile = int(math.ceil(length / denominator))
        overlap = int(round(tile * overlap_fraction))
    else:
        overlap = int(overlap_pixels)
        tile = int(math.ceil((length + (count - 1) * overlap) / count))
    overlap = max(0, min(overlap, tile - 1))
    return min(tile, length), overlap


def auto_tile_schedule(args, shape_zyx: tuple[int, int, int]) -> tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]:
    rows, cols = choose_xy_grid(args.tiles, shape_zyx)
    overlap_arg = str(args.overlap).strip()
    overlap_fraction: float | None = None
    overlap_px = (0, 0, 0)
    if "," in overlap_arg:
        overlap_px = parse_overlap(overlap_arg)
    else:
        value = float(overlap_arg)
        if value < 1.0:
            overlap_fraction = value
        else:
            overlap_px = (0, int(round(value)), int(round(value)))
    z_size, y_size, x_size = shape_zyx
    tile_y, overlap_y = auto_axis_tile(y_size, rows, overlap_fraction, overlap_px[1])
    tile_x, overlap_x = auto_axis_tile(x_size, cols, overlap_fraction, overlap_px[2])
    return (z_size, tile_y, tile_x), (0, overlap_y, overlap_x), (1, rows, cols)


def tile_schedule(args, shape_zyx: tuple[int, int, int]) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    if args.tiles > 0:
        tile_size, overlap, grid = auto_tile_schedule(args, shape_zyx)
        print(
            f"[tile-schedule] auto tiles={args.tiles} grid_zyx={grid} "
            f"tile_size_zyx={tile_size} overlap_zyx={overlap}",
            flush=True,
        )
        return tile_size, overlap
    tile_size = parse_triplet(args.tile, "tile")
    overlap = parse_overlap(args.overlap)
    if args.full_z:
        tile_size = (shape_zyx[0], tile_size[1], tile_size[2])
        overlap = (0, overlap[1], overlap[2])
    return tile_size, overlap


def parse_source_map(value: str) -> dict[int, Path]:
    mapping: dict[int, Path] = {}
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Invalid --moving-source-map entry: {item}")
        channel_s, path_s = item.split("=", 1)
        path = Path(path_s.strip().strip('"')).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Moving source override does not exist for channel {channel_s}: {path}")
        if path.suffix.lower() not in {".ims", ".tif", ".tiff"}:
            raise ValueError(f"Moving source override must be .ims, .tif, or .tiff: {path}")
        mapping[int(channel_s.strip())] = path
    return mapping


def tile_to_dict(tile: Tile) -> dict[str, int]:
    return {"z0": tile.z0, "z1": tile.z1, "y0": tile.y0, "y1": tile.y1, "x0": tile.x0, "x1": tile.x1}


def tile_axis_segments(length: int, tile: int, starts: list[int]) -> list[dict[str, int]]:
    return [{"start": start, "end": min(start + tile, length), "size": min(start + tile, length) - start} for start in starts]


def adjacent_overlaps(segments: list[dict[str, int]]) -> list[dict[str, int]]:
    overlaps = []
    for index, (left, right) in enumerate(zip(segments, segments[1:]), start=1):
        start = max(left["start"], right["start"])
        end = min(left["end"], right["end"])
        if end > start:
            overlaps.append({"between": index, "start": start, "end": end, "size": end - start})
    return overlaps


def write_tile_map_svg(
    svg_path: Path,
    spec: SourceSpec,
    tile_size: tuple[int, int, int],
    requested_overlap: tuple[int, int, int],
    tiles: list[Tile],
    x_overlaps: list[dict[str, int]],
    y_overlaps: list[dict[str, int]],
) -> None:
    z_size, y_size, x_size = spec.shape_zyx
    max_w = 1500.0
    max_h = 760.0
    scale = min(max_w / x_size, max_h / y_size)
    image_w = x_size * scale
    image_h = y_size * scale
    margin_x = 48.0
    image_y = 150.0
    canvas_w = image_w + margin_x * 2
    canvas_h = image_y + image_h + 160.0 + min(len(tiles), 12) * 18.0

    def sx(value: int | float) -> float:
        return margin_x + float(value) * scale

    def sy(value: int | float) -> float:
        return image_y + float(value) * scale

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{canvas_w:.0f}" height="{canvas_h:.0f}" viewBox="0 0 {canvas_w:.0f} {canvas_h:.0f}">',
        "<style>"
        "text{font-family:Arial,sans-serif;fill:#1f2933}"
        ".small{font-size:14px}.label{font-size:16px;font-weight:700}.title{font-size:24px;font-weight:700}"
        ".mono{font-family:Consolas,Menlo,monospace;font-size:13px}"
        "</style>",
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="48" y="42" class="title">CyANTs Tile Map</text>',
        f'<text x="48" y="70" class="small">source: {html.escape(str(spec.path))}</text>',
        (
            f'<text x="48" y="94" class="small">shape Z,Y,X={spec.shape_zyx} | '
            f'tile Z,Y,X={tile_size} | requested overlap Z,Y,X={requested_overlap} | tiles={len(tiles)}</text>'
        ),
        f'<text x="48" y="118" class="small">spacing X,Y,Z={spec.spacing_xyz}</text>',
        f'<rect x="{margin_x:.2f}" y="{image_y:.2f}" width="{image_w:.2f}" height="{image_h:.2f}" fill="#ffffff" stroke="#111827" stroke-width="2"/>',
    ]

    for overlap in x_overlaps:
        x = sx(overlap["start"])
        w = overlap["size"] * scale
        lines.append(
            f'<rect x="{x:.2f}" y="{image_y:.2f}" width="{w:.2f}" height="{image_h:.2f}" fill="#f97316" opacity="0.22"/>'
        )
    for overlap in y_overlaps:
        y = sy(overlap["start"])
        h = overlap["size"] * scale
        lines.append(
            f'<rect x="{margin_x:.2f}" y="{y:.2f}" width="{image_w:.2f}" height="{h:.2f}" fill="#14b8a6" opacity="0.22"/>'
        )

    colors = ["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#0891b2", "#ea580c", "#4f46e5", "#be123c"]
    for index, tile in enumerate(tiles, start=1):
        color = colors[(index - 1) % len(colors)]
        x = sx(tile.x0)
        y = sy(tile.y0)
        w = (tile.x1 - tile.x0) * scale
        h = (tile.y1 - tile.y0) * scale
        lines.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'fill="{color}" opacity="0.12" stroke="#111827" stroke-width="6"/>'
        )
        lines.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'fill="none" stroke="{color}" stroke-width="3"/>'
        )
        lines.append(
            f'<text x="{x + w / 2:.2f}" y="{y + h / 2:.2f}" text-anchor="middle" class="label" fill="{color}">tile {index}</text>'
        )

    legend_y = image_y + image_h + 34
    x_overlap_text = ", ".join(f'{item["start"]}-{item["end"]} ({item["size"]} px)' for item in x_overlaps) or "none"
    y_overlap_text = ", ".join(f'{item["start"]}-{item["end"]} ({item["size"]} px)' for item in y_overlaps) or "none"
    lines.extend(
        [
            f'<rect x="48" y="{legend_y - 16:.2f}" width="18" height="18" fill="#f97316" opacity="0.35"/>',
            f'<text x="74" y="{legend_y:.2f}" class="small">X overlap bands: {html.escape(x_overlap_text)}</text>',
            f'<rect x="48" y="{legend_y + 16:.2f}" width="18" height="18" fill="#14b8a6" opacity="0.35"/>',
            f'<text x="74" y="{legend_y + 32:.2f}" class="small">Y overlap bands: {html.escape(y_overlap_text)}</text>',
        ]
    )

    coord_y = legend_y + 64
    for index, tile in enumerate(tiles[:12], start=1):
        lines.append(
            f'<text x="48" y="{coord_y + (index - 1) * 18:.2f}" class="mono">'
            f'tile {index:02d}: x={tile.x0}-{tile.x1}, y={tile.y0}-{tile.y1}, z={tile.z0}-{tile.z1}'
            "</text>"
        )
    if len(tiles) > 12:
        lines.append(f'<text x="48" y="{coord_y + 12 * 18:.2f}" class="mono">... remaining tile coordinates are in the JSON manifest</text>')

    lines.append("</svg>")
    svg_path.write_text("\n".join(lines), encoding="utf-8")


def running_headless() -> bool:
    if sys.platform.startswith("win") or sys.platform == "darwin":
        return False
    return not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def write_tile_map_only(args, out_dir: Path) -> dict:
    spec = fixed_source_spec(args) if args.fixed_ims else source_spec(args, Path(args.ims).expanduser().resolve(), args.ref_ch)
    tile_size, overlap = tile_schedule(args, spec.shape_zyx)
    tiles = generate_tiles(spec.shape_zyx, tile_size, overlap)
    if args.max_tiles > 0:
        tiles = tiles[: args.max_tiles]
    x_starts = sorted({tile.x0 for tile in tiles})
    y_starts = sorted({tile.y0 for tile in tiles})
    z_starts = sorted({tile.z0 for tile in tiles})
    x_segments = tile_axis_segments(spec.shape_zyx[2], tile_size[2], x_starts)
    y_segments = tile_axis_segments(spec.shape_zyx[1], tile_size[1], y_starts)
    z_segments = tile_axis_segments(spec.shape_zyx[0], tile_size[0], z_starts)
    x_overlaps = adjacent_overlaps(x_segments)
    y_overlaps = adjacent_overlaps(y_segments)
    z_overlaps = adjacent_overlaps(z_segments)

    map_dir = out_dir / "tile_map"
    map_dir.mkdir(parents=True, exist_ok=True)
    svg_path = map_dir / "tile_map_xy.svg"
    json_path = map_dir / "tile_map_xy.json"
    write_tile_map_svg(svg_path, spec, tile_size, overlap, tiles, x_overlaps, y_overlaps)
    manifest = {
        "source": str(spec.path),
        "source_kind": spec.kind,
        "channel": spec.channel,
        "shape_zyx": spec.shape_zyx,
        "spacing_xyz": spec.spacing_xyz,
        "tile_size_zyx": tile_size,
        "overlap_zyx": overlap,
        "tile_count": len(tiles),
        "full_z_tiles": args.full_z,
        "x_segments": x_segments,
        "y_segments": y_segments,
        "z_segments": z_segments,
        "x_overlaps": x_overlaps,
        "y_overlaps": y_overlaps,
        "z_overlaps": z_overlaps,
        "tiles": [tile_to_dict(tile) for tile in tiles],
        "svg": str(svg_path),
    }
    json_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[tile-map] source shape_zyx={spec.shape_zyx}", flush=True)
    print(f"[tile-map] tile_size_zyx={tile_size} overlap_zyx={overlap} tile_count={len(tiles)}", flush=True)
    print(f"[tile-map] wrote SVG: {svg_path}", flush=True)
    print(f"[tile-map] wrote JSON: {json_path}", flush=True)
    should_open = not running_headless()
    if should_open:
        opened = webbrowser.open(svg_path.as_uri())
        if not opened:
            print(f"[tile-map] browser did not report success; open manually: {svg_path}", flush=True)
    elif running_headless():
        print(f"[tile-map] headless session detected; not opening browser. Open manually: {svg_path}", flush=True)
    return manifest


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
            raise ValueError(
                f"Cannot read multi-dimensional TIFF with non-singleton axis {axis}: "
                f"shape={arr.shape} axes={axes}"
            )
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Register an .ims channel to a reference channel using a global transform on whole-volume "
            "proxies followed by tiled QuickSyN on 4x tile proxies, then apply tile transforms back "
            "to full-resolution tiles and overlap-blend them into a BigTIFF."
        )
    )
    parser.add_argument("--ims", "--input-ims", required=True, help="Moving/source multi-channel .ims file")
    parser.add_argument(
        "--fixed-ims",
        default="",
        help=(
            "Optional fixed/reference .ims file. When supplied, the script runs cycle-to-reference "
            "streaming mode: estimate tiled transforms on DAPI, then apply them to --ch channels."
        ),
    )
    parser.add_argument("--ref-ch", "--reference-channel", type=int, default=0, help="Reference channel. Default: 0")
    parser.add_argument("--reg-ch", type=int, default=0, help="Moving DAPI channel used to estimate cycle transforms in --fixed-ims mode. Default: 0")
    parser.add_argument(
        "--ch",
        "--channels",
        dest="channels",
        required=True,
        type=parse_channels,
        help="Moving channel(s) to output. In --fixed-ims mode, transforms are estimated once from --reg-ch and applied to these channels.",
    )
    parser.add_argument("--out", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--spacing", type=parse_spacing_override, default=None, help="Override spacing, for example 0.711")
    parser.add_argument("--ims-resolution-level", type=int, default=0, help=".ims resolution level. Default: 0")
    parser.add_argument("--ims-timepoint", type=int, default=0, help=".ims timepoint. Default: 0")
    parser.add_argument("--ims-axis-order", choices=["auto", "zyx", "xyz"], default="auto", help="Default: auto")
    parser.add_argument("--ims-dataset-template", default="", help="Optional dataset template with {level},{timepoint},{channel}")
    parser.add_argument("--t", "--threads", dest="threads", type=int, default=32, help="ANTs/ITK CPU threads. Default: 32")
    parser.add_argument("--pi", "--progress-interval", dest="progress_interval", type=float, default=30.0, help="Heartbeat seconds. Default: 30")
    parser.add_argument("--trsaa-ds", "--global-ds", dest="trsaa_ds", type=int, default=16, help="Whole-volume global registration downsample stride. Default: 16")
    parser.add_argument(
        "--global-tx",
        default="TRSAA",
        help="Whole-volume global transform type. Use QuickSyN as shorthand for antsRegistrationSyNQuick[s]. Default: TRSAA.",
    )
    parser.add_argument(
        "--global-only",
        action="store_true",
        help="Estimate global registration on the downsampled whole-volume DAPI proxy, write preview/QC output, and stop before tiled full-resolution application.",
    )
    parser.add_argument(
        "--global-qc",
        "--pre-tile-qc",
        dest="global_qc",
        action="store_true",
        help="During a full tiled run, write/open the downsampled global-registration QC before starting tiles.",
    )
    parser.add_argument(
        "--tile-map-only",
        "--map",
        dest="tile_map_only",
        action="store_true",
        help="Write an XY tile/overlap map from image metadata and exit before loading volumes or running registration.",
    )
    parser.add_argument(
        "--qc",
        action="store_true",
        help="After stitched DAPI output, write downsampled before/after TIFF previews, center-slice QC PNG, and similarity metrics.",
    )
    parser.add_argument("--qc-ds", type=int, default=16, help="Downsample stride for stitched DAPI QC output. Default: 16")
    parser.add_argument("--open-qc", action="store_true", help="Write and open the center-slice QC PNG when available.")
    parser.add_argument("--qc-max-panel-side", type=int, default=900, help="Maximum side of each QC panel. Default: 900")
    parser.add_argument("--tile-ds", type=float, default=4.0, help="Tile QuickSyN downsample factor. Default: 4")
    parser.add_argument(
        "--tile",
        "--tile-size",
        default="96,2048,2048",
        help="Full-resolution tile size z,y,x. Default: 96,2048,2048",
    )
    parser.add_argument(
        "--tiles",
        type=int,
        default=0,
        help=(
            "Automatic full-Z XY tiling: specify total tile count and compute the grid, tile size, "
            "and overlap from image metadata. Example: --tiles 2 --overlap 0.25."
        ),
    )
    parser.add_argument(
        "--overlap",
        default="24,512,512",
        help=(
            "Manual mode: full-resolution tile overlap z,y,x. With --tiles, use a fraction like 0.25 "
            "or pixel overlap z,y,x. Default: 24,512,512."
        ),
    )
    parser.add_argument(
        "--full-z",
        action="store_true",
        help="Make every tile span the full Z depth and disable Z overlap; tile only across XY.",
    )
    parser.add_argument("--mxy", "--margin-xy", dest="margin_xy", type=int, default=512, help="Moving read margin XY for global apply. Default: 512")
    parser.add_argument("--mz", "--margin-z", dest="margin_z", type=int, default=16, help="Moving read margin Z for global apply. Default: 16")
    parser.add_argument(
        "--syn-tx",
        default="SyNOnly",
        help=(
            "ANTsPy local transform preset. Default: SyNOnly, because the global transform has already "
            "placed each tile and local rigid/affine initialization can lose overlap."
        ),
    )
    parser.add_argument(
        "--no-local-refine",
        action="store_true",
        help=(
            "Skip per-tile local QuickSyN/SyN and apply only the whole-volume global transform through "
            "the tiled writer. Useful for diagnosing tile seam artifacts."
        ),
    )
    parser.add_argument("--tmp-dir", default="", help="Scratch directory. Default: <out>/scratch")
    parser.add_argument("--keep-temp", action="store_true", help="Keep tile transform scratch files")
    parser.add_argument(
        "--ao",
        "--apply-only",
        dest="apply_only",
        action="store_true",
        help="Cycle-to-reference mode only: reuse existing global and tile transforms, then stream-apply to --ch channels.",
    )
    parser.add_argument(
        "--moving-source-map",
        "--source-map",
        default="",
        help=(
            "Cycle-to-reference mode only: optional semicolon-separated channel=path overrides for moving channels. "
            "Each override can be .ims, .tif, or .tiff."
        ),
    )
    parser.add_argument(
        "--blend",
        choices=["memmap", "ram"],
        default="memmap",
        help="Overlap blending buffer. memmap uses scratch disk; ram is faster but needs about 8 bytes per voxel. Default: memmap.",
    )
    parser.add_argument(
        "--ram-limit-gb",
        type=float,
        default=1500.0,
        help="Safety limit for --blend ram, in GiB. Default: 1500.",
    )
    parser.add_argument(
        "--skip-bg-tiles",
        "--skip-background-tiles",
        action="store_true",
        help=(
            "Cycle-to-reference mode only: skip tiles with little fixed-reference DAPI signal entirely, "
            "leaving those output regions zero for every channel. Off by default."
        ),
    )
    parser.add_argument(
        "--bg-threshold-fraction",
        type=float,
        default=0.02,
        help="Background threshold=min+fraction*(max-min) estimated from reference DAPI proxy. Default: 0.02.",
    )
    parser.add_argument(
        "--min-foreground-fraction",
        type=float,
        default=0.01,
        help="Skip a tile when its fixed DAPI foreground fraction is below this value. Default: 0.01.",
    )
    parser.add_argument(
        "--bg-check-stride",
        type=int,
        default=4,
        help="Subsampling stride used when measuring tile foreground. Default: 4.",
    )
    parser.add_argument("--max-tiles", type=int, default=0, help="Debug limit for number of tiles. Default: all")
    parser.add_argument(
        "--global-transform-dir",
        default="",
        help="Directory containing chN_to_ch0_wholeims_fwdtransforms.json. Default: <out>/transforms",
    )
    parser.add_argument("--force-global", action="store_true", help="Recompute the global transform even if a matching manifest exists")
    parser.add_argument(
        "--uint16-scaling",
        choices=["clip"],
        default="clip",
        help="Output conversion. Currently only clip is supported for streamed tiled output.",
    )
    parser.add_argument(
        "--tiff-prefix",
        "--tif-prefix",
        default="registered_fullres",
        help=(
            "Cycle-to-reference mode only: common final TIFF prefix before Imaris channel marker _C###. "
            "Default: registered_fullres."
        ),
    )
    parser.add_argument(
        "--final-output-dir",
        "--final-out",
        default="",
        help=(
            "Cycle-to-reference mode only: directory for final registered TIFF stacks. "
            "Default: --out. Use one shared folder with --channel-offset to assemble multiple cycles for Imaris."
        ),
    )
    parser.add_argument(
        "--channel-offset",
        "--co",
        type=int,
        default=0,
        help=(
            "Cycle-to-reference mode only: first final _C### channel number for this cycle. "
            "Selected --ch channels are packed consecutively from this value, so the next cycle's offset "
            "is this offset plus the number of exported channels. Default: 0."
        ),
    )
    return parser


def channel_spec(args, channel: int):
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
        Path(args.ims).expanduser().resolve(),
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
    return spec


def source_spec(args, path: Path, channel: int) -> SourceSpec:
    path = path.expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Source does not exist for channel {channel}: {path}")
    suffix = path.suffix.lower()
    if suffix in TIFF_SUFFIXES:
        shape_zyx = tiff_stack_shape_zyx(path)
        spacing_xyz = args.spacing if args.spacing is not None else (1.0, 1.0, 1.0)
        with tifffile.TiffFile(str(path)) as tif:
            series = tif.series[0]
            axes = str(getattr(series, "axes", ""))
            dataset_key = f"TIFF series 0 axes={axes or 'unknown'} pages={len(tif.pages)}"
        return SourceSpec(
            path=path,
            kind="tiff",
            channel=channel,
            dataset_key=dataset_key,
            axis_order="zyx",
            shape_zyx=shape_zyx,
            spacing_xyz=spacing_xyz,
        )
    if suffix != ".ims":
        raise ValueError(f"Unsupported source for channel {channel}: {path}. Use .ims, .tif, or .tiff.")
    ims_spec = build_ims_source_spec(args, path, channel)
    return SourceSpec(
        path=path,
        kind="ims",
        channel=channel,
        dataset_key=ims_spec.dataset_key,
        axis_order=ims_spec.axis_order,
        shape_zyx=ims_spec.shape_zyx,
        spacing_xyz=ims_spec.spacing_xyz,
        ims_spec=ims_spec,
    )


def build_ims_source_spec(args, path: Path, channel: int):
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
        path,
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
    return spec


def moving_source_spec(args, channel: int, source_map: dict[int, Path] | None = None) -> SourceSpec:
    source_map = source_map or {}
    path = source_map.get(channel, Path(args.ims).expanduser().resolve())
    return source_spec(args, path, channel)


def fixed_source_spec(args) -> SourceSpec:
    return source_spec(args, Path(args.fixed_ims).expanduser().resolve(), args.ref_ch)


def read_tiff_tile_from_source(spec: SourceSpec, tile: Tile) -> sitk.Image:
    z_size, y_size, x_size = spec.shape_zyx
    if tile.x0 < 0 or tile.y0 < 0 or tile.z0 < 0 or tile.x1 > x_size or tile.y1 > y_size or tile.z1 > z_size:
        raise ValueError(
            f"Tile {tile} extends outside TIFF source bounds xyz-size=({x_size}, {y_size}, {z_size}): {spec.path}"
        )
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


def read_tiff_downsampled_from_source(spec: SourceSpec, stride_zyx: tuple[int, int, int]) -> sitk.Image:
    sz = max(1, int(stride_zyx[0]))
    sy = max(1, int(stride_zyx[1]))
    sx = max(1, int(stride_zyx[2]))
    with tifffile.TiffFile(str(spec.path)) as tif:
        if len(tif.pages) > 1 and len(tif.pages[0].shape) == 2:
            planes = [np.asarray(tif.pages[z].asarray())[::sy, ::sx] for z in range(0, spec.shape_zyx[0], sz)]
            arr = np.stack(planes, axis=0)
        else:
            series = tif.series[0]
            arr_zyx = tiff_series_array_to_zyx(np.asarray(series.asarray()), str(getattr(series, "axes", "")))
            arr = arr_zyx[::sz, ::sy, ::sx]
    img = sitk.GetImageFromArray(np.ascontiguousarray(arr))
    img.SetSpacing((spec.spacing_xyz[0] * sx, spec.spacing_xyz[1] * sy, spec.spacing_xyz[2] * sz))
    return img


def read_source_tile(spec, tile: Tile) -> sitk.Image:
    if isinstance(spec, SourceSpec):
        if spec.kind == "ims":
            if spec.ims_spec is None:
                raise ValueError(f"Missing .ims metadata for source: {spec.path}")
            return read_ims_tile_from_spec(spec.ims_spec, tile)
        if spec.kind == "tiff":
            return read_tiff_tile_from_source(spec, tile)
        raise ValueError(f"Unsupported source kind: {spec.kind}")
    return read_ims_tile_from_spec(spec, tile)


def read_source_downsampled(spec, stride_zyx: tuple[int, int, int]) -> sitk.Image:
    if isinstance(spec, SourceSpec):
        if spec.kind == "ims":
            if spec.ims_spec is None:
                raise ValueError(f"Missing .ims metadata for source: {spec.path}")
            return read_ims_downsampled_from_spec(spec.ims_spec, stride_zyx)
        if spec.kind == "tiff":
            return read_tiff_downsampled_from_source(spec, stride_zyx)
        raise ValueError(f"Unsupported source kind: {spec.kind}")
    return read_ims_downsampled_from_spec(spec, stride_zyx)


def sitk_to_ants(img: sitk.Image):
    import ants

    return ants.from_numpy(
        sitk.GetArrayFromImage(img).T,
        origin=img.GetOrigin(),
        spacing=img.GetSpacing(),
        direction=np.asarray(img.GetDirection()).reshape((3, 3)),
    )


def ants_to_zyx_float(img) -> np.ndarray:
    return img.numpy().T.astype(np.float32, copy=False)


def ants_to_uint16_zyx(img) -> np.ndarray:
    arr = img.numpy().T
    return np.rint(np.clip(arr, 0, 65535)).astype(np.uint16, copy=False)


def downsample_ants(img, factor: float, label: str, progress_interval: float):
    import ants

    if factor <= 1:
        return img
    spacing = tuple(float(v) * factor for v in img.spacing)
    print(f"[tiled-syn] downsampling {label} by spacing factor {factor:g}", flush=True)
    with progress_heartbeat(f"downsample {label}", progress_interval):
        return ants.resample_image(img, spacing, use_voxels=False, interp_type=0)


def domain_chunk_ants(spec, tile: Tile):
    arr = np.zeros(tile.shape, dtype=np.float32)
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spec.spacing_xyz)
    img.SetOrigin(
        (
            tile.x0 * spec.spacing_xyz[0],
            tile.y0 * spec.spacing_xyz[1],
            tile.z0 * spec.spacing_xyz[2],
        )
    )
    return sitk_to_ants(img)


def padded_tile(tile: Tile, spec, margin_xy: int, margin_z: int) -> Tile:
    z_size, y_size, x_size = spec.shape_zyx
    return Tile(
        z0=max(0, tile.z0 - margin_z),
        z1=min(z_size, tile.z1 + margin_z),
        y0=max(0, tile.y0 - margin_xy),
        y1=min(y_size, tile.y1 + margin_xy),
        x0=max(0, tile.x0 - margin_xy),
        x1=min(x_size, tile.x1 + margin_xy),
    )


def global_manifest_path(tx_dir: Path, channel: int, ref_ch: int) -> Path:
    return tx_dir / f"ch{channel}_to_ch{ref_ch}_wholeims_fwdtransforms.json"


def resolved_global_transform(value: str) -> str:
    if value.strip().lower() in {"quicksyn", "quick-syn"}:
        return "antsRegistrationSyNQuick[s]"
    return value


def global_transform_label(value: str) -> str:
    return "QuickSyN" if resolved_global_transform(value) == "antsRegistrationSyNQuick[s]" else value


def load_or_estimate_global_transform(args, ref_spec, moving_spec, channel: int, tx_dir: Path) -> list[str]:
    import ants

    tx_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = global_manifest_path(tx_dir, channel, args.ref_ch)
    stride = max(1, int(args.trsaa_ds))
    transform_type = resolved_global_transform(args.global_tx)
    transform_label = global_transform_label(args.global_tx)
    if manifest_path.exists() and not args.force_global:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        transforms = payload.get("fwdtransforms")
        if transforms and payload.get("transform", "TRSAA") == transform_type and int(payload.get("downsample_factor", stride)) == stride:
            print(f"[tiled-syn] channel {channel}: reusing global {transform_label} manifest {manifest_path}", flush=True)
            return list(transforms)
        print(
            f"[tiled-syn] channel {channel}: existing global manifest does not match "
            f"requested transform={transform_label} ds={stride}; recomputing",
            flush=True,
        )
    if args.apply_only:
        raise FileNotFoundError(f"--apply-only requested but matching global transform manifest is missing: {manifest_path}")

    stride_zyx = (stride, stride, stride)
    print(f"[tiled-syn] channel {channel}: estimating global {transform_label} at stride_zyx={stride_zyx}", flush=True)
    with progress_heartbeat("load reference global proxy", args.progress_interval):
        fixed_sitk = read_source_downsampled(ref_spec, stride_zyx)
    with progress_heartbeat(f"load channel {channel} global proxy", args.progress_interval):
        moving_sitk = read_source_downsampled(moving_spec, stride_zyx)
    fixed = sitk_to_ants(fixed_sitk)
    moving = sitk_to_ants(moving_sitk)
    del fixed_sitk, moving_sitk
    gc.collect()
    outprefix = tx_dir / f"ch{channel}_to_ch{args.ref_ch}_wholeims_"
    with progress_heartbeat(f"global {transform_label} channel {channel}", args.progress_interval):
        tx = ants.registration(
            fixed=fixed,
            moving=moving,
            type_of_transform=transform_type,
            outprefix=str(outprefix),
            aff_metric="mattes",
            syn_metric="mattes",
            singleprecision=True,
            verbose=True,
        )
    payload = {
        "input_ims": str(Path(args.ims).expanduser().resolve()),
        "reference_channel": args.ref_ch,
        "moving_channel": channel,
        "transform": transform_type,
        "downsample_factor": stride,
        "stride_zyx": stride_zyx,
        "fwdtransforms": tx["fwdtransforms"],
    }
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    del fixed, moving
    gc.collect()
    print(f"[tiled-syn] channel {channel}: wrote global {transform_label} manifest {manifest_path}", flush=True)
    return list(tx["fwdtransforms"])


def write_preview_tiff(img, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    arr = ants_to_uint16_zyx(img)
    tifffile.imwrite(str(output_path), arr, photometric="minisblack", bigtiff=True, metadata=None)
    print(f"[global-preview] wrote downsampled TIFF: {output_path}", flush=True)
    del arr


def write_global_preview_qc(
    args,
    fixed_spec,
    moving_spec,
    transforms: list[str],
    out_dir: Path,
    timings: list[dict] | None = None,
) -> dict:
    import ants

    stride = max(1, int(args.trsaa_ds))
    stride_zyx = (stride, stride, stride)
    transform_label = global_transform_label(args.global_tx)

    stage_start = time.time()
    with progress_heartbeat("load fixed global preview proxy", args.progress_interval):
        fixed = sitk_to_ants(read_source_downsampled(fixed_spec, stride_zyx))
    with progress_heartbeat("load moving global preview proxy", args.progress_interval):
        moving = sitk_to_ants(read_source_downsampled(moving_spec, stride_zyx))
    with progress_heartbeat("place unregistered moving proxy in fixed grid for QC", args.progress_interval):
        moving_before = ants.resample_image_to_target(moving, fixed, interp_type="linear")
    with progress_heartbeat("apply global preview transform", args.progress_interval):
        warped = ants.apply_transforms(
            fixed=fixed,
            moving=moving,
            transformlist=transforms,
            interpolator="linear",
            singleprecision=True,
        )
    if timings is not None:
        record_timing(timings, "apply transform to downsampled global preview", stage_start, channel=args.reg_ch)

    preview_dir = out_dir / "global_preview"
    qc_dir = preview_dir / "qc"
    prefix = f"ch{args.reg_ch}_to_fixed_ch{args.ref_ch}_global_{transform_label}_ds{stride}"
    fixed_path = preview_dir / f"{prefix}_fixed.tif"
    moving_path = preview_dir / f"{prefix}_before.tif"
    warped_path = preview_dir / f"{prefix}_registered.tif"
    stage_start = time.time()
    write_preview_tiff(fixed, fixed_path)
    write_preview_tiff(moving_before, moving_path)
    write_preview_tiff(warped, warped_path)
    export_center_slice_qc(
        fixed,
        moving_before,
        warped,
        qc_dir,
        prefix,
        args.open_qc,
        args.qc_max_panel_side,
    )
    if timings is not None:
        record_timing(timings, "write global preview TIFFs and QC", stage_start, channel=args.reg_ch)

    metrics = {
        "computed_on": f"whole_volume_proxy_ds{stride}",
        "transform": resolved_global_transform(args.global_tx),
        "interpretation": "higher ncc/nmi is better; lower nrmse is better",
        "before": compute_similarity_metrics_from_images(fixed, moving_before),
        "after": compute_similarity_metrics_from_images(fixed, warped),
    }
    metrics_path = preview_dir / f"{prefix}_similarity_metrics.json"
    write_metrics_json(metrics_path, metrics)
    result = {
        "mode": "global_preview",
        "fixed_ims": str(Path(args.fixed_ims).expanduser().resolve()),
        "moving_ims": str(Path(args.ims).expanduser().resolve()),
        "fixed_channel": args.ref_ch,
        "moving_registration_channel": args.reg_ch,
        "transform": resolved_global_transform(args.global_tx),
        "downsample_factor": stride,
        "stride_zyx": stride_zyx,
        "fwdtransforms": transforms,
        "fixed_preview_tif": str(fixed_path),
        "moving_before_preview_tif": str(moving_path),
        "registered_preview_tif": str(warped_path),
        "qc_png": str(qc_dir / f"{prefix}_center_slices_overlay.png"),
        "metrics_json": str(metrics_path),
        "timings": timings,
    }
    del fixed, moving, moving_before, warped
    gc.collect()
    return result


def run_global_preview(args, out_dir: Path, global_tx_dir: Path) -> dict:
    if not args.fixed_ims:
        raise ValueError("--global-only requires --fixed-ims for cycle-to-reference registration")
    if args.apply_only:
        raise ValueError("--global-only cannot be combined with --apply-only")
    if args.reg_ch not in args.channels:
        raise ValueError("--reg-ch must be included in --ch for --global-only output")

    fixed_spec = fixed_source_spec(args)
    source_map = parse_source_map(args.moving_source_map) if args.moving_source_map else {}
    moving_spec = moving_source_spec(args, args.reg_ch, source_map)
    transform_label = global_transform_label(args.global_tx)
    timings: list[dict] = []
    stage_start = time.time()
    transforms = load_or_estimate_global_transform(args, fixed_spec, moving_spec, args.reg_ch, global_tx_dir)
    record_timing(timings, f"global {transform_label} estimate/reuse", stage_start, channel=args.reg_ch)
    return write_global_preview_qc(args, fixed_spec, moving_spec, transforms, out_dir, timings)


def write_tiled_dapi_qc(
    args,
    fixed_spec,
    moving_spec,
    registered_output_path: Path,
    out_dir: Path,
    channel: int,
) -> dict:
    import ants

    stride = max(1, int(args.qc_ds))
    stride_zyx = (stride, stride, stride)
    registered_spec = source_spec(args, registered_output_path, channel)
    qc_dir = out_dir / "qc"
    prefix = f"ch{channel}_to_fixed_ch{args.ref_ch}_tiled_registered_ds{stride}"
    print(f"[qc] reading stitched DAPI output at stride_zyx={stride_zyx}", flush=True)
    with progress_heartbeat("load fixed DAPI QC proxy", args.progress_interval):
        fixed = sitk_to_ants(read_source_downsampled(fixed_spec, stride_zyx))
    with progress_heartbeat("load moving DAPI before-QC proxy", args.progress_interval):
        moving = sitk_to_ants(read_source_downsampled(moving_spec, stride_zyx))
    with progress_heartbeat("place unregistered moving DAPI in fixed QC grid", args.progress_interval):
        moving_before = ants.resample_image_to_target(moving, fixed, interp_type="linear")
    with progress_heartbeat("load registered stitched DAPI QC proxy", args.progress_interval):
        registered = sitk_to_ants(read_source_downsampled(registered_spec, stride_zyx))

    fixed_path = qc_dir / f"{prefix}_fixed.tif"
    before_path = qc_dir / f"{prefix}_before.tif"
    registered_path = qc_dir / f"{prefix}_registered.tif"
    write_preview_tiff(fixed, fixed_path)
    write_preview_tiff(moving_before, before_path)
    write_preview_tiff(registered, registered_path)
    export_center_slice_qc(
        fixed,
        moving_before,
        registered,
        qc_dir,
        prefix,
        args.open_qc,
        args.qc_max_panel_side,
    )
    metrics = {
        "computed_on": f"stitched_fullres_output_proxy_ds{stride}",
        "registered_output_tif": str(registered_output_path),
        "interpretation": "higher ncc/nmi is better; lower nrmse is better",
        "before": compute_similarity_metrics_from_images(fixed, moving_before),
        "after": compute_similarity_metrics_from_images(fixed, registered),
    }
    metrics_path = qc_dir / f"{prefix}_similarity_metrics.json"
    write_metrics_json(metrics_path, metrics)
    result = {
        "downsample_factor": stride,
        "fixed_preview_tif": str(fixed_path),
        "moving_before_preview_tif": str(before_path),
        "registered_preview_tif": str(registered_path),
        "overlay_png": str(qc_dir / f"{prefix}_center_slices_overlay.png"),
        "metrics_json": str(metrics_path),
    }
    del fixed, moving, moving_before, registered
    gc.collect()
    return result


def globally_align_tile(ref_spec, moving_spec, tile: Tile, global_transforms: list[str], margin_xy: int, margin_z: int, progress_interval: float):
    import ants

    fixed_domain = domain_chunk_ants(ref_spec, tile)
    moving_tile = padded_tile(tile, moving_spec, margin_xy, margin_z)
    with progress_heartbeat("load moving padded tile for global apply", progress_interval):
        moving_sitk = read_source_tile(moving_spec, moving_tile)
    moving = sitk_to_ants(moving_sitk)
    del moving_sitk
    with progress_heartbeat("apply global transform to tile", progress_interval):
        aligned = ants.apply_transforms(
            fixed=fixed_domain,
            moving=moving,
            transformlist=global_transforms,
            interpolator="linear",
            singleprecision=True,
        )
    del moving
    return fixed_domain, aligned


def create_blend_buffers(
    shape_zyx: tuple[int, int, int],
    channel_scratch: Path,
    blend_mode: str,
    ram_limit_gb: float,
    channel_label: str,
) -> tuple[np.ndarray, np.ndarray]:
    required_gib = float(np.prod(shape_zyx) * 8 / 1024**3)
    print(
        f"[tiled-syn] {channel_label}: blend buffers need about {required_gib:.1f} GiB "
        f"for accum+weight using {blend_mode}",
        flush=True,
    )
    if blend_mode == "ram":
        if required_gib > ram_limit_gb:
            raise MemoryError(
                f"--blend ram would need about {required_gib:.1f} GiB, above --ram-limit-gb {ram_limit_gb:.1f}. "
                "Use --blend memmap or raise the limit."
            )
        return (
            np.zeros(shape_zyx, dtype=np.float32),
            np.zeros(shape_zyx, dtype=np.float32),
        )

    channel_scratch.mkdir(parents=True, exist_ok=True)
    acc_path = channel_scratch / "accum.float32.dat"
    weight_path = channel_scratch / "weight.float32.dat"
    accum = np.memmap(acc_path, mode="w+", dtype=np.float32, shape=shape_zyx)
    weight = np.memmap(weight_path, mode="w+", dtype=np.float32, shape=shape_zyx)
    accum[:] = 0.0
    weight[:] = 0.0
    return accum, weight


def flush_blend_buffers(accum: np.ndarray, weight: np.ndarray) -> None:
    if isinstance(accum, np.memmap):
        accum.flush()
    if isinstance(weight, np.memmap):
        weight.flush()


def record_timing(timings: list[dict], stage: str, start: float, **extra) -> float:
    elapsed = time.time() - start
    entry = {"stage": stage, "seconds": elapsed, "duration": format_duration(elapsed)}
    entry.update(extra)
    timings.append(entry)
    print(f"[timing] {stage}: {format_duration(elapsed)}", flush=True)
    return elapsed


def write_output_tiff_from_buffers(accum: np.ndarray, weight: np.ndarray, output_path: Path, progress_interval: float) -> None:
    z_size = accum.shape[0]
    start = time.time()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[tiled-syn] writing stitched BigTIFF: {output_path}", flush=True)
    with tifffile.TiffWriter(output_path, bigtiff=True) as writer:
        for z in range(z_size):
            plane = accum[z] / np.maximum(weight[z], np.float32(1e-6))
            plane_u16 = np.rint(np.clip(plane, 0, 65535)).astype(np.uint16, copy=False)
            writer.write(plane_u16, photometric="minisblack", metadata=None, contiguous=True)
            if progress_interval > 0 and (time.time() - start) >= progress_interval:
                # Reset start only for heartbeat style output.
                print(f"[tiled-syn] writing z={z + 1}/{z_size}", flush=True)
                start = time.time()
    print(f"[tiled-syn] wrote stitched BigTIFF: {output_path}", flush=True)


def run_channel(args, ref_spec, channel: int, out_dir: Path, scratch_root: Path, global_tx_dir: Path) -> dict:
    import ants

    moving_spec = channel_spec(args, channel)
    if moving_spec.shape_zyx != ref_spec.shape_zyx:
        raise ValueError(f"Channel {channel} shape mismatch: ref={ref_spec.shape_zyx}, moving={moving_spec.shape_zyx}")

    global_transforms = load_or_estimate_global_transform(args, ref_spec, moving_spec, channel, global_tx_dir)

    tile_size, overlap = tile_schedule(args, ref_spec.shape_zyx)
    tiles = generate_tiles(ref_spec.shape_zyx, tile_size, overlap)
    if args.max_tiles > 0:
        tiles = tiles[: args.max_tiles]
    print(
        f"[tiled-syn] channel {channel}: {len(tiles)} tiles shape_zyx={ref_spec.shape_zyx} "
        f"tile={tile_size} overlap={overlap} tile_ds={args.tile_ds}",
        flush=True,
    )

    channel_scratch = scratch_root / f"ch{channel}"
    accum, weight = create_blend_buffers(
        ref_spec.shape_zyx,
        channel_scratch,
        args.blend,
        args.ram_limit_gb,
        f"channel {channel}",
    )

    weights_cache: dict[tuple[int, int, int], np.ndarray] = {}
    transforms_dir = out_dir / "tile_quicksyn_transforms" / f"ch{channel}"
    transforms_dir.mkdir(parents=True, exist_ok=True)
    start = time.time()

    for idx, tile in enumerate(tiles, start=1):
        tile_start = time.time()
        print(f"[tiled-syn] channel {channel}: tile {idx}/{len(tiles)} {tile}", flush=True)
        with progress_heartbeat(f"load fixed tile {idx}", args.progress_interval):
            fixed_sitk = read_source_tile(ref_spec, tile)
        fixed = sitk_to_ants(fixed_sitk)
        del fixed_sitk

        _, moving_global = globally_align_tile(
            ref_spec,
            moving_spec,
            tile,
            global_transforms,
            args.margin_xy,
            args.margin_z,
            args.progress_interval,
        )

        outprefix = transforms_dir / f"tile_{idx:04d}_"
        fixed_ds = None
        moving_ds = None
        if args.no_local_refine:
            tx = {"fwdtransforms": []}
            refined = moving_global
            print(f"[tiled-syn] channel {channel}: tile {idx}: --no-local-refine, using global transform only", flush=True)
        else:
            fixed_ds = downsample_ants(fixed, args.tile_ds, f"fixed tile {idx}", args.progress_interval)
            moving_ds = downsample_ants(moving_global, args.tile_ds, f"moving tile {idx}", args.progress_interval)
            with progress_heartbeat(f"QuickSyN tile {idx}", args.progress_interval):
                tx = ants.registration(
                    fixed=fixed_ds,
                    moving=moving_ds,
                    type_of_transform=args.syn_tx,
                    outprefix=str(outprefix),
                    aff_metric="mattes",
                    syn_metric="mattes",
                    singleprecision=True,
                    verbose=True,
                )
            with progress_heartbeat(f"apply QuickSyN tile {idx}", args.progress_interval):
                refined = ants.apply_transforms(
                    fixed=fixed,
                    moving=moving_global,
                    transformlist=tx["fwdtransforms"],
                    interpolator="linear",
                    singleprecision=True,
                )

        arr = ants_to_zyx_float(refined)
        if arr.shape != tile.shape:
            raise ValueError(f"Tile {idx} output shape mismatch: expected {tile.shape}, got {arr.shape}")
        tile_weights = weights_cache.get(tile.shape)
        if tile_weights is None:
            tile_weights = raised_cosine_weights(tile.shape, overlap)
            weights_cache[tile.shape] = tile_weights
        accum[tile.z0 : tile.z1, tile.y0 : tile.y1, tile.x0 : tile.x1] += arr * tile_weights
        weight[tile.z0 : tile.z1, tile.y0 : tile.y1, tile.x0 : tile.x1] += tile_weights
        flush_blend_buffers(accum, weight)

        if not args.keep_temp:
            for path in transforms_dir.glob(f"tile_{idx:04d}_*"):
                try:
                    path.unlink()
                except OSError:
                    pass

        del fixed, moving_global, refined, arr
        if fixed_ds is not None:
            del fixed_ds
        if moving_ds is not None:
            del moving_ds
        gc.collect()
        elapsed = time.time() - start
        eta = (elapsed / idx) * (len(tiles) - idx)
        print(
            f"[tiled-syn] channel {channel}: tile {idx}/{len(tiles)} done in "
            f"{format_duration(time.time() - tile_start)}; elapsed={format_duration(elapsed)} eta={format_duration(eta)}",
            flush=True,
        )

    output_path = out_dir / f"ch{channel}_to_ch{args.ref_ch}_wholeims_trsaa_tiled_quicksyn_fullres.tif"
    write_output_tiff_from_buffers(accum, weight, output_path, args.progress_interval)
    result = {
        "channel": channel,
        "output_tif": str(output_path),
        "global_transforms": global_transforms,
        "tile_size_zyx": tile_size,
        "overlap_zyx": overlap,
        "tile_downsample_factor": args.tile_ds,
        "no_local_refine": args.no_local_refine,
        "tile_count": len(tiles),
    }
    (out_dir / f"ch{channel}_tiled_quicksyn_manifest.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def tile_transform_manifest_path(transforms_dir: Path, tile_index: int) -> Path:
    return transforms_dir / f"tile_{tile_index:04d}_fwdtransforms.json"


def save_tile_transform_manifest(
    transforms_dir: Path,
    tile_index: int,
    tile: Tile,
    transforms: list[str],
    skipped: bool = False,
    reason: str = "",
    output_skipped: bool = False,
    foreground_fraction: float | None = None,
) -> None:
    path = tile_transform_manifest_path(transforms_dir, tile_index)
    payload = {
        "tile_index": tile_index,
        "tile": tile_to_dict(tile),
        "fwdtransforms": transforms,
        "local_refinement_skipped": skipped,
        "skip_reason": reason,
        "output_skipped": output_skipped,
        "fixed_dapi_foreground_fraction": foreground_fraction,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def load_tile_manifest(transforms_dir: Path, tile_index: int) -> dict:
    path = tile_transform_manifest_path(transforms_dir, tile_index)
    if not path.exists():
        raise FileNotFoundError(f"Missing tile QuickSyN transform manifest: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_tile_transforms(payload: dict, tile_index: int) -> list[str]:
    transforms = payload.get("fwdtransforms", [])
    if payload.get("local_refinement_skipped"):
        print(
            f"[cycle-stream] tile {tile_index}: local refinement was skipped; applying global transform only "
            f"({payload.get('skip_reason', 'no reason recorded')})",
            flush=True,
        )
        return []
    if not transforms:
        raise ValueError(f"Tile transform manifest has no fwdtransforms: {path}")
    return list(transforms)


def foreground_fraction(arr: np.ndarray, threshold: float, stride: int) -> float:
    stride = max(1, int(stride))
    sampled = arr[::stride, ::stride, ::stride]
    return float(np.mean(sampled > threshold))


def estimate_background_threshold(args, fixed_spec) -> float:
    stride = max(1, int(args.trsaa_ds))
    with progress_heartbeat("load reference DAPI background proxy", args.progress_interval):
        proxy = read_source_downsampled(fixed_spec, (stride, stride, stride))
    arr = sitk.GetArrayViewFromImage(proxy)
    finite = np.asarray(arr)[np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError("Fixed DAPI background proxy contains no finite voxels")
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    threshold = vmin + float(args.bg_threshold_fraction) * (vmax - vmin)
    print(
        f"[cycle-stream] background skipping enabled: fixed DAPI proxy range=({vmin:.6g}, {vmax:.6g}) "
        f"threshold={threshold:.6g} min_foreground_fraction={args.min_foreground_fraction:g}",
        flush=True,
    )
    return threshold


def tile_has_foreground(img, label: str, tile_index: int) -> bool:
    arr = img.numpy()
    finite = arr[np.isfinite(arr)]
    nonzero = int(np.count_nonzero(finite)) if finite.size else 0
    minimum = max(100, int(arr.size * 1e-5))
    usable = bool(finite.size and nonzero >= minimum and float(np.max(finite)) > float(np.min(finite)))
    if not usable:
        print(
            f"[cycle-stream] tile {tile_index}: {label} lacks usable foreground "
            f"(nonzero={nonzero}/{arr.size}); local refinement will be skipped",
            flush=True,
        )
    del arr, finite
    return usable


def cycle_output_channel(channels: list[int], channel: int, channel_offset: int) -> int:
    try:
        selected_index = channels.index(channel)
    except ValueError as exc:
        raise ValueError(f"Channel {channel} is not included in selected output channels {channels}") from exc
    output_channel = channel_offset + selected_index
    if output_channel < 0:
        raise ValueError(f"Final TIFF channel number must be non-negative, got offset {channel_offset} + index {selected_index}")
    return output_channel


def cycle_output_path(output_dir: Path, channels: list[int], channel: int, channel_offset: int, tiff_prefix: str) -> Path:
    output_channel = cycle_output_channel(channels, channel, channel_offset)
    return output_dir / f"{tiff_prefix}_C{output_channel:03d}.tif"


def write_registered_channel_from_buffers(
    args,
    out_dir: Path,
    scratch_root: Path,
    channel: int,
    ref_ch: int,
    reg_ch: int,
    accum: np.ndarray,
    weight: np.ndarray,
) -> Path:
    final_output_dir = Path(args.final_output_dir).expanduser().resolve() if args.final_output_dir else out_dir
    output_path = cycle_output_path(final_output_dir, args.channels, channel, args.channel_offset, args.tiff_prefix)
    write_output_tiff_from_buffers(accum, weight, output_path, args.progress_interval)
    flush_blend_buffers(accum, weight)
    del accum, weight
    gc.collect()
    return output_path


def estimate_tiles_and_write_registration_channel(
    args,
    fixed_spec,
    moving_reg_spec,
    tiles: list[Tile],
    overlap: tuple[int, int, int],
    global_transforms: list[str],
    transforms_dir: Path,
    out_dir: Path,
    scratch_root: Path,
    timings: list[dict],
    background_threshold: float | None,
) -> dict:
    import ants

    transforms_dir.mkdir(parents=True, exist_ok=True)
    channel = args.reg_ch
    channel_scratch = scratch_root / f"cycle_reg_ch{channel}"
    accum, weight = create_blend_buffers(
        fixed_spec.shape_zyx,
        channel_scratch,
        args.blend,
        args.ram_limit_gb,
        f"cycle registration channel {channel}",
    )
    weights_cache: dict[tuple[int, int, int], np.ndarray] = {}
    start = time.time()

    for idx, tile in enumerate(tiles, start=1):
        tile_start = time.time()
        print(f"[cycle-stream] DAPI tile {idx}/{len(tiles)} {tile}", flush=True)
        stage_start = time.time()
        with progress_heartbeat(f"load fixed DAPI tile {idx}", args.progress_interval):
            fixed_sitk = read_source_tile(fixed_spec, tile)
        if background_threshold is not None:
            fixed_arr = sitk.GetArrayViewFromImage(fixed_sitk)
            fg_fraction = foreground_fraction(fixed_arr, background_threshold, args.bg_check_stride)
            print(
                f"[cycle-stream] DAPI tile {idx}: fixed foreground fraction={fg_fraction:.6g}",
                flush=True,
            )
            if fg_fraction < args.min_foreground_fraction:
                skip_reason = (
                    f"fixed DAPI foreground fraction {fg_fraction:.6g} below "
                    f"{args.min_foreground_fraction:.6g}"
                )
                save_tile_transform_manifest(
                    transforms_dir,
                    idx,
                    tile,
                    [],
                    skipped=True,
                    reason=skip_reason,
                    output_skipped=True,
                    foreground_fraction=fg_fraction,
                )
                del fixed_sitk
                record_timing(timings, f"DAPI tile {idx} skipped background", stage_start, channel=channel, tile_index=idx)
                elapsed = time.time() - start
                eta = (elapsed / idx) * (len(tiles) - idx)
                print(
                    f"[cycle-stream] DAPI tile {idx}/{len(tiles)} skipped as background; "
                    f"elapsed={format_duration(elapsed)} eta={format_duration(eta)}",
                    flush=True,
                )
                continue
        fixed = sitk_to_ants(fixed_sitk)
        del fixed_sitk
        record_timing(timings, f"DAPI tile {idx} load fixed tile", stage_start, channel=channel, tile_index=idx)

        stage_start = time.time()
        _, moving_global = globally_align_tile(
            fixed_spec,
            moving_reg_spec,
            tile,
            global_transforms,
            args.margin_xy,
            args.margin_z,
            args.progress_interval,
        )
        record_timing(timings, f"DAPI tile {idx} load/apply global transform", stage_start, channel=channel, tile_index=idx)

        outprefix = transforms_dir / f"tile_{idx:04d}_"
        local_transforms: list[str] = []
        skip_reason = ""
        fixed_ds = None
        moving_ds = None
        if args.no_local_refine:
            skip_reason = "--no-local-refine requested"
            print(f"[cycle-stream] DAPI tile {idx}: --no-local-refine, using global transform only", flush=True)
        else:
            stage_start = time.time()
            fixed_ds = downsample_ants(fixed, args.tile_ds, f"fixed DAPI tile {idx}", args.progress_interval)
            moving_ds = downsample_ants(moving_global, args.tile_ds, f"moving DAPI tile {idx}", args.progress_interval)
            record_timing(timings, f"DAPI tile {idx} downsample for QuickSyN", stage_start, channel=channel, tile_index=idx)
            if not tile_has_foreground(fixed_ds, "fixed DAPI", idx) or not tile_has_foreground(moving_ds, "moving DAPI", idx):
                skip_reason = "insufficient foreground in downsampled DAPI tile"
            else:
                stage_start = time.time()
                try:
                    with progress_heartbeat(f"local SyN DAPI tile {idx}", args.progress_interval):
                        tx = ants.registration(
                            fixed=fixed_ds,
                            moving=moving_ds,
                            type_of_transform=args.syn_tx,
                            initial_transform="Identity",
                            outprefix=str(outprefix),
                            aff_metric="mattes",
                            syn_metric="mattes",
                            singleprecision=True,
                            verbose=True,
                        )
                    local_transforms = list(tx["fwdtransforms"])
                except RuntimeError as exc:
                    skip_reason = f"ANTs local refinement failed: {exc}"
                    print(
                        f"[cycle-stream] WARNING DAPI tile {idx}: local refinement failed; "
                        "falling back to global transform only",
                        flush=True,
                    )
                record_timing(timings, f"DAPI tile {idx} local SyN registration", stage_start, channel=channel, tile_index=idx)
        save_tile_transform_manifest(
            transforms_dir,
            idx,
            tile,
            local_transforms,
            skipped=not bool(local_transforms),
            reason=skip_reason,
            output_skipped=False,
        )
        if local_transforms:
            stage_start = time.time()
            with progress_heartbeat(f"apply DAPI tile local SyN {idx}", args.progress_interval):
                refined = ants.apply_transforms(
                    fixed=fixed,
                    moving=moving_global,
                    transformlist=local_transforms,
                    interpolator="linear",
                    singleprecision=True,
                )
            record_timing(timings, f"DAPI tile {idx} apply local SyN", stage_start, channel=channel, tile_index=idx)
        else:
            refined = moving_global
            print(f"[cycle-stream] DAPI tile {idx}: using global transform output without local refinement", flush=True)

        stage_start = time.time()
        arr = ants_to_zyx_float(refined)
        if arr.shape != tile.shape:
            raise ValueError(f"DAPI tile {idx} output shape mismatch: expected {tile.shape}, got {arr.shape}")
        tile_weights = weights_cache.get(tile.shape)
        if tile_weights is None:
            tile_weights = raised_cosine_weights(tile.shape, overlap)
            weights_cache[tile.shape] = tile_weights
        accum[tile.z0 : tile.z1, tile.y0 : tile.y1, tile.x0 : tile.x1] += arr * tile_weights
        weight[tile.z0 : tile.z1, tile.y0 : tile.y1, tile.x0 : tile.x1] += tile_weights
        flush_blend_buffers(accum, weight)
        record_timing(timings, f"DAPI tile {idx} blend", stage_start, channel=channel, tile_index=idx)

        del fixed, moving_global, refined, arr
        if fixed_ds is not None:
            del fixed_ds
        if moving_ds is not None:
            del moving_ds
        gc.collect()
        elapsed = time.time() - start
        eta = (elapsed / idx) * (len(tiles) - idx)
        print(
            f"[cycle-stream] DAPI tile {idx}/{len(tiles)} done in {format_duration(time.time() - tile_start)}; "
            f"elapsed={format_duration(elapsed)} eta={format_duration(eta)}",
            flush=True,
        )
        record_timing(timings, f"DAPI tile {idx} total", tile_start, channel=channel, tile_index=idx)

    stage_start = time.time()
    output_path = write_registered_channel_from_buffers(
        args,
        out_dir,
        scratch_root,
        channel,
        args.ref_ch,
        args.reg_ch,
        accum,
        weight,
    )
    record_timing(timings, f"channel {channel} write registered DAPI TIFF", stage_start, channel=channel)
    result = {
        "channel": channel,
        "role": "registration_channel",
        "output_tif": str(output_path),
        "tile_transform_dir": str(transforms_dir),
    }
    if args.qc:
        stage_start = time.time()
        result["qc"] = write_tiled_dapi_qc(args, fixed_spec, moving_reg_spec, output_path, out_dir, channel)
        record_timing(timings, f"channel {channel} stitched DAPI QC", stage_start, channel=channel)
    return result


def apply_saved_tile_transforms_to_channel(
    args,
    fixed_spec,
    moving_spec,
    channel: int,
    tiles: list[Tile],
    overlap: tuple[int, int, int],
    global_transforms: list[str],
    transforms_dir: Path,
    out_dir: Path,
    scratch_root: Path,
    timings: list[dict],
) -> dict:
    import ants

    channel_scratch = scratch_root / f"cycle_ch{channel}"
    accum, weight = create_blend_buffers(
        fixed_spec.shape_zyx,
        channel_scratch,
        args.blend,
        args.ram_limit_gb,
        f"cycle channel {channel}",
    )
    weights_cache: dict[tuple[int, int, int], np.ndarray] = {}
    start = time.time()

    for idx, tile in enumerate(tiles, start=1):
        tile_start = time.time()
        print(f"[cycle-stream] channel {channel}: apply tile {idx}/{len(tiles)} {tile}", flush=True)
        tile_manifest = load_tile_manifest(transforms_dir, idx)
        if tile_manifest.get("output_skipped"):
            print(
                f"[cycle-stream] channel {channel}: tile {idx}/{len(tiles)} skipped because "
                f"fixed DAPI tile was background ({tile_manifest.get('skip_reason', '')})",
                flush=True,
            )
            record_timing(timings, f"channel {channel} tile {idx} skipped background", tile_start, channel=channel, tile_index=idx)
            continue
        stage_start = time.time()
        fixed_domain, moving_global = globally_align_tile(
            fixed_spec,
            moving_spec,
            tile,
            global_transforms,
            args.margin_xy,
            args.margin_z,
            args.progress_interval,
        )
        record_timing(
            timings,
            f"channel {channel} tile {idx} load/apply global transform",
            stage_start,
            channel=channel,
            tile_index=idx,
        )
        local_transforms = load_tile_transforms(tile_manifest, idx)
        if local_transforms:
            stage_start = time.time()
            with progress_heartbeat(f"apply ch{channel} tile {idx} local SyN", args.progress_interval):
                refined = ants.apply_transforms(
                    fixed=fixed_domain,
                    moving=moving_global,
                    transformlist=local_transforms,
                    interpolator="linear",
                    singleprecision=True,
                )
            record_timing(timings, f"channel {channel} tile {idx} apply local SyN", stage_start, channel=channel, tile_index=idx)
        else:
            refined = moving_global

        stage_start = time.time()
        arr = ants_to_zyx_float(refined)
        if arr.shape != tile.shape:
            raise ValueError(f"Channel {channel} tile {idx} output shape mismatch: expected {tile.shape}, got {arr.shape}")
        tile_weights = weights_cache.get(tile.shape)
        if tile_weights is None:
            tile_weights = raised_cosine_weights(tile.shape, overlap)
            weights_cache[tile.shape] = tile_weights
        accum[tile.z0 : tile.z1, tile.y0 : tile.y1, tile.x0 : tile.x1] += arr * tile_weights
        weight[tile.z0 : tile.z1, tile.y0 : tile.y1, tile.x0 : tile.x1] += tile_weights
        flush_blend_buffers(accum, weight)
        record_timing(timings, f"channel {channel} tile {idx} blend", stage_start, channel=channel, tile_index=idx)

        del fixed_domain, moving_global, refined, arr
        gc.collect()
        elapsed = time.time() - start
        eta = (elapsed / idx) * (len(tiles) - idx)
        print(
            f"[cycle-stream] channel {channel}: tile {idx}/{len(tiles)} done in "
            f"{format_duration(time.time() - tile_start)}; elapsed={format_duration(elapsed)} eta={format_duration(eta)}",
            flush=True,
        )
        record_timing(timings, f"channel {channel} tile {idx} total", tile_start, channel=channel, tile_index=idx)

    stage_start = time.time()
    output_path = write_registered_channel_from_buffers(
        args,
        out_dir,
        scratch_root,
        channel,
        args.ref_ch,
        args.reg_ch,
        accum,
        weight,
    )
    record_timing(timings, f"channel {channel} write registered TIFF", stage_start, channel=channel)
    return {
        "channel": channel,
        "role": "applied_channel",
        "input_source": str(moving_spec.path) if isinstance(moving_spec, SourceSpec) else "",
        "input_source_kind": moving_spec.kind if isinstance(moving_spec, SourceSpec) else "ims",
        "output_tif": str(output_path),
    }


def run_cycle_to_fixed_streaming(args, out_dir: Path, scratch_root: Path, global_tx_dir: Path) -> dict:
    if not (0.0 <= args.bg_threshold_fraction <= 1.0):
        raise ValueError("--bg-threshold-fraction must be in [0, 1]")
    if not (0.0 <= args.min_foreground_fraction <= 1.0):
        raise ValueError("--min-foreground-fraction must be in [0, 1]")
    if args.bg_check_stride < 1:
        raise ValueError("--bg-check-stride must be >= 1")
    if args.channel_offset < 0:
        raise ValueError("--channel-offset must be >= 0")
    if not args.apply_only and args.reg_ch not in args.channels:
        raise ValueError("--reg-ch must be included in --ch when estimating registration so its registered TIFF has an output channel")
    fixed_spec = fixed_source_spec(args)
    source_map = parse_source_map(args.moving_source_map) if args.moving_source_map else {}
    moving_reg_spec = moving_source_spec(args, args.reg_ch, source_map)
    timings: list[dict] = []
    stage_start = time.time()
    global_transforms = load_or_estimate_global_transform(args, fixed_spec, moving_reg_spec, args.reg_ch, global_tx_dir)
    record_timing(timings, f"global {global_transform_label(args.global_tx)} estimate/reuse", stage_start, channel=args.reg_ch)
    global_qc_result = None
    if args.global_qc and not args.apply_only:
        stage_start = time.time()
        print("[cycle-stream] writing global-registration QC before tiled refinement", flush=True)
        global_qc_result = write_global_preview_qc(args, fixed_spec, moving_reg_spec, global_transforms, out_dir, timings)
        record_timing(timings, "global QC before tiling total", stage_start, channel=args.reg_ch)
    background_threshold = None
    if args.skip_bg_tiles and not args.apply_only:
        stage_start = time.time()
        background_threshold = estimate_background_threshold(args, fixed_spec)
        record_timing(timings, "estimate background threshold", stage_start)

    tile_size, overlap = tile_schedule(args, fixed_spec.shape_zyx)
    tiles = generate_tiles(fixed_spec.shape_zyx, tile_size, overlap)
    if args.max_tiles > 0:
        tiles = tiles[: args.max_tiles]
    print(
        f"[cycle-stream] fixed={fixed_spec.path} fixed_ch={args.ref_ch} shape_zyx={fixed_spec.shape_zyx}",
        flush=True,
    )
    print(
        f"[cycle-stream] moving={Path(args.ims).expanduser().resolve()} reg_ch={args.reg_ch} "
        f"shape_zyx={moving_reg_spec.shape_zyx}",
        flush=True,
    )
    print(
        f"[cycle-stream] tiles={len(tiles)} tile={tile_size} overlap={overlap} tile_ds={args.tile_ds} "
        f"channels={args.channels} blend={args.blend} no_local_refine={args.no_local_refine} "
        f"skip_bg_tiles={args.skip_bg_tiles}",
        flush=True,
    )
    final_output_dir = Path(args.final_output_dir).expanduser().resolve() if args.final_output_dir else out_dir
    output_channel_map = {channel: cycle_output_channel(args.channels, channel, args.channel_offset) for channel in args.channels}
    print(
        f"[cycle-stream] final TIFF directory={final_output_dir} prefix={args.tiff_prefix} "
        f"channel_offset={args.channel_offset} output_channel_map={output_channel_map}",
        flush=True,
    )
    if source_map:
        for channel, path in sorted(source_map.items()):
            print(f"[cycle-stream] moving channel {channel} source override: {path}", flush=True)

    transforms_dir = out_dir / "tile_quicksyn_transforms" / f"dapi_ch{args.reg_ch}_to_fixed_ch{args.ref_ch}"
    results: list[dict] = []
    if args.apply_only:
        print("[cycle-stream] apply-only mode: reusing global and tile QuickSyN transforms", flush=True)
    else:
        stage_start = time.time()
        results.append(
            estimate_tiles_and_write_registration_channel(
                args,
                fixed_spec,
                moving_reg_spec,
                tiles,
                overlap,
                global_transforms,
                transforms_dir,
                out_dir,
                scratch_root,
                timings,
                background_threshold,
            )
        )
        record_timing(timings, f"registration channel {args.reg_ch} total", stage_start, channel=args.reg_ch)

    for channel in args.channels:
        if channel == args.reg_ch and not args.apply_only:
            print(f"[cycle-stream] channel {channel}: already written during DAPI registration pass", flush=True)
            continue
        moving_spec = moving_source_spec(args, channel, source_map)
        stage_start = time.time()
        results.append(
            apply_saved_tile_transforms_to_channel(
                args,
                fixed_spec,
                moving_spec,
                channel,
                tiles,
                overlap,
                global_transforms,
                transforms_dir,
                out_dir,
                scratch_root,
                timings,
            )
        )
        record_timing(timings, f"channel {channel} total", stage_start, channel=channel)

    return {
        "mode": "cycle_to_fixed_streaming",
        "fixed_ims": str(Path(args.fixed_ims).expanduser().resolve()),
        "moving_ims": str(Path(args.ims).expanduser().resolve()),
        "fixed_channel": args.ref_ch,
        "moving_registration_channel": args.reg_ch,
        "channels": args.channels,
        "global_transforms": global_transforms,
        "tile_transform_dir": str(transforms_dir),
        "tile_size_zyx": tile_size,
        "overlap_zyx": overlap,
        "tile_downsample_factor": args.tile_ds,
        "no_local_refine": args.no_local_refine,
        "full_z_tiles": args.full_z,
        "tile_count": len(tiles),
        "blend": args.blend,
        "final_output_dir": str(final_output_dir),
        "tiff_prefix": args.tiff_prefix,
        "channel_offset": args.channel_offset,
        "output_channel_map": output_channel_map,
        "skip_background_tiles": args.skip_bg_tiles,
        "background_threshold": background_threshold,
        "background_threshold_fraction": args.bg_threshold_fraction,
        "min_foreground_fraction": args.min_foreground_fraction,
        "global_qc_before_tiling": global_qc_result,
        "timings": timings,
        "results": results,
    }


def main() -> int:
    args = build_parser().parse_args()
    if args.open_qc:
        args.qc = True
    if args.qc_ds < 1:
        raise ValueError("--qc-ds must be >= 1")

    out_dir = Path(args.out).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.tile_map_only:
        write_tile_map_only(args, out_dir)
        return 0

    configure_ants_runtime(args.threads)
    import ants  # noqa: F401

    scratch_root = Path(args.tmp_dir).expanduser().resolve() if args.tmp_dir else out_dir / "scratch"
    scratch_root.mkdir(parents=True, exist_ok=True)
    global_tx_dir = Path(args.global_transform_dir).expanduser().resolve() if args.global_transform_dir else out_dir / "transforms"

    if args.global_only:
        run_start = time.time()
        result = run_global_preview(args, out_dir, global_tx_dir)
        total_seconds = time.time() - run_start
        result["elapsed_seconds"] = total_seconds
        result["elapsed"] = format_duration(total_seconds)
        timing_path = out_dir / "stage_timing_global_preview.json"
        timing_path.write_text(
            json.dumps(
                {
                    "timings": result.get("timings", []),
                    "total_seconds": total_seconds,
                    "total": result["elapsed"],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        result["timing_log"] = str(timing_path)
        manifest_path = out_dir / "global_preview_manifest.json"
        manifest_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(f"[global-preview] wrote timing log: {timing_path}", flush=True)
        print(f"[global-preview] wrote manifest: {manifest_path}", flush=True)
        print(f"[global-preview] done in {result['elapsed']}", flush=True)
        return 0

    if args.fixed_ims:
        run_start = time.time()
        result = run_cycle_to_fixed_streaming(args, out_dir, scratch_root, global_tx_dir)
        total_seconds = time.time() - run_start
        result["elapsed_seconds"] = total_seconds
        result["elapsed"] = format_duration(total_seconds)
        timing_path = out_dir / "stage_timing_cycle_tiled_quicksyn.json"
        timing_path.write_text(
            json.dumps(
                {
                    "timings": result.get("timings", []),
                    "total_seconds": total_seconds,
                    "total": result["elapsed"],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        result["timing_log"] = str(timing_path)
        manifest_path = out_dir / "cycle_tiled_quicksyn_manifest.json"
        manifest_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(f"[cycle-stream] wrote timing log: {timing_path}", flush=True)
        print(f"[cycle-stream] wrote manifest: {manifest_path}", flush=True)
        print(f"[cycle-stream] done in {result['elapsed']}", flush=True)
        return 0

    if args.apply_only:
        raise ValueError("--apply-only requires --fixed-ims cycle-to-reference mode")
    if args.moving_source_map:
        raise ValueError("--moving-source-map requires --fixed-ims cycle-to-reference mode")

    ref_spec = channel_spec(args, args.ref_ch)
    print(
        f"[tiled-syn] reference channel={args.ref_ch} dataset={ref_spec.dataset_key} "
        f"shape_zyx={ref_spec.shape_zyx} spacing_xyz={ref_spec.spacing_xyz}",
        flush=True,
    )
    print(
        "[tiled-syn] NOTE: this workflow is disk-heavy. It uses full-volume float32 "
        "accumulation and weight scratch arrays for overlap blending.",
        flush=True,
    )

    run_start = time.time()
    results = []
    for channel in args.channels:
        results.append(run_channel(args, ref_spec, channel, out_dir, scratch_root, global_tx_dir))
    manifest = {
        "input_ims": str(Path(args.ims).expanduser().resolve()),
        "reference_channel": args.ref_ch,
        "channels": args.channels,
        "results": results,
        "elapsed": format_duration(time.time() - run_start),
    }
    manifest_path = out_dir / "tiled_quicksyn_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[tiled-syn] wrote manifest: {manifest_path}", flush=True)
    print(f"[tiled-syn] done in {manifest['elapsed']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
