#!/usr/bin/env python3
"""Pad/crop a TIFF or Imaris .ims stack to a requested ZYX shape."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile

from ants_roi_quicksyn import format_duration
from cyants_io import Tile, build_ims_volume_spec, parse_spacing_override, read_ims_tile_from_spec


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pad or crop a TIFF stack or .ims channel to an exact Z,Y,X shape.")
    parser.add_argument("--input", "--in", dest="input_path", required=True, help="Input TIFF stack or .ims file")
    parser.add_argument("--output", "--out", required=True, help="Output padded/cropped BigTIFF stack")
    parser.add_argument(
        "--shape",
        required=True,
        help="Target shape as z,y,x, for example 1568,2304,12800",
    )
    parser.add_argument("--pad-value", type=int, default=0, help="Padding value. Default: 0")
    parser.add_argument("--chunk-z", "--cz", type=int, default=8, help="Z planes per read/write chunk. Default: 8")
    parser.add_argument("--ims-channel", "--ch", type=int, default=0, help=".ims channel to read. Default: 0")
    parser.add_argument("--ims-resolution-level", type=int, default=0, help=".ims resolution level. Default: 0")
    parser.add_argument("--ims-timepoint", type=int, default=0, help=".ims timepoint. Default: 0")
    parser.add_argument(
        "--ims-dataset-path",
        "--ims-dset",
        default="",
        help="Optional explicit .ims HDF5 dataset path",
    )
    parser.add_argument(
        "--ims-axis-order",
        choices=["auto", "zyx", "xyz"],
        default="auto",
        help=".ims dataset axis order. Default: auto",
    )
    parser.add_argument("--spacing", type=parse_spacing_override, default=None, help="Optional .ims spacing override")
    return parser


def parse_shape(value: str) -> tuple[int, int, int]:
    parts = [int(v.strip()) for v in value.split(",")]
    if len(parts) != 3 or any(v <= 0 for v in parts):
        raise ValueError("--shape must be z,y,x with positive integers")
    return tuple(parts)


def write_padded_tiff_from_tiff(
    input_path: Path,
    output_path: Path,
    target_shape: tuple[int, int, int],
    chunk_z: int,
    pad_value: int,
    start: float,
) -> None:
    target_z, target_y, target_x = target_shape
    with tifffile.TiffFile(input_path) as tif:
        in_z = len(tif.pages)
        first = tif.pages[0].asarray()
        in_y, in_x = first.shape
        dtype = first.dtype
        print(
            f"[pad-tiff] input TIFF={input_path} shape_zyx=({in_z}, {in_y}, {in_x}) dtype={dtype}",
            flush=True,
        )
        print(
            f"[pad-tiff] output={output_path} target_zyx=({target_z}, {target_y}, {target_x}) "
            f"chunk_z={chunk_z} pad_value={pad_value}",
            flush=True,
        )
        with tifffile.TiffWriter(output_path, bigtiff=True) as writer:
            for z0 in range(0, target_z, chunk_z):
                z1 = min(target_z, z0 + chunk_z)
                for z in range(z0, z1):
                    out = np.full((target_y, target_x), pad_value, dtype=dtype)
                    if z < in_z:
                        plane = tif.pages[z].asarray()
                        copy_y = min(target_y, plane.shape[0])
                        copy_x = min(target_x, plane.shape[1])
                        out[:copy_y, :copy_x] = plane[:copy_y, :copy_x]
                    writer.write(out, photometric="minisblack", metadata=None, contiguous=True)
                print_progress("pad-tiff", z0, z1, target_z, start)


def write_padded_tiff_from_ims(
    args,
    input_path: Path,
    output_path: Path,
    target_shape: tuple[int, int, int],
    chunk_z: int,
    pad_value: int,
    start: float,
) -> None:
    target_z, target_y, target_x = target_shape
    spec = build_ims_volume_spec(
        input_path,
        res_level=args.ims_resolution_level,
        timepoint=args.ims_timepoint,
        channel=args.ims_channel,
        dataset_path=args.ims_dataset_path,
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
    in_z, in_y, in_x = spec.shape_zyx
    print(
        f"[pad-tiff] input IMS={input_path} channel={args.ims_channel} dataset={spec.dataset_key} "
        f"shape_zyx=({in_z}, {in_y}, {in_x})",
        flush=True,
    )
    print(
        f"[pad-tiff] output={output_path} target_zyx=({target_z}, {target_y}, {target_x}) "
        f"chunk_z={chunk_z} pad_value={pad_value}",
        flush=True,
    )
    dtype = None
    with tifffile.TiffWriter(output_path, bigtiff=True) as writer:
        for z0 in range(0, target_z, chunk_z):
            z1 = min(target_z, z0 + chunk_z)
            out = None
            read_z0 = z0
            read_z1 = min(z1, in_z)
            copy_y = min(target_y, in_y)
            copy_x = min(target_x, in_x)
            if read_z0 < read_z1 and copy_y > 0 and copy_x > 0:
                tile = Tile(z0=read_z0, z1=read_z1, y0=0, y1=copy_y, x0=0, x1=copy_x)
                img = read_ims_tile_from_spec(spec, tile)
                arr = sitk.GetArrayFromImage(img)
                if dtype is None:
                    dtype = arr.dtype
                out = np.full((z1 - z0, target_y, target_x), pad_value, dtype=dtype)
                out[: arr.shape[0], : arr.shape[1], : arr.shape[2]] = arr
            if out is None:
                if dtype is None:
                    dtype = np.uint16
                out = np.full((z1 - z0, target_y, target_x), pad_value, dtype=dtype)
            for plane in out:
                writer.write(plane, photometric="minisblack", metadata=None, contiguous=True)
            print_progress("pad-ims", z0, z1, target_z, start)


def print_progress(label: str, z0: int, z1: int, target_z: int, start: float) -> None:
    elapsed = time.time() - start
    done = z1
    eta = (elapsed / done) * (target_z - done) if done else 0.0
    print(
        f"[{label}] wrote z {z0}:{z1} / {target_z}; "
        f"elapsed={format_duration(elapsed)} eta={format_duration(eta)}",
        flush=True,
    )


def main() -> int:
    args = build_parser().parse_args()
    input_path = Path(args.input_path).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    target_z, target_y, target_x = parse_shape(args.shape)
    chunk_z = max(1, int(args.chunk_z))
    output_path.parent.mkdir(parents=True, exist_ok=True)

    start = time.time()
    if input_path.suffix.lower() == ".ims":
        write_padded_tiff_from_ims(args, input_path, output_path, (target_z, target_y, target_x), chunk_z, args.pad_value, start)
    else:
        write_padded_tiff_from_tiff(input_path, output_path, (target_z, target_y, target_x), chunk_z, args.pad_value, start)
    print(f"[pad-tiff] done in {format_duration(time.time() - start)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
