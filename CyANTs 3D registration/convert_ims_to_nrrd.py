#!/usr/bin/env python3
# Developed by Alex Wong
# Cite: Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in
# Preprint: https://doi.org/10.64898/2026.05.17.725158
# Registration method: ANTsX/ANTsPy - https://github.com/ANTsX/ANTsPy

"""Convert one Imaris .ims channel to a full-resolution NRRD volume."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import SimpleITK as sitk

from cyants_io import Tile, build_ims_volume_spec, parse_spacing_override, read_ims_tile_from_spec


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert one channel from an Imaris .ims volume to NRRD.")
    parser.add_argument("--ims", required=True, help="Input .ims path")
    parser.add_argument("--out", required=True, help="Output .nrrd path")
    parser.add_argument("--ch", "--channel", type=int, default=0, help="Channel index. Default: 0")
    parser.add_argument("--level", "--ims-resolution-level", type=int, default=0, help="Resolution level. Default: 0")
    parser.add_argument("--timepoint", "--ims-timepoint", type=int, default=0, help="Timepoint index. Default: 0")
    parser.add_argument("--ims-dataset-path", default="", help="Optional explicit Imaris HDF5 dataset path")
    parser.add_argument("--ims-axis-order", choices=["auto", "zyx", "xyz"], default="auto", help="Dataset axis order. Default: auto")
    parser.add_argument("--spacing", type=parse_spacing_override, default=None, help="Override output spacing, for example 0.711")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    input_path = Path(args.ims).expanduser().resolve()
    output_path = Path(args.out).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input .ims does not exist: {input_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    start = time.time()
    spec = build_ims_volume_spec(
        input_path,
        res_level=args.level,
        timepoint=args.timepoint,
        channel=args.ch,
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
    print(
        f"[convert] source={input_path} dataset={spec.dataset_key} "
        f"shape_zyx={spec.shape_zyx} spacing_xyz={spec.spacing_xyz}",
        flush=True,
    )
    print("[convert] loading full channel volume into memory", flush=True)
    full_tile = Tile(z0=0, z1=spec.shape_zyx[0], y0=0, y1=spec.shape_zyx[1], x0=0, x1=spec.shape_zyx[2])
    img = read_ims_tile_from_spec(spec, full_tile)
    print(f"[convert] writing NRRD: {output_path}", flush=True)
    sitk.WriteImage(img, str(output_path), True)
    verify = sitk.ReadImage(str(output_path))
    print(
        f"[convert] done in {time.time() - start:.1f}s; "
        f"size_xyz={verify.GetSize()} spacing_xyz={verify.GetSpacing()}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
