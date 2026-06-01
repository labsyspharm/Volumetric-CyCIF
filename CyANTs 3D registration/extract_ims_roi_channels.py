#!/usr/bin/env python3
"""Extract one ROI box from multiple channels in an Imaris .ims file to NRRD."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Tuple

import SimpleITK as sitk
from roifile import ImagejRoi, roiread

from cyants_io import (
    Tile,
    build_ims_volume_spec,
    parse_spacing_override,
    read_ims_tile_from_spec,
    stem_for_output,
)


def parse_roi(value: str) -> Tuple[int, int, int, int, int, int]:
    """Parse x0,y0,z0,x1,y1,z1 with x1/y1/z1 exclusive."""
    if Path(value).exists():
        payload = json.loads(Path(value).read_text(encoding="utf-8"))
        if all(k in payload for k in ("x0", "y0", "z0", "x1", "y1", "z1")):
            vals = tuple(int(payload[k]) for k in ("x0", "y0", "z0", "x1", "y1", "z1"))
        elif all(k in payload for k in ("x", "y", "z", "size_x", "size_y", "size_z")):
            x, y, z = (int(payload[k]) for k in ("x", "y", "z"))
            sx, sy, sz = (int(payload[k]) for k in ("size_x", "size_y", "size_z"))
            vals = (x, y, z, x + sx, y + sy, z + sz)
        else:
            raise ValueError("ROI JSON must contain x0,y0,z0,x1,y1,z1 or x,y,z,size_x,size_y,size_z")
    else:
        parts = [p.strip() for p in value.split(",") if p.strip()]
        if len(parts) != 6:
            raise argparse.ArgumentTypeError("--roi must be x0,y0,z0,x1,y1,z1 or a ROI JSON path")
        vals = tuple(int(p) for p in parts)
    x0, y0, z0, x1, y1, z1 = vals
    if x1 <= x0 or y1 <= y0 or z1 <= z0:
        raise argparse.ArgumentTypeError("ROI end indices must be greater than start indices")
    if min(vals) < 0:
        raise argparse.ArgumentTypeError("ROI indices must be non-negative")
    return vals


def parse_z_range(value: str) -> Tuple[int, int]:
    parts = [p.strip() for p in value.split(",") if p.strip()]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("--z-range must be first_slice,last_slice")
    first_slice, last_slice = (int(p) for p in parts)
    if first_slice < 1 or last_slice < first_slice:
        raise argparse.ArgumentTypeError("--z-range uses ImageJ 1-based slices and must satisfy last >= first >= 1")
    return (first_slice - 1, last_slice)


def roi_position_z0(roi: ImagejRoi) -> int:
    position = int(getattr(roi, "position", 0) or 0)
    z_position = int(getattr(roi, "z_position", 0) or 0)
    if z_position > 0:
        return z_position - 1
    if position > 0:
        return position - 1
    return -1


def parse_imagej_roi(
    path: Path,
    z_range: Tuple[int, int] | None,
    default_z_range: Tuple[int, int] | None = None,
) -> Tuple[int, int, int, int, int, int]:
    rois = roiread(str(path))
    if isinstance(rois, ImagejRoi):
        roi_list = [rois]
    else:
        roi_list = list(rois)
    if not roi_list:
        raise ValueError(f"No ROIs found in {path}")

    lefts = []
    rights = []
    tops = []
    bottoms = []
    zs = []
    for roi in roi_list:
        left = int(getattr(roi, "left"))
        top = int(getattr(roi, "top"))
        right = int(getattr(roi, "right"))
        bottom = int(getattr(roi, "bottom"))
        lefts.append(left)
        rights.append(right)
        tops.append(top)
        bottoms.append(bottom)
        z0 = roi_position_z0(roi)
        if z0 >= 0:
            zs.append(z0)

    x0 = min(lefts)
    x1 = max(rights)
    y0 = min(tops)
    y1 = max(bottoms)
    if z_range is not None:
        z0, z1 = z_range
    elif zs:
        z0 = min(zs)
        z1 = max(zs) + 1
    elif default_z_range is not None:
        z0, z1 = default_z_range
    else:
        raise ValueError(
            "The ImageJ ROI file did not contain z positions. "
            "Pass --z-range first_slice,last_slice, using ImageJ 1-based slice numbers."
        )
    return validate_roi((x0, y0, z0, x1, y1, z1))


def validate_roi(vals: Tuple[int, int, int, int, int, int]) -> Tuple[int, int, int, int, int, int]:
    x0, y0, z0, x1, y1, z1 = vals
    if x1 <= x0 or y1 <= y0 or z1 <= z0:
        raise argparse.ArgumentTypeError("ROI end indices must be greater than start indices")
    if min(vals) < 0:
        raise argparse.ArgumentTypeError("ROI indices must be non-negative")
    return vals


def parse_channels(value: str) -> List[int]:
    channels = set()
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            start_s, end_s = item.split("-", 1)
            start, end = int(start_s), int(end_s)
            if end < start:
                raise argparse.ArgumentTypeError(f"Invalid channel range: {item}")
            channels.update(range(start, end + 1))
        else:
            channels.add(int(item))
    if not channels or min(channels) < 0:
        raise argparse.ArgumentTypeError("--channels must contain non-negative channel indices")
    return sorted(channels)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract a voxel ROI from selected channels of a .ims file and write channel NRRDs."
    )
    parser.add_argument("--input-ims", required=True, help="Source Imaris .ims file")
    parser.add_argument(
        "--roi",
        type=parse_roi,
        help="ROI as x0,y0,z0,x1,y1,z1 with end indices exclusive, or path to ROI JSON",
    )
    parser.add_argument(
        "--roi-imagej",
        default="",
        help="Optional ImageJ/Fiji .roi or ROI .zip file. Its bounding box is used as the ROI.",
    )
    parser.add_argument(
        "--z-range",
        type=parse_z_range,
        default=None,
        help=(
            "Optional ImageJ 1-based z slice range as first_slice,last_slice. "
            "If omitted for --roi-imagej, the full .ims z depth is used."
        ),
    )
    parser.add_argument("--channels", required=True, type=parse_channels, help="Channels, for example 0,1,2 or 0-4")
    parser.add_argument("--output-dir", required=True, help="Directory for ROI channel NRRDs")
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
    return parser


def main() -> int:
    args = build_parser().parse_args()
    input_path = Path(args.input_ims).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    first_channel = args.channels[0]
    first_dataset_path = (
        args.ims_dataset_template.format(
            level=args.ims_resolution_level,
            timepoint=args.ims_timepoint,
            channel=first_channel,
        )
        if args.ims_dataset_template
        else ""
    )
    first_spec = build_ims_volume_spec(
        path=input_path,
        res_level=args.ims_resolution_level,
        timepoint=args.ims_timepoint,
        channel=first_channel,
        dataset_path=first_dataset_path,
        axis_order=args.ims_axis_order,
    )
    full_z_range = (0, first_spec.shape_zyx[0])

    if args.roi_imagej:
        x0, y0, z0, x1, y1, z1 = parse_imagej_roi(
            Path(args.roi_imagej).expanduser().resolve(),
            args.z_range,
            default_z_range=full_z_range,
        )
    elif args.roi:
        x0, y0, z0, x1, y1, z1 = args.roi
    else:
        raise ValueError("Provide either --roi or --roi-imagej")
    tile = Tile(z0=z0, z1=z1, y0=y0, y1=y1, x0=x0, x1=x1)
    prefix = args.prefix or stem_for_output(input_path)

    roi_manifest = {
        "input_ims": str(input_path),
        "roi_imagej": str(Path(args.roi_imagej).expanduser().resolve()) if args.roi_imagej else "",
        "z_defaulted_to_full_depth": bool(args.roi_imagej and args.z_range is None),
        "roi_xyz_exclusive": {"x0": x0, "y0": y0, "z0": z0, "x1": x1, "y1": y1, "z1": z1},
        "channels": args.channels,
        "outputs": [],
    }

    for channel in args.channels:
        dataset_path = (
            args.ims_dataset_template.format(
                level=args.ims_resolution_level,
                timepoint=args.ims_timepoint,
                channel=channel,
            )
            if args.ims_dataset_template
            else ""
        )
        if channel == first_channel and dataset_path == first_dataset_path:
            spec = first_spec
        else:
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
        img = read_ims_tile_from_spec(spec, tile)
        out_path = output_dir / f"{prefix}_ch{channel}_roi.nrrd"
        sitk.WriteImage(img, str(out_path))
        print(f"Wrote channel {channel} ROI: {out_path}", flush=True)
        roi_manifest["outputs"].append({"channel": channel, "path": str(out_path)})

    manifest_path = output_dir / f"{prefix}_roi_manifest.json"
    manifest_path.write_text(json.dumps(roi_manifest, indent=2), encoding="utf-8")
    print(f"Wrote ROI manifest: {manifest_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
