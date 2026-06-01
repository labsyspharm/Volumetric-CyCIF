#!/usr/bin/env python3
"""Estimate whole-.ims intra-cycle channel alignment transforms."""

from __future__ import annotations

import argparse
import gc
import json
import time
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile

from ants_roi_quicksyn import (
    ants_to_sitk_uint16,
    configure_ants_runtime,
    export_center_slice_qc,
    format_duration,
    progress_heartbeat,
)
from extract_ims_roi_channels import parse_channels
from cyants_io import (
    Tile,
    build_ims_volume_spec,
    parse_spacing_override,
    read_ims_downsampled_from_spec,
    read_ims_tile_from_spec,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate channel-to-channel intra-cycle transforms from the whole .ims field. "
            "This reads downsampled whole-volume proxies, so it is separate from ROI extraction."
        )
    )
    parser.add_argument("--input-ims", "--ims", required=True, help="Source multi-channel .ims file")
    parser.add_argument("--reference-ims", "--ref-ims", default="", help="Optional reference .ims path. Default: use --ims.")
    parser.add_argument(
        "--moving-source-map",
        "--source-map",
        default="",
        help="Optional semicolon-separated moving channel path overrides, for example channel=path.",
    )
    parser.add_argument("--reference-channel", "--ref-ch", type=int, default=0, help="Reference channel. Default: 0")
    parser.add_argument("--channels", "--ch", required=True, type=parse_channels, help="Channels to align, for example 1,3")
    parser.add_argument("--output-dir", "--out", required=True, help="Directory for transforms and QC")
    parser.add_argument("--transform", "--tx", default="TRSAA", help="ANTs transform type. Default: TRSAA")
    parser.add_argument(
        "--apply-only",
        "--ao",
        action="store_true",
        help="Reuse existing whole-.ims transform manifests and only write requested TIFF outputs.",
    )
    parser.add_argument(
        "--downsample-factor",
        "--ds",
        type=int,
        default=8,
        help="Whole-volume stride used for transform estimation. Default: 8",
    )
    parser.add_argument("--ims-resolution-level", type=int, default=0, help=".ims resolution level. Default: 0")
    parser.add_argument("--ims-timepoint", type=int, default=0, help=".ims timepoint. Default: 0")
    parser.add_argument(
        "--ims-dataset-template",
        default="",
        help="Optional HDF5 dataset template with {level}, {timepoint}, and {channel}",
    )
    parser.add_argument(
        "--ims-axis-order",
        choices=["auto", "zyx", "xyz"],
        default="auto",
        help="Axis order for .ims dataset interpretation. Default: auto",
    )
    parser.add_argument("--spacing", type=parse_spacing_override, default=None, help="Override spacing, e.g. 0.711")
    parser.add_argument("--threads", "--t", type=int, default=32, help="CPU threads for ANTs/ITK. Default: 32")
    parser.add_argument(
        "--progress-interval",
        "--pi",
        type=float,
        default=30.0,
        help="Heartbeat interval in seconds. Use 0 to disable. Default: 30.",
    )
    parser.add_argument("--no-qc", dest="qc", action="store_false", help="Disable center-slice QC PNG output")
    parser.set_defaults(qc=True)
    parser.add_argument("--open-qc", action="store_true", help="Open QC PNGs after writing them")
    parser.add_argument("--qc-max-panel-side", type=int, default=900, help="Maximum QC panel side. Default: 900")
    parser.add_argument(
        "--write-aligned-tif",
        "--tif",
        action="store_true",
        help=(
            "Write the aligned downsampled whole-.ims proxy for each moving channel as a TIFF stack. "
            "This is downsampled by --ds, not a full-resolution 200 GB volume."
        ),
    )
    parser.add_argument(
        "--aligned-tif-dir",
        default="",
        help="Directory for --write-aligned-tif outputs. Default: <output-dir>/aligned_tif",
    )
    parser.add_argument(
        "--write-fullres-aligned-tif",
        "--fullres-tif",
        "--full-tif",
        action="store_true",
        help=(
            "Write full-resolution aligned TIFF stacks by streaming the .ims in Z chunks. "
            "This can be very slow and creates very large BigTIFF files."
        ),
    )
    parser.add_argument(
        "--fullres-tif-dir",
        default="",
        help="Directory for --write-fullres-aligned-tif outputs. Default: <output-dir>/fullres_aligned_tif",
    )
    parser.add_argument(
        "--fullres-chunk-z",
        "--cz",
        type=int,
        default=8,
        help="Number of full-resolution Z planes to resample/write per chunk. Default: 8.",
    )
    parser.add_argument(
        "--fullres-margin-xy",
        "--mxy",
        type=int,
        default=512,
        help="XY moving-image margin in voxels around each output chunk. Default: 512.",
    )
    parser.add_argument(
        "--fullres-margin-z",
        "--mz",
        type=int,
        default=16,
        help="Z moving-image margin in voxels around each output chunk. Default: 16.",
    )
    parser.add_argument(
        "--uint16-scaling",
        choices=["minmax", "robust", "clip"],
        default="clip",
        help="How to convert aligned float output to uint16 TIFF. Default: clip.",
    )
    return parser


def parse_source_map(value: str) -> dict[int, Path]:
    mapping: dict[int, Path] = {}
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Invalid --source-map entry: {item}")
        channel_s, path_s = item.split("=", 1)
        path = Path(path_s.strip().strip('"')).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Source override does not exist for channel {channel_s}: {path}")
        if path.suffix.lower() != ".ims":
            raise ValueError(f"Intracycle source override must be a .ims file: {path}")
        mapping[int(channel_s.strip())] = path
    return mapping


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
    source_map = getattr(args, "source_map_parsed", {})
    if channel == args.reference_channel:
        input_path = Path(args.reference_ims or args.input_ims).expanduser().resolve()
    else:
        input_path = source_map.get(channel, Path(args.input_ims).expanduser().resolve())
    spec = None
    last_error = None
    for attempt in range(1, 4):
        try:
            spec = build_ims_volume_spec(
                input_path,
                res_level=args.ims_resolution_level,
                timepoint=args.ims_timepoint,
                channel=channel,
                dataset_path=dataset_path,
                axis_order=args.ims_axis_order,
            )
            break
        except OSError as exc:
            last_error = exc
            print(
                f"[whole-ims] WARNING failed to open .ims for channel {channel} "
                f"(attempt {attempt}/3): {exc}",
                flush=True,
            )
            time.sleep(3)
    if spec is None:
        raise OSError(
            f"Could not open .ims after retries: {input_path}. "
            "If the error says 'file signature not found', check that this path is the real .ims file "
            "and that the network drive is not serving a shortcut/partial file."
        ) from last_error
    if args.spacing is not None:
        spec = type(spec)(
            path=spec.path,
            dataset_key=spec.dataset_key,
            axis_order=spec.axis_order,
            shape_zyx=spec.shape_zyx,
            spacing_xyz=args.spacing,
        )
    return spec


def sitk_to_ants(img: sitk.Image):
    import ants

    return ants.from_numpy(
        sitk.GetArrayFromImage(img).T,
        origin=img.GetOrigin(),
        spacing=img.GetSpacing(),
        direction=np.asarray(img.GetDirection()).reshape((3, 3)),
    )


def write_aligned_tif(img, path: Path, uint16_scaling: str, progress_interval: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[whole-ims] writing aligned downsampled TIFF: {path}", flush=True)
    with progress_heartbeat("whole-ims write aligned tif", progress_interval):
        sitk_img = ants_to_sitk_uint16(img, uint16_scaling)
        arr = sitk.GetArrayFromImage(sitk_img)
        tifffile.imwrite(path, arr, bigtiff=True, imagej=True, metadata={"axes": "ZYX"})
    print(f"[whole-ims] wrote aligned downsampled TIFF: {path}", flush=True)


def domain_chunk_ants(spec, tile: Tile):
    arr = np.zeros((tile.z1 - tile.z0, tile.y1 - tile.y0, tile.x1 - tile.x0), dtype=np.float32)
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


def padded_chunk_tile(base_tile: Tile, spec, margin_xy: int, margin_z: int) -> Tile:
    z_size, y_size, x_size = spec.shape_zyx
    return Tile(
        z0=max(0, base_tile.z0 - margin_z),
        z1=min(z_size, base_tile.z1 + margin_z),
        y0=max(0, base_tile.y0 - margin_xy),
        y1=min(y_size, base_tile.y1 + margin_xy),
        x0=max(0, base_tile.x0 - margin_xy),
        x1=min(x_size, base_tile.x1 + margin_xy),
    )


def ants_to_uint16_zyx(img, scaling: str) -> np.ndarray:
    sitk_img = ants_to_sitk_uint16(img, scaling)
    return sitk.GetArrayFromImage(sitk_img)


def write_fullres_aligned_tif(
    channel: int,
    reference_spec,
    moving_spec,
    transformlist: list[str],
    output_path: Path,
    chunk_z: int,
    margin_xy: int,
    margin_z: int,
    uint16_scaling: str,
    progress_interval: float,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    z_size, y_size, x_size = reference_spec.shape_zyx
    chunk_z = max(1, int(chunk_z))
    margin_xy = max(0, int(margin_xy))
    margin_z = max(0, int(margin_z))
    n_chunks = (z_size + chunk_z - 1) // chunk_z
    start = time.time()
    print(
        f"[whole-ims] channel {channel}: writing full-resolution aligned BigTIFF: {output_path}",
        flush=True,
    )
    print(
        f"[whole-ims] channel {channel}: fullres shape_zyx={(z_size, y_size, x_size)} "
        f"chunk_z={chunk_z} chunks={n_chunks} margin_xy={margin_xy} margin_z={margin_z}",
        flush=True,
    )
    with tifffile.TiffWriter(output_path, bigtiff=True) as writer:
        for chunk_idx, z0 in enumerate(range(0, z_size, chunk_z), start=1):
            z1 = min(z_size, z0 + chunk_z)
            fixed_tile = Tile(z0=z0, z1=z1, y0=0, y1=y_size, x0=0, x1=x_size)
            moving_tile = padded_chunk_tile(fixed_tile, moving_spec, margin_xy, margin_z)
            print(
                f"[whole-ims] channel {channel}: chunk {chunk_idx}/{n_chunks} "
                f"fixed_tile={fixed_tile} moving_tile={moving_tile}",
                flush=True,
            )
            with progress_heartbeat(f"whole-ims fullres load ch{channel} chunk {chunk_idx}", progress_interval):
                moving_sitk = read_ims_tile_from_spec(moving_spec, moving_tile)
            with progress_heartbeat(f"whole-ims fullres apply ch{channel} chunk {chunk_idx}", progress_interval):
                fixed_domain = domain_chunk_ants(reference_spec, fixed_tile)
                moving_chunk = sitk_to_ants(moving_sitk)
                aligned = __import__("ants").apply_transforms(
                    fixed=fixed_domain,
                    moving=moving_chunk,
                    transformlist=transformlist,
                    interpolator="linear",
                    singleprecision=True,
                )
            with progress_heartbeat(f"whole-ims fullres write ch{channel} chunk {chunk_idx}", progress_interval):
                arr = ants_to_uint16_zyx(aligned, uint16_scaling)
                for plane in arr:
                    writer.write(plane, photometric="minisblack", metadata=None, contiguous=True)
            del moving_sitk, fixed_domain, moving_chunk, aligned, arr
            gc.collect()
            elapsed = time.time() - start
            eta = (elapsed / chunk_idx) * (n_chunks - chunk_idx) if chunk_idx else 0.0
            print(
                f"[whole-ims] channel {channel}: chunk {chunk_idx}/{n_chunks} done; "
                f"elapsed={format_duration(elapsed)} eta={format_duration(eta)}",
                flush=True,
            )
    print(
        f"[whole-ims] channel {channel}: wrote full-resolution aligned BigTIFF in "
        f"{format_duration(time.time() - start)}: {output_path}",
        flush=True,
    )


def transform_manifest_path(tx_dir: Path, channel: int, reference_channel: int) -> Path:
    return tx_dir / f"ch{channel}_to_ch{reference_channel}_wholeims_fwdtransforms.json"


def load_existing_transform(tx_dir: Path, channel: int, reference_channel: int) -> tuple[list[str], dict]:
    path = transform_manifest_path(tx_dir, channel, reference_channel)
    if not path.exists():
        raise FileNotFoundError(f"Existing transform manifest not found for --apply-only: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    transforms = payload.get("fwdtransforms")
    if not transforms:
        raise ValueError(f"Transform manifest does not contain fwdtransforms: {path}")
    print(f"[whole-ims] channel {channel}: reusing transform manifest: {path}", flush=True)
    return list(transforms), payload


def main() -> int:
    args = build_parser().parse_args()
    args.source_map_parsed = parse_source_map(args.moving_source_map)
    configure_ants_runtime(args.threads)
    import ants

    out_dir = Path(args.output_dir).expanduser().resolve()
    tx_dir = out_dir / "transforms"
    qc_dir = out_dir / "qc"
    tif_dir = Path(args.aligned_tif_dir).expanduser().resolve() if args.aligned_tif_dir else out_dir / "aligned_tif"
    fullres_tif_dir = (
        Path(args.fullres_tif_dir).expanduser().resolve()
        if args.fullres_tif_dir
        else out_dir / "fullres_aligned_tif"
    )
    tx_dir.mkdir(parents=True, exist_ok=True)
    if args.qc:
        qc_dir.mkdir(parents=True, exist_ok=True)
    if args.write_aligned_tif:
        tif_dir.mkdir(parents=True, exist_ok=True)
    if args.write_fullres_aligned_tif:
        fullres_tif_dir.mkdir(parents=True, exist_ok=True)

    stride = max(1, int(args.downsample_factor))
    stride_zyx = (stride, stride, stride)
    start = time.time()

    print("[whole-ims] starting whole-.ims intra-cycle alignment", flush=True)
    print(f"[whole-ims] input={Path(args.input_ims).expanduser().resolve()}", flush=True)
    if args.reference_ims:
        print(f"[whole-ims] reference override={Path(args.reference_ims).expanduser().resolve()}", flush=True)
    if args.source_map_parsed:
        print(
            "[whole-ims] moving source overrides="
            + ", ".join(f"ch{channel}:{path}" for channel, path in sorted(args.source_map_parsed.items())),
            flush=True,
        )
    print(f"[whole-ims] reference channel={args.reference_channel} channels={args.channels}", flush=True)
    print(f"[whole-ims] transform={args.transform} stride_zyx={stride_zyx} threads={args.threads}", flush=True)
    print(f"[whole-ims] output dir={out_dir}", flush=True)
    if args.write_aligned_tif:
        print(
            f"[whole-ims] aligned TIFF output enabled: {tif_dir} "
            f"(downsampled by {stride}x, uint16_scaling={args.uint16_scaling})",
            flush=True,
        )
    if args.write_fullres_aligned_tif:
        print(
            f"[whole-ims] full-resolution aligned TIFF output enabled: {fullres_tif_dir} "
            f"chunk_z={args.fullres_chunk_z} margin_xy={args.fullres_margin_xy} "
            f"margin_z={args.fullres_margin_z} uint16_scaling={args.uint16_scaling}",
            flush=True,
        )
    if args.apply_only:
        print("[whole-ims] apply-only mode: reusing existing transforms; skipping registration", flush=True)
    if args.apply_only and not (args.write_aligned_tif or args.write_fullres_aligned_tif):
        raise ValueError("--apply-only needs --tif or --full-tif so there is an output to write")
    skip_apply_only_qc = args.apply_only and args.write_fullres_aligned_tif and not args.write_aligned_tif
    if skip_apply_only_qc and args.qc:
        print("[whole-ims] apply-only full-res export: skipping downsampled QC reads; add --tif for preview TIFF/QC", flush=True)
        args.qc = False

    ref_spec = channel_spec(args, args.reference_channel)
    print(
        f"[whole-ims] reference dataset={ref_spec.dataset_key} "
        f"shape_zyx={ref_spec.shape_zyx} spacing_xyz={ref_spec.spacing_xyz}",
        flush=True,
    )
    fixed = None
    if not args.apply_only or args.write_aligned_tif or args.qc:
        with progress_heartbeat("whole-ims load reference downsampled", args.progress_interval):
            ref_sitk = read_ims_downsampled_from_spec(ref_spec, stride_zyx)
        fixed = sitk_to_ants(ref_sitk)
        del ref_sitk
        gc.collect()

    manifests = {}
    for channel in args.channels:
        if channel == args.reference_channel:
            print(f"[whole-ims] skipping reference channel {channel}", flush=True)
            continue
        channel_start = time.time()
        moving_spec = channel_spec(args, channel)
        print(
            f"[whole-ims] channel {channel}: dataset={moving_spec.dataset_key} "
            f"shape_zyx={moving_spec.shape_zyx} spacing_xyz={moving_spec.spacing_xyz}",
            flush=True,
        )
        moving = None
        if args.apply_only:
            transform_json = transform_manifest_path(tx_dir, channel, args.reference_channel)
            transformlist, payload = load_existing_transform(tx_dir, channel, args.reference_channel)
        else:
            with progress_heartbeat(f"whole-ims load ch{channel} downsampled", args.progress_interval):
                moving_sitk = read_ims_downsampled_from_spec(moving_spec, stride_zyx)
            moving = sitk_to_ants(moving_sitk)
            del moving_sitk
            gc.collect()

            outprefix = tx_dir / f"ch{channel}_to_ch{args.reference_channel}_wholeims_"
            print(f"[whole-ims] channel {channel}: running ANTs {args.transform}", flush=True)
            with progress_heartbeat(f"whole-ims register ch{channel}", args.progress_interval):
                tx = ants.registration(
                    fixed=fixed,
                    moving=moving,
                    type_of_transform=args.transform,
                    outprefix=str(outprefix),
                    aff_metric="mattes",
                    singleprecision=True,
                    verbose=True,
                )
            transformlist = tx["fwdtransforms"]
            transform_json = transform_manifest_path(tx_dir, channel, args.reference_channel)
            payload = {
                "input_ims": str(Path(args.input_ims).expanduser().resolve()),
                "reference_ims": str(Path(args.reference_ims or args.input_ims).expanduser().resolve()),
                "moving_source_ims": str(moving_spec.path),
                "reference_channel": args.reference_channel,
                "moving_channel": channel,
                "transform": args.transform,
                "downsample_factor": stride,
                "stride_zyx": stride_zyx,
                "fwdtransforms": transformlist,
            }
            transform_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        manifests[str(channel)] = str(transform_json)
        print(f"[whole-ims] channel {channel}: wrote transform manifest: {transform_json}", flush=True)

        aligned_tif_path = ""
        fullres_tif_path = ""
        if args.qc or args.write_aligned_tif:
            if moving is None:
                with progress_heartbeat(f"whole-ims load ch{channel} downsampled", args.progress_interval):
                    moving_sitk = read_ims_downsampled_from_spec(moving_spec, stride_zyx)
                moving = sitk_to_ants(moving_sitk)
                del moving_sitk
                gc.collect()
            if args.qc:
                print(f"[whole-ims] channel {channel}: writing downsampled before/after QC", flush=True)
            if args.write_aligned_tif:
                print(f"[whole-ims] channel {channel}: applying transform for aligned TIFF", flush=True)
            with progress_heartbeat(f"whole-ims qc apply ch{channel}", args.progress_interval):
                aligned = ants.apply_transforms(
                    fixed=fixed,
                    moving=moving,
                    transformlist=transformlist,
                    interpolator="linear",
                    singleprecision=True,
                )
            if args.qc:
                export_center_slice_qc(
                    fixed,
                    moving,
                    aligned,
                    qc_dir,
                    f"whole_ims_ch{channel}_to_ch{args.reference_channel}",
                    args.open_qc,
                    args.qc_max_panel_side,
                )
            if args.write_aligned_tif:
                tif_path = tif_dir / f"ch{channel}_to_ch{args.reference_channel}_wholeims_ds{stride}_aligned.tif"
                write_aligned_tif(aligned, tif_path, args.uint16_scaling, args.progress_interval)
                aligned_tif_path = str(tif_path)
            del aligned

        if aligned_tif_path:
            payload["aligned_downsampled_tif"] = aligned_tif_path
            transform_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        if args.write_fullres_aligned_tif:
            fullres_path = fullres_tif_dir / f"ch{channel}_to_ch{args.reference_channel}_wholeims_fullres_aligned.tif"
            write_fullres_aligned_tif(
                channel,
                ref_spec,
                moving_spec,
                transformlist,
                fullres_path,
                args.fullres_chunk_z,
                args.fullres_margin_xy,
                args.fullres_margin_z,
                args.uint16_scaling,
                args.progress_interval,
            )
            fullres_tif_path = str(fullres_path)
            payload["aligned_fullres_tif"] = fullres_tif_path
            transform_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        if moving is not None:
            del moving
        gc.collect()
        print(f"[whole-ims] channel {channel}: finished in {format_duration(time.time() - channel_start)}", flush=True)

    manifest_path = out_dir / "whole_ims_intracycle_manifest.json"
    manifest_path.write_text(json.dumps({"transforms": manifests}, indent=2), encoding="utf-8")
    print(f"[whole-ims] wrote run manifest: {manifest_path}", flush=True)
    print(f"[whole-ims] done in {format_duration(time.time() - start)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
