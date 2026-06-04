#!/usr/bin/env python3
"""Register one CyCIF cycle ROI directly from an IMS/TIFF source and XY CSV."""

from __future__ import annotations

import argparse
import gc
import json
import shutil
import time
from pathlib import Path

import numpy as np
import SimpleITK as sitk

from ants_roi_quicksyn import (
    compute_downsampled_metrics_from_images,
    configure_ants_runtime,
    describe_image,
    downsample_by_spacing,
    export_center_slice_qc,
    format_duration,
    load_transform_manifest,
    progress_heartbeat,
    write_metrics_json,
)
from apply_ims_roi_channels import (
    build_channel_spec,
    parse_input_ims_map,
    parse_xy_coordinates_csv,
    print_ants_stats,
    print_sitk_stats,
    read_channel_tile_from_spec,
    tile_from_roi,
    write_registered_outputs,
)
from extract_ims_roi_channels import parse_channels, parse_z_range
from cyants_io import parse_spacing_override


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Stream a CSV-defined ROI directly from one moving cycle .ims/.tif source, "
            "register moving DAPI to an existing reference ROI crop, and apply the transform "
            "to selected channels without writing raw ROI intermediates."
        )
    )
    parser.add_argument("--project-root", "--proj", required=True, help="Registration project root containing reference/ and cycles/")
    parser.add_argument("--cycle-id", "--cycle", required=True, help="Cycle folder id under cycles/, for example cycle_003")
    parser.add_argument("--moving-source", "--ims", required=True, help="Moving cycle .ims source, usually containing all channels")
    parser.add_argument("--roi-csv", required=True, help="XY_Coordinates.csv for the moving cycle; copied into the standard cycle ROI folder if needed")
    parser.add_argument(
        "--fixed-crop",
        default="",
        help=(
            "Cached reference DAPI ROI crop (.nrrd, .tif/.tiff, .nii, or .nii.gz), or a fixed/reference .ims source "
            "when --fixed-roi-csv is supplied. Default: infer one non-registered crop from reference/roi"
        ),
    )
    parser.add_argument(
        "--fixed-roi-csv",
        default="",
        help=(
            "Reference/template ROI CSV used when --fixed-crop points to a full fixed .ims. "
            "The fixed ROI is streamed directly from the .ims using these full-image X,Y coordinates."
        ),
    )
    parser.add_argument("--channels", "--ch", type=parse_channels, default=parse_channels("0-3"), help="Channels to register/output, for example 0-3. Default: 0-3")
    parser.add_argument("--dapi-channel", "--dapi", type=int, default=0, help="Moving DAPI channel used to estimate transforms. Default: 0")
    parser.add_argument("--fixed-dapi-channel", type=int, default=0, help="Fixed/reference DAPI channel when --fixed-crop is a .ims source. Default: 0")
    parser.add_argument(
        "--source-map",
        default="",
        help="Optional channel=path overrides separated by semicolons; values can be .ims, .tif, or .tiff.",
    )
    parser.add_argument("--run-id", "--run", default="run_0001", help="Registration run folder id. Default: run_0001")
    parser.add_argument("--spacing", type=parse_spacing_override, default=(0.711, 0.711, 0.711), help="Spacing override, for example 0.711. Default: 0.711")
    parser.add_argument("--z-range", type=parse_z_range, default=None, help="Optional 1-based z slice range. Default: full source depth")
    parser.add_argument("--ims-resolution-level", type=int, default=0, help=".ims resolution level. Default: 0")
    parser.add_argument("--ims-timepoint", type=int, default=0, help=".ims timepoint. Default: 0")
    parser.add_argument("--ims-dataset-template", default="", help="Optional .ims dataset template with {level}, {timepoint}, {channel}")
    parser.add_argument("--ims-axis-order", choices=["auto", "zyx", "xyz"], default="auto", help="Default: auto")
    parser.add_argument("--downsample-factor", "--ds", type=float, default=4.0, help="DAPI registration downsample factor. Default: 4")
    parser.add_argument("--type-of-transform", "--tx", default="antsRegistrationSyNQuick[s]", help="ANTsPy DAPI transform preset. Default: antsRegistrationSyNQuick[s]")
    parser.add_argument("--threads", "--t", type=int, default=32, help="CPU threads for ANTs/ITK. Default: 32")
    parser.add_argument("--progress-interval", "--pi", type=float, default=30.0, help="Heartbeat seconds for long stages. Default: 30")
    parser.add_argument("--interpolator", default="linear", help="Interpolator for intensity channels. Default: linear")
    parser.add_argument("--registered-format", "--reg-format", choices=["nrrd", "tif", "both"], default="tif", help="Registered output format. Default: tif")
    parser.add_argument(
        "--tiff-prefix",
        "--tif-prefix",
        default="registered",
        help="Common TIFF filename prefix before Imaris channel marker _C###. Default: registered.",
    )
    parser.add_argument(
        "--channel-offset",
        "--co",
        type=int,
        default=0,
        help="First final _C### channel number for this ROI cycle; selected --ch channels are packed consecutively. Default: 0",
    )
    parser.add_argument("--series-id", default="", help=argparse.SUPPRESS)
    parser.add_argument("--output-pixeltype", choices=["uint16", "float32"], default="uint16", help="Registered output pixel type. Default: uint16")
    parser.add_argument("--uint16-scaling", "--u16", choices=["clip", "minmax", "robust"], default="clip", help="uint16 conversion. Default: clip to preserve intensity range")
    parser.add_argument("--save-raw-roi", action="store_true", help="Also write unregistered extracted ROI NRRDs for debugging. Default: keep ROI crops in memory only")
    parser.add_argument("--cache-fixed-roi", action="store_true", help="When --fixed-crop is a full .ims, also write the extracted fixed ROI under reference/roi for reuse")
    parser.add_argument(
        "--cache-fixed-roi-format",
        choices=["nrrd", "tif", "tiff", "nii", "nii.gz"],
        default="nrrd",
        help="Format for --cache-fixed-roi when extracting the fixed ROI from a full .ims. Default: nrrd",
    )
    parser.add_argument("--apply-only", "--ao", action="store_true", help="Reuse existing DAPI transforms and only extract/apply requested channels")
    parser.add_argument("--no-qc", dest="qc", action="store_false", default=True, help="Disable DAPI center-slice QC PNG output")
    parser.add_argument("--open-qc", action="store_true", help="Open DAPI QC PNG after it is written")
    parser.add_argument("--qc-max-panel-side", type=int, default=900, help="Maximum side of each QC panel. Default: 900")
    parser.add_argument("--verbose", action="store_true", help="Enable ANTs registration diagnostic output")
    return parser


def record_timing(timings: list[dict], stage: str, start: float, **extra) -> float:
    elapsed = time.time() - start
    entry = {"stage": stage, "seconds": elapsed, "duration": format_duration(elapsed)}
    entry.update(extra)
    timings.append(entry)
    print(f"[timing] {stage}: {format_duration(elapsed)}", flush=True)
    return elapsed


def image_suffix(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".nii.gz"):
        return ".nii.gz"
    return path.suffix.lower()


def infer_fixed_crop(project_root: Path) -> Path:
    roi_dir = project_root / "reference" / "roi"
    candidates = []
    for path in sorted(path for pattern in ("*.nrrd", "*.tif", "*.tiff", "*.nii", "*.nii.gz") for path in roi_dir.glob(pattern)):
        name = path.name.lower()
        if "registered" in name or "_roi_raw" in name:
            continue
        candidates.append(path)
    if len(candidates) != 1:
        raise ValueError(
            f"Could not infer one reference fixed crop under {roi_dir}; found {candidates}. "
            "Pass --fixed-crop explicitly."
        )
    return candidates[0].resolve()


def standard_layout(project_root: Path, cycle_id: str, run_id: str) -> dict[str, Path]:
    cycle_root = project_root / "cycles" / cycle_id
    paths = {
        "cycle_root": cycle_root,
        "roi_dir": cycle_root / "roi",
        "transform_dir": cycle_root / "registration" / run_id / "quicksyn_transforms",
        "output_dir": cycle_root / "registration" / run_id / "outputs",
        "qc_dir": cycle_root / "registration" / run_id / "qc",
        "log_dir": cycle_root / "registration" / run_id / "logs",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def place_in_fixed_roi_frame(img: sitk.Image, fixed) -> sitk.Image:
    img.SetSpacing(tuple(float(v) for v in fixed.spacing))
    img.SetOrigin(tuple(float(v) for v in fixed.origin))
    img.SetDirection(tuple(float(v) for v in np.asarray(fixed.direction).ravel()))
    return img


def sitk_to_ants(img: sitk.Image):
    import ants

    return ants.from_numpy(
        sitk.GetArrayFromImage(img).T,
        origin=img.GetOrigin(),
        spacing=img.GetSpacing(),
        direction=np.asarray(img.GetDirection()).reshape((3, 3)),
    )


def channel_source_path(args, channel: int, source_map: dict[int, Path]) -> Path:
    return source_map.get(channel, Path(args.moving_source).expanduser().resolve())


def extract_channel_in_memory(args, channel: int, source_map: dict[int, Path], tile, fixed):
    source_path = channel_source_path(args, channel, source_map)
    spec = build_channel_spec(args, source_path, channel)
    print(
        f"[roi-cycle] extracting channel {channel} directly from {spec.kind} source: {source_path} "
        f"shape_zyx={spec.shape_zyx} spacing_xyz={spec.spacing_xyz}",
        flush=True,
    )
    with progress_heartbeat(f"extract channel {channel} ROI in memory", args.progress_interval):
        sitk_img = read_channel_tile_from_spec(spec, tile)
    sitk_img = place_in_fixed_roi_frame(sitk_img, fixed)
    print_sitk_stats(f"[roi-cycle] extracted channel {channel} ROI", sitk_img)
    return sitk_to_ants(sitk_img), source_path


def load_fixed_roi(args, fixed_crop_path: Path, project_root: Path, progress_interval: float):
    import ants

    suffix = image_suffix(fixed_crop_path)
    if suffix != ".ims":
        with progress_heartbeat("load fixed ROI", progress_interval):
            return ants.image_read(str(fixed_crop_path)), str(fixed_crop_path), None

    fixed_roi_csv = Path(args.fixed_roi_csv).expanduser().resolve() if args.fixed_roi_csv else None
    if fixed_roi_csv is None or not fixed_roi_csv.exists():
        raise ValueError(
            "--fixed-crop points to a full .ims file, so --fixed-roi-csv is required. "
            "Pass the Cycle0/reference XY_Coordinates.csv with --fixed-roi-csv, "
            "or choose a cached fixed crop .nrrd/.tif/.nii instead of the full .ims."
        )

    fixed_spec = build_channel_spec(args, fixed_crop_path, args.fixed_dapi_channel)
    fixed_roi = parse_xy_coordinates_csv(fixed_roi_csv, args.z_range, (0, fixed_spec.shape_zyx[0]))
    fixed_tile = tile_from_roi(fixed_roi)
    print(
        f"[roi-cycle] fixed crop is a full .ims; streaming fixed DAPI ROI from channel {args.fixed_dapi_channel}: "
        f"{fixed_crop_path}",
        flush=True,
    )
    print(f"[roi-cycle] fixed ROI CSV={fixed_roi_csv}", flush=True)
    print(f"[roi-cycle] fixed ROI xyz exclusive={fixed_roi}", flush=True)
    with progress_heartbeat("extract fixed ROI from .ims", progress_interval):
        sitk_img = read_channel_tile_from_spec(fixed_spec, fixed_tile)
    print_sitk_stats("[roi-cycle] fixed ROI extracted from .ims", sitk_img)
    fixed = sitk_to_ants(sitk_img)

    fixed_manifest_path = str(fixed_crop_path)
    if args.cache_fixed_roi:
        reference_roi_dir = project_root / "reference" / "roi"
        reference_roi_dir.mkdir(parents=True, exist_ok=True)
        cache_format = args.cache_fixed_roi_format
        extension = ".nii.gz" if cache_format == "nii.gz" else f".{cache_format}"
        cached_path = reference_roi_dir / f"{fixed_crop_path.stem}_ch{args.fixed_dapi_channel}_roi_from_csv{extension}"
        if not cached_path.exists():
            print(f"[roi-cycle] caching extracted fixed ROI for reuse: {cached_path}", flush=True)
            with progress_heartbeat(f"write cached fixed ROI {cache_format}", progress_interval):
                sitk.WriteImage(sitk_img, str(cached_path))
        else:
            print(f"[roi-cycle] cached fixed ROI already exists: {cached_path}", flush=True)
        fixed_manifest_path = str(cached_path)
    return fixed, fixed_manifest_path, fixed_roi


def output_channel_number(channels: list[int], channel: int, channel_offset: int) -> int:
    try:
        selected_index = channels.index(channel)
    except ValueError as exc:
        raise ValueError(f"Channel {channel} is not included in selected output channels {channels}") from exc
    output_channel = channel_offset + selected_index
    if output_channel < 0:
        raise ValueError(f"Final TIFF channel number must be non-negative, got offset {channel_offset} + index {selected_index}")
    return output_channel


def registered_cycle_output_paths(
    output_dir: Path,
    cycle_id: str,
    channels: list[int],
    channel: int,
    channel_offset: int,
    tiff_prefix: str,
    registered_format: str,
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    if registered_format in {"nrrd", "both"}:
        paths["nrrd"] = output_dir / f"{cycle_id}_ch{channel}_roi_registered.nrrd"
    if registered_format in {"tif", "both"}:
        paths["tif"] = output_dir / f"{tiff_prefix}_C{output_channel_number(channels, channel, channel_offset):03d}.tif"
    return paths


def main() -> int:
    args = build_parser().parse_args()
    if args.channel_offset < 0:
        raise ValueError("--channel-offset must be >= 0")
    if args.dapi_channel not in args.channels:
        raise ValueError("--dapi-channel must be included in --ch")
    if args.series_id:
        print(
            "[roi-cycle] --series-id is retained only for command compatibility and is ignored; "
            "TIFF outputs use Imaris channel tokens such as registered_C000.tif.",
            flush=True,
        )
    configure_ants_runtime(args.threads)
    import ants

    start_total = time.time()
    timings: list[dict] = []
    project_root = Path(args.project_root).expanduser().resolve()
    layout = standard_layout(project_root, args.cycle_id, args.run_id)
    fixed_crop_path = Path(args.fixed_crop).expanduser().resolve() if args.fixed_crop else infer_fixed_crop(project_root)

    input_csv = Path(args.roi_csv).expanduser().resolve()
    if not input_csv.exists():
        raise FileNotFoundError(f"ROI CSV not found: {input_csv}")
    saved_csv = layout["roi_dir"] / "XY_Coordinates.csv"
    if input_csv != saved_csv.resolve():
        shutil.copy2(input_csv, saved_csv)
        print(f"[roi-cycle] copied ROI CSV into project structure: {saved_csv}", flush=True)
    else:
        print(f"[roi-cycle] using ROI CSV already in project structure: {saved_csv}", flush=True)

    source_map = parse_input_ims_map(args.source_map) if args.source_map else {}
    dapi_path = channel_source_path(args, args.dapi_channel, source_map)
    dapi_spec = build_channel_spec(args, dapi_path, args.dapi_channel)
    roi = parse_xy_coordinates_csv(saved_csv, args.z_range, (0, dapi_spec.shape_zyx[0]))
    tile = tile_from_roi(roi)
    print(f"[roi-cycle] cycle={args.cycle_id} run={args.run_id}", flush=True)
    print(f"[roi-cycle] fixed crop={fixed_crop_path}", flush=True)
    print(f"[roi-cycle] moving DAPI source={dapi_path}", flush=True)
    print(f"[roi-cycle] ROI xyz exclusive={roi}; raw ROI intermediates={'enabled' if args.save_raw_roi else 'disabled'}", flush=True)
    output_channel_map = {channel: output_channel_number(args.channels, channel, args.channel_offset) for channel in args.channels}
    print(f"[roi-cycle] output_channel_map={output_channel_map}", flush=True)
    print(f"[roi-cycle] transform dir={layout['transform_dir']}", flush=True)
    print(f"[roi-cycle] output dir={layout['output_dir']}", flush=True)

    stage_start = time.time()
    fixed, fixed_crop_manifest_path, fixed_roi = load_fixed_roi(args, fixed_crop_path, project_root, args.progress_interval)
    describe_image("[roi-cycle] fixed ROI", fixed)
    record_timing(timings, "load/extract fixed ROI", stage_start)

    manifest = None
    dapi_moving = None
    if not args.apply_only:
        stage_start = time.time()
        dapi_moving, dapi_source_path = extract_channel_in_memory(args, args.dapi_channel, source_map, tile, fixed)
        describe_image("[roi-cycle] moving DAPI ROI", dapi_moving)
        record_timing(timings, "extract moving DAPI ROI from source", stage_start)

        stage_start = time.time()
        fixed_ds = downsample_by_spacing(fixed, args.downsample_factor)
        moving_ds = downsample_by_spacing(dapi_moving, args.downsample_factor)
        describe_image("[roi-cycle] fixed downsampled", fixed_ds)
        describe_image("[roi-cycle] moving downsampled", moving_ds)
        record_timing(timings, "downsample DAPI ROIs in memory", stage_start)

        stage_start = time.time()
        outprefix = str(layout["transform_dir"] / "roi_quicksyn_")
        print(f"[roi-cycle] registering DAPI with {args.type_of_transform}; outprefix={outprefix}", flush=True)
        with progress_heartbeat("register DAPI ROI", args.progress_interval):
            tx = ants.registration(
                fixed=fixed_ds,
                moving=moving_ds,
                type_of_transform=args.type_of_transform,
                outprefix=outprefix,
                aff_metric="mattes",
                syn_metric="mattes",
                singleprecision=True,
                verbose=args.verbose,
            )
        manifest = {
            "mode": "direct_csv_roi_from_source",
            "cycle_id": args.cycle_id,
            "fixed_crop": fixed_crop_manifest_path,
            "fixed_source": str(fixed_crop_path),
            "fixed_roi_csv": str(Path(args.fixed_roi_csv).expanduser().resolve()) if args.fixed_roi_csv else "",
            "fixed_roi_xyz_exclusive": fixed_roi,
            "moving_source": str(dapi_source_path),
            "roi_csv": str(saved_csv),
            "roi_xyz_exclusive": roi,
            "dapi_channel": args.dapi_channel,
            "downsample_factor": args.downsample_factor,
            "type_of_transform": args.type_of_transform,
            "fwdtransforms": tx["fwdtransforms"],
            "invtransforms": tx["invtransforms"],
        }
        (layout["transform_dir"] / "fwdtransforms.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        record_timing(timings, "register downsampled DAPI ROI", stage_start)

        stage_start = time.time()
        metrics = compute_downsampled_metrics_from_images(fixed_ds, moving_ds, manifest["fwdtransforms"], args.progress_interval)
        write_metrics_json(layout["transform_dir"] / "roi_similarity_metrics.json", metrics)
        record_timing(timings, "compute DAPI similarity metrics", stage_start)
        del fixed_ds, moving_ds, tx
        gc.collect()
    else:
        manifest = load_transform_manifest(layout["transform_dir"])
        print(f"[roi-cycle] apply-only: reusing transforms from {layout['transform_dir']}", flush=True)

    outputs: list[dict] = []
    for channel in args.channels:
        channel_start = time.time()
        if channel == args.dapi_channel and dapi_moving is not None:
            moving = dapi_moving
            source_path = dapi_path
        else:
            stage_start = time.time()
            moving, source_path = extract_channel_in_memory(args, channel, source_map, tile, fixed)
            record_timing(timings, f"extract channel {channel} ROI from source", stage_start, channel=channel)

        print_ants_stats(f"[roi-cycle] channel {channel} input ROI", moving)
        if args.save_raw_roi:
            raw_path = layout["roi_dir"] / f"{args.cycle_id}_ch{channel}_roi_raw.nrrd"
            with progress_heartbeat(f"write raw channel {channel} ROI", args.progress_interval):
                ants.image_write(moving, str(raw_path))
            raw_path_text = str(raw_path)
        else:
            raw_path_text = ""

        stage_start = time.time()
        with progress_heartbeat(f"apply channel {channel} transform", args.progress_interval):
            registered = ants.apply_transforms(
                fixed=fixed,
                moving=moving,
                transformlist=manifest["fwdtransforms"],
                interpolator=args.interpolator,
                singleprecision=True,
            )
        record_timing(timings, f"apply transform channel {channel}", stage_start, channel=channel)
        print_ants_stats(f"[roi-cycle] registered channel {channel}", registered)

        if args.qc and channel == args.dapi_channel:
            stage_start = time.time()
            export_center_slice_qc(
                fixed,
                moving,
                registered,
                layout["qc_dir"],
                f"{args.cycle_id}_dapi",
                args.open_qc,
                args.qc_max_panel_side,
            )
            record_timing(timings, "write DAPI QC images", stage_start)

        output_paths = registered_cycle_output_paths(
            layout["output_dir"],
            args.cycle_id,
            args.channels,
            channel,
            args.channel_offset,
            args.tiff_prefix,
            args.registered_format,
        )
        stage_start = time.time()
        write_registered_outputs(
            registered,
            output_paths,
            args.output_pixeltype,
            args.progress_interval,
            args.uint16_scaling,
        )
        record_timing(timings, f"write registered channel {channel}", stage_start, channel=channel)
        outputs.append(
            {
                "channel": channel,
                "output_channel": output_channel_map[channel],
                "input_source": str(source_path),
                "raw_roi_written": raw_path_text,
                "registered_nrrd": str(output_paths.get("nrrd", "")),
                "registered_tif": str(output_paths.get("tif", "")),
                "elapsed": format_duration(time.time() - channel_start),
            }
        )
        if channel == args.dapi_channel:
            dapi_moving = None
        del moving, registered
        gc.collect()

    if dapi_moving is not None:
        del dapi_moving
        gc.collect()

    total_seconds = time.time() - start_total
    record_timing(timings, "total", start_total)
    timing_path = layout["log_dir"] / "stage_timing_direct_roi_cycle.json"
    timing_path.write_text(json.dumps({"timings": timings, "total_seconds": total_seconds}, indent=2), encoding="utf-8")
    run_manifest = {
        "cycle_id": args.cycle_id,
        "moving_source": str(Path(args.moving_source).expanduser().resolve()),
        "saved_roi_csv": str(saved_csv),
        "fixed_crop": fixed_crop_manifest_path,
        "fixed_source": str(fixed_crop_path),
        "fixed_roi_csv": str(Path(args.fixed_roi_csv).expanduser().resolve()) if args.fixed_roi_csv else "",
        "fixed_roi_xyz_exclusive": fixed_roi,
        "roi_xyz_exclusive": roi,
        "channel_offset": args.channel_offset,
        "output_channel_map": output_channel_map,
        "raw_roi_intermediates_written": bool(args.save_raw_roi),
        "transform_dir": str(layout["transform_dir"]),
        "timing_log": str(timing_path),
        "total_elapsed": format_duration(total_seconds),
        "outputs": outputs,
    }
    manifest_path = layout["log_dir"] / "direct_roi_cycle_manifest.json"
    manifest_path.write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")
    print(f"[roi-cycle] wrote timing log: {timing_path}", flush=True)
    print(f"[roi-cycle] wrote run manifest: {manifest_path}", flush=True)
    print(f"[roi-cycle] complete in {format_duration(total_seconds)}", flush=True)
    del fixed
    gc.collect()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
