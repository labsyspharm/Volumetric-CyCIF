#!/usr/bin/env python3
# Developed by Alex Wong
# Cite: Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in
# Preprint: https://doi.org/10.64898/2026.05.17.725158
# Registration method: ANTsX/ANTsPy - https://github.com/ANTsX/ANTsPy

"""Register already-cropped full-resolution NRRDs with ANTsPy quick SyN."""

from __future__ import annotations

import argparse
import gc
import json
import math
import os
import threading
import time
import webbrowser
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import SimpleITK as sitk
try:
    from PIL import Image, ImageDraw
except ImportError:  # pragma: no cover - optional runtime dependency
    Image = None
    ImageDraw = None

ants = None


def configure_ants_runtime(threads: int) -> None:
    global ants
    if threads < 0:
        raise ValueError("--threads must be >= 0")
    if threads:
        thread_value = str(int(threads))
        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = thread_value
        os.environ["OMP_NUM_THREADS"] = thread_value
        os.environ["OPENBLAS_NUM_THREADS"] = thread_value
        os.environ["MKL_NUM_THREADS"] = thread_value
        os.environ["NUMEXPR_NUM_THREADS"] = thread_value
        print(f"[runtime] requested {thread_value} CPU thread(s) for ANTs/ITK", flush=True)
    else:
        current = os.environ.get("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS", "ANTs default")
        print(f"[runtime] CPU threads={current}", flush=True)

    import ants as ants_module

    ants = ants_module


def downsample_by_spacing(img, factor: float):
    new_spacing = tuple(float(s) * factor for s in img.spacing)
    return ants.resample_image(img, new_spacing, use_voxels=False, interp_type=0)


def format_duration(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, sec = divmod(int(round(seconds)), 60)
    if minutes < 60:
        return f"{minutes}m{sec:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m{sec:02d}s"


@contextmanager
def progress_heartbeat(label: str, interval_seconds: float):
    if interval_seconds <= 0:
        yield
        return

    stop_event = threading.Event()
    start = time.time()

    def worker() -> None:
        while not stop_event.wait(interval_seconds):
            print(
                f"[{label}] still running; elapsed={format_duration(time.time() - start)}",
                flush=True,
            )

    thread = threading.Thread(target=worker, daemon=True)
    thread.start()
    try:
        yield
    finally:
        stop_event.set()
        thread.join(timeout=1.0)


def describe_image(label: str, img) -> None:
    print(
        f"{label}: shape={tuple(int(v) for v in img.shape)} "
        f"spacing={tuple(float(v) for v in img.spacing)} pixeltype={img.pixeltype}",
        flush=True,
    )


def ants_to_sitk_uint16(img, scaling: str) -> "sitk.Image":
    arr = img.numpy()
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        arr_u16 = np.zeros(arr.shape, dtype=np.uint16)
    elif scaling == "clip":
        arr_u16 = np.rint(np.clip(arr, 0, 65535)).astype(np.uint16, copy=False)
    elif scaling in {"minmax", "robust"}:
        if scaling == "robust":
            lo, hi = np.percentile(finite, (0.1, 99.9))
        else:
            lo, hi = float(np.min(finite)), float(np.max(finite))
        if hi <= lo:
            arr_u16 = np.zeros(arr.shape, dtype=np.uint16)
        else:
            scaled = np.clip((arr - lo) / (hi - lo), 0.0, 1.0)
            arr_u16 = np.rint(scaled * 65535).astype(np.uint16, copy=False)
        print(f"[write] uint16 scaling={scaling} input_range=({lo:.6g}, {hi:.6g})", flush=True)
    else:
        raise ValueError(f"Unsupported uint16 scaling mode: {scaling}")
    sitk_img = sitk.GetImageFromArray(arr_u16.T)
    sitk_img.SetSpacing(tuple(float(v) for v in img.spacing))
    sitk_img.SetOrigin(tuple(float(v) for v in img.origin))
    sitk_img.SetDirection(tuple(float(v) for v in np.asarray(img.direction).ravel()))
    return sitk_img


def write_ants_image(
    img,
    output_path: Path,
    pixeltype: str,
    progress_interval: float,
    uint16_scaling: str = "minmax",
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if pixeltype == "uint16":
        with progress_heartbeat("write uint16 nrrd", progress_interval):
            sitk.WriteImage(ants_to_sitk_uint16(img, uint16_scaling), str(output_path))
    elif pixeltype == "float32":
        with progress_heartbeat("write float32 nrrd", progress_interval):
            ants.image_write(img, str(output_path))
    else:
        raise ValueError(f"Unsupported output pixel type: {pixeltype}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run rigid+affine+quick SyN on downsampled ROI crops, then apply the "
            "resulting transforms to the full-resolution moving ROI crop."
        )
    )
    parser.add_argument(
        "--mode",
        choices=["all", "staged", "memory", "downsample", "register", "apply"],
        default="memory",
        help=(
            "Workflow mode. 'staged' writes downsampled inputs, registers them, releases memory, "
            "then applies transforms to full-res crops. 'memory' keeps downsampled inputs in RAM "
            "and does not write them to disk. Default: memory."
        ),
    )
    parser.add_argument("--fixed-crop", default="", help="Full-resolution fixed ROI crop NRRD")
    parser.add_argument("--moving-crop", default="", help="Full-resolution moving ROI crop NRRD")
    parser.add_argument(
        "--output",
        default="",
        help=(
            "Output registered full-resolution moving ROI NRRD. Default: "
            "<registration-run>/outputs/<moving-crop-stem>_registered.nrrd"
        ),
    )
    parser.add_argument("--fixed-downsampled", default="", help="Downsampled fixed ROI NRRD path")
    parser.add_argument("--moving-downsampled", default="", help="Downsampled moving ROI NRRD path")
    parser.add_argument(
        "--transform-dir",
        required=True,
        help="Directory for ANTs transform files and fwdtransforms.json",
    )
    parser.add_argument(
        "--downsample-factor",
        type=float,
        default=4.0,
        help="Spacing multiplier for registration inputs. Default: 4",
    )
    parser.add_argument(
        "--type-of-transform",
        default="antsRegistrationSyNQuick[s]",
        help="ANTsPy registration transform preset. Default: antsRegistrationSyNQuick[s]",
    )
    parser.add_argument(
        "--interpolator",
        default="linear",
        help="Interpolator for applying transforms to the full-resolution crop. Default: linear",
    )
    parser.add_argument(
        "--output-pixeltype",
        choices=["uint16", "float32"],
        default="uint16",
        help="Pixel type for registered full-resolution NRRD outputs. Default: uint16",
    )
    parser.add_argument(
        "--uint16-scaling",
        choices=["minmax", "robust", "clip"],
        default="minmax",
        help=(
            "How to convert transformed float images to uint16. 'minmax' rescales the full "
            "registered range to 0-65535 and avoids hard clipping; 'robust' uses 0.1-99.9 "
            "percentiles; 'clip' preserves numeric values but clips outside 0-65535. Default: minmax."
        ),
    )
    parser.add_argument(
        "--metrics-output",
        default="",
        help="Optional JSON path for before/after ROI similarity metrics. Default: <transform-dir>/roi_similarity_metrics.json",
    )
    parser.add_argument(
        "--qc-output-dir",
        default="",
        help="Directory for center-slice QC PNGs. Default: <registration-run>/qc",
    )
    parser.add_argument(
        "--no-qc",
        dest="qc",
        action="store_false",
        default=True,
        help="Disable center-slice QC PNG output.",
    )
    parser.add_argument(
        "--open-qc",
        action="store_true",
        help="Open the center-slice QC PNG after writing it.",
    )
    parser.add_argument(
        "--qc-max-panel-side",
        type=int,
        default=900,
        help="Maximum width or height for each QC slice panel. Default: 900",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=0,
        help=(
            "Number of CPU threads for ANTs/ITK. Use 0 to keep the ANTs default or any "
            "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS value already set in the environment. Default: 0."
        ),
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=30.0,
        help=(
            "Print an elapsed-time progress heartbeat every N seconds during long ANTs operations. "
            "Use 0 to disable. Default: 30."
        ),
    )
    parser.add_argument("--verbose", action="store_true", help="Print ANTs registration output")
    return parser


def assert_same_dimension(fixed, moving) -> None:
    if fixed.dimension != moving.dimension:
        raise ValueError(f"Fixed/moving dimension mismatch: {fixed.dimension} vs {moving.dimension}")
    if fixed.dimension != 3:
        raise ValueError(f"Expected 3D ROI crops, got dimension={fixed.dimension}")


def require_path(value: str, flag: str) -> Path:
    if not value:
        raise ValueError(f"{flag} is required for this mode")
    return Path(value).expanduser().resolve()


def default_registered_output(moving_crop_path: Path, transform_dir: Path) -> Path:
    transform_name = transform_dir.name.lower()
    if transform_name == "transforms" or transform_name.endswith("_transforms"):
        run_dir = transform_dir.parent
    else:
        run_dir = transform_dir
    return run_dir / "outputs" / f"{moving_crop_path.stem}_registered.nrrd"


def registration_run_dir(transform_dir: Path) -> Path:
    transform_name = transform_dir.name.lower()
    if transform_name == "transforms" or transform_name.endswith("_transforms"):
        return transform_dir.parent
    return transform_dir


def resolve_output_path(value: str, moving_crop_path: Path, transform_dir: Path) -> Path:
    if value:
        return Path(value).expanduser().resolve()
    output_path = default_registered_output(moving_crop_path, transform_dir).resolve()
    print(f"[output] no --output supplied; using default: {output_path}", flush=True)
    return output_path


def resolve_qc_dir(value: str, transform_dir: Path) -> Path:
    if value:
        return Path(value).expanduser().resolve()
    return (registration_run_dir(transform_dir) / "qc").resolve()


def robust_uint8(slice_2d: np.ndarray) -> np.ndarray:
    arr = np.asarray(slice_2d, dtype=np.float32)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return np.zeros(arr.shape, dtype=np.uint8)
    lo, hi = np.percentile(finite, (0.5, 99.5))
    if hi <= lo:
        hi = float(np.max(finite))
        lo = float(np.min(finite))
    if hi <= lo:
        return np.zeros(arr.shape, dtype=np.uint8)
    arr = np.clip((arr - lo) / (hi - lo), 0.0, 1.0)
    return (arr * 255).astype(np.uint8)


def overlay_magenta_cyan(fixed_slice: np.ndarray, moving_slice: np.ndarray) -> "Image.Image":
    fixed_u8 = robust_uint8(fixed_slice)
    moving_u8 = robust_uint8(moving_slice)
    h = min(fixed_u8.shape[0], moving_u8.shape[0])
    w = min(fixed_u8.shape[1], moving_u8.shape[1])
    fixed_u8 = fixed_u8[:h, :w]
    moving_u8 = moving_u8[:h, :w]
    rgb = np.zeros((h, w, 3), dtype=np.uint8)
    rgb[..., 0] = np.maximum(rgb[..., 0], fixed_u8)
    rgb[..., 2] = np.maximum(rgb[..., 2], fixed_u8)
    rgb[..., 1] = np.maximum(rgb[..., 1], moving_u8)
    rgb[..., 2] = np.maximum(rgb[..., 2], moving_u8)
    return Image.fromarray(rgb)


def center_slice(arr: np.ndarray, axis: int) -> np.ndarray:
    idx = int(arr.shape[axis] // 2)
    if axis == 0:
        slc = arr[idx, :, :]
    elif axis == 1:
        slc = arr[:, idx, :]
    else:
        slc = arr[:, :, idx]
    return np.asarray(slc).T


def fit_panel(img: "Image.Image", max_side: int) -> "Image.Image":
    if max_side <= 0:
        return img
    out = img.copy()
    resampling = getattr(getattr(Image, "Resampling", Image), "BILINEAR")
    out.thumbnail((max_side, max_side), resampling)
    return out


def label_panel(img: "Image.Image", label: str) -> "Image.Image":
    header_h = 28
    out = Image.new("RGB", (img.width, img.height + header_h), (20, 20, 20))
    out.paste(img.convert("RGB"), (0, header_h))
    draw = ImageDraw.Draw(out)
    draw.text((8, 7), label, fill=(240, 240, 240))
    return out


def export_center_slice_qc(
    fixed_img,
    moving_before_img,
    moving_after_img,
    qc_dir: Path,
    moving_label: str,
    open_qc: bool,
    max_panel_side: int,
) -> None:
    if Image is None or ImageDraw is None:
        print("[qc] Pillow is not installed; skipping center-slice QC PNG", flush=True)
        return

    start = time.time()
    qc_dir.mkdir(parents=True, exist_ok=True)
    print("[qc] building XY/XZ/YZ center-slice before/after overlay PNG", flush=True)
    fixed_arr = np.asarray(fixed_img.numpy())
    before_arr = np.asarray(moving_before_img.numpy())
    after_arr = np.asarray(moving_after_img.numpy())

    panels = []
    for plane_label, axis in (("XY", 2), ("XZ", 1), ("YZ", 0)):
        before = overlay_magenta_cyan(
            center_slice(fixed_arr, axis),
            center_slice(before_arr, axis),
        )
        after = overlay_magenta_cyan(
            center_slice(fixed_arr, axis),
            center_slice(after_arr, axis),
        )
        panels.append(
            (
                label_panel(fit_panel(before, max_panel_side), f"{plane_label} center before"),
                label_panel(fit_panel(after, max_panel_side), f"{plane_label} center after"),
            )
        )

    gap = 12
    width = max(left.width + gap + right.width for left, right in panels)
    height = sum(max(left.height, right.height) for left, right in panels) + gap * (len(panels) - 1)
    canvas = Image.new("RGB", (width, height), (10, 10, 10))
    y = 0
    for left, right in panels:
        canvas.paste(left, (0, y))
        canvas.paste(right, (left.width + gap, y))
        y += max(left.height, right.height) + gap

    qc_path = qc_dir / f"{moving_label}_center_slices_overlay.png"
    canvas.save(qc_path)
    print(f"[qc] wrote center-slice QC: {qc_path}", flush=True)
    if open_qc:
        open_path = str(qc_path)
        if os.name == "nt":
            os.startfile(open_path)  # type: ignore[attr-defined]
        else:
            webbrowser.open(qc_path.as_uri())
        print(f"[qc] opened QC image: {qc_path}", flush=True)

    del fixed_arr, before_arr, after_arr, panels, canvas
    gc.collect()
    print(f"[qc] finished center-slice QC in {format_duration(time.time() - start)}", flush=True)


def downsample_file(input_path: Path, output_path: Path, factor: float, progress_interval: float) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    print(f"[downsample] loading full-resolution ROI: {input_path}", flush=True)
    with progress_heartbeat("downsample load full-res", progress_interval):
        img = ants.image_read(str(input_path))
    describe_image("[downsample] full-resolution ROI", img)
    print(f"[downsample] resampling by spacing factor {factor:g}", flush=True)
    with progress_heartbeat("downsample resample", progress_interval):
        ds = downsample_by_spacing(img, factor)
    describe_image("[downsample] downsampled ROI", ds)
    with progress_heartbeat("downsample write", progress_interval):
        ants.image_write(ds, str(output_path))
    print(f"[downsample] wrote downsampled ROI: {output_path}", flush=True)
    del img, ds
    gc.collect()
    print(
        f"[downsample] released full-resolution ROI from memory in {format_duration(time.time() - start)}",
        flush=True,
    )


def downsample_file_to_memory(input_path: Path, factor: float, label: str, progress_interval: float):
    start = time.time()
    print(f"[memory] loading full-resolution {label} ROI: {input_path}", flush=True)
    with progress_heartbeat(f"memory load {label}", progress_interval):
        img = ants.image_read(str(input_path))
    describe_image(f"[memory] {label} full-resolution ROI", img)
    print(f"[memory] downsampling {label} ROI by spacing factor {factor:g}", flush=True)
    with progress_heartbeat(f"memory downsample {label}", progress_interval):
        ds = downsample_by_spacing(img, factor)
    describe_image(f"[memory] {label} downsampled ROI", ds)
    del img
    gc.collect()
    print(
        f"[memory] released full-resolution {label} ROI; kept only downsampled image "
        f"in {format_duration(time.time() - start)}",
        flush=True,
    )
    return ds


def compute_ncc(a: np.ndarray, b: np.ndarray) -> float:
    a = a.astype(np.float32, copy=False)
    b = b.astype(np.float32, copy=False)
    am = float(np.mean(a))
    bm = float(np.mean(b))
    da = a - am
    db = b - bm
    denom = math.sqrt(float(np.sum(da * da)) * float(np.sum(db * db)))
    if denom < 1e-12:
        return 0.0
    return float(np.sum(da * db) / denom)


def compute_nrmse(a: np.ndarray, b: np.ndarray) -> float:
    a = a.astype(np.float32, copy=False)
    b = b.astype(np.float32, copy=False)
    rmse = float(np.sqrt(np.mean((a - b) ** 2)))
    dynamic = float(np.max(a) - np.min(a))
    return float(rmse / max(dynamic, 1e-6))


def compute_nmi(a: np.ndarray, b: np.ndarray, bins: int = 64) -> float:
    a = a.astype(np.float32, copy=False).ravel()
    b = b.astype(np.float32, copy=False).ravel()
    hist2d, _, _ = np.histogram2d(a, b, bins=bins)
    total = float(np.sum(hist2d))
    if total <= 0:
        return 0.0

    pxy = hist2d / total
    px = np.sum(pxy, axis=1)
    py = np.sum(pxy, axis=0)
    eps = 1e-12
    hx = -float(np.sum(px * np.log(px + eps)))
    hy = -float(np.sum(py * np.log(py + eps)))
    hxy = -float(np.sum(pxy * np.log(pxy + eps)))
    if hxy <= eps:
        return 0.0
    return float((hx + hy) / hxy)


def compute_similarity_metrics_from_images(fixed_img, moving_img) -> dict:
    fixed_arr = fixed_img.numpy().astype(np.float32, copy=False)
    moving_arr = moving_img.numpy().astype(np.float32, copy=False)
    min_shape = tuple(min(a, b) for a, b in zip(fixed_arr.shape, moving_arr.shape))
    slices = tuple(slice(0, s) for s in min_shape)
    fa = fixed_arr[slices]
    ma = moving_arr[slices]
    metrics = {
        "ncc": compute_ncc(fa, ma),
        "nmi": compute_nmi(fa, ma, bins=64),
        "nrmse": compute_nrmse(fa, ma),
    }
    del fixed_arr, moving_arr, fa, ma
    gc.collect()
    return metrics


def compute_downsampled_before_after_metrics(
    fixed_downsampled_path: Path,
    moving_downsampled_path: Path,
    transform_dir: Path,
) -> dict:
    manifest = load_transform_manifest(transform_dir)
    fixed_ds = ants.image_read(str(fixed_downsampled_path))
    moving_ds = ants.image_read(str(moving_downsampled_path))
    assert_same_dimension(fixed_ds, moving_ds)
    warped_ds = ants.apply_transforms(
        fixed=fixed_ds,
        moving=moving_ds,
        transformlist=manifest["fwdtransforms"],
        interpolator="linear",
        singleprecision=True,
    )
    metrics = {
        "computed_on": "downsampled_roi",
        "fixed_downsampled": str(fixed_downsampled_path),
        "moving_downsampled": str(moving_downsampled_path),
        "interpretation": "higher ncc/nmi is better; lower nrmse is better",
        "before": compute_similarity_metrics_from_images(fixed_ds, moving_ds),
        "after": compute_similarity_metrics_from_images(fixed_ds, warped_ds),
    }
    del fixed_ds, moving_ds, warped_ds
    gc.collect()
    return metrics


def write_metrics_json(metrics_output_path: Path, metrics: dict) -> None:
    metrics_output_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_output_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print("Similarity metrics (higher NCC/NMI better, lower NRMSE better):", flush=True)
    for label in ("before", "after"):
        vals = metrics[label]
        print(
            f"  {label:6s} NCC={vals['ncc']:.6f} NMI={vals['nmi']:.6f} NRMSE={vals['nrmse']:.6f}",
            flush=True,
        )
    print(f"Wrote similarity metrics: {metrics_output_path}", flush=True)


def register_downsampled(
    fixed_downsampled_path: Path,
    moving_downsampled_path: Path,
    transform_dir: Path,
    type_of_transform: str,
    verbose: bool,
    manifest_updates: dict,
    progress_interval: float,
) -> dict:
    start = time.time()
    print(f"[register] loading downsampled fixed ROI: {fixed_downsampled_path}", flush=True)
    with progress_heartbeat("register load fixed downsampled", progress_interval):
        fixed_ds = ants.image_read(str(fixed_downsampled_path))
    describe_image("[register] fixed downsampled", fixed_ds)
    print(f"[register] loading downsampled moving ROI: {moving_downsampled_path}", flush=True)
    with progress_heartbeat("register load moving downsampled", progress_interval):
        moving_ds = ants.image_read(str(moving_downsampled_path))
    describe_image("[register] moving downsampled", moving_ds)
    assert_same_dimension(fixed_ds, moving_ds)

    outprefix = str(transform_dir / "roi_quicksyn_")
    print(f"[register] starting ANTs {type_of_transform}; outprefix={outprefix}", flush=True)
    with progress_heartbeat("register ANTs", progress_interval):
        tx = ants.registration(
            fixed=fixed_ds,
            moving=moving_ds,
            type_of_transform=type_of_transform,
            outprefix=outprefix,
            aff_metric="mattes",
            syn_metric="mattes",
            singleprecision=True,
            verbose=verbose,
        )

    manifest = {
        **manifest_updates,
        "fixed_downsampled": str(fixed_downsampled_path),
        "moving_downsampled": str(moving_downsampled_path),
        "type_of_transform": type_of_transform,
        "fwdtransforms": tx["fwdtransforms"],
        "invtransforms": tx["invtransforms"],
    }
    (transform_dir / "fwdtransforms.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[register] wrote transforms: {transform_dir}", flush=True)
    del fixed_ds, moving_ds, tx
    gc.collect()
    print(f"[register] released downsampled registration images in {format_duration(time.time() - start)}", flush=True)
    return manifest


def register_downsampled_images(
    fixed_ds,
    moving_ds,
    transform_dir: Path,
    type_of_transform: str,
    verbose: bool,
    manifest_updates: dict,
    progress_interval: float,
) -> dict:
    start = time.time()
    describe_image("[register] fixed downsampled", fixed_ds)
    describe_image("[register] moving downsampled", moving_ds)
    assert_same_dimension(fixed_ds, moving_ds)

    outprefix = str(transform_dir / "roi_quicksyn_")
    print(f"[register] starting ANTs {type_of_transform}; outprefix={outprefix}", flush=True)
    with progress_heartbeat("register ANTs", progress_interval):
        tx = ants.registration(
            fixed=fixed_ds,
            moving=moving_ds,
            type_of_transform=type_of_transform,
            outprefix=outprefix,
            aff_metric="mattes",
            syn_metric="mattes",
            singleprecision=True,
            verbose=verbose,
        )

    manifest = {
        **manifest_updates,
        "fixed_downsampled": "in_memory",
        "moving_downsampled": "in_memory",
        "type_of_transform": type_of_transform,
        "fwdtransforms": tx["fwdtransforms"],
        "invtransforms": tx["invtransforms"],
    }
    (transform_dir / "fwdtransforms.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(
        f"[register] wrote transforms after {format_duration(time.time() - start)}: {transform_dir}",
        flush=True,
    )
    del tx
    gc.collect()
    return manifest


def compute_downsampled_metrics_from_images(
    fixed_ds,
    moving_ds,
    transformlist: list[str],
    progress_interval: float,
) -> dict:
    print("[metrics] computing before/after metrics on in-memory downsampled ROIs", flush=True)
    with progress_heartbeat("metrics warp downsampled", progress_interval):
        warped_ds = ants.apply_transforms(
            fixed=fixed_ds,
            moving=moving_ds,
            transformlist=transformlist,
            interpolator="linear",
            singleprecision=True,
        )
    metrics = {
        "computed_on": "downsampled_roi",
        "fixed_downsampled": "in_memory",
        "moving_downsampled": "in_memory",
        "interpretation": "higher ncc/nmi is better; lower nrmse is better",
        "before": compute_similarity_metrics_from_images(fixed_ds, moving_ds),
        "after": compute_similarity_metrics_from_images(fixed_ds, warped_ds),
    }
    del warped_ds
    gc.collect()
    return metrics


def load_transform_manifest(transform_dir: Path) -> dict:
    manifest_path = transform_dir / "fwdtransforms.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Transform manifest not found: {manifest_path}")
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def apply_to_fullres(
    fixed_crop_path: Path,
    moving_crop_path: Path,
    output_path: Path,
    transform_dir: Path,
    interpolator: str,
    output_pixeltype: str,
    uint16_scaling: str,
    qc_dir,
    open_qc: bool,
    qc_max_panel_side: int,
    progress_interval: float,
) -> None:
    start = time.time()
    manifest = load_transform_manifest(transform_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[apply] loading full-resolution fixed reference ROI: {fixed_crop_path}", flush=True)
    with progress_heartbeat("apply load fixed full-res", progress_interval):
        fixed_crop = ants.image_read(str(fixed_crop_path))
    describe_image("[apply] fixed full-resolution ROI", fixed_crop)
    print(f"[apply] loading full-resolution moving ROI: {moving_crop_path}", flush=True)
    with progress_heartbeat("apply load moving full-res", progress_interval):
        moving_crop = ants.image_read(str(moving_crop_path))
    describe_image("[apply] moving full-resolution ROI", moving_crop)
    assert_same_dimension(fixed_crop, moving_crop)
    print(f"[apply] applying {len(manifest['fwdtransforms'])} transform(s)", flush=True)
    with progress_heartbeat("apply transforms full-res", progress_interval):
        warped_full_crop = ants.apply_transforms(
            fixed=fixed_crop,
            moving=moving_crop,
            transformlist=manifest["fwdtransforms"],
            interpolator=interpolator,
            singleprecision=True,
        )
    describe_image("[apply] registered full-resolution ROI", warped_full_crop)
    if qc_dir is not None:
        export_center_slice_qc(
            fixed_crop,
            moving_crop,
            warped_full_crop,
            qc_dir,
            moving_crop_path.stem,
            open_qc,
            qc_max_panel_side,
        )
    write_ants_image(warped_full_crop, output_path, output_pixeltype, progress_interval, uint16_scaling)
    print(f"[apply] wrote registered full-resolution ROI: {output_path}", flush=True)
    del fixed_crop, moving_crop, warped_full_crop
    gc.collect()
    print(f"[apply] released full-resolution apply images in {format_duration(time.time() - start)}", flush=True)


def maybe_write_metrics(
    fixed_downsampled_path: Path,
    moving_downsampled_path: Path,
    transform_dir: Path,
    metrics_output_path: Path,
    allow_missing: bool = False,
) -> None:
    if allow_missing and (not fixed_downsampled_path.exists() or not moving_downsampled_path.exists()):
        print(
            "[metrics] skipping downsampled metrics because staged downsampled ROI files do not exist",
            flush=True,
        )
        return
    print("[metrics] computing before/after metrics on downsampled ROIs", flush=True)
    metrics = compute_downsampled_before_after_metrics(
        fixed_downsampled_path,
        moving_downsampled_path,
        transform_dir,
    )
    write_metrics_json(metrics_output_path, metrics)


def main() -> int:
    args = build_parser().parse_args()
    configure_ants_runtime(args.threads)
    print(
        f"[start] mode={args.mode} downsample_factor={args.downsample_factor:g} "
        f"transform={args.type_of_transform}",
        flush=True,
    )

    transform_dir = Path(args.transform_dir).expanduser().resolve()
    transform_dir.mkdir(parents=True, exist_ok=True)
    metrics_output_path = (
        Path(args.metrics_output).expanduser().resolve()
        if args.metrics_output
        else transform_dir / "roi_similarity_metrics.json"
    )
    qc_dir = resolve_qc_dir(args.qc_output_dir, transform_dir) if args.qc else None
    if qc_dir is not None:
        print(f"[qc] center-slice QC output dir: {qc_dir}", flush=True)

    fixed_ds_path = (
        Path(args.fixed_downsampled).expanduser().resolve()
        if args.fixed_downsampled
        else transform_dir / "fixed_roi_downsampled.nrrd"
    )
    moving_ds_path = (
        Path(args.moving_downsampled).expanduser().resolve()
        if args.moving_downsampled
        else transform_dir / "moving_roi_downsampled.nrrd"
    )

    if args.mode in {"all", "staged", "memory", "downsample"}:
        fixed_crop_path = require_path(args.fixed_crop, "--fixed-crop")
        moving_crop_path = require_path(args.moving_crop, "--moving-crop")
        if args.mode == "memory":
            output_path = resolve_output_path(args.output, moving_crop_path, transform_dir)
            fixed_ds = downsample_file_to_memory(
                fixed_crop_path,
                args.downsample_factor,
                "fixed",
                args.progress_interval,
            )
            moving_ds = downsample_file_to_memory(
                moving_crop_path,
                args.downsample_factor,
                "moving",
                args.progress_interval,
            )
            manifest = register_downsampled_images(
                fixed_ds,
                moving_ds,
                transform_dir,
                args.type_of_transform,
                args.verbose,
                {
                    "mode": args.mode,
                    "fixed_crop": str(fixed_crop_path),
                    "moving_crop": str(moving_crop_path),
                    "output": str(output_path),
                    "downsample_factor": args.downsample_factor,
                    "threads": args.threads,
                    "output_pixeltype": args.output_pixeltype,
                    "uint16_scaling": args.uint16_scaling,
                    "qc_dir": str(qc_dir) if qc_dir is not None else "",
                    "progress_interval": args.progress_interval,
                },
                args.progress_interval,
            )
            metrics = compute_downsampled_metrics_from_images(
                fixed_ds,
                moving_ds,
                manifest["fwdtransforms"],
                args.progress_interval,
            )
            write_metrics_json(metrics_output_path, metrics)
            del fixed_ds, moving_ds
            gc.collect()
            print("[memory] released downsampled registration images before full-resolution apply", flush=True)
            apply_to_fullres(
                fixed_crop_path,
                moving_crop_path,
                output_path,
                transform_dir,
                args.interpolator,
                args.output_pixeltype,
                args.uint16_scaling,
                qc_dir,
                args.open_qc,
                args.qc_max_panel_side,
                args.progress_interval,
            )
            return 0

        if args.mode == "all":
            output_path = resolve_output_path(args.output, moving_crop_path, transform_dir)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with progress_heartbeat("all load fixed full-res", args.progress_interval):
                fixed_crop = ants.image_read(str(fixed_crop_path))
            with progress_heartbeat("all load moving full-res", args.progress_interval):
                moving_crop = ants.image_read(str(moving_crop_path))
            assert_same_dimension(fixed_crop, moving_crop)
            with progress_heartbeat("all downsample fixed", args.progress_interval):
                fixed_ds = downsample_by_spacing(fixed_crop, args.downsample_factor)
            with progress_heartbeat("all downsample moving", args.progress_interval):
                moving_ds = downsample_by_spacing(moving_crop, args.downsample_factor)
            outprefix = str(transform_dir / "roi_quicksyn_")
            with progress_heartbeat("all register ANTs", args.progress_interval):
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
                "mode": args.mode,
                "fixed_crop": str(fixed_crop_path),
                "moving_crop": str(moving_crop_path),
                "output": str(output_path),
                "downsample_factor": args.downsample_factor,
                "threads": args.threads,
                "output_pixeltype": args.output_pixeltype,
                "uint16_scaling": args.uint16_scaling,
                "qc_dir": str(qc_dir) if qc_dir is not None else "",
                "progress_interval": args.progress_interval,
                "type_of_transform": args.type_of_transform,
                "fwdtransforms": tx["fwdtransforms"],
                "invtransforms": tx["invtransforms"],
            }
            (transform_dir / "fwdtransforms.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            with progress_heartbeat("all apply transforms full-res", args.progress_interval):
                warped_full_crop = ants.apply_transforms(
                    fixed=fixed_crop,
                    moving=moving_crop,
                    transformlist=tx["fwdtransforms"],
                    interpolator=args.interpolator,
                    singleprecision=True,
                )
            if qc_dir is not None:
                export_center_slice_qc(
                    fixed_crop,
                    moving_crop,
                    warped_full_crop,
                    qc_dir,
                    moving_crop_path.stem,
                    args.open_qc,
                    args.qc_max_panel_side,
                )
            write_ants_image(
                warped_full_crop,
                output_path,
                args.output_pixeltype,
                args.progress_interval,
                args.uint16_scaling,
            )
            with progress_heartbeat("all warp downsampled metrics", args.progress_interval):
                warped_ds = ants.apply_transforms(
                    fixed=fixed_ds,
                    moving=moving_ds,
                    transformlist=tx["fwdtransforms"],
                    interpolator="linear",
                    singleprecision=True,
                )
            metrics = {
                "computed_on": "downsampled_roi",
                "fixed_downsampled": "in_memory",
                "moving_downsampled": "in_memory",
                "interpretation": "higher ncc/nmi is better; lower nrmse is better",
                "before": compute_similarity_metrics_from_images(fixed_ds, moving_ds),
                "after": compute_similarity_metrics_from_images(fixed_ds, warped_ds),
            }
            write_metrics_json(metrics_output_path, metrics)
            print(f"Wrote registered full-resolution ROI: {output_path}", flush=True)
            print(f"Wrote transforms: {transform_dir}", flush=True)
            del warped_ds
            return 0

        downsample_file(fixed_crop_path, fixed_ds_path, args.downsample_factor, args.progress_interval)
        downsample_file(moving_crop_path, moving_ds_path, args.downsample_factor, args.progress_interval)
        if args.mode == "downsample":
            return 0

    if args.mode in {"staged", "register"}:
        manifest_updates = {
            "mode": args.mode,
            "downsample_factor": args.downsample_factor,
            "threads": args.threads,
            "output_pixeltype": args.output_pixeltype,
            "uint16_scaling": args.uint16_scaling,
            "qc_dir": str(qc_dir) if qc_dir is not None else "",
            "progress_interval": args.progress_interval,
        }
        if args.fixed_crop:
            manifest_updates["fixed_crop"] = str(Path(args.fixed_crop).expanduser().resolve())
        if args.moving_crop:
            manifest_updates["moving_crop"] = str(Path(args.moving_crop).expanduser().resolve())
        register_downsampled(
            fixed_ds_path,
            moving_ds_path,
            transform_dir,
            args.type_of_transform,
            args.verbose,
            manifest_updates,
            args.progress_interval,
        )
        if args.mode == "register":
            maybe_write_metrics(fixed_ds_path, moving_ds_path, transform_dir, metrics_output_path)
            return 0

    if args.mode in {"staged", "apply"}:
        fixed_crop_path = require_path(args.fixed_crop, "--fixed-crop")
        moving_crop_path = require_path(args.moving_crop, "--moving-crop")
        output_path = resolve_output_path(args.output, moving_crop_path, transform_dir)
        apply_to_fullres(
            fixed_crop_path,
            moving_crop_path,
            output_path,
            transform_dir,
            args.interpolator,
            args.output_pixeltype,
            args.uint16_scaling,
            qc_dir,
            args.open_qc,
            args.qc_max_panel_side,
            args.progress_interval,
        )
        maybe_write_metrics(
            fixed_ds_path,
            moving_ds_path,
            transform_dir,
            metrics_output_path,
            allow_missing=args.mode == "apply",
        )
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
