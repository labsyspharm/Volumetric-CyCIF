#!/usr/bin/env python3
"""Shared I/O and tiling utilities for CyANTs microscopy workflows."""

from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import SimpleITK as sitk

try:
    import h5py
except ImportError:  # Optional until an Imaris .ims source is read.
    h5py = None


@dataclass(frozen=True)
class Tile:
    z0: int
    z1: int
    y0: int
    y1: int
    x0: int
    x1: int

    @property
    def shape(self) -> Tuple[int, int, int]:
        return (self.z1 - self.z0, self.y1 - self.y0, self.x1 - self.x0)


@dataclass(frozen=True)
class ImsVolumeSpec:
    path: Path
    dataset_key: str
    axis_order: str
    shape_zyx: Tuple[int, int, int]
    spacing_xyz: Tuple[float, float, float]


def path_suffix(path: Path) -> str:
    lower = path.name.lower()
    if lower.endswith(".nii.gz"):
        return ".nii.gz"
    return path.suffix.lower()


def stem_for_output(path: Path) -> str:
    suffix = path_suffix(path)
    name = path.name
    if suffix and name.lower().endswith(suffix):
        return name[: -len(suffix)]
    return path.stem


def parse_spacing_override(value: str) -> Optional[Tuple[float, float, float]]:
    value = value.strip()
    if not value:
        return None
    parts = [p.strip() for p in value.split(",") if p.strip()]
    try:
        if len(parts) == 1:
            spacing = float(parts[0])
            if spacing <= 0:
                raise ValueError
            return (spacing, spacing, spacing)
        if len(parts) == 3:
            spacing_xyz = tuple(float(p) for p in parts)
            if any(v <= 0 for v in spacing_xyz):
                raise ValueError
            return spacing_xyz
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--spacing must be one positive float or x,y,z positive floats") from exc
    raise argparse.ArgumentTypeError("--spacing must be one float or x,y,z")


def compute_starts(length: int, tile: int, overlap: int) -> List[int]:
    if tile >= length:
        return [0]
    stride = tile - overlap
    if stride <= 0:
        raise ValueError("tile size must be greater than overlap for each axis")
    starts = list(range(0, max(length - tile + 1, 1), stride))
    last = length - tile
    if starts[-1] != last:
        starts.append(last)
    return starts


def generate_tiles(
    volume_shape_zyx: Tuple[int, int, int],
    tile_size_zyx: Tuple[int, int, int],
    overlap_zyx: Tuple[int, int, int],
) -> List[Tile]:
    z, y, x = volume_shape_zyx
    tz, ty, tx = tile_size_zyx
    oz, oy, ox = overlap_zyx
    z_starts = compute_starts(z, min(tz, z), min(oz, max(z - 1, 0)))
    y_starts = compute_starts(y, min(ty, y), min(oy, max(y - 1, 0)))
    x_starts = compute_starts(x, min(tx, x), min(ox, max(x - 1, 0)))
    return [
        Tile(
            z0=z0,
            z1=min(z0 + tz, z),
            y0=y0,
            y1=min(y0 + ty, y),
            x0=x0,
            x1=min(x0 + tx, x),
        )
        for z0 in z_starts
        for y0 in y_starts
        for x0 in x_starts
    ]


def raised_cosine_weights(shape_zyx: Tuple[int, int, int], overlap_zyx: Tuple[int, int, int]) -> np.ndarray:
    def axis_weights(length: int, overlap: int) -> np.ndarray:
        if overlap <= 0 or length <= 1:
            return np.ones(length, dtype=np.float32)
        edge = min(overlap, length // 2)
        if edge <= 0:
            return np.ones(length, dtype=np.float32)
        weights = np.ones(length, dtype=np.float32)
        ramp = np.linspace(0.0, np.pi, edge, endpoint=False, dtype=np.float32)
        edge_weights = 0.5 * (1.0 - np.cos(ramp))
        weights[:edge] = edge_weights
        weights[-edge:] = edge_weights[::-1]
        return np.clip(weights, 1e-4, 1.0)

    wz = axis_weights(shape_zyx[0], overlap_zyx[0])[:, None, None]
    wy = axis_weights(shape_zyx[1], overlap_zyx[1])[None, :, None]
    wx = axis_weights(shape_zyx[2], overlap_zyx[2])[None, None, :]
    return (wz * wy * wx).astype(np.float32)


def _extract_ims_spacing_xyz(h5f: "h5py.File", size_zyx: Tuple[int, int, int]) -> Tuple[float, float, float]:
    attrs = h5f.attrs
    needed = {"ExtMin0", "ExtMin1", "ExtMin2", "ExtMax0", "ExtMax1", "ExtMax2"}
    if not needed.issubset(set(attrs.keys())):
        return (1.0, 1.0, 1.0)
    ext_min = [float(attrs[f"ExtMin{i}"]) for i in range(3)]
    ext_max = [float(attrs[f"ExtMax{i}"]) for i in range(3)]
    size_xyz = (size_zyx[2], size_zyx[1], size_zyx[0])
    spacing = []
    for axis, size in enumerate(size_xyz):
        spacing.append(1.0 if size <= 1 else abs(ext_max[axis] - ext_min[axis]) / float(size - 1))
    return (float(spacing[0]), float(spacing[1]), float(spacing[2]))


def _discover_ims_dataset(
    h5f: "h5py.File", res_level: int, timepoint: int, channel: int
) -> "h5py.Dataset":
    candidates = [
        f"DataSet/ResolutionLevel {res_level}/TimePoint {timepoint}/Channel {channel}/Data",
        f"/DataSet/ResolutionLevel {res_level}/TimePoint {timepoint}/Channel {channel}/Data",
    ]
    for candidate in candidates:
        if candidate in h5f and isinstance(h5f[candidate], h5py.Dataset):
            return h5f[candidate]

    found: List[h5py.Dataset] = []

    def visit(_: str, obj: object) -> None:
        if isinstance(obj, h5py.Dataset) and obj.ndim >= 3:
            found.append(obj)

    h5f.visititems(visit)
    if not found:
        raise ValueError("No 3D datasets found in .ims file")
    return sorted(
        found,
        key=lambda dataset: (
            "Data" not in dataset.name,
            f"ResolutionLevel {res_level}" not in dataset.name,
            f"TimePoint {timepoint}" not in dataset.name,
            f"Channel {channel}" not in dataset.name,
            dataset.size,
        ),
    )[0]


def _shape_to_zyx(dataset_shape: Tuple[int, int, int], axis_order: str) -> Tuple[int, int, int]:
    if axis_order == "zyx":
        return dataset_shape
    if axis_order == "xyz":
        return (dataset_shape[2], dataset_shape[1], dataset_shape[0])
    raise ValueError(f"Unsupported --ims-axis-order value: {axis_order}")


def _parse_ims_dataset_indices(dataset_key: str) -> Optional[Tuple[int, int, int]]:
    match = re.search(
        r"ResolutionLevel\s+(\d+).*?TimePoint\s+(\d+).*?Channel\s+(\d+).*?/Data$",
        dataset_key.replace("\\", "/"),
    )
    if not match:
        return None
    return (int(match.group(1)), int(match.group(2)), int(match.group(3)))


def _available_ims_resolution_levels(h5f: "h5py.File") -> List[int]:
    dataset_group = h5f.get("DataSet")
    if dataset_group is None or not hasattr(dataset_group, "keys"):
        return []
    levels: List[int] = []
    for key in dataset_group.keys():
        match = re.fullmatch(r"ResolutionLevel\s+(\d+)", str(key))
        if match:
            levels.append(int(match.group(1)))
    return sorted(levels)


def available_ims_resolution_levels(path: Path) -> List[int]:
    """Return available Imaris pyramid resolution level indices for an .ims file."""
    if h5py is None:
        raise RuntimeError("Reading .ims requires h5py. Install with: pip install h5py")
    with h5py.File(str(path), "r") as h5f:
        return _available_ims_resolution_levels(h5f)


def _stride_expected_shape(shape_zyx: Tuple[int, int, int], stride_zyx: Tuple[int, int, int]) -> Tuple[int, int, int]:
    return tuple((size + max(1, int(stride)) - 1) // max(1, int(stride)) for size, stride in zip(shape_zyx, stride_zyx))


def _matching_pyramid_dataset(
    h5f: "h5py.File",
    spec: ImsVolumeSpec,
    stride_zyx: Tuple[int, int, int],
) -> Optional[Tuple["h5py.Dataset", Tuple[int, int, int], Tuple[float, float, float], int]]:
    if os.environ.get("CYANTS_DISABLE_IMS_PYRAMID", "").strip() in {"1", "true", "TRUE", "yes", "YES"}:
        return None
    if all(max(1, int(stride)) == 1 for stride in stride_zyx):
        return None
    parsed = _parse_ims_dataset_indices(spec.dataset_key)
    if parsed is None:
        return None
    current_level, timepoint, channel = parsed
    expected_shape = _stride_expected_shape(spec.shape_zyx, stride_zyx)
    candidates = []
    for level in _available_ims_resolution_levels(h5f):
        if level == current_level:
            continue
        try:
            dataset = _discover_ims_dataset(h5f, res_level=level, timepoint=timepoint, channel=channel)
        except Exception:
            continue
        dataset_shape = tuple(int(size) for size in dataset.shape)
        if len(dataset_shape) != 3:
            continue
        shape_zyx = _shape_to_zyx(dataset_shape, spec.axis_order)
        shape_delta = tuple(abs(actual - expected) for actual, expected in zip(shape_zyx, expected_shape))
        if any(delta > 1 for delta in shape_delta):
            continue
        effective_scale_zyx = tuple(
            float(base) / float(max(1, reduced)) for base, reduced in zip(spec.shape_zyx, shape_zyx)
        )
        stride_error = sum(abs(effective - max(1, int(stride))) for effective, stride in zip(effective_scale_zyx, stride_zyx))
        candidates.append((sum(shape_delta), stride_error, level, dataset, shape_zyx, effective_scale_zyx))
    if not candidates:
        return None
    _, _, level, dataset, shape_zyx, effective_scale_zyx = min(candidates, key=lambda item: (item[0], item[1], item[2]))
    return dataset, shape_zyx, effective_scale_zyx, level


def _ims_dataset_to_zyx_array(dataset: "h5py.Dataset", axis_order: str) -> np.ndarray:
    if axis_order == "zyx":
        return dataset[:, :, :]
    if axis_order == "xyz":
        return np.transpose(dataset[:, :, :], (2, 1, 0))
    raise ValueError(f"Unsupported axis_order in spec: {axis_order}")


def _pyramid_spacing_xyz(
    spec: ImsVolumeSpec, effective_scale_zyx: Tuple[float, float, float]
) -> Tuple[float, float, float]:
    return (
        spec.spacing_xyz[0] * effective_scale_zyx[2],
        spec.spacing_xyz[1] * effective_scale_zyx[1],
        spec.spacing_xyz[2] * effective_scale_zyx[0],
    )


def read_ims_lowest_resolution_from_spec(
    spec: ImsVolumeSpec,
) -> Tuple[sitk.Image, int, Tuple[float, float, float]]:
    """Read the coarsest available Imaris pyramid level for fast previews.

    Returns the image, the selected resolution level, and the effective z,y,x
    scale relative to the supplied full-resolution spec.
    """
    if h5py is None:
        raise RuntimeError("Reading .ims requires h5py. Install with: pip install h5py")
    parsed = _parse_ims_dataset_indices(spec.dataset_key)
    if parsed is None:
        raise ValueError(f"Cannot infer .ims resolution/time/channel from dataset path: {spec.dataset_key}")
    current_level, timepoint, channel = parsed
    with h5py.File(str(spec.path), "r") as h5f:
        candidates = []
        for level in _available_ims_resolution_levels(h5f):
            if level == current_level:
                continue
            try:
                dataset = _discover_ims_dataset(h5f, res_level=level, timepoint=timepoint, channel=channel)
            except Exception:
                continue
            dataset_shape = tuple(int(size) for size in dataset.shape)
            if len(dataset_shape) != 3:
                continue
            shape_zyx = _shape_to_zyx(dataset_shape, spec.axis_order)
            effective_scale_zyx = tuple(
                float(base) / float(max(1, reduced)) for base, reduced in zip(spec.shape_zyx, shape_zyx)
            )
            candidates.append((dataset.size, -level, level, dataset, shape_zyx, effective_scale_zyx))
        if not candidates:
            raise ValueError(f"No lower .ims pyramid levels were found for {spec.path}")
        _size, _neg_level, level, dataset, shape_zyx, effective_scale_zyx = min(candidates, key=lambda item: (item[0], item[1]))
        arr = _ims_dataset_to_zyx_array(dataset, spec.axis_order)
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(_pyramid_spacing_xyz(spec, effective_scale_zyx))
    print(
        f"[ims-pyramid] using coarsest ResolutionLevel {level} shape_zyx={shape_zyx} "
        f"effective_scale_zyx=({effective_scale_zyx[0]:.3g}, {effective_scale_zyx[1]:.3g}, "
        f"{effective_scale_zyx[2]:.3g}) for fast preview: {spec.path}",
        flush=True,
    )
    return img, level, effective_scale_zyx


def build_ims_volume_spec(
    path: Path,
    res_level: int,
    timepoint: int,
    channel: int,
    dataset_path: str,
    axis_order: str,
) -> ImsVolumeSpec:
    if h5py is None:
        raise RuntimeError("Reading .ims requires h5py. Install with: pip install h5py")
    with h5py.File(str(path), "r") as h5f:
        if dataset_path:
            if dataset_path not in h5f:
                raise ValueError(f"Dataset path not found in .ims: {dataset_path}")
            dataset = h5f[dataset_path]
            if not isinstance(dataset, h5py.Dataset):
                raise ValueError(f"Path exists but is not a dataset: {dataset_path}")
        else:
            dataset = _discover_ims_dataset(h5f, res_level=res_level, timepoint=timepoint, channel=channel)
        dataset_shape = tuple(int(size) for size in dataset.shape)
        if len(dataset_shape) != 3:
            raise ValueError(f"Expected 3D dataset in .ims, got shape={dataset_shape}")
        resolved_axis = "zyx" if axis_order == "auto" else axis_order
        shape_zyx = _shape_to_zyx(dataset_shape, resolved_axis)
        return ImsVolumeSpec(
            path=path,
            dataset_key=dataset.name,
            axis_order=resolved_axis,
            shape_zyx=shape_zyx,
            spacing_xyz=_extract_ims_spacing_xyz(h5f, size_zyx=shape_zyx),
        )


def read_ims_tile_from_spec(spec: ImsVolumeSpec, tile: Tile) -> sitk.Image:
    if h5py is None:
        raise RuntimeError("Reading .ims requires h5py. Install with: pip install h5py")
    with h5py.File(str(spec.path), "r") as h5f:
        dataset = h5f[spec.dataset_key]
        if not isinstance(dataset, h5py.Dataset):
            raise ValueError(f"Dataset missing in .ims file: {spec.dataset_key}")
        if spec.axis_order == "zyx":
            arr = dataset[tile.z0 : tile.z1, tile.y0 : tile.y1, tile.x0 : tile.x1]
        elif spec.axis_order == "xyz":
            arr = np.transpose(
                dataset[tile.x0 : tile.x1, tile.y0 : tile.y1, tile.z0 : tile.z1],
                (2, 1, 0),
            )
        else:
            raise ValueError(f"Unsupported axis_order in spec: {spec.axis_order}")
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spec.spacing_xyz)
    img.SetOrigin(
        (
            tile.x0 * spec.spacing_xyz[0],
            tile.y0 * spec.spacing_xyz[1],
            tile.z0 * spec.spacing_xyz[2],
        )
    )
    return img


def read_ims_downsampled_from_spec(spec: ImsVolumeSpec, stride_zyx: Tuple[int, int, int]) -> sitk.Image:
    if h5py is None:
        raise RuntimeError("Reading .ims requires h5py. Install with: pip install h5py")
    sz, sy, sx = (max(1, int(stride)) for stride in stride_zyx)
    spacing_xyz = (spec.spacing_xyz[0] * sx, spec.spacing_xyz[1] * sy, spec.spacing_xyz[2] * sz)
    with h5py.File(str(spec.path), "r") as h5f:
        pyramid_match = _matching_pyramid_dataset(h5f, spec, (sz, sy, sx))
        if pyramid_match is not None:
            dataset, shape_zyx, effective_scale_zyx, level = pyramid_match
            source_label = (
                f"ResolutionLevel {level} shape_zyx={shape_zyx} "
                f"effective_scale_zyx=({effective_scale_zyx[0]:.3g}, "
                f"{effective_scale_zyx[1]:.3g}, {effective_scale_zyx[2]:.3g})"
            )
            spacing_xyz = _pyramid_spacing_xyz(spec, effective_scale_zyx)
            arr = _ims_dataset_to_zyx_array(dataset, spec.axis_order)
            print(
                f"[ims-pyramid] using {source_label} for requested stride_zyx={(sz, sy, sx)}: {spec.path}",
                flush=True,
            )
        else:
            dataset = h5f[spec.dataset_key]
            if not isinstance(dataset, h5py.Dataset):
                raise ValueError(f"Dataset missing in .ims file: {spec.dataset_key}")
            if spec.axis_order == "zyx":
                arr = dataset[::sz, ::sy, ::sx]
            elif spec.axis_order == "xyz":
                arr = np.transpose(dataset[::sx, ::sy, ::sz], (2, 1, 0))
            else:
                raise ValueError(f"Unsupported axis_order in spec: {spec.axis_order}")
            if (sz, sy, sx) != (1, 1, 1):
                print(
                    f"[ims-pyramid] no matching .ims pyramid level for requested stride_zyx={(sz, sy, sx)}; "
                    f"using direct strided read: {spec.path}",
                    flush=True,
                )
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spacing_xyz)
    return img
