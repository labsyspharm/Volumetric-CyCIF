#!/usr/bin/env python3
# Developed by Alex Wong
# Cite: Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in
# Preprint: https://doi.org/10.64898/2026.05.17.725158
# Registration method: ANTsX/ANTsPy - https://github.com/ANTsX/ANTsPy

"""Export one full-resolution .ims channel to an uncompressed BigTIFF."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import h5py
import numpy as np
import tifffile
from tqdm import tqdm


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ims", required=True, help="Source Imaris .ims file")
    parser.add_argument("--ch", type=int, required=True, help="Channel index to export")
    parser.add_argument("--out", required=True, help="Output BigTIFF path")
    parser.add_argument("--level", type=int, default=0, help="Imaris resolution level. Default: 0")
    parser.add_argument("--timepoint", type=int, default=0, help="Imaris timepoint. Default: 0")
    parser.add_argument(
        "--block-z",
        "--bz",
        type=int,
        default=128,
        help="Z planes read per block. Larger is usually faster if RAM permits. Default: 128",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing output TIFF")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    source = Path(args.ims).expanduser().resolve()
    output = Path(args.out).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists; pass --overwrite to replace it: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    dataset_key = f"DataSet/ResolutionLevel {args.level}/TimePoint {args.timepoint}/Channel {args.ch}/Data"
    block_z = max(1, int(args.block_z))
    started = time.time()
    with h5py.File(source, "r") as handle:
        if dataset_key not in handle:
            raise KeyError(f"Dataset not found: /{dataset_key}")
        dataset = handle[dataset_key]
        if dataset.ndim != 3:
            raise ValueError(f"Expected a 3D ZYX dataset, found shape={dataset.shape}")
        z_size, y_size, x_size = (int(v) for v in dataset.shape)
        block_z = min(block_z, z_size)
        buffer = np.empty((block_z, y_size, x_size), dtype=dataset.dtype)
        gib = dataset.size * dataset.dtype.itemsize / 1024**3
        print(
            f"[export] source=/{dataset_key} shape_zyx={dataset.shape} dtype={dataset.dtype} "
            f"size={gib:.2f} GiB block_z={block_z}",
            flush=True,
        )
        print(f"[export] output={output}", flush=True)

        with tifffile.TiffWriter(output, bigtiff=True) as writer:
            with tqdm(total=z_size, unit="plane", desc=f"export ch{args.ch}", dynamic_ncols=True) as progress:
                for z0 in range(0, z_size, block_z):
                    z1 = min(z_size, z0 + block_z)
                    view = buffer[: z1 - z0]
                    dataset.read_direct(view, source_sel=np.s_[z0:z1, :, :])
                    writer.write(
                        view,
                        photometric="minisblack",
                        compression=None,
                        contiguous=True,
                        metadata=None,
                    )
                    progress.update(z1 - z0)

    elapsed = time.time() - started
    print(f"[export] finished in {elapsed / 60:.1f} min: {output}", flush=True)
    with tifffile.TiffFile(output) as tif:
        print(f"[export] verification shape={tif.series[0].shape} axes={tif.series[0].axes}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
