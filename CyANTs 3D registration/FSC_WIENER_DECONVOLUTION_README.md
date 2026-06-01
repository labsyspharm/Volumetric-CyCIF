# FSC/SFSC-Guided Wiener Deconvolution

`fsc_wiener_deconvolution.py` runs chunked 3D Wiener deconvolution for lightsheet volumes. It can read Imaris/HDF5 volumes or a folder of serial TIFF slices, estimate an effective Gaussian PSF from sectioned FSC, use one GPU for small volumes, automatically tile large volumes, and optionally write QC figures.

## Quick Start

Minimal stitched Imaris input:

```bash
python /Users/alexwong/Documents/Codex/registration/fsc_wiener_deconvolution.py \
  --input /path/to/sample.ims \
  --output /path/to/ch0_deconvolved.tif \
  --ims-channel 0 \
  --gpu-memory-gb 11 \
  --overwrite
```

Minimal stitched serial TIFF folder input:

```bash
python /Users/alexwong/Documents/Codex/registration/fsc_wiener_deconvolution.py \
  --input /path/to/stitched_tiff_folder \
  --output /path/to/ch0_deconvolved.tif \
  --tiff-channel-token _ch0 \
  --spacing-xyz-um 0.2,0.2,1.0 \
  --gpu-memory-gb 11 \
  --overwrite
```

Recommended 8-GPU run with QC:

```bash
python /Users/alexwong/Documents/Codex/registration/fsc_wiener_deconvolution.py \
  --input /path/to/stitched_tiff_folder \
  --input-stage post-stitch \
  --tiff-channel-token _ch0 \
  --spacing-xyz-um 0.2,0.2,1.0 \
  --output /path/to/ch0_deconvolved.tif \
  --qc-dir /path/to/ch0_deconvolved_qc \
  --gpu-memory-gb 11 \
  --tiling auto \
  --gpus 0,1,2,3,4,5,6,7 \
  --workers 8 \
  --psf-mode global \
  --wiener-k 0.005 \
  --output-dtype uint16 \
  --clip-negative \
  --verbose \
  --overwrite
```

For pre-stitch deskewed blocks, run one block at a time:

```bash
python /Users/alexwong/Documents/Codex/registration/fsc_wiener_deconvolution.py \
  --input /path/to/deskewed_block_ch0 \
  --input-stage pre-stitch \
  --tiff-channel-token _ch0 \
  --spacing-xyz-um 0.2,0.2,1.0 \
  --output /path/to/deconvolved_blocks/block001_ch0_deconvolved.tif \
  --gpu-memory-gb 11 \
  --tiling auto \
  --overwrite
```

## What Defaults Do

With only required input/output arguments, the script uses:

- `--tiling auto`: one full-volume GPU job if it fits, tiled processing if not.
- `--device auto`: CUDA if available.
- `--gpus 0 --workers 1`: one GPU by default.
- `--psf-mode global`: estimate one PSF from a central crop.
- `--wiener-k 0.005`: moderate Wiener regularization.
- `--output-dtype float32`: preserve deconvolved dynamic range.
- `--input-stage post-stitch`: assume a stitched volume unless told otherwise.

For your 8-GPU PC, usually add only:

```bash
--gpu-memory-gb 11 --gpus 0,1,2,3,4,5,6,7 --workers 8
```

## Input Modes

Imaris/HDF5:

```bash
--input sample.ims
--ims-resolution-level 0
--ims-timepoint 0
--ims-channel 0
```

Explicit HDF5 dataset path if auto-discovery fails:

```bash
--ims-dataset-path "/DataSet/ResolutionLevel 0/TimePoint 0/Channel 0/Data"
```

Serial TIFF folder:

```bash
--input /path/to/slices
--tiff-channel-token _ch0
--spacing-xyz-um 0.2,0.2,1.0
```

The TIFF channel token filters filenames. For example, `_ch0` selects files whose names contain `_ch0`.

## Output Modes

HDF5 output:

```bash
--output /path/to/deconvolved.h5
```

The output dataset is named `deconvolved`.

BigTIFF output:

```bash
--output /path/to/deconvolved.tif
```

For TIFF output, the script stages results in temporary HDF5, then streams slices into a BigTIFF stack. Keep the staging file with:

```bash
--keep-staging-h5
```

Or choose its location:

```bash
--staging-h5 /fast_disk/tmp_deconvolution.h5
```

## QC Outputs

Enable QC:

```bash
--qc-dir /path/to/qc
```

Optional custom QC planes:

```bash
--qc-slices 120,800,900
```

This means `z,y,x`. `auto` uses center planes.

QC files:

- `before_after_xy_*.png`
- `before_after_xz_*.png`
- `before_after_yz_*.png`
- `tile_psf_timing.png`
- `qc_summary.json`

If `matplotlib` is not installed, deconvolution still runs and QC figures are skipped.

## Important Choices

### `--input-stage`

Use `post-stitch` when the input is a stitched volume. Use `pre-stitch` when the input is one deskewed acquisition block before stitching. The deconvolution math is the same, but this label is stored in metadata and changes the QC interpretation.

Recommended:

- If you still have deskewed pre-stitch blocks, deconvolve each block before stitching.
- If you only have the stitched `.ims` or serial TIFFs, use `post-stitch`.

### `--psf-mode`

- `global`: estimate one PSF from a central crop. Recommended first run.
- `local`: estimate a PSF per processing tile and fall back to global if unstable.
- `fixed`: use `--fixed-fwhm-um x,y,z`.

Start with `global`. Try `local` only after checking QC and if blur clearly varies across depth or position.

### `--tiling`

- `auto`: process full volume on one GPU if estimated memory fits; otherwise tile.
- `on`: always tile.
- `off`: force one full-volume GPU job.

For 11 GB GPUs:

```bash
--gpu-memory-gb 11 --tiling auto
```

The script prints the memory estimate and decision.

## Full Flag Reference

### Required

| Flag | Default | Meaning |
| --- | --- | --- |
| `--input` | required | `.ims`, `.h5`, `.hdf5`, or folder of serial TIFF slices. |
| `--output` | required | `.h5`, `.hdf5`, `.tif`, or `.tiff`. |

### Input identity

| Flag | Default | Meaning |
| --- | --- | --- |
| `--input-stage` | `post-stitch` | `post-stitch`, `pre-stitch`, or `unknown`. |
| `--ims-resolution-level` | `0` | Imaris resolution level. |
| `--ims-timepoint` | `0` | Imaris timepoint. |
| `--ims-channel` | `0` | Imaris channel index. |
| `--ims-dataset-path` | empty | Explicit HDF5 dataset path. |
| `--ims-axis-order` | `auto` | `auto`, `zyx`, or `xyz`. |
| `--tiff-channel-token` | empty | Filename token for TIFF channel selection, e.g. `_ch0`. |
| `--validate-tiff-slices` | off | Check every TIFF slice shape and dtype before processing. |
| `--spacing-xyz-um` | empty | Voxel size as `x,y,z` in microns. Important for TIFF input. |

### Output

| Flag | Default | Meaning |
| --- | --- | --- |
| `--overwrite` | off | Replace existing output. |
| `--staging-h5` | auto next to output | Temporary HDF5 path for TIFF output. |
| `--keep-staging-h5` | off | Keep the staging HDF5 after TIFF writing. |
| `--output-dtype` | `float32` | Output dtype: `float32`, `uint16`, `same`, etc. |
| `--compression` | `lzf` | HDF5 compression: `lzf`, `gzip`, or `none`. |

### GPU and tiling

| Flag | Default | Meaning |
| --- | --- | --- |
| `--device` | `auto` | `auto`, `cuda`, or `cpu`. |
| `--gpus` | `0` | GPU IDs, e.g. `0,1,2,3,4,5,6,7`, or `auto`. |
| `--workers` | `1` | Parallel tile workers. Usually one per selected GPU. |
| `--tiling` | `auto` | `auto`, `on`, or `off`. |
| `--gpu-memory-gb` | `auto` | GPU memory per worker. Use `11` for 11 GB cards. |
| `--memory-safety-factor` | `0.65` | Fraction of GPU memory allowed for auto full-volume mode. |
| `--memory-bytes-per-voxel` | `96` | Conservative FFT memory estimate. Increase if auto mode OOMs. |
| `--tile-size` | `96,256,256` | Core tile size as `z,y,x`. |
| `--overlap` | `24,64,64` | Halo/overlap as `z,y,x`. |

### PSF and deconvolution

| Flag | Default | Meaning |
| --- | --- | --- |
| `--psf-mode` | `global` | `global`, `local`, or `fixed`. |
| `--fixed-fwhm-um` | empty | Fixed Gaussian PSF FWHM as `x,y,z` in microns. |
| `--section-angle-deg` | `15` | Sectioned FSC wedge half-angle. |
| `--fsc-bins` | `64` | Number of frequency bins for FSC. |
| `--min-fsc-points` | `32` | Minimum Fourier voxels per FSC bin. |
| `--threshold-mode` | `one-bit` | `one-bit` or `fixed`. |
| `--threshold-snr-e` | `0.5` | Expected SNR for one-bit FSC threshold. |
| `--fixed-threshold` | `1/7` | Fixed FSC threshold when using `--threshold-mode fixed`. |
| `--wiener-k` | `0.005` | Wiener regularization. Lower is sharper/noisier. |
| `--background-percentile` | `0` | Subtract percentile background before deconvolution. |
| `--clip-negative` | off | Clip restored negative values to zero. |

### Foreground and robustness

| Flag | Default | Meaning |
| --- | --- | --- |
| `--min-foreground-fraction` | `0.005` | Minimum foreground fraction for FSC estimation. |
| `--foreground-threshold-fraction` | `0.02` | Foreground threshold as fraction of local intensity range. |

### Progress and QC

| Flag | Default | Meaning |
| --- | --- | --- |
| `--verbose` | off | Print extra QC details and PSF fallback reasons. |
| `--quiet` | off | Suppress per-tile progress lines. |
| `--progress-every` | `1` | Print progress every N completed tiles. |
| `--qc-dir` | empty | Directory for QC PNGs and JSON. |
| `--qc-slices` | `auto` | QC planes as `z,y,x`, or `auto`. |

## Suggested First Pass

Use global PSF and QC first:

```bash
--psf-mode global --qc-dir /path/to/qc --verbose
```

Inspect the before/after panels and `tile_psf_timing.png`. If blur varies strongly by depth or region, test:

```bash
--psf-mode local
```

on a smaller crop or representative subset before running the entire dataset.

