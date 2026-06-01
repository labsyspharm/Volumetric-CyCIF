#!/usr/bin/env python3
"""Download stitched TIFF slices from SFTP and convert them to a multi-channel .ims."""

from __future__ import annotations

import argparse
import os
import posixpath
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import tifffile

try:
    import paramiko
except ImportError:  # pragma: no cover - surfaced as a friendly runtime error
    paramiko = None


DEFAULT_REMOTE_BASE = "/data/alw749/reg"
DEFAULT_CHANNEL_DIR_GLOB = "Ex_*_deskewed_stitched"


@dataclass(frozen=True)
class ChannelVolume:
    name: str
    directory: Path
    slices: tuple[Path, ...]


def parse_triplet(value: str) -> tuple[float, float, float]:
    parts = [part.strip() for part in value.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("Expected x,y,z")
    try:
        triplet = tuple(float(part) for part in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Voxel sizes must be numeric") from exc
    if any(component <= 0 for component in triplet):
        raise argparse.ArgumentTypeError("Voxel sizes must be > 0")
    return triplet  # type: ignore[return-value]


def natural_sort_key(path: Path) -> list[object]:
    return [int(token) if token.isdigit() else token.lower() for token in re.split(r"(\d+)", path.name)]


def resolve_dataset_name(args: argparse.Namespace) -> str:
    if args.dataset_name:
        return args.dataset_name
    if args.local_dataset_dir:
        return Path(args.local_dataset_dir).resolve().name
    raise ValueError("--dataset-name is required unless --local-dataset-dir is provided")


def resolve_output_stem(args: argparse.Namespace, dataset_name: str) -> str:
    if args.output_name:
        return Path(args.output_name).stem
    return dataset_name


def ensure_writer_binary(path_arg: str | None) -> Path:
    candidates: list[Path] = []
    if path_arg:
        candidates.append(Path(path_arg).expanduser())

    which_hit = shutil.which("parallelimariswriter")
    if which_hit:
        candidates.append(Path(which_hit))

    for candidate in candidates:
        resolved = candidate.expanduser()
        if resolved.exists() and os.access(resolved, os.X_OK):
            return resolved.resolve()

    search_paths = "\n".join(f"  - {candidate}" for candidate in candidates)
    raise FileNotFoundError(
        "Could not find an executable `parallelimariswriter` binary.\n"
        "Build or install the writer first, then pass it with --writer-binary.\n"
        f"Searched:\n{search_paths}"
    )


def sftp_connect(args: argparse.Namespace) -> "paramiko.SFTPClient":
    if paramiko is None:
        raise RuntimeError("paramiko is required for SFTP downloads. Install with `pip install -r requirements.txt`.")

    password = os.environ.get(args.sftp_password_env) if args.sftp_password_env else None
    passphrase = os.environ.get(args.sftp_passphrase_env) if args.sftp_passphrase_env else None

    client = paramiko.SSHClient()
    client.load_system_host_keys()
    if args.strict_host_key_checking:
        client.set_missing_host_key_policy(paramiko.RejectPolicy())
    else:
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    connect_kwargs = {
        "hostname": args.sftp_host,
        "port": args.sftp_port,
        "username": args.sftp_user,
        "password": password,
        "passphrase": passphrase,
        "look_for_keys": True,
        "allow_agent": True,
    }
    if args.sftp_key:
        connect_kwargs["key_filename"] = str(Path(args.sftp_key).expanduser())

    client.connect(**connect_kwargs)
    sftp = client.open_sftp()
    setattr(sftp, "_ssh_client", client)
    return sftp


def close_sftp(sftp: "paramiko.SFTPClient | None") -> None:
    if sftp is None:
        return
    ssh_client = getattr(sftp, "_ssh_client", None)
    sftp.close()
    if ssh_client is not None:
        ssh_client.close()


def iter_remote_files(sftp: "paramiko.SFTPClient", remote_root: str) -> Iterable[str]:
    entries = sorted(sftp.listdir_attr(remote_root), key=lambda entry: entry.filename)
    for entry in entries:
        remote_path = posixpath.join(remote_root, entry.filename)
        if stat.S_ISDIR(entry.st_mode):
            yield from iter_remote_files(sftp, remote_path)
        elif stat.S_ISREG(entry.st_mode):
            yield remote_path


def download_remote_tree(sftp: "paramiko.SFTPClient", remote_root: str, local_root: Path) -> Path:
    try:
        attrs = sftp.stat(remote_root)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Remote dataset not found: {remote_root}") from exc
    if not stat.S_ISDIR(attrs.st_mode):
        raise NotADirectoryError(f"Remote dataset path is not a directory: {remote_root}")

    dataset_dir = local_root / Path(remote_root.rstrip("/")).name
    dataset_dir.mkdir(parents=True, exist_ok=True)

    remote_files = list(iter_remote_files(sftp, remote_root))
    if not remote_files:
        raise FileNotFoundError(f"No files found under remote dataset: {remote_root}")

    total = len(remote_files)
    for index, remote_file in enumerate(remote_files, start=1):
        relative = posixpath.relpath(remote_file, remote_root)
        local_file = dataset_dir / Path(relative)
        local_file.parent.mkdir(parents=True, exist_ok=True)
        print(f"[download {index}/{total}] {remote_file} -> {local_file}", flush=True)
        sftp.get(remote_file, str(local_file))
    return dataset_dir


def discover_channels(dataset_dir: Path, channel_dir_glob: str) -> list[ChannelVolume]:
    candidate_dirs = sorted(
        [path for path in dataset_dir.glob(channel_dir_glob) if path.is_dir()],
        key=lambda path: natural_sort_key(path),
    )
    if not candidate_dirs:
        candidate_dirs = sorted(
            [path for path in dataset_dir.iterdir() if path.is_dir() and any(path.glob("*.tif*"))],
            key=lambda path: natural_sort_key(path),
        )

    channels: list[ChannelVolume] = []
    for directory in candidate_dirs:
        slices = sorted(
            [*directory.glob("*.tif"), *directory.glob("*.tiff"), *directory.glob("*.TIF"), *directory.glob("*.TIFF")],
            key=natural_sort_key,
        )
        if not slices:
            continue
        channels.append(ChannelVolume(name=directory.name, directory=directory, slices=tuple(slices)))

    if not channels:
        raise FileNotFoundError(f"No channel folders with TIFF slices found under {dataset_dir}")

    slice_counts = {len(channel.slices) for channel in channels}
    if len(slice_counts) != 1:
        details = ", ".join(f"{channel.name}={len(channel.slices)}" for channel in channels)
        raise ValueError(f"Channel folders do not contain the same number of TIFF slices: {details}")

    return channels


def validate_channel_slices(channels: Sequence[ChannelVolume]) -> tuple[tuple[int, int], str]:
    reference_shape: tuple[int, int] | None = None
    reference_dtype: str | None = None

    for channel in channels:
        image = tifffile.imread(channel.slices[0])
        if image.ndim != 2:
            raise ValueError(f"Expected 2D TIFF slices, got shape {image.shape} in {channel.slices[0]}")
        shape = (int(image.shape[0]), int(image.shape[1]))
        dtype = str(image.dtype)
        if reference_shape is None:
            reference_shape = shape
            reference_dtype = dtype
        elif shape != reference_shape or dtype != reference_dtype:
            raise ValueError(
                "All channels must share a single slice shape and dtype. "
                f"Mismatch at {channel.slices[0]}: shape={shape}, dtype={dtype}, "
                f"expected shape={reference_shape}, dtype={reference_dtype}"
            )

    if reference_shape is None or reference_dtype is None:
        raise RuntimeError("No TIFF slices were validated")

    return reference_shape, reference_dtype


def assemble_channel_stack(channel: ChannelVolume, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    stage_name = re.sub(r"[^A-Za-z0-9._-]+", "_", channel.name)
    output_path = output_dir / f"{stage_name}.tif"
    print(f"[assemble] {channel.name} -> {output_path}", flush=True)

    expected_shape: tuple[int, int] | None = None
    expected_dtype: str | None = None
    with tifffile.TiffWriter(output_path, bigtiff=True) as writer:
        for slice_index, slice_path in enumerate(channel.slices):
            image = tifffile.imread(slice_path)
            if image.ndim != 2:
                raise ValueError(f"Expected 2D TIFF slices, got shape {image.shape} in {slice_path}")
            shape = (int(image.shape[0]), int(image.shape[1]))
            dtype = str(image.dtype)
            if expected_shape is None:
                expected_shape = shape
                expected_dtype = dtype
            elif shape != expected_shape or dtype != expected_dtype:
                raise ValueError(
                    f"Channel {channel.name} contains mismatched slice data at {slice_path}: "
                    f"shape={shape}, dtype={dtype}, expected shape={expected_shape}, dtype={expected_dtype}"
                )
            writer.write(
                image,
                contiguous=True,
                compression=None,
                metadata=None,
                photometric="minisblack",
            )
            if slice_index == 0 or (slice_index + 1) == len(channel.slices) or (slice_index + 1) % 250 == 0:
                print(
                    f"  [assemble {channel.name}] wrote slice {slice_index + 1}/{len(channel.slices)}",
                    flush=True,
                )
    return output_path


def run_parallel_imaris_writer(
    writer_binary: Path,
    staged_tiffs: Sequence[Path],
    output_dir: Path,
    output_stem: str,
    voxel_size_xyz: tuple[float, float, float],
    block_size_xyz: str | None,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    patterns = [stage_path.stem for stage_path in staged_tiffs]
    cmd = [
        str(writer_binary),
        "-F",
        str(staged_tiffs[0].parent),
        "-P",
        ",".join(patterns),
        "-c",
        str(len(staged_tiffs)),
        "-t",
        "1",
        "-o",
        str(output_dir),
        "-n",
        output_stem,
        "-r",
        "tiff",
        "-v",
        ",".join(f"{value:g}" for value in voxel_size_xyz),
    ]
    if block_size_xyz:
        cmd.extend(["-b", block_size_xyz])

    print("[writer] " + " ".join(cmd), flush=True)
    run = subprocess.run(cmd, check=False)
    if run.returncode != 0:
        raise RuntimeError(f"Parallel Imaris Writer failed with exit code {run.returncode}")

    output_path = output_dir / f"{output_stem}.ims"
    if not output_path.exists():
        raise FileNotFoundError(f"Writer completed without producing {output_path}")
    return output_path


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download stitched TIFF slices from SFTP and convert them into a multi-channel .ims file."
    )
    parser.add_argument("--dataset-name", help="Dataset name under the remote base path")
    parser.add_argument("--remote-base", default=DEFAULT_REMOTE_BASE, help=f"Remote dataset base (default: {DEFAULT_REMOTE_BASE})")
    parser.add_argument("--local-dataset-dir", help="Use an already-downloaded local dataset directory instead of SFTP")
    parser.add_argument("--output-dir", required=True, help="Directory that will receive the final .ims output")
    parser.add_argument("--output-name", help="Output file stem (default: dataset name)")
    parser.add_argument("--temp-root", help="Temporary workspace root (default: system temp dir)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep the downloaded and staged TIFF data")
    parser.add_argument(
        "--channel-dir-glob",
        default=DEFAULT_CHANNEL_DIR_GLOB,
        help=f"Glob used to find channel directories (default: {DEFAULT_CHANNEL_DIR_GLOB})",
    )
    parser.add_argument("--writer-binary", help="Path to the compiled parallelimariswriter executable")
    parser.add_argument("--voxel-size-um", default="1.0,1.0,1.0", help="Voxel size as x,y,z in microns")
    parser.add_argument("--block-size-xyz", help="Optional writer block size override as x,y,z")

    parser.add_argument("--sftp-host", default=os.environ.get("SFTP_HOST"), help="SFTP host")
    parser.add_argument("--sftp-port", type=int, default=int(os.environ.get("SFTP_PORT", "22")), help="SFTP port")
    parser.add_argument("--sftp-user", default=os.environ.get("SFTP_USER"), help="SFTP username")
    parser.add_argument("--sftp-key", default=os.environ.get("SFTP_KEY"), help="Path to a private key file")
    parser.add_argument(
        "--sftp-password-env",
        default="SFTP_PASSWORD",
        help="Environment variable that stores the SFTP password (default: SFTP_PASSWORD)",
    )
    parser.add_argument(
        "--sftp-passphrase-env",
        default="SFTP_PASSPHRASE",
        help="Environment variable that stores the private key passphrase (default: SFTP_PASSPHRASE)",
    )
    parser.add_argument(
        "--strict-host-key-checking",
        action="store_true",
        help="Reject unknown SFTP host keys instead of auto-adding them",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    dataset_name = resolve_dataset_name(args)
    output_stem = resolve_output_stem(args, dataset_name)
    voxel_size_xyz = parse_triplet(args.voxel_size_um)
    writer_binary = ensure_writer_binary(args.writer_binary)

    temp_parent = Path(args.temp_root).expanduser() if args.temp_root else None
    temp_dir = Path(tempfile.mkdtemp(prefix="ims_convert_", dir=temp_parent))
    sftp = None

    try:
        if args.local_dataset_dir:
            dataset_dir = Path(args.local_dataset_dir).expanduser().resolve()
            if not dataset_dir.is_dir():
                raise NotADirectoryError(f"--local-dataset-dir is not a directory: {dataset_dir}")
        else:
            missing = [name for name in ("sftp_host", "sftp_user") if not getattr(args, name)]
            if missing:
                parser.error(
                    "SFTP mode requires "
                    + ", ".join(f"--{name.replace('_', '-')}" for name in missing)
                    + " (or matching environment variables)"
                )
            remote_dataset = posixpath.join(args.remote_base.rstrip("/"), dataset_name)
            print(f"[sftp] downloading {remote_dataset}", flush=True)
            sftp = sftp_connect(args)
            dataset_dir = download_remote_tree(sftp, remote_dataset, temp_dir / "download")

        channels = discover_channels(dataset_dir, args.channel_dir_glob)
        slice_shape, dtype = validate_channel_slices(channels)
        print(
            f"[dataset] {dataset_dir} | channels={len(channels)} | z_slices={len(channels[0].slices)} "
            f"| yx={slice_shape} | dtype={dtype}",
            flush=True,
        )

        stage_dir = temp_dir / "staged_volumes"
        staged_tiffs = [assemble_channel_stack(channel, stage_dir) for channel in channels]
        ims_path = run_parallel_imaris_writer(
            writer_binary=writer_binary,
            staged_tiffs=staged_tiffs,
            output_dir=Path(args.output_dir).expanduser().resolve(),
            output_stem=output_stem,
            voxel_size_xyz=voxel_size_xyz,
            block_size_xyz=args.block_size_xyz,
        )
        print(f"[done] wrote {ims_path}", flush=True)
        return 0
    finally:
        close_sftp(sftp)
        if args.keep_temp:
            print(f"[temp] kept {temp_dir}", flush=True)
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
