#!/usr/bin/env python3
"""Optional desktop launcher for CyANTs command-line workflows.

The GUI intentionally stays thin: it collects common paths/options, shows the
exact command, then runs the existing headless scripts in a subprocess.
"""

from __future__ import annotations

import json
import os
import queue
import shlex
import subprocess
import sys
import threading
import contextlib
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox
import tkinter as tk
from tkinter import ttk

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD  # type: ignore
except Exception:  # pragma: no cover - optional GUI nicety
    DND_FILES = None
    TkinterDnD = None


APP_TITLE = "CyANTs Launcher"
THEME = {
    "bg": "#07111f",
    "panel": "#0c1f33",
    "panel_alt": "#102a43",
    "border": "#1d4f75",
    "text": "#e6f7ff",
    "muted": "#9bc7d9",
    "accent": "#26d6d0",
    "accent_blue": "#1aa7ff",
    "field": "#f7fbff",
    "field_text": "#07111f",
    "log_bg": "#04101d",
    "log_fg": "#c8f3ff",
}
DEFAULT_SPACING = "0.711"
DEFAULT_THREADS = "32"
DEFAULT_PROGRESS_INTERVAL = "30"
COMMAND_FORMAT_AUTO = "Auto platform"
COMMAND_FORMAT_WINDOWS = "Windows CMD"
COMMAND_FORMAT_POSIX = "macOS/Linux/HPC"
COMMAND_FORMAT_LOCAL = "Exact local run"
RUNNER_FLAG = "--cyants-run-script"
RUNNABLE_SCRIPTS = {
    "vcycif_ROI.py",
    "ants_ims_tiled_quicksyn.py",
    "ants_ims_intracycle_whole.py",
}


@dataclass
class Field:
    label: str
    var: tk.StringVar
    browse: str | None = None
    width: int = 72


def repo_script(name: str) -> str:
    if getattr(sys, "frozen", False):
        return str(Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent)).joinpath(name))
    return str(Path(__file__).resolve().with_name(name))


def resource_path(relative_path: str) -> Path:
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent)).joinpath(relative_path)
    return Path(__file__).resolve().parent.joinpath(relative_path)


def packaged_command(script_name: str, args: list[str]) -> list[str]:
    return [sys.executable, RUNNER_FLAG, script_name, *args]


def local_script_command(script_name: str, args: list[str]) -> list[str]:
    if getattr(sys, "frozen", False):
        return packaged_command(script_name, args)
    return [sys.executable, repo_script(script_name), *args]


def quote_command(command: list[str], command_format: str = COMMAND_FORMAT_AUTO) -> str:
    if command_format == COMMAND_FORMAT_LOCAL:
        command_format = COMMAND_FORMAT_AUTO
    if command_format == COMMAND_FORMAT_WINDOWS or (command_format == COMMAND_FORMAT_AUTO and os.name == "nt"):
        return subprocess.list2cmdline(command)
    return shlex.join(command)


def portable_command(command: list[str]) -> list[str]:
    """Return a repo-root command suitable for copying to another machine."""
    if len(command) < 2:
        return command
    if len(command) >= 4 and command[1] == RUNNER_FLAG:
        return ["python", command[2], *command[3:]]
    script_name = Path(command[1]).name
    return ["python", script_name, *command[2:]]


def parse_bool(value: tk.BooleanVar, flag: str) -> list[str]:
    return [flag] if value.get() else []


def parse_roi_csv_points(path: Path) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        first_line = next((line for line in sample.splitlines() if line.strip()), "")
        has_header = all(token in first_line.lower() for token in ("x", "y"))
        if has_header:
            reader = csv.DictReader(handle)
            for row in reader:
                lowered = {str(key).strip().lower(): value for key, value in row.items() if key is not None}
                x_value = lowered.get("x") or lowered.get("x coordinate") or lowered.get("xcoordinate")
                y_value = lowered.get("y") or lowered.get("y coordinate") or lowered.get("ycoordinate")
                if x_value is None or y_value is None:
                    continue
                points.append((float(x_value), float(y_value)))
        else:
            reader = csv.reader(handle)
            for row in reader:
                values = [value.strip() for value in row if value.strip()]
                if len(values) >= 2:
                    try:
                        points.append((float(values[0]), float(values[1])))
                    except ValueError:
                        continue
    if len(points) < 3:
        raise ValueError(f"ROI CSV needs at least three X,Y points: {path}")
    return points


def preview_roi_to_full_points(
    roi: tuple[float, float, float, float] | list[float],
    preview_shape_yx: tuple[int, int],
    full_shape_zyx: tuple[int, int, int],
) -> list[tuple[int, int]]:
    """Convert a preview-space ROI box into full-resolution original .ims X/Y points."""
    preview_y, preview_x = preview_shape_yx
    _full_z, full_y, full_x = full_shape_zyx
    if preview_x <= 0 or preview_y <= 0 or full_x <= 0 or full_y <= 0:
        raise ValueError("Preview and full-resolution shapes must be positive")
    x0, y0, x1, y1 = roi
    x0, x1 = sorted((float(x0), float(x1)))
    y0, y1 = sorted((float(y0), float(y1)))
    x0 = max(0.0, min(float(preview_x - 1), x0))
    x1 = max(0.0, min(float(preview_x - 1), x1))
    y0 = max(0.0, min(float(preview_y - 1), y0))
    y1 = max(0.0, min(float(preview_y - 1), y1))

    def scale_x(value: float) -> int:
        return max(0, min(full_x - 1, int(round(value * full_x / preview_x))))

    def scale_y(value: float) -> int:
        return max(0, min(full_y - 1, int(round(value * full_y / preview_y))))

    return [
        (scale_x(x0), scale_y(y0)),
        (scale_x(x1), scale_y(y0)),
        (scale_x(x1), scale_y(y1)),
        (scale_x(x0), scale_y(y1)),
    ]


def choose_xy_grid(tile_count: int, shape_zyx: tuple[int, int, int]) -> tuple[int, int]:
    if tile_count < 1:
        raise ValueError("Tiles must be >= 1")
    _z_size, y_size, x_size = shape_zyx
    aspect = max(float(x_size) / max(float(y_size), 1.0), 1e-9)
    candidates = []
    for rows in range(1, int(math.sqrt(tile_count)) + 1):
        if tile_count % rows:
            continue
        cols = tile_count // rows
        grid_aspect = cols / rows
        score = abs(math.log(max(grid_aspect, 1e-9) / aspect))
        candidates.append((score, rows, cols))
    _score, rows, cols = min(candidates)
    return rows, cols


def auto_axis_tile(length: int, count: int, overlap_fraction: float) -> tuple[int, int]:
    if count <= 1:
        return length, 0
    denominator = count - (count - 1) * overlap_fraction
    tile = int(math.ceil(length / denominator))
    overlap = int(round(tile * overlap_fraction))
    overlap = max(0, min(overlap, max(tile - 1, 0)))
    return min(tile, length), overlap


def compute_starts(length: int, tile: int, overlap: int) -> list[int]:
    if tile >= length:
        return [0]
    stride = tile - overlap
    if stride <= 0:
        raise ValueError("Tile size must be greater than overlap.")
    starts = list(range(0, max(length - tile + 1, 1), stride))
    last = length - tile
    if starts[-1] != last:
        starts.append(last)
    return starts


def compute_auto_tiles(
    shape_zyx: tuple[int, int, int],
    tile_count: int,
    overlap_fraction: float,
) -> tuple[tuple[int, int, int], tuple[int, int, int], list[tuple[int, int, int, int, int, int]], tuple[int, int]]:
    rows, cols = choose_xy_grid(tile_count, shape_zyx)
    z_size, y_size, x_size = shape_zyx
    tile_y, overlap_y = auto_axis_tile(y_size, rows, overlap_fraction)
    tile_x, overlap_x = auto_axis_tile(x_size, cols, overlap_fraction)
    y_starts = compute_starts(y_size, tile_y, overlap_y)
    x_starts = compute_starts(x_size, tile_x, overlap_x)
    tiles = [
        (0, z_size, y0, min(y0 + tile_y, y_size), x0, min(x0 + tile_x, x_size))
        for y0 in y_starts
        for x0 in x_starts
    ][:tile_count]
    return (z_size, tile_y, tile_x), (0, overlap_y, overlap_x), tiles, (rows, cols)


class CommandRunner:
    def __init__(self, log: tk.Text, on_done):
        self.log = log
        self.on_done = on_done
        self.process: subprocess.Popen[str] | None = None
        self.output_queue: queue.Queue[str | None] = queue.Queue()

    def running(self) -> bool:
        return self.process is not None and self.process.poll() is None

    def append(self, text: str) -> None:
        self.log.insert("end", text)
        self.log.see("end")

    def start(self, command: list[str]) -> None:
        if self.running():
            messagebox.showinfo(APP_TITLE, "A command is already running.")
            return
        self.append("\n[gui] running:\n" + quote_command(command) + "\n\n")
        self.process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        threading.Thread(target=self._reader, daemon=True).start()
        self.log.after(100, self._drain)

    def stop(self) -> None:
        if self.running() and self.process is not None:
            self.append("\n[gui] terminating process...\n")
            self.process.terminate()

    def _reader(self) -> None:
        assert self.process is not None
        assert self.process.stdout is not None
        for line in self.process.stdout:
            self.output_queue.put(line)
        return_code = self.process.wait()
        self.output_queue.put(f"\n[gui] process exited with code {return_code}\n")
        self.output_queue.put(None)

    def _drain(self) -> None:
        done = False
        while True:
            try:
                item = self.output_queue.get_nowait()
            except queue.Empty:
                break
            if item is None:
                done = True
            else:
                self.append(item)
        if done:
            self.on_done()
        elif self.running():
            self.log.after(100, self._drain)


class CyAntsGui:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title(APP_TITLE)
        self.root.geometry("1120x780")
        self.root.minsize(920, 640)
        self._icon_image = None
        self._header_icon = None

        self.runner: CommandRunner | None = None
        self.status = tk.StringVar(value="Ready")

        self._configure_icon()
        self._configure_style()
        self._build_vars()
        self._build_layout()
        self._refresh_command()

    def _configure_icon(self) -> None:
        icon_path = resource_path("assets/cyants_icon_256.png")
        if not icon_path.exists():
            return
        try:
            self._icon_image = tk.PhotoImage(file=str(icon_path))
            self.root.iconphoto(True, self._icon_image)
        except tk.TclError:
            self._icon_image = None

    def _apply_window_icon(self, window: tk.Toplevel) -> None:
        if self._icon_image is None:
            return
        with contextlib.suppress(tk.TclError):
            window.iconphoto(False, self._icon_image)

    def _configure_style(self) -> None:
        self.root.configure(bg=THEME["bg"])
        style = ttk.Style(self.root)
        with contextlib.suppress(tk.TclError):
            style.theme_use("clam")
        style.configure(".", background=THEME["bg"], foreground=THEME["text"], fieldbackground=THEME["field"])
        style.configure("TFrame", background=THEME["bg"])
        style.configure("Panel.TFrame", background=THEME["panel"])
        style.configure("TLabel", background=THEME["bg"], foreground=THEME["text"])
        style.configure("Panel.TLabel", background=THEME["panel"], foreground=THEME["text"])
        style.configure("PanelMuted.TLabel", background=THEME["panel"], foreground=THEME["muted"])
        style.configure("Title.TLabel", background=THEME["panel"], foreground=THEME["accent"], font=("TkDefaultFont", 18, "bold"))
        style.configure("Muted.TLabel", background=THEME["bg"], foreground=THEME["muted"])
        style.configure(
            "TLabelframe",
            background=THEME["bg"],
            foreground=THEME["text"],
            bordercolor=THEME["border"],
            relief="solid",
        )
        style.configure("TLabelframe.Label", background=THEME["bg"], foreground=THEME["accent"])
        style.configure(
            "TButton",
            background=THEME["panel_alt"],
            foreground=THEME["text"],
            bordercolor=THEME["border"],
            focusthickness=2,
            focuscolor=THEME["accent"],
            padding=(10, 5),
        )
        style.map(
            "TButton",
            background=[("active", THEME["accent_blue"]), ("pressed", THEME["accent"])],
            foreground=[("active", "#ffffff"), ("pressed", "#001827")],
        )
        style.configure(
            "Accent.TButton",
            background=THEME["accent"],
            foreground="#001827",
            bordercolor=THEME["accent"],
            padding=(14, 6),
        )
        style.map("Accent.TButton", background=[("active", THEME["accent_blue"]), ("pressed", THEME["accent"])])
        style.configure("TEntry", fieldbackground=THEME["field"], foreground=THEME["field_text"], insertcolor=THEME["field_text"])
        style.configure("TCombobox", fieldbackground=THEME["field"], foreground=THEME["field_text"])
        style.configure("TCheckbutton", background=THEME["bg"], foreground=THEME["text"])
        style.map("TCheckbutton", background=[("active", THEME["bg"])], foreground=[("active", THEME["accent"])])
        style.configure("TNotebook", background=THEME["bg"], borderwidth=0)
        style.configure(
            "TNotebook.Tab",
            background=THEME["panel"],
            foreground=THEME["muted"],
            padding=(14, 7),
        )
        style.map(
            "TNotebook.Tab",
            background=[("selected", THEME["panel_alt"]), ("active", THEME["panel_alt"])],
            foreground=[("selected", THEME["accent"]), ("active", THEME["text"])],
        )
        style.configure("Horizontal.TScale", background=THEME["bg"], troughcolor=THEME["panel_alt"])

    def _build_vars(self) -> None:
        self.common_project_root = tk.StringVar()
        self.common_spacing = tk.StringVar(value=DEFAULT_SPACING)
        self.common_threads = tk.StringVar(value=DEFAULT_THREADS)
        self.common_progress = tk.StringVar(value=DEFAULT_PROGRESS_INTERVAL)
        default_format = COMMAND_FORMAT_WINDOWS if os.name == "nt" else COMMAND_FORMAT_POSIX
        self.command_format = tk.StringVar(value=default_format)

        self.roi_cycle = tk.StringVar(value="cycle_005")
        self.roi_ims = tk.StringVar()
        self.roi_csv = tk.StringVar()
        self.roi_fixed_crop = tk.StringVar()
        self.roi_channels = tk.StringVar(value="0-3")
        self.roi_channel_offset = tk.StringVar(value="0")
        self.roi_ds = tk.StringVar(value="4")
        self.roi_prefix = tk.StringVar(value="cycle_005_roi_registered")
        self.roi_apply_only = tk.BooleanVar(value=False)
        self.roi_open_qc = tk.BooleanVar(value=True)
        self.roi_save_raw = tk.BooleanVar(value=False)

        self.whole_fixed_ims = tk.StringVar()
        self.whole_ims = tk.StringVar()
        self.whole_out = tk.StringVar()
        self.whole_final_out = tk.StringVar()
        self.whole_channels = tk.StringVar(value="0-3")
        self.whole_channel_offset = tk.StringVar(value="0")
        self.whole_source_map = tk.StringVar()
        self.whole_global_tx = tk.StringVar(value="TRSAA")
        self.whole_global_ds = tk.StringVar(value="16")
        self.whole_tiles = tk.StringVar(value="3")
        self.whole_overlap = tk.StringVar(value="0.25")
        self.whole_tile_ds = tk.StringVar(value="6")
        self.whole_syn_tx = tk.StringVar(value="SyNOnly")
        self.whole_blend = tk.StringVar(value="memmap")
        self.whole_ram = tk.StringVar(value="1500")
        self.whole_map_only = tk.BooleanVar(value=False)
        self.whole_global_qc = tk.BooleanVar(value=True)
        self.whole_qc = tk.BooleanVar(value=True)
        self.whole_open_qc = tk.BooleanVar(value=True)
        self.whole_no_local = tk.BooleanVar(value=False)
        self.whole_apply_only = tk.BooleanVar(value=False)

        self.intra_ims = tk.StringVar()
        self.intra_out = tk.StringVar()
        self.intra_ref_ch = tk.StringVar(value="0")
        self.intra_channels = tk.StringVar(value="1,3")
        self.intra_transform = tk.StringVar(value="TRSAA")
        self.intra_ds = tk.StringVar(value="16")
        self.intra_full_tif = tk.BooleanVar(value=True)
        self.intra_apply_only = tk.BooleanVar(value=False)
        self.intra_open_qc = tk.BooleanVar(value=True)

    def _build_layout(self) -> None:
        outer = ttk.Frame(self.root, padding=10)
        outer.pack(fill="both", expand=True)

        header = ttk.Frame(outer, style="Panel.TFrame", padding=(10, 8))
        header.pack(fill="x", pady=(0, 10))
        icon_path = resource_path("assets/cyants_icon_256.png")
        if icon_path.exists():
            try:
                self._header_icon = tk.PhotoImage(file=str(icon_path)).subsample(4, 4)
                ttk.Label(header, image=self._header_icon, style="Panel.TLabel").pack(side="left", padx=(0, 10))
            except tk.TclError:
                self._header_icon = None
        title_block = ttk.Frame(header, style="Panel.TFrame")
        title_block.pack(side="left", fill="x", expand=True)
        ttk.Label(title_block, text="CyANTs", style="Title.TLabel").pack(anchor="w")
        ttk.Label(title_block, text="Registration launcher", style="PanelMuted.TLabel").pack(anchor="w")

        common = ttk.LabelFrame(outer, text="Common")
        common.pack(fill="x", pady=(0, 8))
        self._path_row(common, Field("Project Reg folder", self.common_project_root, "dir"), 0)
        self._entry_row(common, "Spacing", self.common_spacing, 1, width=12, column=1)
        self._entry_row(common, "Threads", self.common_threads, 1, width=8, column=3)
        self._entry_row(common, "Progress interval", self.common_progress, 1, width=8, column=5)

        self.notebook = ttk.Notebook(outer)
        self.notebook.pack(fill="both", expand=True)
        self.notebook.bind("<<NotebookTabChanged>>", lambda _event: self._refresh_command())

        self.roi_tab = ttk.Frame(self.notebook, padding=8)
        self.whole_tab = ttk.Frame(self.notebook, padding=8)
        self.intra_tab = ttk.Frame(self.notebook, padding=8)
        self.notebook.add(self.roi_tab, text="ROI Registration")
        self.notebook.add(self.whole_tab, text="Whole Volume")
        self.notebook.add(self.intra_tab, text="Intracycle")

        self._build_roi_tab()
        self._build_whole_tab()
        self._build_intra_tab()

        command_frame = ttk.LabelFrame(outer, text="Generated command")
        command_frame.pack(fill="x", pady=(8, 8))
        command_options = ttk.Frame(command_frame)
        command_options.pack(fill="x", padx=6, pady=(6, 0))
        ttk.Label(command_options, text="Command format").pack(side="left", padx=(0, 6))
        command_format = ttk.Combobox(
            command_options,
            textvariable=self.command_format,
            values=(COMMAND_FORMAT_WINDOWS, COMMAND_FORMAT_POSIX, COMMAND_FORMAT_LOCAL),
            width=18,
            state="readonly",
        )
        command_format.pack(side="left")
        command_format.bind("<<ComboboxSelected>>", lambda _event: self._refresh_command())
        ttk.Label(
            command_options,
            text="Copy/export can use portable python script names; Run always uses this local GUI environment.",
        ).pack(side="left", padx=(10, 0))
        self.command_text = tk.Text(
            command_frame,
            height=4,
            wrap="word",
            bg=THEME["log_bg"],
            fg=THEME["log_fg"],
            insertbackground=THEME["accent"],
            selectbackground=THEME["border"],
            relief="flat",
            padx=8,
            pady=6,
        )
        self.command_text.pack(fill="x", padx=6, pady=6)

        buttons = ttk.Frame(outer)
        buttons.pack(fill="x", pady=(0, 8))
        ttk.Button(buttons, text="Refresh Command", command=self._refresh_command).pack(side="left")
        ttk.Button(buttons, text="Copy Command", command=self._copy_command).pack(side="left", padx=6)
        ttk.Button(buttons, text="Save Command", command=self._save_command).pack(side="left")
        ttk.Button(buttons, text="Save Config", command=self._save_config).pack(side="left")
        ttk.Button(buttons, text="Load Config", command=self._load_config).pack(side="left", padx=6)
        self.run_button = ttk.Button(buttons, text="Run", command=self._run, style="Accent.TButton")
        self.run_button.pack(side="right")
        self.stop_button = ttk.Button(buttons, text="Stop", command=self._stop)
        self.stop_button.pack(side="right", padx=6)

        log_frame = ttk.LabelFrame(outer, text="Log")
        log_frame.pack(fill="both", expand=True)
        self.log = tk.Text(
            log_frame,
            wrap="word",
            bg=THEME["log_bg"],
            fg=THEME["log_fg"],
            insertbackground=THEME["accent"],
            selectbackground=THEME["border"],
            relief="flat",
            padx=8,
            pady=6,
        )
        self.log.pack(fill="both", expand=True, padx=6, pady=6)
        self.runner = CommandRunner(self.log, self._run_done)

        status = ttk.Label(outer, textvariable=self.status, anchor="w")
        status.pack(fill="x")
        if TkinterDnD is None:
            self.status.set("Ready. Install tkinterdnd2 for drag/drop; Browse buttons work without it.")

    def _build_roi_tab(self) -> None:
        rows = [
            Field("Cycle id", self.roi_cycle),
            Field("Moving .ims", self.roi_ims, "file"),
            Field("ROI CSV", self.roi_csv, "file"),
            Field("Fixed reference crop", self.roi_fixed_crop, "file"),
            Field("Channels", self.roi_channels),
            Field("Channel offset", self.roi_channel_offset),
            Field("Downsample", self.roi_ds),
            Field("TIFF prefix", self.roi_prefix),
        ]
        for i, field in enumerate(rows):
            if field.browse:
                self._path_row(self.roi_tab, field, i)
            else:
                self._entry_row(self.roi_tab, field.label, field.var, i)
        self._check_row(
            self.roi_tab,
            len(rows),
            (("Apply only", self.roi_apply_only), ("Open QC", self.roi_open_qc), ("Save raw ROI", self.roi_save_raw)),
        )
        preview = ttk.Frame(self.roi_tab)
        preview.grid(row=len(rows) + 1, column=0, columnspan=3, sticky="w", pady=(4, 3))
        ttk.Button(preview, text="Preview ROI Fast", command=self._preview_roi).pack(side="left")
        ttk.Label(preview, text="Loads the coarsest .ims pyramid level and overlays the ROI CSV.").pack(side="left", padx=(10, 0))

    def _build_whole_tab(self) -> None:
        rows = [
            Field("Fixed Cycle0 .ims", self.whole_fixed_ims, "file"),
            Field("Moving cycle .ims", self.whole_ims, "file"),
            Field("Cycle output folder", self.whole_out, "dir"),
            Field("Final TIFF folder", self.whole_final_out, "dir"),
            Field("Channels", self.whole_channels),
            Field("Channel offset", self.whole_channel_offset),
            Field("Source map", self.whole_source_map),
            Field("Global transform", self.whole_global_tx),
            Field("Global downsample", self.whole_global_ds),
            Field("Tiles", self.whole_tiles),
            Field("Overlap", self.whole_overlap),
            Field("Tile downsample", self.whole_tile_ds),
            Field("Tile transform", self.whole_syn_tx),
            Field("Blend", self.whole_blend),
            Field("RAM limit GB", self.whole_ram),
        ]
        for i, field in enumerate(rows):
            if field.browse:
                self._path_row(self.whole_tab, field, i)
            else:
                self._entry_row(self.whole_tab, field.label, field.var, i)
        self._check_row(
            self.whole_tab,
            len(rows),
            (
                ("Map only", self.whole_map_only),
                ("Global QC", self.whole_global_qc),
                ("Stitched QC", self.whole_qc),
                ("Open QC", self.whole_open_qc),
                ("No local refine", self.whole_no_local),
                ("Apply only", self.whole_apply_only),
            ),
        )
        preview = ttk.Frame(self.whole_tab)
        preview.grid(row=len(rows) + 1, column=0, columnspan=3, sticky="w", pady=(4, 3))
        ttk.Button(preview, text="Preview Tiles", command=self._preview_tiles).pack(side="left")
        ttk.Label(preview, text="Draws the auto full-Z XY tile grid from .ims metadata.").pack(side="left", padx=(10, 0))

    def _build_intra_tab(self) -> None:
        rows = [
            Field("Source .ims", self.intra_ims, "file"),
            Field("Output folder", self.intra_out, "dir"),
            Field("Reference channel", self.intra_ref_ch),
            Field("Channels to align", self.intra_channels),
            Field("Transform", self.intra_transform),
            Field("Downsample", self.intra_ds),
        ]
        for i, field in enumerate(rows):
            if field.browse:
                self._path_row(self.intra_tab, field, i)
            else:
                self._entry_row(self.intra_tab, field.label, field.var, i)
        self._check_row(
            self.intra_tab,
            len(rows),
            (("Full-res TIFF", self.intra_full_tif), ("Apply only", self.intra_apply_only), ("Open QC", self.intra_open_qc)),
        )

    def _entry_row(self, parent, label: str, var: tk.StringVar, row: int, width: int = 72, column: int = 0) -> ttk.Entry:
        ttk.Label(parent, text=label).grid(row=row, column=column, sticky="w", padx=(0, 6), pady=3)
        entry = ttk.Entry(parent, textvariable=var, width=width)
        entry.grid(row=row, column=column + 1, sticky="ew", padx=(0, 6), pady=3)
        parent.columnconfigure(column + 1, weight=1)
        var.trace_add("write", lambda *_args: self._refresh_command())
        self._enable_drop(entry, var)
        return entry

    def _path_row(self, parent, field: Field, row: int) -> None:
        entry = self._entry_row(parent, field.label, field.var, row, field.width)
        ttk.Button(parent, text="Browse", command=lambda: self._browse(field)).grid(row=row, column=2, sticky="e", pady=3)
        self._enable_drop(entry, field.var)

    def _check_row(self, parent, row: int, checks: tuple[tuple[str, tk.BooleanVar], ...]) -> None:
        frame = ttk.Frame(parent)
        frame.grid(row=row, column=0, columnspan=3, sticky="w", pady=(8, 3))
        for label, var in checks:
            ttk.Checkbutton(frame, text=label, variable=var, command=self._refresh_command).pack(side="left", padx=(0, 12))

    def _enable_drop(self, widget: ttk.Entry, var: tk.StringVar) -> None:
        if DND_FILES is None:
            return
        widget.drop_target_register(DND_FILES)
        widget.dnd_bind("<<Drop>>", lambda event: self._drop_path(event, var))

    def _drop_path(self, event, var: tk.StringVar) -> None:
        value = event.data.strip()
        if value.startswith("{") and value.endswith("}"):
            value = value[1:-1]
        var.set(value)

    def _browse(self, field: Field) -> None:
        if field.browse == "dir":
            value = filedialog.askdirectory()
        else:
            value = filedialog.askopenfilename()
        if value:
            field.var.set(value)

    def _active_mode(self) -> str:
        return self.notebook.tab(self.notebook.select(), "text")

    def _base_common(self) -> list[str]:
        return ["--spacing", self.common_spacing.get(), "--t", self.common_threads.get(), "--pi", self.common_progress.get()]

    def _roi_command(self) -> list[str]:
        args = [
            "--proj",
            self.common_project_root.get(),
            "--cycle",
            self.roi_cycle.get(),
            "--ims",
            self.roi_ims.get(),
            "--roi-csv",
            self.roi_csv.get(),
            "--fixed-crop",
            self.roi_fixed_crop.get(),
            "--ch",
            self.roi_channels.get(),
            "--co",
            self.roi_channel_offset.get(),
            "--ds",
            self.roi_ds.get(),
            "--reg-format",
            "tif",
            "--u16",
            "clip",
            "--tif-prefix",
            self.roi_prefix.get(),
        ]
        args.extend(self._base_common())
        args.extend(parse_bool(self.roi_open_qc, "--open-qc"))
        args.extend(parse_bool(self.roi_apply_only, "--ao"))
        args.extend(parse_bool(self.roi_save_raw, "--save-raw-roi"))
        return local_script_command("vcycif_ROI.py", args)

    def _whole_command(self) -> list[str]:
        args = [
            "--fixed-ims",
            self.whole_fixed_ims.get(),
            "--ims",
            self.whole_ims.get(),
            "--ref-ch",
            "0",
            "--reg-ch",
            "0",
            "--ch",
            self.whole_channels.get(),
            "--out",
            self.whole_out.get(),
            "--co",
            self.whole_channel_offset.get(),
            "--global-tx",
            self.whole_global_tx.get(),
            "--global-ds",
            self.whole_global_ds.get(),
            "--tiles",
            self.whole_tiles.get(),
            "--overlap",
            self.whole_overlap.get(),
            "--tile-ds",
            self.whole_tile_ds.get(),
            "--mxy",
            "512",
            "--mz",
            "16",
            "--syn-tx",
            self.whole_syn_tx.get(),
            "--blend",
            self.whole_blend.get(),
            "--ram-limit-gb",
            self.whole_ram.get(),
        ]
        if self.whole_final_out.get().strip():
            args.extend(["--final-out", self.whole_final_out.get()])
        if self.whole_source_map.get().strip():
            args.extend(["--source-map", self.whole_source_map.get()])
        args.extend(self._base_common())
        args.extend(parse_bool(self.whole_map_only, "--map"))
        args.extend(parse_bool(self.whole_global_qc, "--global-qc"))
        args.extend(parse_bool(self.whole_qc, "--qc"))
        args.extend(parse_bool(self.whole_open_qc, "--open-qc"))
        args.extend(parse_bool(self.whole_no_local, "--no-local-refine"))
        args.extend(parse_bool(self.whole_apply_only, "--ao"))
        return local_script_command("ants_ims_tiled_quicksyn.py", args)

    def _intra_command(self) -> list[str]:
        args = [
            "--ims",
            self.intra_ims.get(),
            "--out",
            self.intra_out.get(),
            "--ref-ch",
            self.intra_ref_ch.get(),
            "--ch",
            self.intra_channels.get(),
            "--tx",
            self.intra_transform.get(),
            "--ds",
            self.intra_ds.get(),
        ]
        args.extend(self._base_common())
        args.extend(parse_bool(self.intra_full_tif, "--full-tif"))
        args.extend(parse_bool(self.intra_apply_only, "--ao"))
        args.extend(parse_bool(self.intra_open_qc, "--open-qc"))
        return local_script_command("ants_ims_intracycle_whole.py", args)

    def _load_ims_spec(self, path_value: str, channel: int = 0):
        from cyants_io import build_ims_volume_spec

        path = Path(path_value).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f".ims file does not exist: {path}")
        spec = build_ims_volume_spec(
            path,
            res_level=0,
            timepoint=0,
            channel=channel,
            dataset_path="",
            axis_order="auto",
        )
        spacing = self.common_spacing.get().strip()
        if spacing:
            from cyants_io import parse_spacing_override

            spacing_xyz = parse_spacing_override(spacing)
            if spacing_xyz is not None:
                spec = type(spec)(
                    path=spec.path,
                    dataset_key=spec.dataset_key,
                    axis_order=spec.axis_order,
                    shape_zyx=spec.shape_zyx,
                    spacing_xyz=spacing_xyz,
                )
        return spec

    def _preview_roi(self) -> None:
        ims_path = self.roi_ims.get().strip()
        roi_path = self.roi_csv.get().strip()
        if not ims_path:
            messagebox.showinfo(APP_TITLE, "Choose a moving .ims first.")
            return
        try:
            from PIL import Image, ImageTk
            import numpy as np
            import SimpleITK as sitk
            from cyants_io import read_ims_downsampled_from_spec, read_ims_lowest_resolution_from_spec

            self.status.set("Loading fastest ROI preview...")
            self.root.update_idletasks()
            spec = self._load_ims_spec(ims_path, channel=0)
            preview_source = "16x direct stride fallback"
            try:
                preview_img, level, effective_scale_zyx = read_ims_lowest_resolution_from_spec(spec)
                preview_source = (
                    f"ResolutionLevel {level}, effective zyx scale "
                    f"({effective_scale_zyx[0]:.3g}, {effective_scale_zyx[1]:.3g}, {effective_scale_zyx[2]:.3g})"
                )
            except Exception as exc:
                print(f"[gui-preview] falling back to 16x direct ROI preview: {exc}", flush=True)
                preview_img = read_ims_downsampled_from_spec(spec, (16, 16, 16))
            arr = sitk.GetArrayFromImage(preview_img)
            lo, hi = np.percentile(arr, [1, 99.7])
            if hi <= lo:
                lo, hi = float(np.min(arr)), float(np.max(arr) or 1.0)
            _full_z, full_y, full_x = spec.shape_zyx
            preview_y, preview_x = arr.shape[1], arr.shape[2]
            points: list[tuple[float, float]] = []
            if roi_path and Path(roi_path).expanduser().exists():
                points = parse_roi_csv_points(Path(roi_path).expanduser())
            scaled = [(x * preview_x / full_x, y * preview_y / full_y) for x, y in points]
            roi_note = f"ROI points={len(points)}" if points else "draw or resize ROI, then Confirm ROI"
            self._show_roi_stack_viewer(
                arr=arr,
                intensity_range=(float(lo), float(hi)),
                roi_points=scaled,
                full_shape_zyx=spec.shape_zyx,
                csv_path=roi_path,
                caption=(
                    f"{Path(ims_path).name} ch0 preview from {preview_source}\n"
                    f"full shape zyx={spec.shape_zyx}; preview shape zyx={arr.shape}; {roi_note}"
                ),
                image_module=Image,
                image_tk_module=ImageTk,
            )
            self.status.set("ROI preview ready")
        except Exception as exc:
            self.status.set("Ready")
            messagebox.showerror(APP_TITLE, f"Could not build ROI preview:\n{exc}")

    def _show_image_popup(self, title: str, image, caption: str, image_tk_module) -> None:
        popup = tk.Toplevel(self.root)
        popup.title(title)
        self._apply_window_icon(popup)
        frame = ttk.Frame(popup, padding=8)
        frame.pack(fill="both", expand=True)
        photo = image_tk_module.PhotoImage(image)
        label = ttk.Label(frame, image=photo)
        label.image = photo
        label.pack(fill="both", expand=True)
        ttk.Label(frame, text=caption, justify="left").pack(fill="x", pady=(8, 0))

    def _show_roi_stack_viewer(
        self,
        arr,
        intensity_range: tuple[float, float],
        roi_points: list[tuple[float, float]],
        full_shape_zyx: tuple[int, int, int],
        csv_path: str,
        caption: str,
        image_module,
        image_tk_module,
    ) -> None:
        import numpy as np

        popup = tk.Toplevel(self.root)
        popup.title("ROI Preview Fast")
        self._apply_window_icon(popup)
        frame = ttk.Frame(popup, padding=8)
        frame.pack(fill="both", expand=True)

        state = {
            "photo": None,
            "max_projection": tk.BooleanVar(value=True),
            "z": tk.IntVar(value=max(0, arr.shape[0] // 2)),
            "display_scale": 1.0,
            "roi": None,
            "drag": None,
        }
        preview_h, preview_w = int(arr.shape[1]), int(arr.shape[2])
        if roi_points:
            xs = [point[0] for point in roi_points]
            ys = [point[1] for point in roi_points]
            state["roi"] = [
                max(0.0, min(xs)),
                max(0.0, min(ys)),
                min(float(preview_w - 1), max(xs)),
                min(float(preview_h - 1), max(ys)),
            ]
        else:
            margin_x = preview_w * 0.30
            margin_y = preview_h * 0.30
            state["roi"] = [margin_x, margin_y, preview_w - margin_x, preview_h - margin_y]

        canvas = tk.Canvas(frame, bg="#0f172a", highlightthickness=1, highlightbackground="#94a3b8")
        canvas.pack(fill="both", expand=True)

        controls = ttk.Frame(frame)
        controls.pack(fill="x", pady=(8, 0))
        ttk.Checkbutton(controls, text="Max projection", variable=state["max_projection"]).pack(side="left")
        ttk.Label(controls, text="Z").pack(side="left", padx=(14, 4))
        z_slider = ttk.Scale(controls, from_=0, to=max(0, arr.shape[0] - 1), orient="horizontal")
        z_slider.pack(side="left", fill="x", expand=True)
        z_label = ttk.Label(controls, width=18)
        z_label.pack(side="left", padx=(8, 0))
        ttk.Button(controls, text="Confirm ROI", command=lambda: confirm_roi()).pack(side="left", padx=(12, 0))
        caption_label = ttk.Label(frame, text=caption, justify="left")
        caption_label.pack(fill="x", pady=(8, 0))
        hint_label = ttk.Label(
            frame,
            text="Drag outside the ROI to draw a new box. Drag inside to move it. Drag corners to resize.",
            justify="left",
        )
        hint_label.pack(fill="x", pady=(4, 0))

        lo, hi = intensity_range
        max_w, max_h = 1050, 720

        def normalize_roi(roi: list[float]) -> list[float]:
            x0, y0, x1, y1 = roi
            x0, x1 = sorted((max(0.0, min(float(preview_w - 1), x0)), max(0.0, min(float(preview_w - 1), x1))))
            y0, y1 = sorted((max(0.0, min(float(preview_h - 1), y0)), max(0.0, min(float(preview_h - 1), y1))))
            return [x0, y0, x1, y1]

        def roi_to_canvas(roi: list[float]) -> tuple[float, float, float, float]:
            scale = float(state["display_scale"])
            x0, y0, x1, y1 = normalize_roi(roi)
            return x0 * scale, y0 * scale, x1 * scale, y1 * scale

        def canvas_to_preview(x_value: float, y_value: float) -> tuple[float, float]:
            scale = max(float(state["display_scale"]), 1e-9)
            x = max(0.0, min(float(preview_w - 1), x_value / scale))
            y = max(0.0, min(float(preview_h - 1), y_value / scale))
            return x, y

        def draw_roi_overlay() -> None:
            canvas.delete("roi")
            roi = state["roi"]
            if roi is None:
                return
            x0, y0, x1, y1 = roi_to_canvas(roi)
            canvas.create_rectangle(x0, y0, x1, y1, outline="#ff4614", width=3, tags="roi")
            canvas.create_rectangle(x0, y0, x1, y1, outline="#ffffff", width=1, tags="roi")
            handle_size = 8
            for hx, hy in ((x0, y0), (x1, y0), (x0, y1), (x1, y1)):
                canvas.create_rectangle(
                    hx - handle_size / 2,
                    hy - handle_size / 2,
                    hx + handle_size / 2,
                    hy + handle_size / 2,
                    fill="#ff4614",
                    outline="#ffffff",
                    tags="roi",
                )
            canvas.tag_raise("roi")

        def render() -> None:
            z_index = int(round(float(z_slider.get())))
            z_index = max(0, min(z_index, arr.shape[0] - 1))
            state["z"].set(z_index)
            if state["max_projection"].get():
                plane = np.max(arr, axis=0).astype(np.float32, copy=False)
                z_label.configure(text=f"max projection ({arr.shape[0]} z)")
            else:
                plane = arr[z_index].astype(np.float32, copy=False)
                z_label.configure(text=f"{z_index + 1}/{arr.shape[0]}")
            gray = np.clip((plane - lo) / max(hi - lo, 1e-6) * 255.0, 0, 255).astype(np.uint8)
            image = image_module.fromarray(gray, mode="L").convert("RGB")
            scale = min(max_w / image.width, max_h / image.height, 1.0)
            display = image
            if scale < 1.0:
                resampling = getattr(getattr(image_module, "Resampling", image_module), "BILINEAR")
                display = image.resize((int(image.width * scale), int(image.height * scale)), resampling)
            state["display_scale"] = scale
            canvas.configure(width=display.width, height=display.height, scrollregion=(0, 0, display.width, display.height))
            photo = image_tk_module.PhotoImage(display)
            state["photo"] = photo
            canvas.delete("preview")
            canvas.create_image(0, 0, image=photo, anchor="nw", tags="preview")
            canvas.tag_lower("preview")
            draw_roi_overlay()

        def on_slider(value: str) -> None:
            if not state["max_projection"].get():
                render()
            else:
                z_label.configure(text=f"max projection ({arr.shape[0]} z)")

        def on_toggle() -> None:
            z_slider.configure(state="disabled" if state["max_projection"].get() else "normal")
            render()

        def hit_test(x_value: float, y_value: float) -> str:
            roi = state["roi"]
            if roi is None:
                return "draw"
            x0, y0, x1, y1 = normalize_roi(roi)
            handle_radius = max(8.0 / max(float(state["display_scale"]), 1e-9), 4.0)
            corners = {
                "resize-nw": (x0, y0),
                "resize-ne": (x1, y0),
                "resize-sw": (x0, y1),
                "resize-se": (x1, y1),
            }
            for mode, (hx, hy) in corners.items():
                if abs(x_value - hx) <= handle_radius and abs(y_value - hy) <= handle_radius:
                    return mode
            if x0 <= x_value <= x1 and y0 <= y_value <= y1:
                return "move"
            return "draw"

        def on_press(event) -> None:
            x_value, y_value = canvas_to_preview(event.x, event.y)
            mode = hit_test(x_value, y_value)
            state["drag"] = {
                "mode": mode,
                "start": (x_value, y_value),
                "roi": list(state["roi"] or [x_value, y_value, x_value, y_value]),
            }
            if mode == "draw":
                state["roi"] = [x_value, y_value, x_value, y_value]
                draw_roi_overlay()

        def on_motion(event) -> None:
            drag = state.get("drag")
            if not drag:
                return
            x_value, y_value = canvas_to_preview(event.x, event.y)
            mode = drag["mode"]
            start_x, start_y = drag["start"]
            original = list(drag["roi"])
            if mode == "move":
                dx = x_value - start_x
                dy = y_value - start_y
                width = original[2] - original[0]
                height = original[3] - original[1]
                x0 = max(0.0, min(float(preview_w - 1) - width, original[0] + dx))
                y0 = max(0.0, min(float(preview_h - 1) - height, original[1] + dy))
                state["roi"] = [x0, y0, x0 + width, y0 + height]
            elif mode == "resize-nw":
                state["roi"] = [x_value, y_value, original[2], original[3]]
            elif mode == "resize-ne":
                state["roi"] = [original[0], y_value, x_value, original[3]]
            elif mode == "resize-sw":
                state["roi"] = [x_value, original[1], original[2], y_value]
            elif mode == "resize-se":
                state["roi"] = [original[0], original[1], x_value, y_value]
            else:
                state["roi"] = [start_x, start_y, x_value, y_value]
            state["roi"] = normalize_roi(state["roi"])
            draw_roi_overlay()

        def on_release(_event) -> None:
            if state.get("roi") is not None:
                state["roi"] = normalize_roi(state["roi"])
                draw_roi_overlay()
            state["drag"] = None

        def confirm_roi() -> None:
            roi = state["roi"]
            if roi is None:
                messagebox.showinfo(APP_TITLE, "Draw an ROI first.")
                return
            x0, y0, x1, y1 = normalize_roi(roi)
            if (x1 - x0) < 1 or (y1 - y0) < 1:
                messagebox.showinfo(APP_TITLE, "Draw a larger ROI before confirming.")
                return
            full_points = preview_roi_to_full_points(
                (x0, y0, x1, y1),
                preview_shape_yx=(preview_h, preview_w),
                full_shape_zyx=full_shape_zyx,
            )
            output_path: Path | None = None
            if csv_path.strip():
                output_path = Path(csv_path).expanduser()
            else:
                chosen = filedialog.asksaveasfilename(
                    title="Save ROI coordinates",
                    defaultextension=".csv",
                    filetypes=[("CSV", "*.csv"), ("All files", "*.*")],
                )
                if chosen:
                    output_path = Path(chosen)
            if output_path is None:
                return
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["X", "Y"])
                writer.writerows(full_points)
            self.roi_csv.set(str(output_path))
            self.status.set(f"Saved ROI CSV: {output_path}")
            caption_label.configure(
                text=(
                    f"{caption}\nSaved ROI CSV: {output_path}\n"
                    f"full-resolution X/Y points: {full_points}"
                )
            )

        canvas.bind("<ButtonPress-1>", on_press)
        canvas.bind("<B1-Motion>", on_motion)
        canvas.bind("<ButtonRelease-1>", on_release)
        z_slider.set(state["z"].get())
        z_slider.configure(command=on_slider)
        state["max_projection"].trace_add("write", lambda *_args: on_toggle())
        on_toggle()

    def _preview_tiles(self) -> None:
        ims_path = self.whole_fixed_ims.get().strip() or self.whole_ims.get().strip()
        if not ims_path:
            messagebox.showinfo(APP_TITLE, "Choose a fixed or moving .ims first.")
            return
        try:
            spec = self._load_ims_spec(ims_path, channel=0)
            tile_count = int(self.whole_tiles.get())
            if tile_count < 1:
                raise ValueError("Tiles must be >= 1 for GUI preview.")
            overlap_fraction = float(self.whole_overlap.get())
            if not (0.0 <= overlap_fraction < 1.0):
                raise ValueError("GUI tile preview expects fractional overlap, for example 0.25.")
            tile_size, overlap, tiles, grid = compute_auto_tiles(spec.shape_zyx, tile_count, overlap_fraction)
            self._show_tile_popup(spec, tile_size, overlap, tiles, grid)
            self.status.set("Tile preview ready")
        except Exception as exc:
            self.status.set("Ready")
            messagebox.showerror(APP_TITLE, f"Could not build tile preview:\n{exc}")

    def _show_tile_popup(
        self,
        spec,
        tile_size: tuple[int, int, int],
        overlap: tuple[int, int, int],
        tiles: list[tuple[int, int, int, int, int, int]],
        grid: tuple[int, int],
    ) -> None:
        z_size, y_size, x_size = spec.shape_zyx
        max_w, max_h = 1050.0, 680.0
        scale = min(max_w / x_size, max_h / y_size)
        canvas_w = int(x_size * scale) + 60
        canvas_h = int(y_size * scale) + 105
        popup = tk.Toplevel(self.root)
        popup.title("Tile Preview")
        self._apply_window_icon(popup)
        frame = ttk.Frame(popup, padding=8)
        frame.pack(fill="both", expand=True)
        canvas = tk.Canvas(
            frame,
            width=canvas_w,
            height=canvas_h,
            bg=THEME["bg"],
            highlightthickness=1,
            highlightbackground=THEME["border"],
        )
        canvas.pack(fill="both", expand=True)
        ox, oy = 30, 28
        image_w = x_size * scale
        image_h = y_size * scale
        canvas.create_rectangle(ox, oy, ox + image_w, oy + image_h, fill=THEME["panel"], outline=THEME["accent"], width=2)
        colors = ["#26d6d0", "#1aa7ff", "#f72585", "#7df9c6", "#ffd166", "#9b5de5", "#00f5d4", "#ff8fab"]
        for index, (_z0, _z1, y0, y1, x0, x1) in enumerate(tiles, start=1):
            color = colors[(index - 1) % len(colors)]
            x = ox + x0 * scale
            y = oy + y0 * scale
            w = (x1 - x0) * scale
            h = (y1 - y0) * scale
            canvas.create_rectangle(x, y, x + w, y + h, outline=color, width=3)
            canvas.create_text(x + w / 2, y + h / 2, text=str(index), fill=color, font=("TkDefaultFont", 14, "bold"))
        caption = (
            f"source={Path(spec.path).name} | shape zyx={spec.shape_zyx} | grid yx={grid} | "
            f"tile zyx={tile_size} | overlap zyx={overlap}"
        )
        canvas.create_text(ox, oy + image_h + 26, text=caption, anchor="w", fill=THEME["text"])
        canvas.create_text(
            ox,
            oy + image_h + 50,
            text="Preview uses metadata only. Tiles span full Z; rectangles show XY extents and overlap.",
            anchor="w",
            fill=THEME["muted"],
        )

    def _current_command(self) -> list[str]:
        mode = self._active_mode()
        if mode == "ROI Registration":
            return self._roi_command()
        if mode == "Whole Volume":
            return self._whole_command()
        return self._intra_command()

    def _display_command(self) -> list[str]:
        command = self._current_command()
        if self.command_format.get() == COMMAND_FORMAT_LOCAL:
            return command
        return portable_command(command)

    def _display_command_text(self) -> str:
        command_format = self.command_format.get()
        command = self._display_command()
        return quote_command(command, command_format)

    def _refresh_command(self) -> None:
        if not hasattr(self, "command_text"):
            return
        self.command_text.delete("1.0", "end")
        self.command_text.insert("1.0", self._display_command_text())

    def _copy_command(self) -> None:
        self._refresh_command()
        text = self.command_text.get("1.0", "end").strip()
        self.root.clipboard_clear()
        self.root.clipboard_append(text)
        self.status.set("Command copied to clipboard")

    def _save_command(self) -> None:
        self._refresh_command()
        default_ext = ".bat" if self.command_format.get() == COMMAND_FORMAT_WINDOWS else ".sh"
        path = filedialog.asksaveasfilename(
            defaultextension=default_ext,
            filetypes=[("Shell scripts", "*.sh"), ("Windows batch", "*.bat"), ("Text", "*.txt"), ("All files", "*.*")],
        )
        if not path:
            return
        command_text = self.command_text.get("1.0", "end").strip()
        Path(path).write_text(command_text + "\n", encoding="utf-8")
        self.status.set(f"Saved command: {path}")

    def _config(self) -> dict[str, object]:
        return {
            "mode": self._active_mode(),
            "common": {
                "project_root": self.common_project_root.get(),
                "spacing": self.common_spacing.get(),
                "threads": self.common_threads.get(),
                "progress_interval": self.common_progress.get(),
                "command_format": self.command_format.get(),
            },
            "roi": {
                "cycle": self.roi_cycle.get(),
                "ims": self.roi_ims.get(),
                "roi_csv": self.roi_csv.get(),
                "fixed_crop": self.roi_fixed_crop.get(),
                "channels": self.roi_channels.get(),
                "channel_offset": self.roi_channel_offset.get(),
                "downsample": self.roi_ds.get(),
                "prefix": self.roi_prefix.get(),
                "apply_only": self.roi_apply_only.get(),
                "open_qc": self.roi_open_qc.get(),
                "save_raw": self.roi_save_raw.get(),
            },
            "whole": {
                "fixed_ims": self.whole_fixed_ims.get(),
                "ims": self.whole_ims.get(),
                "out": self.whole_out.get(),
                "final_out": self.whole_final_out.get(),
                "channels": self.whole_channels.get(),
                "channel_offset": self.whole_channel_offset.get(),
                "source_map": self.whole_source_map.get(),
                "global_tx": self.whole_global_tx.get(),
                "global_ds": self.whole_global_ds.get(),
                "tiles": self.whole_tiles.get(),
                "overlap": self.whole_overlap.get(),
                "tile_ds": self.whole_tile_ds.get(),
                "syn_tx": self.whole_syn_tx.get(),
                "blend": self.whole_blend.get(),
                "ram": self.whole_ram.get(),
                "map_only": self.whole_map_only.get(),
                "global_qc": self.whole_global_qc.get(),
                "qc": self.whole_qc.get(),
                "open_qc": self.whole_open_qc.get(),
                "no_local_refine": self.whole_no_local.get(),
                "apply_only": self.whole_apply_only.get(),
            },
            "intracycle": {
                "ims": self.intra_ims.get(),
                "out": self.intra_out.get(),
                "ref_ch": self.intra_ref_ch.get(),
                "channels": self.intra_channels.get(),
                "transform": self.intra_transform.get(),
                "downsample": self.intra_ds.get(),
                "full_tif": self.intra_full_tif.get(),
                "apply_only": self.intra_apply_only.get(),
                "open_qc": self.intra_open_qc.get(),
            },
        }

    def _save_config(self) -> None:
        path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON", "*.json"), ("All files", "*.*")])
        if not path:
            return
        Path(path).write_text(json.dumps(self._config(), indent=2), encoding="utf-8")
        self.status.set(f"Saved config: {path}")

    def _load_config(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("JSON", "*.json"), ("All files", "*.*")])
        if not path:
            return
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        self._apply_config(data)
        self.status.set(f"Loaded config: {path}")
        self._refresh_command()

    def _set_many(self, mapping: dict[str, object], setters: dict[str, tk.Variable]) -> None:
        for key, var in setters.items():
            if key in mapping:
                var.set(mapping[key])

    def _apply_config(self, data: dict[str, object]) -> None:
        common = data.get("common", {})
        if isinstance(common, dict):
            self._set_many(
                common,
                {
                    "project_root": self.common_project_root,
                    "spacing": self.common_spacing,
                    "threads": self.common_threads,
                    "progress_interval": self.common_progress,
                    "command_format": self.command_format,
                },
            )
        roi = data.get("roi", {})
        if isinstance(roi, dict):
            self._set_many(
                roi,
                {
                    "cycle": self.roi_cycle,
                    "ims": self.roi_ims,
                    "roi_csv": self.roi_csv,
                    "fixed_crop": self.roi_fixed_crop,
                    "channels": self.roi_channels,
                    "channel_offset": self.roi_channel_offset,
                    "downsample": self.roi_ds,
                    "prefix": self.roi_prefix,
                    "apply_only": self.roi_apply_only,
                    "open_qc": self.roi_open_qc,
                    "save_raw": self.roi_save_raw,
                },
            )
        whole = data.get("whole", {})
        if isinstance(whole, dict):
            self._set_many(
                whole,
                {
                    "fixed_ims": self.whole_fixed_ims,
                    "ims": self.whole_ims,
                    "out": self.whole_out,
                    "final_out": self.whole_final_out,
                    "channels": self.whole_channels,
                    "channel_offset": self.whole_channel_offset,
                    "source_map": self.whole_source_map,
                    "global_tx": self.whole_global_tx,
                    "global_ds": self.whole_global_ds,
                    "tiles": self.whole_tiles,
                    "overlap": self.whole_overlap,
                    "tile_ds": self.whole_tile_ds,
                    "syn_tx": self.whole_syn_tx,
                    "blend": self.whole_blend,
                    "ram": self.whole_ram,
                    "map_only": self.whole_map_only,
                    "global_qc": self.whole_global_qc,
                    "qc": self.whole_qc,
                    "open_qc": self.whole_open_qc,
                    "no_local_refine": self.whole_no_local,
                    "apply_only": self.whole_apply_only,
                },
            )
        intracycle = data.get("intracycle", {})
        if isinstance(intracycle, dict):
            self._set_many(
                intracycle,
                {
                    "ims": self.intra_ims,
                    "out": self.intra_out,
                    "ref_ch": self.intra_ref_ch,
                    "channels": self.intra_channels,
                    "transform": self.intra_transform,
                    "downsample": self.intra_ds,
                    "full_tif": self.intra_full_tif,
                    "apply_only": self.intra_apply_only,
                    "open_qc": self.intra_open_qc,
                },
            )
        mode = data.get("mode")
        if isinstance(mode, str):
            for index in range(self.notebook.index("end")):
                if self.notebook.tab(index, "text") == mode:
                    self.notebook.select(index)
                    break

    def _run(self) -> None:
        self._refresh_command()
        if self.runner is None:
            return
        self.status.set("Running")
        self.run_button.configure(state="disabled")
        self.runner.start(self._current_command())

    def _stop(self) -> None:
        if self.runner is not None:
            self.runner.stop()

    def _run_done(self) -> None:
        self.run_button.configure(state="normal")
        self.status.set("Ready")


def run_bundled_script(argv: list[str]) -> int:
    if len(argv) < 2 or argv[0] != RUNNER_FLAG:
        raise ValueError(f"Expected {RUNNER_FLAG} <script.py> [args...]")
    script_name = Path(argv[1]).name
    if script_name not in RUNNABLE_SCRIPTS:
        raise ValueError(f"Unsupported bundled script: {script_name}")
    script_path = Path(repo_script(script_name))
    if not script_path.exists():
        raise FileNotFoundError(f"Bundled script not found: {script_path}")
    old_argv = sys.argv[:]
    old_path = sys.path[:]
    sys.argv = [str(script_path), *argv[2:]]
    if str(script_path.parent) not in sys.path:
        sys.path.insert(0, str(script_path.parent))
    module_globals = {"__name__": "__main__", "__file__": str(script_path)}
    try:
        code = compile(script_path.read_text(encoding="utf-8"), str(script_path), "exec")
        with contextlib.redirect_stdout(sys.stdout), contextlib.redirect_stderr(sys.stderr):
            exec(code, module_globals)
    except SystemExit as exc:
        if isinstance(exc.code, int):
            return exc.code
        if exc.code is None:
            return 0
        print(exc.code, file=sys.stderr)
        return 1
    finally:
        sys.argv = old_argv
        sys.path = old_path
    return 0


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == RUNNER_FLAG:
        return run_bundled_script(sys.argv[1:])
    if TkinterDnD is not None:
        root = TkinterDnD.Tk()
    else:
        root = tk.Tk()
    CyAntsGui(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
