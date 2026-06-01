# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller build spec for the optional CyANTs Windows GUI executable."""

from pathlib import Path

from PyInstaller.utils.hooks import collect_submodules


ROOT = Path(SPECPATH).resolve().parents[1]
PIPELINE_SCRIPTS = [
    "vcycif_ROI.py",
    "ants_ims_tiled_quicksyn.py",
    "ants_ims_intracycle_whole.py",
    "ants_roi_quicksyn.py",
    "apply_ims_roi_channels.py",
    "extract_ims_roi_channels.py",
    "cyants_io.py",
    "pad_tiff_stack.py",
]

datas = [(str(ROOT / script), ".") for script in PIPELINE_SCRIPTS]
datas += [
    (str(ROOT / "assets" / "cyants_icon.png"), "assets"),
    (str(ROOT / "assets" / "cyants_icon_256.png"), "assets"),
]
hiddenimports = []
for package in (
    "ants",
    "h5py",
    "numpy",
    "PIL",
    "roifile",
    "SimpleITK",
    "tifffile",
):
    hiddenimports += collect_submodules(package)


a = Analysis(
    ["cyants_gui.py"],
    pathex=[str(ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="CyANTs",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=str(ROOT / "packaging" / "windows" / "assets" / "cyants_icon.ico"),
)
