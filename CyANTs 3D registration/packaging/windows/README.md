# CyANTs Windows Packaging

This folder contains the optional Windows GUI executable and installer build assets.

## Install The Provided Download

Use `CyANTs_Setup.exe` when available. It is the from-scratch Windows installer for users who do not already have Python or conda. It installs CyANTs for the current user under `%LOCALAPPDATA%\Programs\CyANTs`, bootstraps Miniforge if conda is missing, creates or repairs the `cyants` conda environment, then launches the GUI from that environment.

1. Download the provided `CyANTs_Setup.exe` from GitHub Actions, a release, or the lab share.
2. Double-click it.
3. Leave `Create or update the CyANTs conda environment now` checked on the final installer page.
4. Leave `Launch CyANTs GUI` checked to open the GUI immediately after setup.
5. Later, launch `CyANTs GUI` from the desktop icon or Start Menu.

If you only have `CyANTs.exe`, place it in a folder such as `C:\Users\Administrator1\Documents\cyants` and double-click it. The bare `.exe` is a portable GUI launcher, not the full dependency installer. For a new shared PC, distribute `CyANTs_Setup.exe`.

Uninstalling `CyANTs_Setup.exe` removes application files and shortcuts only. It intentionally leaves Miniforge, the `cyants` conda environment, project folders, registered TIFFs, `.ims` sources, and `Reg` output data untouched.
