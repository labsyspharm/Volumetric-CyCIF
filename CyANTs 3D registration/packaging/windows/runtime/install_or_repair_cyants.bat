@echo off
setlocal EnableExtensions

set "APP_DIR=%~1"
if "%APP_DIR%"=="" set "APP_DIR=%~dp0"
for %%P in ("%APP_DIR%") do set "APP_DIR=%%~fP"

set "ENV_NAME=%~2"
if "%ENV_NAME%"=="" set "ENV_NAME=cyants"

set "NO_PAUSE=%~3"

echo [CyANTs] Install/repair starting
echo [CyANTs] App directory: %APP_DIR%
echo [CyANTs] Conda environment: %ENV_NAME%

set "CONDA_EXE="
for %%P in (
  "%USERPROFILE%\miniforge3\Scripts\conda.exe"
  "%LOCALAPPDATA%\miniforge3\Scripts\conda.exe"
  "%USERPROFILE%\mambaforge\Scripts\conda.exe"
  "%LOCALAPPDATA%\mambaforge\Scripts\conda.exe"
  "%USERPROFILE%\anaconda3\Scripts\conda.exe"
  "%LOCALAPPDATA%\anaconda3\Scripts\conda.exe"
  "C:\ProgramData\anaconda3\Scripts\conda.exe"
) do (
  if exist "%%~P" (
    set "CONDA_EXE=%%~P"
    goto :have_conda
  )
)

where conda >nul 2>nul
if not errorlevel 1 (
  for /f "delims=" %%P in ('where conda 2^>nul') do (
    set "CONDA_EXE=%%P"
    goto :have_conda
  )
)

echo [CyANTs] No conda installation found. Installing Miniforge for this user...
set "MINIFORGE_DIR=%USERPROFILE%\miniforge3"
set "MINIFORGE_EXE=%TEMP%\Miniforge3-Windows-x86_64.exe"
powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "$ProgressPreference='SilentlyContinue'; Invoke-WebRequest -Uri 'https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe' -OutFile $env:MINIFORGE_EXE"
if errorlevel 1 (
  echo [CyANTs] Failed to download Miniforge.
  echo [CyANTs] Install Miniforge manually, then rerun this repair shortcut.
  exit /b 1
)

start /wait "" "%MINIFORGE_EXE%" /S /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /D=%MINIFORGE_DIR%
if errorlevel 1 (
  echo [CyANTs] Miniforge installer failed.
  exit /b 1
)
set "CONDA_EXE=%MINIFORGE_DIR%\Scripts\conda.exe"

:have_conda
echo [CyANTs] Using conda: %CONDA_EXE%
for %%P in ("%CONDA_EXE%") do set "CONDA_BAT=%%~dpP..\condabin\conda.bat"
if not exist "%CONDA_BAT%" set "CONDA_BAT=%CONDA_EXE%"

call "%CONDA_BAT%" env list | findstr /R /C:"^%ENV_NAME%[ ]" >nul
if errorlevel 1 (
  echo [CyANTs] Creating environment "%ENV_NAME%"...
  call "%CONDA_BAT%" env create -n "%ENV_NAME%" -f "%APP_DIR%\environment.yml"
) else (
  echo [CyANTs] Updating environment "%ENV_NAME%"...
  call "%CONDA_BAT%" env update -n "%ENV_NAME%" -f "%APP_DIR%\environment.yml" --prune
)
if errorlevel 1 (
  echo [CyANTs] Conda environment setup failed.
  exit /b 1
)

call "%CONDA_BAT%" activate "%ENV_NAME%"
if errorlevel 1 (
  echo [CyANTs] Environment activation failed.
  exit /b 1
)

python -m pip install -r "%APP_DIR%\requirements.txt"
if errorlevel 1 (
  echo [CyANTs] Pip dependency installation failed.
  exit /b 1
)

python -c "import ants, h5py, SimpleITK, tifffile; print('[CyANTs] Python environment ready')"
if errorlevel 1 (
  echo [CyANTs] Dependency verification failed.
  exit /b 1
)

echo [CyANTs] Install/repair complete.
if /I not "%NO_PAUSE%"=="--no-pause" pause
exit /b 0
