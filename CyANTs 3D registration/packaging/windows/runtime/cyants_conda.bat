@echo off
set "ENV_NAME=%~1"
if "%ENV_NAME%"=="" set "ENV_NAME=cyants"

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
    goto :found_conda
  )
)

where conda >nul 2>nul
if not errorlevel 1 (
  for /f "delims=" %%P in ('where conda 2^>nul') do (
    set "CONDA_EXE=%%P"
    goto :found_conda
  )
)

echo [CyANTs] Could not find conda.
echo [CyANTs] Run Install/Repair CyANTs from the Start Menu first.
exit /b 1

:found_conda
for %%P in ("%CONDA_EXE%") do set "CONDA_BAT=%%~dpP..\condabin\conda.bat"
if not exist "%CONDA_BAT%" set "CONDA_BAT=%CONDA_EXE%"

call "%CONDA_BAT%" activate "%ENV_NAME%"
if errorlevel 1 (
  echo [CyANTs] Could not activate environment "%ENV_NAME%".
  echo [CyANTs] Run Install/Repair CyANTs from the Start Menu first.
  exit /b 1
)
exit /b 0
