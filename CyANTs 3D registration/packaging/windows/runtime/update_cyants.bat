@echo off
setlocal

set "APP_DIR=%~dp0"

echo [CyANTs] Updating installed files from Git, if available...
cd /d "%APP_DIR%" || exit /b 1
where git >nul 2>nul
if errorlevel 1 (
  echo [CyANTs] Git was not found on PATH. Skipping git pull.
) else (
  git pull
)

call "%APP_DIR%install_or_repair_cyants.bat" "%APP_DIR%" cyants

endlocal
