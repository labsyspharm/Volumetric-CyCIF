@echo off
setlocal

set "APP_DIR=%~dp0"
set "ENV_NAME=cyants"

call "%APP_DIR%cyants_conda.bat" "%ENV_NAME%" || exit /b 1

if exist "%APP_DIR%cyants_gui.py" (
  start "" python "%APP_DIR%cyants_gui.py"
  exit /b 0
)

if exist "%APP_DIR%CyANTs.exe" (
  start "" "%APP_DIR%CyANTs.exe"
  exit /b 0
)

echo [CyANTs] Could not find cyants_gui.py or CyANTs.exe in %APP_DIR%
pause
exit /b 1

endlocal
