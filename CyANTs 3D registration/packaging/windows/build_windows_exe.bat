@echo off
setlocal

set "BUILD_DIR=%~dp0"
for %%P in ("%BUILD_DIR%..\..") do set "REPO_ROOT=%%~fP"

echo [CyANTs] Building Windows GUI executable from %REPO_ROOT%
echo [CyANTs] This uses the currently activated conda environment.

pushd "%REPO_ROOT%" || exit /b 1
python -c "import sys; print('[CyANTs] Python:', sys.executable)" || exit /b 1
python -m pip install --upgrade pyinstaller tkinterdnd2 || exit /b 1
python -m PyInstaller --clean --noconfirm --distpath "%REPO_ROOT%\dist" --workpath "%REPO_ROOT%\build\pyinstaller" "%BUILD_DIR%cyants_gui.spec" || exit /b 1
popd

echo.
echo [CyANTs] Done.
echo [CyANTs] Executable:
echo   %REPO_ROOT%\dist\CyANTs.exe
echo.
endlocal
