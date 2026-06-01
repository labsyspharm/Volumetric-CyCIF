@echo off
setlocal

set "BUILD_DIR=%~dp0"
for %%P in ("%BUILD_DIR%..\..") do set "REPO_ROOT=%%~fP"

echo [CyANTs] Building Windows executable first...
call "%BUILD_DIR%build_windows_exe.bat" || exit /b 1

set "ISCC_EXE="
for %%P in (
  "%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe"
  "%ProgramFiles%\Inno Setup 6\ISCC.exe"
  "%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe"
) do (
  if exist "%%~P" (
    set "ISCC_EXE=%%~P"
    goto :found_iscc
  )
)

where ISCC.exe >nul 2>nul
if not errorlevel 1 (
  for /f "delims=" %%P in ('where ISCC.exe 2^>nul') do (
    set "ISCC_EXE=%%P"
    goto :found_iscc
  )
)

echo [CyANTs] Inno Setup compiler was not found.
echo [CyANTs] Install Inno Setup 6, then rerun this command:
echo   https://jrsoftware.org/isdl.php
exit /b 1

:found_iscc
echo [CyANTs] Using Inno Setup compiler: %ISCC_EXE%
"%ISCC_EXE%" "%BUILD_DIR%installer\CyANTs_Setup.iss" || exit /b 1

echo.
echo [CyANTs] Done.
echo [CyANTs] Installer:
echo   %REPO_ROOT%\dist\CyANTs_Setup.exe
echo.
endlocal
