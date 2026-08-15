@echo off
REM ============================================================================
REM  Build fieldrrs.exe  -  a standalone Windows executable.
REM
REM  Run this ONCE on any Windows machine that has Python and an internet
REM  connection. The resulting dist\fieldrrs.exe is self-contained: copy it to
REM  the field tablet and it runs with NO Python installed there at all.
REM
REM  This has to run on Windows. PyInstaller cannot cross-compile, so an exe
REM  cannot be produced from Linux or macOS.
REM ============================================================================
setlocal
cd /d "%~dp0"

echo.
echo === 1/4  locating Python ===
set PY=
where py >nul 2>nul && set PY=py -3
if "%PY%"=="" ( where python >nul 2>nul && set PY=python )
if "%PY%"=="" (
    echo.
    echo  Python was not found.
    echo  Install it from https://www.python.org/downloads/windows/
    echo  and TICK "Add python.exe to PATH", then run this again.
    echo.
    pause
    exit /b 1
)
%PY% --version

echo.
echo === 2/4  self-test before building ===
%PY% tests\test_fieldrrs.py
if %errorlevel% neq 0 (
    echo.
    echo  TESTS FAILED. Not building an exe from a broken tree.
    pause
    exit /b 1
)

echo.
echo === 3/4  installing PyInstaller (needs internet, one time) ===
%PY% -m pip install --upgrade pyinstaller
if %errorlevel% neq 0 (
    echo.
    echo  Could not install PyInstaller. If this machine has no internet, run
    echo  this script on one that does and copy dist\fieldrrs.exe across.
    pause
    exit /b 1
)

echo.
echo === 4/4  building ===
%PY% -m PyInstaller --clean --noconfirm fieldrrs.spec
if %errorlevel% neq 0 (
    echo.
    echo  Build failed. The output above says why.
    pause
    exit /b 1
)

echo.
echo ============================================================================
echo  DONE.   dist\fieldrrs.exe
echo.
echo  Copy that single file anywhere. It needs no Python.
echo  Take FIELD_CARD.pdf with you as well.
echo.
echo  TEST IT NOW, before the field:
echo    1. %PY% make_demo_data.py
echo    2. run dist\fieldrrs.exe
echo    3. LOAD WATER  -^> demo_scans\station1_water.sed
echo       LOAD SKY    -^> demo_scans\station1_sky.sed
echo    4. COMPUTE  -^>  expect Rrs(443) = 0.00294 and no warnings
echo ============================================================================
echo.
pause
