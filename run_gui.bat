@echo off
REM ============================================================
REM  fieldrrs - double-click this file to start the field GUI.
REM  Needs only a standard python.org install. No pip, no internet.
REM ============================================================
cd /d "%~dp0"

REM Try the Windows launcher first, then a plain python on PATH.
where py >nul 2>nul
if %errorlevel%==0 (
    py -3 -m fieldrrs
    goto done
)

where python >nul 2>nul
if %errorlevel%==0 (
    python -m fieldrrs
    goto done
)

echo.
echo ============================================================
echo  Python was not found on this tablet.
echo.
echo  Install it BEFORE going into the field:
echo    1. https://www.python.org/downloads/windows/
echo    2. Run the installer
echo    3. TICK the box "Add python.exe to PATH"
echo    4. Finish, then double-click this file again
echo.
echo  Nothing else needs installing. No pip, no internet.
echo ============================================================
echo.
pause
exit /b 1

:done
if %errorlevel% neq 0 (
    echo.
    echo The GUI exited with an error - the message is above.
    pause
)
