@echo off
setlocal
cd /d "%~dp0"
set PORT=8123

where python >nul 2>nul
if %errorlevel%==0 (set PY=python) else (set PY=py -3)

echo.
echo   Warehouse Area Measure
echo   ----------------------------------------
echo   Opening http://localhost:%PORT%
echo.
echo   A small server window will appear, minimised.
echo   Close that window when you are done measuring.
echo.

start "Warehouse Area Measure - server (close to stop)" /min %PY% -m http.server %PORT%
timeout /t 2 /nobreak >nul
start "" http://localhost:%PORT%

timeout /t 3 /nobreak >nul
endlocal
