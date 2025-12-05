@echo off
REM BioHabs Dashboard Launcher for Windows

echo ========================================
echo   BioHabs Interactive Dashboard
echo ========================================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found. Please install Python 3.8+
    pause
    exit /b 1
)

echo [1/3] Checking dependencies...
pip show flask >nul 2>&1
if errorlevel 1 (
    echo Installing Python dependencies...
    pip install -r requirements.txt
)

echo [2/3] Starting Flask backend...
start "BioHabs Backend" python backend\app.py

timeout /t 3 >nul

echo [3/3] Opening dashboard in browser...
start http://localhost:5000

echo.
echo ========================================
echo   Dashboard is running!
echo   Backend: http://localhost:5000
echo   Press Ctrl+C in backend window to stop
echo ========================================
echo.

pause
