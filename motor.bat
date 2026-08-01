@echo off
setlocal

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0motor.ps1" %*
set "MOTOR_EXIT_CODE=%ERRORLEVEL%"

endlocal & exit /b %MOTOR_EXIT_CODE%
