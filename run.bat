@echo off
setlocal

set ENV_NAME=motor
set CONDA_EXE=D:\miniconda3\Scripts\conda.exe

if exist "%CONDA_EXE%" (
    "%CONDA_EXE%" run -n %ENV_NAME% python design_pmsm.py --save-csv --save-png --out output
) else (
    python design_pmsm.py --save-csv --save-png --out output
)

endlocal
