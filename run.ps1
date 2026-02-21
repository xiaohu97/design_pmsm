param(
    [string]$EnvName = "motor"
)

$ErrorActionPreference = "Stop"
$CondaExe = "D:\miniconda3\Scripts\conda.exe"

if (Test-Path $CondaExe) {
    & $CondaExe run -n $EnvName python .\design_pmsm.py --save-csv --save-png --out .\output
} else {
    python .\design_pmsm.py --save-csv --save-png --out .\output
}
