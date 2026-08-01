param(
    [string]$EnvName = "motor",
    [string]$PythonExe = ""
)

$Entry = Join-Path $PSScriptRoot "motor.ps1"
& $Entry quick -CondaEnv $EnvName -PythonExe $PythonExe
exit $LASTEXITCODE
