[CmdletBinding()]
param(
    [Parameter(Position = 0)]
    [ValidateSet("quick", "femm-validate", "test")]
    [string]$Mode = "quick",

    [string]$PythonExe,
    [string]$CondaExe,
    [string]$CondaEnv = "motor",
    [string]$Config,
    [string]$Out,
    [Nullable[double]]$Torque,
    [Nullable[double]]$MaxRpm,
    [Nullable[int]]$Points,
    [Nullable[double]]$Vdc,
    [Nullable[double]]$CurrentLimit,
    [switch]$NoPng,

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ToolArgs
)

$ErrorActionPreference = "Stop"
$repoRoot = $PSScriptRoot
$invariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

if ([string]::IsNullOrWhiteSpace($Config)) {
    $Config = Join-Path $repoRoot "motor_config.json"
}

function Resolve-Application {
    param([Parameter(Mandatory = $true)][string]$NameOrPath)

    if (Test-Path -LiteralPath $NameOrPath -PathType Leaf) {
        return (Resolve-Path -LiteralPath $NameOrPath).Path
    }

    $command = Get-Command -Name $NameOrPath -CommandType Application -ErrorAction SilentlyContinue |
        Select-Object -First 1
    if ($null -ne $command) {
        return $command.Source
    }
    return $null
}

function Test-DirectPython {
    param([Parameter(Mandatory = $true)][string]$Executable)

    & $Executable -c "import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)" *> $null
    return $LASTEXITCODE -eq 0
}

function Test-CondaPython {
    param(
        [Parameter(Mandatory = $true)][string]$Executable,
        [Parameter(Mandatory = $true)][string]$EnvironmentName
    )

    & $Executable run -n $EnvironmentName python -c "import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)" *> $null
    return $LASTEXITCODE -eq 0
}

$runnerKind = $null
$runnerPath = $null

if (-not [string]::IsNullOrWhiteSpace($PythonExe)) {
    $candidate = Resolve-Application $PythonExe
    if ($null -eq $candidate -or -not (Test-DirectPython $candidate)) {
        throw "-PythonExe does not point to a working Python 3.10+ interpreter: $PythonExe"
    }
    $runnerKind = "python"
    $runnerPath = $candidate
}

if ($null -eq $runnerPath) {
    $venvPython = Join-Path $repoRoot ".venv\Scripts\python.exe"
    if ((Test-Path -LiteralPath $venvPython -PathType Leaf) -and (Test-DirectPython $venvPython)) {
        $runnerKind = "python"
        $runnerPath = $venvPython
    }
}

if ($null -eq $runnerPath) {
    $condaCandidates = @()
    if (-not [string]::IsNullOrWhiteSpace($CondaExe)) {
        $condaCandidates += $CondaExe
    }
    if (-not [string]::IsNullOrWhiteSpace($env:MOTOR_CONDA_EXE)) {
        $condaCandidates += $env:MOTOR_CONDA_EXE
    }
    $pathConda = Get-Command -Name "conda" -CommandType Application -ErrorAction SilentlyContinue |
        Select-Object -First 1
    if ($null -ne $pathConda) {
        $condaCandidates += $pathConda.Source
    }
    $condaCandidates += "D:\miniconda3\Scripts\conda.exe"

    foreach ($condaCandidate in ($condaCandidates | Select-Object -Unique)) {
        $candidate = Resolve-Application $condaCandidate
        if ($null -ne $candidate -and (Test-CondaPython $candidate $CondaEnv)) {
            $runnerKind = "conda"
            $runnerPath = $candidate
            break
        }
    }
}

if ($null -eq $runnerPath) {
    $candidate = Resolve-Application "python"
    if ($null -ne $candidate -and (Test-DirectPython $candidate)) {
        $runnerKind = "python"
        $runnerPath = $candidate
    }
}

if ($null -eq $runnerPath) {
    throw (
        "No working Python 3.10+ interpreter was found. Use -PythonExe, create " +
        "'.venv', or set MOTOR_CONDA_EXE for a Conda installation containing the '$CondaEnv' environment."
    )
}

if ($Mode -ne "test" -and -not (Test-Path -LiteralPath $Config -PathType Leaf)) {
    throw "Motor configuration does not exist: $Config"
}

switch ($Mode) {
    "quick" {
        $quickOut = if ([string]::IsNullOrWhiteSpace($Out)) {
            Join-Path $repoRoot "output_quick"
        } else {
            $Out
        }
        $pythonArgs = @(
            (Join-Path $repoRoot "design_pmsm.py"),
            "--config", $Config,
            "--out", $quickOut
        )
        if ($null -ne $Torque) {
            $pythonArgs += @("--torque", $Torque.ToString($invariantCulture))
        }
        if ($null -ne $MaxRpm) {
            $pythonArgs += @("--max-rpm", $MaxRpm.ToString($invariantCulture))
        }
        if ($null -ne $Points) {
            $pythonArgs += @("--points", $Points.ToString($invariantCulture))
        }
        if ($null -ne $Vdc) {
            $pythonArgs += @("--vdc", $Vdc.ToString($invariantCulture))
        }
        if ($null -ne $CurrentLimit) {
            $pythonArgs += @(
                "--current-limit",
                $CurrentLimit.ToString($invariantCulture)
            )
        }
        if ($NoPng) { $pythonArgs += "--no-png" }
        $pythonArgs += @($ToolArgs)
    }
    "femm-validate" {
        $femmOut = if ([string]::IsNullOrWhiteSpace($Out)) {
            Join-Path $repoRoot "output_femm_validation"
        } else {
            $Out
        }
        $pythonArgs = @(
            (Join-Path $repoRoot "femm_spm_template.py"),
            "--config", $Config,
            "--mesh-level", "coarse",
            "--validation-iq", "1",
            "--validation-delta-current", "1",
            "--validation-steps", "3",
            "--out", $femmOut
        ) + @($ToolArgs) + @("--analysis", "validate")
    }
    "test" {
        $pythonArgs = @(
            "-m", "unittest", "discover",
            "-s", (Join-Path $repoRoot "tests"),
            "-v"
        ) + @($ToolArgs)
    }
}

if ($runnerKind -eq "conda") {
    Write-Host "Using Conda environment '$CondaEnv': $runnerPath"
    & $runnerPath run --no-capture-output -n $CondaEnv python @pythonArgs
} else {
    Write-Host "Using Python: $runnerPath"
    & $runnerPath @pythonArgs
}

exit $LASTEXITCODE
