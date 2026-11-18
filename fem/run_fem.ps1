param([switch]$Quick)
$ErrorActionPreference = 'Stop'
$TaskRoot = Split-Path -Parent $PSScriptRoot
$TaskPython = Join-Path $TaskRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $TaskPython)) {
    $TaskPython = (Get-Command python -ErrorAction Stop).Source
}
if ($Quick) { & $TaskPython (Join-Path $PSScriptRoot 'run_fem.py') --quick }
else { & $TaskPython (Join-Path $PSScriptRoot 'run_fem.py') }
if ($LASTEXITCODE -ne 0) { throw "FEM run failed: $LASTEXITCODE" }
