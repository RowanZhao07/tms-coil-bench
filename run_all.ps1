param([switch]$UseGPU, [switch]$SkipFEM)
$ErrorActionPreference = 'Stop'
$projectDir = $PSScriptRoot
$matlabSource = (Join-Path $projectDir 'matlab').Replace('\','/').Replace("'","''")
$gpuLiteral = if ($UseGPU) { 'true' } else { 'false' }
$matlabCommand = "addpath('$matlabSource'); run_study($gpuLiteral); benchmark_device; run_sensitivity; run_spatial_convergence($gpuLiteral);"
& matlab -batch $matlabCommand -logfile (Join-Path $projectDir 'results\matlab_batch_latest.log')
if ($LASTEXITCODE -ne 0) { throw 'MATLAB experiment failed. Read results/matlab_batch_latest.log.' }
$projectPython = Join-Path $projectDir '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $projectPython)) { throw 'Run setup_experiment.ps1 first.' }
& $projectPython (Join-Path $projectDir 'simulation\compare_magpylib.py')
if ($LASTEXITCODE -ne 0) { throw 'Magpylib comparison failed.' }
& $projectPython (Join-Path $projectDir 'simulation\export_gmsh_views.py')
if ($LASTEXITCODE -ne 0) { throw 'Gmsh visualization export failed.' }
if (-not $SkipFEM) {
    & (Join-Path $projectDir 'fem\run_fem.ps1')
    if ($LASTEXITCODE -ne 0) { throw 'FEM experiment failed.' }
    & $projectPython (Join-Path $projectDir 'fem\run_multiturn_case.py')
    if ($LASTEXITCODE -ne 0) { throw 'Multi-turn FEM experiment failed.' }
}
if (Test-Path -LiteralPath (Join-Path $projectDir 'fem\results\axis_comparison.csv')) {
    & matlab -batch "addpath('$matlabSource'); make_validation_figure;"
    if ($LASTEXITCODE -ne 0) { throw 'Validation figure export failed.' }
}
Write-Host 'Experiments complete. Open README.md and results/.'
