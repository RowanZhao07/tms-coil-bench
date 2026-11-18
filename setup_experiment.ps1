param([switch]$DownloadFEMTools)
$ErrorActionPreference = 'Stop'
$projectDir = $PSScriptRoot
$envDir = Join-Path $projectDir '.venv'
if (-not (Test-Path -LiteralPath (Join-Path $envDir 'Scripts\python.exe'))) {
    & python -m venv $envDir
    if ($LASTEXITCODE -ne 0) { throw 'Could not create the virtual environment. Use Python 3.13 for the tested dependency set.' }
}
$projectPython = Join-Path $envDir 'Scripts\python.exe'
& $projectPython -m pip install -r (Join-Path $projectDir 'requirements.txt')
if ($LASTEXITCODE -ne 0) { throw 'Python dependency installation failed.' }
if ($DownloadFEMTools) {
    $toolDir = Join-Path $projectDir 'tools'
    New-Item -ItemType Directory -Force -Path $toolDir | Out-Null
    $archives = @(
        @{ Name='gmsh-4.15.2-Windows64.zip'; Url='https://gmsh.info/bin/Windows/gmsh-4.15.2-Windows64.zip' },
        @{ Name='getdp-3.5.0-Windows64c.zip'; Url='https://getdp.info/bin/Windows/getdp-3.5.0-Windows64c.zip' }
    )
    foreach ($archive in $archives) {
        $destination = Join-Path $toolDir $archive.Name
        if (-not (Test-Path -LiteralPath $destination)) {
            Invoke-WebRequest -Uri $archive.Url -OutFile $destination
        }
        Expand-Archive -LiteralPath $destination -DestinationPath $toolDir -Force
    }
}
Write-Host 'Setup complete. Run .\run_all.ps1'
