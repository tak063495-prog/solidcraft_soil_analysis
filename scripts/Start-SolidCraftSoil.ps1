param(
    [string]$Design,
    [int]$Port = 7894,
    [string]$SolidCraftPath,
    [string]$SolverPath
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
if (-not $SolidCraftPath) {
    $SolidCraftPath = Join-Path $projectRoot 'external\solidcraft-soil-fem\target\debug\solidcraft.exe'
}
$SolidCraftPath = (Resolve-Path -LiteralPath $SolidCraftPath).Path
if (-not $SolverPath) { $SolverPath = Join-Path $projectRoot 'target\release\soil-fem-3d.exe' }
$SolverPath = (Resolve-Path -LiteralPath $SolverPath).Path
if (-not (Test-Path -LiteralPath $SolverPath -PathType Leaf)) { throw 'SolverPath must name an installed solver file' }
foreach ($relative in @('.venv\Scripts\python.exe', 'tools\native_job.py', 'target\release\soil-fem-3d.exe')) {
    if (-not (Test-Path -LiteralPath (Join-Path $projectRoot $relative) -PathType Leaf)) {
        throw "Missing FEM environment file: $relative"
    }
}
if ($Port -lt 1024 -or $Port -gt 65535) { throw 'Choose a control port between 1024 and 65535' }
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port is already in use" }
$arguments = @('--control', [string]$Port)
if ($Design) {
    $designPath = (Resolve-Path -LiteralPath $Design).Path
    if ($designPath.Contains('"')) { throw 'A design path cannot contain a quote' }
    $arguments = @(('"' + $designPath + '"')) + $arguments
}
# Interactive CAD requested by the user. Only this child's environment is changed.
$cadProcess = Start-Process -FilePath $SolidCraftPath -ArgumentList $arguments -WorkingDirectory $projectRoot `
    -Environment @{ SOLIDCRAFT_SOIL_FEM_ROOT = $projectRoot; SOLIDCRAFT_SOIL_FEM_SOLVER = $SolverPath } -PassThru
Write-Host "SolidCraft soil host PID=$($cadProcess.Id), control=127.0.0.1:$Port"
Write-Host 'Open Soil Analysis > Analysis to start a job in a new output folder.'
