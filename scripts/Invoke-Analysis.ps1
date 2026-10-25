param(
    [string]$Design,
    [switch]$Live,
    [string]$Connect = '127.0.0.1:7878',
    [double]$MeshSizeM = 2.5,
    [string]$CaseFile
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$cli = Join-Path $projectRoot 'tools\solidcraft-0.20.0\solidcraft-0.20.0-windows-x64-portable\solidcraft-cli.exe'
if ($Live -and $Design) { throw 'Choose -Live or -Design' }
if (-not $CaseFile) { $CaseFile = Join-Path $projectRoot 'examples\elastic.json' }
if (-not $Design -and -not $Live) { $Design = Join-Path $projectRoot 'examples\block.json' }
if ($Design) { $Design = (Resolve-Path -LiteralPath $Design).Path }
$CaseFile = (Resolve-Path -LiteralPath $CaseFile).Path
$runName = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$runDir = Join-Path $projectRoot "runs\$runName"
New-Item -ItemType Directory -Path $runDir -Force | Out-Null
Push-Location -LiteralPath $projectRoot
try {
    Copy-Item -LiteralPath $CaseFile -Destination "$runDir\case.json"
    if ($Live) {
        & $python tools\solidcraft_bridge.py --cli $cli --connect $Connect --output-dir $runDir
        if ($LASTEXITCODE -ne 0) { throw 'SolidCraft live export failed' }
    } else {
        & $cli run --in $Design --out "$runDir\model.step" --save "$runDir\model.solidcraft" --json-out "$runDir\solidcraft-measure.json"
        if ($LASTEXITCODE -ne 0) { throw 'SolidCraft export failed' }
    }
    & $python tools\mesh_step.py "$runDir\model.step" "$runDir\mesh.json" --mesh-size-m $MeshSizeM
    if ($LASTEXITCODE -ne 0) { throw 'Gmsh meshing failed' }
    cargo run --locked --release -- validate "$runDir\mesh.json"
    if ($LASTEXITCODE -ne 0) { throw 'Rust mesh validation failed' }
    cargo run --locked --release -- solve-elastic "$runDir\mesh.json" "$runDir\case.json" "$runDir\solution.json"
    if ($LASTEXITCODE -ne 0) { throw 'Rust elasticity solve failed' }
    & $python tools\export_vtu.py "$runDir\mesh.json" "$runDir\solution.json" "$runDir\result.vtu"
    if ($LASTEXITCODE -ne 0) { throw 'VTU export failed' }
    Write-Host "Completed (linear elastic bootstrap): $runDir"
} finally { Pop-Location }
