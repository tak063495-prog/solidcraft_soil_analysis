param([string]$Design, [int]$Port = 7878)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$app = Join-Path $projectRoot 'tools\solidcraft-0.20.0\solidcraft-0.20.0-windows-x64-portable\solidcraft.exe'
if (-not $Design) { $Design = Join-Path $projectRoot 'tools\bootstrap-block.solidcraft' }
$Design = (Resolve-Path -LiteralPath $Design).Path
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port is already in use" }
# Interactive CAD window requested for modelling; bind the control server to loopback.
$cadProcess = Start-Process -FilePath $app -ArgumentList @(('"' + $Design + '"'), '--control', $Port) -PassThru
Write-Host "SolidCraft PID=$($cadProcess.Id), control=127.0.0.1:$Port"
Write-Host "Single body elasticity: .\scripts\Invoke-Analysis.ps1 -Live -Connect 127.0.0.1:$Port"
Write-Host "Layered K0/MCC state: .\scripts\Invoke-Ground.ps1 -Live -Connect 127.0.0.1:$Port"
