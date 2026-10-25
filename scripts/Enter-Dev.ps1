$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$projectRoot\.venv\Scripts;$projectRoot\tools\solidcraft-0.20.0\solidcraft-0.20.0-windows-x64-portable;$env:PATH"
Set-Location -LiteralPath $projectRoot
Write-Host "Soil FEM 3D: $projectRoot"
Write-Host 'cargo test --locked | .\scripts\Invoke-Analysis.ps1 | .\scripts\Start-SolidCraft.ps1'
