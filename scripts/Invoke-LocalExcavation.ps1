param(
    [string]$Design,
    [string]$GroundFile,
    [string]$StageFile,
    [switch]$Live,
    [string]$Connect = '127.0.0.1:7878',
    [double]$MeshSizeM = 2.0,
    [ValidateSet(1,2)][int]$ElementOrder = 2,
    [ValidateSet(3,4,5)][int]$QuadratureOrder = 3,
    [string]$SolverPath,
    [string]$OutputDirectory
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
if (-not $Design -and -not $Live) { $Design = Join-Path $projectRoot 'examples\local-excavation.json' }
if (-not $GroundFile) { $GroundFile = Join-Path $projectRoot 'examples\ground-local.json' }
if (-not $StageFile) { $StageFile = Join-Path $projectRoot 'examples\stages-local-sparse.json' }
$arguments = @{ GroundFile = $GroundFile; StageFile = $StageFile; MeshSizeM = $MeshSizeM; Connect = $Connect; ElementOrder=$ElementOrder; QuadratureOrder=$QuadratureOrder }
if ($Design) { $arguments.Design = $Design }
if ($Live) { $arguments.Live = $true }
if ($OutputDirectory) { $arguments.OutputDirectory = $OutputDirectory }
if ($SolverPath) { $arguments.SolverPath = $SolverPath }
& "$PSScriptRoot\Invoke-Stages.ps1" @arguments
