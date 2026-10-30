param(
    [Parameter(Mandatory = $true)][string]$EnvironmentPath,
    [Parameter(Mandatory = $true)][string]$Wheelhouse,
    [string]$Python = 'python'
)
$ErrorActionPreference = 'Stop'
$verificationRoot = Split-Path -Parent $PSScriptRoot
$verificationTarget = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($EnvironmentPath)
if (Test-Path -LiteralPath $verificationTarget) { throw 'Verification environment must be a new path; existing files are preserved' }
if (-not (Test-Path -LiteralPath ([IO.Path]::GetDirectoryName($verificationTarget)) -PathType Container)) { throw 'Environment parent directory must exist' }
$verificationWheelhouse = (Resolve-Path -LiteralPath $Wheelhouse).Path
$verificationLock = Join-Path $verificationRoot 'requirements-verification-win-cp313.lock.txt'
$verificationInventoryPath = Join-Path $verificationRoot 'tools\verification_wheels_windows.json'
$verificationInventory = Get-Content -LiteralPath $verificationInventoryPath -Raw | ConvertFrom-Json
if ($verificationInventory.schema -ne 'soil-fem-verification-wheels/1' -or $verificationInventory.python -ne 'cp313' -or $verificationInventory.platform -ne 'win_amd64') { throw 'Unsupported verification wheel inventory' }
$verificationRows = @($verificationInventory.wheels.PSObject.Properties)
if ($verificationRows.Count -ne 9) { throw 'Expected all nine pinned verification wheels' }
$verificationLockRows = @{}
foreach ($verificationLine in (Get-Content -LiteralPath $verificationLock)) {
    if ($verificationLine.Trim().Length -eq 0 -or $verificationLine.Trim().StartsWith('#')) { continue }
    if ($verificationLine -cnotmatch '^([A-Za-z0-9_.-]+)==([A-Za-z0-9_.+!-]+) --hash=sha256:([0-9a-f]{64})$') { throw 'Unsupported or unpinned verification requirement' }
    $verificationName = $Matches[1].ToLowerInvariant() -replace '[-_.]+','-'
    if ($verificationLockRows.ContainsKey($verificationName)) { throw 'Duplicate verification requirement' }
    $verificationLockRows[$verificationName] = @{version=$Matches[2];sha256=$Matches[3]}
}
if ($verificationLockRows.Count -ne $verificationRows.Count) { throw 'Verification lock and wheel inventory differ' }
$verificationSeen = @{}
foreach ($verificationRow in $verificationRows) {
    if ([IO.Path]::GetFileName($verificationRow.Name) -cne $verificationRow.Name -or -not $verificationRow.Name.EndsWith('.whl')) { throw 'Invalid wheel inventory filename' }
    $verificationName = $verificationRow.Value.name.ToLowerInvariant() -replace '[-_.]+','-'
    if (-not $verificationLockRows.ContainsKey($verificationName)) { throw 'Wheel missing from verification lock' }
    if ($verificationLockRows[$verificationName].version -cne $verificationRow.Value.version -or $verificationLockRows[$verificationName].sha256 -cne $verificationRow.Value.sha256) { throw 'Wheel and locked requirement identity differ' }
    $verificationFile = Join-Path $verificationWheelhouse $verificationRow.Name
    if ((Get-Item -LiteralPath $verificationFile).Length -ne $verificationRow.Value.bytes) { throw "Verification wheel size differs: $($verificationRow.Name)" }
    if ((Get-FileHash -LiteralPath $verificationFile -Algorithm SHA256).Hash.ToLowerInvariant() -cne $verificationRow.Value.sha256) { throw "Verification wheel hash differs: $($verificationRow.Name)" }
    $verificationSeen[$verificationRow.Name] = $verificationRow.Value.sha256
}
foreach ($verificationFile in (Get-ChildItem -LiteralPath $verificationWheelhouse -File -Filter '*.whl')) {
    if (-not $verificationSeen.ContainsKey($verificationFile.Name)) { throw 'Wheelhouse contains an unlisted wheel' }
}
$verificationPythonCommand = Get-Command $Python -CommandType Application
$verificationBasePython = $verificationPythonCommand.Source
$verificationBaseJson = & $verificationBasePython -B -I -X utf8 -c 'import json,platform,struct,sys;print(json.dumps(dict(version=platform.python_version(),platform=sys.platform,bits=struct.calcsize(''P'')*8,implementation=sys.implementation.name)))'
if ($LASTEXITCODE -ne 0) { throw 'Verification Python identity query failed' }
$verificationBase = $verificationBaseJson | ConvertFrom-Json
if ($verificationBase.version -ne '3.13.14' -or $verificationBase.platform -ne 'win32' -or $verificationBase.bits -ne 64 -or $verificationBase.implementation -ne 'cpython') { throw 'This pinned verification runtime requires CPython 3.13.14 Windows x64' }
$verificationSourceHashes = @{}
foreach ($verificationSource in @($PSCommandPath,$verificationLock,$verificationInventoryPath)) { $verificationSourceHashes[$verificationSource] = (Get-FileHash -LiteralPath $verificationSource -Algorithm SHA256).Hash.ToLowerInvariant() }
& $verificationBasePython -B -I -X utf8 -m venv $verificationTarget
if ($LASTEXITCODE -ne 0) { throw 'Verification venv creation failed; preserve the directory for inspection' }
$verificationPython = Join-Path $verificationTarget 'Scripts\python.exe'
$verificationPipCode = @'
import os,runpy,sys
for name in list(os.environ):
    if name.upper().startswith('PIP_'):del os.environ[name]
# pip's configuration loader treats this as disabling every configuration file.
os.environ['PIP_CONFIG_FILE']=os.devnull
sys.argv=['pip',*sys.argv[1:]]
runpy.run_module('pip',run_name='__main__')
'@
$verificationPipReport = Join-Path $verificationTarget 'pip-install-report.json'
& $verificationPython -B -I -X utf8 -c $verificationPipCode install --disable-pip-version-check --no-cache-dir --no-index --find-links $verificationWheelhouse --require-hashes --report $verificationPipReport -r $verificationLock *> (Join-Path $verificationTarget 'pip-install.log')
if ($LASTEXITCODE -ne 0) { throw 'Pinned offline verification dependency installation failed; inspect pip-install.log' }
& $verificationPython -B -I -X utf8 -c $verificationPipCode check *> (Join-Path $verificationTarget 'pip-check.log')
if ($LASTEXITCODE -ne 0) { throw 'Verification dependency consistency failed' }
$verificationCheckCode = @'
from pathlib import Path
import importlib,importlib.metadata,json,platform,re,site,sys
target,inventory_path,report_path=map(Path,sys.argv[1:])
assert Path(sys.prefix).resolve()==target.resolve() and sys.flags.isolated and not site.ENABLE_USER_SITE
assert platform.python_version()=='3.13.14' and sys.platform=='win32'
inventory=json.loads(inventory_path.read_text('utf-8-sig'))
normalize=lambda name:re.sub('[-_.]+','-',name).lower()
expected={normalize(row['name']):row for row in inventory['wheels'].values()}
report=json.loads(report_path.read_text('utf-8'))
assert len(report['install'])==len(expected)==9
installed=set()
for row in report['install']:
    name=normalize(row['metadata']['name']);assert name not in installed;installed.add(name)
    assert row['metadata']['version']==expected[name]['version']
    assert row['download_info']['url'].startswith('file:')
    assert row['download_info']['archive_info']['hashes']['sha256']==expected[name]['sha256']
assert installed==set(expected)
modules={'gmsh':'gmsh','ijson':'ijson','markdown-it-py':'markdown_it','mdurl':'mdurl','meshio':'meshio','numpy':'numpy','pygments':'pygments','rich':'rich','scipy':'scipy'}
packages={}
for name,row in expected.items():
    module=importlib.import_module(modules[name]);path=Path(module.__file__).resolve()
    assert path.is_relative_to(target.resolve())
    version=importlib.metadata.version(row['name']);assert version==row['version']
    packages[name]=dict(version=version,path=str(path))
print(json.dumps(dict(passed=True,python=platform.python_version(),prefix=sys.prefix,
    isolated=True,user_site_disabled=True,all9_imports_from_new_env=True,
    all9_install_sources_local_and_hash_pinned=True,packages=packages,same_host_only=True)))
'@
$verificationRuntimeJson = & $verificationPython -B -I -X utf8 -c $verificationCheckCode $verificationTarget $verificationInventoryPath $verificationPipReport
if ($LASTEXITCODE -ne 0) { throw 'Verification environment isolation or installed artifact identity failed' }
$verificationRuntime = $verificationRuntimeJson | ConvertFrom-Json
foreach ($verificationSource in $verificationSourceHashes.Keys) {
    if ((Get-FileHash -LiteralPath $verificationSource -Algorithm SHA256).Hash.ToLowerInvariant() -cne $verificationSourceHashes[$verificationSource]) { throw 'Verification setup source changed during installation' }
}
foreach ($verificationRow in $verificationRows) {
    if ((Get-FileHash -LiteralPath (Join-Path $verificationWheelhouse $verificationRow.Name) -Algorithm SHA256).Hash.ToLowerInvariant() -cne $verificationRow.Value.sha256) { throw 'Verification wheel changed during installation' }
}
$verificationRecord = [ordered]@{passed=$true;utc=[DateTime]::UtcNow.ToString('o');powershell_version=$PSVersionTable.PSVersion.ToString();runtime=$verificationRuntime;source_sha256=$verificationSourceHashes;wheel_sha256=$verificationSeen;pip_report_sha256=(Get-FileHash -LiteralPath $verificationPipReport -Algorithm SHA256).Hash.ToLowerInvariant();offline=$true;pip_environment_and_configuration_disabled_in_children=$true;existing_environment_packages_not_modified=$true;new_PC_or_other_GPU_not_claimed=$true}
[IO.File]::WriteAllText((Join-Path $verificationTarget 'verification-environment.json'),($verificationRecord | ConvertTo-Json -Depth 10),[Text.UTF8Encoding]::new($false))
Write-Output ($verificationRecord | ConvertTo-Json -Depth 10)
