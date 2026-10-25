param([switch]$WithReferenceSources)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
Push-Location -LiteralPath $projectRoot
try {
    if (-not (Get-Command rustup -ErrorAction SilentlyContinue)) {
        $downloadDir = Join-Path $projectRoot 'tools\downloads'
        New-Item -ItemType Directory -Path $downloadDir -Force | Out-Null
        Invoke-WebRequest -Uri 'https://win.rustup.rs/x86_64' -OutFile "$downloadDir\rustup-init.exe"
        & "$downloadDir\rustup-init.exe" -y --profile minimal --default-toolchain none
        if ($LASTEXITCODE -ne 0) { throw 'rustup installation failed' }
    }
    $vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
    if (-not (Test-Path -LiteralPath $vswhere)) { throw 'Install Visual Studio Build Tools with Desktop development with C++ first' }
    $msvc = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $msvc) { throw 'MSVC x64 build tools missing' }
    rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
    if ($LASTEXITCODE -ne 0) { throw 'Rust toolchain installation failed' }
    if (-not (Test-Path -LiteralPath .venv\Scripts\python.exe)) {
        python -m venv .venv
        if ($LASTEXITCODE -ne 0) { throw 'Python venv creation failed; Python 3.13 recommended' }
    }
    & .\.venv\Scripts\python.exe -m pip install --disable-pip-version-check -r requirements.lock.txt
    if ($LASTEXITCODE -ne 0) { throw 'Python dependency installation failed' }
    $sourceLock = Get-Content -LiteralPath sources.lock.json -Raw | ConvertFrom-Json
    $packageDir = Join-Path $projectRoot 'tools\solidcraft-0.20.0'
    $cli = Join-Path $packageDir 'solidcraft-0.20.0-windows-x64-portable\solidcraft-cli.exe'
    if (-not (Test-Path -LiteralPath $cli)) {
        New-Item -ItemType Directory -Path tools\downloads -Force | Out-Null
        $archive = Join-Path $projectRoot 'tools\downloads\solidcraft-0.20.0-windows-x64-portable.zip'
        if (-not (Test-Path -LiteralPath $archive)) { Invoke-WebRequest -Uri $sourceLock.solidcraft_binary.url -OutFile $archive }
        $actual = (Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLower()
        if ($actual -ne $sourceLock.solidcraft_binary.sha256) { throw 'SolidCraft SHA256 mismatch' }
        if (Test-Path -LiteralPath $packageDir) { throw 'Incomplete SolidCraft directory exists; inspect it before retrying' }
        Expand-Archive -LiteralPath $archive -DestinationPath $packageDir
    }
    if ($WithReferenceSources) {
        New-Item -ItemType Directory -Path external -Force | Out-Null
        foreach ($repo in $sourceLock.repositories) {
            $name = [IO.Path]::GetFileNameWithoutExtension($repo.url)
            $checkout = Join-Path $projectRoot "external\$name"
            if (Test-Path -LiteralPath $checkout) {
                $head = git -C $checkout rev-parse HEAD
                if ($head -ne $repo.commit) { throw "Existing checkout differs: $checkout" }
            } else {
                git clone --depth 1 --filter=blob:none --sparse --no-checkout $repo.url $checkout
                if ($LASTEXITCODE -ne 0) { throw 'Reference clone failed' }
                git -C $checkout sparse-checkout set --no-cone '/README.md' '/AGENTS.md' '/Cargo.toml' '/docs/*.md' '/src/vba/*' '/native/*.c' '/native/*.h'
                if ($LASTEXITCODE -ne 0) { throw 'Sparse reference selection failed' }
                git -C $checkout fetch --depth 1 origin $repo.commit
                if ($LASTEXITCODE -ne 0) { throw 'Pinned reference fetch failed' }
                git -C $checkout checkout --detach $repo.commit
                if ($LASTEXITCODE -ne 0) { throw 'Reference checkout failed' }
            }
        }
    }
    cargo build --locked --release
    if ($LASTEXITCODE -ne 0) { throw 'Rust build failed' }
    if (-not (Test-Path -LiteralPath tools\bootstrap-block.solidcraft)) {
        & $cli run examples\block.json --save tools\bootstrap-block.solidcraft --out tools\bootstrap-block.step --json-out tools\bootstrap-block-measure.json
        if ($LASTEXITCODE -ne 0) { throw 'SolidCraft bootstrap model failed' }
    }
    if (-not (Test-Path -LiteralPath tools\bootstrap-layers.solidcraft)) {
        & $cli run examples\layers.json --save tools\bootstrap-layers.solidcraft --json-out tools\bootstrap-layers-measure.json
        if ($LASTEXITCODE -ne 0) { throw 'SolidCraft layered model failed' }
    }
    Write-Host 'Environment ready. Run .\scripts\Invoke-Analysis.ps1'
} finally { Pop-Location }
