$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $RepoRoot
$Python = if ($env:PYTHON) { $env:PYTHON } else { "python" }
# La configuración de pytest.ini (raíz) ya limita las pruebas a las de Track 1 y los esquemas compartidos.
& $Python -m pytest -v
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
