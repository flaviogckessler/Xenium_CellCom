$ErrorActionPreference = 'Stop'
$bookDirectory = $PSScriptRoot
$temporaryDirectory = Join-Path $bookDirectory '.render-tmp'
New-Item -ItemType Directory -Path $temporaryDirectory -Force | Out-Null

$quartoCommand = Get-Command quarto -ErrorAction SilentlyContinue
$quartoCandidates = @(
    $(if ($quartoCommand) { $quartoCommand.Source }),
    "$env:ProgramFiles\Quarto\bin\quarto.exe",
    "$env:ProgramFiles\RStudio\resources\app\bin\quarto\bin\quarto.exe",
    "$env:LOCALAPPDATA\Programs\Quarto\bin\quarto.exe"
)
$quarto = $quartoCandidates | Where-Object { $_ -and (Test-Path -LiteralPath $_) } | Select-Object -First 1
if (-not $quarto) { throw 'Install Quarto, or use an RStudio installation that includes Quarto.' }

$environmentNames = @('TEMP', 'TMP', 'TMPDIR', 'QUARTO_R', 'R_LIBS_USER', 'LC_ALL', 'LC_CTYPE')
$savedEnvironment = @{}
foreach ($name in $environmentNames) {
    $savedEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
}

try {
    # Pandoc needs a temporary path without URL fragment characters such as #.
    $env:TEMP = $temporaryDirectory
    $env:TMP = $temporaryDirectory
    $env:TMPDIR = $temporaryDirectory
    $rInstallations = Get-ChildItem -LiteralPath "$env:ProgramFiles\R" -Directory -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -match '^R-\d+\.\d+\.\d+$' } |
        Sort-Object { [version]($_.Name -replace '^R-', '') } -Descending
    if (-not $env:QUARTO_R -and $rInstallations) {
        $selectedR = $rInstallations | Select-Object -First 1
        $env:QUARTO_R = Join-Path $selectedR.FullName 'bin'
        $versionParts = [version]($selectedR.Name -replace '^R-', '')
        $userLibrary = Join-Path $env:USERPROFILE (
            'AppData\Local\R\win-library\{0}.{1}' -f $versionParts.Major, $versionParts.Minor
        )
        if (-not $env:R_LIBS_USER -and (Test-Path -LiteralPath $userLibrary)) {
            $env:R_LIBS_USER = $userLibrary
        }
    }
    if ($env:LC_ALL -eq 'C.UTF-8') { $env:LC_ALL = 'English_United States.utf8' }
    if ($env:LC_CTYPE -eq 'C.UTF-8') { $env:LC_CTYPE = 'English_United States.utf8' }

    Push-Location -LiteralPath $bookDirectory
    try {
        & $quarto render 'Xenium Cell Com book.qmd' --to html 2>&1 |
            Tee-Object -FilePath (Join-Path $bookDirectory 'render.log')
        if ($LASTEXITCODE -ne 0) { throw "Quarto render failed with exit code $LASTEXITCODE." }
    } finally {
        Pop-Location
    }
} finally {
    foreach ($name in $environmentNames) {
        [Environment]::SetEnvironmentVariable($name, $savedEnvironment[$name], 'Process')
    }
}
