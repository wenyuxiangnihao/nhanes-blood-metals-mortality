# =====================================================================
#  publish_to_github.ps1
#  Publish the NHANES blood-metal analysis code to GitHub and prepare
#  the version tag that Zenodo archives.
#
#  Usage (run inside the unzipped repository folder):
#      powershell -ExecutionPolicy Bypass -File .\publish_to_github.ps1 -GitHubUser <yourname>
#      powershell -ExecutionPolicy Bypass -File .\publish_to_github.ps1 -GitHubUser <yourname> -Push
#
#  Step 1 (no -Push): checks the folder, replaces the placeholder URL,
#                     creates the git repository, the first commit and tag v1.0.0.
#  Step 2 (-Push)   : pushes to the empty GitHub repository you created.
# =====================================================================
param(
    [Parameter(Mandatory = $true)][string]$GitHubUser,
    [string]$RepoName = "nhanes-blood-metals-mortality",
    [string]$CommitName  = "Yuxiang Wen",
    [string]$CommitEmail = "wenyuxiang@yangtzeu.edu.cn",
    [switch]$Push
)

$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent $MyInvocation.MyCommand.Path)
$URL = "https://github.com/$GitHubUser/$RepoName"

Write-Host "== 1. Environment check ==" -ForegroundColor Cyan
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw "git was not found in PATH. Install Git (or add D:\Doubao\tools\Git\cmd to PATH) and retry."
}
Write-Host ("git  : " + (git --version))

Write-Host "== 2. Data-leak check (the repository must contain code only) ==" -ForegroundColor Cyan
$bad = Get-ChildItem -Recurse -File | Where-Object {
    $_.Extension -in ".xpt", ".XPT", ".rds", ".RDS", ".dat", ".tif", ".eps" -or
    ($_.Extension -eq ".csv" -and $_.Name -ne "data_dictionary.csv")
}
if ($bad) {
    Write-Host "STOP - data-like files found. Remove them before publishing:" -ForegroundColor Red
    $bad | ForEach-Object { Write-Host ("   " + $_.FullName) }
    throw "Aborted: the repository must not contain participant data."
}
Write-Host "OK - only data_dictionary.csv is present." -ForegroundColor Green

Write-Host "== 3. Replace the placeholder repository URL ==" -ForegroundColor Cyan
foreach ($f in "CITATION.cff", ".zenodo.json", "README.md") {
    if (Test-Path $f) {
        $t = [IO.File]::ReadAllText((Resolve-Path $f))
        $t = $t.Replace("https://github.com/your-org/nhanes-blood-metals-mortality", $URL)
        [IO.File]::WriteAllText((Resolve-Path $f), $t, (New-Object Text.UTF8Encoding($false)))
        Write-Host "   updated $f"
    }
}

Write-Host "== 4. Create the local repository ==" -ForegroundColor Cyan
if (-not (Test-Path ".git")) {
    git init | Out-Null
    git branch -M main
    git config user.name  $CommitName
    git config user.email $CommitEmail
    Write-Host "   git initialised on branch main"
} else {
    Write-Host "   .git already exists - keeping it"
}
git add .
$staged = (git status --porcelain).Count
if ($staged -gt 0) {
    git commit -m "NHANES blood metal mixtures and mortality: analysis code for the manuscript (v1.0.0)" | Out-Null
    Write-Host "   commit created ($staged paths)"
} else {
    Write-Host "   nothing new to commit"
}
if (-not (git tag --list "v1.0.0")) {
    git tag -a v1.0.0 -m "Code for: Opposing directions in blood metal mixtures and mortality in US adults (NHANES 1999-2018)"
    Write-Host "   tag v1.0.0 created"
}

Write-Host "== 5. Push ==" -ForegroundColor Cyan
if ($Push) {
    $remotes = git remote
    if ($remotes -notcontains "origin") { git remote add origin "$URL.git" }
    git push -u origin main
    git push origin v1.0.0
    Write-Host "Pushed to $URL" -ForegroundColor Green
} else {
    Write-Host "Not pushed yet (no -Push). Create the EMPTY repository first:" -ForegroundColor Yellow
    Write-Host "   1) https://github.com/new"
    Write-Host "      Name: $RepoName   Visibility: Public   (do NOT add README/.gitignore/Licence)"
    Write-Host "   2) Re-run this script with -Push:"
    Write-Host "      powershell -ExecutionPolicy Bypass -File .\publish_to_github.ps1 -GitHubUser $GitHubUser -Push"
}

Write-Host ""
Write-Host "== 6. Zenodo archive (do this once) ==" -ForegroundColor Cyan
Write-Host "   1) https://zenodo.org  -> sign in with GitHub"
Write-Host "   2) Settings -> GitHub -> switch ON the repository '$RepoName'"
Write-Host "   3) GitHub -> Releases -> Draft a new release -> choose tag v1.0.0 -> Publish release"
Write-Host "   4) Back in Zenodo: the new record shows a DOI '10.5281/zenodo.XXXXXXX'"
Write-Host "      (use the Concept DOI - it always points to the latest version)"
Write-Host ""
Write-Host "== 7. Fill the DOI into the manuscript ==" -ForegroundColor Cyan
Write-Host "   Code availability statement to paste into the journal system:"
Write-Host ""
Write-Host "   Code availability. The complete analysis code (data download, cleaning,"
Write-Host "   statistical analysis and figure generation) is openly available at"
Write-Host "   $URL and archived on Zenodo (DOI: 10.5281/zenodo.XXXXXXX)."
Write-Host ""
Write-Host "   Send the repository URL and the DOI back to ima.copilot and the manuscript,"
Write-Host "   cover letter and submission files will be updated for you in one step."
