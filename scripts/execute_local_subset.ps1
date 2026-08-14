param(
  [string]$ProjectRoot = 'D:\Xiaonan\CODEX_projects\Yanan_Xenium',
  [string]$RunLabel = 'local_notebook_test'
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Join-Path $ProjectRoot 'adipose_analysis\YNH_Xenium_scWAT'
$InputRoot = Join-Path $ProjectRoot 'adipose_analysis\subset_input\adipose_data'
$RunRoot = Join-Path $ProjectRoot "adipose_analysis\scwat_qc_outputs\$RunLabel"
$TempRoot = Join-Path $ProjectRoot 'adipose_analysis\tmp'
$ExecutedRoot = Join-Path $RunRoot 'executed_notebooks'
$Python = 'C:\Users\Xiaonan_Wang\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$Rscript = 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe'

foreach ($Path in @($RepoRoot, $InputRoot)) {
  if (-not (Test-Path -LiteralPath $Path)) { throw "Required path is absent: $Path" }
}
foreach ($Path in @($RunRoot, $TempRoot, $ExecutedRoot)) {
  New-Item -ItemType Directory -Path $Path -Force | Out-Null
}
$env:TMPDIR = $TempRoot
$env:TEMP = $TempRoot
$env:TMP = $TempRoot

& $Rscript (Join-Path $RepoRoot 'tests\test_source.R')
if ($LASTEXITCODE -ne 0) { throw 'Reusable R tests failed.' }

$Renderer = Join-Path $RepoRoot 'scripts\render_notebooks.py'
$Executor = Join-Path $RepoRoot 'scripts\execute_r_notebook.R'
$SectionTemplate = Join-Path $RepoRoot 'notebooks\01_section_phase0_2_QC.ipynb'

foreach ($Region in @('Region_1', 'Region_2', 'Region_3', 'Region_4')) {
  $Injected = Join-Path $ExecutedRoot "$Region.input.ipynb"
  $Executed = Join-Path $ExecutedRoot "$Region.executed.ipynb"
  & $Python $Renderer --inject $SectionTemplate $Injected `
    --set "PROJECT_ROOT=$($ProjectRoot -replace '\\','/')" `
    --set "PIPELINE_REPO=$($RepoRoot -replace '\\','/')" `
    --set "INPUT_ROOT=$($InputRoot -replace '\\','/')" `
    --set "REGION_ID=$Region" --set "RUN_LABEL=$RunLabel" `
    --set 'METADATA_PATH=' --set 'EXPECTED_SECTION_COUNT=4L' `
    --set 'SEED=20260814L' --set 'STRICT_MODE=FALSE'
  if ($LASTEXITCODE -ne 0) { throw "Parameter injection failed for $Region." }
  & $Rscript $Executor $Injected $Executed
  if ($LASTEXITCODE -ne 0) { throw "Notebook execution failed for $Region." }
  $SectionOutput = Join-Path $RunRoot "sections\$Region"
  & $Rscript -e "source('$($RepoRoot -replace '\\','/')/R/source.R'); stopifnot(validate_section_artifacts('$($SectionOutput -replace '\\','/')', '$Region')); q <- read.delim(gzfile('$($SectionOutput -replace '\\','/')/cell_qc_metadata.tsv.gz')); stopifnot(nrow(q) == 500L)"
  if ($LASTEXITCODE -ne 0) { throw "Artifact validation failed for $Region." }
}

$SummaryTemplate = Join-Path $RepoRoot 'notebooks\02_slide_QC_summary.ipynb'
$SummaryInjected = Join-Path $ExecutedRoot 'slide_summary.input.ipynb'
$SummaryExecuted = Join-Path $ExecutedRoot 'slide_summary.executed.ipynb'
& $Python $Renderer --inject $SummaryTemplate $SummaryInjected `
  --set "PROJECT_ROOT=$($ProjectRoot -replace '\\','/')" `
  --set "PIPELINE_REPO=$($RepoRoot -replace '\\','/')" `
  --set "RUN_LABEL=$RunLabel" --set 'EXPECTED_SECTION_COUNT=4L'
if ($LASTEXITCODE -ne 0) { throw 'Summary parameter injection failed.' }
& $Rscript $Executor $SummaryInjected $SummaryExecuted
if ($LASTEXITCODE -ne 0) { throw 'Summary notebook execution failed.' }

& $Rscript -e "x <- readRDS('$($RunRoot -replace '\\','/')/slide_summary/slide_qc_summary.rds'); stopifnot(nrow(x[['data']][['cell_metadata']]) == 2000L, length(unique(x[['data']][['cell_metadata']][['region_id']])) == 4L)"
if ($LASTEXITCODE -ne 0) { throw 'Slide summary reload validation failed.' }
Write-Output "Local four-section notebook validation complete: $RunRoot"
