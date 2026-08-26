param(
  [string]$ProjectRoot = 'D:\Xiaonan\CODEX_projects\Yanan_Xenium',
  [string]$RunLabel = 'local_colon_parity_qc',
  [string]$RepoRoot = '',
  [ValidateSet('ALL', 'SUMMARY', 'Region_1', 'Region_2', 'Region_3', 'Region_4')]
  [string]$RegionId = 'ALL'
)

$ErrorActionPreference = 'Stop'
if (-not $RepoRoot) { $RepoRoot = Join-Path $ProjectRoot 'adipose_analysis\YNH_Xenium_scWAT' }
$InputRoot = Join-Path $ProjectRoot 'adipose_analysis\subset_input\adipose_data'
$RunRoot = Join-Path $ProjectRoot "adipose_analysis\scwat_qc_outputs\$RunLabel"
$TempRoot = Join-Path $ProjectRoot "adipose_analysis\tmp\$RunLabel"
$ExecutedRoot = Join-Path $RunRoot 'executed_notebooks'
$Rscript = 'D:\Programs\R-4.6.1\bin\Rscript.exe'
$RLibrary = 'D:\Programs\R_library'

foreach ($Path in @($RepoRoot, $InputRoot, $Rscript, $RLibrary)) {
  if (-not (Test-Path -LiteralPath $Path)) { throw "Required path is absent: $Path" }
}
foreach ($Path in @($RunRoot, $TempRoot, $ExecutedRoot)) {
  New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

$env:TMPDIR = $TempRoot
$env:TEMP = $TempRoot
$env:TMP = $TempRoot
$env:R_LIBS_USER = $RLibrary
$env:SCWAT_QC_MODE = 'LOCAL_SUBSET'
$env:SCWAT_PROJECT_ROOT = $ProjectRoot
$env:SCWAT_PIPELINE_REPO = $RepoRoot
$env:SCWAT_INPUT_ROOT = $InputRoot
$env:SCWAT_RUN_LABEL = $RunLabel
$env:SCWAT_RUN_ROOT = $RunRoot
$env:SCWAT_TEMP_ROOT = $TempRoot

$Executor = Join-Path $RepoRoot 'scripts\execute_r_notebook.R'
$Regions = if ($RegionId -eq 'ALL') {
  @('Region_1', 'Region_2', 'Region_3', 'Region_4')
} elseif ($RegionId -eq 'SUMMARY') {
  @()
} else {
  @($RegionId)
}

foreach ($Region in $Regions) {
  $Notebook = Join-Path $RepoRoot "notebooks\01_QC_$($Region -replace '_','').ipynb"
  $Executed = Join-Path $ExecutedRoot "$Region.executed.ipynb"
  & $Rscript $Executor $Notebook $Executed
  if ($LASTEXITCODE -ne 0) { throw "Notebook execution failed for $Region." }
  $SectionOutput = Join-Path $RunRoot "sections\$Region"
  & $Rscript -e "source('$($RepoRoot -replace '\\','/')/R/source.R'); stopifnot(validate_scwat_region_qc_bundle('$($SectionOutput -replace '\\','/')', '$Region')); x <- readRDS('$($SectionOutput -replace '\\','/')/$Region.phase0_2_qc.rds'); expected <- with(x[['cells']], nFeature_Xenium > 5 & nFeature_Xenium < 200 & nCount_Xenium > 10 & nCount_Xenium < 1000); stopifnot(identical(as.logical(x[['cells']][['primary_include']]), expected), isTRUE(x[['raw_counts_preserved']]), x[['cells_deleted']] == 0L)"
  if ($LASTEXITCODE -ne 0) { throw "Artifact validation failed for $Region." }
}

if ($RegionId -in @('ALL', 'SUMMARY')) {
  $SummaryNotebook = Join-Path $RepoRoot 'notebooks\02_slide_QC_summary.ipynb'
  $SummaryExecuted = Join-Path $ExecutedRoot 'slide_summary.executed.ipynb'
  & $Rscript $Executor $SummaryNotebook $SummaryExecuted
  if ($LASTEXITCODE -ne 0) { throw 'Summary notebook execution failed.' }
  $SummaryOutput = Join-Path $RunRoot 'slide_summary'
  & $Rscript -e "source('$($RepoRoot -replace '\\','/')/R/source.R'); stopifnot(validate_scwat_slide_qc_bundle('$($SummaryOutput -replace '\\','/')')); x <- readRDS('$($SummaryOutput -replace '\\','/')/scwat_slide_qc.rds'); stopifnot(identical(x[['slide_data']][['coverage']][['region_id']], expected_scwat_regions()), nrow(x[['mouse_summary']][['mouse_summary']]) == 2L)"
  if ($LASTEXITCODE -ne 0) { throw 'Slide summary reload validation failed.' }
}

Write-Output "Local Colon-parity notebook validation complete for ${RegionId}: $RunRoot"
