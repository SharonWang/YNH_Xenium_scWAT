param(
  [string]$ProjectRoot = 'D:\Xiaonan\CODEX_projects\Yanan_Xenium',
  [string]$RunLabel = 'local_extended_qc_test',
  [string]$MetadataPath = '',
  [string]$RepoRoot = '',
  [ValidateSet('ALL', 'SUMMARY', 'Region_1', 'Region_2', 'Region_3', 'Region_4')]
  [string]$RegionId = 'ALL'
)

$ErrorActionPreference = 'Stop'
if (-not $RepoRoot) { $RepoRoot = Join-Path $ProjectRoot 'adipose_analysis\YNH_Xenium_scWAT' }
if (-not $MetadataPath) { $MetadataPath = Join-Path $RepoRoot 'config\scwat_sample_manifest.tsv' }
$InputRoot = Join-Path $ProjectRoot 'adipose_analysis\subset_input\adipose_data'
$RunRoot = Join-Path $ProjectRoot "adipose_analysis\scwat_qc_outputs\$RunLabel"
$TempRoot = Join-Path $ProjectRoot 'adipose_analysis\tmp'
$ExecutedRoot = Join-Path $RunRoot 'executed_notebooks'
$Python = 'C:\Users\Xiaonan_Wang\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$Rscript = 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe'

foreach ($Path in @($RepoRoot, $InputRoot, $MetadataPath)) {
  if (-not (Test-Path -LiteralPath $Path)) { throw "Required path is absent: $Path" }
}
foreach ($Path in @($RunRoot, $TempRoot, $ExecutedRoot)) {
  New-Item -ItemType Directory -Path $Path -Force | Out-Null
}
$env:TMPDIR = $TempRoot
$env:TEMP = $TempRoot
$env:TMP = $TempRoot
$env:PYTHONDONTWRITEBYTECODE = '1'

& $Rscript (Join-Path $RepoRoot 'tests\test_source.R')
if ($LASTEXITCODE -ne 0) { throw 'Reusable R tests failed.' }
& $Rscript (Join-Path $RepoRoot 'tests\test_extended_qc.R')
if ($LASTEXITCODE -ne 0) { throw 'Extended QC R tests failed.' }

$Renderer = Join-Path $RepoRoot 'scripts\render_notebooks.py'
$Executor = Join-Path $RepoRoot 'scripts\execute_r_notebook.R'
$Regions = if ($RegionId -eq 'ALL') { @('Region_1', 'Region_2', 'Region_3', 'Region_4') } elseif ($RegionId -eq 'SUMMARY') { @() } else { @($RegionId) }
$RunSummary = $RegionId -in @('ALL', 'SUMMARY')

foreach ($Region in $Regions) {
  $SectionTemplate = (& $Python $Renderer --region-notebook $Region).Trim()
  if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $SectionTemplate)) { throw "Named notebook resolution failed for $Region." }
  $Injected = Join-Path $ExecutedRoot "$Region.input.ipynb"
  $Executed = Join-Path $ExecutedRoot "$Region.executed.ipynb"
  & $Python $Renderer --inject $SectionTemplate $Injected `
    --set "PROJECT_ROOT=$($ProjectRoot -replace '\\','/')" `
    --set "PIPELINE_REPO=$($RepoRoot -replace '\\','/')" `
    --set "INPUT_ROOT=$($InputRoot -replace '\\','/')" `
    --set "REGION_ID=$Region" --set "RUN_LABEL=$RunLabel" `
    --set "METADATA_PATH=$($MetadataPath -replace '\\','/')" --set 'EXPECTED_SECTION_COUNT=4L' `
    --set 'SEED=20260814L' --set 'STRICT_MODE=FALSE' `
    --set 'EXTENDED_QC_MODE=LOCAL_SUBSET' `
    --set "EXTENDED_QC_CONFIG_PATH=$(($RepoRoot -replace '\\','/') + '/config/extended_qc_defaults.tsv')"
  if ($LASTEXITCODE -ne 0) { throw "Parameter injection failed for $Region." }
  & $Rscript $Executor $Injected $Executed
  if ($LASTEXITCODE -ne 0) { throw "Notebook execution failed for $Region." }
  $SectionOutput = Join-Path $RunRoot "sections\$Region"
  & $Rscript -e "source('$($RepoRoot -replace '\\','/')/R/source.R'); stopifnot(validate_section_artifacts('$($SectionOutput -replace '\\','/')', '$Region')); stopifnot(validate_extended_section_artifacts('$($SectionOutput -replace '\\','/')', '$Region', 'LOCAL_SUBSET')); q <- read.delim(gzfile('$($SectionOutput -replace '\\','/')/cell_qc_metadata.tsv.gz')); g <- read.delim('$($SectionOutput -replace '\\','/')/section_readiness_gates.tsv'); e <- read_extended_section_artifacts('$($SectionOutput -replace '\\','/')', '$Region', 'LOCAL_SUBSET'); stopifnot(nrow(q) == 500L, all(q[['metadata_status']] == 'VERIFIED_USER_SUPPLIED'), g[['status']][g[['gate']] == 'metadata'] == 'PASS', all(e[['gene_quality']][['transcript_status']] == 'NOT_RUN_LOCAL_SUBSET'), e[['status']][['cells_deleted']] == 0L, nrow(e[['spatial_cells']]) == 500L)"
  if ($LASTEXITCODE -ne 0) { throw "Artifact validation failed for $Region." }
}

if ($RunSummary) {
  $SummaryTemplate = Join-Path $RepoRoot 'notebooks\02_slide_QC_summary.ipynb'
  $SummaryInjected = Join-Path $ExecutedRoot 'slide_summary.input.ipynb'
  $SummaryExecuted = Join-Path $ExecutedRoot 'slide_summary.executed.ipynb'
  & $Python $Renderer --inject $SummaryTemplate $SummaryInjected `
    --set "PROJECT_ROOT=$($ProjectRoot -replace '\\','/')" `
    --set "PIPELINE_REPO=$($RepoRoot -replace '\\','/')" `
    --set "RUN_LABEL=$RunLabel" --set 'EXPECTED_SECTION_COUNT=4L' `
    --set "METADATA_PATH=$($MetadataPath -replace '\\','/')" `
    --set "EXTENDED_QC_CONFIG_PATH=$(($RepoRoot -replace '\\','/') + '/config/extended_qc_defaults.tsv')" `
    --set "SUBSET_REFERENCE_PATH=$(($RepoRoot -replace '\\','/') + '/config/subset_qc_reference.tsv')"
  if ($LASTEXITCODE -ne 0) { throw 'Summary parameter injection failed.' }
  & $Rscript $Executor $SummaryInjected $SummaryExecuted
  if ($LASTEXITCODE -ne 0) { throw 'Summary notebook execution failed.' }

  & $Rscript -e "source('$($RepoRoot -replace '\\','/')/R/source.R'); x <- readRDS('$($RunRoot -replace '\\','/')/slide_summary/slide_qc_summary.rds'); s <- read.delim('$($RunRoot -replace '\\','/')/slide_summary/combined_qc_summary.tsv', check.names=FALSE); r <- read.delim('$($RunRoot -replace '\\','/')/slide_summary/combined_readiness.tsv', check.names=FALSE); e <- read.delim('$($RunRoot -replace '\\','/')/slide_summary/extended_slide_qc_status.tsv', check.names=FALSE); rank <- read.delim('$($RunRoot -replace '\\','/')/slide_summary/subset_full_qc_ranking.tsv', check.names=FALSE); cand <- read.delim('$($RunRoot -replace '\\','/')/slide_summary/candidate_cycle_affected_genes.tsv', check.names=FALSE); stopifnot(nrow(x[['data']][['cell_metadata']]) == 2000L, length(unique(x[['data']][['cell_metadata']][['region_id']])) == 4L, identical(s[['core_qc_pass']], c(493L,500L,497L,499L)), identical(s[['review_flagged']], c(37L,16L,27L,11L)), identical(r[['status']], c('HOLD','HOLD','PASS','HOLD')), all(s[['cells_deleted']] == 0L), validate_extended_slide_qc_artifacts('$($RunRoot -replace '\\','/')', stop_on_error=TRUE), e[['subset_full_status']] == 'NOT_RUN_LOCAL_SUBSET', all(rank[['comparison_status']] == 'NOT_RUN_LOCAL_SUBSET'), all(cand[['candidate_status']] == 'CANDIDATE_NOT_CONFIRMED'), all(cand[['exact_cycle_status']] == 'REQUIRES_10X_DIAGNOSTICS'))"
  if ($LASTEXITCODE -ne 0) { throw 'Slide summary reload validation failed.' }
}
Write-Output "Local notebook validation complete for target ${RegionId}: $RunRoot"
