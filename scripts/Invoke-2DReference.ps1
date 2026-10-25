param(
  [Parameter(Mandatory=$true)][string]$SourceWorkbook,
  [Parameter(Mandatory=$true)][string]$Fixture,
  [Parameter(Mandatory=$true)][string]$OutDirectory
)
$ErrorActionPreference='Stop'
$application=$null;$book=$null;$previousSecurity=$null
$source=(Resolve-Path -LiteralPath $SourceWorkbook).Path
$fixturePath=(Resolve-Path -LiteralPath $Fixture).Path
$sourceDll=Join-Path (Split-Path -Parent $source) 'soilfem_native_x64.dll'
$sourceHash=(Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash
$dllHash=(Get-FileHash -LiteralPath $sourceDll -Algorithm SHA256).Hash
$fixtureHash=(Get-FileHash -LiteralPath $fixturePath -Algorithm SHA256).Hash
$runnerHash=(Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash
if(Get-Process EXCEL -ErrorAction SilentlyContinue){throw 'An Excel instance already exists; no workbook was opened.'}
if(Test-Path -LiteralPath $OutDirectory){throw 'Output directory already exists.'}
[void](New-Item -ItemType Directory -Path $OutDirectory)
$output=(Resolve-Path -LiteralPath $OutDirectory).Path
$copy=Join-Path $output 'case.xlsm'
Copy-Item -LiteralPath $source -Destination $copy
Copy-Item -LiteralPath $sourceDll -Destination (Join-Path $output 'soilfem_native_x64.dll')
$data=Get-Content -LiteralPath $fixturePath -Raw -Encoding UTF8 | ConvertFrom-Json
function Write-Json([string]$Path,$Value){
  [IO.File]::WriteAllText($Path,($Value|ConvertTo-Json -Depth 30),[Text.UTF8Encoding]::new($false))
}
function Code-Fingerprint($Workbook){
  $result=[ordered]@{}
  foreach($component in $Workbook.VBProject.VBComponents){
    $count=$component.CodeModule.CountOfLines
    $text=if($count -gt 0){$component.CodeModule.Lines(1,$count)}else{''}
    $hash=[Convert]::ToHexString([Security.Cryptography.SHA256]::HashData([Text.Encoding]::UTF8.GetBytes($text))).ToLowerInvariant()
    $result[$component.Name]=[ordered]@{lines=$count;sha256=$hash;component_type=[int]$component.Type}
  }
  return $result
}
function Put-Rows([string]$Name,$Rows,[int]$Columns){
  $sheet=$book.Worksheets.Item($Name);$sheet.Unprotect()
  $values=[object[,]]::new($Rows.Count,$Columns)
  for($r=0;$r -lt $Rows.Count;$r++){
    if($Rows[$r].Count -ne $Columns){throw "Invalid fixture column count in $Name"}
    for($c=0;$c -lt $Columns;$c++){$values[$r,$c]=$Rows[$r][$c]}
  }
  $sheet.Range($sheet.Cells(2,1),$sheet.Cells($Rows.Count+1,$Columns)).Value2=$values
}
function Read-Rows($Matrix){
  $rows=[Collections.Generic.List[object]]::new()
  for($r=$Matrix.GetLowerBound(0);$r -le $Matrix.GetUpperBound(0);$r++){
    $row=[object[]]::new($Matrix.GetLength(1))
    for($c=$Matrix.GetLowerBound(1);$c -le $Matrix.GetUpperBound(1);$c++){$row[$c-$Matrix.GetLowerBound(1)]=$Matrix.GetValue($r,$c)}
    $rows.Add($row)
  }
  return ,$rows.ToArray()
}
try{
  $application=New-Object -ComObject Excel.Application
  $application.Visible=$false;$application.EnableEvents=$false;$application.DisplayAlerts=$false
  $previousSecurity=$application.AutomationSecurity
  # Instance-scoped execution of this pinned, dedicated test copy. No Office
  # registry, Trust Center setting or existing workbook is changed.
  $application.AutomationSecurity=1
  $book=$application.Workbooks.Open($copy,0,$false)
  $prefix="'"+$book.Name+"'!"
  $before=Code-Fingerprint $book
  Write-Json (Join-Path $output 'code-before.json') $before
  [void]$application.Run($prefix+'FEMWriteBuildStamp',$false)
  [void]$application.Run($prefix+'SetP0SilentMode',$true)
  foreach($name in '材料データ','節点データ','要素データ','ステージ','載荷','接合','結果節点','結果要素','ステージ結果','結果増分'){
    if($name -eq '結果増分'){
      try{$sheet=$book.Worksheets.Item($name)}catch{continue}
    }else{$sheet=$book.Worksheets.Item($name)}
    $sheet.Unprotect()
    $last=[Math]::Max(2000,$sheet.UsedRange.Rows.Count)
    $sheet.Range($sheet.Cells(2,1),$sheet.Cells($last,[Math]::Max(72,$sheet.UsedRange.Columns.Count))).ClearContents()
  }
  $shape=$book.Worksheets.Item('形状入力');$shape.Unprotect()
  $shape.Range('B8:D2000').ClearContents();$shape.Range('G8:J2000').ClearContents()
  Put-Rows '材料データ' $data.materials 13
  Put-Rows '節点データ' $data.nodes 9
  Put-Rows '要素データ' $data.elements 10
  Put-Rows 'ステージ' $data.stages 6
  Put-Rows '載荷' $data.loads 10
  foreach($property in $data.settings.PSObject.Properties){
    if($property.Value -is [string]){[void]$application.Run($prefix+'FEMWriteTextSetting',$property.Name,$property.Value)}
    else{[void]$application.Run($prefix+'FEMWriteNumericSetting',$property.Name,[double]$property.Value)}
  }
  [void]$application.Run($prefix+'FEMInvalidateInputCache')
  [void]$application.Run($prefix+'FEMInvalidateSettingCache')
  $book.Save()
  $timer=[Diagnostics.Stopwatch]::StartNew()
  # Published program entry point. No VBE command bars, UI clicks or source
  # patches are used; VBA compiles its own required code normally on Run.
  [void]$application.Run($prefix+'P0_RunAnalysis')
  $timer.Stop()
  $after=Code-Fingerprint $book
  Write-Json (Join-Path $output 'code-after.json') $after
  $programChanges=@();$emptyDocumentChanges=@()
  foreach($name in (@($before.Keys)+@($after.Keys)|Sort-Object -Unique)){
    if($before.Contains($name) -and $after.Contains($name)){
      if($before[$name].sha256 -cne $after[$name].sha256 -or $before[$name].component_type -ne $after[$name].component_type){$programChanges+=,$name}
    }else{
      $entry=if($before.Contains($name)){$before[$name]}else{$after[$name]}
      # Creating an output worksheet also creates an empty document module.
      # It adds no VBA code. Keep all pre-existing/nonempty program text exact.
      if($entry.lines -eq 0 -and $entry.component_type -eq 100){$emptyDocumentChanges+=,$name}
      else{$programChanges+=,$name}
    }
  }
  Write-Json (Join-Path $output 'code-comparison.json') ([ordered]@{program_changes=$programChanges;empty_document_module_changes=$emptyDocumentChanges})
  if($programChanges.Count -gt 0){throw ('VBA program source changed: '+($programChanges -join ', '))}
  $tables=[ordered]@{};$absentOptional=@()
  foreach($entry in @(@('節点データ',9),@('要素データ',10),@('材料データ',13),@('ステージ',6),@('載荷',10),@('ステージ結果',20),@('結果節点',12),@('結果要素',72),@('結果増分',20),@('診断',5),@('設定',5))){
    if($entry[0] -eq '結果増分'){
      try{$sheet=$book.Worksheets.Item($entry[0])}catch{$absentOptional+=,$entry[0];continue}
    }else{$sheet=$book.Worksheets.Item($entry[0])}
    $last=$sheet.Cells($sheet.Rows.Count,1).End(-4162).Row
    $tables[$entry[0]]=Read-Rows $sheet.Range($sheet.Cells(1,1),$sheet.Cells([Math]::Max(2,$last),[int]$entry[1])).Value2
  }
  $status=[string]$book.Worksheets.Item('診断').Cells(2,2).Value2
  $record=[ordered]@{status=$status;elapsed_seconds=$timer.Elapsed.TotalSeconds;excel_version=$application.Version;
    source_workbook_sha256=$sourceHash;source_dll_sha256=$dllHash;fixture_sha256=$fixtureHash;runner_sha256=$runnerHash;
    absent_optional_output_sheets=$absentOptional;
    application_code_unchanged=$true;vba_modules=$after;tables=$tables;
    empty_document_module_changes=$emptyDocumentChanges;
    scope='Actual pinned 2D workbook execution through P0_RunAnalysis; no UI/source modification or numerical-setting relaxation'}
  Write-Json (Join-Path $output 'result.json') $record
  $book.Save()
  if((Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash -cne $sourceHash -or
     (Get-FileHash -LiteralPath $sourceDll -Algorithm SHA256).Hash -cne $dllHash -or
     (Get-FileHash -LiteralPath $fixturePath -Algorithm SHA256).Hash -cne $fixtureHash -or
     (Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash -cne $runnerHash){throw 'Reference input or runner bytes changed.'}
  if($status -ne 'PASS'){throw "2D reference analysis did not pass: $status"}
  Write-Output "2D reference PASS: $($data.nodes.Count) nodes, $($data.elements.Count) Q8 elements"
}catch{
  Write-Json (Join-Path $output 'failure.json') ([ordered]@{error=$_.Exception.Message;source_workbook_sha256=$sourceHash;fixture_sha256=$fixtureHash;runner_sha256=$runnerHash})
  throw
}finally{
  if($book){try{$book.Close($false)}catch{}}
  if($application){
    if($null -ne $previousSecurity){$application.AutomationSecurity=$previousSecurity}
    if($application.Workbooks.Count -eq 0){$application.Quit()}
    else{$application.EnableEvents=$true;$application.DisplayAlerts=$true;$application.Visible=$true}
    [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($application)
  }
  [GC]::Collect();[GC]::WaitForPendingFinalizers()
}
