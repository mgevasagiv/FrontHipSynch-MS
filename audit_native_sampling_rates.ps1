# audit_native_sampling_rates.ps1
# Scan E:\MAZE_1.0\data\iEEG\raw_orig for native acquisition sampling rates.
# No MATLAB required. Paste the printed table back into the Cursor chat.
#
# Usage (PowerShell):
#   cd C:\Users\mgeva\Documents\GitHub\Kamin_iEEG_MAZE_codes
#   powershell -ExecutionPolicy Bypass -File .\audit_native_sampling_rates.ps1

$ErrorActionPreference = 'Continue'
$Root = 'E:\MAZE_1.0\data\iEEG\raw_orig'
$Subjects = @('s01','s06','s08','s09','s12','s14','s15','s16','s17','s18')

function Get-NlxSamplingFrequency([string]$Path) {
    try {
        $fs = [System.IO.File]::OpenRead($Path)
        $buf = New-Object byte[] 16384
        [void]$fs.Read($buf, 0, 16384)
        $fs.Close()
        $text = [System.Text.Encoding]::ASCII.GetString($buf)
        if ($text -match '-SamplingFrequency\s+([0-9.]+)') {
            return [double]$Matches[1]
        }
    } catch {}
    return $null
}

function Get-EdfSamplingFrequency([string]$Path) {
    try {
        $fs = [System.IO.File]::OpenRead($Path)
        $hdr = New-Object byte[] 256
        [void]$fs.Read($hdr, 0, 256)
        $text = [System.Text.Encoding]::ASCII.GetString($hdr)
        $nchan = [int]($text.Substring(252,4).Trim())
        $dur = [double]($text.Substring(244,8).Trim())
        $skip = 256 + $nchan*(16+80+8+8+8+8+8+80)
        [void]$fs.Seek($skip, 'Begin')
        $sampBuf = New-Object byte[] 8
        [void]$fs.Read($sampBuf, 0, 8)
        $fs.Close()
        $ns = [double]([System.Text.Encoding]::ASCII.GetString($sampBuf).Trim())
        if ($dur -gt 0 -and $ns -gt 0) { return ($ns / $dur) }
    } catch {}
    return $null
}

Write-Host "Scanning: $Root"
if (-not (Test-Path $Root)) {
    Write-Error "Path not found: $Root"
    exit 1
}

Write-Host "`nTop-level folders:"
Get-ChildItem $Root -Directory | ForEach-Object { Write-Host ("  " + $_.Name) }

$rows = @()
foreach ($subj in $Subjects) {
    $dirs = Get-ChildItem $Root -Directory -Recurse -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -match ("(?i)^" + [regex]::Escape($subj) + "$") -or $_.Name -match ("(?i)" + [regex]::Escape($subj)) } |
        Select-Object -ExpandProperty FullName -Unique

    # Also check neuralynx\<subj>
    $nlx = Join-Path $Root "neuralynx\$subj"
    if (Test-Path $nlx) { $dirs = @($dirs + $nlx) | Select-Object -Unique }

    $rates = @()
    $evidence = @()
    foreach ($d in $dirs) {
        Get-ChildItem $d -Recurse -Include *.ncs,*.NCS -ErrorAction SilentlyContinue |
            Select-Object -First 5 | ForEach-Object {
                $r = Get-NlxSamplingFrequency $_.FullName
                if ($null -ne $r) { $rates += $r; $evidence += ("ncs:" + $_.Name + "=" + $r) }
            }
        Get-ChildItem $d -Recurse -Include *.edf,*.EDF -ErrorAction SilentlyContinue |
            Select-Object -First 3 | ForEach-Object {
                $r = Get-EdfSamplingFrequency $_.FullName
                if ($null -ne $r) { $rates += $r; $evidence += ("edf:" + $_.Name + "=" + $r) }
            }
    }

    $native = $null
    if ($rates.Count -gt 0) {
        $grouped = $rates | Group-Object | Sort-Object Count -Descending
        $native = [double]$grouped[0].Name
    }

    $flag = 'unknown'
    if ($null -ne $native) {
        if ($native -ge 500 -and $native -le 520) { $flag = 'impacted_native_~512Hz' }
        elseif ($native -ge 900) { $flag = 'not_impacted_native_>=1kHz' }
        else { $flag = 'check_manually' }
    }

    $rows += [pscustomobject]@{
        subj = $subj
        native_fs_Hz = $native
        reviewer_AA_concern = $flag
        raw_folders = ($dirs -join ' | ')
        evidence = ($evidence -join '; ')
    }
}

Write-Host "`n=== Per-subject native sampling rates ==="
$rows | Format-Table -AutoSize | Out-String | Write-Host

$impacted = @($rows | Where-Object { $_.reviewer_AA_concern -eq 'impacted_native_~512Hz' })
$ok = @($rows | Where-Object { $_.reviewer_AA_concern -eq 'not_impacted_native_>=1kHz' })
$unk = @($rows | Where-Object { $_.reviewer_AA_concern -in @('unknown','check_manually') })

Write-Host "=== Summary (N=$($Subjects.Count)) ==="
Write-Host ("Native ~512 Hz (reviewer AA concern): {0}" -f $impacted.Count)
Write-Host ("Native >= ~1 kHz: {0}" -f $ok.Count)
Write-Host ("Unknown / manual check: {0}" -f $unk.Count)
if ($impacted.Count -gt 0) { Write-Host ("Impacted: " + (($impacted | ForEach-Object subj) -join ', ')) }
if ($ok.Count -gt 0) { Write-Host ("Higher-rate: " + (($ok | ForEach-Object subj) -join ', ')) }

$outDir = 'E:\MAZE_1.0\analysis\sampling_rate_audit'
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
$csv = Join-Path $outDir 'sampling_rate_audit_raw_orig.csv'
$rows | Export-Csv -NoTypeInformation -Path $csv
Write-Host "`nWrote $csv"
Write-Host "Paste this table back into the Cursor chat."
