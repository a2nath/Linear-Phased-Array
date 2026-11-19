param([ValidateSet('Configure', 'Release', 'Debug', 'Help', 'HelpDebug', 'HelpRelease', 'SmokeDebug', 'SmokeRelease', 'ArchitecturesDebug', 'ArchitecturesRelease', 'ArgsDebug', 'ArgsRelease', 'CpuDebug', 'CpuRelease', 'RunCpuDebug', 'RunCpuRelease', 'RunCudaDebug', 'RunCudaRelease')]
    [string]$Mode = 'CpuDebug',
    [string]$BuildDir = '')

$root = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($BuildDir)) {
    $directoryName = if ($Mode -like '*Cpu*') { 'build-cpu-windows' } else { 'build-cuda-windows' }
    $BuildDir = Join-Path $root 'build' $directoryName
}
$cmakeCommand = Get-Command cmake -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
$cmake = if ($cmakeCommand) { $cmakeCommand.Source } else { 'C:\Program Files\CMake\bin\cmake.exe' }
if (-not (Test-Path -LiteralPath $cmake -PathType Leaf)) {
    throw 'CMake was not found. Install CMake or the Visual Studio CMake tools and add cmake to PATH.'
}
$toolkit = Join-Path $PSScriptRoot 'cuda-toolkit-12.9-windows'
$build = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($BuildDir)
$log = Join-Path $BuildDir "last_run.log"

if ($Mode -eq 'Configure') {
    & $cmake --version
    & "$toolkit\bin\nvcc.exe" --version
    $arguments = @('-S', $root, '-B', $build, '-G', 'Visual Studio 17 2022',
        '-A', 'x64', '-T', "cuda=$toolkit", '-DPHASEARRAY_BUILD_CUDA=ON',
        '-DBUILD_SHARED_LIBS=OFF', "-DCUDAToolkit_ROOT=$toolkit")
    $sources = @{
        'SFML' = 'sfml'; 'IMGUI' = 'imgui'; 'IMGUI-SFML' = 'imgui-sfml';
        'ARGPARSE' = 'argparse'; 'SPDLOG' = 'spdlog'; 'JSONCPP' = 'jsoncpp';
        'RAPIDJSON' = 'rapidjson'
    }
    foreach ($dependency in $sources.GetEnumerator()) {
        $source = Join-Path $PSScriptRoot "_deps\$($dependency.Value)-src"
        $arguments += "-DFETCHCONTENT_SOURCE_DIR_$($dependency.Key)=$source"
    }
    & $cmake @arguments 2>&1 | Out-File -Encoding utf8 $log
    $result = $LASTEXITCODE
} elseif ($Mode -like 'Cpu*') {
    $configuration = $Mode.Substring(3)
    # Configure every CPU build so a fresh clone needs only this script.
    & $cmake -S $root -B $build -G 'Visual Studio 17 2022' -A x64 `
        '-DPHASEARRAY_BUILD_CUDA=OFF' '-DBUILD_SHARED_LIBS=OFF' 2>&1 |
        ForEach-Object { $_ | Out-Host; $_ } | Out-File -Encoding utf8 $log
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $cmake --build $build --target PhaseArray --config $configuration --parallel 8 2>&1 |
        ForEach-Object { $_ | Out-Host; $_ } | Out-File -Append -Encoding utf8 $log
    $result = $LASTEXITCODE
} elseif ($Mode -like 'Run*') {
    $cpu = $Mode -like 'RunCpu*'
    $configuration = if ($cpu) { $Mode.Substring(6) } else { $Mode.Substring(7) }
    $target = if ($cpu) { 'PhaseArray' } else { 'CudaSolution' }
    $executable = Join-Path $build "bin\$configuration\$target.exe"
    $runDirectory = Join-Path $PSScriptRoot "actual-windows-$target-$configuration"
    New-Item -ItemType Directory -Force -Path $runDirectory | Out-Null
    $fontDirectory = Join-Path $runDirectory 'font'
    New-Item -ItemType Directory -Force -Path $fontDirectory | Out-Null
    Copy-Item -Path (Join-Path $root 'font\OpenSans-Light.ttf') -Destination $fontDirectory -Force
    $stdout = Join-Path $runDirectory 'stdout.log'
    $stderr = Join-Path $runDirectory 'stderr.log'
    $arguments = @('--timeslot', '3', '--frequency', '1900000000.0', '--bandwidth', '20000000.0',
        '--symbolrate', '3840000.0', '--blockpersymbol', '768', '--height', '200',
        '--ms_grx', '-2', '--system_noise', '5', '--mobile_stations', '15',
        '--base_stations', '3', '--slimit', '24.0', '--base_station_theta_c', '75.0,90.0,105.0',
        '--base_station_location', '250,200 500,180 750,200',
        '--mobile_station_location', '250,500 350,500 450,500 550,500 650,500 750,500 250,600 350,600 450,600 550,600 650,600 750,600 450,325 500,325 550,325',
        '--base_station_antenna_counts', '5,3,5', '--base_station_power_range_dBm', '-30,30',
        '--base_station_scan_angle_range_deg', '-90,90', '--antenna_spacing', '0.3,0.4,0.3',
        '--antenna_dims', '0.2,0.4', '--antenna_txpower', '30.0,29.0,28.0',
        '--scan_angle', '40.0,65.0,-40.0', '--ms_selection', '7,4,10',
        '--output_dir', $runDirectory)
    if ($cpu) { $arguments += '--plot' }
    # These fixed arguments contain no embedded quotes or trailing backslashes.
    $argumentLine = ($arguments | ForEach-Object { '"' + $_ + '"' }) -join ' '
    Add-Type -TypeDefinition @'
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;
public static class PhaseArrayVerificationWindows {
    public delegate bool EnumCallback(IntPtr window, IntPtr parameter);
    [DllImport("user32.dll")] static extern bool EnumWindows(EnumCallback callback, IntPtr parameter);
    [DllImport("user32.dll")] static extern uint GetWindowThreadProcessId(IntPtr window, out uint processId);
    [DllImport("user32.dll", CharSet=CharSet.Unicode)] static extern int GetWindowText(IntPtr window, StringBuilder title, int capacity);
    [DllImport("user32.dll")] static extern bool PostMessage(IntPtr window, uint message, IntPtr wParam, IntPtr lParam);
    public static string[] Inspect(int processId, bool closeSimulationWindow) {
        var titles = new List<string>();
        EnumWindows((window, parameter) => {
            uint owner;
            GetWindowThreadProcessId(window, out owner);
            if (owner == processId) {
                var title = new StringBuilder(512);
                GetWindowText(window, title, title.Capacity);
                string text = title.ToString();
                if (text.Length > 0) titles.Add(text);
                if (closeSimulationWindow && text == "SFML Grid Plot")
                    PostMessage(window, 0x0010, IntPtr.Zero, IntPtr.Zero);
            }
            return true;
        }, IntPtr.Zero);
        return titles.ToArray();
    }
}
'@
    $runStarted = Get-Date
    $testProcess = Start-Process -FilePath $executable -ArgumentList $argumentLine `
        -WorkingDirectory $runDirectory -RedirectStandardOutput $stdout `
        -RedirectStandardError $stderr -PassThru
    $null = $testProcess.Handle
    $watch = [System.Diagnostics.Stopwatch]::StartNew()
    $titles = @()
    while (-not $testProcess.HasExited -and $watch.Elapsed.TotalSeconds -lt 20) {
        if (-not $cpu) {
            $titles += [PhaseArrayVerificationWindows]::Inspect($testProcess.Id, $watch.Elapsed.TotalSeconds -gt 3)
        }
        Start-Sleep -Milliseconds 200
        $testProcess.Refresh()
    }
    $timedOut = -not $testProcess.HasExited
    if ($timedOut) { Stop-Process -Id $testProcess.Id -Force }
    $testProcess.WaitForExit()
    $result = if ($timedOut) { 124 } else { $testProcess.ExitCode }
    $textOutput = Get-Content -Raw $stdout
    $record = [ordered]@{
        executable = $executable; configuration = $configuration; target = $target;
        arguments = $arguments; exit_code = $result; timed_out = $timedOut;
        simulation_result = $textOutput -match 'timeslot:.*sinr:';
        png_count = @(Get-ChildItem $runDirectory -Filter '*.png' | Where-Object { $_.LastWriteTime -ge $runStarted }).Count;
        observed_windows = @($titles | Select-Object -Unique);
        stdout = $stdout; stderr = $stderr;
        binary_sha256 = (Get-FileHash -Algorithm SHA256 $executable).Hash
    }
    $record | ConvertTo-Json -Depth 4 | Out-File -Encoding utf8 (Join-Path $runDirectory 'result.json')
    $record | ConvertTo-Json -Depth 4
} elseif ($Mode -like 'Help*') {
    $configurations = if ($Mode -eq 'Help') { @('Debug', 'Release') } else { @($Mode.Substring(4)) }
    foreach ($configuration in $configurations) {
        & "$build\bin\$configuration\CudaSolution.exe" --help
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    }
    $result = 0
} elseif ($Mode -like 'Args*') {
    $configuration = $Mode.Substring(4)
    $smokeSource = Join-Path $PSScriptRoot 'verify-cuda-smoke-src'
    $smokeBuild = Join-Path $PSScriptRoot 'verify-cuda-smoke-windows'
    & $cmake -S $smokeSource -B $smokeBuild -G 'Visual Studio 17 2022' -A x64 `
        "-DCMAKE_CONFIGURATION_TYPES=$configuration" "-DSMOKE_CONFIGURATION=$configuration" 2>&1 |
        Out-File -Encoding utf8 $log
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $cmake --build $smokeBuild --target PhaseArrayArgsSmoke --config $configuration --parallel 8 2>&1 |
        Out-File -Append -Encoding utf8 $log
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    $output = Join-Path $PSScriptRoot "verify-args-$configuration-output"
    & "$smokeBuild\bin\$configuration\PhaseArrayArgsSmoke.exe" `
        --frequency 1900000000 --bandwidth 20000000 --symbolrate 3840000 --blockpersymbol 768 `
        --height 200 --ms_grx -2 --system_noise 5 --base_stations 2 --mobile_stations 3 --timeslots 2 `
        --base_station_theta_c '90,90' --base_station_location '2,2 4,2' `
        --mobile_station_location '6,6 5,7 7,5' --base_station_antenna_counts '3,3' `
        --antenna_dims '0.2,0.4' --antenna_txpower '30,20 10,0' --scan_angle '30,-30 60,0' `
        --ms_selection '1,2 1,3' --area '16,16' --output_dir $output 2>&1 |
        Tee-Object -FilePath (Join-Path $PSScriptRoot "verify-cuda-windows-$Mode-output.log")
    $result = $LASTEXITCODE
} elseif ($Mode -like 'Smoke*') {
    $configuration = $Mode.Substring(5)
    $smokeSource = Join-Path $PSScriptRoot 'verify-cuda-smoke-src'
    $smokeBuild = Join-Path $PSScriptRoot 'verify-cuda-smoke-windows'
    & $cmake -S $smokeSource -B $smokeBuild -G 'Visual Studio 17 2022' -A x64 `
        "-DCMAKE_CONFIGURATION_TYPES=$configuration" "-DSMOKE_CONFIGURATION=$configuration" 2>&1 |
        Out-File -Encoding utf8 $log
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $cmake --build $smokeBuild --config $configuration --parallel 8 2>&1 |
        Out-File -Append -Encoding utf8 $log
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & "$smokeBuild\bin\$configuration\PhaseArrayCudaSmoke.exe" 2>&1 |
        Tee-Object -FilePath (Join-Path $PSScriptRoot "verify-cuda-windows-$Mode-output.log")
    $result = $LASTEXITCODE
} elseif ($Mode -like 'Architectures*') {
    $configuration = $Mode.Substring(13)
    & "$toolkit\bin\cuobjdump.exe" --list-elf "$build\bin\$configuration\CudaSolution.exe" 2>&1 |
        Tee-Object -FilePath $log
    $result = $LASTEXITCODE
} else {
    & $cmake --build $build --target CudaSolution --config $Mode --parallel 8 2>&1 |
        Out-File -Encoding utf8 $log
    $result = $LASTEXITCODE
}
Write-Output "$Mode exit code: $result; log: $log"
exit $result
