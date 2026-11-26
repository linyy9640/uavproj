# 编译当前 .cpp（自动识别源码编码），输出同名 .exe
# 供 .vscode/tasks.json 的 "build: 编译当前文件" 任务调用
param(
    [Parameter(Mandatory = $true)]
    [string]$Source
)

$gpp = 'E:\DEVC\Dev-Cpp\MinGW64\bin\g++.exe'

if (-not (Test-Path -LiteralPath $Source)) {
    Write-Host "[build] 找不到源文件: $Source" -ForegroundColor Red
    exit 1
}
if (-not (Test-Path -LiteralPath $gpp)) {
    Write-Host "[build] 找不到编译器: $gpp" -ForegroundColor Red
    exit 1
}

# 编码识别：能被严格 UTF-8 解码就按 UTF-8 处理，否则按系统 ANSI(GBK) 处理。
# 两种情况都把执行字符集设为 GBK，这样中文在系统控制台(代码页 936)里显示正常。
$bytes  = [System.IO.File]::ReadAllBytes($Source)
$strict = New-Object System.Text.UTF8Encoding($false, $true)
$isUtf8 = $true
try { $null = $strict.GetString($bytes) } catch { $isUtf8 = $false }

if ($isUtf8) {
    $charset = @('-finput-charset=UTF-8', '-fexec-charset=GBK')
    $encName = 'UTF-8'
} else {
    $charset = @('-finput-charset=GBK', '-fexec-charset=GBK')
    $encName = 'GBK/ANSI'
}

$output   = [System.IO.Path]::ChangeExtension($Source, '.exe')
$compArgs = @('-std=c++11', '-g', '-Wall') + $charset + @($Source, '-o', $output)

Write-Host "[build] 源码编码: $encName" -ForegroundColor DarkGray
Write-Host "[build] $gpp $($compArgs -join ' ')" -ForegroundColor DarkGray

& $gpp @compArgs
$code = $LASTEXITCODE

if ($code -eq 0) {
    Write-Host "[build] 编译成功 -> $output" -ForegroundColor Green
} else {
    Write-Host "[build] 编译失败，退出码 $code" -ForegroundColor Red
}
exit $code
