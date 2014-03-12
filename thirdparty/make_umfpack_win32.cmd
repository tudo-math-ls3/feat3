@echo off

rem ===========================================================================
rem Check for Visual Studio version
echo Checking for Visual Studio installation
if "%VS120COMNTOOLS%" neq "" set VSPATH=%VS120COMNTOOLS%..\..

rem Ensure that we have a path to devenv.exe
if "%VSPATH%" == "" goto novs

rem Call VC batch script
call "%VSPATH%\VC\vcvarsall.bat" x86

rem ===========================================================================

call ./make_umfpack_internal.cmd dbg x86
call ./make_umfpack_internal.cmd opt x86

goto end

rem ===========================================================================
:novs

echo.
echo ERROR: No compatible Visual Studio installation found
echo.
pause
goto end

rem ===========================================================================
:end