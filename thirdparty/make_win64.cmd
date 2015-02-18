@echo off

rem ===========================================================================
rem Check for Visual Studio version
echo Checking for Visual Studio installation...
if "%VS120COMNTOOLS%" neq "" set VSPATH=%VS120COMNTOOLS%..\..

rem Ensure that we have a path to devenv.exe
if "%VSPATH%" == "" goto novs

rem Call VC batch script
call "%VSPATH%\VC\vcvarsall.bat" amd64

rem ===========================================================================
echo **************************************************************************
if exist "./ALGLIB" (
  call ./vc12_internal/make_alglib_internal.cmd dbg x64
  call ./vc12_internal/make_alglib_internal.cmd opt x64
) else (
  echo ALGLIB not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./SuiteSparse" (
  call ./vc12_internal/make_suite_sparse_internal.cmd dbg x64
  call ./vc12_internal/make_suite_sparse_internal.cmd opt x64
) else (
  echo SuiteSparse not found; skipping...
  echo.
)

echo **************************************************************************

rem ===========================================================================
goto end
:novs

echo.
echo ERROR: No compatible Visual Studio installation found
echo.
pause
goto end

rem ===========================================================================
:end