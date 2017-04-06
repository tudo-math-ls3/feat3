@echo off

rem ===========================================================================
rem Check for Visual Studio version
echo Checking for Visual Studio installation...
if "%VS140COMNTOOLS%" neq "" (
  echo Found Visual Studio 2015 installation
  set VSVER=14
) else if "%VS120COMNTOOLS%" neq "" (
  echo Found Visual Studio 2013 installation
  set VSVER=12
)

rem Ensure that we have a path to devenv.exe
if "%VSVER%" == "" goto novs

rem Call VC batch script
call "%%VS%VSVER%0COMNTOOLS%%..\..\VC\vcvarsall.bat" amd64

rem ===========================================================================
echo **************************************************************************
if exist "./ALGLIB" (
  call ./vc_internal/make_alglib_vc%VSVER%.cmd dbg x64
  call ./vc_internal/make_alglib_vc%VSVER%.cmd opt x64
) else (
  echo ALGLIB not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./SuiteSparse" (
  call ./vc_internal/make_umfpack_vc%VSVER%.cmd dbg x64
  call ./vc_internal/make_umfpack_vc%VSVER%.cmd opt x64
) else (
  echo SuiteSparse not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./parmetis" (
  call ./vc_internal/make_parmetis_vc%VSVER%.cmd dbg x64
  call ./vc_internal/make_parmetis_vc%VSVER%.cmd opt x64
) else (
  echo ParMETIS not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./fparser" (
  call ./vc_internal/make_fparser_vc%VSVER%.cmd dbg x64
  call ./vc_internal/make_fparser_vc%VSVER%.cmd opt x64
) else (
  echo fparser not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./zlib" (
  call ./vc_internal/make_zlib_vc%VSVER%.cmd dbg x64
  call ./vc_internal/make_zlib_vc%VSVER%.cmd opt x64
) else (
  echo zlib not found; skipping...
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