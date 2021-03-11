@echo off
@rem FEAT3: Finite Element Analysis Toolbox, Version 3
@rem Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
@rem FEAT3 is released under the GNU General Public License version 3,
@rem see the file 'copyright.txt' in the top level directory for details.

rem ===========================================================================
rem Check for Visual Studio version
echo Checking for Visual Studio installation...
if "%VS16PATH%" neq "" (
  echo Found Visual Studio 2019 installation
  set VSVER=16
  call "%VS16PATH%\Common7\Tools\VsDevCmd.bat" -arch=amd64 -host_arch=amd64
) else if "%VS140COMNTOOLS%" neq "" (
  echo Found Visual Studio 2015 installation
  set VSVER=14
  call "%VS140COMNTOOLS%..\..\VC\vcvarsall.bat" amd64
) else if "%VS120COMNTOOLS%" neq "" (
  echo Found Visual Studio 2013 installation
  set VSVER=12
  call "%VS120COMNTOOLS%..\..\VC\vcvarsall.bat" amd64
)

rem Ensure that we have a path to devenv.exe
if "%VSVER%" == "" goto novs

rem Create output directories if necessary
if not exist "..\lib" mkdir "..\lib"
if not exist "..\obj" mkdir "..\obj"

rem ===========================================================================
echo **************************************************************************
if exist "./ALGLIB" (
  call ./vc_internal/make_alglib.cmd vc%VSVER% dbg x64
  call ./vc_internal/make_alglib.cmd vc%VSVER% opt x64
) else (
  echo ALGLIB not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./SuiteSparse" (
  call ./vc_internal/make_umfpack.cmd vc%VSVER% dbg x64
  call ./vc_internal/make_umfpack.cmd vc%VSVER% opt x64
) else (
  echo SuiteSparse not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./parmetis" (
  call ./vc_internal/make_parmetis.cmd vc%VSVER% dbg x64
  call ./vc_internal/make_parmetis.cmd vc%VSVER% opt x64
) else (
  echo ParMETIS not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./fparser" (
  call ./vc_internal/make_fparser.cmd vc%VSVER% dbg x64
  call ./vc_internal/make_fparser.cmd vc%VSVER% opt x64
) else (
  echo fparser not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./zlib" (
  call ./vc_internal/make_zlib.cmd vc%VSVER% dbg x64
  call ./vc_internal/make_zlib.cmd vc%VSVER% opt x64
) else (
  echo zlib not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./zfp" (
  call ./vc_internal/make_zfp.cmd vc%VSVER% dbg x64
  call ./vc_internal/make_zfp.cmd vc%VSVER% opt x64
) else (
  echo zfp not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./triangle" (
  call ./vc_internal/make_triangle.cmd vc%VSVER% dbg x64
  call ./vc_internal/make_triangle.cmd vc%VSVER% opt x64
) else (
  echo triangle not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./hypre" (
  call ./vc_internal/make_hypre.cmd vc%VSVER% dbg x64 serial
  call ./vc_internal/make_hypre.cmd vc%VSVER% opt x64 serial
  if not "%MSMPI_INC%" == "" (
    call ./vc_internal/make_hypre.cmd vc%VSVER% dbg x64 mpi
    call ./vc_internal/make_hypre.cmd vc%VSVER% opt x64 mpi
  ) else (
    echo MPI not found; skipping MPI-based hypre library...
  )
) else (
  echo hypre not found; skipping...
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