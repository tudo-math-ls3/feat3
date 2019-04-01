@echo off
@rem FEAT3: Finite Element Analysis Toolbox, Version 3
@rem Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
@rem FEAT3 is released under the GNU General Public License version 3,
@rem see the file 'copyright.txt' in the top level directory for details.

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
call "%%VS%VSVER%0COMNTOOLS%%..\..\VC\vcvarsall.bat" x86

rem ===========================================================================
echo **************************************************************************
if exist "./ALGLIB" (
  call ./vc_internal/make_alglib_vc%VSVER%.cmd dbg x86
  call ./vc_internal/make_alglib_vc%VSVER%.cmd opt x86
) else (
  echo ALGLIB not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./SuiteSparse" (
  call ./vc_internal/make_umfpack_vc%VSVER%.cmd dbg x86
  call ./vc_internal/make_umfpack_vc%VSVER%.cmd opt x86
) else (
  echo SuiteSparse not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./parmetis" (
  call ./vc_internal/make_parmetis_vc%VSVER%.cmd dbg x86
  call ./vc_internal/make_parmetis_vc%VSVER%.cmd opt x86
) else (
  echo ParMETIS not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./fparser" (
  call ./vc_internal/make_fparser_vc%VSVER%.cmd dbg x86
  call ./vc_internal/make_fparser_vc%VSVER%.cmd opt x86
) else (
  echo fparser not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./zlib" (
  call ./vc_internal/make_zlib_vc%VSVER%.cmd dbg x86
  call ./vc_internal/make_zlib_vc%VSVER%.cmd opt x86
) else (
  echo zlib not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./triangle" (
  call ./vc_internal/make_triangle_vc%VSVER%.cmd dbg x86
  call ./vc_internal/make_triangle_vc%VSVER%.cmd opt x86
) else (
  echo triangle not found; skipping...
  echo.
)

rem ===========================================================================
echo **************************************************************************
if exist "./hypre" (
  call ./vc_internal/make_hypre_vc%VSVER%.cmd dbg x86 serial
  call ./vc_internal/make_hypre_vc%VSVER%.cmd opt x86 serial
  if not "%MSMPI_INC%" == "" (
    call ./vc_internal/make_hypre_vc%VSVER%.cmd dbg x86 mpi
    call ./vc_internal/make_hypre_vc%VSVER%.cmd opt x86 mpi
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