@echo off
@rem FEAT3: Finite Element Analysis Toolbox, Version 3
@rem Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
@rem FEAT3 is released under the GNU General Public License version 3,
@rem see the file 'copyright.txt' in the top level directory for details.

if "%3" == "" goto errcall

rem Check build mode
if "%2" == "dbg" set CXXFLAGS=/Od /RTC1 /MDd
if "%2" == "opt" set CXXFLAGS=/MP /Gy /O2 /Ob2 /Oi /MD
if "%3" == "x64" set LIBFLAGS=/MACHINE:X64 /NOLOGO
if "%3" == "x86" set LIBFLAGS=/MACHINE:X86 /NOLOGO
goto build

:errcall
echo.
echo ERROR: Do not execute this script directly.
echo        Execute 'make_win32.cmd' or 'make_win64.cmd' instead
echo.
goto end

rem ===============================================================================================
:build
set OBJPATH=..\obj\triangle.%1-%2-%3

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /TC /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4996" /wd"4267" /wd"4244" /wd"4311" /wd"4312"
set CXXFLAGS=%CXXFLAGS% /DANSI_DECLARATORS /DTRILIBRARY /DNO_TIMER /DCPU86
set CXXFLAGS=%CXXFLAGS% /I"./triangle"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/triangle.%1-%2-%3/triangle.%1-%2-%3.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/triangle.%1-%2-%3/triangle.%1-%2-%3.pch"

echo Compiling Triangle Sources...
cl %CXXFLAGS% ./triangle/triangle.c    /Fo"%OBJPATH%/triangle.obj"

rem ===============================================================================================
:link

echo.
echo Linking 'triangle.%1-%2-%3'
lib %LIBFLAGS% /OUT:"..\lib\triangle.%1-%2-%3.lib" ..\obj\triangle.%1-%2-%3\*.obj
echo.

rem ===============================================================================================
:end