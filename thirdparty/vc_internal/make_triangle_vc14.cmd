@echo off

if "%2" == "" goto errcall

rem Check build mode
if "%1" == "dbg" set CXXFLAGS=/Gm /Od /RTC1 /MDd
if "%1" == "opt" set CXXFLAGS=/MP /Gy /Gm- /O2 /Ob2 /Oi /MD
if "%2" == "x64" set LIBFLAGS=/MACHINE:X64
if "%2" == "x86" set LIBFLAGS=/MACHINE:X86
goto build

:errcall
echo.
echo ERROR: Do not execute this script directly.
echo        Execute 'make_win32.cmd' or 'make_win64.cmd' instead
echo.
goto end

rem ===============================================================================================
:build
set OBJPATH=..\obj\triangle.vc14-%1-%2

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
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/triangle.vc14-%1-%2/triangle.vc14-%1-%2.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/triangle.vc14-%1-%2/triangle.vc14-%1-%2.pch"

echo Compiling zlib Sources...
cl %CXXFLAGS% ./triangle/triangle.c    /Fo"%OBJPATH%/triangle.obj"

rem ===============================================================================================
:link
set LIBFLAGS=%LIBFLAGS% /NOLOGO

echo.
echo Linking 'triangle.vc14-%1-%2'
lib %LIBFLAGS% /OUT:"..\lib\triangle.vc14-%1-%2.lib" ..\obj\triangle.vc14-%1-%2\*.obj
echo.

rem ===============================================================================================
:end