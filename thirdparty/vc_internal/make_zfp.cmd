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
set OBJPATH=..\obj\zfp.%1-%2-%3

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /TC /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4996" /wd"4267" /wd"4146"
set CXXFLAGS=%CXXFLAGS% /I"./zfp/include" /I"./zfp/src"
set CXXFLAGS=%CXXFLAGS% /DZFP_INT64="long long" /DZFP_INT64_SUFFIX="ll"
set CXXFLAGS=%CXXFLAGS% /DZFP_UINT64="unsigned long long" -DZFP_UINT64_SUFFIX="ull"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/zfp.%1-%2-%3/zfp.%1-%2-%3.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/zfp.%1-%2-%3/zfp.%1-%2-%3.pch"

echo Compiling zfp Sources...
cl %CXXFLAGS% ./zfp/src/bitstream.c /Fo"%OBJPATH%/bitstream.obj"
cl %CXXFLAGS% ./zfp/src/decode1d.c  /Fo"%OBJPATH%/decode1d.obj"
cl %CXXFLAGS% ./zfp/src/decode1f.c  /Fo"%OBJPATH%/decode1f.obj"
cl %CXXFLAGS% ./zfp/src/decode1i.c  /Fo"%OBJPATH%/decode1i.obj"
cl %CXXFLAGS% ./zfp/src/decode1l.c  /Fo"%OBJPATH%/decode1l.obj"
cl %CXXFLAGS% ./zfp/src/decode2d.c  /Fo"%OBJPATH%/decode2d.obj"
cl %CXXFLAGS% ./zfp/src/decode2f.c  /Fo"%OBJPATH%/decode2f.obj"
cl %CXXFLAGS% ./zfp/src/decode2i.c  /Fo"%OBJPATH%/decode2i.obj"
cl %CXXFLAGS% ./zfp/src/decode2l.c  /Fo"%OBJPATH%/decode2l.obj"
cl %CXXFLAGS% ./zfp/src/decode3d.c  /Fo"%OBJPATH%/decode3d.obj"
cl %CXXFLAGS% ./zfp/src/decode3f.c  /Fo"%OBJPATH%/decode3f.obj"
cl %CXXFLAGS% ./zfp/src/decode3i.c  /Fo"%OBJPATH%/decode3i.obj"
cl %CXXFLAGS% ./zfp/src/decode3l.c  /Fo"%OBJPATH%/decode3l.obj"
cl %CXXFLAGS% ./zfp/src/decode4d.c  /Fo"%OBJPATH%/decode4d.obj"
cl %CXXFLAGS% ./zfp/src/decode4f.c  /Fo"%OBJPATH%/decode4f.obj"
cl %CXXFLAGS% ./zfp/src/decode4i.c  /Fo"%OBJPATH%/decode4i.obj"
cl %CXXFLAGS% ./zfp/src/decode4l.c  /Fo"%OBJPATH%/decode4l.obj"
cl %CXXFLAGS% ./zfp/src/encode1d.c  /Fo"%OBJPATH%/encode1d.obj"
cl %CXXFLAGS% ./zfp/src/encode1f.c  /Fo"%OBJPATH%/encode1f.obj"
cl %CXXFLAGS% ./zfp/src/encode1i.c  /Fo"%OBJPATH%/encode1i.obj"
cl %CXXFLAGS% ./zfp/src/encode1l.c  /Fo"%OBJPATH%/encode1l.obj"
cl %CXXFLAGS% ./zfp/src/encode2d.c  /Fo"%OBJPATH%/encode2d.obj"
cl %CXXFLAGS% ./zfp/src/encode2f.c  /Fo"%OBJPATH%/encode2f.obj"
cl %CXXFLAGS% ./zfp/src/encode2i.c  /Fo"%OBJPATH%/encode2i.obj"
cl %CXXFLAGS% ./zfp/src/encode2l.c  /Fo"%OBJPATH%/encode2l.obj"
cl %CXXFLAGS% ./zfp/src/encode3d.c  /Fo"%OBJPATH%/encode3d.obj"
cl %CXXFLAGS% ./zfp/src/encode3f.c  /Fo"%OBJPATH%/encode3f.obj"
cl %CXXFLAGS% ./zfp/src/encode3i.c  /Fo"%OBJPATH%/encode3i.obj"
cl %CXXFLAGS% ./zfp/src/encode3l.c  /Fo"%OBJPATH%/encode3l.obj"
cl %CXXFLAGS% ./zfp/src/encode4d.c  /Fo"%OBJPATH%/encode4d.obj"
cl %CXXFLAGS% ./zfp/src/encode4f.c  /Fo"%OBJPATH%/encode4f.obj"
cl %CXXFLAGS% ./zfp/src/encode4i.c  /Fo"%OBJPATH%/encode4i.obj"
cl %CXXFLAGS% ./zfp/src/encode4l.c  /Fo"%OBJPATH%/encode4l.obj"
cl %CXXFLAGS% ./zfp/src/zfp.c       /Fo"%OBJPATH%/zfp.obj"

rem ===============================================================================================
:link

echo.
echo Linking 'zfp.%1-%2-%3'
lib %LIBFLAGS% /OUT:"..\lib\zfp.%1-%2-%3.lib" ..\obj\zfp.%1-%2-%3\*.obj
echo.

rem ===============================================================================================
:end