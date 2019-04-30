@echo off
@rem FEAT3: Finite Element Analysis Toolbox, Version 3
@rem Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
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
set OBJPATH=..\obj\zlib.%1-%2-%3

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /TC /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4996" /wd"4267"
set CXXFLAGS=%CXXFLAGS% /I"./zlib"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/zlib.%1-%2-%3/zlib.%1-%2-%3.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/zlib.%1-%2-%3/zlib.%1-%2-%3.pch"

echo Compiling zlib Sources...
cl %CXXFLAGS% ./zlib/adler32.c    /Fo"%OBJPATH%/adler32.obj"
cl %CXXFLAGS% ./zlib/compress.c   /Fo"%OBJPATH%/compress.obj"
cl %CXXFLAGS% ./zlib/crc32.c      /Fo"%OBJPATH%/crc32.obj"
cl %CXXFLAGS% ./zlib/deflate.c    /Fo"%OBJPATH%/deflate.obj"
cl %CXXFLAGS% ./zlib/gzclose.c    /Fo"%OBJPATH%/gzclose.obj"
cl %CXXFLAGS% ./zlib/gzlib.c      /Fo"%OBJPATH%/gzlib.obj"
cl %CXXFLAGS% ./zlib/gzread.c     /Fo"%OBJPATH%/gzread.obj"
cl %CXXFLAGS% ./zlib/gzwrite.c    /Fo"%OBJPATH%/gzwrite.obj"
cl %CXXFLAGS% ./zlib/infback.c    /Fo"%OBJPATH%/infback.obj"
cl %CXXFLAGS% ./zlib/inffast.c    /Fo"%OBJPATH%/inffast.obj"
cl %CXXFLAGS% ./zlib/inflate.c    /Fo"%OBJPATH%/inflate.obj"
cl %CXXFLAGS% ./zlib/inftrees.c   /Fo"%OBJPATH%/inftrees.obj"
cl %CXXFLAGS% ./zlib/trees.c      /Fo"%OBJPATH%/trees.obj"
cl %CXXFLAGS% ./zlib/uncompr.c    /Fo"%OBJPATH%/uncompr.obj"
cl %CXXFLAGS% ./zlib/zutil.c      /Fo"%OBJPATH%/zutil.obj"

rem ===============================================================================================
:link

echo.
echo Linking 'zlib.%1-%2-%3'
lib %LIBFLAGS% /OUT:"..\lib\zlib.%1-%2-%3.lib" ..\obj\zlib.%1-%2-%3\*.obj
echo.

rem ===============================================================================================
:end