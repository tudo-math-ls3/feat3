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
set OBJPATH=..\obj\zlib.vc14-%1-%2

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /TC /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4996" /wd"4267"
set CXXFLAGS=%CXXFLAGS% /I"./zlib"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/zlib.vc14-%1-%2/zlib.vc14-%1-%2.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/zlib.vc14-%1-%2/zlib.vc14-%1-%2.pch"

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
set LIBFLAGS=%LIBFLAGS% /NOLOGO

echo.
echo Linking 'zlib.vc14-%1-%2'
lib %LIBFLAGS% /OUT:"..\lib\zlib.vc14-%1-%2.lib" ..\obj\zlib.vc14-%1-%2\*.obj
echo.

rem ===============================================================================================
:end