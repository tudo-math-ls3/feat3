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
set OBJPATH=..\obj\alglib.vc14-%1-%2

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /TP /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4244"
set CXXFLAGS=%CXXFLAGS% /I"./ALGLIB/cpp/src"
set CXXFLAGS=%CXXFLAGS% /D "NBLAS" /D "NCHOLMOD" /D "_MBCS"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/alglib.vc14-%1-%2/alglib.vc14-%1-%2.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/alglib.vc14-%1-%2/alglib.vc14-%1-%2.pch"

echo Compiling ALGLIB Sources...
cl %CXXFLAGS% ./ALGLIB/cpp/src/alglibinternal.cpp   /Fo"%OBJPATH%/alglibinternal.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/alglibmisc.cpp       /Fo"%OBJPATH%/alglibmisc.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/ap.cpp               /Fo"%OBJPATH%/ap.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/dataanalysis.cpp     /Fo"%OBJPATH%/dataanalysis.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/diffequations.cpp    /Fo"%OBJPATH%/diffequations.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/fasttransforms.cpp   /Fo"%OBJPATH%/fasttransforms.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/integration.cpp      /Fo"%OBJPATH%/integration.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/interpolation.cpp    /Fo"%OBJPATH%/interpolation.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/linalg.cpp           /Fo"%OBJPATH%/linalg.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/optimization.cpp     /Fo"%OBJPATH%/optimization.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/solvers.cpp          /Fo"%OBJPATH%/solvers.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/specialfunctions.cpp /Fo"%OBJPATH%/specialfunctions.obj"
cl %CXXFLAGS% ./ALGLIB/cpp/src/statistics.cpp       /Fo"%OBJPATH%/statistics.obj"

rem ===============================================================================================
:link
set LIBFLAGS=%LIBFLAGS% /NOLOGO

echo.
echo Linking 'alglib.vc14-%1-%2'
lib %LIBFLAGS% /OUT:"..\lib\alglib.vc14-%1-%2.lib" ..\obj\alglib.vc14-%1-%2\*.obj
echo.

rem ===============================================================================================
:end