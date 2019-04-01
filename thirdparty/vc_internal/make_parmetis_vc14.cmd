@echo off
@rem FEAT3: Finite Element Analysis Toolbox, Version 3
@rem Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
@rem FEAT3 is released under the GNU General Public License version 3,
@rem see the file 'copyright.txt' in the top level directory for details.

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
set OBJPATH1=..\obj\gklib.vc14-%1-%2
set OBJPATH2=..\obj\metis.vc14-%1-%2
set OBJPATH3=..\obj\parmetis.vc14-%1-%2

if exist %OBJPATH1% (
  del %OBJPATH1%\*.obj
) else (
  mkdir %OBJPATH1%
)

if exist %OBJPATH2% (
  del %OBJPATH2%\*.obj
) else (
  mkdir %OBJPATH2%
)

if exist %OBJPATH3% (
  del %OBJPATH3%\*.obj
) else (
  mkdir %OBJPATH3%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4028" /wd"4018" /wd"4101" /wd"4102" /wd"4244" /wd"4267" /wd"4305" /wd"4996"
set CXXFLAGS=%CXXFLAGS% /TC /D "USE_GKREGEX"

rem Set GKLib compiler flags
set CXXFLAGS1=%CXXFLAGS% /I"./parmetis/metis/GKlib"
set CXXFLAGS1=%CXXFLAGS1% /Fd"../obj/gklib.vc14-%1-%2/gklib.vc14-%1-%2.pdb"
set CXXFLAGS1=%CXXFLAGS1% /Fp"../obj/gklib.vc14-%1-%2/gklib.vc14-%1-%2.pch"

rem Set METIS compiler flags
set CXXFLAGS2=%CXXFLAGS% /I"./parmetis/metis/include" /I"./parmetis/metis/libmetis" /I"./parmetis/metis/GKlib"
set CXXFLAGS2=%CXXFLAGS2% /Fd"../obj/metis.vc14-%1-%2/metis.vc14-%1-%2.pdb"
set CXXFLAGS2=%CXXFLAGS2% /Fp"../obj/metis.vc14-%1-%2/metis.vc14-%1-%2.pch"

rem Set ParMETIS compiler flags
set CXXFLAGS3=%CXXFLAGS%  /I"./parmetis/include" /I"./parmetis/libparmetis"
set CXXFLAGS3=%CXXFLAGS3% /I"./parmetis/metis/include" /I"./parmetis/metis/libmetis"
set CXXFLAGS3=%CXXFLAGS3% /I"./parmetis/metis/GKlib" /I"%MSMPI_INC% "
set CXXFLAGS3=%CXXFLAGS3% /Fd"../obj/parmetis.vc14-%1-%2/parmetis.vc14-%1-%2.pdb"
set CXXFLAGS3=%CXXFLAGS3% /Fp"../obj/parmetis.vc14-%1-%2/parmetis.vc14-%1-%2.pch"

echo.
echo Compiling GKlib Sources...
cl %CXXFLAGS1% ./parmetis/metis/GKlib/b64.c         /Fo"%OBJPATH1%/b64.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/blas.c        /Fo"%OBJPATH1%/blas.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/csr.c         /Fo"%OBJPATH1%/csr.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/error.c       /Fo"%OBJPATH1%/error.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/evaluate.c    /Fo"%OBJPATH1%/evaluate.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/fkvkselect.c  /Fo"%OBJPATH1%/fkvkselect.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/fs.c          /Fo"%OBJPATH1%/fs.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/getopt.c      /Fo"%OBJPATH1%/getopt.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/gkregex.c     /Fo"%OBJPATH1%/gkregex.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/graph.c       /Fo"%OBJPATH1%/graph.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/htable.c      /Fo"%OBJPATH1%/htable.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/io.c          /Fo"%OBJPATH1%/io.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/itemsets.c    /Fo"%OBJPATH1%/itemsets.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/mcore.c       /Fo"%OBJPATH1%/mcore.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/memory.c      /Fo"%OBJPATH1%/memory.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/omp.c         /Fo"%OBJPATH1%/omp.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/pdb.c         /Fo"%OBJPATH1%/pdb.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/pqueue.c      /Fo"%OBJPATH1%/pqueue.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/random.c      /Fo"%OBJPATH1%/random.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/rw.c          /Fo"%OBJPATH1%/rw.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/seq.c         /Fo"%OBJPATH1%/seq.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/sort.c        /Fo"%OBJPATH1%/sort.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/string.c      /Fo"%OBJPATH1%/string.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/timers.c      /Fo"%OBJPATH1%/timers.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/tokenizer.c   /Fo"%OBJPATH1%/tokenizer.obj
cl %CXXFLAGS1% ./parmetis/metis/GKlib/util.c        /Fo"%OBJPATH1%/util.obj

echo.
echo Compiling METIS Sources...
cl %CXXFLAGS2% ./parmetis/metis/libmetis/auxapi.c       /Fo"%OBJPATH2%/auxapi.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/balance.c      /Fo"%OBJPATH2%/balance.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/bucketsort.c   /Fo"%OBJPATH2%/bucketsort.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/checkgraph.c   /Fo"%OBJPATH2%/checkgraph.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/coarsen.c      /Fo"%OBJPATH2%/coarsen.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/compress.c     /Fo"%OBJPATH2%/compress.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/contig.c       /Fo"%OBJPATH2%/contig.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/debug.c        /Fo"%OBJPATH2%/debug.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/fm.c           /Fo"%OBJPATH2%/fm.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/fortran.c      /Fo"%OBJPATH2%/fortran.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/frename.c      /Fo"%OBJPATH2%/frename.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/gklib.c        /Fo"%OBJPATH2%/gklib.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/graph.c        /Fo"%OBJPATH2%/graph.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/initpart.c     /Fo"%OBJPATH2%/initpart.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/kmetis.c       /Fo"%OBJPATH2%/kmetis.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/kwayfm.c       /Fo"%OBJPATH2%/kwayfm.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/kwayrefine.c   /Fo"%OBJPATH2%/kwayrefine.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/mcutil.c       /Fo"%OBJPATH2%/mcutil.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/mesh.c         /Fo"%OBJPATH2%/mesh.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/meshpart.c     /Fo"%OBJPATH2%/meshpart.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/minconn.c      /Fo"%OBJPATH2%/minconn.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/mincover.c     /Fo"%OBJPATH2%/mincover.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/mmd.c          /Fo"%OBJPATH2%/mmd.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/ometis.c       /Fo"%OBJPATH2%/ometis.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/options.c      /Fo"%OBJPATH2%/options.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/parmetis.c     /Fo"%OBJPATH2%/parmetis.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/pmetis.c       /Fo"%OBJPATH2%/pmetis.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/refine.c       /Fo"%OBJPATH2%/refine.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/separator.c    /Fo"%OBJPATH2%/separator.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/sfm.c          /Fo"%OBJPATH2%/sfm.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/srefine.c      /Fo"%OBJPATH2%/srefine.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/stat.c         /Fo"%OBJPATH2%/stat.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/timing.c       /Fo"%OBJPATH2%/timing.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/util.c         /Fo"%OBJPATH2%/util.obj
cl %CXXFLAGS2% ./parmetis/metis/libmetis/wspace.c       /Fo"%OBJPATH2%/wspace.obj

echo.
echo Compiling ParMETIS Sources...
cl %CXXFLAGS3% ./parmetis/libparmetis/akwayfm.c         /Fo"%OBJPATH3%/akwayfm.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/ametis.c          /Fo"%OBJPATH3%/ametis.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/balancemylink.c   /Fo"%OBJPATH3%/balancemylink.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/comm.c            /Fo"%OBJPATH3%/comm.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/csrmatch.c        /Fo"%OBJPATH3%/csrmatch.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/ctrl.c            /Fo"%OBJPATH3%/ctrl.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/debug.c           /Fo"%OBJPATH3%/debug.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/diffutil.c        /Fo"%OBJPATH3%/diffutil.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/frename.c         /Fo"%OBJPATH3%/frename.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/gkmetis.c         /Fo"%OBJPATH3%/gkmetis.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/gkmpi.c           /Fo"%OBJPATH3%/gkmpi.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/graph.c           /Fo"%OBJPATH3%/graph.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/initbalance.c     /Fo"%OBJPATH3%/initbalance.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/initmsection.c    /Fo"%OBJPATH3%/initmsection.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/initpart.c        /Fo"%OBJPATH3%/initpart.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/kmetis.c          /Fo"%OBJPATH3%/kmetis.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/kwayrefine.c      /Fo"%OBJPATH3%/kwayrefine.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/match.c           /Fo"%OBJPATH3%/match.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/mdiffusion.c      /Fo"%OBJPATH3%/mdiffusion.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/mesh.c            /Fo"%OBJPATH3%/mesh.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/mmetis.c          /Fo"%OBJPATH3%/mmetis.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/move.c            /Fo"%OBJPATH3%/move.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/msetup.c          /Fo"%OBJPATH3%/msetup.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/node_refine.c     /Fo"%OBJPATH3%/node_refine.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/ometis.c          /Fo"%OBJPATH3%/ometis.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/pspases.c         /Fo"%OBJPATH3%/pspases.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/redomylink.c      /Fo"%OBJPATH3%/redomylink.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/remap.c           /Fo"%OBJPATH3%/remap.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/renumber.c        /Fo"%OBJPATH3%/renumber.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/rmetis.c          /Fo"%OBJPATH3%/rmetis.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/selectq.c         /Fo"%OBJPATH3%/selectq.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/serial.c          /Fo"%OBJPATH3%/serial.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/stat.c            /Fo"%OBJPATH3%/stat.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/timer.c           /Fo"%OBJPATH3%/timer.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/util.c            /Fo"%OBJPATH3%/util.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/wave.c            /Fo"%OBJPATH3%/wave.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/weird.c           /Fo"%OBJPATH3%/weird.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/wspace.c          /Fo"%OBJPATH3%/wspace.obj
cl %CXXFLAGS3% ./parmetis/libparmetis/xyzpart.c         /Fo"%OBJPATH3%/xyzpart.obj

rem ===============================================================================================
:link
set LIBFLAGS=%LIBFLAGS% /NOLOGO

echo.
echo Linking 'metis.vc14-%1-%2'
lib %LIBFLAGS% /OUT:"..\lib\metis.vc14-%1-%2.lib" ..\obj\gklib.vc14-%1-%2\*.obj ..\obj\metis.vc14-%1-%2\*.obj

echo.
echo Linking 'parmetis.vc14-%1-%2'
lib %LIBFLAGS% /OUT:"..\lib\parmetis.vc14-%1-%2.lib" ..\obj\parmetis.vc14-%1-%2\*.obj
echo.

rem ===============================================================================================
:end