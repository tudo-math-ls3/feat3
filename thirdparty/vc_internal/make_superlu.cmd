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
set OBJPATH=..\obj\superlu-dist.%1-%2-%3

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4996" /wd"4005" /wd"4013" /wd"4068" /wd"4101" /wd"4146" /wd"4244" /wd"4267" /wd"4305" /wd"4334" /wd"4477" /wd"4715"
set CXXFLAGS=%CXXFLAGS% /I"%MSMPI_INC% " /I"./SuperLU_DIST/CBLAS" /I"./SuperLU_DIST/SRC"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/superlu-dist.%1-%2-%3/superlu-dist.%1-%2-%3.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/superlu-dist.%1-%2-%3/superlu-dist.%1-%2-%3.pch"

echo.
echo Compiling SuperLU CBLAS Sources...
cl %CXXFLAGS% ./SuperLU_DIST/CBLAS/dgemm.c /Fo"%OBJPATH%/dgemm.obj"
cl %CXXFLAGS% ./SuperLU_DIST/CBLAS/dgemv.c /Fo"%OBJPATH%/dgemv.obj"
cl %CXXFLAGS% ./SuperLU_DIST/CBLAS/dger.c  /Fo"%OBJPATH%/dger.obj"
cl %CXXFLAGS% ./SuperLU_DIST/CBLAS/dtrsm.c /Fo"%OBJPATH%/dtrsm.obj"
cl %CXXFLAGS% ./SuperLU_DIST/CBLAS/dtrsv.c /Fo"%OBJPATH%/dtrsv.obj"
cl %CXXFLAGS% ./SuperLU_DIST/CBLAS/daxpy.c /Fo"%OBJPATH%/daxpy.obj"
cl %CXXFLAGS% ./SuperLU_DIST/CBLAS/dscal.c /Fo"%OBJPATH%/dscal.obj"
cl %CXXFLAGS% ./SuperLU_DIST/CBLAS/input_error_dist.c /Fo"%OBJPATH%/input_error_dist.obj"

echo.
echo Compiling SuperLU Sources...
cl %CXXFLAGS% ./SuperLU_DIST/SRC/TreeInterface.cpp          /Fo"%OBJPATH%/TreeInterface.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/sp_ienv.c                  /Fo"%OBJPATH%/sp_ienv.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/etree.c                    /Fo"%OBJPATH%/etree.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/sp_colorder.c              /Fo"%OBJPATH%/sp_colorder.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/get_perm_c.c               /Fo"%OBJPATH%/get_perm_c.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/colamd.c                   /Fo"%OBJPATH%/colamd.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/mmd.c                      /Fo"%OBJPATH%/mmd.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/comm.c                     /Fo"%OBJPATH%/comm.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/memory.c                   /Fo"%OBJPATH%/memory.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/util.c                     /Fo"%OBJPATH%/util.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/superlu_grid.c             /Fo"%OBJPATH%/superlu_grid.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pxerr_dist.c               /Fo"%OBJPATH%/pxerr_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/superlu_timer.c            /Fo"%OBJPATH%/superlu_timer.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/symbfact.c                 /Fo"%OBJPATH%/symbfact.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/psymbfact.c                /Fo"%OBJPATH%/psymbfact.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/psymbfact_util.c           /Fo"%OBJPATH%/psymbfact_util.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/get_perm_c_parmetis.c      /Fo"%OBJPATH%/get_perm_c_parmetis.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/mc64ad_dist.c              /Fo"%OBJPATH%/mc64ad_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/xerr_dist.c                /Fo"%OBJPATH%/xerr_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/smach_dist.c               /Fo"%OBJPATH%/smach_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dmach_dist.c               /Fo"%OBJPATH%/dmach_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/superlu_dist_version.c     /Fo"%OBJPATH%/superlu_dist_version.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/superlu_grid3d.c           /Fo"%OBJPATH%/superlu_grid3d.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/supernodal_etree.c         /Fo"%OBJPATH%/supernodal_etree.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/supernodalForest.c         /Fo"%OBJPATH%/supernodalForest.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/trfAux.c                   /Fo"%OBJPATH%/trfAux.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/communication_aux.c        /Fo"%OBJPATH%/communication_aux.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/treeFactorization.c        /Fo"%OBJPATH%/treeFactorization.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/sec_structs.c              /Fo"%OBJPATH%/sec_structs.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dlangs_dist.c              /Fo"%OBJPATH%/dlangs_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dgsequ_dist.c              /Fo"%OBJPATH%/dgsequ_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dlaqgs_dist.c              /Fo"%OBJPATH%/dlaqgs_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dutil_dist.c               /Fo"%OBJPATH%/dutil_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dmemory_dist.c             /Fo"%OBJPATH%/dmemory_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dmyblas2_dist.c            /Fo"%OBJPATH%/dmyblas2_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dsp_blas2_dist.c           /Fo"%OBJPATH%/dsp_blas2_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dsp_blas3_dist.c           /Fo"%OBJPATH%/dsp_blas3_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgssvx.c                  /Fo"%OBJPATH%/pdgssvx.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgssvx_ABglobal.c         /Fo"%OBJPATH%/pdgssvx_ABglobal.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dreadhb.c                  /Fo"%OBJPATH%/dreadhb.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dreadrb.c                  /Fo"%OBJPATH%/dreadrb.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dreadtriple.c              /Fo"%OBJPATH%/dreadtriple.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dreadMM.c                  /Fo"%OBJPATH%/dreadMM.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dbinary_io.c               /Fo"%OBJPATH%/dbinary_io.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgsequ.c                  /Fo"%OBJPATH%/pdgsequ.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdlaqgs.c                  /Fo"%OBJPATH%/pdlaqgs.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dldperm_dist.c             /Fo"%OBJPATH%/dldperm_dist.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdlangs.c                  /Fo"%OBJPATH%/pdlangs.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdutil.c                   /Fo"%OBJPATH%/pdutil.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdsymbfact_distdata.c      /Fo"%OBJPATH%/pdsymbfact_distdata.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/ddistribute.c              /Fo"%OBJPATH%/ddistribute.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pddistribute.c             /Fo"%OBJPATH%/pddistribute.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgstrf.c                  /Fo"%OBJPATH%/pdgstrf.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dstatic_schedule.c         /Fo"%OBJPATH%/dstatic_schedule.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgstrf2.c                 /Fo"%OBJPATH%/pdgstrf2.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdGetDiagU.c               /Fo"%OBJPATH%/pdGetDiagU.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgstrs.c                  /Fo"%OBJPATH%/pdgstrs.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgstrs1.c                 /Fo"%OBJPATH%/pdgstrs1.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgstrs_lsum.c             /Fo"%OBJPATH%/pdgstrs_lsum.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgstrs_Bglobal.c          /Fo"%OBJPATH%/pdgstrs_Bglobal.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgsrfs.c                  /Fo"%OBJPATH%/pdgsrfs.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgsmv.c                   /Fo"%OBJPATH%/pdgsmv.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgsrfs_ABXglobal.c        /Fo"%OBJPATH%/pdgsrfs_ABXglobal.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgsmv_AXglobal.c          /Fo"%OBJPATH%/pdgsmv_AXglobal.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dreadtriple_noheader.c     /Fo"%OBJPATH%/dreadtriple_noheader.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dsuperlu_blas.c            /Fo"%OBJPATH%/dsuperlu_blas.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgssvx3d.c                /Fo"%OBJPATH%/pdgssvx3d.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pdgstrf3d.c                /Fo"%OBJPATH%/pdgstrf3d.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dtreeFactorization.c       /Fo"%OBJPATH%/dtreeFactorization.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dscatter3d.c               /Fo"%OBJPATH%/dscatter3d.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dgather.c                  /Fo"%OBJPATH%/dgather.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/pd3dcomm.c                 /Fo"%OBJPATH%/pd3dcomm.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dtrfAux.c                  /Fo"%OBJPATH%/dtrfAux.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dcommunication_aux.c       /Fo"%OBJPATH%/dcommunication_aux.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dtrfCommWrapper.c          /Fo"%OBJPATH%/dtrfCommWrapper.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dnrformat_loc3d.c          /Fo"%OBJPATH%/dnrformat_loc3d.obj"
cl %CXXFLAGS% ./SuperLU_DIST/SRC/dtreeFactorizationGPU.c    /Fo"%OBJPATH%/dtreeFactorizationGPU.obj"

rem ===============================================================================================
:link
set LIBFLAGS=%LIBFLAGS% /NOLOGO

echo.
echo Linking 'superlu-dist.%1-%2-%3'
lib %LIBFLAGS% /OUT:"..\lib\superlu-dist.%1-%2-%3.lib" ..\obj\superlu-dist.%1-%2-%3\*.obj
echo.

rem ===============================================================================================
:end