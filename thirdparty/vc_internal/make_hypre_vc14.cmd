@echo off

if "%3" == "" goto errcall

rem Check build mode
if "%1" == "dbg" set CXXFLAGS=/FS /Od /RTC1 /MDd
if "%1" == "opt" set CXXFLAGS=/MP /Gy /Gm- /O2 /Ob2 /Oi /MD
if "%2" == "x64" set LIBFLAGS=/MACHINE:X64
if "%2" == "x86" set LIBFLAGS=/MACHINE:X86
if "%3" == "mpi" set CXXFLAGS=%CXXFLAGS% /DHYPRE_HAVE_MPI /I"%MSMPI_INC% "
if "%3" == "serial" set CXXFLAGS=%CXXFLAGS% /DHYPRE_SEQUENTIAL
goto build

:errcall
echo.
echo ERROR: Do not execute this script directly.
echo        Execute 'make_win32.objmd' or 'make_win64.objmd' instead
echo.
goto end

rem ===============================================================================================
:build
set OBJPATH=..\obj\hypre.vc14-%1-%3-%2

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4028" /wd"4244" /wd"4293" /wd"4305" /wd"4996"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/hypre.vc14-%1-%3-%2/hypre.vc14-%1-%2.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/hypre.vc14-%1-%3-%2/hypre.vc14-%1-%2.pch"
set CXXFLAGS=%CXXFLAGS% /I"./hypre/src"

echo.
echo Compiling hypre utilities Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/utilities/amg_linklist.c          /Fo"%OBJPATH%/amg_linklist.obj
cl %CXXFLAGS1% ./hypre/src/utilities/binsearch.c             /Fo"%OBJPATH%/binsearch.obj
cl %CXXFLAGS1% ./hypre/src/utilities/exchange_data.c         /Fo"%OBJPATH%/exchange_data.obj
cl %CXXFLAGS1% ./hypre/src/utilities/fortran_matrix.c        /Fo"%OBJPATH%/fortran_matrix.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_ap.c              /Fo"%OBJPATH%/hypre_ap.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_complex.c         /Fo"%OBJPATH%/hypre_complex.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_memory.c          /Fo"%OBJPATH%/hypre_memory.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_printf.c          /Fo"%OBJPATH%/hypre_printf.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_qsort.c           /Fo"%OBJPATH%/hypre_qsort.obj
cl %CXXFLAGS1% ./hypre/src/utilities/memory_dmalloc.c        /Fo"%OBJPATH%/memory_dmalloc.obj
cl %CXXFLAGS1% ./hypre/src/utilities/mpistubs.c              /Fo"%OBJPATH%/mpistubs.obj
cl %CXXFLAGS1% ./hypre/src/utilities/qsplit.c                /Fo"%OBJPATH%/qsplit.obj
cl %CXXFLAGS1% ./hypre/src/utilities/random.c                /Fo"%OBJPATH%/random.obj
cl %CXXFLAGS1% ./hypre/src/utilities/threading.c             /Fo"%OBJPATH%/threading.obj
cl %CXXFLAGS1% ./hypre/src/utilities/timer.c                 /Fo"%OBJPATH%/timer.obj
cl %CXXFLAGS1% ./hypre/src/utilities/timing.c                /Fo"%OBJPATH%/timing.obj
cl %CXXFLAGS1% ./hypre/src/utilities/umalloc_local.c         /Fo"%OBJPATH%/umalloc_local.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_mpi_comm_f2c.c    /Fo"%OBJPATH%/hypre_mpi_comm_f2c.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_error.c           /Fo"%OBJPATH%/hypre_error.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_prefix_sum.c      /Fo"%OBJPATH%/hypre_prefix_sum.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_merge_sort.c      /Fo"%OBJPATH%/hypre_merge_sort.obj
cl %CXXFLAGS1% ./hypre/src/utilities/hypre_hopscotch_hash.c  /Fo"%OBJPATH%/hypre_hopscotch_hash.obj

echo.
echo Compiling hypre blas Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/blas" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/blas/blas_utils.c  /Fo"%OBJPATH%/blas_utils.obj
cl %CXXFLAGS1% ./hypre/src/blas/dasum.c       /Fo"%OBJPATH%/dasum.obj
cl %CXXFLAGS1% ./hypre/src/blas/daxpy.c       /Fo"%OBJPATH%/daxpy.obj
cl %CXXFLAGS1% ./hypre/src/blas/dcopy.c       /Fo"%OBJPATH%/dcopy.obj
cl %CXXFLAGS1% ./hypre/src/blas/ddot.c        /Fo"%OBJPATH%/ddot.obj
cl %CXXFLAGS1% ./hypre/src/blas/dgemm.c       /Fo"%OBJPATH%/dgemm.obj
cl %CXXFLAGS1% ./hypre/src/blas/dgemv.c       /Fo"%OBJPATH%/dgemv.obj
cl %CXXFLAGS1% ./hypre/src/blas/dger.c        /Fo"%OBJPATH%/dger.obj
cl %CXXFLAGS1% ./hypre/src/blas/dnrm2.c       /Fo"%OBJPATH%/dnrm2.obj
cl %CXXFLAGS1% ./hypre/src/blas/drot.c        /Fo"%OBJPATH%/drot.obj
cl %CXXFLAGS1% ./hypre/src/blas/dscal.c       /Fo"%OBJPATH%/dscal.obj
cl %CXXFLAGS1% ./hypre/src/blas/dswap.c       /Fo"%OBJPATH%/dswap.obj
cl %CXXFLAGS1% ./hypre/src/blas/dsymm.c       /Fo"%OBJPATH%/dsymm.obj
cl %CXXFLAGS1% ./hypre/src/blas/dsymv.c       /Fo"%OBJPATH%/dsymv.obj
cl %CXXFLAGS1% ./hypre/src/blas/dsyr2.c       /Fo"%OBJPATH%/dsyr2.obj
cl %CXXFLAGS1% ./hypre/src/blas/dsyr2k.c      /Fo"%OBJPATH%/dsyr2k.obj
cl %CXXFLAGS1% ./hypre/src/blas/dsyrk.c       /Fo"%OBJPATH%/dsyrk.obj
cl %CXXFLAGS1% ./hypre/src/blas/dtrmm.c       /Fo"%OBJPATH%/dtrmm.obj
cl %CXXFLAGS1% ./hypre/src/blas/dtrmv.c       /Fo"%OBJPATH%/dtrmv.obj
cl %CXXFLAGS1% ./hypre/src/blas/dtrsm.c       /Fo"%OBJPATH%/dtrsm.obj
cl %CXXFLAGS1% ./hypre/src/blas/dtrsv.c       /Fo"%OBJPATH%/dtrsv.obj
cl %CXXFLAGS1% ./hypre/src/blas/idamax.c      /Fo"%OBJPATH%/idamax.obj

echo.
echo Compiling hypre lapack Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/blas" /I"./hypre/src/lapack" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/lapack/dbdsqr.c       /Fo"%OBJPATH%/dbdsqr.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgebd2.c       /Fo"%OBJPATH%/dgebd2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgebrd.c       /Fo"%OBJPATH%/dgebrd.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgelq2.c       /Fo"%OBJPATH%/dgelq2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgelqf.c       /Fo"%OBJPATH%/dgelqf.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgels.c        /Fo"%OBJPATH%/dgels.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgeqr2.c       /Fo"%OBJPATH%/dgeqr2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgeqrf.c       /Fo"%OBJPATH%/dgeqrf.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgesvd.c       /Fo"%OBJPATH%/dgesvd.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgetrf.c       /Fo"%OBJPATH%/dgetrf.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgetrs.c       /Fo"%OBJPATH%/dgetrs.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dgetf2.c       /Fo"%OBJPATH%/dgetf2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlabad.c       /Fo"%OBJPATH%/dlabad.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlabrd.c       /Fo"%OBJPATH%/dlabrd.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlacpy.c       /Fo"%OBJPATH%/dlacpy.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlae2.c        /Fo"%OBJPATH%/dlae2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlaev2.c       /Fo"%OBJPATH%/dlaev2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlange.c       /Fo"%OBJPATH%/dlange.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlanst.c       /Fo"%OBJPATH%/dlanst.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlansy.c       /Fo"%OBJPATH%/dlansy.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlapy2.c       /Fo"%OBJPATH%/dlapy2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlarfb.c       /Fo"%OBJPATH%/dlarfb.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlarf.c        /Fo"%OBJPATH%/dlarf.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlarfg.c       /Fo"%OBJPATH%/dlarfg.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlarft.c       /Fo"%OBJPATH%/dlarft.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlartg.c       /Fo"%OBJPATH%/dlartg.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlas2.c        /Fo"%OBJPATH%/dlas2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlascl.c       /Fo"%OBJPATH%/dlascl.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlaset.c       /Fo"%OBJPATH%/dlaset.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasq1.c       /Fo"%OBJPATH%/dlasq1.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasq2.c       /Fo"%OBJPATH%/dlasq2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasq3.c       /Fo"%OBJPATH%/dlasq3.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasq4.c       /Fo"%OBJPATH%/dlasq4.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasq5.c       /Fo"%OBJPATH%/dlasq5.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasq6.c       /Fo"%OBJPATH%/dlasq6.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasr.c        /Fo"%OBJPATH%/dlasr.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasrt.c       /Fo"%OBJPATH%/dlasrt.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlassq.c       /Fo"%OBJPATH%/dlassq.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlaswp.c       /Fo"%OBJPATH%/dlaswp.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlasv2.c       /Fo"%OBJPATH%/dlasv2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlatrd.c       /Fo"%OBJPATH%/dlatrd.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorg2l.c       /Fo"%OBJPATH%/dorg2l.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorg2r.c       /Fo"%OBJPATH%/dorg2r.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorgbr.c       /Fo"%OBJPATH%/dorgbr.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorgl2.c       /Fo"%OBJPATH%/dorgl2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorglq.c       /Fo"%OBJPATH%/dorglq.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorgql.c       /Fo"%OBJPATH%/dorgql.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorgqr.c       /Fo"%OBJPATH%/dorgqr.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorgtr.c       /Fo"%OBJPATH%/dorgtr.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorm2r.c       /Fo"%OBJPATH%/dorm2r.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dormbr.c       /Fo"%OBJPATH%/dormbr.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dorml2.c       /Fo"%OBJPATH%/dorml2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dormlq.c       /Fo"%OBJPATH%/dormlq.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dormqr.c       /Fo"%OBJPATH%/dormqr.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dpotf2.c       /Fo"%OBJPATH%/dpotf2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dpotrf.c       /Fo"%OBJPATH%/dpotrf.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dpotrs.c       /Fo"%OBJPATH%/dpotrs.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dsteqr.c       /Fo"%OBJPATH%/dsteqr.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dsterf.c       /Fo"%OBJPATH%/dsterf.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dsyev.c        /Fo"%OBJPATH%/dsyev.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dsygs2.c       /Fo"%OBJPATH%/dsygs2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dsygst.c       /Fo"%OBJPATH%/dsygst.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dsygv.c        /Fo"%OBJPATH%/dsygv.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dsytd2.c       /Fo"%OBJPATH%/dsytd2.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dsytrd.c       /Fo"%OBJPATH%/dsytrd.obj
cl %CXXFLAGS1% ./hypre/src/lapack/ieeeck.c       /Fo"%OBJPATH%/ieeeck.obj
cl %CXXFLAGS1% ./hypre/src/lapack/ilaenv.c       /Fo"%OBJPATH%/ilaenv.obj
cl %CXXFLAGS1% ./hypre/src/lapack/lapack_utils.c /Fo"%OBJPATH%/lapack_utils.obj
cl %CXXFLAGS1% ./hypre/src/lapack/lsame.c        /Fo"%OBJPATH%/lsame.obj
cl %CXXFLAGS1% ./hypre/src/lapack/xerbla.c       /Fo"%OBJPATH%/xerbla.obj
cl %CXXFLAGS1% ./hypre/src/lapack/dlamch.c       /Fo"%OBJPATH%/dlamch.obj

echo.
echo Compiling hypre multivector Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/multivector" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/multivector/multivector.c       /Fo"%OBJPATH%/multivector.obj
cl %CXXFLAGS1% ./hypre/src/multivector/temp_multivector.c  /Fo"%OBJPATH%/temp_multivector.obj

echo.
echo Compiling hypre krylov Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/krylov" /I"./hypre/src/multivector" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/krylov/bicgstab.c        /Fo"%OBJPATH%/bicgstab.obj
cl %CXXFLAGS1% ./hypre/src/krylov/cgnr.c            /Fo"%OBJPATH%/cgnr.obj
cl %CXXFLAGS1% ./hypre/src/krylov/gmres.c           /Fo"%OBJPATH%/gmres.obj
cl %CXXFLAGS1% ./hypre/src/krylov/flexgmres.c       /Fo"%OBJPATH%/flexgmres.obj
cl %CXXFLAGS1% ./hypre/src/krylov/lgmres.c          /Fo"%OBJPATH%/lgmres.obj
cl %CXXFLAGS1% ./hypre/src/krylov/HYPRE_bicgstab.c  /Fo"%OBJPATH%/HYPRE_bicgstab.obj
cl %CXXFLAGS1% ./hypre/src/krylov/HYPRE_cgnr.c      /Fo"%OBJPATH%/HYPRE_cgnr.obj
cl %CXXFLAGS1% ./hypre/src/krylov/HYPRE_gmres.c     /Fo"%OBJPATH%/HYPRE_gmres.obj
cl %CXXFLAGS1% ./hypre/src/krylov/HYPRE_lgmres.c    /Fo"%OBJPATH%/HYPRE_lgmres.obj
cl %CXXFLAGS1% ./hypre/src/krylov/HYPRE_flexgmres.c /Fo"%OBJPATH%/HYPRE_flexgmres.obj
cl %CXXFLAGS1% ./hypre/src/krylov/HYPRE_pcg.c       /Fo"%OBJPATH%/HYPRE_pcg.obj
cl %CXXFLAGS1% ./hypre/src/krylov/pcg.c             /Fo"%OBJPATH%/pcg.obj
cl %CXXFLAGS1% ./hypre/src/krylov/HYPRE_lobpcg.c    /Fo"%OBJPATH%/HYPRE_lobpcg.obj
cl %CXXFLAGS1% ./hypre/src/krylov/lobpcg.c          /Fo"%OBJPATH%/lobpcg.obj

echo.
echo Compiling hypre seq_mv Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/seq_mv" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/seq_mv/csr_matop.c                /Fo"%OBJPATH%/csr_matop.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/csr_matrix.c               /Fo"%OBJPATH%/csr_matrix.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/csr_matvec.c               /Fo"%OBJPATH%/csr_matvec.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/genpart.c                  /Fo"%OBJPATH%/genpart.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/HYPRE_csr_matrix.c         /Fo"%OBJPATH%/HYPRE_csr_matrix.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/HYPRE_mapped_matrix.c      /Fo"%OBJPATH%/HYPRE_mapped_matrix.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/HYPRE_multiblock_matrix.c  /Fo"%OBJPATH%/HYPRE_multiblock_matrix.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/HYPRE_vector.c             /Fo"%OBJPATH%/HYPRE_vector.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/mapped_matrix.c            /Fo"%OBJPATH%/mapped_matrix.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/multiblock_matrix.c        /Fo"%OBJPATH%/multiblock_matrix.obj
cl %CXXFLAGS1% ./hypre/src/seq_mv/vector.c                   /Fo"%OBJPATH%/vector.obj

echo.
echo Compiling hypre parcsr_mv Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/parcsr_mv" /I"./hypre/src/seq_mv" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/communicationT.c           /Fo"%OBJPATH%/communicationT.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/HYPRE_parcsr_matrix.c      /Fo"%OBJPATH%/HYPRE_parcsr_matrix.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/HYPRE_parcsr_vector.c      /Fo"%OBJPATH%/HYPRE_parcsr_vector.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/new_commpkg.c              /Fo"%OBJPATH%/new_commpkg.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/numbers.c                  /Fo"%OBJPATH%/numbers.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_aat.c              /Fo"%OBJPATH%/par_csr_aat.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_assumed_part.c     /Fo"%OBJPATH%/par_csr_assumed_part.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_bool_matop.c       /Fo"%OBJPATH%/par_csr_bool_matop.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_bool_matrix.c      /Fo"%OBJPATH%/par_csr_bool_matrix.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_communication.c    /Fo"%OBJPATH%/par_csr_communication.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_matop.c            /Fo"%OBJPATH%/par_csr_matop.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_matrix.c           /Fo"%OBJPATH%/par_csr_matrix.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_matop_marked.c     /Fo"%OBJPATH%/par_csr_matop_marked.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_csr_matvec.c           /Fo"%OBJPATH%/par_csr_matvec.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_vector.c               /Fo"%OBJPATH%/par_vector.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_mv/par_make_system.c          /Fo"%OBJPATH%/par_make_system.obj

echo.
echo Compiling hypre parcsr_block_mv Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/krylov"  /I"./hypre/src/IJ_mv" /I"./hypre/src/parcsr_ls" /I"./hypre/src/parcsr_mv" /I"./hypre/src/seq_mv" /I"./hypre/src/multivector" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/csr_block_matrix.c                  /Fo"%OBJPATH%/csr_block_matrix.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/csr_block_matvec.c                  /Fo"%OBJPATH%/csr_block_matvec.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/par_csr_block_matrix.c              /Fo"%OBJPATH%/par_csr_block_matrix.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/par_csr_block_matvec.c              /Fo"%OBJPATH%/par_csr_block_matvec.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/par_csr_block_comm.c                /Fo"%OBJPATH%/par_csr_block_comm.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/par_csr_block_rap.c                 /Fo"%OBJPATH%/par_csr_block_rap.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/par_csr_block_rap_communication.c   /Fo"%OBJPATH%/par_csr_block_rap_communication.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/par_csr_block_interp.c              /Fo"%OBJPATH%/par_csr_block_interp.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/par_csr_block_relax.c               /Fo"%OBJPATH%/par_csr_block_relax.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_block_mv/par_block_nodal_systems.c           /Fo"%OBJPATH%/par_block_nodal_systems.obj

echo.
echo Compiling hypre distributed_matrix Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/distributed_matrix" /I"./hypre/src/parcsr_mv" /I"./hypre/src/seq_mv" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/distributed_matrix/distributed_matrix.c          /Fo"%OBJPATH%/distributed_matrix.obj
cl %CXXFLAGS1% ./hypre/src/distributed_matrix/HYPRE_distributed_matrix.c    /Fo"%OBJPATH%/HYPRE_distributed_matrix.obj
cl %CXXFLAGS1% ./hypre/src/distributed_matrix/distributed_matrix_ISIS.c     /Fo"%OBJPATH%/distributed_matrix_ISIS.obj
cl %CXXFLAGS1% ./hypre/src/distributed_matrix/distributed_matrix_PETSc.c    /Fo"%OBJPATH%/distributed_matrix_PETSc.obj
cl %CXXFLAGS1% ./hypre/src/distributed_matrix/distributed_matrix_parcsr.c   /Fo"%OBJPATH%/distributed_matrix_parcsr.obj

echo.
echo Compiling hypre IJ_mv Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/IJ_mv" /I"./hypre/src/parcsr_ls" /I"./hypre/src/struct_mv" /I"./hypre/src/parcsr_mv" /I"./hypre/src/seq_mv" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/IJ_mv/aux_parcsr_matrix.c    /Fo"%OBJPATH%/aux_parcsr_matrix.obj
cl %CXXFLAGS1% ./hypre/src/IJ_mv/aux_par_vector.c       /Fo"%OBJPATH%/aux_par_vector.obj
cl %CXXFLAGS1% ./hypre/src/IJ_mv/HYPRE_IJMatrix.c       /Fo"%OBJPATH%/HYPRE_IJMatrix.obj
cl %CXXFLAGS1% ./hypre/src/IJ_mv/HYPRE_IJVector.c       /Fo"%OBJPATH%/HYPRE_IJVector.obj
cl %CXXFLAGS1% ./hypre/src/IJ_mv/IJ_assumed_part.c      /Fo"%OBJPATH%/IJ_assumed_part.obj
cl %CXXFLAGS1% ./hypre/src/IJ_mv/IJMatrix.c             /Fo"%OBJPATH%/IJMatrix.obj
cl %CXXFLAGS1% ./hypre/src/IJ_mv/IJMatrix_parcsr.c      /Fo"%OBJPATH%/IJMatrix_parcsr.obj
cl %CXXFLAGS1% ./hypre/src/IJ_mv/IJVector.c             /Fo"%OBJPATH%/IJVector.obj
cl %CXXFLAGS1% ./hypre/src/IJ_mv/IJVector_parcsr.c      /Fo"%OBJPATH%/IJVector_parcsr.obj

echo.
echo Compiling hypre matrix_matrix Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/matrix_matrix" /I"./hypre/src/distributed_matrix" /I"./hypre/src/parcsr_mv" /I"./hypre/src/seq_mv" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/matrix_matrix/HYPRE_ConvertParCSRMatrixToDistributedMatrix.c  /Fo"%OBJPATH%/HYPRE_ConvertParCSRMatrixToDistributedMatrix.obj
cl %CXXFLAGS1% ./hypre/src/matrix_matrix/HYPRE_ConvertPETScMatrixToDistributedMatrix.c   /Fo"%OBJPATH%/HYPRE_ConvertPETScMatrixToDistributedMatrix.obj

echo.
echo Compiling hypre distributed_ls/pilut Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/distributed_ls/pilut" /I"./hypre/src/distributed_matrix" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/comm.c                                  /Fo"%OBJPATH%/comm.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/debug.c                                 /Fo"%OBJPATH%/debug.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/distributed_qsort.c                     /Fo"%OBJPATH%/distributed_qsort.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/distributed_qsort_si.c                  /Fo"%OBJPATH%/distributed_qsort_si.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/HYPRE_DistributedMatrixPilutSolver.c    /Fo"%OBJPATH%/HYPRE_DistributedMatrixPilutSolver.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/ilut.c                                  /Fo"%OBJPATH%/ilut.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/parilut.c                               /Fo"%OBJPATH%/parilut.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/parutil.c                               /Fo"%OBJPATH%/parutil.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/pblas1.c                                /Fo"%OBJPATH%/pblas1.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/serilut.c                               /Fo"%OBJPATH%/serilut.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/trifactor.c                             /Fo"%OBJPATH%/trifactor.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/pilut/util.c                                  /Fo"%OBJPATH%/util.obj

echo.
echo Compiling hypre distributed_ls/ParaSails Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/distributed_ls/ParaSails" /I"./hypre/src/distributed_matrix" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/ConjGrad.c         /Fo"%OBJPATH%/ConjGrad.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/DiagScale.c        /Fo"%OBJPATH%/DiagScale.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/FGmres.c           /Fo"%OBJPATH%/FGmres.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/Hash.c             /Fo"%OBJPATH%/Hash.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/hypre_ParaSails.c  /Fo"%OBJPATH%/hypre_ParaSails.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/LoadBal.c          /Fo"%OBJPATH%/LoadBal.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/Matrix.c           /Fo"%OBJPATH%/Matrix.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/Mem.c              /Fo"%OBJPATH%/Mem.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/Numbering.c        /Fo"%OBJPATH%/Numbering.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/OrderStat.c        /Fo"%OBJPATH%/OrderStat.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/ParaSails.c        /Fo"%OBJPATH%/ParaSails.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/PrunedRows.c       /Fo"%OBJPATH%/PrunedRows.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/RowPatt.c          /Fo"%OBJPATH%/RowPatt.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/ParaSails/StoredRows.c       /Fo"%OBJPATH%/StoredRows.obj

echo.
echo Compiling hypre distributed_ls/Euclid Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/distributed_ls/Euclid" /I"./hypre/src/parcsr_mv" /I"./hypre/src/seq_mv" /I"./hypre/src/distributed_matrix" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/blas_dh.c            /Fo"%OBJPATH%/blas_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Euclid_apply.c       /Fo"%OBJPATH%/Euclid_apply.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Euclid_dh.c          /Fo"%OBJPATH%/Euclid_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/ExternalRows_dh.c    /Fo"%OBJPATH%/ExternalRows_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Factor_dh.c          /Fo"%OBJPATH%/Factor_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/getRow_dh.c          /Fo"%OBJPATH%/getRow_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/globalObjects.c      /Fo"%OBJPATH%/globalObjects.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Hash_dh.c            /Fo"%OBJPATH%/Hash_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Hash_i_dh.c          /Fo"%OBJPATH%/Hash_i_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/ilu_mpi_bj.c         /Fo"%OBJPATH%/ilu_mpi_bj.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/ilu_mpi_pilu.c       /Fo"%OBJPATH%/ilu_mpi_pilu.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/ilu_seq.c            /Fo"%OBJPATH%/ilu_seq.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/io_dh.c              /Fo"%OBJPATH%/io_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/krylov_dh.c          /Fo"%OBJPATH%/krylov_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Mat_dh.c             /Fo"%OBJPATH%/Mat_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/mat_dh_private.c     /Fo"%OBJPATH%/mat_dh_private.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/MatGenFD.c           /Fo"%OBJPATH%/MatGenFD.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Mem_dh.c             /Fo"%OBJPATH%/Mem_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Numbering_dh.c       /Fo"%OBJPATH%/Numbering_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Parser_dh.c          /Fo"%OBJPATH%/Parser_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/shellSort_dh.c       /Fo"%OBJPATH%/shellSort_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/sig_dh.c             /Fo"%OBJPATH%/sig_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/SortedList_dh.c      /Fo"%OBJPATH%/SortedList_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/SortedSet_dh.c       /Fo"%OBJPATH%/SortedSet_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/SubdomainGraph_dh.c  /Fo"%OBJPATH%/SubdomainGraph_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/TimeLog_dh.c         /Fo"%OBJPATH%/TimeLog_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Timer_dh.c           /Fo"%OBJPATH%/Timer_dh.obj
cl %CXXFLAGS1% ./hypre/src/distributed_ls/Euclid/Vec_dh.c             /Fo"%OBJPATH%/Vec_dh.obj

echo.
echo Compiling hypre parcsr_ls Sources...
set CXXFLAGS1=%CXXFLAGS% /I"./hypre/src/parcsr_ls" /I"./hypre/src/multivector" /I"./hypre/src/krylov" /I"./hypre/src/distributed_matrix" /I"./hypre/src/matrix_matrix" /I"./hypre/src/IJ_mv" /I"./hypre/src/parcsr_block_mv" /I"./hypre/src/parcsr_mv" /I"./hypre/src/seq_mv" /I"./hypre/src/utilities"
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/amg_hybrid.c                    /Fo"%OBJPATH%/amg_hybrid.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/aux_interp.c                    /Fo"%OBJPATH%/aux_interp.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/gen_redcs_mat.c                 /Fo"%OBJPATH%/gen_redcs_mat.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_amg.c              /Fo"%OBJPATH%/HYPRE_parcsr_amg.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_bicgstab.c         /Fo"%OBJPATH%/HYPRE_parcsr_bicgstab.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_block.c            /Fo"%OBJPATH%/HYPRE_parcsr_block.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_cgnr.c             /Fo"%OBJPATH%/HYPRE_parcsr_cgnr.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_Euclid.c           /Fo"%OBJPATH%/HYPRE_parcsr_Euclid.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_gmres.c            /Fo"%OBJPATH%/HYPRE_parcsr_gmres.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_flexgmres.c        /Fo"%OBJPATH%/HYPRE_parcsr_flexgmres.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_lgmres.c           /Fo"%OBJPATH%/HYPRE_parcsr_lgmres.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_hybrid.c           /Fo"%OBJPATH%/HYPRE_parcsr_hybrid.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_int.c              /Fo"%OBJPATH%/HYPRE_parcsr_int.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_ParaSails.c        /Fo"%OBJPATH%/HYPRE_parcsr_ParaSails.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_pcg.c              /Fo"%OBJPATH%/HYPRE_parcsr_pcg.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_pilut.c            /Fo"%OBJPATH%/HYPRE_parcsr_pilut.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_parcsr_schwarz.c          /Fo"%OBJPATH%/HYPRE_parcsr_schwarz.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_ams.c                     /Fo"%OBJPATH%/HYPRE_ams.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_ads.c                     /Fo"%OBJPATH%/HYPRE_ads.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/HYPRE_ame.c                     /Fo"%OBJPATH%/HYPRE_ame.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_amg.c                       /Fo"%OBJPATH%/par_amg.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_amg_setup.c                 /Fo"%OBJPATH%/par_amg_setup.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_amg_solve.c                 /Fo"%OBJPATH%/par_amg_solve.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_amg_solveT.c                /Fo"%OBJPATH%/par_amg_solveT.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_cg_relax_wt.c               /Fo"%OBJPATH%/par_cg_relax_wt.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_coarsen.c                   /Fo"%OBJPATH%/par_coarsen.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_cgc_coarsen.c               /Fo"%OBJPATH%/par_cgc_coarsen.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_coarse_parms.c              /Fo"%OBJPATH%/par_coarse_parms.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_coordinates.c               /Fo"%OBJPATH%/par_coordinates.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_cr.c                        /Fo"%OBJPATH%/par_cr.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_cycle.c                     /Fo"%OBJPATH%/par_cycle.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_add_cycle.c                 /Fo"%OBJPATH%/par_add_cycle.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_difconv.c                   /Fo"%OBJPATH%/par_difconv.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_gsmg.c                      /Fo"%OBJPATH%/par_gsmg.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_indepset.c                  /Fo"%OBJPATH%/par_indepset.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_interp.c                    /Fo"%OBJPATH%/par_interp.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_jacobi_interp.c             /Fo"%OBJPATH%/par_jacobi_interp.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_multi_interp.c              /Fo"%OBJPATH%/par_multi_interp.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_laplace_27pt.c              /Fo"%OBJPATH%/par_laplace_27pt.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_laplace_9pt.c               /Fo"%OBJPATH%/par_laplace_9pt.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_laplace.c                   /Fo"%OBJPATH%/par_laplace.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_lr_interp.c                 /Fo"%OBJPATH%/par_lr_interp.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_nongalerkin.c               /Fo"%OBJPATH%/par_nongalerkin.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_nodal_systems.c             /Fo"%OBJPATH%/par_nodal_systems.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_rap.c                       /Fo"%OBJPATH%/par_rap.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_rap_communication.c         /Fo"%OBJPATH%/par_rap_communication.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_rotate_7pt.c                /Fo"%OBJPATH%/par_rotate_7pt.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_vardifconv.c                /Fo"%OBJPATH%/par_vardifconv.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_relax.c                     /Fo"%OBJPATH%/par_relax.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_relax_more.c                /Fo"%OBJPATH%/par_relax_more.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_relax_interface.c           /Fo"%OBJPATH%/par_relax_interface.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_scaled_matnorm.c            /Fo"%OBJPATH%/par_scaled_matnorm.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_schwarz.c                   /Fo"%OBJPATH%/par_schwarz.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_stats.c                     /Fo"%OBJPATH%/par_stats.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_strength.c                  /Fo"%OBJPATH%/par_strength.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_sv_interp.c                 /Fo"%OBJPATH%/par_sv_interp.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/par_sv_interp_ln.c              /Fo"%OBJPATH%/par_sv_interp_ln.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/partial.c                       /Fo"%OBJPATH%/partial.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/pcg_par.c                       /Fo"%OBJPATH%/pcg_par.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/schwarz.c                       /Fo"%OBJPATH%/schwarz.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/block_tridiag.c                 /Fo"%OBJPATH%/block_tridiag.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/ams.c                           /Fo"%OBJPATH%/ams.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/ads.c                           /Fo"%OBJPATH%/ads.obj
cl %CXXFLAGS1% ./hypre/src/parcsr_ls/ame.c                           /Fo"%OBJPATH%/ame.obj

rem ===============================================================================================
:link
set LIBFLAGS=%LIBFLAGS% /NOLOGO

echo.
echo Linking 'hypre.vc14-%1-%3-%2'
lib %LIBFLAGS% /OUT:"..\lib\hypre.vc14-%1-%3-%2.lib" ..\obj\hypre.vc14-%1-%3-%2\*.obj

rem ===============================================================================================
:end