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
echo        Execute 'make_umfpack_win32.cmd' or 'make_umfpack_win64.cmd' instead
echo.
goto end

rem ===============================================================================================
:build
set OBJPATH=..\obj\umfpack.vc12-%1-%2

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /TC /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo
set CXXFLAGS=%CXXFLAGS% /wd"4068" /wd"4101" /wd"4244" /wd"4267" /wd"4996"
set CXXFLAGS=%CXXFLAGS% /I"./UMFPACK/Include" /I"./AMD/Include" /I"./SuiteSparse_config"
set CXXFLAGS=%CXXFLAGS% /D "NBLAS" /D "NCHOLMOD" /D "_MBCS"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/umfpack.vc12-%1-%2/umfpack.vc12-%1-%2.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/umfpack.vc12-%1-%2/umfpack.vc12-%1-%2.pch"

echo Compiling SuiteSparse_config Sources...
cl.exe %CXXFLAGS% ./SuiteSparse_config/SuiteSparse_config.c /Fo"%OBJPATH%/SuiteSparse_config.obj"
cl.exe %CXXFLAGS% ./SuiteSparse_config/xerbla/xerbla.c /Fo"%OBJPATH%/xerbla.obj"

echo.
echo Compiling AMD Sources...
cl %CXXFLAGS% ./AMD/Source/amd_global.c /Fo"%OBJPATH%/amd_global.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_aat.c /Fo"%OBJPATH%/amd_i_aat.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_1.c /Fo"%OBJPATH%/amd_i_1.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_2.c /Fo"%OBJPATH%/amd_i_2.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_dump.c /Fo"%OBJPATH%/amd_i_dump.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_postorder.c /Fo"%OBJPATH%/amd_i_postorder.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_post_tree.c /Fo"%OBJPATH%/amd_i_post_tree.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_defaults.c /Fo"%OBJPATH%/amd_i_defaults.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_order.c /Fo"%OBJPATH%/amd_i_order.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_control.c /Fo"%OBJPATH%/amd_i_control.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_info.c /Fo"%OBJPATH%/amd_i_info.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_valid.c /Fo"%OBJPATH%/amd_i_valid.obj"
cl %CXXFLAGS% /D "DINT" ./AMD/Source/amd_preprocess.c /Fo"%OBJPATH%/amd_i_preprocess.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_aat.c /Fo"%OBJPATH%/amd_l_aat.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_1.c /Fo"%OBJPATH%/amd_l_1.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_2.c /Fo"%OBJPATH%/amd_l_2.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_dump.c /Fo"%OBJPATH%/amd_l_dump.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_postorder.c /Fo"%OBJPATH%/amd_l_postorder.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_post_tree.c /Fo"%OBJPATH%/amd_l_post_tree.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_defaults.c /Fo"%OBJPATH%/amd_l_defaults.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_order.c /Fo"%OBJPATH%/amd_l_order.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_control.c /Fo"%OBJPATH%/amd_l_control.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_info.c /Fo"%OBJPATH%/amd_l_info.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_valid.c /Fo"%OBJPATH%/amd_l_valid.obj"
cl %CXXFLAGS% /D "DLONG" ./AMD/Source/amd_preprocess.c /Fo"%OBJPATH%/amd_l_preprocess.obj"

echo.
echo Compiling UMFPACK Sources...
cl %CXXFLAGS% ./UMFPACK/Source/umfpack_global.c /Fo"%OBJPATH%/umfpack_global.obj"
cl %CXXFLAGS% ./UMFPACK/Source/umfpack_timer.c /Fo"%OBJPATH%/umfpack_timer.obj"
cl %CXXFLAGS% ./UMFPACK/Source/umfpack_tictoc.c /Fo"%OBJPATH%/umfpack_tictoc.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_analyze.c /Fo"%OBJPATH%/umf_i_analyze.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_apply_order.c /Fo"%OBJPATH%/umf_i_apply_order.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_colamd.c /Fo"%OBJPATH%/umf_i_colamd.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_cholmod.c /Fo"%OBJPATH%/umf_i_cholmod.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_free.c /Fo"%OBJPATH%/umf_i_free.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_fsize.c /Fo"%OBJPATH%/umf_i_fsize.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_is_permutation.c /Fo"%OBJPATH%/umf_i_is_permutation.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_malloc.c /Fo"%OBJPATH%/umf_i_malloc.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_realloc.c /Fo"%OBJPATH%/umf_i_realloc.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_report_perm.c /Fo"%OBJPATH%/umf_i_report_perm.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_singletons.c /Fo"%OBJPATH%/umf_i_singletons.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_analyze.c /Fo"%OBJPATH%/umf_l_analyze.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_apply_order.c /Fo"%OBJPATH%/umf_l_apply_order.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_colamd.c /Fo"%OBJPATH%/umf_l_colamd.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_cholmod.c /Fo"%OBJPATH%/umf_l_cholmod.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_free.c /Fo"%OBJPATH%/umf_l_free.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_fsize.c /Fo"%OBJPATH%/umf_l_fsize.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_is_permutation.c /Fo"%OBJPATH%/umf_l_is_permutation.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_malloc.c /Fo"%OBJPATH%/umf_l_malloc.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_realloc.c /Fo"%OBJPATH%/umf_l_realloc.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_report_perm.c /Fo"%OBJPATH%/umf_l_report_perm.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_singletons.c /Fo"%OBJPATH%/umf_l_singletons.obj"
cl %CXXFLAGS% /D "DINT" /D "CONJUGATE_SOLVE" ./UMFPACK/Source/umf_ltsolve.c /Fo"%OBJPATH%/umf_di_lhsolve.obj"
cl %CXXFLAGS% /D "DINT" /D "CONJUGATE_SOLVE" ./UMFPACK/Source/umf_utsolve.c /Fo"%OBJPATH%/umf_di_uhsolve.obj"
cl %CXXFLAGS% /D "DINT" /D "DO_MAP" ./UMFPACK/Source/umf_triplet.c /Fo"%OBJPATH%/umf_di_triplet_map_nox.obj"
cl %CXXFLAGS% /D "DINT" /D "DO_VALUES" ./UMFPACK/Source/umf_triplet.c /Fo"%OBJPATH%/umf_di_triplet_nomap_x.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_triplet.c /Fo"%OBJPATH%/umf_di_triplet_nomap_nox.obj"
cl %CXXFLAGS% /D "DINT" /D "DO_MAP" /D "DO_VALUES" ./UMFPACK/Source/umf_triplet.c /Fo"%OBJPATH%/umf_di_triplet_map_x.obj"
cl %CXXFLAGS% /D "DINT" /D "FIXQ" ./UMFPACK/Source/umf_assemble.c /Fo"%OBJPATH%/umf_di_assemble_fixq.obj"
cl %CXXFLAGS% /D "DINT" /D "DROP" ./UMFPACK/Source/umf_store_lu.c /Fo"%OBJPATH%/umf_di_store_lu_drop.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_assemble.c /Fo"%OBJPATH%/umf_di_assemble.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_blas3_update.c /Fo"%OBJPATH%/umf_di_blas3_update.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_build_tuples.c /Fo"%OBJPATH%/umf_di_build_tuples.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_create_element.c /Fo"%OBJPATH%/umf_di_create_element.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_dump.c /Fo"%OBJPATH%/umf_di_dump.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_extend_front.c /Fo"%OBJPATH%/umf_di_extend_front.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_garbage_collection.c /Fo"%OBJPATH%/umf_di_garbage_collection.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_get_memory.c /Fo"%OBJPATH%/umf_di_get_memory.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_init_front.c /Fo"%OBJPATH%/umf_di_init_front.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_kernel.c /Fo"%OBJPATH%/umf_di_kernel.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_kernel_init.c /Fo"%OBJPATH%/umf_di_kernel_init.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_kernel_wrapup.c /Fo"%OBJPATH%/umf_di_kernel_wrapup.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_local_search.c /Fo"%OBJPATH%/umf_di_local_search.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_lsolve.c /Fo"%OBJPATH%/umf_di_lsolve.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_ltsolve.c /Fo"%OBJPATH%/umf_di_ltsolve.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_mem_alloc_element.c /Fo"%OBJPATH%/umf_di_mem_alloc_element.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_mem_alloc_head_block.c /Fo"%OBJPATH%/umf_di_mem_alloc_head_block.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_mem_alloc_tail_block.c /Fo"%OBJPATH%/umf_di_mem_alloc_tail_block.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_mem_free_tail_block.c /Fo"%OBJPATH%/umf_di_mem_free_tail_block.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_mem_init_memoryspace.c /Fo"%OBJPATH%/umf_di_mem_init_memoryspace.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_report_vector.c /Fo"%OBJPATH%/umf_di_report_vector.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_row_search.c /Fo"%OBJPATH%/umf_di_row_search.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_scale_column.c /Fo"%OBJPATH%/umf_di_scale_column.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_set_stats.c /Fo"%OBJPATH%/umf_di_set_stats.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_solve.c /Fo"%OBJPATH%/umf_di_solve.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_symbolic_usage.c /Fo"%OBJPATH%/umf_di_symbolic_usage.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_transpose.c /Fo"%OBJPATH%/umf_di_transpose.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_tuple_lengths.c /Fo"%OBJPATH%/umf_di_tuple_lengths.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_usolve.c /Fo"%OBJPATH%/umf_di_usolve.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_utsolve.c /Fo"%OBJPATH%/umf_di_utsolve.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_valid_numeric.c /Fo"%OBJPATH%/umf_di_valid_numeric.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_valid_symbolic.c /Fo"%OBJPATH%/umf_di_valid_symbolic.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_grow_front.c /Fo"%OBJPATH%/umf_di_grow_front.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_start_front.c /Fo"%OBJPATH%/umf_di_start_front.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_store_lu.c /Fo"%OBJPATH%/umf_di_store_lu.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umf_scale.c /Fo"%OBJPATH%/umf_di_scale.obj"
cl %CXXFLAGS% /D "DINT" /D "WSOLVE" ./UMFPACK/Source/umfpack_solve.c /Fo"%OBJPATH%/umfpack_di_wsolve.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_col_to_triplet.c /Fo"%OBJPATH%/umfpack_di_col_to_triplet.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_defaults.c /Fo"%OBJPATH%/umfpack_di_defaults.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_free_numeric.c /Fo"%OBJPATH%/umfpack_di_free_numeric.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_free_symbolic.c /Fo"%OBJPATH%/umfpack_di_free_symbolic.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_get_numeric.c /Fo"%OBJPATH%/umfpack_di_get_numeric.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_get_lunz.c /Fo"%OBJPATH%/umfpack_di_get_lunz.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_get_symbolic.c /Fo"%OBJPATH%/umfpack_di_get_symbolic.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_get_determinant.c /Fo"%OBJPATH%/umfpack_di_get_determinant.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_numeric.c /Fo"%OBJPATH%/umfpack_di_numeric.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_qsymbolic.c /Fo"%OBJPATH%/umfpack_di_qsymbolic.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_control.c /Fo"%OBJPATH%/umfpack_di_report_control.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_info.c /Fo"%OBJPATH%/umfpack_di_report_info.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_matrix.c /Fo"%OBJPATH%/umfpack_di_report_matrix.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_numeric.c /Fo"%OBJPATH%/umfpack_di_report_numeric.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_perm.c /Fo"%OBJPATH%/umfpack_di_report_perm.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_status.c /Fo"%OBJPATH%/umfpack_di_report_status.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_symbolic.c /Fo"%OBJPATH%/umfpack_di_report_symbolic.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_triplet.c /Fo"%OBJPATH%/umfpack_di_report_triplet.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_report_vector.c /Fo"%OBJPATH%/umfpack_di_report_vector.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_solve.c /Fo"%OBJPATH%/umfpack_di_solve.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_symbolic.c /Fo"%OBJPATH%/umfpack_di_symbolic.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_transpose.c /Fo"%OBJPATH%/umfpack_di_transpose.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_triplet_to_col.c /Fo"%OBJPATH%/umfpack_di_triplet_to_col.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_scale.c /Fo"%OBJPATH%/umfpack_di_scale.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_load_numeric.c /Fo"%OBJPATH%/umfpack_di_load_numeric.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_save_numeric.c /Fo"%OBJPATH%/umfpack_di_save_numeric.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_load_symbolic.c /Fo"%OBJPATH%/umfpack_di_load_symbolic.obj"
cl %CXXFLAGS% /D "DINT" ./UMFPACK/Source/umfpack_save_symbolic.c /Fo"%OBJPATH%/umfpack_di_save_symbolic.obj"
cl %CXXFLAGS% /D "DLONG" /D "CONJUGATE_SOLVE" ./UMFPACK/Source/umf_ltsolve.c /Fo"%OBJPATH%/umf_dl_lhsolve.obj"
cl %CXXFLAGS% /D "DLONG" /D "CONJUGATE_SOLVE" ./UMFPACK/Source/umf_utsolve.c /Fo"%OBJPATH%/umf_dl_uhsolve.obj"
cl %CXXFLAGS% /D "DLONG" /D "DO_MAP" ./UMFPACK/Source/umf_triplet.c /Fo"%OBJPATH%/umf_dl_triplet_map_nox.obj"
cl %CXXFLAGS% /D "DLONG" /D "DO_VALUES" ./UMFPACK/Source/umf_triplet.c /Fo"%OBJPATH%/umf_dl_triplet_nomap_x.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_triplet.c /Fo"%OBJPATH%/umf_dl_triplet_nomap_nox.obj"
cl %CXXFLAGS% /D "DLONG" /D "DO_MAP" /D "DO_VALUES" ./UMFPACK/Source/umf_triplet.c /Fo"%OBJPATH%/umf_dl_triplet_map_x.obj"
cl %CXXFLAGS% /D "DLONG" /D "FIXQ" ./UMFPACK/Source/umf_assemble.c /Fo"%OBJPATH%/umf_dl_assemble_fixq.obj"
cl %CXXFLAGS% /D "DLONG" /D "DROP" ./UMFPACK/Source/umf_store_lu.c /Fo"%OBJPATH%/umf_dl_store_lu_drop.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_assemble.c /Fo"%OBJPATH%/umf_dl_assemble.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_blas3_update.c /Fo"%OBJPATH%/umf_dl_blas3_update.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_build_tuples.c /Fo"%OBJPATH%/umf_dl_build_tuples.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_create_element.c /Fo"%OBJPATH%/umf_dl_create_element.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_dump.c /Fo"%OBJPATH%/umf_dl_dump.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_extend_front.c /Fo"%OBJPATH%/umf_dl_extend_front.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_garbage_collection.c /Fo"%OBJPATH%/umf_dl_garbage_collection.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_get_memory.c /Fo"%OBJPATH%/umf_dl_get_memory.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_init_front.c /Fo"%OBJPATH%/umf_dl_init_front.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_kernel.c /Fo"%OBJPATH%/umf_dl_kernel.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_kernel_init.c /Fo"%OBJPATH%/umf_dl_kernel_init.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_kernel_wrapup.c /Fo"%OBJPATH%/umf_dl_kernel_wrapup.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_local_search.c /Fo"%OBJPATH%/umf_dl_local_search.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_lsolve.c /Fo"%OBJPATH%/umf_dl_lsolve.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_ltsolve.c /Fo"%OBJPATH%/umf_dl_ltsolve.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_mem_alloc_element.c /Fo"%OBJPATH%/umf_dl_mem_alloc_element.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_mem_alloc_head_block.c /Fo"%OBJPATH%/umf_dl_mem_alloc_head_block.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_mem_alloc_tail_block.c /Fo"%OBJPATH%/umf_dl_mem_alloc_tail_block.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_mem_free_tail_block.c /Fo"%OBJPATH%/umf_dl_mem_free_tail_block.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_mem_init_memoryspace.c /Fo"%OBJPATH%/umf_dl_mem_init_memoryspace.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_report_vector.c /Fo"%OBJPATH%/umf_dl_report_vector.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_row_search.c /Fo"%OBJPATH%/umf_dl_row_search.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_scale_column.c /Fo"%OBJPATH%/umf_dl_scale_column.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_set_stats.c /Fo"%OBJPATH%/umf_dl_set_stats.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_solve.c /Fo"%OBJPATH%/umf_dl_solve.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_symbolic_usage.c /Fo"%OBJPATH%/umf_dl_symbolic_usage.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_transpose.c /Fo"%OBJPATH%/umf_dl_transpose.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_tuple_lengths.c /Fo"%OBJPATH%/umf_dl_tuple_lengths.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_usolve.c /Fo"%OBJPATH%/umf_dl_usolve.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_utsolve.c /Fo"%OBJPATH%/umf_dl_utsolve.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_valid_numeric.c /Fo"%OBJPATH%/umf_dl_valid_numeric.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_valid_symbolic.c /Fo"%OBJPATH%/umf_dl_valid_symbolic.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_grow_front.c /Fo"%OBJPATH%/umf_dl_grow_front.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_start_front.c /Fo"%OBJPATH%/umf_dl_start_front.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_store_lu.c /Fo"%OBJPATH%/umf_dl_store_lu.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umf_scale.c /Fo"%OBJPATH%/umf_dl_scale.obj"
cl %CXXFLAGS% /D "DLONG" /D "WSOLVE" ./UMFPACK/Source/umfpack_solve.c /Fo"%OBJPATH%/umfpack_dl_wsolve.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_col_to_triplet.c /Fo"%OBJPATH%/umfpack_dl_col_to_triplet.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_defaults.c /Fo"%OBJPATH%/umfpack_dl_defaults.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_free_numeric.c /Fo"%OBJPATH%/umfpack_dl_free_numeric.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_free_symbolic.c /Fo"%OBJPATH%/umfpack_dl_free_symbolic.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_get_numeric.c /Fo"%OBJPATH%/umfpack_dl_get_numeric.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_get_lunz.c /Fo"%OBJPATH%/umfpack_dl_get_lunz.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_get_symbolic.c /Fo"%OBJPATH%/umfpack_dl_get_symbolic.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_get_determinant.c /Fo"%OBJPATH%/umfpack_dl_get_determinant.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_numeric.c /Fo"%OBJPATH%/umfpack_dl_numeric.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_qsymbolic.c /Fo"%OBJPATH%/umfpack_dl_qsymbolic.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_control.c /Fo"%OBJPATH%/umfpack_dl_report_control.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_info.c /Fo"%OBJPATH%/umfpack_dl_report_info.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_matrix.c /Fo"%OBJPATH%/umfpack_dl_report_matrix.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_numeric.c /Fo"%OBJPATH%/umfpack_dl_report_numeric.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_perm.c /Fo"%OBJPATH%/umfpack_dl_report_perm.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_status.c /Fo"%OBJPATH%/umfpack_dl_report_status.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_symbolic.c /Fo"%OBJPATH%/umfpack_dl_report_symbolic.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_triplet.c /Fo"%OBJPATH%/umfpack_dl_report_triplet.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_report_vector.c /Fo"%OBJPATH%/umfpack_dl_report_vector.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_solve.c /Fo"%OBJPATH%/umfpack_dl_solve.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_symbolic.c /Fo"%OBJPATH%/umfpack_dl_symbolic.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_transpose.c /Fo"%OBJPATH%/umfpack_dl_transpose.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_triplet_to_col.c /Fo"%OBJPATH%/umfpack_dl_triplet_to_col.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_scale.c /Fo"%OBJPATH%/umfpack_dl_scale.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_load_numeric.c /Fo"%OBJPATH%/umfpack_dl_load_numeric.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_save_numeric.c /Fo"%OBJPATH%/umfpack_dl_save_numeric.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_load_symbolic.c /Fo"%OBJPATH%/umfpack_dl_load_symbolic.obj"
cl %CXXFLAGS% /D "DLONG" ./UMFPACK/Source/umfpack_save_symbolic.c /Fo"%OBJPATH%/umfpack_dl_save_symbolic.obj"

rem ===============================================================================================
:link
set LIBFLAGS=%LIBFLAGS% /NOLOGO

echo.
echo Linking 'umfpack.vc12-%1-%2'
lib %LIBFLAGS% /OUT:"..\lib\umfpack.vc12-%1-%2.lib" ..\obj\umfpack.vc12-%1-%2\*.obj
echo.

rem ===============================================================================================
:end