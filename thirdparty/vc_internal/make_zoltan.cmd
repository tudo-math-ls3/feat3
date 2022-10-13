@echo off
@rem FEAT3: Finite Element Analysis Toolbox, Version 3
@rem Copyright (C) 2010 - 2022 by Stefan Turek & the FEAT group
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
set OBJPATH=..\obj\zoltan.%1-%2-%3

if exist %OBJPATH% (
  del %OBJPATH%\*.obj
) else (
  mkdir %OBJPATH%
)

rem Set Compiler flags
set CXXFLAGS=%CXXFLAGS% /c /GS /Gd /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo /TC
set CXXFLAGS=%CXXFLAGS% /wd"4028" /wd"4018" /wd"4090" /wd"4101" /wd"4102" /wd"4113" /wd"4244" /wd"4267" /wd"4305" /wd"4311" /wd"4477" /wd"4716" /wd"4996"
set CXXFLAGS=%CXXFLAGS% /I"%MSMPI_INC% " /I"./Zoltan/src" /I"./Zoltan/src/include" /I"./Zoltan/src/all" /I"./Zoltan/src/Utilities/shared" /I"./Zoltan/src/Utilities/Timer" /I"./Zoltan/src/zz" /I"./Zoltan/src/coloring" /I"./Zoltan/src/graph" /I"./Zoltan/src/ha" /I"./Zoltan/src/hier" /I"./Zoltan/src/hsfc" /I"./Zoltan/src/lb" /I"./Zoltan/src/matrix" /I"./Zoltan/src/order" /I"./Zoltan/src/par" /I"./Zoltan/src/params" /I"./Zoltan/src/phg" /I"./Zoltan/src/rcb" /I"./Zoltan/src/reftree" /I"./Zoltan/src/simple" /I"./Zoltan/src/timer" /I"./Zoltan/src/tpls"
set CXXFLAGS=%CXXFLAGS% /Fd"../obj/zoltan.%1-%2-%3/zoltan.%1-%2-%3.pdb"
set CXXFLAGS=%CXXFLAGS% /Fp"../obj/zoltan.%1-%2-%3/zoltan.%1-%2-%3.pch"

if "%MSMPI_INC%" == "" (
  echo.
  echo MPI not found; skipping Zoltan library...
  echo.
  goto end
)

echo.
echo Compiling Zoltan Sources...
cl %CXXFLAGS% ./Zoltan/src/all/all_allo.c /Fo"%OBJPATH%/all_allo.obj"
cl %CXXFLAGS% ./Zoltan/src/coloring/coloring.c /Fo"%OBJPATH%/coloring.obj"
cl %CXXFLAGS% ./Zoltan/src/coloring/color_test.c /Fo"%OBJPATH%/color_test.obj"
cl %CXXFLAGS% ./Zoltan/src/coloring/bucket.c /Fo"%OBJPATH%/bucket.obj"
cl %CXXFLAGS% ./Zoltan/src/coloring/g2l_hash.c /Fo"%OBJPATH%/g2l_hash.obj"
cl %CXXFLAGS% ./Zoltan/src/graph/graph.c /Fo"%OBJPATH%/graph.obj"
cl %CXXFLAGS% ./Zoltan/src/ha/divide_machine.c /Fo"%OBJPATH%/divide_machine.obj"
cl %CXXFLAGS% ./Zoltan/src/ha/get_processor_name.c /Fo"%OBJPATH%/get_processor_name.obj"
cl %CXXFLAGS% ./Zoltan/src/ha/ha_ovis.c /Fo"%OBJPATH%/ha_ovis.obj"
cl %CXXFLAGS% ./Zoltan/src/hier/hier.c /Fo"%OBJPATH%/hier.obj"
cl %CXXFLAGS% ./Zoltan/src/hier/hier_free_struct.c /Fo"%OBJPATH%/hier_free_struct.obj"
cl %CXXFLAGS% ./Zoltan/src/hsfc/hsfc_box_assign.c /Fo"%OBJPATH%/hsfc_box_assign.obj"
cl %CXXFLAGS% ./Zoltan/src/hsfc/hsfc.c /Fo"%OBJPATH%/hsfc.obj"
cl %CXXFLAGS% ./Zoltan/src/hsfc/hsfc_hilbert.c /Fo"%OBJPATH%/hsfc_hilbert.obj"
cl %CXXFLAGS% ./Zoltan/src/hsfc/hsfc_point_assign.c /Fo"%OBJPATH%/hsfc_point_assign.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_balance.c /Fo"%OBJPATH%/lb_balance.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_box_assign.c /Fo"%OBJPATH%/lb_box_assign.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_copy.c /Fo"%OBJPATH%/lb_copy.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_eval.c /Fo"%OBJPATH%/lb_eval.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_free.c /Fo"%OBJPATH%/lb_free.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_init.c /Fo"%OBJPATH%/lb_init.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_invert.c /Fo"%OBJPATH%/lb_invert.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_migrate.c /Fo"%OBJPATH%/lb_migrate.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_part2proc.c /Fo"%OBJPATH%/lb_part2proc.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_point_assign.c /Fo"%OBJPATH%/lb_point_assign.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_remap.c /Fo"%OBJPATH%/lb_remap.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_set_fn.c /Fo"%OBJPATH%/lb_set_fn.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_set_method.c /Fo"%OBJPATH%/lb_set_method.obj"
cl %CXXFLAGS% ./Zoltan/src/lb/lb_set_part_sizes.c /Fo"%OBJPATH%/lb_set_part_sizes.obj"
cl %CXXFLAGS% ./Zoltan/src/matrix/matrix_build.c /Fo"%OBJPATH%/matrix_build.obj"
cl %CXXFLAGS% ./Zoltan/src/matrix/matrix_distribute.c /Fo"%OBJPATH%/matrix_distribute.obj"
cl %CXXFLAGS% ./Zoltan/src/matrix/matrix_operations.c /Fo"%OBJPATH%/matrix_operations.obj"
cl %CXXFLAGS% ./Zoltan/src/matrix/matrix_sym.c /Fo"%OBJPATH%/matrix_sym.obj"
cl %CXXFLAGS% ./Zoltan/src/matrix/matrix_utils.c /Fo"%OBJPATH%/matrix_utils.obj"
cl %CXXFLAGS% ./Zoltan/src/order/order.c /Fo"%OBJPATH%/order.obj"
cl %CXXFLAGS% ./Zoltan/src/order/order_struct.c /Fo"%OBJPATH%/order_struct.obj"
cl %CXXFLAGS% ./Zoltan/src/order/order_tools.c /Fo"%OBJPATH%/order_tools.obj"
cl %CXXFLAGS% ./Zoltan/src/order/hsfcOrder.c /Fo"%OBJPATH%/hsfcOrder.obj"
cl %CXXFLAGS% ./Zoltan/src/order/perm.c /Fo"%OBJPATH%/perm.obj"
cl %CXXFLAGS% ./Zoltan/src/par/par_average.c /Fo"%OBJPATH%/par_average.obj"
cl %CXXFLAGS% ./Zoltan/src/par/par_bisect.c /Fo"%OBJPATH%/par_bisect.obj"
cl %CXXFLAGS% ./Zoltan/src/par/par_median.c /Fo"%OBJPATH%/par_median.obj"
cl %CXXFLAGS% ./Zoltan/src/par/par_median_randomized.c /Fo"%OBJPATH%/par_median_randomized.obj"
cl %CXXFLAGS% ./Zoltan/src/par/par_stats.c /Fo"%OBJPATH%/par_stats.obj"
cl %CXXFLAGS% ./Zoltan/src/par/par_sync.c /Fo"%OBJPATH%/par_sync.obj"
cl %CXXFLAGS% ./Zoltan/src/par/par_tflops_special.c /Fo"%OBJPATH%/par_tflops_special.obj"
cl %CXXFLAGS% ./Zoltan/src/params/assign_param_vals.c /Fo"%OBJPATH%/assign_param_vals.obj"
cl %CXXFLAGS% ./Zoltan/src/params/bind_param.c /Fo"%OBJPATH%/bind_param.obj"
cl %CXXFLAGS% ./Zoltan/src/params/check_param.c /Fo"%OBJPATH%/check_param.obj"
cl %CXXFLAGS% ./Zoltan/src/params/free_params.c /Fo"%OBJPATH%/free_params.obj"
cl %CXXFLAGS% ./Zoltan/src/params/key_params.c /Fo"%OBJPATH%/key_params.obj"
cl %CXXFLAGS% ./Zoltan/src/params/print_params.c /Fo"%OBJPATH%/print_params.obj"
cl %CXXFLAGS% ./Zoltan/src/params/set_param.c /Fo"%OBJPATH%/set_param.obj"
cl %CXXFLAGS% ./Zoltan/src/tpls/build_graph.c /Fo"%OBJPATH%/build_graph.obj"
cl %CXXFLAGS% ./Zoltan/src/tpls/postprocessing.c /Fo"%OBJPATH%/postprocessing.obj"
cl %CXXFLAGS% ./Zoltan/src/tpls/preprocessing.c /Fo"%OBJPATH%/preprocessing.obj"
cl %CXXFLAGS% ./Zoltan/src/tpls/scatter_graph.c /Fo"%OBJPATH%/scatter_graph.obj"
cl %CXXFLAGS% ./Zoltan/src/tpls/third_library.c /Fo"%OBJPATH%/third_library.obj"
cl %CXXFLAGS% ./Zoltan/src/tpls/verify_graph.c /Fo"%OBJPATH%/verify_graph.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_build.c /Fo"%OBJPATH%/phg_build.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_build_calls.c /Fo"%OBJPATH%/phg_build_calls.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg.c /Fo"%OBJPATH%/phg.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_lookup.c /Fo"%OBJPATH%/phg_lookup.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_verbose.c /Fo"%OBJPATH%/phg_verbose.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_coarse.c /Fo"%OBJPATH%/phg_coarse.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_comm.c /Fo"%OBJPATH%/phg_comm.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_distrib.c /Fo"%OBJPATH%/phg_distrib.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_gather.c /Fo"%OBJPATH%/phg_gather.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_hypergraph.c /Fo"%OBJPATH%/phg_hypergraph.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_match.c /Fo"%OBJPATH%/phg_match.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_order.c /Fo"%OBJPATH%/phg_order.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_parkway.c /Fo"%OBJPATH%/phg_parkway.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_patoh.c /Fo"%OBJPATH%/phg_patoh.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_plot.c /Fo"%OBJPATH%/phg_plot.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_rdivide.c /Fo"%OBJPATH%/phg_rdivide.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_refinement.c /Fo"%OBJPATH%/phg_refinement.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_scale.c /Fo"%OBJPATH%/phg_scale.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_serialpartition.c /Fo"%OBJPATH%/phg_serialpartition.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_util.c /Fo"%OBJPATH%/phg_util.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_tree.c /Fo"%OBJPATH%/phg_tree.obj"
cl %CXXFLAGS% ./Zoltan/src/phg/phg_Vcycle.c /Fo"%OBJPATH%/phg_Vcycle.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/box_assign.c /Fo"%OBJPATH%/box_assign.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/create_proc_list.c /Fo"%OBJPATH%/create_proc_list.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/inertial1d.c /Fo"%OBJPATH%/inertial1d.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/inertial2d.c /Fo"%OBJPATH%/inertial2d.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/inertial3d.c /Fo"%OBJPATH%/inertial3d.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/point_assign.c /Fo"%OBJPATH%/point_assign.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/rcb_box.c /Fo"%OBJPATH%/rcb_box.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/rcb.c /Fo"%OBJPATH%/rcb.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/rcb_util.c /Fo"%OBJPATH%/rcb_util.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/rib.c /Fo"%OBJPATH%/rib.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/rib_util.c /Fo"%OBJPATH%/rib_util.obj"
cl %CXXFLAGS% ./Zoltan/src/rcb/shared.c /Fo"%OBJPATH%/shared.obj"
cl %CXXFLAGS% ./Zoltan/src/reftree/reftree_build.c /Fo"%OBJPATH%/reftree_build.obj"
cl %CXXFLAGS% ./Zoltan/src/reftree/reftree_coarse_path.c /Fo"%OBJPATH%/reftree_coarse_path.obj"
cl %CXXFLAGS% ./Zoltan/src/reftree/reftree_hash.c /Fo"%OBJPATH%/reftree_hash.obj"
cl %CXXFLAGS% ./Zoltan/src/reftree/reftree_part.c /Fo"%OBJPATH%/reftree_part.obj"
cl %CXXFLAGS% ./Zoltan/src/simple/block.c /Fo"%OBJPATH%/block.obj"
cl %CXXFLAGS% ./Zoltan/src/simple/cyclic.c /Fo"%OBJPATH%/cyclic.obj"
cl %CXXFLAGS% ./Zoltan/src/simple/random.c /Fo"%OBJPATH%/random.obj"
cl %CXXFLAGS% ./Zoltan/src/timer/timer_params.c /Fo"%OBJPATH%/timer_params.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_exchange_sizes.c /Fo"%OBJPATH%/comm_exchange_sizes.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_invert_map.c /Fo"%OBJPATH%/comm_invert_map.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_do.c /Fo"%OBJPATH%/comm_do.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_do_reverse.c /Fo"%OBJPATH%/comm_do_reverse.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_info.c /Fo"%OBJPATH%/comm_info.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_create.c /Fo"%OBJPATH%/comm_create.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_resize.c /Fo"%OBJPATH%/comm_resize.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_sort_ints.c /Fo"%OBJPATH%/comm_sort_ints.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_destroy.c /Fo"%OBJPATH%/comm_destroy.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Communication/comm_invert_plan.c /Fo"%OBJPATH%/comm_invert_plan.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Timer/zoltan_timer.c /Fo"%OBJPATH%/zoltan_timer.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Timer/timer.c /Fo"%OBJPATH%/timer.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Memory.c /Fo"%OBJPATH%/DD_Memory.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Find.c /Fo"%OBJPATH%/DD_Find.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Destroy.c /Fo"%OBJPATH%/DD_Destroy.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn3.c /Fo"%OBJPATH%/DD_Set_Neighbor_Hash_Fn3.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Remove.c /Fo"%OBJPATH%/DD_Remove.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Create.c /Fo"%OBJPATH%/DD_Create.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Update.c /Fo"%OBJPATH%/DD_Update.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Stats.c /Fo"%OBJPATH%/DD_Stats.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Hash2.c /Fo"%OBJPATH%/DD_Hash2.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Print.c /Fo"%OBJPATH%/DD_Print.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn2.c /Fo"%OBJPATH%/DD_Set_Neighbor_Hash_Fn2.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Set_Hash_Fn.c /Fo"%OBJPATH%/DD_Set_Hash_Fn.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn1.c /Fo"%OBJPATH%/DD_Set_Neighbor_Hash_Fn1.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/Memory/mem.c /Fo"%OBJPATH%/mem.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/shared/zoltan_align.c /Fo"%OBJPATH%/zoltan_align.obj"
cl %CXXFLAGS% ./Zoltan/src/Utilities/shared/zoltan_id.c /Fo"%OBJPATH%/zoltan_id.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_coord.c /Fo"%OBJPATH%/zz_coord.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_gen_files.c /Fo"%OBJPATH%/zz_gen_files.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_hash.c /Fo"%OBJPATH%/zz_hash.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/murmur3.c /Fo"%OBJPATH%/murmur3.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_map.c /Fo"%OBJPATH%/zz_map.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_heap.c /Fo"%OBJPATH%/zz_heap.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_init.c /Fo"%OBJPATH%/zz_init.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_obj_list.c /Fo"%OBJPATH%/zz_obj_list.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_rand.c /Fo"%OBJPATH%/zz_rand.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_set_fn.c /Fo"%OBJPATH%/zz_set_fn.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_sort.c /Fo"%OBJPATH%/zz_sort.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_struct.c /Fo"%OBJPATH%/zz_struct.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_back_trace.c /Fo"%OBJPATH%/zz_back_trace.obj"
cl %CXXFLAGS% ./Zoltan/src/zz/zz_util.c /Fo"%OBJPATH%/zz_util.obj"


echo.
echo Linking 'zoltan.%1-%2-%3'
lib %LIBFLAGS% /OUT:"..\lib\zoltan.%1-%2-%3.lib" ..\obj\zoltan.%1-%2-%3\*.obj
echo.

rem ===============================================================================================
:end