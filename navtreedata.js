/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "FEAT", "index.html", [
    [ "FEAT - Finite Element Analysis Toolbox", "index.html", [
      [ "Basic Information", "index.html#main_sec_basic", null ],
      [ "Build-System", "index.html#main_sec_buildsys", null ],
      [ "Mesh related Topics", "index.html#main_sec_mesh", null ],
      [ "Kernel Subsystems", "index.html#main_sec_subsystems", null ],
      [ "Details about specific topics:", "index.html#main_sec_details", null ]
    ] ],
    [ "Analytic Lambda Expression Functions", "analytic_lambda_function.html", [
      [ "General Information", "analytic_lambda_function.html#ana_lambda_general", null ],
      [ "Lambda Expression Interface", "analytic_lambda_function.html#ana_lambda_format", [
        [ "Additional Parameters in Lambda Expressions", "analytic_lambda_function.html#ana_lambda_params", null ],
        [ "Lambda Expression Thread Safety", "analytic_lambda_function.html#ana_lambda_threads", null ]
      ] ],
      [ "Lambda Function Object Creation", "analytic_lambda_function.html#ana_lambda_create", [
        [ "Richardson Extrapolation of Derivatives", "analytic_lambda_function.html#ana_lambda_derive_extrapol", null ],
        [ "Scalar Examples", "analytic_lambda_function.html#ana_lambda_examples_scalar", null ],
        [ "Vector-Valued Examples", "analytic_lambda_function.html#ana_lambda_examples_vector", null ]
      ] ]
    ] ],
    [ "Coding Conventions", "coding_conventions.html", [
      [ "1. General Conventions", "coding_conventions.html#codeconv_general_conventions", null ],
      [ "2. Naming Conventions", "coding_conventions.html#codeconv_naming_conventions", [
        [ "2.1 Camel-Case", "coding_conventions.html#codeconv_nameconv_camelcase", null ],
        [ "2.2 Lower Snake-Case", "coding_conventions.html#codeconv_nameconv_lower_snakecase", null ],
        [ "2.3 Upper Snake-Case", "coding_conventions.html#codeconv_nameconv_upper_snakecase", null ]
      ] ],
      [ "2. Indentation, Spacing and Line Breaks", "coding_conventions.html#codeconv_indent_spacing", null ],
      [ "3. Preprocessor Directives and Macros", "coding_conventions.html#codeconv_preproc", null ],
      [ "4. Functions", "coding_conventions.html#codeconv_functions", null ],
      [ "5. Classes, Structs and Unions", "coding_conventions.html#codeconv_classes", null ],
      [ "6. Templates", "coding_conventions.html#codeconv_templates", null ],
      [ "7. Namespaces", "coding_conventions.html#codeconv_namespace", null ],
      [ "8. Miscellaneous", "coding_conventions.html#codeconv_misc", null ]
    ] ],
    [ "Frequently Asked Questions", "faq_page.html", [
      [ "Generic FAQ", "faq_page.html#faq_generic", [
        [ "Why are some arrays indexed using the 'int' type and others using the 'FEAT::Index' type?", "faq_page.html#faq_gen_int_vs_index", null ],
        [ "Some source files in the kernel subdirectories have \"eickt\" in their filename. What are these files for?", "faq_page.html#faq_gen_what_is_eickt", null ],
        [ "Should I use the --eickt option for configure or not?", "faq_page.html#faq_gen_use_eickt_or_not", null ]
      ] ],
      [ "LAFEM FAQ", "faq_page.html#faq_lafem", [
        [ "Why can't I use the assignment operator '=' to copy LAFEM containers?", "faq_page.html#faq_lafem_assign_op", null ],
        [ "What is the difference between the clone(), copy() and convert() member functions?", "faq_page.html#faq_lafem_clone_convert_copy", null ]
      ] ]
    ] ],
    [ "INI Data File Format", "ini_format.html", [
      [ "Comments", "ini_format.html#sec_comment", null ],
      [ "Properties", "ini_format.html#sec_property", null ],
      [ "Sections", "ini_format.html#sec_section", null ],
      [ "Examples", "ini_format.html#sec_example", null ]
    ] ],
    [ "LAFEM design description", "lafem_design.html", [
      [ "Overview", "lafem_design.html#lafem_sec_overview", null ],
      [ "Container properties", "lafem_design.html#lafem_sec_container_properties", null ],
      [ "LA Operations", "lafem_design.html#lafem_sec_ops", null ],
      [ "Data handling", "lafem_design.html#lafem_sec_data_handling", null ],
      [ "Container reset", "lafem_design.html#lafem_sec_container_reset", null ],
      [ "Sparse Layout", "lafem_design.html#lafem_sec_sparse_layout", null ],
      [ "Common container typedefs", "lafem_design.html#lafem_sec_typedefs", null ],
      [ "Common typename abbreviations", "lafem_design.html#lafem_sec_template_abbreviation", null ],
      [ "Common Backends", "lafem_design.html#lafem_sec_backends", null ]
    ] ],
    [ "LIKWID marker API with FEAT", "likwid_for_feat.html", [
      [ "LIKWID General Information", "likwid_for_feat.html#LIKWID_descr", null ],
      [ "Prerequisites for using LIWWID", "likwid_for_feat.html#LIKWID_prereques", null ],
      [ "Using the LIKWID marker API", "likwid_for_feat.html#LIKWID_use", null ]
    ] ],
    [ "FEAT Mesh File Format", "mesh_file_format.html", [
      [ "Basic Information", "mesh_file_format.html#meshfile_basic", [
        [ "Mesh Terminology", "mesh_file_format.html#meshfile_basic_terminology", null ],
        [ "Mesh Topology", "mesh_file_format.html#meshfile_basic_topology", [
          [ "Triangle Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_tria", null ],
          [ "Quadrilateral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_quad", null ],
          [ "Tetrahedral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_tetra", null ],
          [ "Hexahedral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_hexa", null ]
        ] ],
        [ "Boundary Parameterization by Charts and Mesh-Parts", "mesh_file_format.html#meshfile_parameterization", [
          [ "Boundary Adaption by Explicit Parameterization", "mesh_file_format.html#meshfile_parameterization_explicit", null ],
          [ "Boundary Adaption by Implicit Parameterization", "mesh_file_format.html#meshfile_parameterization_implicit", null ]
        ] ]
      ] ],
      [ "Overall Mesh File Structure", "mesh_file_format.html#meshfile_file_structure", [
        [ "Mesh Type Specification", "mesh_file_format.html#meshfile_meshtype", null ]
      ] ],
      [ "Mesh File Root XML Node Description", "mesh_file_format.html#meshfile_root", [
        [ "Examples", "mesh_file_format.html#meshfile_root_examples", null ]
      ] ],
      [ "Info XML Node Description", "mesh_file_format.html#meshfile_info", null ],
      [ "Chart XML Node Description", "mesh_file_format.html#meshfile_charts", [
        [ "2D Circle Chart Description", "mesh_file_format.html#meshfile_chart_circle", [
          [ "2D Circle Parameterization", "mesh_file_format.html#meshfile_chart_circle_param", null ],
          [ "2D Circle Chart Examples", "mesh_file_format.html#meshfile_chart_circle_examples", null ]
        ] ],
        [ "3D Sphere Chart Description", "mesh_file_format.html#meshfile_chart_sphere", [
          [ "3D Sphere Chart Examples", "mesh_file_format.html#meshfile_chart_sphere_examples", null ]
        ] ],
        [ "2D Bezier Chart Description", "mesh_file_format.html#meshfile_chart_bezier", [
          [ "Bezier Points", "mesh_file_format.html#meshfile_chart_bezier_points", null ],
          [ "Bezier Parameters", "mesh_file_format.html#meshfile_chart_bezier_params", null ],
          [ "2D Bezier Chart Examples", "mesh_file_format.html#meshfile_chart_bezier_examples", null ]
        ] ],
        [ "3D SurfaceMesh Chart Description", "mesh_file_format.html#meshfile_chart_surfmesh", [
          [ "3D SurfaceMesh Chart Examples", "mesh_file_format.html#meshfile_chart_surfmesh_examples", null ]
        ] ],
        [ "3D Extrude Chart Description", "mesh_file_format.html#meshfile_chart_extrude", [
          [ "Transformation", "mesh_file_format.html#meshfile_chart_extrude_transform", null ],
          [ "3D Extrude Chart Examples", "mesh_file_format.html#meshfile_chart_extrude_examples", null ]
        ] ],
        [ "3D CGAL based SurfaceMesh Chart Description", "mesh_file_format.html#meshfile_chart_CGAL_surfacemesh", null ]
      ] ],
      [ "Mesh XML Node Description", "mesh_file_format.html#meshfile_mesh", [
        [ "Vertices Sub-Node Description", "mesh_file_format.html#meshfile_mesh_vertices", null ],
        [ "Topology Sub-Node Description", "mesh_file_format.html#meshfile_mesh_topology", null ],
        [ "Mesh Examples", "mesh_file_format.html#meshfile_mesh_examples", [
          [ "2D Quadrilateral Unit-Square", "mesh_file_format.html#meshfile_mesh_examples_2d_unit_square_quad", null ],
          [ "2D Triangle Unit-Circle", "mesh_file_format.html#meshfile_mesh_examples_2d_unit_circle_tria", null ],
          [ "3D Hexahedral Unit-Cube", "mesh_file_format.html#meshfile_mesh_examples_3d_cube_hexa", null ]
        ] ]
      ] ],
      [ "MeshPart XML Node Description", "mesh_file_format.html#meshfile_meshpart", [
        [ "Mapping Sub-Node Description", "mesh_file_format.html#meshfile_meshpart_mapping", null ],
        [ "Topology Description", "mesh_file_format.html#meshfile_meshpart_topology", null ],
        [ "Attribute Description", "mesh_file_format.html#meshfile_meshpart_attribute", null ],
        [ "MeshPart Examples", "mesh_file_format.html#meshfile_meshpart_examples", [
          [ "2D Unit-Square Top Edge Mesh-Part", "mesh_file_format.html#meshfile_meshpart_example_unit_square_simple", null ],
          [ "2D Unit-Circle Boundary Mesh-Part", "mesh_file_format.html#meshfile_meshpart_example_unit_circle", null ]
        ] ]
      ] ],
      [ "Partition XML Node Description", "mesh_file_format.html#meshfile_partition", [
        [ "Patch Sub-Node Description", "mesh_file_format.html#meshfile_partition_patch", null ],
        [ "Partition Examples", "mesh_file_format.html#meshfile_partition_examples", [
          [ "2D Unit-Circle Partitions", "mesh_file_format.html#meshfile_partition_examples_unit_circle_tria", null ],
          [ "2D Unit-Square Partitions", "mesh_file_format.html#meshfile_partition_examples_unit_square_quad", null ]
        ] ]
      ] ]
    ] ],
    [ "Mesh Tools Overview", "tools_meshtools.html", [
      [ "Mesh Generator Scripts", "tools_meshtools.html#meshtools_meshgen", [
        [ "2d_quad_circle.py", "tools_meshtools.html#meshtools_2d_quad_circle", null ],
        [ "2d_quad_onion.py", "tools_meshtools.html#meshtools_2d_quad_onion", null ],
        [ "2d_quad_rect.py", "tools_meshtools.html#meshtools_2d_quad_rect", null ],
        [ "2d_quad_ring.py", "tools_meshtools.html#meshtools_2d_quad_ring", null ],
        [ "2d_tria_rect.py", "tools_meshtools.html#meshtools_2d_tria_rect", null ],
        [ "3d_hexa_cube.py", "tools_meshtools.html#meshtools_3d_hexa_cube", null ]
      ] ],
      [ "The mesh-extruder tool", "tools_meshtools.html#meshtools_mesh_extruder", null ],
      [ "The mesh-indexer tool", "tools_meshtools.html#meshtools_mesh_indexer", null ],
      [ "The mesh-partitioner tool", "tools_meshtools.html#meshtools_mesh_partitioner", null ],
      [ "The mesh-validator tool", "tools_meshtools.html#meshtools_mesh_validator", null ],
      [ "The mesh2vtk tool", "tools_meshtools.html#meshtools_mesh2vtk", null ],
      [ "The mesh2eps tool", "tools_meshtools.html#meshtools_mesh2eps", null ],
      [ "The mesh2svg tool", "tools_meshtools.html#meshtools_mesh2svg", null ],
      [ "The mesh2tri tool", "tools_meshtools.html#meshtools_mesh2tri", null ],
      [ "The tri2mesh tool", "tools_meshtools.html#meshtools_tri2mesh", null ]
    ] ],
    [ "Description of the mesh optimization tools", "meshopt_design.html", [
      [ "Basic overview", "meshopt_design.html#basic_overview", null ],
      [ "Workflow", "meshopt_design.html#workflow", null ],
      [ "Implemented classes", "meshopt_design.html#imp_classes", null ]
    ] ],
    [ "Description of the geometric MultiGrid solver", "multigrid_design.html", [
      [ "Multigrid Cycles", "multigrid_design.html#multigrid_cycles", null ],
      [ "W-Cycle Implementation Details", "multigrid_design.html#multigrid_wcycle_impl", null ]
    ] ],
    [ "Syntax overview of ParsedScalar/VectorFunction", "parsed_function.html", [
      [ "Numeric literals", "parsed_function.html#parsed_function_numeric_literals", null ],
      [ "Identifier names", "parsed_function.html#parsed_function_identifier_names", null ],
      [ "The function string syntax", "parsed_function.html#parsed_function_function_string_syntax", null ],
      [ "Inline variables", "parsed_function.html#parsed_function_inline_variables", null ]
    ] ],
    [ "Preprocessor Macros", "preproc_macros.html", [
      [ "Compiler-Detection Macros", "preproc_macros.html#preproc_sec_comp_detect", [
        [ "FEAT_COMPILER", "preproc_macros.html#ppm_feat_compiler", null ],
        [ "FEAT_COMPILER_GNU", "preproc_macros.html#ppm_feat_compiler_gnu", null ],
        [ "FEAT_COMPILER_CLANG", "preproc_macros.html#ppm_feat_compiler_clang", null ],
        [ "FEAT_COMPILER_CRAY", "preproc_macros.html#ppm_feat_compiler_cray", null ],
        [ "FEAT_COMPILER_INTEL", "preproc_macros.html#ppm_feat_compiler_intel", null ],
        [ "FEAT_COMPILER_INTEL_ONEAPI", "preproc_macros.html#ppm_feat_compiler_intel_oneapi", null ],
        [ "FEAT_COMPILER_MICROSOFT", "preproc_macros.html#ppm_feat_compiler_microsoft", null ]
      ] ],
      [ "Backend / Library Macros", "preproc_macros.html#preproc_sec_backends", [
        [ "FEAT_HAVE_ALGLIB", "preproc_macros.html#ppm_feat_have_alglib", null ],
        [ "FEAT_HAVE_CGAL", "preproc_macros.html#ppm_feat_have_cgal", null ],
        [ "FEAT_HAVE_BOOST", "preproc_macros.html#ppm_feat_have_boost", null ],
        [ "FEAT_HAVE_DEATH_HANDLER", "preproc_macros.html#ppm_feat_have_deathhandler", null ],
        [ "FEAT_HAVE_CUDA", "preproc_macros.html#ppm_feat_have_cuda", null ],
        [ "FEAT_HAVE_CUDSS", "preproc_macros.html#ppm_feat_have_cudss", null ],
        [ "FEAT_HAVE_FPARSER", "preproc_macros.html#ppm_feat_have_fparser", null ],
        [ "FEAT_HAVE_FLOATX", "preproc_macros.html#ppm_feat_have_floatx", null ],
        [ "FEAT_HAVE_HALFMATH", "preproc_macros.html#ppm_feat_have_halfmath", null ],
        [ "FEAT_HAVE_HYPRE", "preproc_macros.html#ppm_feat_have_hypre", null ],
        [ "FEAT_HAVE_METIS", "preproc_macros.html#ppm_feat_have_metis", null ],
        [ "FEAT_HAVE_MKL", "preproc_macros.html#ppm_feat_have_mkl", null ],
        [ "FEAT_HAVE_MPI", "preproc_macros.html#ppm_feat_have_mpi", null ],
        [ "FEAT_HAVE_OMP", "preproc_macros.html#ppm_feat_have_omp", null ],
        [ "FEAT_HAVE_PARMETIS", "preproc_macros.html#ppm_feat_have_parmetis", null ],
        [ "FEAT_HAVE_QUADMATH", "preproc_macros.html#ppm_feat_have_quadmath", null ],
        [ "FEAT_HAVE_SUPERLU_DIST", "preproc_macros.html#ppm_feat_have_superlu_dist", null ],
        [ "FEAT_HAVE_TRIANGLE", "preproc_macros.html#ppm_feat_have_triangle", null ],
        [ "FEAT_HAVE_TRILINOS", "preproc_macros.html#ppm_feat_have_trilinos", null ],
        [ "FEAT_HAVE_UMFPACK", "preproc_macros.html#ppm_feat_have_umfpack", null ],
        [ "FEAT_HAVE_ZFP", "preproc_macros.html#ppm_feat_have_zfp", null ],
        [ "FEAT_HAVE_ZLIB", "preproc_macros.html#ppm_feat_have_zlib", null ],
        [ "FEAT_HAVE_ZOLTAN", "preproc_macros.html#ppm_feat_have_zoltan", null ]
      ] ],
      [ "User-Defineable Macros", "preproc_macros.html#preproc_sec_user_def", [
        [ "FEAT_EICKT", "preproc_macros.html#ppm_feat_eickt", null ],
        [ "FEAT_INDEX_U32", "preproc_macros.html#ppm_feat_index_ull", null ],
        [ "FEAT_NO_CONFIG", "preproc_macros.html#ppm_feat_no_config", null ],
        [ "FEAT_OVERRIDE_MPI_OPS", "preproc_macros.html#ppm_feat_ovr_mpi_ops", null ],
        [ "FEAT_MPI_THREAD_MULTIPLE", "preproc_macros.html#ppm_feat_mpi_thread_multiple", null ],
        [ "FEAT_USE_MKL_SPARSE_EXECUTOR", "preproc_macros.html#ppm_feat_mkl_sparse_executor", null ],
        [ "FEAT_UNROLL_BANDED", "preproc_macros.html#ppm_feat_unroll_banded", null ]
      ] ],
      [ "Code-Generation Control Macros", "preproc_macros.html#preproc_sec_codegen", [
        [ "FEAT_PRAGMA_IVDEP", "preproc_macros.html#ppm_feat_pragma_ivdep", null ],
        [ "FEAT_PRAGMA_OMP", "preproc_macros.html#ppm_feat_pragma_omp", null ],
        [ "FORCE_INLINE", "preproc_macros.html#ppm_feat_force_inline", null ],
        [ "NOINLINE", "preproc_macros.html#ppm_feat_noinline", null ]
      ] ],
      [ "Miscellaneous Macros", "preproc_macros.html#preproc_sec_misc", [
        [ "FEAT_DISABLE_WARNINGS / FEAT_RESTORE_WARNINGS", "preproc_macros.html#ppm_feat_warnings", null ],
        [ "FEAT_GIT_SHA1", "preproc_macros.html#ppm_feat_git_sha1", null ],
        [ "FEAT_F128C(x)", "preproc_macros.html#ppm_feat_f128c", null ]
      ] ],
      [ "Build-System Macros", "preproc_macros.html#preproc_sec_build_system", [
        [ "FEAT_DEBUG_MODE", "preproc_macros.html#ppm_feat_debug_mode", null ],
        [ "VISUAL_STUDIO", "preproc_macros.html#ppm_visual_studio", null ],
        [ "FEAT_TESTING_VC", "preproc_macros.html#ppm_feat_testing_vc", null ],
        [ "FEAT_SOURCE_DIR", "preproc_macros.html#ppm_feat_source_dir", null ],
        [ "FEAT_BINARY_DIR", "preproc_macros.html#ppm_feat_binary_dir", null ],
        [ "FEAT_BUILD_ID", "preproc_macros.html#ppm_feat_build_id", null ],
        [ "FEAT_CPU_TYPE", "preproc_macros.html#ppm_feat_cpu_type", null ],
        [ "FEAT_SYSTEM_HOST_COMPILER", "preproc_macros.html#ppm_feat_system_host_compiler", null ],
        [ "FEAT_CUDA_HOST_COMPILER", "preproc_macros.html#ppm_feat_cuda_host_compiler", null ],
        [ "CMAKE_CXX_COMPILER_ID", "preproc_macros.html#ppm_cmake_cxx_compiler_id", null ],
        [ "FEAT_COMPILER_ID", "preproc_macros.html#ppm_feat_compiler_id", null ],
        [ "CMAKE_CXX_COMPILER_VERSION", "preproc_macros.html#ppm_cmake_cxx_compiler_version", null ],
        [ "CMAKE_CXX_COMPILER", "preproc_macros.html#ppm_cmake_cxx_compiler", null ],
        [ "FEAT_USE_COMPILER_WRAPPER", "preproc_macros.html#ppm_feat_compiler_wrapper", null ],
        [ "CMAKE_CXX_COMPILER_ARG1", "preproc_macros.html#ppm_cmake_cxx_compiler_arg1", null ],
        [ "CMAKE_CXX_COMPILER_ARG1_PATH", "preproc_macros.html#ppm_cmake_cxx_compiler_arg1_path", null ],
        [ "CMAKE_CXX_FLAGS", "preproc_macros.html#ppm_cmake_cxx_flags", null ],
        [ "FEAT_CUDA_NVCC_EXECUTABLE", "preproc_macros.html#ppm_cuda_nvcc_executable", null ],
        [ "FEAT_CUDA_NVCC_FLAGS", "preproc_macros.html#ppm_cuda_nvcc_flags", null ],
        [ "FEAT_CUDA_VERSION", "preproc_macros.html#ppm_cuda_version", null ],
        [ "FEAT_CUDA_VERSION_MAJOR", "preproc_macros.html#ppm_cuda_version_major", null ],
        [ "FEAT_CUDA_ARCH", "preproc_macros.html#ppm_cuda_arch", null ],
        [ "MPI_CXX_COMPILER", "preproc_macros.html#ppm_mpi_cxx_compiler", null ],
        [ "MPIEXEC", "preproc_macros.html#ppm_mpiexec", null ],
        [ "CMAKE_MPI_VERSION", "preproc_macros.html#ppm_cmake_mpi_version", null ],
        [ "FEAT_HOSTNAME", "preproc_macros.html#ppm_feat_hostname", null ],
        [ "CMAKE_COMMAND", "preproc_macros.html#ppm_cmake_command", null ],
        [ "CMAKE_GENERATOR", "preproc_macros.html#ppm_cmake_generator", null ]
      ] ]
    ] ],
    [ "Resident vs Transient Reference / Pointer Parameters", "resident_vs_transient.html", [
      [ "Transient Reference/Pointer Parameters", "resident_vs_transient.html#res_vs_tran_transient", null ],
      [ "Resident Reference/Pointer Parameters", "resident_vs_transient.html#res_vs_tran_resident", null ]
    ] ],
    [ "Solver configuration via PropertyMaps", "solver_configuration.html", null ],
    [ "Tool overview", "tools.html", [
      [ "Matrix Info", "tools.html#tools_sec_matrix_info", null ],
      [ "File conversion", "tools.html#tools_sec_io", null ],
      [ "FEM tools", "tools.html#tools_fem", null ],
      [ "Solver tools", "tools.html#tools_solver", [
        [ "Tool for debugging of nonlinear optimizers", "tools.html#tools_subsec_dbg_nlopt", null ]
      ] ]
    ] ],
    [ "The nonlinear optimizer debugging tool", "tools_dbg_nlopt.html", null ],
    [ "The matrix-info Tool", "tools_matrix_info.html", null ],
    [ "FEAT for Visual Studio 2022 on Windows systems", "feat_for_win_vs17.html", [
      [ "Prerequisites", "feat_for_win_vs17.html#win_vs17_prereqs", [
        [ "Microsoft Visual Studio 2022", "feat_for_win_vs17.html#win_vs17_prereqs_vs", [
          [ "Setting the Visual Studio path environment variable", "feat_for_win_vs17.html#win_vs17_pathvar", null ]
        ] ],
        [ "Python 3 Interpreter", "feat_for_win_vs17.html#win_vs17_prereqs_python", null ],
        [ "Microsoft MPI", "feat_for_win_vs17.html#win_vs17_prereqs_mpi", null ],
        [ "Nvidia CUDA", "feat_for_win_vs17.html#win_vs17_prereqs_cuda", null ],
        [ "Intel oneAPI MKL", "feat_for_win_vs17.html#win_vs17_prereqs_imkl", null ]
      ] ],
      [ "Additional Third-Party Libraries", "feat_for_win_vs17.html#win_vs17_thirdparty", null ],
      [ "Creating Visual Studio Project and Solution files", "feat_for_win_vs17.html#win_vs17_vcgen_tool", [
        [ "Creating Kernel Symlink", "feat_for_win_vs17.html#win_vs17_buildsys_vc_symlink", null ],
        [ "Generating Empty Application Projects", "feat_for_win_vs17.html#win_vs17_vcgen_use_empty", null ],
        [ "Generating Simple Application Projects", "feat_for_win_vs17.html#win_vs17_vcgen_use_app_simple", null ],
        [ "Generating Kernel Unit-Test Projects", "feat_for_win_vs17.html#win_vs17_vcgen_use_test", null ],
        [ "The kernel Project file", "feat_for_win_vs17.html#win_vs17_kernel_project", null ]
      ] ],
      [ "Build Configurations and Platforms", "feat_for_win_vs17.html#win_vs17_build_config", null ]
    ] ],
    [ "FEAT Voxel Map File Format", "voxel_map_file_format.html", [
      [ "Basic Information about Voxel Maps", "voxel_map_file_format.html#voxel_map_file_basic", [
        [ "Voxel Map Strides", "voxel_map_file_format.html#voxel_map_strides", null ],
        [ "Voxel Map Indices", "voxel_map_file_format.html#voxel_map_file_indices", null ],
        [ "Voxel Coordinates", "voxel_map_file_format.html#voxel_map_file_coords", null ]
      ] ],
      [ "Voxel Map File Overview", "voxel_map_file_format.html#voxel_map_file_overview", [
        [ "Voxel Map File Header", "voxel_map_file_format.html#voxel_map_file_header", null ],
        [ "Voxel Map File Compression Blocks", "voxel_map_file_format.html#voxel_map_file_compression_blocks", null ]
      ] ]
    ] ],
    [ "FEAT for Visual Studio Code on Unix and Windows Systems", "feat_for_vscode.html", [
      [ "Prerequisites for using FEAT with VScode", "feat_for_vscode.html#vscode_prereques", [
        [ "Prerequisites under unix", "feat_for_vscode.html#unix_vscode_prereques", [
          [ "Special case: Modul based enviroment", "feat_for_vscode.html#unix_vscode_greenhouse", null ]
        ] ],
        [ "Prerequisites under Windows", "feat_for_vscode.html#win_vscode_prereques", [
          [ "Microsoft MPI", "feat_for_vscode.html#win_vscode_prereqs_mpi", null ],
          [ "Nvidia CUDA", "feat_for_vscode.html#win_vscode_prereqs_cuda", null ],
          [ "Intel MKL", "feat_for_vscode.html#win_vscode_prereqs_imkl", null ]
        ] ]
      ] ],
      [ "VScode setup", "feat_for_vscode.html#vscode_setup", null ],
      [ "Configuration Setup", "feat_for_vscode.html#vscode_configuration_setup", null ],
      [ "Working with CMake Tools", "feat_for_vscode.html#vscode_cmake_tools", [
        [ "Configuring FEAT", "feat_for_vscode.html#vscode_cmake_configure", null ],
        [ "Building FEAT", "feat_for_vscode.html#vscode_cmake_build", null ]
      ] ],
      [ "Running FEAT", "feat_for_vscode.html#vscode_cmake_run", null ],
      [ "Working with VScode", "feat_for_vscode.html#vscode_editor", [
        [ "Problems with file discovery", "feat_for_vscode.html#intellisense_filediscovery", null ]
      ] ],
      [ "Using clangd as Language Server", "feat_for_vscode.html#clangd", [
        [ "Setup of clangd", "feat_for_vscode.html#clangd_setup", null ],
        [ "Configuring clangd", "feat_for_vscode.html#clangd_config", null ]
      ] ]
    ] ],
    [ "Todo List", "todo.html", null ],
    [ "Platform depenendend code branches", "platformswitches.html", null ],
    [ "Deprecated List", "deprecated.html", null ],
    [ "Compiler Hacks", "compilerhacks.html", null ],
    [ "Bibliography", "citelist.html", null ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", "namespacemembers_dup" ],
        [ "Functions", "namespacemembers_func.html", "namespacemembers_func" ],
        [ "Variables", "namespacemembers_vars.html", null ],
        [ "Typedefs", "namespacemembers_type.html", null ],
        [ "Enumerations", "namespacemembers_enum.html", null ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", "functions_vars" ],
        [ "Typedefs", "functions_type.html", "functions_type" ],
        [ "Enumerations", "functions_enum.html", null ],
        [ "Related Functions", "functions_rela.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"2d__p1__unrolled_8hpp_source.html",
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html#a7c1c46728c31d39ef595ee3af6580d1f",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html#a667be493fe7fba4d29e225e06cbf5ed4",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol.html",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow_1_1_evaluator.html#a1b624e7d52058a5599660145723d4a45",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_r_h_s2_d.html#aa9afcaf0dcdf11b8f926ae3eddb60033",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_tensor_static.html#a374f763afb91441eb549057949377b1a",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_plane_distance_function_s_d.html",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#afe9068b858f19c48889ba2d13f95a8a4",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#ad3e817886424347711fcfccc4a335f48",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits2.html#a48b1c02d6de9ba7bbea3b418545c1385",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_1_1_evaluator.html#a3445d9c7997262f2e2692c1c5aa769a8",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#afcecbec9a3d1cefb70e672ae895198a2",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a34e3074c9cd4919881367127328809cc",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a34e3074c9cd4919881367127328809cc",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job.html#acc317e0d3cbf48f1ab560754cd4f9e70",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator_1_1_evaluator.html#aa8e4e186b36a7adeb4e2c26ccbf73819",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_functional_1_1_evaluator.html#ac4762c437ba1f46e019dcf9703a75dc5",
"class_f_e_a_t_1_1_assembly_1_1_discrete_function_integral_job_1_1_task.html#ad69b2ff7702e6027a0e1ab74caab0987",
"class_f_e_a_t_1_1_assembly_1_1_force_functional_assembly_job_1_1_task.html",
"class_f_e_a_t_1_1_assembly_1_1_linear_functional_assembly_job_1_1_task.html#ab807fc92cc3192026f1eb0cf1962a626",
"class_f_e_a_t_1_1_assembly_1_1_symbolic_assembler.html#a06564ada37e37d7bf7f169c16b3dd76a",
"class_f_e_a_t_1_1_control_1_1_blocked_unit_filter_system_level.html",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#ac3ce6d36f32fa36f83c5f71f8d2ba38b",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#a4b12d9f60d937a3a5a7e3ce97e27005e",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#af4fef940cc5509f16c2f7acad08e870b",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_system_level.html#ae6619a23cbd3c598f03bdfa38a787e0e",
"class_f_e_a_t_1_1_control_1_1_scalar_combined_system_level.html#a0c231ecf0ca4b2f50026e0cf881237da",
"class_f_e_a_t_1_1_control_1_1_stokes_blocked_unit_velo_none_pres_system_level.html",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_driver_factory_3_01_driver___00_01true_01_4.html#a6add187fd7044dd193a04cdd7c78485e",
"class_f_e_a_t_1_1_dist_1_1_comm.html#a2167a338365aa0a2c4149a9edbf86939",
"class_f_e_a_t_1_1_dist_1_1_status.html#a52171bf1aa1c5afee9c7c44139c28b4c",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier_chart_parser.html#a4b38a70774328011a33830a6c1c02a2d",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#ae87dd290af29f0a0d77c6e69943b15d6",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_surface_mesh.html#a521e739ac61761a63e30ab712b83ccf2",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#ac2307143ce8b8b95e9b09ad6e6321e8f",
"class_f_e_a_t_1_1_geometry_1_1_hit_test_factory.html#a700cc5bff089eb6b005f06d371a982cc",
"class_f_e_a_t_1_1_geometry_1_1_mesh_distortion.html#a92acf0f2c5235da619b658ce1ab50c56",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node_1_1_mesh_part_node_bin.html#a73d042560a4be7e6e6d6ac8ace77b1ef",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#a92c6cf62de2cf7aa41885d70eae66e63",
"class_f_e_a_t_1_1_geometry_1_1_patch_halo_factory_3_01_geometry_1_1_conformal_mesh_3_01_shape___e75c5f271492b3ed0db46ab700c23399.html#abc348e938663ae60e29fddabbdf33b4f",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#a99e8f23a1d0d74c880494580fa076a4d",
"class_f_e_a_t_1_1_geometry_1_1_struct_index_set.html#af31a9b5b569cb242a74891cce473978e",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#a4d48e32cf90368ac14c0c541e74f9c19",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_matrix.html#a2cd30c273fb5308b2a8c4079a90afee1",
"class_f_e_a_t_1_1_global_1_1_matrix.html#a5f7098ae1dd3e547ade4ab6c6fbab720",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar5e6f49c15ccb43840cb894a3375f33c3.html#a357612ed1bb93ffd3625f20137527a40",
"class_f_e_a_t_1_1_global_1_1_synch_scalar_ticket.html#a1630a6104e9b1fa7f7c61d5d06f2f5f4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#a47b718987b62807efb0461015eb18b56",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a23a25e3270a4a0570987d74db1a14817",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a86d0a968a4baec5e597e5fd587439c9b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#a694adc8ae93a07712499d961870189f5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_meta_vector_test_base.html#ae3c2f906f0942c48f4e2b6b4dd181a85",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_d2.html#a3a10660f8f390ac0bed5b39424097572",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#a0d16a14bda7fc7b740b87c6711866b7c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#acc6c844ac68960651eff8bc5a9d75848",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_serial_config.html#a1f25769a2f353a9e6da96da684670879",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#aa0f8a8c78b7e7b226dd476bd9fe924d3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#aa60b2a7f663f102b2fcda559d448176d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#af6adb2d2d65dd0333357ad1f93058c48",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a332cf155ce040c3dade8ef73b0a96335",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a478bbe02aba354da414f5c4b2014a8c0",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_test_matrix_factory.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix.html#a2d440fdd517a165ce379dda46d62eb4d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a926a84217711ab1c22d67fd41e03185a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a9060d0e4e26b6ca72132ec9504a8ec90",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_pow_of_dist.html#a65f1dd685478575a31617909d5e0c62e",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a956e7651d354c5f1018d2eaee9d3ceff",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#a7f2eb50fc47420f4d0025450f2c8156d",
"class_f_e_a_t_1_1_simple_arg_parser.html#a2f3fe2c2f05b25e13544a7d4814ac1e1",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a97bbba38684a4f7352179bd4519872fb",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#aa35ece4c3d8cdb530df482a45c7d6029",
"class_f_e_a_t_1_1_solver_1_1_approximate_hessian_precond.html",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a1f17cb192fab20b3591eafdafd7c5f6e",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_convert_precond.html#ae551eab5ff44d8976ffc77881119b63b",
"class_f_e_a_t_1_1_solver_1_1_diagonal_precond.html#a2be9f4c549b31693fb60a242b2deca00",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_maf2bb55a50c59e69123dfd500b5d8c62b.html#af1b378f730f61d25f3d448f39652b75f",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#ad0e0aa027e766929e7e741f6384f887a",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#ad8716d5245471de11330a2f693d6f4f2",
"class_f_e_a_t_1_1_solver_1_1_generic_umfpack.html#a803cd1653739576cfab7090523f66063",
"class_f_e_a_t_1_1_solver_1_1_hypre_solver_base.html#a9f5d7f67b85c27d9995b31e4604699a8",
"class_f_e_a_t_1_1_solver_1_1_i_l_u_precond_with_backend_3_01_preferred_backend_1_1generic_00_01_785a361efe2523a8ba3c26ca3bd778e8.html#a182fd5d48e8773614562d08c0119d912",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html",
"class_f_e_a_t_1_1_solver_1_1_m_k_l_d_s_s.html#a75eae5bb5307e73b15a06a1f5bc07b6b",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#add5c03ded9d5568e14ec60cdfb1572b8",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_level_base.html#a8c5095a4609a90d158513c0b6d43b789",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#aba7a35edf5c4e1ea1643a3847ef38cf4",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#ac44abdc6decc99909473d0ae2541f08b",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#abd5ba6ffbb459e9be5d0acb50e27ec97",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#ac2ed9d06b8fc5632e175864d6b0aec6c",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#ad586a1b68477198d76f4e794eaa83b5d",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#a43bbf72e681f56ef9a6a96af37265a76",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#aa7cf18c4f6e62acfd0c3cf3285c00b96",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a1f17cb192fab20b3591eafdafd7c5f6e",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#a27d93b800c82dae6cf39f5ffc9a94015",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#a26458def0ef8654e3fff8f444df8108f",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a7d583ff48c7678f1cffa30fdd48315a7",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#ad5a06be09b1b40ab2c5b7cceb8274ab8",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#a4c11796c53851c2e9f599168917ef6ec",
"class_f_e_a_t_1_1_solver_1_1_saddle_umfpack_mean.html#a1fc13a753cfd7a6ba9b4d4dc893493f7",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a9694e5d8b7db0c434b6f8f7890bdd5fc",
"class_f_e_a_t_1_1_solver_1_1_umfpack.html#af5ea439102e6a6b6db2d3bb5b434af34",
"class_f_e_a_t_1_1_solver_1_1_voxel_ama_vanka.html#a76aca114b5b79bb9db75c116ae04b1ed",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_element_1_1_dof_assignment.html",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_element.html#a8b224fb72515bc1bdc61dafccad1c546",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_evaluator_3_01_space___00_01_trafo_evaluator_cd584e2ae53d2b4306ff680d7265a0ec.html#a928c1fc25e4c56643b88d1a7510cb40f",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___0b6ac4939673a6cd90495f3edd5624cc3.html#a0e7f8c0fa4765783dc8aa9bd85deb982",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_964e26ad63d9177b3f747659af96be6e.html#ae3968cef6969d7edfd3e29680c987fb6",
"class_f_e_a_t_1_1_space_1_1_dof_mapping_base.html#a7f2b5590e436c6a2f6de8d3082d30fdf",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_element.html#a7a00a22acd058302b50e2e3d05768e5c",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercube_a7111ca547438b901c926f188c55fe99.html#abc7087cde691d01f94b0084c96966334",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_se3129d8853eb9690234f254237b7f1ec.html#a1f5b8556e9fa1ee98993eeebeb5337c3",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_node_functional.html#abc7087cde691d01f94b0084c96966334",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sb3aa60112044c986d3e79f1b1403fdf7.html",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_evaluator.html",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_2be5aeec58e5b7d19e38f7fb4e901a36.html#abbeace51bb9218c5002d90ce9f7ba59b",
"class_f_e_a_t_1_1_statistics.html#ab2ab9f302b49ea12380d12547aea5449",
"class_f_e_a_t_1_1_thread_fence.html#ae8e2106f34f6f5be326fd1747a28656b",
"class_f_e_a_t_1_1_tiny_1_1_vector.html#a6476874775adb0e63c7ed826f3838fca",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_011_01_4_01_4.html#a44b64fa989b225892d2e869ac480c707",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_012_01_4_01_4.html#a973a38ad683103552dade007a17cb255",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_6f4e5b7f5f9d6e2d04743ee3976d1657.html#a638be8b19ba18b9b933b485e9a209847",
"class_f_e_a_t_1_1_trafo_1_1_mapping_base.html#a226bc9a13bd705855ae2b7935fa59401",
"class_f_e_a_t_1_1_voxel_assembly_1_1_voxel_burgers_assembler_3_01_q2_standard_hyper_cube_3_01dimcd264f458d76378c7602bad0426c2704.html#a19f9d6de05ee74ac73ce6e2b8ac5d2e0",
"class_f_e_a_t_1_1_xml_1_1_syntax_error.html#aad4387e6feb2068a99620dc8a8ced21f",
"function_8hpp_source.html",
"likwid__marker_8hpp_source.html",
"namespace_f_e_a_t_1_1_control_1_1_domain.html#a26f4622eb76d3da5d2d2e40650ec744ca1b7d5726533ab525a8760351e9b5e415",
"namespace_f_e_a_t_1_1_shape.html#a13f39565eb0b9391f364f92fd10eb4bd",
"namespace_f_e_a_t_1_1_voxel_assembly.html#a637d205e4e3d706dcc1d6a5fb619c0ed",
"row__norm__generic_8hpp_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#a7258cd7a7e2e67eefd55115862b942fe",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html#a195637fb1d5ed7b4e1422d275679cdb3",
"struct_f_e_a_t_1_1_space_1_1_argyris_1_1_node_functional_3_01_space___00_011_00_01_data_type___01_4_1_1_value.html",
"struct_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_node_functional_3_01_space___00_01_shape_type___00_01e7800677649f5ff35e76794b0e5e8cf0.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_eval_policy.html#a08461895b388d9937fd61ac2aa1d6e37",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01short_01_4.html#a1712c2a5a21fd7507f66abcee40e3a26"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';