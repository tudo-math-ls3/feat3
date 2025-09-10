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
    [ "FEATS's Thirdparty Dependencies", "feat_tpls.html", [
      [ "Reference", "feat_tpls.html#Reference", null ],
      [ "Enabling thirdparty dependencies", "feat_tpls.html#enabling_tpls", null ],
      [ "Resolving thirdparty dependencies", "feat_tpls.html#including_tpls", [
        [ "Finding thirdparty dependencies", "feat_tpls.html#including_tpls_finding", null ],
        [ "Downloading thirdparty dependencies", "feat_tpls.html#including_tpls_dowloading", null ]
      ] ],
      [ "For FEAT developers", "feat_tpls.html#developing_feat", null ]
    ] ],
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
          [ "Special case: Module based enviroment", "feat_for_vscode.html#unix_vscode_greenhouse", null ]
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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html#a3e2f70d9769a220323eacfc07a5b219c",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html#a0745dcc8cb89511d3c281b4992ce2565",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_frankes_function_1_1_evaluator.html#a209a38a4c0310ba9bbed6a6c4402c97e",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow_1_1_evaluator.html",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_r_h_s2_d.html#a605d243593f2d8541b9370566296c305",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_tensor_static.html",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_inverse_distance_function_1_1_evaluator.html#ae4df19352cefc54ac4335428ac187eda",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#ae5a4e6789890d4b833c0021860181c33",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#aa17f28e8cbbe153fab1d9d3f2e6a0699",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits2.html#af735b9f13a5ee9411fb85b4b2b442b2c",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_task_base.html#ad1ffa1d8dd3067c503578fdace46faed",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#aeb81ae5049b5cd987c5782c183ea3dbd",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#af0e331577ab72f707e03979b9fd4115a",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#a9c3ea52ca671041ee3ebcebafda9a13f",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_force_functional_1_1_evaluator.html#ab8571c66e8d93996d3f2045b9a6dc1a6",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_normal_gradient_trial_operator_blocked_1_1_evaluator.html#a11beb9bb5d154eb3761026e3f3dd0925",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a8096ae609313c8bac0e9b8812563ad91",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_basic_vector_task_c_r_t_p.html#a2bced2b7e1c71aabe159a7557d094e51",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_discrete_function_integral_job.html#ae69644f00a07fd63929de43090a1060b",
"class_f_e_a_t_1_1_assembly_1_1_function_integral_info.html",
"class_f_e_a_t_1_1_assembly_1_1_slip_filter_assembler.html#a4c6ab9d4fc71eaf1a92b56f8a4eff3ec",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_matrix_task_c_r_t_p1.html#a3ff6b62efe7a1c0be3df3fa4c2c1708d",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_task_base2.html#a35c73286bd7103e41a2449b1ad486ad9",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_bilinear_operator_apply_vector_job1_1_1_task.html#ae5a44dd9ec9a2aeedac0318a8818ab93",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_discrete_function_integral_job_1_1_task.html#a3ff6b62efe7a1c0be3df3fa4c2c1708d",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_mass_matrix_job_1_1_task.html#a896dfe2496917efb5e7d1d74a6e26751",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_task_base1.html",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_stokes_body_force_assembly_job_1_1_task.html#aa7e3ae883cb3ca99ca3f21b648f28cf5",
"class_f_e_a_t_1_1_control_1_1_blocked_combined_system_level.html#a33daee45bd7f59e999af8ef9e108f258",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#a800d3c72b7195f11cf496daf02e4ae71",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#a208fd72cd4a598bbf31e6452e107f912",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#aa68124254f24405dd2288c0aae3689a3",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_system_level.html#a8618adc26f27180a24d54f33a4cde70d",
"class_f_e_a_t_1_1_control_1_1_scalar_basic_system_level.html#a1fd0d52031673e65105c9895088fbc0c",
"class_f_e_a_t_1_1_control_1_1_stokes_blocked_slip_unit_velo_mean_pres_system_level.html#ab13df046ed96bbd03b4ca40e3aeb5904",
"class_f_e_a_t_1_1_cubature_1_1_refine_factory.html",
"class_f_e_a_t_1_1_cubature_1_1_trapezoidal_driver_3_01_shape_1_1_simplex_3_01dim___01_4_01_4.html#a914bb935ae6dda79ed9a4434ec91df57",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#a9789ad04018b1345170995052e567275",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_mesh.html#a2848f6a64116a9d0db177dd3e0b3860b",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a5bebdfed21494d0b74603ed724c6d601",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_chart_c_r_t_p.html#aed7ebf57c13d6fe2b19177dea31ee588",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#a8cc9395304ca6033d55d0fc1eafff38d",
"class_f_e_a_t_1_1_geometry_1_1_chart_hit_test_factory.html#aa019f738bc569ba8abb67bb6398d2ef7",
"class_f_e_a_t_1_1_geometry_1_1_export_v_t_k.html#af740fc82d811cb97e5456084c9ed05a5",
"class_f_e_a_t_1_1_geometry_1_1_intern_1_1_adaptive_mesh_storage.html#a37439b790ec5fe6377aa1d1a789f7bba",
"class_f_e_a_t_1_1_geometry_1_1_mesh_distortion.html#adc33d6694ec6457399f1e184900fad4f",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node_linker_error.html",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#a28d7c5907981ef38c4a0192a1622d9f7",
"class_f_e_a_t_1_1_geometry_1_1_partition.html#a232e09c505a6e1d72c88c4a5991a2ba6",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#a1c3fd615b09f5f30b14518eaa46ce02d",
"class_f_e_a_t_1_1_geometry_1_1_standard_refinery_3_01_mesh_part_3_01_parent_mesh___01_4_01_4.html#a3f9be236f37bf207a9212a308440163e",
"class_f_e_a_t_1_1_geometry_1_1_sun_zhao_ma_data.html#aaf6d1ff74afeb770c3e2984be936ab4a",
"class_f_e_a_t_1_1_geometry_1_1_umbrella_smoother.html#aa25b46f410c9df6f5ab6f0a2ca31bc04",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#ae713d2e252aef7d9efbe43677315baf3",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_vector.html#ae889bac9ab7946098ad674c2cef58974",
"class_f_e_a_t_1_1_global_1_1_muxer.html#a856bfd983d139d5439df8c4dc5a61f6f",
"class_f_e_a_t_1_1_global_1_1_splitter.html#a488ba393a93d662c7dc2a0494b83283f",
"class_f_e_a_t_1_1_global_1_1_transfer.html#afac360cc74545f7cdca892d3aae8736d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a6466f8f82437e54a1ee25eb142a70cee",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#ada2b8bd276f75ab9df2a22cb52f48e7c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_chain.html#a5e78e76f6644369d5c5c909f6c4daf94",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter.html#a70e4d91539875674f5f6786bd77e9470",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#a6dc9b5b7c2af872d327732d5ccb1453f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html#ab02dfcf70747f4b71d205f90f2662ba4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_full_matrix.html#adc5f4b467c22ead5fc1a24e639585cb4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#aaa8d347869f563a30524bedf8373ffa1",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a17a23310907ad8282f0e93f309551e2f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a119f63d4e3ff2c01eb77f29139e82635",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a37fa6870946a1dbe0774074952555123",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#a86415f91eb3b42033341fa003af22742",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#aabfe5da81fcd1564e731e36a576d481d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#a0cd2ba5d16debab9ef084cd18c27fd68",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_diag_matrix.html#a394e10315e566336fd77dce4ea3407dd",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix_row.html#ab1c12b1d919d0868ebc0b796be3df299",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter_blocked.html#a318772496d3749014b2f4c3a19884ccc",
"class_f_e_a_t_1_1_memory_usage.html#acb8b5721a409a571c1c3532316a6beeb",
"class_f_e_a_t_1_1_meshopt_1_1_du_dv_functional.html#ab8933ac2b533b57f0ae8e9ced72daf6e",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function.html#a59d1fa97f2b7df0580135fa474bf00f6",
"class_f_e_a_t_1_1_property_map.html#a075a35a04f549de653b983d4313c7225",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a00adcb4d7171bac9a8b415dc286954a1",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a0810de389dc450bea5cb6917b66ca051",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka.html#a0a4d513a3b85b917939eb25e7a2d9516",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#a1feb792955098cf117423f18dfd0db8f",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a5b0b165b5aba580f1b68bd3c2a49933a",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_cu_solver_l_u.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#a388fea3b3a144ccc59b6cdf59f709cbd",
"class_f_e_a_t_1_1_solver_1_1_euclid_precond.html#a8c34d07a55cd09eaafef49269f140536",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a031ba8078fe8170cc0c79bcc6a900341",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#a40db3688d9794b0334a0226994b17c05",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#a27d93b800c82dae6cf39f5ffc9a94015",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#a16f527ebfa7b2e516f4e8cd72a38c4d7",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a2daf610e7977a9e96ff84a6e027293d0",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a1d005a8fa11c118322cb35c7f1297155",
"class_f_e_a_t_1_1_solver_1_1_matrix_precond.html#ad0f65bbe800bcef97fbd1e9de0079d53",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_level_std.html#a93d89e7941fe4c31d059845a1501a451",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#ae142305afb0934862700c46e18620f66",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#af8e655f0cdc31c3b245587b58a14d5b3",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#aec87d011de4859420af0f06559fbef2e",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#aeebdc6fb0b2e5072cb5c0c52205ff58b",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#a77f32fe5475099d5ae331eb7b18d9280",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#aec76c75b2d2582ef2cf0ed928b583314",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a5e276e3ab4c84a834377997024bf746e",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#a6fd0f61a9458e3b6c32b4c502fcdc7f6",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#a6fd0f61a9458e3b6c32b4c502fcdc7f6",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#ac59ca39a0b808cbab77a4b8b028abb7b",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#a0fd9bce681b27352bb6a37c2ef10e300",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#a77f32fe5475099d5ae331eb7b18d9280",
"class_f_e_a_t_1_1_solver_1_1_scale_precond.html#a6863ca5947832845346c0055e0d92353",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#abfeec9e7a424dee82dccb22be9f43fd4",
"class_f_e_a_t_1_1_solver_1_1_uzawa_precond.html#a23fe1659d7117d9dd8818900cf592b0a",
"class_f_e_a_t_1_1_solver_1_1_voxel_ama_vanka.html#af5f950a7e9a34eadc1b73e63927af189",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_045e90cd69ef865eaf5b37b181a99332.html#a4b1ac4e614fad291b9dd1b93491c8342",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_evaluator_3_01_space___00_01_trafo_evaluator__802e55c907078015a390ad2e449f9e5b.html#a2644f5cb97bfdc7ce8f437a8c9176b4b",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_node_functional_3_01_space___00_01_shape_1_1_d215266554069458711a15ecd1479e21.html#a12dd52ec74b75aff5d9471e1c3572485",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___0b6ac4939673a6cd90495f3edd5624cc3.html#a9c6e5d8443bd6ac814fa42e8fc698099",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_b2c6b42d6f9b9e25d469d3575a373a61.html#ab1c4da11e6e9a9a48a9cb6de2f52a72f",
"class_f_e_a_t_1_1_space_1_1_dof_mapping_single_entity.html#a0e8c50272c4dcf6f23bfc44bb26c9f7e",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_element_1_1_evaluator.html#af00f51b518f4ffb3043df3a4abb731de",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_element.html#a8981c83488d337c11302106fa96d02d9",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_se3129d8853eb9690234f254237b7f1ec.html#ae4d17e132c111387fbe169819fc1c2aa",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_element.html#a8d82bce6274b32899746656ddb66160c",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sb3aa60112044c986d3e79f1b1403fdf7.html#ae4d17e132c111387fbe169819fc1c2aa",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s88a2ce84347f190b629ca7916f4787f8.html#ab70b19b6962179079f8a2d38ddddcd28",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_3da8886277cafb8343dcd8258f6373f1.html#a4d3ca23df7379587ca10bb9087e583ac",
"class_f_e_a_t_1_1_stop_watch.html",
"class_f_e_a_t_1_1_tiny_1_1_matrix.html#a1c34fdf3c13f041e4ce8ff4a1c4c245e",
"class_f_e_a_t_1_1_tiny_1_1_vector.html#ac437fd2ef50881a1608e5f79e34b0b9e",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_011_01_4_01_4.html#abd605a87e6152b8167a01c21a7d1f4bb",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a2b61e7d246e324460ea95d351ed4471e",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_6f4e5b7f5f9d6e2d04743ee3976d1657.html#ae59246c38c8f44464b8171d7e43e3d69",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator.html#a357cba752a903cb7f4f6727315e4ea59",
"class_f_e_a_t_1_1_voxel_assembly_1_1_voxel_burgers_assembler_3_01_q2_standard_hyper_cube_3_01dimcd264f458d76378c7602bad0426c2704.html#a7eda66687b41623799694f676ca15f93",
"coding_conventions.html#codeconv_misc",
"functions_func.html",
"lumping__generic_8hpp_source.html",
"namespace_f_e_a_t_1_1_control_1_1_domain.html#a26f4622eb76d3da5d2d2e40650ec744c",
"namespace_f_e_a_t_1_1_pack.html#a7c69c0366e113fd42901ba9bcea55d4dab992b60ea8f82c69c5c7581c1ae96080",
"namespace_f_e_a_t_1_1_tiny.html#a79a96f7f8b812ae88911f5a1ac3c4baf",
"preproc_macros.html#ppm_visual_studio",
"struct_f_e_a_t_1_1_assembly_1_1_scalar_error_info.html#a5d143c00d25f626301f20ad8025462da",
"struct_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere_traits.html#aa130188fcee55b4be01861ed7a3a7132",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html#aa7700e2529d379aa801194002d0600c3",
"struct_f_e_a_t_1_1_space_1_1_argyris_1_1_node_functional_3_01_space___00_010_00_01_data_type___01_4_1_1_value.html",
"struct_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_node_functional_3_01_space___00_01_shape_type___00_0196f09c74bf71ba3bcb4729c098016440.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_eval_policy.html",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01short_01_4.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';