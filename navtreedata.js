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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive.html#a67f96d76def386406d59a39512a06761",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_frankes_function_1_1_evaluator.html#a11f180a7f977049207ea1e0979e8149b",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow.html#aaf22fee3705228e0ae511a18259d9e2c",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_pres2_d_1_1_evaluator.html#ac4ad534e51aa598e3cede929f47d3915",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_taylor_green_vortex_velo2_d_1_1_evaluator.html#a9634dc077cdba105ce70844a09e90e8e",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_inverse_distance_function_1_1_evaluator.html#a67df40393f52e0362f9d63f1a9052e89",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#a8f1ef61ad6383304720f6b55c4ecc0e4",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a37457e9f74d8b5b3ee346602b8097aff",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits1.html#a7a5873ee37868be4044eb4813b812b5c",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#af19d13d192d42ceec47f03db6e14bc65",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a25a4188596a0fd8552e17da078005e36",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_assembly_job_base.html#a0cd20e62f56caf3d8f2e6da35ad2ae9a",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_matrix_assembly_job.html#a3630dde9b818165f948842dbab26e331",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_vector_assembly_job_1_1_task.html#aa480d81a445d3936c72f812b256e9c1c",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_scalar_matrix_assembly_job_1_1_task.html#a75f1fccb2290c861e63210141d95bbed",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a56faffd7256b9d2b897a3ee00f631328",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job.html#adfb1add4c6a0ef5730d58fc1c417e598",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator.html#ab8b5dccda304af814ff1fbda7a8ba66a",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_functional.html#ad8a040955136a818eda272c3768c827d",
"class_f_e_a_t_1_1_assembly_1_1_discrete_cell_projector.html",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_basic_matrix_task_c_r_t_p1.html#a764b7ab7d09ab969256190ae146a762a",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_bilinear_operator_matrix_job1_1_1_task.html#a887caab2010f6165562c7dc7b1aa66a0",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_force_functional_vector_job_1_1_task.html#ae38716258fdf9f6bbfbf7b249d30b833",
"class_f_e_a_t_1_1_assembly_1_1_linear_functional_1_1_evaluator.html#a63997f4f05078c674fe46aa007ae5de8",
"class_f_e_a_t_1_1_assembly_1_1_stokes_f_b_m_assembler.html#a0ae2ae18a4d6237b14a8a5791a1da429",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_matrix_task_c_r_t_p1.html#a950d578cb6bfdb43e5a6dcc789523605",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_task_base2.html#ac95ca9809432c1da68dca77e7a213d0d",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_bilinear_operator_matrix_job1_1_1_task.html#a47949ff39f179d72b985099ad0792ce3",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_discrete_function_integral_job_1_1_task.html#a9c384a8d33215f36b13a291ed6f87c59",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_mass_matrix_job_1_1_task.html#ae50c47be9234530ba390dfd25256db59",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_task_base1.html#abb2644a5bea1a9875360c0761daba65d",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_stokes_vector_analysis_task_c_r_t_p.html#a0535ebea4561af3a70391a24065e934f",
"class_f_e_a_t_1_1_control_1_1_blocked_unit_filter_system_level.html",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#acaaef6f1e579cc13a14ee455e1f7774e",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#a48c57a3b04a4af44e7d3fd710ad09785",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#adf85be4249212b934df605b6b2a484da",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_system_level.html#ad9489e469b3eabe96587dbfecb3362a6",
"class_f_e_a_t_1_1_control_1_1_scalar_basic_system_level.html#aed3c7548f8c24dc41c012b20e3755802",
"class_f_e_a_t_1_1_control_1_1_stokes_blocked_system_level.html#a2c72ce340611881422f31b48eb3f1b2e",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_driver_factory_3_01_driver___00_01false_01_4.html#ab6d81ad1b08f35b98df9ad79f3331a21",
"class_f_e_a_t_1_1_dist_1_1_comm.html#a1472ad3b2745a0a7aca762068b39a9d9",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#af6e17a1fc1be4ef74af643aba3e2ddc1",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_mesh.html#a7fc8aa693d506f8059454c26314f8b6a",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#aa4ead2d485a4c8048a16218249c9e8a9",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#a47580cb7576d19ff11916e3dcca2cc55",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#ae9bc655c31ce3aac351586a347a3222e",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#a031cb5054b3fd5f767847c71ff4bf2a2",
"class_f_e_a_t_1_1_geometry_1_1_factory_3_01_conformal_mesh_3_01_shape___00_01num__coords___00_01_coord_type___01_4_01_4.html#a32a12cc3d66729adb2d1a2b82cb7a887",
"class_f_e_a_t_1_1_geometry_1_1_intern_1_1_adaptive_mesh_storage.html#a8562ffdfd3c597de754f761ce8a5f6ae",
"class_f_e_a_t_1_1_geometry_1_1_mesh_extruder_3_01_conformal_mesh_3_01_shape_1_1_hypercube_3_012_d55947478d0079b597f439c47758a098.html#a112204ffda7b0d95e0b6c6b0383f3677",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part.html#a0f8acb17161ebd7596977a8ded829aca",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#a48bf531fd2b533d0499f3875d884ffe9",
"class_f_e_a_t_1_1_geometry_1_1_partition.html#a4ca2b3db4768a461663794d782cd521b",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#a37d2690edd188b7da1081ef716f037be",
"class_f_e_a_t_1_1_geometry_1_1_standard_refinery_3_01_mesh_part_3_01_parent_mesh___01_4_01_4.html#a84ce9c8de04f1f03655fc196b892804c",
"class_f_e_a_t_1_1_geometry_1_1_sun_zhao_ma_expansion_data.html",
"class_f_e_a_t_1_1_geometry_1_1_unit_cube_patch_generator.html",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#afdad224c3b4221b7c718eef25bd45a9b",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_system.html#ae5320c9c0136a70257dd765c3e1226a5",
"class_f_e_a_t_1_1_global_1_1_muxer.html#a41a4e05af9cd6be92674ba1d73f6c2eb",
"class_f_e_a_t_1_1_global_1_1_splitter.html#a09d1406888b33c985bc2950b6c546796",
"class_f_e_a_t_1_1_global_1_1_transfer.html#abcad3df455e97eb9b058f32b966a206a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a44130fc02d29b0b3560b101191962e5c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#aa2c36f6aadece3f8451a3a3cae8a15fd",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#af5dc00e2137f8c37f5462b3423f7bd0d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#a4c8588083209fc7a5dd44fb29ebba84f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html#a82dd6f1f97f59ef3740cd9d142824bb5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_full_matrix.html#a837e3f5c3bc0f1c45981f1c1d3a4ac9d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#a6a5266ffee69b98876cde386197f46b6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_slip_filter.html#afdf76eedf179ae4589847f886d5a6101",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#aea7cd62e111b9a08e7f9f009c641e120",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#af4875a200ee0fa3f2b473780b63666c7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#a35f74a32d4565b6ce185d665918440ff",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a74c9279f1115605cd6ac24be25d4bdef",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a958343c0926a9f67dffc5997c3c05b31",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_transfer.html#a1e90eb773a357faefd7c6baa50b09584",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix.html#ac91cafb5a7a9499f5d52b18961e36b3b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter.html#a3c3ae498f41f107c70283b4fae7697b5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#ad8ba52e668a73138f2fd8bbb9238ed65",
"class_f_e_a_t_1_1_meshopt_1_1_du_dv_functional.html#a447d737515d86497d9614357d57d183d",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#ad30a3161633b42d49a92df6d31b00ad7",
"class_f_e_a_t_1_1_meshopt_1_1_rumpf_functional.html#a800ba01856798ca8df8de933c42464c9",
"class_f_e_a_t_1_1_solver_1_1_a_d_p_solver_base.html#a394e71e99c87e99a5eda16c302df725e",
"class_f_e_a_t_1_1_solver_1_1_a_d_p_solver_base_3_01_global_1_1_p_m_d_c_d_s_c_matrix_3_01_matrix_040b7e1caa6985ce9572fb7c7833a8b9.html#af2e2fe610ff0c3b63fe23b73ed0316b7",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#afaab86a331f64d630f134f3c6027b5f7",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka.html",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#a1f17cb192fab20b3591eafdafd7c5f6e",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#ae1a798867103a1ab77313c5f429626ca",
"class_f_e_a_t_1_1_solver_1_1_convert_precond.html#aa04fb62591a5379e3ca4088c0775d6a9",
"class_f_e_a_t_1_1_solver_1_1_direct_sparse_solver.html#a12f5142bae88b96eed8fdfea372c9342",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#aa49f8a3ff5c644668d3c83c5e84bd7e2",
"class_f_e_a_t_1_1_solver_1_1_euclid_precond.html#afdc8520e21ad255c8c577b9ec7893c6a",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a1678798c09d463bd121bf52812719ba9",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#a27d93b800c82dae6cf39f5ffc9a94015",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#a5e276e3ab4c84a834377997024bf746e",
"class_f_e_a_t_1_1_solver_1_1_hypre_solver_base.html#afdc8520e21ad255c8c577b9ec7893c6a",
"class_f_e_a_t_1_1_solver_1_1_i_l_u_precond_with_backend_3_01_preferred_backend_1_1generic_00_01_aed6c8770199c660f2b23c3052158c26.html#acd49cd3547318e90690c227c56a29d89",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a12edcf32734ef6628a3320fe7653007f",
"class_f_e_a_t_1_1_solver_1_1_m_k_l_d_s_s.html#aefd52f3a0eb0b20c62dc780734722771",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#aeebdc6fb0b2e5072cb5c0c52205ff58b",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_level_std.html",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#ac93194a978f521b1a3e4ce546a3cd757",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#ad586a1b68477198d76f4e794eaa83b5d",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#acf7a6a13468eae67e3446f004032bf21",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#ad586a1b68477198d76f4e794eaa83b5d",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#aec76c75b2d2582ef2cf0ed928b583314",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#a576f7cf3e1c847ac5611e477bd551b73",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#abfeec9e7a424dee82dccb22be9f43fd4",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a288adb7c820ebbb82a0a8b07f178db13",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html#ae6bedcbf11e2f8710f3c0a9b577a4e35",
"class_f_e_a_t_1_1_solver_1_1_precon_wrapper.html#a70018a5e4ca2cd3dfc0cac9db3df051e",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a576f7cf3e1c847ac5611e477bd551b73",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#ab12c4d5da695f0c49d2ebd82eeccffcd",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_s_s_o_r_precond.html#af513f0b3a23a450609db8e887de91c1e",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a6e0074b3a336238364169fb1a0fb0701",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#ab934585a0eb28b50149e32ff156685bc",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#aa80981f7205f164f6f7578c6e8fc3d44",
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_node_functional_3_01_space___00_010_00_01_data_type___01_4.html",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_82866b3cc738de4ae4935135fac5e732.html#ab70b19b6962179079f8a2d38ddddcd28",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_element.html#a8b224fb72515bc1bdc61dafccad1c546",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___001b011599bed78529504db50cb49456d.html#ac25e50167de479a49d61c8fa0ddbd179",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element.html#a8b224fb72515bc1bdc61dafccad1c546",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_base.html#a8b994c0cb6ba94e5475b21cea15768a0",
"class_f_e_a_t_1_1_space_1_1_element_base_1_1_dof_assignment.html",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#ad5886f4911cac95096215def29085034",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#a6db7e5c3954c49646d2823108329675d",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_element_1_1_dof_assignment.html#a219f4619cd071ad0f25275f84aecafb6",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html#a6bfdd412851516ebd9e42a508b51edb6",
"class_f_e_a_t_1_1_space_1_1_node_functional_base.html#ad42b3539fd9aefad4ad90a6f93b72fba",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_element.html#ae0abe6963c2356d0f7450af8544ef1ab",
"class_f_e_a_t_1_1_space_1_1_standard_vector_eval_traits.html#a97dfebf83204b64b549294d9b31e25df",
"class_f_e_a_t_1_1_syntax_error.html#a76bcfc0b1c2ed8cc4d368d67c2c5c193",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#a37fbaec87d0ddc1f7f9f9266fda53405",
"class_f_e_a_t_1_1_trafo_1_1_evaluator_base.html#acdb7bfb197ee64e4d83a5901a0aaa456",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#a3712e4ae60ddf267304ff196ce8de8e3",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#a04385656ab1f3e606550f78430ae4b23",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#a4d4f718b82c4369a953bf9b2f21564b7",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a97164a88685009b2dfcd657bc4e65679",
"class_f_e_a_t_1_1_xml_1_1_content_error.html#a7aee2ab09e1e1840e107a63579f9db63",
"cuda__util_8cu_source.html",
"functions_type_g.html",
"mesh_file_format.html#meshfile_chart_bezier_examples",
"namespace_f_e_a_t_1_1_dist.html#aed6f8b30f1c9c545ac8f0f1fae41a5d3",
"namespace_f_e_a_t_1_1_solver.html#a2865a0f7572372119b5ab0a71b357e28",
"namespace_f_e_a_t_1_1_voxel_assembly.html#a637d205e4e3d706dcc1d6a5fb619c0ed",
"random_8hpp_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_scalar_error_info.html#ac25e6ad602e2839a14fd2be4c49dda5b",
"struct_f_e_a_t_1_1_geometry_1_1_conformal_mesh_1_1_index_set.html#aff6a33e9af6db7e2e498a6cd743cd1a8",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html#af6d3939eede3073c49191a7276cf1066",
"struct_f_e_a_t_1_1_space_1_1_bernstein2_1_1_dof_traits_3_01_shape_1_1_hypercube_3_01shape__dim___01_4_00_01face__dim___01_4.html#a75bcea2a2ec25e3d2e18b0e15b52cfaa",
"struct_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_dof_traits.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_eval_policy.html#a7c1585285806b11a4f6215bd073517b1",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01short_01_4.html#af1ce789c9e07b4db55a27bbdc0fec213"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';