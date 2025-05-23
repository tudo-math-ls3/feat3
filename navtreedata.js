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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html#a29f5261d5491c3fc190caa5a7e6eb7f8",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_frankes_function_1_1_evaluator.html#a1eb375034ddde8778d98406e71725128",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow.html#adf290f410d3d0978ddaa37c587128768",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_r_h_s2_d.html#a4e6c4a413e4e6078ac48ea09aa23992e",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_taylor_green_vortex_velo2_d_1_1_evaluator.html#ae4df19352cefc54ac4335428ac187eda",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_inverse_distance_function_1_1_evaluator.html#ad7a88326f8fd783170fa1bf5cea987ff",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#aceff256971abdac8875f1efea819d6e7",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a9fbe5c00df64b0b5fc137bb44365f154",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits2.html#a03c3ede75963845a3dd9702f62f5b224",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_1_1_evaluator.html",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#ad74b37015de3ba9088143fdb0c175a1a",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a25a4188596a0fd8552e17da078005e36",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a16cf785893294825c0949ea32b1ec969",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job.html#aa6598c1a041c7e2d45ec8962469cb4d2",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator_1_1_evaluator.html#a4dfc1807680a2ab83c090f87c3283ce7",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_functional_1_1_evaluator.html#a98c728b49b2ba064f46fff286b9cc0f0",
"class_f_e_a_t_1_1_assembly_1_1_discrete_function_integral_job_1_1_task.html#aa4746f1277b9f518283bd84c3338cb8d",
"class_f_e_a_t_1_1_assembly_1_1_force_functional_assembly_job.html#ab92890acd137fa134026349564c5c58a",
"class_f_e_a_t_1_1_assembly_1_1_linear_functional_assembly_job_1_1_task.html#a76b20a47c59a1a614c5f8434bce71934",
"class_f_e_a_t_1_1_assembly_1_1_stokes_f_b_m_assembler.html#afa1b0161c50ccdfc15b9e35a4d642f30",
"class_f_e_a_t_1_1_control_1_1_blocked_combined_system_level.html#ae375e8ccc7302013a235e80d1ab6aa8d",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#abaa5c7803a9d1ab3fc22a9a891ed68a9",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#a426bf2533f64c2ebafc4112b72626347",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#adf85be4249212b934df605b6b2a484da",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_system_level.html#ad9489e469b3eabe96587dbfecb3362a6",
"class_f_e_a_t_1_1_control_1_1_scalar_basic_system_level.html#aed3c7548f8c24dc41c012b20e3755802",
"class_f_e_a_t_1_1_control_1_1_stokes_blocked_system_level.html#a2c72ce340611881422f31b48eb3f1b2e",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_driver_factory_3_01_driver___00_01false_01_4.html#ab6d81ad1b08f35b98df9ad79f3331a21",
"class_f_e_a_t_1_1_dist_1_1_comm.html#a1472ad3b2745a0a7aca762068b39a9d9",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#af6e17a1fc1be4ef74af643aba3e2ddc1",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_mesh.html#a7fc8aa693d506f8059454c26314f8b6a",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#aad0b74b5f1bbbb9888cdf5388d476dd7",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#a4c5036b7c5368e34ae94a6bb33468093",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#aea82c135894a81131dcc68316cac37a6",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#a492cc6c48fd1a5797000cff811799b8e",
"class_f_e_a_t_1_1_geometry_1_1_global_masked_boundary_factory.html#a16e12db98573f4717cc9ac329a71aff2",
"class_f_e_a_t_1_1_geometry_1_1_intern_1_1_parti_iterative_individual.html#a25dad35b5a056a554241ded9fc84f491",
"class_f_e_a_t_1_1_geometry_1_1_mesh_file_reader.html#a88297dfd82fbfda44b44d3e7cf32b195",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part.html#ae75ff5e5a33f114cbb9d26f958722b21",
"class_f_e_a_t_1_1_geometry_1_1_parti_iterative_3_01_conformal_mesh_3_01_shape___00_01num__coords___00_01_coord___01_4_01_4.html#a4d81a35eda3a71acae7e737238c94da0",
"class_f_e_a_t_1_1_geometry_1_1_patch_mesh_part_splitter_3_01_geometry_1_1_conformal_mesh_3_01_shba0a87c4bb4c4d2c925ee5e6ab30bd83.html#a585e54ea2b3ad870d18ebb12ffe1e515",
"class_f_e_a_t_1_1_geometry_1_1_schneiders_data.html#a7f0601539a3f8e8f7c21e214e403d071",
"class_f_e_a_t_1_1_geometry_1_1_struct_unit_cube_factory.html",
"class_f_e_a_t_1_1_geometry_1_1_template_builder.html#aba6adc71b5a28fe9659b96c62fefbee6",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#a635bc81c98aca4f29d7e67cb943f57ab",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_matrix.html#a6bc3f8a151741f223c2e74c81ed1f056",
"class_f_e_a_t_1_1_global_1_1_matrix.html#aafc7cee7eceb0956e1f29b2a6816d18c",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar5e6f49c15ccb43840cb894a3375f33c3.html#a425ab19b989dcd0832de43fa63f04c28",
"class_f_e_a_t_1_1_global_1_1_synch_scalar_ticket.html#aa5f2e6d6a455ae40dcb6910d64c487d2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#a7915be74ff0aba8d97b530cf87c663b6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a3d2a2536488a06809e3b435c0d8435a4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a96a306b2187e5def6102928d6516e92a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#a7e1173858100512c2a6922e0dacadc03",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_none_filter.html#a4c3decf2cdc29901916fa38bb11db43a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_d2.html#ae21822b403259f3ccde707cbf33d1559",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#a413a8c59b519fb93b962d8bbad1d3c19",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#adec6314c4974dd31804f803b70569d37",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_slip_filter.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#aacece407ebb7a34807eb1f317d459b69",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#ab6ce85e7ab759a914f20985aa13b79e2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a43b3061774550efb2221a851f690393e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a719dfbaf1b0d869d46d710093e6daf94",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_test_matrix_factory.html#a5c1319a992be7d06de192abf840bf804",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix.html#a5ab97585d4307608f9bec397aea0b721",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#aaf518e5245006923971cabaa7d929d9d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a9e02b51ad88f3c4dc81c6077108b3302",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_pow_of_dist.html#af45d9df55a1627e83d69185146501f22",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#aa59f2fbec36c0c37482bb4d13afb693b",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#aaad404cd571a207f3c855e30acc559aa",
"class_f_e_a_t_1_1_simple_arg_parser.html#a99f09de8dec9e2cd253508bf975459b8",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#aa7281f8e3498245d411c5829b9d0cec0",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_approximate_hessian_precond.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#a2908f6e22afbf7c0de75b88d872f6885",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#adeda7bbb6c2498a790d4bee34052b0c2",
"class_f_e_a_t_1_1_solver_1_1_descent.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_maf2bb55a50c59e69123dfd500b5d8c62b.html#a11ab45df003dfac0adf117ab3cdcd015",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a8965935265c8ba7609b4461b0a5672e0",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#aa95d79033c0912710a988a61e7312a15",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_hessian_precond.html#ac687f31fa9389b4c0c0ec211bc262e0c",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#af415e6a1a3687fef310f8b2012a6b040",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#af5be9eafbcdd5d8d6b43d4781698f198",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#ad8026d89fbcbd8d65ffdc207d15b06f6",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#aaec597375935042fdeee5ff49d9e5dba",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#aad34701b722b596d0e0140f5c799d2d5",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a87caaf864220f9cb8add2b0eb90e2ad2",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a87caaf864220f9cb8add2b0eb90e2ad2",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a9694e5d8b7db0c434b6f8f7890bdd5fc",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a973d3942df838d68f63666265e73b025",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#aeebdc6fb0b2e5072cb5c0c52205ff58b",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#ad0e0aa027e766929e7e741f6384f887a",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html#a90679ca33476caf229d7f21d97270eba",
"class_f_e_a_t_1_1_solver_1_1_polynomial_precond.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a43bbf72e681f56ef9a6a96af37265a76",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#aa49a59066c84ccc87b0249d8f59d99f3",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#af415e6a1a3687fef310f8b2012a6b040",
"class_f_e_a_t_1_1_solver_1_1_s_s_o_r_precond.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a5e276e3ab4c84a834377997024bf746e",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#ad94d705e3aefd0e7709db1ef698de19a",
"class_f_e_a_t_1_1_solver_1_1_vanka_factor_error.html#ac922497e9e4a4a2e340e7173bb4d4780",
"class_f_e_a_t_1_1_space_1_1_basis_data_reduced.html#a087ac5432b656a7349650a3e7bcd874c",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercub4124820b6987228b00243ae6979c030a.html",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_element_1_1_evaluator.html#a66e38aca4ff839852987f450a76e70e6",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___03bb1a3d92225206453ee71c00db5aa7a.html#a3c7283f42e823561cd3cc5da70bb5108",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element_1_1_evaluator.html",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_single_entity.html#a3b9da3c3b96ce741c7f022c784ed68e1",
"class_f_e_a_t_1_1_space_1_1_evaluator_base.html#a24159b15f72a5584d653894370369e3a",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spd08b4bf16ff97c315f2f45da544757eb.html#a89576bb728b70079c6be450b671dd945",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sbede9c86e9ef3e51ccbac5494ff49399.html#a1f5b8556e9fa1ee98993eeebeb5337c3",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_evaluator.html#a4b1ac4e614fad291b9dd1b93491c8342",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s5e9945e9db3d72bec06849dc7e957367.html#a1f5b8556e9fa1ee98993eeebeb5337c3",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_element.html#a57d3800a4007b3aea0b184440d3df44b",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_2be5aeec58e5b7d19e38f7fb4e901a36.html#a0ad9c5786977c9dd1fe7dabc17e88d1a",
"class_f_e_a_t_1_1_statistics.html#a3dd21536839595cd7f03d71ba46065eb",
"class_f_e_a_t_1_1_test_system_1_1_test_list_1_1_delete_on_destruction.html",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#ac885c405d542d655a38f0f8d0e1c1a3e",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping_data.html",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#ae378db672d86ce801a6bbd96d0b3b029",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#abd66371dcb2b499cdba2dac74e59d997",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_mapping.html",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_mapping.html#ab029a0459681a10c64549ca9a280ea64",
"class_f_e_a_t_1_1_xml_1_1_scanner.html#a1e25aa44e6ed481fbc0f5061f52c2cd2",
"empty__voxel_8cpp_source.html",
"kernel_2lafem_2base_8hpp.html#a82e983d75b2911bc2cd014a541d50711a6adf97f83acf6453d4a6a4b1070f3754",
"namespace_f_e_a_t.html#ad574925e947c823619c356d5d12b0815aa41eb1b7074026b30ca99409df3704a0",
"namespace_f_e_a_t_1_1_math.html#a61139d6ac642cd3298e36e4ad43fe277",
"namespace_f_e_a_t_1_1_solver.html#adab9653730c87172104d5178417d8fa1",
"preproc_macros.html#ppm_cmake_mpi_version",
"struct_f_e_a_t_1_1_assembly_1_1_discrete_function_integral_3_01_l_a_f_e_m_1_1_dense_vector_block282ab9ff04649d6c90f61fc974e9149f.html",
"struct_f_e_a_t_1_1_geometry_1_1_adaption_stats.html#a3fdc1c91631cc2a7cc587ae8aed6abbc",
"struct_f_e_a_t_1_1_geometry_1_1_subdivision_level_tuple.html#a0f42b6005b711fed79efdb4bcb6f178c",
"struct_f_e_a_t_1_1_solver_1_1_optimization_test_traits.html",
"struct_f_e_a_t_1_1_space_1_1_lagrange3_1_1_dof_traits_3_01_shape_1_1_simplex_3_01shape__dim___01_4_00_012_01_4.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___00d0928204f7a6d17533bb19e7a72add33.html#ae2312efefedc551cea24f44386fe1e5e",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01int_01_4.html#ae0b90fe096e3b4445e5a6fe2a3b1b6a5",
"vanka_8hpp_source.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';