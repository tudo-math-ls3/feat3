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
"class_f_e_a_t_1_1_analytic_1_1_auto_derive.html#a67fdcf3e6dfffbb044d802d6ab13093f",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_frankes_function_1_1_evaluator.html#a1770796225f7fe17fd349e5b0fc148fc",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow.html#ac43934be653cb23ff0ba83de2f67d357",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_pres2_d_1_1_evaluator.html#ae4df19352cefc54ac4335428ac187eda",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_taylor_green_vortex_velo2_d_1_1_evaluator.html#ab6fd0a7ee32d7ad4aff743e337a10a28",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_inverse_distance_function_1_1_evaluator.html#a713f3de35e7abcba2d97e81ad687845a",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#aac360f552a5009cf9e86667c88198826",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a4c9916f051e50f1f48678db79548a675",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits1.html#a8486ea3feaecf5278fd1e6cf92edf975",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_task_base.html",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a3c4b8515f1afe1b4036d12169ced1f38",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_assembly_job_base.html#a46dabf841d82815a815ea464ff078f5f",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_matrix_assembly_job.html#a89f43b87d30e6920925dbe75dcf22ac4",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_vector_assembly_job_1_1_task.html#acede4391e02e7de869cc6597c7631b44",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_scalar_matrix_assembly_job_1_1_task.html#ab8cf5a81feb521eb9a559b05d2d0501a",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a8e6d7a03ac19e89bd11b490fd99c81ff",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#a2f57e23777c8d214b414ee41407f5887",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator_1_1_evaluator.html#ad0521bf4b7169457b4b75249ce3fe46f",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_functional_1_1_evaluator.html#ac6dcdb4b3b28a8e517f3a49eba70474e",
"class_f_e_a_t_1_1_assembly_1_1_discrete_evaluator.html#ac4d5435f58b216a7c66a5a9604579582",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_basic_matrix_task_c_r_t_p1.html#af92684cd70fdb18ac96e502d41051a16",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_bilinear_operator_matrix_job2.html#a66d0b040a7fe0b90ccaa1d25e7ab7fb7",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_job_1_1_task.html#ac46bb91fd5d6dc1421bb283a6ad174ab",
"class_f_e_a_t_1_1_assembly_1_1_oldroyd_assembler.html#ac4e6ebad11d10177b5f4018aaae836c1",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembler.html#a82a98629c15d565ef5f51b3f70d9d1f5",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_task_base1.html#a0bcd7951a7742b0321872ce5d856e5e2",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_bilinear_operator_apply_vector_job1_1_1_task.html#a273daf6dcb8161c4a705e7f87ca3ce32",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_bilinear_operator_matrix_job2_1_1_task.html#a9a5a8a94a30c33e370bf3627d2f7335f",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_job_1_1_task.html#a1e0014a4d2e03af4b6dfd21368f50646",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_stabilization_matrix_job_1_1_task.html#a46d32c2822c84d5496ce0a387fb3dcca",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_stokes_body_force_assembly_job_1_1_task.html",
"class_f_e_a_t_1_1_assembly_1_1_velocity_analyser.html#a4e7e4ab0d8f47466fc70365053df30a4",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_hierarch_unit_cube_domain_control2.html#abaa5c7803a9d1ab3fc22a9a891ed68a9",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control_base.html#afc0f046709c34a9d79407ca29b28dd9f",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_level_wrapper.html#a0ec069bf4f63a020c39d1310c2a58502",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_control_base.html#a83d805b11c78229a3be5137e13722469",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_quadratic_system_level.html#a3a7ece2bbe00811b932bf9f74d91bbb6",
"class_f_e_a_t_1_1_control_1_1_scalar_unit_filter_system_level.html#a530bc0001ac286d083b3df677f0114cc",
"class_f_e_a_t_1_1_cubature_1_1_dunavant_driver_3_01_shape_1_1_simplex_3_012_01_4_01_4.html#a8fe0b70f58b273631c40d512e4cd3235",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_trapezoidal_driver.html#a03e297ab5b6ded4744152d085898c910",
"class_f_e_a_t_1_1_dist_1_1_comm.html#afb8e0aad6059b40d0ef1334cb2c15bb0",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_export_v_t_u_3_01_geometry_1_1_adaptive_mesh_3_01_templab17a45e1da2bf1be268079986731df4a.html#a69e3dd1c273f6a57631afd3037a1696d",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_neighbors.html#aac33e630fec64ef342f4e9825e45114c",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_chart_base.html#a59dbc0be992c90afde38bc9549481f5e",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_extrude.html#ae755c482803a84d6913755d8a754a1ce",
"class_f_e_a_t_1_1_geometry_1_1_attribute_set.html#aa2d5d001ffa3eb39e10d36390f7232eb",
"class_f_e_a_t_1_1_geometry_1_1_distributed_umbrella_smoother.html#a52a9a7d42e9427ad362399036a372d6c",
"class_f_e_a_t_1_1_geometry_1_1_index_set.html#a33220cbd7f779c4afef88ddc6bb44fd6",
"class_f_e_a_t_1_1_geometry_1_1_masked_boundary_factory.html#a66452a19a648e5b10e94b2ad58f199f4",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node.html#a04fb633874e522f61988fd5861732daa",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part_node.html#ab80198336a632da79ca33a6a868635a7",
"class_f_e_a_t_1_1_geometry_1_1_parti_par_m_e_t_i_s.html#a42009bbed414e9b289818bcefd73c793",
"class_f_e_a_t_1_1_geometry_1_1_patch_mesh_part_splitter_3_01_geometry_1_1_conformal_mesh_3_01_shba0a87c4bb4c4d2c925ee5e6ab30bd83.html#afa8dbf0ff8f4a958e136aa1d4caa8ff0",
"class_f_e_a_t_1_1_geometry_1_1_shape_convert_factory.html#a47275a6de0f82edd19056874042c5cd8",
"class_f_e_a_t_1_1_geometry_1_1_structured_mesh.html#a6223dbd07c6841d3381693f0a95888a3",
"class_f_e_a_t_1_1_geometry_1_1_template_search_space.html#a4b63b7850ff783c9192a6431ad946e9f",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#a63e6a9dbc780ce63e3016c5c4c682163",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_matrix.html#a703724161bd18cd18e6c8e822699ea47",
"class_f_e_a_t_1_1_global_1_1_matrix.html#aafc7cee7eceb0956e1f29b2a6816d18c",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar5e6f49c15ccb43840cb894a3375f33c3.html#a425ab19b989dcd0832de43fa63f04c28",
"class_f_e_a_t_1_1_global_1_1_synch_scalar_ticket.html#aa5f2e6d6a455ae40dcb6910d64c487d2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#a750b1210fab448b6154bfb22fcb36d98",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a34097d4a4f075e6096f8c0ea53996cfe",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a8aef4998758d1f2304fd2acca1076265",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#a6f3b4a5abd557187e9af23019968e8b7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_none_filter.html#a1910da2f3826e6bd5537c23021eb2aa3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_d2.html#a74a819eb6b929a479eaef30f52d4a3a0",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#a0d16a14bda7fc7b740b87c6711866b7c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#ac0ee2821cef88fbdd6e1a7b6732b7f95",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_saddle_point_matrix.html#afd36b623b8348c867b73c23cc30639da",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a9ef67b4c24ac885b2876763175a5f555",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a9e02b51ad88f3c4dc81c6077108b3302",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#adfded131472fdef6233c2957e724cb9e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a262c21e61815e101720371833fde4539",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a10f550bc09571b4accdb8ce14530c6b1",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#adc3c50680c59df206bdec4057f914dff",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_filter.html#af89174272556cef2470125428f90ab53",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a55e68321bded65673581ade0db918fc6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a76bc69841d76f29d5462baed6aa5f097",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_default.html#a27f13ff9194d8de5cd9f4da05e38c6fb",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a87e2ea0763273b9e6a1c63030b564368",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#a275ac94fbaed9620654896bd8d0c4100",
"class_f_e_a_t_1_1_runtime.html#a831515a57ea208a209c3086dd2fe5621",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a73654a64fe9f26443c5c12b460960983",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_analytic_function_operator.html#a6e309d5709bd7ef8c8d2583330e884d7",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#aa86287c6d142b78ac494420116688d66",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_descent.html#a82061dfdc0313a9da082736c7c48559e",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_global_1_1_matrix_3_01_local_matrix___00_12c7f8d5d38cd40f20a5510074419a06.html#a28cfc6782adf797385af1c3b7fc4d75a",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a78f9d130dfce27ae46f7dc79045e0850",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#ad0e0aa027e766929e7e741f6384f887a",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#abfeec9e7a424dee82dccb22be9f43fd4",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#a06c481c67eb3c7d71e77137062633968",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a64b8b83c3cff4d963bea8b04d9f3800a",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a5b0b165b5aba580f1b68bd3c2a49933a",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a5a3a1e9b5bae10f47d0a764a6f99908e",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a69fab2263b0c35bf03aa77d95e3a5917",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a5e276e3ab4c84a834377997024bf746e",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#abb072260d9458116ae6aa240d5685fb3",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a24df490c6a8a8000cec5ada30ac0c920",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#afca304c1b40eb84c276dd0510308bb34",
"class_f_e_a_t_1_1_solver_1_1_polynomial_precond.html#a111e1611874c524c6f6d157f16e3b01c",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a0ab9bef4ea22df32065a43399a990561",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a6950d64c2b9db6d1b76eff702a2bffa1",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#abfcf290e6fcb166ed8869d38f7d972fc",
"class_f_e_a_t_1_1_solver_1_1_s_o_r_precond.html#ad5975b047133ffcf21e1c7d0a88843e3",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a2c590983a3dcf459985b357ca1399543",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#a37367c98f42d457dbfc34fb36bc8da7a",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#aa00114dcc63934436ab46e03a00a5823",
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_node_functional.html#abc7087cde691d01f94b0084c96966334",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_82866b3cc738de4ae4935135fac5e732.html#aa8684b1a76f905215d1f7de7b7ac527f",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_element.html#a468ff0d4f362e01570e0e465014afcc3",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___001b011599bed78529504db50cb49456d.html#ab2cc1ebfc0e21abdac8f91b793a856fe",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element.html#a6cceeb32a988c6eaff5dd4a1d69c81ca",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_node_functional_3_01_space___00_010_00_01_variant_ab4cc4159982446cf6be9b262da9c709.html#abc7087cde691d01f94b0084c96966334",
"class_f_e_a_t_1_1_space_1_1_element_base.html#ac0a978ce5dd25bd74b99be6479891052",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#ab28a5371dee9efc07044f69a4e944535",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#a3c7283f42e823561cd3cc5da70bb5108",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_element.html#af250c0bc81fe3e1e64cc45cd2bef2dc6",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html#a3c7283f42e823561cd3cc5da70bb5108",
"class_f_e_a_t_1_1_space_1_1_node_functional_base.html#aa6cb95c632710f9e50e675a723bc5532",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_element.html#ac2407fb73711543919c560d477f2fb60",
"class_f_e_a_t_1_1_space_1_1_standard_vector_eval_traits.html#a21daec59b1df06832a9a6140b4f6a74b",
"class_f_e_a_t_1_1_syntax_error.html",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#a247b01215190e0845912001afc3c3b6f",
"class_f_e_a_t_1_1_trafo_1_1_evaluator_base.html#acdb7bfb197ee64e4d83a5901a0aaa456",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#a3712e4ae60ddf267304ff196ce8de8e3",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#a04385656ab1f3e606550f78430ae4b23",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#a4d4f718b82c4369a953bf9b2f21564b7",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a97164a88685009b2dfcd657bc4e65679",
"class_f_e_a_t_1_1_xml_1_1_content_error.html#a7aee2ab09e1e1840e107a63579f9db63",
"cuda__util_8cu_source.html",
"functions_type_j.html",
"mesh_file_format.html#meshfile_chart_circle",
"namespace_f_e_a_t_1_1_geometry.html#a05aca5dec690a05c85a05dea815e8176",
"namespace_f_e_a_t_1_1_solver.html#a2f9eb79ed3fa1defdf955da765d6ce8e",
"namespace_f_e_a_t_1_1_voxel_assembly_1_1_kernel.html",
"rumpf__functional_8hpp_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#a8e78f41a6f18578d8ec5b75f65d82e9a",
"struct_f_e_a_t_1_1_geometry_1_1_index_tuple.html#a13824dd90d90a27b43b08eb9557d2500",
"struct_f_e_a_t_1_1_l_a_f_e_m_1_1_arch_1_1_min_abs_index.html",
"struct_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_dof_traits_3_01_shape_1_1_hypercube_3_01shape__dim___01_4_00_010_01_4.html#af309ee5d36291d86171a31cda6aad6d3",
"struct_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01644f37e7e95ab6c4e55eccb88d5d5c51.html",
"struct_f_e_a_t_1_1_type_1_1_traits.html#a2b68ae86ec8c47e87ff3ac12bdf149dc",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_assembly_mapping_data.html#ae661ed72bb123eaad0ed782d45740a4d"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';