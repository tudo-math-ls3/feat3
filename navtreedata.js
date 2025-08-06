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
"class_f_e_a_t_1_1_assembly_1_1_discrete_function_integral_job_1_1_task.html#a724793c484805faa2d129713b570e97a",
"class_f_e_a_t_1_1_assembly_1_1_error_function_integral_job_1_1_task.html#ae7808db8cc580311ba7ee949cab0feb6",
"class_f_e_a_t_1_1_assembly_1_1_linear_functional_assembly_job_1_1_task.html#a2125f4b5965dc0a444114322344463c9",
"class_f_e_a_t_1_1_assembly_1_1_stokes_f_b_m_assembler.html#aeda51b65644952ecab2fd9e411efb8f1",
"class_f_e_a_t_1_1_control_1_1_blocked_combined_system_level.html#aa0fea8e531b5841b4281c802e991c156",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#aa6418bba40405c43dea84f0c99a04a4f",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#a376650a155aae2fc37f254dc30d5bf38",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#acc4f1f332f2c38f05bea2e9de64fa850",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_system_level.html#abc7e5f93ffcd2477664ec64fdce921b8",
"class_f_e_a_t_1_1_control_1_1_scalar_basic_system_level.html#ac1fea29a946fa5b9ea3b60bee1a16c9c",
"class_f_e_a_t_1_1_control_1_1_stokes_blocked_system_level.html",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_driver_base.html#ad3e90f9feb48734b4c4ab797eb3abd45",
"class_f_e_a_t_1_1_dist_1_1_comm.html#a071bffc0b9c220e8c4d9d0b8bc18f51d",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#add2683f64364cab03fffef5b5d2bae4a",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_mesh.html#a6eb78cf6dc7317003d1d66fce696f210",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a8b234feb91c3fc1e4b24b8fbc432ecc5",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#a2b037a7a14247a41a9b68ba4b2168951",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#ad72b359245ac85712bd6a6a411b888ba",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#a2f3158ab506ca75ae1bec39ede7e89a1",
"class_f_e_a_t_1_1_geometry_1_1_factory_3_01_structured_mesh_3_01shape__dim___00_01num__coords___00_01_coord___01_4_01_4.html#a49f01a67b8bc3ec7d5751cf20dc10104",
"class_f_e_a_t_1_1_geometry_1_1_intern_1_1_mesh_layer_3_01_template_set___00_01_mesh_shape___00_01_vertex_type___00_010_01_4.html#a68f267a37ed51131dde3923a40df8b1e",
"class_f_e_a_t_1_1_geometry_1_1_mesh_file_reader.html#a7151af9ea20b00d88dd8d662e6517366",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part.html#aca67498132474564bfd9e1a40dca7f71",
"class_f_e_a_t_1_1_geometry_1_1_parti_iterative.html",
"class_f_e_a_t_1_1_geometry_1_1_patch_mesh_part_splitter_3_01_geometry_1_1_conformal_mesh_3_01_shba0a87c4bb4c4d2c925ee5e6ab30bd83.html#a0de827a32394b0b9fa158ce3c24eaedc",
"class_f_e_a_t_1_1_geometry_1_1_schneiders_data.html#a31f967c44567ba29e594637112ccb3f6",
"class_f_e_a_t_1_1_geometry_1_1_struct_index_set_1_1_index_tuple_type.html#a3c3d3c0724e330b84eed80e236d8def9",
"class_f_e_a_t_1_1_geometry_1_1_template_builder.html#a7de96a5ba69c9e4f7194ebd7e5f6e743",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#a4d9ad15b11066968a42745e775d5b6b3",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_matrix.html#a3b25821f54a168a6aca6ebd797027be3",
"class_f_e_a_t_1_1_global_1_1_matrix.html#a855fb33f78eb70bbb127a8d4ae639cdf",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar5e6f49c15ccb43840cb894a3375f33c3.html#a3708b609a7ec545f104ba387a2bacd67",
"class_f_e_a_t_1_1_global_1_1_synch_scalar_ticket.html#a4db4ca744ac008b986c0ecc60db3ac47",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#a5839f37f6c1a3f41f714b68fd1514659",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a23c27253672b9a3900681acf8b6ce090",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a8908f5beea3784008d7d41d7bbebe815",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#a6adc39c851a24c861ab53b7037ec8728",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_none_filter.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_d2.html#a41d8fba193c740ec7786257b797b74a2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#a0f705074f86604abf1fc6955cc8d2850",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#accb79a551092f52fbc5613bd72a44744",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_serial_config.html#a4899ed29ca313b2aa863342ef61b827d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#aa1eb8ec95680a75b4572adf3165525f4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#aabfe5da81fcd1564e731e36a576d481d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#afe13d80e7fde8dd9296455addf1b6a4c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a34c1f6862697f1377d4d59ad6c3e50a6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a47b718987b62807efb0461015eb18b56",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_test_matrix_factory.html#a015ef8ebd71416585d6c5dae8801cd05",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix.html#a3947607b17928861b1ed6aa7ef6046fa",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a93ec9424f64afc350b9163f16892c5b3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a91980ce64f289928b47776a351d82244",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_pow_of_dist.html#a72c2c749c59be63b5098a56249b7367c",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a96989f24f88a4d725d78f26570136ea4",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#a7fd670002dbdc6aa041caf1b3bef984d",
"class_f_e_a_t_1_1_simple_arg_parser.html#a2f3fe2c2f05b25e13544a7d4814ac1e1",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a973d3942df838d68f63666265e73b025",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#aa2b2f7492bf2c6b7d9f740b95e99969f",
"class_f_e_a_t_1_1_solver_1_1_approximate_hessian_precond.html",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#ad355136ba5cd50b1f2ece6658ab1c0c8",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#ad355136ba5cd50b1f2ece6658ab1c0c8",
"class_f_e_a_t_1_1_solver_1_1_descent.html#ab562b0e476df821799ffd4151bf146bb",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_global_1_1_matrix_3_01_local_matrix___00_12c7f8d5d38cd40f20a5510074419a06.html#aeae87053092ba231f4af3a38131c539b",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a7647edeb0b1945e4a74e8d46e95a9c4c",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#aa2b2f7492bf2c6b7d9f740b95e99969f",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#ad4ef2095e84619a29292a37f078a877b",
"class_f_e_a_t_1_1_solver_1_1_hessian_precond.html#a0e97356074a327f5495808cbe24e43ce",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#aea4133c0018531fdfb8f23c90c97667f",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#aec76c75b2d2582ef2cf0ed928b583314",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#aa7cf18c4f6e62acfd0c3cf3285c00b96",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#a66c1cfdd52da05579dd92bec255803c2",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a77f32fe5475099d5ae331eb7b18d9280",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a7d583ff48c7678f1cffa30fdd48315a7",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a76a7e3315fcb44c2a3ac466ae1ee1337",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a7fc0cca82272a85f64e9a2110366f6fe",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a787225aa0743f3e1edc22f1dc39eebbe",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#ae6fc9035510853c34e062754732c17ba",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a5079e20c7708d1409a0798b7a226e9de",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#ab562b0e476df821799ffd4151bf146bb",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html#a6b319348e75d0163e9273ef052ae2d37",
"class_f_e_a_t_1_1_solver_1_1_polynomial_precond.html#a9e5d0c401afaf45dce4308642894a1e0",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a2ae68da1a782ff371c7288c3ae9a2922",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a84aa5df455caeb5545849abc0bb36fcb",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#ae314aed268f854678d35272cd5e41340",
"class_f_e_a_t_1_1_solver_1_1_s_o_r_precond_with_backend_3_01_preferred_backend_1_1generic_00_01_b057fb2a8d1f4c64f3e3aed2663fd10f.html#a0dde8161d8fcf1f885aade5892f786b5",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a4fdec11243202ce55a49e67d259cab2a",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#ab52cfa7a8be01f137c4366528ce7cbfa",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#aec7767fbc67979953cd9cc0f61efd015",
"class_f_e_a_t_1_1_space_1_1_basis_data.html#a5f6bb7d6750afe4064df0b9d32ee90aa",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_node_functional.html#aa6cb95c632710f9e50e675a723bc5532",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_element.html#af250c0bc81fe3e1e64cc45cd2bef2dc6",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___03bb1a3d92225206453ee71c00db5aa7a.html",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element.html#ae7c5f6ff90f5bbd5746bc18e7169ed82",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_null.html#a6707ee4f731a79db45187493aa48ab23",
"class_f_e_a_t_1_1_space_1_1_eval_data_reduced.html#a394d63a90ff721b7b41b93c9cfb65343",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spd08b4bf16ff97c315f2f45da544757eb.html#a4b1ac4e614fad291b9dd1b93491c8342",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#ad5886f4911cac95096215def29085034",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_evaluator.html#a1f5b8556e9fa1ee98993eeebeb5337c3",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html#ae4d17e132c111387fbe169819fc1c2aa",
"class_f_e_a_t_1_1_space_1_1_node_functional_null.html#af3b7b7f553156a6574c71c886f760a7c",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_element_1_1_node_functional.html#aae17cbffdd31a51654fa388d051b074d",
"class_f_e_a_t_1_1_statistics.html#a2de96f5e1bb466d92fd7c8bbeab409b5",
"class_f_e_a_t_1_1_test_system_1_1_test_list.html#aa0104f23b5e79c31dc56972a44dfccb4",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#aa3142caf7396109cdd211a1eb164fb56",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping.html#ab5798a3ebcbf4e11bc60ed65f08d011b",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#ac44d2d37d62bb78ca2bfabf63a8eaac9",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#a8d7e27c639581fef4896e6704d403162",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#aea22a7d5ba302d8e09c0df5fe18d05bf",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_mapping.html#a44f1f4be1cfbfa5a2aea869ae489fa2c",
"class_f_e_a_t_1_1_xml_1_1_markup_parser.html#ae38534833f6bda453dce182e4441c14a",
"dudv__functional_8hpp_source.html",
"isoparam_2mapping_8hpp_source.html",
"namespace_f_e_a_t.html#abc175499504e9be922533b3a3627ccaca0b731d2ebd143daed5111a68cffef1b3",
"namespace_f_e_a_t_1_1_math.html#a2af61b8fc64b284759c9f400c552a42a",
"namespace_f_e_a_t_1_1_solver.html#ac6a1659b386401b23df34866c47bb406",
"preproc_macros.html",
"struct_f_e_a_t_1_1_assembly_1_1_analytic_function_integral.html",
"struct_f_e_a_t_1_1_dist_1_1_operation.html#aa1634df4b3c359a1f09cb3c755bec538",
"struct_f_e_a_t_1_1_geometry_1_1_refinement_template_3_01_shape___00_01_shape_1_1_vertex_01_4.html#af3467d088ddc0d86793f61aab14fff11",
"struct_f_e_a_t_1_1_shape_1_1_simplex.html#a1a34154d62efd9f36488778493c4bcdc",
"struct_f_e_a_t_1_1_space_1_1_lagrange3_1_1_dof_traits_3_01_shape_1_1_hypercube_3_01shape__dim___01_4_00_010_01_4.html#aba1b9aa65fb05d33f552aecab132c1ea",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___00d0928204f7a6d17533bb19e7a72add33.html#a869e2b043c8a5d8a89023ec5ebe152ad",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01int_01_4.html",
"unit__filter__blocked__generic-eickt_8cpp_source.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';