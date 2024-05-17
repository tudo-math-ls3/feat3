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
    [ "Algebraic Multigrid", "amg.html", [
      [ "AMG Introduction", "amg.html#amg_introduction", null ],
      [ "Integrating AMG in your application", "amg.html#amg_integration", null ]
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
    [ "Feature-requests and TO-DO-List", "feature_requests.html", [
      [ "Project coordinates from fine-grid onto the coarse-grid", "feature_requests.html#grid_coord_projection", null ],
      [ "Document the differences betweeen a type-0-Vector and a type-1-Vector. Also explain the effect of the sync_X-function on the shared", "feature_requests.html#documentation_type0_type1_sync", null ]
    ] ],
    [ "INI Data File Format", "ini_format.html", [
      [ "Comments", "ini_format.html#sec_comment", null ],
      [ "Properties", "ini_format.html#sec_property", null ],
      [ "Sections", "ini_format.html#sec_section", null ],
      [ "Examples", "ini_format.html#sec_example", null ]
    ] ],
    [ "Known Issues", "known_issues.html", null ],
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
    [ "<c>LIKWID</c> marker API with FEAT", "likwid_for_feat.html", [
      [ "LIKWID General Information", "likwid_for_feat.html#LIKWID_descr", null ],
      [ "Prerequisites for using LIWWID", "likwid_for_feat.html#LIKWID_prereques", null ],
      [ "Using the LIKWID marker API", "likwid_for_feat.html#LIKWID_use", null ]
    ] ],
    [ "FEAT Mesh File Format", "mesh_file_format.html", [
      [ "Basic Information", "mesh_file_format.html#meshfile_basic", [
        [ "Overall File Structure", "mesh_file_format.html#meshfile_file_structure", null ],
        [ "Terminology", "mesh_file_format.html#meshfile_basic_terminology", null ],
        [ "Encoding", "mesh_file_format.html#meshfile_basic_encoding", null ],
        [ "Parser limitations", "mesh_file_format.html#meshfile_basic_parser_info", null ]
      ] ],
      [ "Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering", [
        [ "Triangle Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_tria", null ],
        [ "Quadrilateral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_quad", null ],
        [ "Tetrahedral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_tetra", null ],
        [ "Hexahedral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_hexa", null ],
        [ "Mesh Type Specification", "mesh_file_format.html#meshfile_meshtype", null ]
      ] ],
      [ "The Root Node", "mesh_file_format.html#meshfile_root", [
        [ "Examples", "mesh_file_format.html#meshfile_root_examples", null ]
      ] ],
      [ "Chart Description", "mesh_file_format.html#meshfile_charts", [
        [ "2D Circle Chart Description", "mesh_file_format.html#meshfile_chart_circle", null ],
        [ "Circle Parameterization", "mesh_file_format.html#meshfile_chart_circle_param", [
          [ "2D Circle Chart Examples", "mesh_file_format.html#meshfile_chart_circle_examples", null ]
        ] ],
        [ "3D Sphere Chart Description", "mesh_file_format.html#meshfile_chart_sphere", [
          [ "3D Sphere Chart Examples", "mesh_file_format.html#meshfile_chart_sphere_examples", null ]
        ] ],
        [ "2D Bezier Chart Description", "mesh_file_format.html#meshfile_chart_bezier", null ],
        [ "Terminology", "mesh_file_format.html#meshfile_terminology", [
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
      [ "Root Mesh Description", "mesh_file_format.html#meshfile_mesh", [
        [ "Vertices Description", "mesh_file_format.html#meshfile_mesh_vertices", null ],
        [ "Topology Description", "mesh_file_format.html#meshfile_mesh_topology", null ],
        [ "RootMesh Examples", "mesh_file_format.html#meshfile_mesh_examples", [
          [ "2D Quadrilateral Unit-Square", "mesh_file_format.html#meshfile_mesh_examples_2d_unit_square_quad", null ],
          [ "2D Triangle Unit-Circle", "mesh_file_format.html#meshfile_mesh_examples_2d_unit_circle_tria", null ],
          [ "3D Hexahedral Unit-Cube", "mesh_file_format.html#meshfile_mesh_examples_3d_cube_hexa", null ]
        ] ]
      ] ],
      [ "Mesh-Part Description", "mesh_file_format.html#meshfile_meshpart", [
        [ "Mapping Description", "mesh_file_format.html#meshfile_meshpart_mapping", null ],
        [ "Topology Description", "mesh_file_format.html#meshfile_meshpart_topology", null ],
        [ "Attribute Description", "mesh_file_format.html#meshfile_meshpart_attribute", null ],
        [ "MeshPart Examples", "mesh_file_format.html#meshfile_meshpart_examples", [
          [ "2D Unit-Square Top Edge Mesh-Part", "mesh_file_format.html#meshfile_meshpart_example_unit_square_simple", null ],
          [ "2D Unit-Circle Boundary Mesh-Part", "mesh_file_format.html#meshfile_meshpart_example_unit_circle", null ]
        ] ]
      ] ],
      [ "Partition Description", "mesh_file_format.html#meshfile_partition", [
        [ "Patch Description", "mesh_file_format.html#meshfile_partition_patch", null ],
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
      [ "The mesh-extruder Tool", "tools_meshtools.html#meshtools_mesh_extruder", null ],
      [ "The mesh-indexer Tool", "tools_meshtools.html#meshtools_mesh_indexer", null ],
      [ "The mesh2vtk Tool", "tools_meshtools.html#meshtools_mesh2vtk", null ],
      [ "The mesh2eps Tool", "tools_meshtools.html#meshtools_mesh2eps", null ],
      [ "The mesh2tri Tool", "tools_meshtools.html#meshtools_mesh2tri", null ],
      [ "The tri2mesh Tool", "tools_meshtools.html#meshtools_tri2mesh", null ]
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
        [ "FEAT_COMPILER_INTEL", "preproc_macros.html#ppm_feat_compiler_intel", null ],
        [ "FEAT_COMPILER_MICROSOFT", "preproc_macros.html#ppm_feat_compiler_microsoft", null ],
        [ "FEAT_COMPILER_OPEN64", "preproc_macros.html#ppm_feat_compiler_open64", null ],
        [ "FEAT_COMPILER_ORACLE", "preproc_macros.html#ppm_feat_compiler_oracle", null ],
        [ "FEAT_COMPILER_PGI", "preproc_macros.html#ppm_feat_compiler_pgi", null ]
      ] ],
      [ "Backend / Library Macros", "preproc_macros.html#preproc_sec_backends", [
        [ "FEAT_HAVE_ALGLIB", "preproc_macros.html#ppm_feat_have_alglib", null ],
        [ "FEAT_HAVE_CGAL", "preproc_macros.html#ppm_feat_have_cgal", null ],
        [ "FEAT_HAVE_BOOST", "preproc_macros.html#ppm_feat_have_boost", null ],
        [ "FEAT_HAVE_DEATH_HANDLER", "preproc_macros.html#ppm_feat_have_deathhandler", null ],
        [ "FEAT_HAVE_CUDA", "preproc_macros.html#ppm_feat_have_cuda", null ],
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
        [ "FEAT_HAVE_ZLIB", "preproc_macros.html#ppm_feat_have_zlib", null ]
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
        [ "Nvidia CUDA v11.8", "feat_for_win_vs17.html#win_vs17_prereqs_cuda", null ],
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
        [ "Enumerator", "functions_eval.html", null ],
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
"class_f_e_a_t_1_1_adjacency_1_1_cuthill_mc_kee.html#a0b672fc2ceece924836f033bf5bffb52",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_ball_cap_function2_d_1_1_evaluator.html#ae4df19352cefc54ac4335428ac187eda",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol_pressure_1_1_evaluator.html#a1eb375034ddde8778d98406e71725128",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_rigid_body_vortex_pres2_d.html#ac87e338e0b9e2dfb144a9ea842172e8f",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_tensor_static.html",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_x_y_plane_rotation_1_1_evaluator.html#a9f9e678d53aae13b5d234fcfd43f1854",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function2_d_1_1_lambda_set.html#a6664577682bdf70d986631f7a071483f",
"class_f_e_a_t_1_1_analytic_1_1_parsed_scalar_function_1_1_evaluator.html#a1eb375034ddde8778d98406e71725128",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits1.html#a9aaa57c8cba82d0686ed0a694262a3d5",
"class_f_e_a_t_1_1_assembly_1_1_basic_vector_assembly_task_c_r_t_p.html#ab807fc92cc3192026f1eb0cf1962a626",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#a34f836bb6de7a3f2304de70d46516716",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job.html#aeeb2a0f5ed21f0ba6fccb54bdeba7459",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_vector_assembly_job_1_1_task.html#aeb81ae5049b5cd987c5782c183ea3dbd",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job.html#a3352b79041bcacc5b6c5ac8a04042799",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_force_functional.html",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_operator_blocked_1_1_evaluator.html#a11beb9bb5d154eb3761026e3f3dd0925",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a54ce4df37dfd4222766a741f5b849499",
"class_f_e_a_t_1_1_assembly_1_1_force_functional_assembly_job_1_1_task.html#acc03c5544df20a41f24be7496909c5d9",
"class_f_e_a_t_1_1_assembly_1_1_oldroyd_assembler.html#ad452b05fc0d9e06ab4207bdb5c49c5e6",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembler.html#abb817e8b27489789b529e66be7f1469f",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#af443adafe2aa4127d2f638c7a9a3b724",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#aca72a6f6fe43030410e410736a1a1c65",
"class_f_e_a_t_1_1_control_1_1_scalar_basic_system_level.html#a9f7285e2c8f77b7274757f478df39b0a",
"class_f_e_a_t_1_1_cubature_1_1_hammer_stroud_d3_driver.html",
"class_f_e_a_t_1_1_cubature_1_1_simplex_scalar_factory_3_01_scalar_driver___00_01false_01_4.html",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#a1c5efdb5b2546346edb00ae520a3c4b5",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a726de2b50f4a191e962d71aba69b2237",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_chart_c_r_t_p.html#afe05bab6765b43567adab8d0dda0abd4",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#a93ba9bfe598632e8c3d671140a882b59",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html",
"class_f_e_a_t_1_1_geometry_1_1_factory_3_01_structured_mesh_3_01shape__dim___00_01num__coords___00_01_coord___01_4_01_4.html#a4350419eea85983cc5239866f36af711",
"class_f_e_a_t_1_1_geometry_1_1_mesh_distortion.html#adc33d6694ec6457399f1e184900fad4f",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node_linker_error.html#a04373352091dcf206f09e86fcc9a910f",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#ade7bd2978b0005b7042040623f9cd399",
"class_f_e_a_t_1_1_geometry_1_1_patch_mesh_factory.html",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#af2367136c5b680432ce6646d88b0a002",
"class_f_e_a_t_1_1_geometry_1_1_structured_mesh.html#a5a5fcf762b9aa0b3417d97380228453b",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_vector.html#a8bb7a7bbb9a4fce788b8c9450034c0ab",
"class_f_e_a_t_1_1_global_1_1_muxer.html#a8f1b20e66a567ae850128287a4d98c90",
"class_f_e_a_t_1_1_global_1_1_splitter.html#a686df90e7739303e563bcf2e647dea77",
"class_f_e_a_t_1_1_global_1_1_vector.html#a2d611cb63a8d5690538a9c4e63236c87",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a927a42091d38f2db9be7944a97c7b2ae",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a06793709fe57f0a870e3e72564ea06c3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_sequence.html#adf6d857bb135564ff9271a29f2d542c5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_meta_matrix_test_base.html#af755ae56bc91c0ceb3894ca0910190f1",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_d2.html#a13eff2fd0eb68d2fbcc0e1a87cd604e2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#a0fd5e738cf0c80888eeb614b8368a3e7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#ad4638ed2020a114071b0593697245fe9",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_slip_filter.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#acbcdba2e6aab53431ecb3e0393fb87b8",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#ae3050afee64fc41ed0ca4d2845fbe026",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#a2a8542d823f9e80ee431e1cdcda6ea4c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a7915be74ff0aba8d97b530cf87c663b6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#ab36d5093e53ccf1cf5a652467de0e225",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_diag_matrix.html#a0aaacfab540d1111896983a53fe1b01c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix_row.html#aa5279ff55102d791ff12aa08a9d90804",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter_blocked.html#a318772496d3749014b2f4c3a19884ccc",
"class_f_e_a_t_1_1_meshopt_1_1_chart_distance_function.html#a00d8c49f24c67dedda8df3f218b5bcf0",
"class_f_e_a_t_1_1_meshopt_1_1_du_dv_functional.html#ad3ad8e77f124e776efbebfe6a9c0e9f9",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function.html#a8c19f208b86039dc4ea790d7def9f2ba",
"class_f_e_a_t_1_1_property_map.html#a132551dc2bd19a32c14ea0bda2bd2ce1",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a1d005a8fa11c118322cb35c7f1297155",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a25362c6f4bb0c5d1deac80d947442303",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka.html#a5ec627d1c43520c8eb7ed436c3fcf065",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#a7dd50fe57adbb7e2de822fa82f66afe1",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#abb072260d9458116ae6aa240d5685fb3",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#a6bc593354df3004878555aecfd280fad",
"class_f_e_a_t_1_1_solver_1_1_expression_restriction.html",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a441d5b2029da644fa997c2b506c2696c",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#a27d93b800c82dae6cf39f5ffc9a94015",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#a26458def0ef8654e3fff8f444df8108f",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a27d93b800c82dae6cf39f5ffc9a94015",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a40db3688d9794b0334a0226994b17c05",
"class_f_e_a_t_1_1_solver_1_1_multi_grid.html#a48cf3105183de2a90678b3d86387252a",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#afaab86a331f64d630f134f3c6027b5f7",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_precond.html#adb60a343f20ba4b9243a7df635657a1c",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_nonlinear_operator_precond_wrapper.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#a3bd60beb7925eaeada45028627f85281",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#a24df490c6a8a8000cec5ada30ac0c920",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#aa86287c6d142b78ac494420116688d66",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#aa86287c6d142b78ac494420116688d66",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#aef06db65ca494feb7fa2fd46459c8fe5",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#a50df7223074bef35a5a7cbd0c8ec9bec",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#ab2e6ef76d11d3a708a85cf56d85d9caf",
"class_f_e_a_t_1_1_solver_1_1_schwarz_precond_3_01_global_1_1_vector_3_01_local_vector___00_01_mic9e5dece21cb7735a0bc194fb3569eb4.html#a6e6988db4e7c2368e682888d9ca623f2",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#aeb397e3e98d6e2edb5e4ee5cb8f73082",
"class_f_e_a_t_1_1_solver_1_1_uzawa_precond.html#af0a4d795e681d81561539d11fbe8e272",
"class_f_e_a_t_1_1_space_1_1_basis_data.html#a5466cd997dfd93399cc0a8d5d0b40795",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_node_functional.html#ad091f3ff7d0c9c412833b3a8446a43ff",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_element_1_1_evaluator.html",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___0d5fd0fd65b9232ba00b88a0e0458fe94.html#a6db7e5c3954c49646d2823108329675d",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_d256e15bbc7b3e4885146b0f6ddde539.html#aa34421fd18b00301822b68b59073481f",
"class_f_e_a_t_1_1_space_1_1_dof_mapping_uniform.html",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#a05bec126e2b37bf1c13ca9c419cb4eeb",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s189823cc2b39e0d1a920e7698210faed.html#a6fa689aeb431d8551cf2977d30edd661",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_node_functional.html#ad091f3ff7d0c9c412833b3a8446a43ff",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s17de294e509ec5be8be361ba5ec7318e.html#a6bfdd412851516ebd9e42a508b51edb6",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_node_functional.html#aa6cb95c632710f9e50e675a723bc5532",
"class_f_e_a_t_1_1_space_1_1_parametric_evaluator.html#a8129c36880171d8d4d4d86ce6814fcc8",
"class_f_e_a_t_1_1_string.html#ac8f44f2b1d8f4c2952538e8639319585",
"class_f_e_a_t_1_1_tiny_1_1_matrix.html#acae5c4812be8730832e6f6597386e876",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping.html",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#a49d3cd9c5a6d5b8c88cd3a061de04970",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#a25dcb1fb84889eabfa3850c1a632faf1",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#a7e20e35337f3f0879235fe4c905dee0f",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#aaaa3c71f9b2bad6c11ff174de7c61602",
"class_f_e_a_t_1_1_xml_1_1_markup_parser.html#ae00f0ff351c6d59e2adfe1daca338dc3",
"empty__lafem_8cpp_source.html",
"kernel_2lafem_2base_8hpp.html#a82e983d75b2911bc2cd014a541d50711ad84a121691fde30486dfe366a7389d59",
"namespace_f_e_a_t.html#aeeeaf7afcb3a42f69b57031db04b2428aaaff6440f63c38bcfdc6816c4d85411b",
"namespace_f_e_a_t_1_1_pack.html#a7c69c0366e113fd42901ba9bcea55d4daacd20a348a5c2b33a193bd82b393b25c",
"namespace_f_e_a_t_1_1_voxel_assembly.html#a637d205e4e3d706dcc1d6a5fb619c0ed",
"runtime_8hpp_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#abd4d61f2de8adccfce479cb5c298d51a",
"struct_f_e_a_t_1_1_l_a_f_e_m_1_1_arch_1_1_intern_1_1_apply_banded_1_1_iteration___right.html",
"struct_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___f7f932e34cb694f42be1d9bf6f7db3ce.html#a47c8fc876e1ccc3341e634564dc195fa",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___0027d78bfa85bd09ff93bcee03217b623c.html#aac6d1844087c35cdcf45816843267a7e",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01signed_01char_01_4.html#aee9f2fe93b3d28b7a9e7406928838940",
"tetris__hexa_8hpp_source.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';