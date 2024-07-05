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
"class_f_e_a_t_1_1_adjacency_1_1_cuthill_mc_kee.html#a9af7891d15ae84ebd0a71930d9a4b449",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_ball_cap_function2_d_1_1_evaluator.html#a1eb375034ddde8778d98406e71725128",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol_rhs_1_1_evaluator.html#a1eb375034ddde8778d98406e71725128",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sin_y_t0.html#acc011af950d314ec9e2f8b16ed52feca",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sphere_sin_bubble_function_1_1_evaluator.html#ae4df19352cefc54ac4335428ac187eda",
"class_f_e_a_t_1_1_analytic_1_1_curl_1_1_evaluator.html#aed244fc095f976a1431e02b2b00ef4b4",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function2_d.html#ab2352cd43ccbb40e7772208a857223e9",
"class_f_e_a_t_1_1_analytic_1_1_parsed_scalar_function.html#a50702e8df25f9028b46d8453bb41dba9",
"class_f_e_a_t_1_1_assembly_1_1_analytic_function_integral_job_1_1_task.html#a8ba97748f6da9e37a06a0702f6deef51",
"class_f_e_a_t_1_1_assembly_1_1_basic_matrix_assembly_task_c_r_t_p2.html#af1091d97ff4a02475c1ac76c7147c480",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembler.html#afa9730fcf687d0c4e012dbabd33e7fb8",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_assembly_task_base.html#af80337cfa1a32793a469d2656f292426",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_vector_assembly_job_1_1_task.html#a9992621baf02dd812a9f30ce1d60a307",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_matrix_assembly_job_1_1_task.html#ab0d4b588e5198fd8569b812139cdfdb5",
"class_f_e_a_t_1_1_assembly_1_1_cell_error_function_integral_job_1_1_task.html#afcf3f8ad8b48a3e423a7d86e81d5c37a",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_identity_operator_blocked_1_1_evaluator.html#ab5e794de6dd40e01fbce22aa7084435e",
"class_f_e_a_t_1_1_assembly_1_1_discrete_cell_projector.html",
"class_f_e_a_t_1_1_assembly_1_1_error_function_integral_job.html#a68690391fdfb74f26c5ce971fe1cb7ee",
"class_f_e_a_t_1_1_assembly_1_1_linear_functional_1_1_evaluator.html#a63997f4f05078c674fe46aa007ae5de8",
"class_f_e_a_t_1_1_assembly_1_1_stokes_f_b_m_assembler.html#a838d858d52d5abb819c7b6d0c552c021",
"class_f_e_a_t_1_1_control_1_1_checkpoint_control.html#a24cc743f59c1c34ba51ae326c2e241b6",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#a81b90dfc6cf05e026a5a2274704376cd",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_system_level.html#ae8147c5c2addd69b35bcc91e2d3c8775",
"class_f_e_a_t_1_1_control_1_1_time_1_1_nvs_bdf_q.html#aa92f4ae9caa234225fd5a38e86ae5137",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_maclaurin_driver.html#af06dfb81097f5ac0f913f74d040c0b1b",
"class_f_e_a_t_1_1_dist_1_1_comm.html#a92a62cb6f962c57bb8581f05ce82a02e",
"class_f_e_a_t_1_1_dist_file_i_o.html#ab2b09f77756278daae19a8f82640b961",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_c_g_a_l_surface_mesh.html#ac6034c009a5e1947156609746276d976",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_extrude.html#a82d358547b19b834dfec22a2e4eb89fc",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_surface_mesh.html#ae52ee9f971dd82d8196c0a2e21644548",
"class_f_e_a_t_1_1_geometry_1_1_export_e_p_s.html#a6fb066fe520f06a190563dbb736c2dbb",
"class_f_e_a_t_1_1_geometry_1_1_intern_1_1_parti_iterative_individual.html",
"class_f_e_a_t_1_1_geometry_1_1_mesh_file_writer.html#a70c15170179ceb543be5d92eaa65e4e0",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part_node.html#a9f4626787d7443c58f2248869920db84",
"class_f_e_a_t_1_1_geometry_1_1_parti_zoltan.html#ab9423ccb8908d0a47e0bb26519dc675a",
"class_f_e_a_t_1_1_geometry_1_1_reference_cell_factory.html#a2306e42600b5d6425e8efdac123cca22",
"class_f_e_a_t_1_1_geometry_1_1_standard_refinery_3_01_mesh_part_3_01_parent_mesh___01_4_01_4.html#afaa4e9cc12bc2e89aaab61a0d14ccd40",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_3_01_l_a_f_e_m_1_1_dense_vector_3_01_d_t___00_01_i_t_914574362a1d637068ac8fa9eb1a3de4.html#aed21eec0a6a21c299c7acb3b6af2cfce",
"class_f_e_a_t_1_1_global_1_1_matrix.html#a536d5ef27a9bfbae92e03ab9f267838f",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar4d3a70ff49e16f0d9ba9d236f0dafdc0.html#a3e43f3bc5757fad55572501af362cbdd",
"class_f_e_a_t_1_1_global_1_1_synch_scalar_ticket.html#aa5f2e6d6a455ae40dcb6910d64c487d2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#a9f5658b930db28e77a5f75e1b4e48879",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a72557b404ef0202f93f901d6faf30b0f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#ade7d7a8fe820d08f0c855a8e81f5fb53",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#adeb4b7feed0f853a4accd76569a914dd",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#aa3256534af8268ec485d428eab08e845",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_diag_matrix.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_mirror.html#a689b6113c73b6043afd45fb924452530",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_saddle_point_matrix.html#a07c29a0f928575eeef6c66540076ed32",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a438fce4e323c32f6053360e43d476ece",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a62d1fc8caae749e4d4fbdde9fd2d66af",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#aa2f8f09ef17204825e8c8b53a1165fed",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a024cc69b9a5b9a820f723497db88ae6d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_factory.html#a17c52b5d04e28206a31dda9587b4873b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#aa0f8a8c78b7e7b226dd476bd9fe924d3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix.html#a1e5a1800d55a6929bbeb82f36d504b8f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a903b2e5e9655cbdc23d4453baa72f995",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a91980ce64f289928b47776a351d82244",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_pow_of_dist.html#ae1d6fb2da827adca004079dd3834f7c9",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#aa07a9b1103093ea07f2a204b7869d722",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#a9ce7bb0307f28b70dddba3553133abd4",
"class_f_e_a_t_1_1_simple_arg_parser.html#a99f09de8dec9e2cd253508bf975459b8",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#aa7281f8e3498245d411c5829b9d0cec0",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_approximate_hessian_precond.html#ac687f31fa9389b4c0c0ec211bc262e0c",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a27d93b800c82dae6cf39f5ffc9a94015",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#a173c0e90d872337278ef536f07fb1ac9",
"class_f_e_a_t_1_1_solver_1_1_cu_solver_l_u.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_maf2bb55a50c59e69123dfd500b5d8c62b.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#ac1b7fdb282872c9e1ef6598f5175d3d5",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#ad0632c8f5c807745f061f9afdda0f83b",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#af5be9eafbcdd5d8d6b43d4781698f198",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#add5c03ded9d5568e14ec60cdfb1572b8",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#add5c03ded9d5568e14ec60cdfb1572b8",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#abfeec9e7a424dee82dccb22be9f43fd4",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#abb072260d9458116ae6aa240d5685fb3",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#adada1c798d6725f2dd0f75acd5822671",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a9920812266b8f49921c9b133dd1451b4",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#aa42f0f9d64cafd591b77644cec9dab53",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a9bfcb3b1ddc260d5c5017f1c729bd85b",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#aa5944f939b4e98706eff002e1c4caa3e",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#aa86287c6d142b78ac494420116688d66",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#a0810de389dc450bea5cb6917b66ca051",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a7647edeb0b1945e4a74e8d46e95a9c4c",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#ae6fc9035510853c34e062754732c17ba",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html#adb7e9889ab00a5710a7e359bc3f82754",
"class_f_e_a_t_1_1_solver_1_1_precon_wrapper.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#ab48c468782b6695d34bee1089a88b785",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#a1feb792955098cf117423f18dfd0db8f",
"class_f_e_a_t_1_1_solver_1_1_s_s_o_r_precond_base.html",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a73654a64fe9f26443c5c12b460960983",
"class_f_e_a_t_1_1_solver_1_1_umfpack.html#a77fb10f30c5b4eaa9aed7f6cfb8b07a2",
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_element.html#a8d82bce6274b32899746656ddb66160c",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_64897ed7485eba824a67fb42fd5ab9c0.html#a2644f5cb97bfdc7ce8f437a8c9176b4b",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_evaluator_3_01_space___00_01_trafo_evaluator__802e55c907078015a390ad2e449f9e5b.html#adc9c506c493349492a7576a9c19fc304",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___03bb1a3d92225206453ee71c00db5aa7a.html#ae4d17e132c111387fbe169819fc1c2aa",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_964e26ad63d9177b3f747659af96be6e.html#a928c1fc25e4c56643b88d1a7510cb40f",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_uniform.html#a8fa6fdc0825b0c04c099ff2150a58724",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_element.html#a1002e67d8130c5f70dfa107f11cfc904",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_node_functional.html#ad091f3ff7d0c9c412833b3a8446a43ff",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sbede9c86e9ef3e51ccbac5494ff49399.html#af4d8faee5ca28cbbdabc180824d218b5",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_node_functional.html",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s5e9945e9db3d72bec06849dc7e957367.html#ad0065e2fbb5b4ca9de85db982abbfcbc",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_element_1_1_evaluator.html#a4fe957dedb7326b088c17b209864ac2e",
"class_f_e_a_t_1_1_statistics.html#ab307f8998aa50e24bfbd462b95d5115d",
"class_f_e_a_t_1_1_thread_fence.html#ad161e6a2064c1db8a28cdedc1c4a3003",
"class_f_e_a_t_1_1_tiny_1_1_vector.html#a78722591e5b03ddd8abaf6b75a0a5318",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_011_01_4_01_4.html#a99aacf4b77b94a83e0558da0601a233f",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_6f4e5b7f5f9d6e2d04743ee3976d1657.html#ae06d7af49337ffa78c0f58646a59b918",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator.html#a17dbe86e0894d5638b0d8a07e7484c7e",
"class_f_e_a_t_1_1_voxel_assembly_1_1_voxel_defo_assembler.html",
"compiler__microsoft_8hpp_source.html",
"functions_o.html",
"mesh_file_format.html#meshfile_basic_encoding",
"namespace_f_e_a_t_1_1_l_a_f_e_m.html#a82e983d75b2911bc2cd014a541d50711a6adf97f83acf6453d4a6a4b1070f3754",
"namespace_f_e_a_t_1_1_solver.html#aa2e0c5a93d63bfa65939ba63207d6f96",
"pcg_8hpp_source.html",
"struct_f_e_a_t_1_1_analytic_1_1_image_1_1_vector.html#a044f79d19f8de09b180178319c081d2a",
"struct_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle_traits.html#a05d9ddccee41308f0988c9938e2ca8d1",
"struct_f_e_a_t_1_1_shape_1_1_face_traits.html#adbcd406bffc596476c67d38dbfb5a384",
"struct_f_e_a_t_1_1_space_1_1_lagrange3_1_1_dof_traits.html#a25a22a28d9328282e53c3588bf4a4a00",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___00f78f81b86ee7ca50fe24a6c9ff685ed0.html#aa8b9794c45b1f771c68fab8aebe08a96",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01long_01long_01_4.html#a1abba6f2325064272569bfe0898a81e4"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';