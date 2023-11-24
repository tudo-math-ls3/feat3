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
        ] ]
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
"class_f_e_a_t_1_1_adjacency_1_1_dynamic_graph.html#a01dfbc54d76dc15f98afeec6da04297f",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_constant_function_1_1_evaluator.html#ae4df19352cefc54ac4335428ac187eda",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol_rhs.html#a784f02c7ddeaad33949653b892db0c31",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_rigid_body_vortex_pres2_d_1_1_evaluator.html#ae4df19352cefc54ac4335428ac187eda",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_tensor_static.html#a7d1aa3bc7bcd5156fafbde65f080c469",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_y_z_plane_parabolic.html#a64fe8c91a3a06815d0a9ece1385aa699",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#a5759200287edd9018266ee196896bf79",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a2f395c33a4e56dfe1e3a490033642dd6",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits2.html#a0cf5f2ce479e80d16d8135798c6b43bf",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_1_1_evaluator.html#a11beb9bb5d154eb3761026e3f3dd0925",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#adae2e75ea23d8fce5804558721dd7465",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a292a09f05426e7838e5b816b616fb5e3",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a292a09f05426e7838e5b816b616fb5e3",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job.html#aabe050d720f74771f8ea87ae63bc28ef",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_force_functional_1_1_evaluator.html#a8443c5f1bcdaf7d6720f83383202eab8",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_operator_blocked_1_1_evaluator.html#ad4e31ab56875fc586a66d64769492f99",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a7be4351570f5479fc5122163cc687f70",
"class_f_e_a_t_1_1_assembly_1_1_function_integral_info.html#a1d14af139749141a2d884e80438c89d6",
"class_f_e_a_t_1_1_assembly_1_1_scalar_discrete_eval_data.html#a24abcac0245d84346918b34ee15c1f3f",
"class_f_e_a_t_1_1_assembly_1_1_vector_discrete_eval_data.html#a993c2db681fb8332b8febbbdb9a934d9",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#a718d6a41ca76f3c54d79939d3a770bb8",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_quadratic_system_level.html#acc909b2745fa2a31aac81c92c1125d72",
"class_f_e_a_t_1_1_cubature_1_1_hammer_stroud_d2_driver.html",
"class_f_e_a_t_1_1_cubature_1_1_silvester_open_driver_3_01_shape_1_1_simplex_3_012_01_4_01_4.html#a4c21ca8ddebd3106f6e702a1b54097e5",
"class_f_e_a_t_1_1_dist_1_1_request.html#af0ec8e9a9ee9da08508002fb34438372",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a5583b8875f60874bbf5f5338fff31162",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#a61f6594a72b9cba3b1dd03074ba8648c",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere_chart_parser.html",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#aa5fc0dae89f687d0ea8ebd36cfa1b1fd",
"class_f_e_a_t_1_1_geometry_1_1_index_set.html#a253b2500e18e24d0b7929565c9807ee6",
"class_f_e_a_t_1_1_geometry_1_1_mesh_file_reader.html#a7304e3c465c3747f4aefbc90732ef830",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part.html#ae75ff5e5a33f114cbb9d26f958722b21",
"class_f_e_a_t_1_1_geometry_1_1_parti_par_m_e_t_i_s.html#a00439f21cf2b4a1b3e2c0b04b3720d6a",
"class_f_e_a_t_1_1_geometry_1_1_patch_part_map.html#a098697351463d6ca843ba3e52c97723d",
"class_f_e_a_t_1_1_geometry_1_1_standard_refinery_3_01_conformal_mesh_3_01_shape___00_01num__coor2b7a1f41c741938e6ddb9e7bf2c48774.html#aa5c40c48353cb8261e16543f3cf10093",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_3_01_l_a_f_e_m_1_1_dense_vector_3_01_d_t___00_01_i_t_914574362a1d637068ac8fa9eb1a3de4.html#a3a938907e73499f340668807dcc3a206",
"class_f_e_a_t_1_1_global_1_1_gate.html#ad617af511bed171cf241422788233de5",
"class_f_e_a_t_1_1_global_1_1_nonlinear_functional.html#ad0c92c6f3d34f120f9c9a36a18fd8245",
"class_f_e_a_t_1_1_global_1_1_synch_matrix.html#a22badc8ad48aa1faa5d47384eff2a434",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#a1005b3492abf7a4a246abdf7fee5fe89",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a1bfb1a46a75ac2c11fd6a6b16d157abb",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a9072241eb92f6dcb09898a4ae6f6edd0",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#a86d0a968a4baec5e597e5fd587439c9b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#a307326ef092b6d4cd0299a6ed594a0cf",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html#a8994516a57c288a94e88ce6a38dfd32f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_full_matrix.html#abee856642c868790e8d576c6be58ea8d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#a957f45e1ac77e72892fdd51440cae9cc",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a17f68b60689b37ca1eaeea7149a6dca7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a319f4512ada6844ad97fdd8dd89d07e7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a70e31f4e5cf0ae6af1ff63bbc235287b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#ac33488b8572ff5da8765f0e0f4eaaf0a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#ae8a894f5239d92a6006a1011ad5564d6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#a77a6db9e56999e4cd7bd665806eb5e36",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_filter.html#a5e1ed3e338abe987f38d36f8625b6302",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a381a42713a363bdbd2b1a30a3c286f2e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a62d1fc8caae749e4d4fbdde9fd2d66af",
"class_f_e_a_t_1_1_meshopt_1_1_chart_distance_function.html#af78222a1deae40d8c6fd0bf49dac4aeb",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a7c6df42b27f15609bdd86dffedfcf1d7",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function_base.html#af6374ccb2d0dc5474c77b2eb5fbcf7e6",
"class_f_e_a_t_1_1_runtime.html",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a76a7e3315fcb44c2a3ac466ae1ee1337",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a82061dfdc0313a9da082736c7c48559e",
"class_f_e_a_t_1_1_solver_1_1_analytic_function_operator.html#a8008cea5666559c4206db819e44ebb00",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#af8e655f0cdc31c3b245587b58a14d5b3",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#a8c34d07a55cd09eaafef49269f140536",
"class_f_e_a_t_1_1_solver_1_1_convert_precond_3_01_global_1_1_vector_3_01_local_vector_outer___00120b74f43424364c71ee3ecd9deab297.html#a25fafb09e4368a228e2c225048f8b94d",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_global_1_1_matrix_3_01_local_matrix___00_12c7f8d5d38cd40f20a5510074419a06.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a8965935265c8ba7609b4461b0a5672e0",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#aa95d79033c0912710a988a61e7312a15",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#aaf199e0c5145dc34413628edc7b94448",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#aa7cf18c4f6e62acfd0c3cf3285c00b96",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#aa0702a006fe8456127b0bb2ea4f013bb",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#aa0702a006fe8456127b0bb2ea4f013bb",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#a2e20a41bbb2460128529b4baf2450afd",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a73654a64fe9f26443c5c12b460960983",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a731e6e902833e07bd086615675b51c31",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a6fd0f61a9458e3b6c32b4c502fcdc7f6",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a76a7e3315fcb44c2a3ac466ae1ee1337",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a73654a64fe9f26443c5c12b460960983",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#ad5a06be09b1b40ab2c5b7cceb8274ab8",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a4832c9683db7a56c1dc72479e852c5f3",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#ab2e6ef76d11d3a708a85cf56d85d9caf",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html#a61697a335687bd2b7f1c3df5563f342d",
"class_f_e_a_t_1_1_solver_1_1_polynomial_precond.html#a85d4787bf4b51cc165f524ace74a7a89",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a288adb7c820ebbb82a0a8b07f178db13",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#ae0ecc3adc98c16bcfe9819a79d669af8",
"class_f_e_a_t_1_1_solver_1_1_s_o_r_precond_with_backend_3_01_preferred_backend_1_1generic_00_01_b057fb2a8d1f4c64f3e3aed2663fd10f.html",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a4fdec11243202ce55a49e67d259cab2a",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#aab005281520844d1c2e4bac4e0354aff",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#aeab098abe3b4dab611831f8a805d212c",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_045e90cd69ef865eaf5b37b181a99332.html#a5be99000b3ed16b5234932176f74d18f",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_evaluator_3_01_space___00_01_trafo_evaluator__802e55c907078015a390ad2e449f9e5b.html#a29bafa5636409dda3697334dc98bdf3c",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___03bb1a3d92225206453ee71c00db5aa7a.html#a24e648bd29f80680fe2028b3fe286673",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element.html#af50c6fb029acb0ad8f20f271c474c4b1",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_null.html#a96876c09d6b238bbf38ff30c82fcbd9b",
"class_f_e_a_t_1_1_space_1_1_evaluator_base.html#ab70b19b6962179079f8a2d38ddddcd28",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spd08b4bf16ff97c315f2f45da544757eb.html#ab9c09073e90efa2d66505fbf8f9b2062",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sbede9c86e9ef3e51ccbac5494ff49399.html#a3c7283f42e823561cd3cc5da70bb5108",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s535ddd64c2ab8dfee2b86e813f56bd8b.html#a6db7e5c3954c49646d2823108329675d",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_node_functional.html#aba5b55179b16f3517fee72d4ce4e0540",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s5e9945e9db3d72bec06849dc7e957367.html#ae4d17e132c111387fbe169819fc1c2aa",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_evaluator.html",
"class_f_e_a_t_1_1_statistics.html#ac09354c26c685beb919b5b3fc76a93bf",
"class_f_e_a_t_1_1_time_stamp.html",
"class_f_e_a_t_1_1_tiny_1_1_vector.html#ac0dc31ec783bdb5de112f748cf94223f",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_011_01_4_01_4.html#abf4777fd65acf5a08c20e0be3b47e27b",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a3682d5e3aa23c276cfba8e6c45d49d83",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_6f4e5b7f5f9d6e2d04743ee3976d1657.html#af2640bae0e23758fadff1879857857f9",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_011_01_4_01_4.html#a2adc79d0234d4cd7d4f32cc45b55d7f4",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#a8d2b2b8130db810ab9561bf348cb0cc7",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#aeb6960cda56f5b7ed5c379496a8e23ec",
"coarse__fine__cell__mapping_8hpp_source.html",
"functions_m.html",
"mesh_file_format.html#meshfile_chart_sphere_examples",
"namespace_f_e_a_t_1_1_l_a_f_e_m.html#ad372dae2e320ef82768327bb8342786fa21d01af167e6874b88ed5f6dc7c4b1e4",
"namespace_f_e_a_t_1_1_solver.html#aace7dd0b283bfeac3bba4392e013ed09",
"preproc_macros.html#ppm_feat_compiler_open64",
"struct_f_e_a_t_1_1_assembly_1_1_f_e_interpolator_3_01_space___00_01_space___01_4.html#ac61c5c93c112d6f1723839eaefd46801",
"struct_f_e_a_t_1_1_geometry_1_1_mesh_part_1_1_target_set.html",
"struct_f_e_a_t_1_1_space_1_1_bernstein2_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercuf2dc830d9bb054392783ac3d465519f0.html",
"struct_f_e_a_t_1_1_type_1_1_traits.html#a43b0fe93cf62e1a1c689e07ac29c877b",
"test__system_8hpp.html#a7585e793a391d301bd4cbdc66f024018"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';