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
"class_f_e_a_t_1_1_adjacency_1_1_dynamic_graph.html",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_constant_function_1_1_evaluator.html#ad4a815291497cde8b9420ed474c64633",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol_rhs.html#a4d9f44a50686dd46d07f3be9a19123a7",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_rigid_body_vortex_pres2_d_1_1_evaluator.html#aa2f90323893b36e0bdb8e117c7c53b6b",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_tensor_static.html#a781c7c6a1e8840664be7d16198cf41af",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_y_z_plane_parabolic.html#a447397b3cd0d0271be9d658228215b07",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#a468df095b44ffe15cd65c2d31d487da2",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a14d9e06c7bb55cda7cfea1d2f1562fcd",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits2.html#a03c3ede75963845a3dd9702f62f5b224",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_1_1_evaluator.html",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#ad74b37015de3ba9088143fdb0c175a1a",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a25a4188596a0fd8552e17da078005e36",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a16cf785893294825c0949ea32b1ec969",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job.html#aa6598c1a041c7e2d45ec8962469cb4d2",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_force_functional_1_1_evaluator.html#a795f3f587dbb52cf614c83547337aff4",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_operator_blocked_1_1_evaluator.html#acd107b01fd1f388a76929a06fd0d1509",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a7bb2dc51f98f92bd55b741ed4026371b",
"class_f_e_a_t_1_1_assembly_1_1_function_integral_info.html#a196296f5ca6d8a87cab8556578b594d7",
"class_f_e_a_t_1_1_assembly_1_1_scalar_discrete_eval_data.html#a0e4976973cb6f2062cba86d5fc2e3efe",
"class_f_e_a_t_1_1_assembly_1_1_vector_discrete_eval_data.html#a931970b05ad9e4bd73e941b3a0dddc25",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#a6e66ba14a3b61ed7f6323c36770460fc",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_quadratic_system_level.html#acbc8e6087e323a35d52be35512d4365c",
"class_f_e_a_t_1_1_cubature_1_1_hammer_stroud_d2_driver.html",
"class_f_e_a_t_1_1_cubature_1_1_silvester_open_driver_3_01_shape_1_1_simplex_3_012_01_4_01_4.html#a4c21ca8ddebd3106f6e702a1b54097e5",
"class_f_e_a_t_1_1_dist_1_1_request.html#af0ec8e9a9ee9da08508002fb34438372",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a5583b8875f60874bbf5f5338fff31162",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_chart_c_r_t_p.html#ac21fc4262e3e1ba20276ab2631de6095",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#a82d358547b19b834dfec22a2e4eb89fc",
"class_f_e_a_t_1_1_geometry_1_1_chart_hit_test_factory.html#aa019f738bc569ba8abb67bb6398d2ef7",
"class_f_e_a_t_1_1_geometry_1_1_factory_3_01_mesh_part_3_01_mesh_type___01_4_01_4.html#ad97c28b35e6866da9c784e5bb051794e",
"class_f_e_a_t_1_1_geometry_1_1_mesh_extruder_3_01_conformal_mesh_3_01_shape_1_1_hypercube_3_012_d55947478d0079b597f439c47758a098.html#a582acf9c62e87287a56cd9bdd5b60242",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part.html#a2f5f806b61407754ca2bd01f9e635aa7",
"class_f_e_a_t_1_1_geometry_1_1_parsed_hit_test_factory.html#ac0f7547977bfb2e9e2691aa62686c12f",
"class_f_e_a_t_1_1_geometry_1_1_patch_mesh_part_factory.html#aed038bcf7b351b7b11f72359c6dab62d",
"class_f_e_a_t_1_1_geometry_1_1_shape_convert_factory_3_01_conformal_mesh_3_01_shape___00_01num__ddbfcb9c037d180b3bd81aab02e7610b.html#ad060758281c2ee20b0d08b91eff01c85",
"class_f_e_a_t_1_1_geometry_1_1_structured_mesh.html#af0e0db5896fd215feb52ae36fe0f2a28",
"class_f_e_a_t_1_1_global_1_1_gate.html#a0e168a2db82e0b71013b90322fad6c2c",
"class_f_e_a_t_1_1_global_1_1_nonlinear_functional.html#a14423566363f4eec07179c7ea151b111",
"class_f_e_a_t_1_1_global_1_1_symmetric_lumped_schur_matrix.html#a12c3313a7174c71985ab979732a925ab",
"class_f_e_a_t_1_1_global_1_1_vector.html#a74598b8ff7154381a22a369997f6aa1f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#abe3a1a4ba76a8f6a1870e895149e96d2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a3ad1c85eff5e55cfbc1bfe529e235d97",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror.html#a6cc66bf79640f5c7e77629068adcf4dd",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_none_filter.html#a2364d67ba141994d53e13a2970bc2326",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_e.html#a72aa8a2b990f5f44f38085c8364c1448",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#aab073914b60167bb30556976c78f479b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#a22588a7124bcf45c9cff7d6bb9daed1f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_slip_filter.html#a8874e4e8640ed973ea70dec7dbcf5210",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#aed9515b87e88602977928308733a52db",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a009d9bd9b95fea502d0fbe10358f73ab",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#a7915be74ff0aba8d97b530cf87c663b6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a9dae85fd4d7b9a2768e6cedcdb6859da",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#af5dc00e2137f8c37f5462b3423f7bd0d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_diag_matrix.html#a76ae4e72b459f51e77198a2899de4c76",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_mirror.html#a0686e72be212194cd0dbb57b04d8520f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter_blocked.html#abb59f993183f0c13910dc80bea1eb96b",
"class_f_e_a_t_1_1_meshopt_1_1_chart_distance_function.html#a3bab4dd99f0ed337269283295273c651",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a1b345cec59df9ccfd4702a7efb8eaefc",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function.html#af4613cd50db05475fb923fb474b86b66",
"class_f_e_a_t_1_1_property_map.html#a8a3bb18a30ad3e366912c14882d9f0e3",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a3b1a9cd30df7bf5a6d7c899d171c954f",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a40db3688d9794b0334a0226994b17c05",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#ab2e6ef76d11d3a708a85cf56d85d9caf",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#ae314aed268f854678d35272cd5e41340",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#ad5a06be09b1b40ab2c5b7cceb8274ab8",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#aabc5092b932e3669792d5e0867a2b76e",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a3a29e25850729825faf35116ecd66b9f",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a6d1131adeb401ad3964fcfd278e5b2c4",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#a6db06343bf3a13f014148a1e831c87cd",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#a44a696c771d3d261646f8de9a353cf38",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a5e276e3ab4c84a834377997024bf746e",
"class_f_e_a_t_1_1_solver_1_1_multi_grid.html#a8d7c5b8b824c6185867a4c090b8c3fb4",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a2f56865d47780501d77a27ed07daf11c",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a26458def0ef8654e3fff8f444df8108f",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a288adb7c820ebbb82a0a8b07f178db13",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a3a018a9fb70c9a5cb45047206c5f4d06",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a259e4ecfc3e51e9603b08e2c56ee3b0f",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#a731e6e902833e07bd086615675b51c31",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#ae0d3746b41b92d62d18e4e289d9ebd18",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#a50df7223074bef35a5a7cbd0c8ec9bec",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#ad62132dc8322e57ea884251664f8b822",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#ad586a1b68477198d76f4e794eaa83b5d",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a26c2873f682cb9f1722ab827d4f04a92",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#a82061dfdc0313a9da082736c7c48559e",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#ae314aed268f854678d35272cd5e41340",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_solver_base.html#a6f6fea1f8a375f49339ccb3e36e83bd5",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#a408d5c6093083db120a7ce801a371bb5",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_element.html#aaeb5cdc5568b733d0d86d55a2ea244af",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_element.html#a3e76ce09ab818c97e0be2b18357ff1a3",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___001b011599bed78529504db50cb49456d.html#a37e882b4a9a1a613ebad1b3244c43c06",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_node_functional.html#ac36f33ce44f2589fa5e008245c42dcb4",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_node_functional.html#ac36f33ce44f2589fa5e008245c42dcb4",
"class_f_e_a_t_1_1_space_1_1_element_base.html#a60113177e57c05c8dc561c93d1862f25",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#a79202b1f14cad175fb7aaeb5eecdaaf8",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#a1f5b8556e9fa1ee98993eeebeb5337c3",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_element.html#abdb421687dbb2dad0f79b33ef88c6c65",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s9a6fe2c5a6ce71c9b4fa308d8ff6eafb.html#afe9a67fa6c5b0994ca917225b7400c3c",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html#acac0922cdbfababeb05ff3a07e27b2da",
"class_f_e_a_t_1_1_space_1_1_node_functional_null.html#ad091f3ff7d0c9c412833b3a8446a43ff",
"class_f_e_a_t_1_1_statistics.html#a3bf61b6af4d00bfcba04af169a4e378e",
"class_f_e_a_t_1_1_test_system_1_1_test_list.html#aaa91bd3a9927b50298dc4777c416c74e",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#ace86d63e53e051cc11fc85c8e43b3f9e",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping_error.html#a04373352091dcf206f09e86fcc9a910f",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_012_01_4_01_4.html#a5be4bbf43479d8f523138fd6394fa8d2",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_6f4e5b7f5f9d6e2d04743ee3976d1657.html#a221c90b05d665aed17157cefab78da35",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_mapping.html#a8cb36075bac1a408d8409603920e670f",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_013_01_4_01_4.html#aaaf5ed4fd3dcc71e4f052170bd011fa7",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a0ef01e7a5c11e94a9fe2ce58d3b8f7d2",
"class_f_e_a_t_1_1_xml_1_1_scanner.html#a63d2f8a859df4d44a2d95757407d7931",
"feat_for_win_vs17.html#win_vs17_prereqs",
"lagrange1_2element_8hpp_source.html",
"namespace_f_e_a_t_1_1_assembly.html#a4d6f3cbc2e767ac015f0e462d6767814",
"namespace_f_e_a_t_1_1_solver.html#a56a1737a447bf0816faf68de552c2064ac14d2f1f2526ae3542b6c4f0f83ca91d",
"nonlinear__functional_8hpp_source.html",
"standard__tria_8hpp_source.html",
"struct_f_e_a_t_1_1_dist_1_1_operation.html",
"struct_f_e_a_t_1_1_shape_1_1_hypercube.html",
"struct_f_e_a_t_1_1_space_1_1_lagrange3_1_1_dof_traits_3_01_shape_1_1_simplex_3_01shape__dim___01_4_00_012_01_4.html",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01char_01_4.html#a93cc263a4de5df4e6d4aaadf0bc61e61"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';