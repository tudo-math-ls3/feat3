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
        [ "FEAT_COMPILER_MICROSOFT", "preproc_macros.html#ppm_feat_compiler_microsoft", null ]
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
      [ "Voxel Map File Format", "voxel_map_file_format.html#voxel_map_file_format", [
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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html#a8ced5640e35e5dd66b6402d016e2a85a",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html#aca18359aabfd7fd8a5e20bb1de3d5ba5",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol_pressure.html#afd9ed0af471602552ea19d3899ac105b",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_rigid_body_vortex_velo2_d.html#a611acaf4638f9cc4d259f5080d784b3d",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_tensor_static.html#a849a033f8c2eb2c6bd74605dc13f486e",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_y_z_plane_parabolic.html#a2ffa0e2790f316194a550da5081eed5c",
"class_f_e_a_t_1_1_analytic_1_1_gradient_1_1_evaluator.html",
"class_f_e_a_t_1_1_analytic_1_1_lambda_vector_function2_d_1_1_lambda_set.html#ae20fa78eaddd6a448d366a654db99ed4",
"class_f_e_a_t_1_1_analytic_1_1_static_function.html#a7228e78d25d9d683af4449549f310396",
"class_f_e_a_t_1_1_assembly_1_1_basic_matrix_assembly_task_c_r_t_p1.html#a4811516bf8862a5133339e7ee31cfe6c",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_matrix_assembly_job2.html#a1868b89af605e9383a76315629f5b923",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_assembly_task_base.html",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_vector_assembly_job.html#a3b10e6340b78da910ce99c28d24ebcd9",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_matrix_assembly_job.html#a42373c206757aea13cf306baf8c69cee",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#ad1ffa1d8dd3067c503578fdace46faed",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_gradient_test_operator_blocked_1_1_evaluator.html#a3445d9c7997262f2e2692c1c5aa769a8",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_strain_rate_tensor_operator_1_1_evaluator.html#ade1dfe8eeffd4a15751521335a2c24f1",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#ad0d72b39e71c14e4d3895dea1cf09080",
"class_f_e_a_t_1_1_assembly_1_1_function_integral_info.html#a497ffb9465ff124997d6ab5c33bf3fd0",
"class_f_e_a_t_1_1_assembly_1_1_scalar_discrete_eval_data.html",
"class_f_e_a_t_1_1_assembly_1_1_unit_filter_assembler.html#a39c56fb5d2145d0ea39f9eb4d2a46459",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control_base.html#a00f555bab5a53e5d0920843655f02541",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_hyperelasticity_functional_control.html#a4e9909e4afc21fbf98fce30912c0459d",
"class_f_e_a_t_1_1_control_1_1_scalar_combined_system_level.html#aed3c7548f8c24dc41c012b20e3755802",
"class_f_e_a_t_1_1_cubature_1_1_hammer_stroud_d5_driver_3_01_shape_1_1_simplex_3_013_01_4_01_4.html#a41db8460ddc3af8ad250dd09fe9c4d81",
"class_f_e_a_t_1_1_cubature_1_1_symmetric_simplex_driver_3_01_shape_1_1_simplex_3_012_01_4_01_4.html#ae5e3aef395954417b29593bcf95e80cd",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#a485087d415ee7c01b42bcd4f31e1588e",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a8b234feb91c3fc1e4b24b8fbc432ecc5",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#a2b037a7a14247a41a9b68ba4b2168951",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#ad72b359245ac85712bd6a6a411b888ba",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#a2f3158ab506ca75ae1bec39ede7e89a1",
"class_f_e_a_t_1_1_geometry_1_1_global_masked_boundary_factory.html#a12f4f81c63436b274c694881954674df",
"class_f_e_a_t_1_1_geometry_1_1_mesh_atlas.html#aaf6c22bf8d1fba810f169e8a564c2f0b",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node.html#ac57ab6ce12d574dc738386d9b852bb86",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#a413b2a2c6035ce5a1673142d76d4cbe9",
"class_f_e_a_t_1_1_geometry_1_1_patch_halo_factory.html",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#aada548a249edd9a417fed9c6daaf9d0a",
"class_f_e_a_t_1_1_geometry_1_1_struct_index_set_1_1_image_iterator.html#a9e46bc18a3c6b582fe37307c011c9502",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#aa3b130254f3bb341d627446138f914c6",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_vector.html#a165ba822f747e9acb09effdd7bef23e1",
"class_f_e_a_t_1_1_global_1_1_muxer.html#a0581479a30d0d0da99ec49e5395ed395",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar4d3a70ff49e16f0d9ba9d236f0dafdc0.html#af3a4e05bd34103541d9a4d98ff571a38",
"class_f_e_a_t_1_1_global_1_1_transfer.html#a6c7af5fe39e8bfeb19903e1a8c3c596f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a5de387b906c5e065448dbdb7703f8dee",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#ad486486693c22a3b456aad94630b2ada",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_chain.html#a591fac466b59c692fe1ed6b3d866d958",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter.html#a69d09187108f4a1a276af6dbbfbe8e2d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#a8d65695f73c9a296ee7ccc4755280307",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html#ae00dc7021f0d34405b7ff58ec6b5d6b3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_mirror.html#a3c8682f128efee6d01a7292bf38cc77c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#af8a5fc749952c24eb84604d48eb74eca",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a3c76c031a06c6405891174efccf5b876",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a4780b9550664ac6ae3f14c2ea7cad700",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a9060d0e4e26b6ca72132ec9504a8ec90",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#ae193114a4761cd2ad7440c2fdb9d700b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#af8904c03021f10e3c892b25bf1a9c51f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#a9060d0e4e26b6ca72132ec9504a8ec90",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_filter.html#adca7971a01f336b59f9a5fd0971d54e9",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a55e68321bded65673581ade0db918fc6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a77a6db9e56999e4cd7bd665806eb5e36",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_default.html#a82192885f0d9bb5d98ab04a47695b6d0",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a8b2690f0081d027191122f2f8a8fddfd",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#a2e1fec781341ce08ad76a469b87c5d6e",
"class_f_e_a_t_1_1_runtime_1_1_scope_guard.html",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a87caaf864220f9cb8add2b0eb90e2ad2",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a973d3942df838d68f63666265e73b025",
"class_f_e_a_t_1_1_solver_1_1_analytic_function_operator.html#aad23ad7ac46dead39909d6f8b067ccb9",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a02fd7ffda14e175b1abf044d706876c8",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#ab8b1f9ade67e6d44803a9b0e9128d94f",
"class_f_e_a_t_1_1_solver_1_1_convert_precond.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_global_1_1_matrix_3_01_local_matrix___00_12c7f8d5d38cd40f20a5510074419a06.html#a36f26920a94a7678c21690b7ec48f7b5",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a8ed25011c486ac360b5705a54aa004e7",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#aa5524e82d3d153e57af1787421c0ec90",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#a87ca75831356cc78d31b03b197fac8f4",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#a731e6e902833e07bd086615675b51c31",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a76a7e3315fcb44c2a3ac466ae1ee1337",
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
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spa96d47ea9782e73b0166b29cca748b7d3.html#a534b5f5f994a226794ca55a4f7f3bcef",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_82866b3cc738de4ae4935135fac5e732.html#a218c56b45a6fb44568f536f496afda81",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_evaluator_3_01_space___00_01_trafo_evaluator__cefb60ca526d73ea9fbcb2c4bda04e64.html#a989106fa0a904e55dc673397b97fb6ef",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___001b011599bed78529504db50cb49456d.html#a3b5673df67ee682095f0264ccb25904f",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_node_functional.html#af3b7b7f553156a6574c71c886f760a7c",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_node_functional.html#af3b7b7f553156a6574c71c886f760a7c",
"class_f_e_a_t_1_1_space_1_1_element_base.html#a437479f3de53315ef4f5349e588bcf15",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#a6bfdd412851516ebd9e42a508b51edb6",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#a0033f9bbad229aef82a64d05d1149d01",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_element.html#a8d82bce6274b32899746656ddb66160c",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_node_functional_3_01_space___00_01_shape_1_1_simplex_35518806867547cd79f3d31eba693e2b8.html#aa6cb95c632710f9e50e675a723bc5532",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_element.html#a79a0ac427b4d7a25aabbc912bc6ab8dc",
"class_f_e_a_t_1_1_space_1_1_standard_scalar_eval_traits.html#a772e25b153101dabfde18682173c2d3c",
"class_f_e_a_t_1_1_string_mapped.html#a69f8735082a1a7abf6c17f08e838d6ba",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#a141d74be9d7caa7f8963fc152de657eb",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping.html#a68e490a0dac77cb70cc9ad09f52fcc08",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#a99aacf4b77b94a83e0558da0601a233f",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#a74c0fe3b32dbb480c017f185d59b2779",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#ac80bdf67ee98f7557dec2d66c6be6346",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_mapping.html#a018fc3467811ca7f0fdbb0611fd6ca4b",
"class_f_e_a_t_1_1_xml_1_1_scanner.html#a3b489fcdbaf1dcc52769df52cab021c1",
"ext__vtk__writer_8hpp_source.html",
"kernel_2lafem_2base_8hpp.html#ad372dae2e320ef82768327bb8342786fa928d0f1285ee7d36c1c2fa1b1b7a164c",
"namespace_f_e_a_t_1_1_adjacency.html#a85fc833108471df98e28d6b06324f3e7aa64cbf11a4a0dd8ce3b57b86d42ee06b",
"namespace_f_e_a_t_1_1_pack.html#a7c69c0366e113fd42901ba9bcea55d4dae203531cddee75a95323e3d6977c4789",
"namespace_f_e_a_t_1_1_type.html",
"row__norm__generic_8hpp_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#a7258cd7a7e2e67eefd55115862b942fe",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html#af74641e4503a92a36b06462d38c8b7c6",
"struct_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_dof_traits_3_01_shape_1_1_hypercube_3_01dim___01_4_00_01dim___01_4.html#a27db0ef70b995859c9d868052bb1dc8a",
"struct_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercu20b7eb172803d372336a0173fae7ea6a.html",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01bool_01_4.html#a4c1d25020d1224116237e25274fc3409",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_space_helper_3_01_q2_standard_f_e_3_01_shape___01_4_00_01_d_t___00_01_i_t___01_4.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';