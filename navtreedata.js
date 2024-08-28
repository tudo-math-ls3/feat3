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
"class_f_e_a_t_1_1_assembly_1_1_function_integral_info.html#a7265c06cfd229ca95a81f0796959348c",
"class_f_e_a_t_1_1_assembly_1_1_scalar_discrete_eval_data.html#a6ce634d2e25ff7890a048fba23940922",
"class_f_e_a_t_1_1_assembly_1_1_unit_filter_assembler.html#a805f3495f573fa5bfaddee50a717b144",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control_base.html#a1509ad37235a2b298961a25caad3395e",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_hyperelasticity_functional_control.html#a769267847e08cecc56cc53a7c6039f40",
"class_f_e_a_t_1_1_control_1_1_scalar_mean_filter_system_level.html#a0c0a1e555bf2a343005289eec432cae9",
"class_f_e_a_t_1_1_cubature_1_1_hammer_stroud_d5_driver_3_01_shape_1_1_simplex_3_013_01_4_01_4.html#ada60bd9977396850624e62d3f4c0e67d",
"class_f_e_a_t_1_1_cubature_1_1_tensor_product_driver_3_01_shape_1_1_hypercube_3_011_01_4_01_4.html",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#a4cc05eb2bc19939aa879db47a2c0ff27",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a93ba9bfe598632e8c3d671140a882b59",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#a2ea8f8c0aa6074b687abc5a2bf390e8c",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#adc63154b796271f060d449593562a1e2",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#a30ec1b3c48e3752e574cef3cbf66a41e",
"class_f_e_a_t_1_1_geometry_1_1_global_masked_boundary_factory.html#a16e12db98573f4717cc9ac329a71aff2",
"class_f_e_a_t_1_1_geometry_1_1_mesh_atlas.html#ab15ef445a62f53e51726601d978dd8b9",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node.html#ace55e8912f1468f724045671d1a3dfc5",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#a48bf531fd2b533d0499f3875d884ffe9",
"class_f_e_a_t_1_1_geometry_1_1_patch_halo_factory_3_01_geometry_1_1_conformal_mesh_3_01_shape___e75c5f271492b3ed0db46ab700c23399.html",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#aae20ce09f918ecb9fd9eb22c6377a64d",
"class_f_e_a_t_1_1_geometry_1_1_struct_index_set_1_1_image_iterator.html#aaf4b16526a160fae7450ba1ec905aaa1",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#aaaedb5285c65492271324cc3cf56dd7c",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_vector.html#a18650e971d3d31726062b5387f242035",
"class_f_e_a_t_1_1_global_1_1_muxer.html#a07305c8dead88e0023cfbc720ee33f7b",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar4d3a70ff49e16f0d9ba9d236f0dafdc0.html#af3cdfc3054594504d9df7909a2af50dd",
"class_f_e_a_t_1_1_global_1_1_transfer.html#a7038d633d8d07eb9712806c6c030d09a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a62d1fc8caae749e4d4fbdde9fd2d66af",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#ada2b8bd276f75ab9df2a22cb52f48e7c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_chain.html#a5e78e76f6644369d5c5c909f6c4daf94",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter.html#a70e4d91539875674f5f6786bd77e9470",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#a937950632dddf68407c4bca5678423c7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html#ae8cb31ae7cb4019e937648ec01c00a2e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_mirror.html#a50c172559ba491b4f64b1539aa4cdb74",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#affb0e3e9e6c38fa0affa4be400c266f5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a3cd9033a51eb7f18d03f6bdd6a8fb144",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a47b718987b62807efb0461015eb18b56",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a958343c0926a9f67dffc5997c3c05b31",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#ae79de30b391f87aa8d10abc4189fdba7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#aff0d5603f6fdb56d4634a28b1ebc339f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#a958343c0926a9f67dffc5997c3c05b31",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_filter.html#aea6ac3243d61d2f8f03792bed9f03581",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a57fdfb40250022d5b847bca837a6db5a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a77af3799df4759e98b5e5a4303969788",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_default.html#ab09cdac7b8b9aa9f2cdfca5123a33918",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a8bb2b72b93e8a5c94d214e9193a12f07",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#a2ea3b5694b5064d3e0cc5a453d04a890",
"class_f_e_a_t_1_1_runtime_1_1_scope_guard.html#a24fbdb75acf6059d7254290936840c42",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a88f7e43b3f7a18196862c64fc3c4d034",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a97f6e79dbe58d08330cec6ba5c0195c3",
"class_f_e_a_t_1_1_solver_1_1_analytic_function_operator.html#aced31cd180995ec979f55c6c2eea407e",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a0628f8a1a0a10403066c86e27ab103ea",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#ab934585a0eb28b50149e32ff156685bc",
"class_f_e_a_t_1_1_solver_1_1_convert_precond.html#a8f13dcde65e5d0d85f60bae8433c4cc5",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_global_1_1_matrix_3_01_local_matrix___00_12c7f8d5d38cd40f20a5510074419a06.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a6950d64c2b9db6d1b76eff702a2bffa1",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#aa7281f8e3498245d411c5829b9d0cec0",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#a73654a64fe9f26443c5c12b460960983",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a77f32fe5475099d5ae331eb7b18d9280",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_multi_grid.html#aa0f03f9448519ec79a105bbc7db3ce95",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a3b1a9cd30df7bf5a6d7c899d171c954f",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a27c5e117b43720cbd96717472af0e517",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a293d87d12e78eb731f41496ed2bf8dec",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a3b1a9cd30df7bf5a6d7c899d171c954f",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a26458def0ef8654e3fff8f444df8108f",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#a73654a64fe9f26443c5c12b460960983",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#ad0e0aa027e766929e7e741f6384f887a",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#ad8026d89fbcbd8d65ffdc207d15b06f6",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#ad5a06be09b1b40ab2c5b7cceb8274ab8",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a27c5e117b43720cbd96717472af0e517",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#a881142b80ffaca87c955550651e82b25",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#ae6fc9035510853c34e062754732c17ba",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a116b3f6b1736900d75d29aa0d7f644df",
"class_f_e_a_t_1_1_solver_1_1_solver_base.html#a70018a5e4ca2cd3dfc0cac9db3df051e",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#a445a63cba717d2c14abc4db6f2d19047",
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spa96d47ea9782e73b0166b29cca748b7d3.html#a5a22e90ba4d55c88908d6d01a31cea6c",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_82866b3cc738de4ae4935135fac5e732.html#a24e648bd29f80680fe2028b3fe286673",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_evaluator_3_01_space___00_01_trafo_evaluator__cefb60ca526d73ea9fbcb2c4bda04e64.html#aa0233c0c9c72793576b89f606d3553f9",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___001b011599bed78529504db50cb49456d.html#a499f5fa0db91640392c8586cee3917bf",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_node_functional_3_01_space___00_01_shape_1_1_hyp14d51a68531283eec3ff37b863183932.html",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_node_functional_3_01_space___00_010_00_01_variant_62f9079e75faf0a6f75d99b10b03a599.html",
"class_f_e_a_t_1_1_space_1_1_element_base.html#a442e1cf4b479140937bd7e44aacfe277",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#a6db7e5c3954c49646d2823108329675d",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_element.html#aa7119f0e4691ad7f99b53be2e0071ba7",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_node_functional_3_01_space___00_01_shape_1_1_vertex_00_01_data_type___01_4.html",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_element.html#a7fbbea2698e84601deaa163443f8f955",
"class_f_e_a_t_1_1_space_1_1_standard_scalar_eval_traits.html#a7a6c2291f10143353b685949e45216fe",
"class_f_e_a_t_1_1_string_mapped.html#a92ac894dccb24fcb0eb990adb7b012d8",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#a19050b771bd8b4a256c0246f80139bd5",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping.html#a7fa0ddc85db0fc846c93ef9755216ce4",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#a9eb11964c43671facf3b5a358bcfc767",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#a7653399d919df8bfcd035819d872e467",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#acdb7bfb197ee64e4d83a5901a0aaa456",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_mapping.html#a045e3535fe6e1f9de7c82d504c5e3cb8",
"class_f_e_a_t_1_1_xml_1_1_scanner.html#a40d656eb15388e7bd7277c36f2afb8c6",
"face__index__mapping_8hpp_source.html",
"kernel_2lafem_2base_8hpp.html#ad372dae2e320ef82768327bb8342786faebd9bec4d70abc789d439c1f136b0538",
"namespace_f_e_a_t_1_1_adjacency.html#a85fc833108471df98e28d6b06324f3e7ae76496363b1788de189d5a6f5d63f80c",
"namespace_f_e_a_t_1_1_pack.html#a7c69c0366e113fd42901ba9bcea55d4dae7e62f6928f76df671b5a0379793fab6",
"namespace_f_e_a_t_1_1_voxel_assembly.html#a2a3b61e1647254610f980a46530d2ea8",
"rumpf__functional-eickt_8cpp_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#a8e78f41a6f18578d8ec5b75f65d82e9a",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html#aff2e00b7b046783b4d60591e71418f1b",
"struct_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_dof_traits_3_01_shape_1_1_hypercube_3_01dim___091_01_4_00_01dim___01_4.html",
"struct_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercu46861c7ddfce2b4bef65304cdf16ecfa.html",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01bool_01_4.html#a616ccef26088f890f0040d4032008258",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_space_helper_3_01_q2_standard_f_e_3_01_shape___01_4_00_01_d_t___00_01_i_t___01_4.html#a103109b1b817235b0f98c90f34fbe67f"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';