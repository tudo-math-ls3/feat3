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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html#a7c1c46728c31d39ef595ee3af6580d1f",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html#abff6abc7e9840424ee6f105fa015dee2",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol_pressure.html#ae8342b7f2e59ce20083d974b06c5084f",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_rigid_body_vortex_velo2_d.html#a2757b0ecfda171d8c3a4b06246398213",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_tensor_static.html#a83d47b894377f151209eb71657cfa175",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_y_z_plane_parabolic.html#a0e8cdb958b827b25e3625ae6a6ebd3aa",
"class_f_e_a_t_1_1_analytic_1_1_gradient.html#af2e048546e4229ad7b23203dedb5763d",
"class_f_e_a_t_1_1_analytic_1_1_lambda_vector_function2_d_1_1_lambda_set.html#ad3ba6d4da7635a9ac5d77b5d18dfadbf",
"class_f_e_a_t_1_1_analytic_1_1_static_function.html#a696f559552227ede6d9e8356b8b3c0a5",
"class_f_e_a_t_1_1_assembly_1_1_basic_matrix_assembly_task_c_r_t_p1.html#a4045342172197b0bb0bfd5fef5738d8b",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_matrix_assembly_job2.html#a1198fbb05df28db7572bfb3ab9e3bb5f",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_task_base.html#af80337cfa1a32793a469d2656f292426",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_vector_assembly_job.html#a34f836bb6de7a3f2304de70d46516716",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_matrix_assembly_job.html#a3b10e6340b78da910ce99c28d24ebcd9",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#ace416de1c80d272e9bfc9c4d5ca0aefa",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_gradient_test_operator_blocked_1_1_evaluator.html#a1d37a2ea952006926e0c71ffd0c6dd83",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_strain_rate_tensor_operator_1_1_evaluator.html#adac0216b546ff948b489bcf43a6735f7",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#acbb6f1a306815a85e5366bda3e1509cf",
"class_f_e_a_t_1_1_assembly_1_1_function_integral_info.html#a4259e2a593b42cb8c390347bd1d2d6d0",
"class_f_e_a_t_1_1_assembly_1_1_rew_projector.html#aa92be8bc60f19e289729348d491954c5",
"class_f_e_a_t_1_1_assembly_1_1_unit_filter_assembler.html#a3277a73270d3b240ec579fec72691ec1",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#afddf5df0ac7908f4b270aa208ed76012",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_hyperelasticity_functional_control.html#a40543baedc948e935edf2779e3422b13",
"class_f_e_a_t_1_1_control_1_1_scalar_combined_system_level.html#ad81bdc8156d4b005ce5e4767598d6939",
"class_f_e_a_t_1_1_cubature_1_1_hammer_stroud_d5_driver_3_01_shape_1_1_simplex_3_013_01_4_01_4.html",
"class_f_e_a_t_1_1_cubature_1_1_symmetric_simplex_driver_3_01_shape_1_1_simplex_3_012_01_4_01_4.html#a8fe0b70f58b273631c40d512e4cd3235",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#a4508c9c82d78f1ce89f74b620a231887",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a88d5e49bd7b8e544c8b8821f837f8a6e",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#a2724665b990d064dc48f63e61f95cd03",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#abf159301aa38e5bf4befc1108a7bae1b",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#a2eaaa4759effda6daa3b44a9163a0415",
"class_f_e_a_t_1_1_geometry_1_1_global_masked_boundary_factory.html",
"class_f_e_a_t_1_1_geometry_1_1_mesh_atlas.html#a9d8799d745fe25f62c0395272740121c",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node.html#ac3c13c8f38d7fcabad8f6e555706548f",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#a3f1cf371ee897a690cc8d61db07225c4",
"class_f_e_a_t_1_1_geometry_1_1_partition_set.html#abbc31164285f934201f591fb9135a509",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#aaaaac3c1c4cb929cad8d7f4ad6ee44c2",
"class_f_e_a_t_1_1_geometry_1_1_struct_index_set_1_1_image_iterator.html#a53eec4cbef7fadef0a5e422adc7b7af2",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#aa2e1037003b2d6a9c14b100318b7c6ca",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_vector.html#a0462b5109a66475a2f600ece2f05cecf",
"class_f_e_a_t_1_1_global_1_1_mean_filter.html#adf5f5ba82048f0efbd8acec8a70c8f61",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar5e6f49c15ccb43840cb894a3375f33c3.html#ab923dbbd9e72776d4073dfa7b8f7e128",
"class_f_e_a_t_1_1_global_1_1_transfer.html#a608c2083b393fc35cc3ab11618696bbc",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a45d3730c71ce42bb9586567acb41638c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#aa960f3412469e8fad2dc681700bce955",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_chain.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter.html#a1c235f628487bac78619cba80c59a75f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#a58f328a9bd265fe665059fc2d71b07a7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html#a8994516a57c288a94e88ce6a38dfd32f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_full_matrix.html#a9e4e6c522f78b17d1890ff769558f85d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#a84d9aa21f0225e4fe352a854f621523b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a011b5e99833fcf98ee2a898fe724aed0",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a22181615f6beb1d9b1844165d269c419",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#a7d0d7ac83684b2e3d144d9cd8d29f166",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a9ebf80cb3b4cff2eb77876df0065d44a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#af5dc00e2137f8c37f5462b3423f7bd0d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_diag_matrix.html#a72c67d3f8d60e6e142e053ac6cd18eab",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix_row.html#afecdc7b41c31bcba1d7ba9babf38f366",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter_blocked.html#aafe11b47befb7332f523d1990b79bd4b",
"class_f_e_a_t_1_1_meshopt_1_1_chart_distance_function.html#a195e4725ec42ea2adc9ec70b6d8ec4e7",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function.html#ac126cb4f277172bd8930cf338015cb5c",
"class_f_e_a_t_1_1_property_map.html#a5c6c5247a02b665f0b0ffc23b1a79eb2",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a288adb7c820ebbb82a0a8b07f178db13",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a3430f79a13866b83205fd94d1e3ea454",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka.html#a9b9d1846417be83357859c2763c89a2f",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#a9ea9dc273f42b64becae5cc2071db416",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#ad355136ba5cd50b1f2ece6658ab1c0c8",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#aa2b2f7492bf2c6b7d9f740b95e99969f",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#a5962d67caa85a29115abe6f04c24adfa",
"class_f_e_a_t_1_1_solver_1_1_expression_call_precond_l.html",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a27ba39cc738b9a4fa6208feb779ca42d",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_hypre_solver_base.html#ae6bedcbf11e2f8710f3c0a9b577a4e35",
"class_f_e_a_t_1_1_solver_1_1_i_l_u_precond_with_backend_3_01_preferred_backend_1_1generic_00_01_aed6c8770199c660f2b23c3052158c26.html#abb820161a5311c363eb00e458f8be9a7",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a1552e7e7182ac79d0c55149fe5347358",
"class_f_e_a_t_1_1_solver_1_1_m_k_l_d_s_s.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#aeeb10e81625c18084050c7cf00d25330",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_level_std.html",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#ac93194a978f521b1a3e4ce546a3cd757",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#ad586a1b68477198d76f4e794eaa83b5d",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#acf7a6a13468eae67e3446f004032bf21",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#ad586a1b68477198d76f4e794eaa83b5d",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#aec76c75b2d2582ef2cf0ed928b583314",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#a5e276e3ab4c84a834377997024bf746e",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#ac59b69cc5fe017dd82287034e6a25a5a",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a2d49a7a1f657c2072f629e67d2820de9",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#a4fdec11243202ce55a49e67d259cab2a",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#a4d6b6c7809c161f59cbaf7cb01149178",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#aa5524e82d3d153e57af1787421c0ec90",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#aec76c75b2d2582ef2cf0ed928b583314",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#a5bfdd86c9fa12e542b4262ea6a219090",
"class_f_e_a_t_1_1_solver_1_1_saddle_umfpack_mean.html#ac1ad9f737ddb43594e23c50dc42bb8c9",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#aa7cf18c4f6e62acfd0c3cf3285c00b96",
"class_f_e_a_t_1_1_solver_1_1_umfpack_mean.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_voxel_ama_vanka.html#aab268be7c22a62d9dc250d42e242a92b",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_045e90cd69ef865eaf5b37b181a99332.html",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_element_1_1_dof_assignment.html#ac3756f73cee50e45d26bb29c71c19f8b",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_node_functional.html#a12dd52ec74b75aff5d9471e1c3572485",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___0b6ac4939673a6cd90495f3edd5624cc3.html#a4c6b7814611a0f215f825559c4872e45",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_b2c6b42d6f9b9e25d469d3575a373a61.html#a29bafa5636409dda3697334dc98bdf3c",
"class_f_e_a_t_1_1_space_1_1_dof_mapping_identity.html#a593d66d5fd29a2a8eb8a800c7a45011d",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_element.html#ac38a8b8a3720601e25e0973b4d8bdff3",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_node_functional_3_01_space___00_01_shape_1_1_simplex_3_b3c8e2911d072d6c232081658e8d85f6.html#a12dd52ec74b75aff5d9471e1c3572485",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_se3129d8853eb9690234f254237b7f1ec.html#a6c3d4541971ee8acd4547988657d94b1",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_node_functional_3_01_space___00_01_shape_1_1_vertex_00_01_data_type___01_4.html",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sb3aa60112044c986d3e79f1b1403fdf7.html#a419e835cb193bf3e549d12b7e5132aed",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s88a2ce84347f190b629ca7916f4787f8.html#a4b1ac4e614fad291b9dd1b93491c8342",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_3da8886277cafb8343dcd8258f6373f1.html#a039d795545bbd7867c97e252a4ed7113",
"class_f_e_a_t_1_1_statistics.html#adb408a83d079531eae9b5ef89c026f98",
"class_f_e_a_t_1_1_time_stamp.html#aa13f6c5955f4b04dc00ff264b9a0f842",
"class_f_e_a_t_1_1_tiny_1_1_vector.html#ae0e56c9cac1ec4b3ded55c256e13a840",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_011_01_4_01_4.html#af0b5c9b9c23d707f271b3b3d7c1e0de9",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a577665c3357b7ebf50ded1a08b0ceaab",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree___00_01_shape_1_1_vertex_01_4.html#a2f25ae5e56324c6744203fb710e0b1ee",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator.html#a63bb4db2d02cb15c4a26d76690d7ad5f",
"class_f_e_a_t_1_1_voxel_assembly_1_1_voxel_poisson_assembler_3_01_q2_standard_hyper_cube_3_01dimb8f4748e519b8a0db5e2dfec2b4d07be.html#a1ddf3b807a7783d8e86ef656bc27c50d",
"component__product__generic_8hpp_source.html",
"functions_type_d.html",
"mesh_file_format.html#meshfile_chart_circle_examples",
"namespace_f_e_a_t_1_1_l_a_f_e_m.html#a82e983d75b2911bc2cd014a541d50711afc9847d479fb14be64841af20cde25c1",
"namespace_f_e_a_t_1_1_solver.html#aa2e0c5a93d63bfa65939ba63207d6f96",
"parti__zoltan_8hpp_source.html",
"struct_f_e_a_t_1_1_analytic_1_1_eval_traits_base.html",
"struct_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier_traits.html#adeec0b78fb0da098603570203269150f",
"struct_f_e_a_t_1_1_meshopt_1_1_rumpf_trafo.html#a1e369be82aadb85474e9481c70f66d3c",
"struct_f_e_a_t_1_1_space_1_1_lagrange1_1_1_dof_traits.html#a723dec5cb52f91dcbc4963abcd284758",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___002b260e771ecb4737ffdd677d6d46ebb5.html#ade71c5704f62d028ffaaca9654e78e2a",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01signed_01long_01_4.html#ad8d1f9f2d042ba7b327e789802bb64c3",
"tools_meshtools.html#meshtools_mesh_extruder"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';