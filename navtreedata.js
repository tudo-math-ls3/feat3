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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html#a3e2f70d9769a220323eacfc07a5b219c",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html#a667be493fe7fba4d29e225e06cbf5ed4",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol_pressure.html",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_rigid_body_vortex_pres2_d_1_1_evaluator.html#a1b624e7d52058a5599660145723d4a45",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_tensor_static.html#a37ccf0f0efa1ef23f56052fcf7431b13",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_x_y_plane_rotation_1_1_evaluator.html#ae4df19352cefc54ac4335428ac187eda",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function1_d.html#abd6fe4a42592fe279a13a11ac4699008",
"class_f_e_a_t_1_1_analytic_1_1_lambda_vector_function3_d.html#adfd28f5670e7d62c4fae9e94ef3b8fed",
"class_f_e_a_t_1_1_analytic_1_1_static_wrapper_function.html",
"class_f_e_a_t_1_1_assembly_1_1_basic_matrix_assembly_task_c_r_t_p1.html#ab5d4a86307142bee347d2c918e0fc735",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_matrix_assembly_job2_1_1_task.html#a28f6f83391139aa0be0fdb4f4bbbdba3",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_assembly_task_base.html#a442dcb2721d2dd367d5ceaeabe122409",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_vector_assembly_job.html#adae2e75ea23d8fce5804558721dd7465",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_matrix_assembly_job.html#afcecbec9a3d1cefb70e672ae895198a2",
"class_f_e_a_t_1_1_assembly_1_1_cell_error_function_integral_job.html#a0b2c06631d2292c1eb2a7fbdc95a990e",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_gradient_trial_operator_blocked_1_1_evaluator.html#a3cb162847584f978a1dd65e702089018",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_stress_divergence_operator_1_1_evaluator.html#ad4e31ab56875fc586a66d64769492f99",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler_1_1_worker.html#a1d06c26574b4d4ce9f9f4e4a7e43b70f",
"class_f_e_a_t_1_1_assembly_1_1_grad_operator_assembler.html",
"class_f_e_a_t_1_1_assembly_1_1_slip_filter_assembler.html#a797efe36cc84563bbd3405b89d187ccc",
"class_f_e_a_t_1_1_backend.html#a76990bdb5dbc08669f823b3c215f93e7",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control_base.html#ad1303539a3cd1760c77f851b0532fa0d",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_control_base.html#a7b191dd7a837ba8e0e86f4074d5c41f7",
"class_f_e_a_t_1_1_control_1_1_scalar_unit_filter_system_level.html#a4083f7e1800aa849ac36a06d051c5d6a",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_driver_base.html#ac423f8206f38140a4f4a664ae870e219",
"class_f_e_a_t_1_1_dist_1_1_comm.html#a065144405c5833ce42addeef33a4e9ab",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#add2c790bdae1965106af1c95d3551d72",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#af5c927a35758584141189fa69075d072",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#aa7a35726ccee0f1824fcfd0bb50afc26",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_surface_mesh.html#a2a9366f9b9d54027c9c4ea28d5b0b2ee",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#aadd250ae63ad8e1855dd486437988738",
"class_f_e_a_t_1_1_geometry_1_1_global_masked_boundary_factory.html#af7adef3d7a9518a7f9c7d169465674bf",
"class_f_e_a_t_1_1_geometry_1_1_mesh_distortion.html#adc33d6694ec6457399f1e184900fad4f",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node_linker_error.html#a76bcfc0b1c2ed8cc4d368d67c2c5c193",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#ade7bd2978b0005b7042040623f9cd399",
"class_f_e_a_t_1_1_geometry_1_1_patch_mesh_factory.html",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#af2367136c5b680432ce6646d88b0a002",
"class_f_e_a_t_1_1_geometry_1_1_structured_mesh.html#a5a5fcf762b9aa0b3417d97380228453b",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#ae517d835f8d148c9cc16f19f8cb61fa7",
"class_f_e_a_t_1_1_global_1_1_filter.html",
"class_f_e_a_t_1_1_global_1_1_muxer.html#aa5baeabafc5e0f2dc2111a9d5aee5364",
"class_f_e_a_t_1_1_global_1_1_splitter.html#a686df90e7739303e563bcf2e647dea77",
"class_f_e_a_t_1_1_global_1_1_vector.html#a16d2f034202beeb857c61ce448ab98b1",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a9060d0e4e26b6ca72132ec9504a8ec90",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_sequence.html#ac0fd1c8692fcb8f31065a494c8fb5a45",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_meta_matrix_test_base.html#aec7e586e14fd149269cc4f2340747ff7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_d.html#aaf6dce4065990914040da0542d1f857a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#a0d16a14bda7fc7b740b87c6711866b7c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#accb79a551092f52fbc5613bd72a44744",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_serial_config.html#ac4c50812d967aafc8366d190634705cd",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#ac096a8e63c8eac202322503a7c0ce868",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#adfded131472fdef6233c2957e724cb9e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#a253bca838ebbda28b193164988dae415",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a74c9279f1115605cd6ac24be25d4bdef",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#aa0f8a8c78b7e7b226dd476bd9fe924d3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_diag_matrix.html#a074b2fbf8efbfced7d43bb84e93713e5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix_row.html#a859adf41078d9c2da8a3368970b40950",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter_blocked.html#a157e29726c65993670f105d62550f46b",
"class_f_e_a_t_1_1_memory_usage.html#a5af240e4280dadc64deda23392939022",
"class_f_e_a_t_1_1_meshopt_1_1_du_dv_functional.html#ab3e4f01c20046d7195c5f96b1c29d612",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function.html#a41c4e59100a9e4cd10908e632a7aed3b",
"class_f_e_a_t_1_1_property_map.html",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a0810de389dc450bea5cb6917b66ca051",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka.html#a0c1b70ddda92eaf4c44cdf34e45f7602",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#a73654a64fe9f26443c5c12b460960983",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#aa7cf18c4f6e62acfd0c3cf3285c00b96",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#a9920812266b8f49921c9b133dd1451b4",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#a51a5b19738ed7d75ce3765d191232d54",
"class_f_e_a_t_1_1_solver_1_1_expression_call_uzawa_s.html",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a2daf610e7977a9e96ff84a6e027293d0",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#a1f17cb192fab20b3591eafdafd7c5f6e",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_invalid_matrix_structure_exception.html#ac922497e9e4a4a2e340e7173bb4d4780",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a1feb792955098cf117423f18dfd0db8f",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a288adb7c820ebbb82a0a8b07f178db13",
"class_f_e_a_t_1_1_solver_1_1_multi_grid.html#a269c4b8a4fb22d37a3c8593007515ce5",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_level_std.html#af7906a2bcdd3bbc2b2919be584c44a3b",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#aeeb10e81625c18084050c7cf00d25330",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_precond.html#a466d07fad209727943213dcc5e41f2bb",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#af8e655f0cdc31c3b245587b58a14d5b3",
"class_f_e_a_t_1_1_solver_1_1_nonlinear_operator_precond_wrapper.html#a1725f2b0500be33f4ccf18d499184a4d",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#a27c5e117b43720cbd96717472af0e517",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#a9920812266b8f49921c9b133dd1451b4",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a7647edeb0b1945e4a74e8d46e95a9c4c",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#a973d3942df838d68f63666265e73b025",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#a96fbaf2b10f359e84f9edc8240199445",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#a43bbf72e681f56ef9a6a96af37265a76",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#aa2b2f7492bf2c6b7d9f740b95e99969f",
"class_f_e_a_t_1_1_solver_1_1_schwarz_precond_3_01_global_1_1_vector_3_01_local_vector___00_01_mic9e5dece21cb7735a0bc194fb3569eb4.html",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#ad8716d5245471de11330a2f693d6f4f2",
"class_f_e_a_t_1_1_solver_1_1_uzawa_precond.html#ac55724acc579a6368efe8efeea2f7f6e",
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_element.html#a8d82bce6274b32899746656ddb66160c",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_64897ed7485eba824a67fb42fd5ab9c0.html#a2644f5cb97bfdc7ce8f437a8c9176b4b",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_evaluator_3_01_space___00_01_trafo_evaluator__802e55c907078015a390ad2e449f9e5b.html#adc9c506c493349492a7576a9c19fc304",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_element.html#adc5863102651ecfd7291ac3bf23e0310",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___0d5fd0fd65b9232ba00b88a0e0458fe94.html#a3ac2e2b361e2462561b68de850437eac",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_d256e15bbc7b3e4885146b0f6ddde539.html#a8e93574d8e24905aa2cad25e55bef35b",
"class_f_e_a_t_1_1_space_1_1_dof_mapping_single_entity_3_01_space___00_010_00_01dofs__per__cell___01_4.html#a593d66d5fd29a2a8eb8a800c7a45011d",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sp5bb7cbd2260293853cf4c05462b7a1c2.html#abeda000eff0270b2a43789f0148fb536",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s189823cc2b39e0d1a920e7698210faed.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sfefdf9344ae32fc63142f60cb7a695a8.html#ad5886f4911cac95096215def29085034",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s17de294e509ec5be8be361ba5ec7318e.html#a1882ec6cd8f12d25d29b3368ab837cb8",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sf226d38d10549efcd3a6fff34f617cde.html#ad9046d38bb6ce0666a023ca375e1f04d",
"class_f_e_a_t_1_1_space_1_1_parametric_evaluator.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_string.html#a7b50b063518d287f2fc6c795f6f57530",
"class_f_e_a_t_1_1_tiny_1_1_matrix.html#a90c43d2c2af4e1b5c20dd405c0f79341",
"class_f_e_a_t_1_1_trafo_1_1_evaluator_base.html#a997be4a421fa226b409f5746b9bf018f",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_012_01_4_01_4.html#add9092926384fb5314de0140b987b0e9",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_mapping_1_1_evaluator.html",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#a1f6c101a92f1382b729c84b77b12a661",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a5e160dfb0c7a1fa8bd9cb5f31206cc88",
"class_f_e_a_t_1_1_xml_1_1_grammar_error.html#a4663f51d4270a50368f9bd0681400447",
"dof__mapping__base_8hpp_source.html",
"ini_format.html#sec_example",
"namespace_f_e_a_t.html#abc175499504e9be922533b3a3627ccac",
"namespace_f_e_a_t_1_1_pack.html#a45291a05866291324a62a68598399a1a",
"namespace_f_e_a_t_1_1_tiny.html#a6ad53ed9caa8145a5086f4d1f1aba469",
"preproc_macros.html#preproc_sec_user_def",
"struct_f_e_a_t_1_1_assembly_1_1_scalar_error_info.html#ae89dab4c35aa6806d800712c5dcf8a04",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html#a7bc9870851acb52fd0879bcc64b7a5cd",
"struct_f_e_a_t_1_1_space_1_1_bernstein2_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercud0df8e3d6307173e0e27fffdaebdea90.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___000b42a86895f896d40a64f335ad2bfb25.html#a454c814e39bb6aee89a3d04258830407",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01double_01_4.html#a2979c4e6c0138ce52f183d70de7c17ef",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_space_helper_3_01_q2_standard_f_e_3_01_shape___01_4_00_01_d_t___00_01_i_t___01_4.html#a78ff1eca0d0f3eebfdc46bd58aa6a112"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';