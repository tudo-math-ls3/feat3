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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html#a7c1c46728c31d39ef595ee3af6580d1f",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html#aa07d2b3527180315ee02981a202667e3",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol.html#acdec0103d79a95cb46b9bc7b1910d2a5",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_polynomial_function1_d.html#a39ad59ea3dd2e4975144af9a24a68700",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_r_h_s2_d_1_1_evaluator.html#a777899b492ca80bad0bf5676982585cb",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_tensor_static.html#a9ab44a0818a8aa7c977973700427637e",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_plane_distance_function_s_d_1_1_evaluator.html#a203af6f7321f123a46df73ff0c0c9e6e",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d_1_1_lambda_set.html#a7b2a2574ac58cd0c478bec8e3b35bab4",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function_1_1_evaluator.html#aded30ae1fc332099109b11884694bf68",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits2.html#aa258c13338ea07088ad4b66d3554e06d",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_assembler.html#a4589b56f2b5e35b9e45052a2f3ec9cbe",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_task_base.html#a305d72e8357bbf4698944acbe7847cc2",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a6d0f5495e4ce7f4577412e83e7f7ff88",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a725222cb930c886a1cc4d50b0d14ebca",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#a09c40fffde21adf120c338152bbf6827",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator_blocked_1_1_evaluator.html#a10a1460da5b6a521c6be97c399a82b8d",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_operator_1_1_evaluator.html#a348d6a5e0dc6e06fe71f59b92f1a51bc",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a047dfbc486082bbf11e5421b4c556a18",
"class_f_e_a_t_1_1_assembly_1_1_force_functional_assembly_job_1_1_task.html#a447bf1198b3bf57ad3a7eb1683b54bdf",
"class_f_e_a_t_1_1_assembly_1_1_matrix_discrete_eval_data.html",
"class_f_e_a_t_1_1_assembly_1_1_symbolic_assembler.html#a7877c84ca56a475a47996934263911dd",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_domain_layer.html",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#af220c1f03fdfe2c396185c576cdf4575",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_quadratic_system_level.html#a10a0672ceb90c26b00b3f55a760c0de9",
"class_f_e_a_t_1_1_cubature_1_1_driver_factory.html",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_newton_cotes_open_driver.html#a9dd57fbb1f2cce703c67ed1a4d97b711",
"class_f_e_a_t_1_1_dist_1_1_comm.html#acd1de9d2509adb85bfcf71b7fd9c6e5d",
"class_f_e_a_t_1_1_file_error.html#a9b1719a6dd6d6b64af5a1b065b0cb950",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_chart_base.html#a7619901767a0842436e3b8c616b6df81",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_extrude.html#af084dcbf0c45a239d8b707454bdfb413",
"class_f_e_a_t_1_1_geometry_1_1_attribute_set.html#af4186f55b3934faa50254bcea9b685b2",
"class_f_e_a_t_1_1_geometry_1_1_export_v_t_k.html#a7fe564c1d629e4c33ede4a20e3f2322d",
"class_f_e_a_t_1_1_geometry_1_1_index_tree.html#a7ab63e7f3af80116f534b34febc71bad",
"class_f_e_a_t_1_1_geometry_1_1_mesh_file_reader.html#ad2942d03a0369e27b8aa1a6ebdd6be63",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part_node.html#a322a377ee1064468c89070ea4ebb9b40",
"class_f_e_a_t_1_1_geometry_1_1_parti_par_m_e_t_i_s.html#ab4950dad3f186a43fde1f873013bdaeb",
"class_f_e_a_t_1_1_geometry_1_1_patch_part_map_holder.html#a1d39b66076754ba56cc9a1430170556f",
"class_f_e_a_t_1_1_geometry_1_1_standard_refinery_3_01_conformal_mesh_3_01_shape___00_01num__coor2b7a1f41c741938e6ddb9e7bf2c48774.html#aac4290103efbb74004a611e3b835d971",
"class_f_e_a_t_1_1_geometry_1_1_voxel_chart_masker.html#a27587a52d8e2b839df1edc71945b208d",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_3_01_l_a_f_e_m_1_1_dense_vector_3_01_d_t___00_01_i_t_914574362a1d637068ac8fa9eb1a3de4.html#a60621ecf3a459c5e9c50da98205d7a9f",
"class_f_e_a_t_1_1_global_1_1_gate.html#ad5ba7d2b6d51f6cecd283514f97e8e9a",
"class_f_e_a_t_1_1_global_1_1_nonlinear_functional.html#a913f884417d42cb2e58ecd732c39b98e",
"class_f_e_a_t_1_1_global_1_1_symmetric_lumped_schur_matrix.html#aa36af6b21a5424b419478cb00c0e9e32",
"class_f_e_a_t_1_1_global_1_1_vector.html#aed028af8d8f1e0524a0b6b9b44d5c313",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#af178f7566cab501705635b5ed22949a3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a64459a40e457cca79923f4deb3567b31",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror.html#aef05a5c35f39f89394d1807ce281287a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_meta_matrix_test_base.html#aca7247c86c1a32944732b6e08060cb79",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_base2.html#a81d2644ebde6a1b0393788fb1ab960fc",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_diag_matrix.html#acbb2f151cc0e05495054a2cecdf80cd2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#a7fad1ab6c9ab67cba13694b217366748",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_saddle_point_matrix.html#ac9a4c080eb891fa38703cd1aac14667b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a8604f1e4fff3ef8c11af4608b8d1b380",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a87f57bf8016ec61114c08dfddf241356",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#ad252fd28aaf3d001e3b4d74be74fd5cc",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a1e818235b7f50a3f0802958c4ca5c664",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a0a94d17901d3eb78142e1d3d4993938b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#ad078e209bce38973280864a8c04af1ab",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_filter.html#acea06ab1f4c077cdf637e4fef6234a6a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a4a3fc278d145d9c1ec38522e1c8f501d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a750b1210fab448b6154bfb22fcb36d98",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_default.html",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a81f425a4409c3acf4c4b602df7c8dd4d",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html",
"class_f_e_a_t_1_1_runtime.html#a831515a57ea208a209c3086dd2fe5621",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a7d583ff48c7678f1cffa30fdd48315a7",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a87caaf864220f9cb8add2b0eb90e2ad2",
"class_f_e_a_t_1_1_solver_1_1_analytic_function_operator.html#a974b96a786422cb4afc4704097fe24dd",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#afaab86a331f64d630f134f3c6027b5f7",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#a90679ca33476caf229d7f21d97270eba",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#afe5cc4723ba8581503c7b41a2e6893c8",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver.html",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#aae4b946889712e9e8afc0f27e7b6f12b",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#ae314aed268f854678d35272cd5e41340",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#ad29a0a21a3e91eaf8b46e12da46d0468",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#ac6850204b9181f10a1221bebcf5a1b6b",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#aada42ef9762a02842dc9f189725d3baf",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a92961d837ffc56d85cff3b820340eb2a",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#a162520eb548bc62b00ae058edee0bbd4",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a6bc6e53489ef8f0699399a419f0ec929",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a64b8b83c3cff4d963bea8b04d9f3800a",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a64b8b83c3cff4d963bea8b04d9f3800a",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a6e7c7ed985749618424e0373981a0c06",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a6950d64c2b9db6d1b76eff702a2bffa1",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#ac6850204b9181f10a1221bebcf5a1b6b",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a27c5e117b43720cbd96717472af0e517",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#a973d3942df838d68f63666265e73b025",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html#a010e5439714453aa3faf23b0b6d27a48",
"class_f_e_a_t_1_1_solver_1_1_polynomial_precond.html#a25e7d4394f26f6198f7b1a1711583bc9",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a118e9cd8518f424af8a625a95f14bb1e",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a6db06343bf3a13f014148a1e831c87cd",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#ac6850204b9181f10a1221bebcf5a1b6b",
"class_f_e_a_t_1_1_solver_1_1_s_o_r_precond.html#add306b84e9c14ef2cb491f105e63df35",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a3b1a9cd30df7bf5a6d7c899d171c954f",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#a418f269b40df7f5772a85fabfb2702f6",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#aa2d039f01373316c0e20b28699f9698a",
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_node_functional.html#af3b7b7f553156a6574c71c886f760a7c",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_82866b3cc738de4ae4935135fac5e732.html#ab70b19b6962179079f8a2d38ddddcd28",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_element.html#a8b224fb72515bc1bdc61dafccad1c546",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___001b011599bed78529504db50cb49456d.html#ac25e50167de479a49d61c8fa0ddbd179",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element.html#a8b224fb72515bc1bdc61dafccad1c546",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_base.html#a8b994c0cb6ba94e5475b21cea15768a0",
"class_f_e_a_t_1_1_space_1_1_element_base_1_1_dof_assignment.html",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#ad5886f4911cac95096215def29085034",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#a6db7e5c3954c49646d2823108329675d",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_element_1_1_dof_assignment.html#a219f4619cd071ad0f25275f84aecafb6",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html#a6bfdd412851516ebd9e42a508b51edb6",
"class_f_e_a_t_1_1_space_1_1_node_functional_base.html#ad42b3539fd9aefad4ad90a6f93b72fba",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_element.html#ae0abe6963c2356d0f7450af8544ef1ab",
"class_f_e_a_t_1_1_space_1_1_standard_vector_eval_traits.html#a97dfebf83204b64b549294d9b31e25df",
"class_f_e_a_t_1_1_syntax_error.html#ac922497e9e4a4a2e340e7173bb4d4780",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#a66eea93e9da4471bb685ae0583e922ef",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping_data.html#a0060219bcf4c278b71b47a31d22ad17c",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#ae94e874057ae9a215d3dac7dc8779771",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#ac4f5132a6a9b4f4b71976b86cd1bbee2",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_mapping.html#a08e99c74983c4fbbbd1d214be96ff1d1",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_mapping.html#ac494c99dfdc48ffe3733ef64b56ce9c8",
"class_f_e_a_t_1_1_xml_1_1_scanner.html#a809944eff677cc0d2d71eddb139e63a4",
"faq_page.html#faq_lafem",
"kernel_2lafem_2base_8hpp.html#ada0e93baecfc2c4933d9c5eb9c96a646ab6328d88ed1bc5247e0844c0b2ef7fe1",
"namespace_f_e_a_t_1_1_analytic.html#a81cd518a48bc53f93615ff4252c9e33b",
"namespace_f_e_a_t_1_1_pack.html#a7c69c0366e113fd42901ba9bcea55d4dadb03acfd680c5352bad66a6020fb77be",
"namespace_f_e_a_t_1_1_trafo_1_1_standard.html#a4753e53e2670882c652bacfe0b6b917d",
"resident_vs_transient.html#res_vs_tran_resident",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#a3f032f8afcbb1804abf4cbc22893c819",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html#ae1f798fcbfeff5dd3e9f0242f9ea1f49",
"struct_f_e_a_t_1_1_space_1_1_bernstein2_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercud0df8e3d6307173e0e27fffdaebdea90.html",
"struct_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_dof_traits_3_01_shape_1_1_hypercube_3_01dim___091_01_4_00_01dim___01_4.html",
"struct_f_e_a_t_1_1_type_1_1_helper.html",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_assembly_mapping_data.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';