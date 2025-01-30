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
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol.html#a06c11ef8fa0e6f36e17a21427272396e",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow_1_1_evaluator.html#a1eb375034ddde8778d98406e71725128",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_r_h_s2_d.html#abc7769f66663a2d02cec7592309456ab",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_tensor_static.html#a42497bb85b16056174e7ad693b3b0f2f",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_plane_distance_function_s_d.html#a19f9c0038e4ca37c510b5b8a6bb6970a",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d_1_1_lambda_set.html",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function_1_1_evaluator.html",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits2.html#a554ef7f691fb352312a385e19f9c55b1",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_1_1_evaluator.html#a3fd6874777715b1cb43f0b67156269ba",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_task_base.html",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a3c4b8515f1afe1b4036d12169ced1f38",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a3c4b8515f1afe1b4036d12169ced1f38",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job.html#adae2e75ea23d8fce5804558721dd7465",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator_1_1_evaluator.html#ad0521bf4b7169457b4b75249ce3fe46f",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_functional_1_1_evaluator.html#ac6dcdb4b3b28a8e517f3a49eba70474e",
"class_f_e_a_t_1_1_assembly_1_1_discrete_function_integral_job_1_1_task.html#af5daa5c60b1f9e91053bcf01d23cb104",
"class_f_e_a_t_1_1_assembly_1_1_force_functional_assembly_job_1_1_task.html#a090631b96cfff3e6586518197463c5eb",
"class_f_e_a_t_1_1_assembly_1_1_linear_functional_assembly_job_1_1_task.html#ac5ec3b8f0f1c66dacbb1c1b18316b751",
"class_f_e_a_t_1_1_assembly_1_1_symbolic_assembler.html#a07779170c583ef1393c4ef18276d3906",
"class_f_e_a_t_1_1_control_1_1_checkpoint_control.html#acc5cdb3b26e223bd5ae6fd455bed5617",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#ac814c290ffdfd7512b61da651bd7c052",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_nonlinear_system_level.html#acc909b2745fa2a31aac81c92c1125d72",
"class_f_e_a_t_1_1_cubature_1_1_auto_alias.html",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_newton_cotes_closed_driver.html#a7eff73d30e34354a0be28a8bf00220ba",
"class_f_e_a_t_1_1_dist_1_1_comm.html#abecd25fafb90350cd5508393addaeff5",
"class_f_e_a_t_1_1_exception.html#a76bcfc0b1c2ed8cc4d368d67c2c5c193",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_chart_base.html#a1ebf2a280d49969ca750db40c6c0e59f",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_extrude.html#abf159301aa38e5bf4befc1108a7bae1b",
"class_f_e_a_t_1_1_geometry_1_1_attribute_set.html#a551f84d3cbd9b013b0d697ce862ae92b",
"class_f_e_a_t_1_1_geometry_1_1_export_v_t_k.html#a464aac7a1552bd3f2aae3e1016615265",
"class_f_e_a_t_1_1_geometry_1_1_index_tree.html#a11eacfe8f837c28079c3bf4dcfa296e9",
"class_f_e_a_t_1_1_geometry_1_1_mesh_file_reader.html#abef75c771edd48fe54f531fe8ccf1999",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part_extruder_factory_3_01_conformal_mesh_3_01_shape_1_1_hyp1e036ec37b875c8064101b4ff69fd188.html",
"class_f_e_a_t_1_1_geometry_1_1_parti_par_m_e_t_i_s.html#a45da9c1c6b4fd72979bb63bd6079cee6",
"class_f_e_a_t_1_1_geometry_1_1_patch_part_map.html#a87718ca3b6157225f37b254b0bdd4779",
"class_f_e_a_t_1_1_geometry_1_1_standard_refinery_3_01_conformal_mesh_3_01_shape___00_01num__coor2b7a1f41c741938e6ddb9e7bf2c48774.html#a7322d8937ba5b10713b6d6834ea65619",
"class_f_e_a_t_1_1_geometry_1_1_voxel_c_g_a_l_masker.html",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map_1_1_write_result.html",
"class_f_e_a_t_1_1_global_1_1_gate.html#a151f690bd1866d13f920f76e40332aba",
"class_f_e_a_t_1_1_global_1_1_muxer.html#ab5aeadc2994978aa8a82769a8a3c544f",
"class_f_e_a_t_1_1_global_1_1_splitter.html#a82f587cf3341c9985021d7f5e75623d3",
"class_f_e_a_t_1_1_global_1_1_vector.html#a2d611cb63a8d5690538a9c4e63236c87",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a9060d0e4e26b6ca72132ec9504a8ec90",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a02f9a68c88154f23ebe3943db9eeb12a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_sequence.html#a432edc213c3651d73c6758b0c31f2c32",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter_blocked.html#a2c5496e7cbee191de7aea9295325b5c2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#abd457828f842fa32fa4090d0b8bc5a36",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_diag_matrix.html#a16bea892bdd2e56b0a8c8551b0737b4d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_mirror.html#a83316953b075233629f14ca1802ff128",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_saddle_point_matrix.html#a0bfc7d866e1de15d94bdfb3fa7d19271",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a3d2a2536488a06809e3b435c0d8435a4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a3c9b732785c183c76eb5429fdae287c0",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a83fec5dca47080be3ee7d4e0fd9129bc",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#abe31b9ab867ff23f476d68b58d164b5d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#ae3720e3fd0f356b3d47a877490f93d5d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#a721d549ae9bd76b39f6926c3c3e6b924",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_diag_matrix.html#ab271e704fb8a18bbff8b045895d442c7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_mirror.html#a553a942a032ec40ad8644a815ad1d484",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter_blocked.html#ad7d3d256d2e68540badc384d5efa49fc",
"class_f_e_a_t_1_1_meshopt_1_1_chart_distance_function.html#a3bab4dd99f0ed337269283295273c651",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a1b345cec59df9ccfd4702a7efb8eaefc",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function.html#af4613cd50db05475fb923fb474b86b66",
"class_f_e_a_t_1_1_property_map.html#a8a3bb18a30ad3e366912c14882d9f0e3",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#a39abb0cee57537441aad338075911f8d",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#a4067f7a2a268e06d5f57f6097d14aaff",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka.html#ac72948c9353e0b61f53f2e2e03d7b27a",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#aae4b946889712e9e8afc0f27e7b6f12b",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#aae4b946889712e9e8afc0f27e7b6f12b",
"class_f_e_a_t_1_1_solver_1_1_descent.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_global_1_1_matrix_3_01_local_matrix___00_12c7f8d5d38cd40f20a5510074419a06.html#a28cfc6782adf797385af1c3b7fc4d75a",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a55827be045f60899c3c9d7cc9df4f310",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a7d583ff48c7678f1cffa30fdd48315a7",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#ad0e0aa027e766929e7e741f6384f887a",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#ac59b69cc5fe017dd82287034e6a25a5a",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a903b89f07e60974b07e6348291fc6741",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#a1536745b0a48dfa72ae7bc51e62e031b",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a69fab2263b0c35bf03aa77d95e3a5917",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a6e0074b3a336238364169fb1a0fb0701",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#ac59b69cc5fe017dd82287034e6a25a5a",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a26458def0ef8654e3fff8f444df8108f",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#a96d7c2134c0f82e5dae7f5058c480e83",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html",
"class_f_e_a_t_1_1_solver_1_1_polynomial_precond.html#a1b635d98d96ea088044b336132594c72",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a0f794ef8c4c58303491aaac0d37f64ec",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a69fab2263b0c35bf03aa77d95e3a5917",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#ac59b69cc5fe017dd82287034e6a25a5a",
"class_f_e_a_t_1_1_solver_1_1_s_o_r_precond.html#ad7085e5ff1f9040f07126f2ef7c333eb",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a3216d57242be2e62ed958a84be4ff2e5",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#a38b619f96dd684b78c67d9794081a0de",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#aa0dd471a743889d2357f7b57fa169dab",
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_node_functional.html#ad091f3ff7d0c9c412833b3a8446a43ff",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_82866b3cc738de4ae4935135fac5e732.html#ab532b05d5e1b0aaf9bbc18d3280f1c93",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_element.html#a7fbbea2698e84601deaa163443f8f955",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___001b011599bed78529504db50cb49456d.html#abeda000eff0270b2a43789f0148fb536",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element.html#a7fbbea2698e84601deaa163443f8f955",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_base.html#a8864b1ef8c6a9cba7a595f8e8f5565aa",
"class_f_e_a_t_1_1_space_1_1_element_base.html#af250c0bc81fe3e1e64cc45cd2bef2dc6",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#abeda000eff0270b2a43789f0148fb536",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#a6bfdd412851516ebd9e42a508b51edb6",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_element_1_1_dof_assignment.html",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html#a68cd1384c7d7438e90a63fbed4a22f0f",
"class_f_e_a_t_1_1_space_1_1_node_functional_base.html#ac36f33ce44f2589fa5e008245c42dcb4",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_element.html#adc5863102651ecfd7291ac3bf23e0310",
"class_f_e_a_t_1_1_space_1_1_standard_vector_eval_traits.html#a803a1f93bf2d379b6c2cd1fd9421489d",
"class_f_e_a_t_1_1_syntax_error.html#ac3251651933f53a6f97b68c9348205c0",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#a5eab64e8c5c5ecb9b0802a842a7de865",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping.html#a1f9df2846ff9e3cd54da4e4ad198e000",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#a62adbfc0c0a7f59f56429b7d4ba82347",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#a48d4f98f064dd604e9fa5e7f61962f8b",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#a997be4a421fa226b409f5746b9bf018f",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#acdb7bfb197ee64e4d83a5901a0aaa456",
"class_f_e_a_t_1_1_xml_1_1_markup_parser.html#a5d72a9406f0a480b890e3f8bf8768560",
"driver__base_8hpp_source.html",
"kahan__summation_8hpp_source.html",
"namespace_f_e_a_t.html#acf92a649a38de9336e758d9ef0880668",
"namespace_f_e_a_t_1_1_pack.html#a094c0e6ec330169ae1210a1cd0fbd8a6",
"namespace_f_e_a_t_1_1_space_1_1_discontinuous.html#a4847f978e054b47c5b104d7a288f0e36",
"preproc_macros.html#ppm_feat_have_omp",
"struct_f_e_a_t_1_1_assembly_1_1_lagrange1_to2_dof_at_entity.html#ab8cbdd086b330cf0a7f065c31798e597",
"struct_f_e_a_t_1_1_geometry_1_1_mesh_part_1_1_index_set.html#a348282f4f2220a6e2911ea1378dde2f9",
"struct_f_e_a_t_1_1_solver_1_1_optimization_test_traits.html#ab28b44ed4b7e5407369946317564c5d3",
"struct_f_e_a_t_1_1_space_1_1_lagrange3_1_1_node_functional_3_01_space___00_01_shape_1_1_vertex_00_01_data_type___01_4_1_1_value.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___00f78f81b86ee7ca50fe24a6c9ff685ed0.html#a22db587d2da3894c981dc5db8c656d80",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01long_01_4.html#a606d66292c853ed12d61578945e5d4ae"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';