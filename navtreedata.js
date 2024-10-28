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
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html#abff6abc7e9840424ee6f105fa015dee2",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_guermond_stokes_sol_1_1_evaluator.html#a104669d0db96ced7ab4fa586bf64c474",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_polynomial_function1_d.html#a422141600866acac3c5f9af315808507",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_r_h_s2_d_1_1_evaluator.html#ae48290222fc151d423fbfc54fced31ae",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_tensor_static.html#ab73266d57c7c4a8d2fed79f372d07f6e",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_plane_distance_function_s_d_1_1_evaluator.html#a422f17a859a529c4eb5ff834d901ed7c",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d_1_1_lambda_set.html#a9ce4e69a127dbe7e6351f71e2a525462",
"class_f_e_a_t_1_1_analytic_1_1_polar_coordinate.html",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits2.html#aa674f7e79397b200282bcd7432dae8a4",
"class_f_e_a_t_1_1_assembly_1_1_bilinear_operator_assembler.html#ab5f5e0dbc00851cb255e4132370daa84",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_task_base.html#a3c4b8515f1afe1b4036d12169ced1f38",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a736cba992ff363a87a3e3d73ed8bc6af",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a763134edee03d0a3ea34bac0d85c42d8",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#a0ff2d3581984c2e3b14fc96f1e9a7e13",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator_blocked_1_1_evaluator.html#a17f47f6a653121b2456f20f8f12943d4",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_operator_1_1_evaluator.html#a3fd6874777715b1cb43f0b67156269ba",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a1925bec919cdbc7be597f9a9b94fda9a",
"class_f_e_a_t_1_1_assembly_1_1_force_functional_assembly_job_1_1_task.html#a6381ea3c134621fb49071ff8a51cf36b",
"class_f_e_a_t_1_1_assembly_1_1_matrix_discrete_eval_data.html#a18a08309c54dc2d2b830f50ef846f7b3",
"class_f_e_a_t_1_1_assembly_1_1_symbolic_assembler.html#a863e661c8bf2c8ee80d3024c144d9a7f",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#a11c869f6821e9da2134cf948b1de3064",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#a28d9ea65fdfe9eed7a13a966049fb347",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_quadratic_system_level.html#a856a6f10d7a1b2cb5aa7bc3319c005c0",
"class_f_e_a_t_1_1_cubature_1_1_driver_factory_3_01_driver___00_01_shape___00_01true_01_4.html#ab3a34dc256e97bfab0ed68ae4a320fd0",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_rule.html#a8cddb8d2e5cd54a166bdfba97b458692",
"class_f_e_a_t_1_1_dist_1_1_comm.html#aea790a032fdd274e62957fb36b08c556",
"class_f_e_a_t_1_1_file_not_found.html#ac922497e9e4a4a2e340e7173bb4d4780",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_chart_c_r_t_p.html#a18b0d28b8e5f862af606fb34bcd5d78e",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_extrude_chart_parser_3_01_mesh___00_01true_01_4.html#a554dbe9f4ff5c3727d27da2516c97d56",
"class_f_e_a_t_1_1_geometry_1_1_boundary_factory.html#ab5c16bbd70f81929b027d54bd045f142",
"class_f_e_a_t_1_1_geometry_1_1_facet_flipper.html",
"class_f_e_a_t_1_1_geometry_1_1_intern_1_1_parti_iterative_individual.html#aa9308c72415e035cc04c22b92ced9006",
"class_f_e_a_t_1_1_geometry_1_1_mesh_file_writer.html#a7f082bc4d67769847d3fed78a5781de1",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part_node.html#aa68cab635e509a13b81d14dbcd69357c",
"class_f_e_a_t_1_1_geometry_1_1_parti_zoltan.html#ac309aabb25b2ed65e473f0b5c9e4885f",
"class_f_e_a_t_1_1_geometry_1_1_polyline_factory.html#adbd5d861df36c60bad67c7fdc6fdfa37",
"class_f_e_a_t_1_1_geometry_1_1_standard_refinery_3_01_mesh_part_3_01_parent_mesh___01_4_01_4.html#a84ce9c8de04f1f03655fc196b892804c",
"class_f_e_a_t_1_1_geometry_1_1_voxel_lambda_masker.html#a31a02cc669ecb0d2ac565881686a1382",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_3_01_l_a_f_e_m_1_1_dense_vector_3_01_d_t___00_01_i_t_914574362a1d637068ac8fa9eb1a3de4.html#ad5d9d6d8df53b839018342529e829acd",
"class_f_e_a_t_1_1_global_1_1_matrix.html#a25d5d0dfb026410f24421dd2757ae72d",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix.html",
"class_f_e_a_t_1_1_global_1_1_synch_matrix.html#a3891c9e7e6fd753412335716d3e53edc",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#a2a8542d823f9e80ee431e1cdcda6ea4c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a1bfb1a46a75ac2c11fd6a6b16d157abb",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a85f4140dfe65a2caa7fa30684e1784b1",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#a692b56906ddae814b44ac56e3ff3dbd5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_meta_vector_test_base.html#aca7247c86c1a32944732b6e08060cb79",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_d2.html#a1a2089d7024c143071c5cffbc7e2c208",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#a06f59965c5178a0e1feac7c1cb9f15cb",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#acadcd429015ba7db1b1fc487730ba75f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_serial_config.html#a0bcb4fbcd7f06f18671b1ca447d8551e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#aa0a3167fe55210236f2229a0da466fa3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#aa60b2a7f663f102b2fcda559d448176d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#af6adb2d2d65dd0333357ad1f93058c48",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a332cf155ce040c3dade8ef73b0a96335",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a47b718987b62807efb0461015eb18b56",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_test_matrix_factory.html#a015ef8ebd71416585d6c5dae8801cd05",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix.html#a3947607b17928861b1ed6aa7ef6046fa",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a93ec9424f64afc350b9163f16892c5b3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a958343c0926a9f67dffc5997c3c05b31",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_pow_of_dist.html#ac42b94423a27803cca0c353e235abc0a",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a9dcb9147addb6798e41d8ac2226c2dab",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#a854acd8e812449adc1e7fa7216c41d80",
"class_f_e_a_t_1_1_simple_arg_parser.html#a92b99ecd71d2ff67374889be651ffb69",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#aa5524e82d3d153e57af1787421c0ec90",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#aa86287c6d142b78ac494420116688d66",
"class_f_e_a_t_1_1_solver_1_1_approximate_hessian_precond.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a27c5e117b43720cbd96717472af0e517",
"class_f_e_a_t_1_1_solver_1_1_c_u_d_s_s.html#a351361de75fc52fecc0cb9b28003f22a",
"class_f_e_a_t_1_1_solver_1_1_convert_precond_3_01_global_1_1_vector_3_01_local_vector_outer___00120b74f43424364c71ee3ecd9deab297.html#a7518e1a16896f5a40f973533839e3257",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_solver_3_01_global_1_1_matrix_3_01_local_matrix___00_12c7f8d5d38cd40f20a5510074419a06.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#a7d583ff48c7678f1cffa30fdd48315a7",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#aa7281f8e3498245d411c5829b9d0cec0",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#ad5a06be09b1b40ab2c5b7cceb8274ab8",
"class_f_e_a_t_1_1_solver_1_1_hessian_precond.html#a466d07fad209727943213dcc5e41f2bb",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#aec87d011de4859420af0f06559fbef2e",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#af2ea6e8c96de45c8408f3dcc59b8f186",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#ad355136ba5cd50b1f2ece6658ab1c0c8",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#aab60dec26fb66c0cd2faeac189467f8f",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#a8ed3f4332277b2340972eaa65f8208a8",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a82061dfdc0313a9da082736c7c48559e",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a89a10069d43959492baab771bd0f183e",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a83c8455eb38065b8af5225a6bd94e899",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a92961d837ffc56d85cff3b820340eb2a",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#aed3c11b1edebf5eadd830e2912ddc941",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#ac93b8a3ec65386e4ba33600ff83c0139",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html#a8c2248c533b58422b9cb1316f349f9f8",
"class_f_e_a_t_1_1_solver_1_1_polynomial_precond.html#ac1fb944402f5fda01412384225e5c180",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a3d2a6e135e4af6af6be0c0aaf3330b3a",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a9920812266b8f49921c9b133dd1451b4",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#aeeb10e81625c18084050c7cf00d25330",
"class_f_e_a_t_1_1_solver_1_1_s_s_o_r_precond.html#a2c7e478addce5d2430d5db30beee8457",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a57e3c7c9e141643dd81fbb6f7bc439e4",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#ad546dfcee26f58c96b032be0955afe8d",
"class_f_e_a_t_1_1_solver_1_1_vanka_factor_error.html#a76bcfc0b1c2ed8cc4d368d67c2c5c193",
"class_f_e_a_t_1_1_space_1_1_basis_data_reduced.html",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_node_functional.html#af3b7b7f553156a6574c71c886f760a7c",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_element_1_1_evaluator.html#a52a02bca8404815741d03242d4074973",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___03bb1a3d92225206453ee71c00db5aa7a.html#a29bafa5636409dda3697334dc98bdf3c",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element_1_1_dof_assignment.html#a05cce9327e135fc6b29af8788d2cfe58",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_single_entity.html",
"class_f_e_a_t_1_1_space_1_1_evaluator_base.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spd08b4bf16ff97c315f2f45da544757eb.html#a79d664e46c384d5808e89e8ae2b9905f",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sbede9c86e9ef3e51ccbac5494ff49399.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_evaluator.html#a40dd6a3cebe11cda4d26fe1ad546e58d",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s5e9945e9db3d72bec06849dc7e957367.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_element.html#a3d583a93a3354fcf41736c3c98b78e46",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_2be5aeec58e5b7d19e38f7fb4e901a36.html#a043a5a49d13cf6c3ac076cd0fb29e80a",
"class_f_e_a_t_1_1_statistics.html#a4599690ce805dae631a95f70ec79e0d4",
"class_f_e_a_t_1_1_test_system_1_1_test_list.html#ae6fa4e251f3031a3bd0f588073233a09",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#ad158c076870ba63056d030289e6b0505",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator.html",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_012_01_4_01_4.html#a67774a6c2ae96e0ff58eca6c1756e7f0",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_6f4e5b7f5f9d6e2d04743ee3976d1657.html#a4220f98b9a2542aa5867d622a405cc15",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_mapping.html#ae9d1a5896315d62563438ef7dcb1bd06",
"class_f_e_a_t_1_1_voxel_assembly_1_1_lagrange_data_handler.html#a9ead7d9de90cb61f7d6bce2d5b572a26",
"class_f_e_a_t_1_1_xml_1_1_scanner.html#afdf2a4cb9f2c35dfefe573f1aa3c18c9",
"feat_for_win_vs17.html#win_vs17_prereqs_vs",
"lagrange3_2node__functional_8hpp_source.html",
"namespace_f_e_a_t_1_1_assembly.html#a28186a335b7da11c0a4a2a3770462bf0add5c07036f2975ff4bce568b6511d3bc",
"namespace_f_e_a_t_1_1_solver.html#a25100c25aa3685c8392ae3cce25f9618",
"namespace_f_e_a_t_1_1_windows.html#ace71d32a51d5f170810bce006dca5be8",
"scale__row__col_8hpp_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#ae779b0255940ed5c448e41a501300eef",
"struct_f_e_a_t_1_1_l_a_f_e_m_1_1_arch_1_1_min_index.html",
"struct_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_dof_traits_3_01_shape_1_1_hypercube_3_01dim___091_01_4_00_01dim___01_4.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___000b42a86895f896d40a64f335ad2bfb25.html#a090cbd87546b2092af578d0187694e50",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01double_01_4.html#a1b54a0535a6d416f0d4805f20c98acf3",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_space_helper_3_01_q2_standard_f_e_3_01_shape___01_4_00_01_d_t___00_01_i_t___01_4.html#a1cd506a341a668ab6a0e67a8703b15af"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';