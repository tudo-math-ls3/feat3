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
        [ "Mesh Terminology", "mesh_file_format.html#meshfile_basic_terminology", null ],
        [ "Mesh Topology", "mesh_file_format.html#meshfile_basic_topology", [
          [ "Triangle Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_tria", null ],
          [ "Quadrilateral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_quad", null ],
          [ "Tetrahedral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_tetra", null ],
          [ "Hexahedral Topology Ordering", "mesh_file_format.html#meshfile_shape_numbering_hexa", null ]
        ] ],
        [ "Boundary Parameterization by Charts and Mesh-Parts", "mesh_file_format.html#meshfile_parameterization", [
          [ "Boundary Adaption by Explicit Parameterization", "mesh_file_format.html#meshfile_parameterization_explicit", null ],
          [ "Boundary Adaption by Implicit Parameterization", "mesh_file_format.html#meshfile_parameterization_implicit", null ]
        ] ]
      ] ],
      [ "Overall Mesh File Structure", "mesh_file_format.html#meshfile_file_structure", [
        [ "Mesh Type Specification", "mesh_file_format.html#meshfile_meshtype", null ]
      ] ],
      [ "Mesh File Root XML Node Description", "mesh_file_format.html#meshfile_root", [
        [ "Examples", "mesh_file_format.html#meshfile_root_examples", null ]
      ] ],
      [ "Info XML Node Description", "mesh_file_format.html#meshfile_info", null ],
      [ "Chart XML Node Description", "mesh_file_format.html#meshfile_charts", [
        [ "2D Circle Chart Description", "mesh_file_format.html#meshfile_chart_circle", [
          [ "2D Circle Parameterization", "mesh_file_format.html#meshfile_chart_circle_param", null ],
          [ "2D Circle Chart Examples", "mesh_file_format.html#meshfile_chart_circle_examples", null ]
        ] ],
        [ "3D Sphere Chart Description", "mesh_file_format.html#meshfile_chart_sphere", [
          [ "3D Sphere Chart Examples", "mesh_file_format.html#meshfile_chart_sphere_examples", null ]
        ] ],
        [ "2D Bezier Chart Description", "mesh_file_format.html#meshfile_chart_bezier", [
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
      [ "Mesh XML Node Description", "mesh_file_format.html#meshfile_mesh", [
        [ "Vertices Sub-Node Description", "mesh_file_format.html#meshfile_mesh_vertices", null ],
        [ "Topology Sub-Node Description", "mesh_file_format.html#meshfile_mesh_topology", null ],
        [ "Mesh Examples", "mesh_file_format.html#meshfile_mesh_examples", [
          [ "2D Quadrilateral Unit-Square", "mesh_file_format.html#meshfile_mesh_examples_2d_unit_square_quad", null ],
          [ "2D Triangle Unit-Circle", "mesh_file_format.html#meshfile_mesh_examples_2d_unit_circle_tria", null ],
          [ "3D Hexahedral Unit-Cube", "mesh_file_format.html#meshfile_mesh_examples_3d_cube_hexa", null ]
        ] ]
      ] ],
      [ "MeshPart XML Node Description", "mesh_file_format.html#meshfile_meshpart", [
        [ "Mapping Sub-Node Description", "mesh_file_format.html#meshfile_meshpart_mapping", null ],
        [ "Topology Description", "mesh_file_format.html#meshfile_meshpart_topology", null ],
        [ "Attribute Description", "mesh_file_format.html#meshfile_meshpart_attribute", null ],
        [ "MeshPart Examples", "mesh_file_format.html#meshfile_meshpart_examples", [
          [ "2D Unit-Square Top Edge Mesh-Part", "mesh_file_format.html#meshfile_meshpart_example_unit_square_simple", null ],
          [ "2D Unit-Circle Boundary Mesh-Part", "mesh_file_format.html#meshfile_meshpart_example_unit_circle", null ]
        ] ]
      ] ],
      [ "Partition XML Node Description", "mesh_file_format.html#meshfile_partition", [
        [ "Patch Sub-Node Description", "mesh_file_format.html#meshfile_partition_patch", null ],
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
      [ "The mesh-extruder tool", "tools_meshtools.html#meshtools_mesh_extruder", null ],
      [ "The mesh-indexer tool", "tools_meshtools.html#meshtools_mesh_indexer", null ],
      [ "The mesh-partitioner tool", "tools_meshtools.html#meshtools_mesh_partitioner", null ],
      [ "The mesh-validator tool", "tools_meshtools.html#meshtools_mesh_validator", null ],
      [ "The mesh2vtk tool", "tools_meshtools.html#meshtools_mesh2vtk", null ],
      [ "The mesh2eps tool", "tools_meshtools.html#meshtools_mesh2eps", null ],
      [ "The mesh2svg tool", "tools_meshtools.html#meshtools_mesh2svg", null ],
      [ "The mesh2tri tool", "tools_meshtools.html#meshtools_mesh2tri", null ],
      [ "The tri2mesh tool", "tools_meshtools.html#meshtools_tri2mesh", null ]
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
"class_f_e_a_t_1_1_control_1_1_blocked_unit_filter_system_level.html#a0e7504a25cdfda8a0e5aba029faca89c",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#acaaef6f1e579cc13a14ee455e1f7774e",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#a4bc22ff565df9dddbbfa73feb4f66e2e",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#afbe934c80a678ad77ebf993f9294c21f",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_system_level.html#ae8147c5c2addd69b35bcc91e2d3c8775",
"class_f_e_a_t_1_1_control_1_1_scalar_combined_system_level.html#a1fd0d52031673e65105c9895088fbc0c",
"class_f_e_a_t_1_1_control_1_1_stokes_blocked_unit_velo_none_pres_system_level.html#a0cd49d826919fd160c3c8a4e98d14d35",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_driver_factory_3_01_driver___00_01true_01_4.html#ae2a54bd283d2777a1fd4893508b0db68",
"class_f_e_a_t_1_1_dist_1_1_comm.html#a3591c6db86690149cbe8e3a947c90260",
"class_f_e_a_t_1_1_dist_1_1_status.html#a696cc999dd1a915a6c132cd639e8ca9a",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier_chart_parser.html#a6b8983ffee7ad8aa0c59be15843bc44e",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#aecd136f6f28b0a470e384fa78b74725f",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_surface_mesh.html#a5a9175bfccb5c69b1f20f8a19f0ac2ca",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#ad45140dfb8189b440315564495f34754",
"class_f_e_a_t_1_1_geometry_1_1_hit_test_factory.html#ad14bf18180b25c14b29f22cbf2911b47",
"class_f_e_a_t_1_1_geometry_1_1_mesh_distortion.html#a9696dbbf4d51b42face6cdb1bf01abb4",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node_linker.html#a1b28d68de30faea078ac79e1d9337aa7",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#a97d3c100e1e6603b84703b15cc2b3044",
"class_f_e_a_t_1_1_geometry_1_1_patch_halo_splitter_3_01_conformal_mesh_3_01_shape___00_01num__coords___00_01_coord___01_4_01_4.html",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#a9f4626787d7443c58f2248869920db84",
"class_f_e_a_t_1_1_geometry_1_1_struct_index_set_1_1_image_iterator.html",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#a4dc6a4317c03622927df4516286a6914",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_matrix.html#a3d16b0cfd267be22c27d60c6944e9618",
"class_f_e_a_t_1_1_global_1_1_matrix.html#a87d7c73a96a31ef42d12f5f978437ff6",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar5e6f49c15ccb43840cb894a3375f33c3.html#a3a2fe09c4e1dffc3db7a471591c51871",
"class_f_e_a_t_1_1_global_1_1_synch_scalar_ticket.html#a5211f595224ce59b596f2ab36e8abfbe",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#a750b1210fab448b6154bfb22fcb36d98",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a3a65fd06d84b2c0b434d402fa5b96a54",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a958343c0926a9f67dffc5997c3c05b31",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#a7915be74ff0aba8d97b530cf87c663b6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_none_filter.html#a383e078c70d8629aa394993ffd2fd1d2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_pointstar_factory_f_d2.html#a8755a9bf0d5c9ec3665de94b37c1a6c5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#a2b89c7badcc9688c8a46670131874ed5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_row_matrix.html#add2555aaf437c71dfcbdf8c20de2c67e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_serial_config.html#ae97769e01d37ef8a6e803b398ddced75",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#aac3e00ad889c37831467302cee76e698",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#ab50cd3cf23b6f440797d57285ebba44c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded_1_1_scatter_axpy.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a42939644a39228d491e172cb109a5f73",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a6f3b4a5abd557187e9af23019968e8b7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_test_matrix_factory.html#a533a8813362f3e577c77826144b9028e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix.html#a5309f0e8960f26c0bef7ea9f71b68cb4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#aa809b8d8f27e23aa3cda04c6f6dfde9b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a9dae85fd4d7b9a2768e6cedcdb6859da",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_pow_of_dist.html#af29521425310bef3eb58f5998081e5f2",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#aa34167f0fd276303de1cf481031fce4e",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#aaa2d7e8a7dfb2f512db0b264054de133",
"class_f_e_a_t_1_1_simple_arg_parser.html#a99f09de8dec9e2cd253508bf975459b8",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#aa7281f8e3498245d411c5829b9d0cec0",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_approximate_hessian_precond.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a27c5e117b43720cbd96717472af0e517",
"class_f_e_a_t_1_1_solver_1_1_c_u_d_s_s.html#a351361de75fc52fecc0cb9b28003f22a",
"class_f_e_a_t_1_1_solver_1_1_convert_precond_3_01_global_1_1_vector_3_01_local_vector_outer___00120b74f43424364c71ee3ecd9deab297.html#a7518e1a16896f5a40f973533839e3257",
"class_f_e_a_t_1_1_solver_1_1_diagonal_precond.html#abc01448f27932c45beb3d53c177dd94b",
"class_f_e_a_t_1_1_solver_1_1_euclid_precond.html#a0b4002ca7fea449e5e7184111d32fced",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#ad8026d89fbcbd8d65ffdc207d15b06f6",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#ae6fc9035510853c34e062754732c17ba",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html",
"class_f_e_a_t_1_1_solver_1_1_hypre_solver_base.html#ad91ccd9edb259b13c5e776b77e3046a3",
"class_f_e_a_t_1_1_solver_1_1_i_l_u_precond_with_backend_3_01_preferred_backend_1_1generic_00_01_aed6c8770199c660f2b23c3052158c26.html#ab2b24531991e7d8f95b6621e4e4241f7",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a116b3f6b1736900d75d29aa0d7f644df",
"class_f_e_a_t_1_1_solver_1_1_m_k_l_d_s_s.html#acd83f64bd02b2ebed8fd32e62a4dd44e",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#aec76c75b2d2582ef2cf0ed928b583314",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_level_base.html#ade70c77e1a4b2e512d7248cfa85920d4",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#ac0e9f4bfffc09a8c2f55e5b2064dc096",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#ad29a0a21a3e91eaf8b46e12da46d0468",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#ac6850204b9181f10a1221bebcf5a1b6b",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#ad0e0aa027e766929e7e741f6384f887a",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#a50df7223074bef35a5a7cbd0c8ec9bec",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#ab562b0e476df821799ffd4151bf146bb",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a27c5e117b43720cbd96717472af0e517",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#a44a696c771d3d261646f8de9a353cf38",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#a43bbf72e681f56ef9a6a96af37265a76",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a9e960432ac8d4397fdcf7a19f83ddc95",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#ae314aed268f854678d35272cd5e41340",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#a50df7223074bef35a5a7cbd0c8ec9bec",
"class_f_e_a_t_1_1_solver_1_1_saddle_umfpack_mean.html#a978ff9c1e2475cea985f898a5e7a2cd4",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#aa5524e82d3d153e57af1787421c0ec90",
"class_f_e_a_t_1_1_solver_1_1_umfpack_mean.html#a62488b56e9de5e61911f25258bfa3785",
"class_f_e_a_t_1_1_solver_1_1_voxel_ama_vanka.html#a9b9d1846417be83357859c2763c89a2f",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_element_1_1_evaluator.html#ab3bbcda7a1f3dafe9be57d4e877e6b27",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_element.html#adc5863102651ecfd7291ac3bf23e0310",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_evaluator_3_01_space___00_01_trafo_evaluator_cd584e2ae53d2b4306ff680d7265a0ec.html#ad5886f4911cac95096215def29085034",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___0b6ac4939673a6cd90495f3edd5624cc3.html#a2644f5cb97bfdc7ce8f437a8c9176b4b",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_b2c6b42d6f9b9e25d469d3575a373a61.html#a1f5b8556e9fa1ee98993eeebeb5337c3",
"class_f_e_a_t_1_1_space_1_1_dof_mapping_base.html#afeedd71694157c54c9e09ca0dbdc38d5",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_element.html#aa8091b8119167113a57d2848d7f165bf",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_node_functional_3_01_space___00_01_shape_1_1_simplex_3_1390b0d61411af5d48632709957ef834.html#a12dd52ec74b75aff5d9471e1c3572485",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_se3129d8853eb9690234f254237b7f1ec.html#a3d51ed7e77717aee822fc7d91ad10b1c",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercubedbe93bfdbe7bb2d3cee5e3353a45cf73.html#aa6cb95c632710f9e50e675a723bc5532",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sb3aa60112044c986d3e79f1b1403fdf7.html#a2644f5cb97bfdc7ce8f437a8c9176b4b",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s88a2ce84347f190b629ca7916f4787f8.html#a24e648bd29f80680fe2028b3fe286673",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_2be5aeec58e5b7d19e38f7fb4e901a36.html#ad668fc46853ed1f0982a26cbf224a392",
"class_f_e_a_t_1_1_statistics.html#ac11e99466627cbc9ac812b614d59c996",
"class_f_e_a_t_1_1_time_stamp.html#a53d3d4fba86b5dd28edc23c21e390004",
"class_f_e_a_t_1_1_tiny_1_1_vector.html#a6c6e547f5e75a739ff28efc736fc8386",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_011_01_4_01_4.html#a633cb22af8831f64b258e86b2b8439f7",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_012_01_4_01_4.html#aa29419993e89dcf964bcedda2e00ea8f",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_6f4e5b7f5f9d6e2d04743ee3976d1657.html#a997be4a421fa226b409f5746b9bf018f",
"class_f_e_a_t_1_1_trafo_1_1_mapping_base.html#a8a1aa58ee0299d3d2b4da01d10021eeb",
"class_f_e_a_t_1_1_voxel_assembly_1_1_voxel_burgers_assembler_3_01_q2_standard_hyper_cube_3_01dimcd264f458d76378c7602bad0426c2704.html#a4d7ca863c6ce5eda30dc70c27c319746",
"classes.html",
"functions_func.html",
"macro__index__mapping_8hpp_source.html",
"namespace_f_e_a_t_1_1_dist.html#a10e8568d8af174896fa7c0553b9ab9e3",
"namespace_f_e_a_t_1_1_solver.html#a289617fe7cdace66d3123dc0472c31d7aa28d8f7ff8168d4aa740054357bcac61",
"namespace_f_e_a_t_1_1_voxel_assembly_1_1_kernel.html#acd63fb3587dd8f322e5795e1f8b91b66",
"scale_8cu_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#acfc2cb61e8367a0961e9991ee9ad435c",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html#aa7700e2529d379aa801194002d0600c3",
"struct_f_e_a_t_1_1_space_1_1_bernstein2_1_1_dof_traits_3_01_shape___00_010_01_4.html",
"struct_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_dof_traits_3_01_shape_1_1_hypercube_3_012_01_4_00_012_01_4.html#a9e1293bc5b64f1cb221b415e8537b6f6",
"struct_f_e_a_t_1_1_trafo_1_1_standard_eval_policy.html#ae017b397bfe82014498177839b71e1e3",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_assembly_cubature_data.html#a374011fc495dc679a9e813d632700cc4"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';