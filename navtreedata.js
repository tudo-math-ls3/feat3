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
    [ "FEAT Domain Controllers", "domain_controllers.html", [
      [ "Concepts", "domain_controllers.html#Concepts", [
        [ "Levels", "domain_controllers.html#domain_controllers_concept_levels", null ],
        [ "Desired Levels", "domain_controllers.html#domain_controllers_concept_desired_levels", null ],
        [ "Domain Layers", "domain_controllers.html#domain_controllers_concept_domain_layers", null ],
        [ "Ancestors", "domain_controllers.html#domain_controllers_concept_ancestors", null ]
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
    [ "FEATS's Thirdparty Dependencies", "feat_tpls.html", [
      [ "Reference", "feat_tpls.html#Reference", null ],
      [ "Enabling thirdparty dependencies", "feat_tpls.html#enabling_tpls", null ],
      [ "Resolving thirdparty dependencies", "feat_tpls.html#including_tpls", [
        [ "Finding thirdparty dependencies", "feat_tpls.html#including_tpls_finding", null ],
        [ "Downloading thirdparty dependencies", "feat_tpls.html#including_tpls_dowloading", null ]
      ] ],
      [ "For FEAT developers", "feat_tpls.html#developing_feat", null ]
    ] ],
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
          [ "Special case: Module based enviroment", "feat_for_vscode.html#unix_vscode_greenhouse", null ]
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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor.html#aee128d3e2e26251079f76083fd4fef3b",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive.html#a5fa797cf4859e4bdf3e6dd8e53bd87e0",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_frankes_function.html#af0811293f12eaa6e925012d166359e27",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow.html#a95d3a8dfab0faf499ab4510c5c34b874",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_pres2_d_1_1_evaluator.html#aa01f2bf540ecc6dafac1839eb657f745",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_taylor_green_vortex_velo2_d_1_1_evaluator.html#a1eb375034ddde8778d98406e71725128",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_inverse_distance_function_1_1_evaluator.html#a1b624e7d52058a5599660145723d4a45",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#a57e30c00f18017dc3b7ee06c5659312a",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a3253004f1e3e3dba964b886a0d6c5ada",
"class_f_e_a_t_1_1_assembly_1_1_analytic_cell_projector.html#a5381f53b131d1e076234fc198fa69a76",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#a3352b79041bcacc5b6c5ac8a04042799",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_diagonal_assembly_job.html#ab2190d9b81db6b3c5aab56334587a70c",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#ade32a12249e4aab9f82d5405fd57ee65",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_assembly_task_base.html#a75f1fccb2290c861e63210141d95bbed",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_matrix_assembly_job_1_1_task.html#a85dc327b6dc8edd06341885a46bf749a",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_scalar_assembly_task_base.html#a2b0267985a827e2b3b30192674698902",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_scalar_vector_assembly_job.html#a3eba53ea815dfbe68d2f5a98e4048317",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_matrix_assembly_job.html#a368cb0f310e20dbcc3eea2ad509930e2",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#ab0d4b588e5198fd8569b812139cdfdb5",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_force_functional_1_1_evaluator.html",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_operator_blocked_1_1_evaluator.html#a4dfc1807680a2ab83c090f87c3283ce7",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a57321922a430a9a10c0876e162496edb",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_basic_matrix_task_c_r_t_p2.html#ade8ccf5ae2a6f3c9568ea2f83e67ed88",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_bilinear_operator_matrix_job2_1_1_task.html#af039394a0128e719bbc85ee625bd6f44",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_linear_functional_vector_job_1_1_task.html#a934ab21ce749c35bb84e5740118bb61a",
"class_f_e_a_t_1_1_assembly_1_1_normal_gradient_surface_integrator_job_1_1_task.html#ae896f403fef42ba2cff2c7cb5d70fc2d",
"class_f_e_a_t_1_1_assembly_1_1_surface_integrator.html#a2ae2ef1ea8f07288f7368495340ee868",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_matrix_task_c_r_t_p2.html#a7e8d087d196933008bc197961eba21aa",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_vector_task_c_r_t_p.html#aa3d2a9ba34b6cb01567c801387d0f348",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_bilinear_operator_matrix_job2.html#a88577a0a637b1fe8fae9771cc4fcad45",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_error_function_integral_job_1_1_task.html#a6c7400e8ad38b360cdac060c452cc2c2",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_matrix_task_c_r_t_p1.html#abe22ca72129ad05a0796a15f358ca0ea",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_linear_functional_vector_job_1_1_task.html#a7257d6635fe4d71e38c33b35edc6b269",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_stokes_vector_analysis_task_c_r_t_p.html#ae27032836f6864125e5d08270654143d",
"class_f_e_a_t_1_1_control_1_1_checkpoint_control.html#aa2a6250c0cf403150d9efb312612fe56",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_c_g_a_l_domain_control.html#a391cbd168f5c2f6cee2d3846d1e0c363",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#a6c4a8d80f62666149e4d8b767abc7623",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#a0e75137c41f89029edec1c51981f7e88",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#a6e66ba14a3b61ed7f6323c36770460fc",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_system_level.html#a3fe803a8ac42dd0f2bcf9c65e1eabb6c",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_quadratic_system_level.html#acbc8e6087e323a35d52be35512d4365c",
"class_f_e_a_t_1_1_control_1_1_stokes_blocked_combined_system_level.html",
"class_f_e_a_t_1_1_cubature_1_1_hammer_stroud_d5_driver_3_01_shape_1_1_simplex_3_013_01_4_01_4.html",
"class_f_e_a_t_1_1_cubature_1_1_symmetric_simplex_driver_3_01_shape_1_1_simplex_3_012_01_4_01_4.html#a8fe0b70f58b273631c40d512e4cd3235",
"class_f_e_a_t_1_1_dist_1_1_request_vector.html#a35aafb969d4823f685ba08706aa92ece",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_index_tuple.html#a24d2621605f732c00e0d3f88304e18d5",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier.html#a163cf8205ab5353bcc2ecd46b5d5186f",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_chart_c_r_t_p.html#a54c780e47bf75e90d5756485341cc3d8",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_sphere.html#a0d4420c71e8c294879edb637c199c457",
"class_f_e_a_t_1_1_geometry_1_1_boundary_factory.html#aedd82de3417d9d11ba64b34cf62dbd8e",
"class_f_e_a_t_1_1_geometry_1_1_distributed_umbrella_smoother.html#abea6edab0d9ce7fe38beb28b106c39c3",
"class_f_e_a_t_1_1_geometry_1_1_index_calculator.html#a75061832865c15a8af1cc498a2da028e",
"class_f_e_a_t_1_1_geometry_1_1_lambda_hit_test_factory_1_1_hit_test_lambda.html#aa909f26e683b9618548d3fecb04847c8",
"class_f_e_a_t_1_1_geometry_1_1_mesh_file_writer.html#a7f082bc4d67769847d3fed78a5781de1",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part_node.html#a322a377ee1064468c89070ea4ebb9b40",
"class_f_e_a_t_1_1_geometry_1_1_parti2_lvl_3_01_conformal_mesh_3_01_shape___00_01num__coords___00_01_coord___01_4_01_4.html#a591e185229cf86325ff6a64d337b26a7",
"class_f_e_a_t_1_1_geometry_1_1_patch_mesh_part_factory.html#affc44c293000e15ca9a3d08d6fdf2e5e",
"class_f_e_a_t_1_1_geometry_1_1_root_mesh_node.html#ae91a13f109ac8630b7a4b02db4e33155",
"class_f_e_a_t_1_1_geometry_1_1_struct_index_set.html#aeb84788f51c202f1317815cb9eb82544",
"class_f_e_a_t_1_1_geometry_1_1_template_builder.html#a2101fb55bb35a72612ebfa93b23eaf19",
"class_f_e_a_t_1_1_geometry_1_1_voxel_lambda_masker.html#a6d9c604b6d05d56b564888dd13b492c0",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti.html#a3643b92d210b374a2b0ac0593797f21a",
"class_f_e_a_t_1_1_global_1_1_gate.html#a6d9a1339ec04f39be6e26b54e57d68ed",
"class_f_e_a_t_1_1_global_1_1_nonlinear_functional.html#a1a4643d2f259640d6098e9038ca2b771",
"class_f_e_a_t_1_1_global_1_1_symmetric_lumped_schur_matrix.html#a082c691f5fe3a4b3ec35bcb9b821d0d3",
"class_f_e_a_t_1_1_global_1_1_vector.html#a72f61774aa21c22ca72d753605b4d8c2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a8668efa558ba27a530c9479b1bafd814",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#af2f4b1c57faf1b0069ea502eee3f4234",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_chain.html#afbef41b5e93781da0314747cb7bc3cec",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter_blocked.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#aa3256534af8268ec485d428eab08e845",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html#ae8cb31ae7cb4019e937648ec01c00a2e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_mirror.html#a1684dbd3fbda1d39402d247d46ef6295",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#adbc6c744b75e9a13d7266030627dbb1d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a1d62f7bf82d9902bcf71c00d06db48e5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a1e818235b7f50a3f0802958c4ca5c664",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a352f79d7473139c3fffca5d6672ada4c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#a81e65b6d916b6690e1a5b583c46846aa",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a9ff3a798b8e01ad5a6bcf83965046658",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#af2c12c4048b0df5249077b10d89602a4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_diag_matrix.html#a074b2fbf8efbfced7d43bb84e93713e5",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix_row.html#a70f92b5b0fe3d6fb1a9d52240e7694fa",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter.html#ad6353c5114e51971444f1774ed49e4a4",
"class_f_e_a_t_1_1_memory_pool.html#abea66f37f6a9a3f3f3be103d4f61500d",
"class_f_e_a_t_1_1_meshopt_1_1_du_dv_functional.html#a993e0d4bb4b7359ce29f0a1b6eda6706",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function.html#a12dfdb1aeed3fe45a12c15cb6b11d7db",
"class_f_e_a_t_1_1_meshopt_1_1_rumpf_functional_base.html#ad32c3a582003b5678ebcfecea805547a",
"class_f_e_a_t_1_1_solver_1_1_a_d_p_solver_base.html#a49bb14b2fdf83babbaf47fc72b2835d9",
"class_f_e_a_t_1_1_solver_1_1_a_d_p_solver_base_3_01_global_1_1_p_m_d_c_d_s_c_matrix_3_01_matrix_040b7e1caa6985ce9572fb7c7833a8b9.html#af2e2fe610ff0c3b63fe23b73ed0316b7",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#afaab86a331f64d630f134f3c6027b5f7",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka.html",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#a1f17cb192fab20b3591eafdafd7c5f6e",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_convert_precond.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_direct_sparse_solver.html#a0dc948cf98b8a0a5715a1844b91c9528",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#a8af6d1f99d37ad98fd2a0bf6695936f4",
"class_f_e_a_t_1_1_solver_1_1_euclid_precond.html#ad87c27c7d905efdeafb574306117ddb9",
"class_f_e_a_t_1_1_solver_1_1_f_r_o_sch_preconditioner.html#a07f5437d9e62afe0be9ead7bac2056ff",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a973d3942df838d68f63666265e73b025",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#ad0e0aa027e766929e7e741f6384f887a",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#af8e655f0cdc31c3b245587b58a14d5b3",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#aa5524e82d3d153e57af1787421c0ec90",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#a8fa021d1ba593ca3ac009c053192d6d4",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a92961d837ffc56d85cff3b820340eb2a",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a6e7c7ed985749618424e0373981a0c06",
"class_f_e_a_t_1_1_solver_1_1_multi_grid.html#ab33422280445f1b8c038045a0cc2bb5e",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a44a696c771d3d261646f8de9a353cf38",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#a3b1a9cd30df7bf5a6d7c899d171c954f",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a3b1a9cd30df7bf5a6d7c899d171c954f",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a44a696c771d3d261646f8de9a353cf38",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#a3b1a9cd30df7bf5a6d7c899d171c954f",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#aec87d011de4859420af0f06559fbef2e",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#a5e276e3ab4c84a834377997024bf746e",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#ad586a1b68477198d76f4e794eaa83b5d",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#a9920812266b8f49921c9b133dd1451b4",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#a973d3942df838d68f63666265e73b025",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#ae6fc9035510853c34e062754732c17ba",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#a44a696c771d3d261646f8de9a353cf38",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#aa5524e82d3d153e57af1787421c0ec90",
"class_f_e_a_t_1_1_solver_1_1_schwarz_precond_3_01_global_1_1_vector_3_01_local_vector___00_01_mic9e5dece21cb7735a0bc194fb3569eb4.html#a3462cb641a929c772c65a0e45f04bef7",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#adb7b00a6751c00889ecb58ebe33f5060",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#a9002bc7d51bcd73fa2322b7a5c9c5bbb",
"class_f_e_a_t_1_1_solver_1_1_umfpack_mean.html#afe8803e26cff48158a7e2df71067f4ed",
"class_f_e_a_t_1_1_solver_1_1_voxel_ama_vanka.html#aef1655cdd37337774e68a34df4aa2d1e",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_045e90cd69ef865eaf5b37b181a99332.html#a29bafa5636409dda3697334dc98bdf3c",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_evaluator_3_01_space___00_01_trafo_evaluator__802e55c907078015a390ad2e449f9e5b.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_node_functional_3_01_space___00_01_shape_1_1_19e50c1472d62363091f828c6d82ab47.html#a12dd52ec74b75aff5d9471e1c3572485",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___0b6ac4939673a6cd90495f3edd5624cc3.html#a7bdda3b50a806b322f4df0ac22b98ba2",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_b2c6b42d6f9b9e25d469d3575a373a61.html#a6db7e5c3954c49646d2823108329675d",
"class_f_e_a_t_1_1_space_1_1_dof_mapping_renderer.html",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_element_1_1_evaluator.html#a0de51bcd1808f8e5465be06791b7c1a3",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_element.html#a3712ee7d71e34ce9d5824ccbaa2ee454",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_se3129d8853eb9690234f254237b7f1ec.html#ad05dbe2dde335dc94ba1a3b1830732ed",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_element.html#a7fbbea2698e84601deaa163443f8f955",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sb3aa60112044c986d3e79f1b1403fdf7.html#ac52f13660e702f543f09ae3a981ad0ea",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s88a2ce84347f190b629ca7916f4787f8.html#a9c746d797c0a8885c685170cec3fc250",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_3da8886277cafb8343dcd8258f6373f1.html#a36ca20d3ee0823eaaf0b82249752ce0c",
"class_f_e_a_t_1_1_statistics.html#af17c4e5e638826881ceed5a0a7882e2b",
"class_f_e_a_t_1_1_tiny_1_1_matrix.html#a0d869085110b7b48ad1de52dc1957ee5",
"class_f_e_a_t_1_1_tiny_1_1_vector.html#a844d999707fb464cf4769ee68a81e4d1",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_011_01_4_01_4.html#a661faef97931ce2f1f1a40701971af1c",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_012_01_4_01_4.html#aa688b36311b83c6bc1511e8ba33e5999",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_6f4e5b7f5f9d6e2d04743ee3976d1657.html#ab6a57d7b1a9e8030f1e047c87f32800b",
"class_f_e_a_t_1_1_trafo_1_1_mapping_base.html#ac7b22c127901c1d733ec33ad6618b3dd",
"class_f_e_a_t_1_1_voxel_assembly_1_1_voxel_burgers_assembler_3_01_q2_standard_hyper_cube_3_01dimcd264f458d76378c7602bad0426c2704.html",
"class_f_e_a_t_1_1_xml_1_1_scanner.html#a40d656eb15388e7bd7277c36f2afb8c6",
"empty__voxel_8cpp_source.html",
"kernel_2lafem_2base_8hpp.html#a82e983d75b2911bc2cd014a541d50711a6996a13fc5a9c8c12e71a3650cd5468c",
"namespace_f_e_a_t.html#abc175499504e9be922533b3a3627ccaca2a09e33185769aa027ad299c564da81b",
"namespace_f_e_a_t_1_1_l_a_f_e_m.html#ada0e93baecfc2c4933d9c5eb9c96a646a34b9b130e68b0a9ab2718496480e3142",
"namespace_f_e_a_t_1_1_solver.html#a911026db26857a46df6a675d2b53258b",
"node__functional__base_8hpp_source.html",
"sparse__matrix__bwrappedcsr_8hpp_source.html",
"struct_f_e_a_t_1_1_control_1_1_stokes_blocked_unit_velo_mean_pres_system_level.html",
"struct_f_e_a_t_1_1_geometry_1_1_intern_1_1_foundation_topology_collector_3_01_adaptive_mesh_type266448c054e9a95e3a2c44088bffa35d.html#aee8bebf41a5e776854fdc03eca792d79",
"struct_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_gather_scatter_helper_3_01_space___00_01_d_t___00_01_i_tf417dc54e181ee7869a07811099e7e1e.html#a3e382bb1c99805b03a282c0adcd6bb5b",
"struct_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___f7f932e34cb694f42be1d9bf6f7db3ce.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___000b42a86895f896d40a64f335ad2bfb25.html#a7b38c6399b62a66d8dbf9235d5182b7f",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01float_01_4.html#a265ee5643518ccf42303d6cded25c4b9",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_space_helper_3_01_q2_standard_f_e_3_01_shape___01_4_00_01_d_t___00_01_i_t___01_4.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';