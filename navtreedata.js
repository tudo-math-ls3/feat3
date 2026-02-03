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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor.html#afdb590ddb07f3922a0dc182f814be000",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive.html#a6494c65c7913f5be7467a7c23e6ac535",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_frankes_function_1_1_evaluator.html",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow.html#aa25084c8acbb4339b3e41fb82327c7dc",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_pres2_d_1_1_evaluator.html#ab9830b01eedcfb335bff302531347911",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_taylor_green_vortex_velo2_d_1_1_evaluator.html#a34c2a7c0b6667b9599fce727f0f1d682",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_inverse_distance_function_1_1_evaluator.html#a1eb375034ddde8778d98406e71725128",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#a5f075c3df417090a0586decbcc12cee0",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a356992efb6443bc6fbc8ca882eeb914d",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits1.html#a687b6ba2998c8fa94f9406ce472ac0fe",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#adf82041f4e71c6b7aad8da35d9928aa5",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a1050035b8c4a9db89fa42fdc7e21509d",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_assembly_job_base.html#a046d3da41aa78a1b4f91ac9b64b6ed30",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_matrix_assembly_job.html#a2fcd58686030ff2ce68f8615e9195767",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_vector_assembly_job_1_1_task.html#a9040fa67d8287413afc76c16ec5d1769",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_scalar_matrix_assembly_job_1_1_task.html#a7493db55e2547eb282b9f3fa4407e390",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#a4b0f648dd3f1c32f602154a4c73a8d5d",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job.html#adf82041f4e71c6b7aad8da35d9928aa5",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator.html#a210eb80ce04ad03aeb8a89c9a0b70ba1",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_functional.html#abc41268a9f01a7b6374db6cc038afbb8",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_trial_derivative_operator_1_1_evaluator.html#afebf92c3395d790ba46d157930b705b0",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_basic_matrix_task_c_r_t_p1.html#a7226c697bddf874056338903feaff3ef",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_bilinear_operator_matrix_job1_1_1_task.html#a806fb501b0edcb97175baa49889e1f6c",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_force_functional_vector_job_1_1_task.html#adefc3297625bbb2c8c97da472604d7c1",
"class_f_e_a_t_1_1_assembly_1_1_linear_functional_1_1_evaluator.html#a463dd4b6ff3c8763ba86a5c12188cfe1",
"class_f_e_a_t_1_1_assembly_1_1_stokes_f_b_m_assembler.html#a0a7627b02595acce735db1f93669059a",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_matrix_task_c_r_t_p1.html#a93b67f1701a21e293fc6666fe942dbb2",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_task_base2.html#ac677d94014929f5437028824f2aa587c",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_bilinear_operator_matrix_job1_1_1_task.html#a3ff6b62efe7a1c0be3df3fa4c2c1708d",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_discrete_function_integral_job_1_1_task.html#a950d578cb6bfdb43e5a6dcc789523605",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_mass_matrix_job_1_1_task.html#ae179891927dd8f8b9deb268df12013aa",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_task_base1.html#ab64f069c26eaaa29aa8c93fd8ff86ff9",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_stokes_vector_analysis_task_c_r_t_p.html#a03ac93098355b1b99e8e3f6cb6ae2bb3",
"class_f_e_a_t_1_1_control_1_1_blocked_combined_system_level.html#af6ea4a66e668984705765e36d0436a0b",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_adaptive_parti_domain_control.html#ade59618d6d10184b7a388518b8e35fc5",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control_base.html#a426bf2533f64c2ebafc4112b72626347",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#aad5134f3db2b8c9152c83b29d1f6852e",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_hyperelasticity_functional_control.html#a974365723f8cdc1c8f1cd45c087eecdc",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_nonlinear_system_level.html#a97be4a1a12940212366b7c8f2676e58d",
"class_f_e_a_t_1_1_control_1_1_scalar_mean_filter_system_level.html#a2c924648ff822a5ec7876caf48a881fe",
"class_f_e_a_t_1_1_control_1_1_time_1_1_nvs_bdf_q.html#aa92f4ae9caa234225fd5a38e86ae5137",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_midpoint_driver.html",
"class_f_e_a_t_1_1_dist_1_1_comm.html#a92a62cb6f962c57bb8581f05ce82a02e",
"class_f_e_a_t_1_1_dist_file_i_o.html#aa4db7630d032b16dbdab3c8429f4342d",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_mesh_layer.html#a135a88ca798d4f64c9ff8af42ce8e2ce",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier_params_parser.html",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle.html#af50eb65b88ae82113fef4dab2495b395",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_surface_mesh.html#a67da106a897714a0e7206a2542e47cd1",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#a70f17e1b619f947342801ff7cadfc0f4",
"class_f_e_a_t_1_1_geometry_1_1_factory_3_01_mesh_part_3_01_mesh_type___01_4_01_4.html#a673021963cdd4baf368c32a808c5d549",
"class_f_e_a_t_1_1_geometry_1_1_intern_1_1_adaptive_mesh_storage.html#af5c49786c9ab0bacc6692e6653a2961f",
"class_f_e_a_t_1_1_geometry_1_1_mesh_extruder_3_01_conformal_mesh_3_01_shape_1_1_hypercube_3_012_d55947478d0079b597f439c47758a098.html#a7cb944d0a26d165fed12a833198533c5",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node_linker_error.html#ac922497e9e4a4a2e340e7173bb4d4780",
"class_f_e_a_t_1_1_geometry_1_1_mesh_permutation.html#a3861274b7b941b7ac2430cd3e1ebc591",
"class_f_e_a_t_1_1_geometry_1_1_partition.html#a3efaa64e126e12a59431d4c86c8ae5b9",
"class_f_e_a_t_1_1_geometry_1_1_refine_factory.html",
"class_f_e_a_t_1_1_geometry_1_1_standard_refinery_3_01_conformal_mesh_3_01_shape___00_01num__coor2b7a1f41c741938e6ddb9e7bf2c48774.html#ae50273415e5cb711c8d90f84be71e244",
"class_f_e_a_t_1_1_geometry_1_1_sun_zhao_ma_data.html#a8a5acb2b8029ae2a224ff0fe3699893a",
"class_f_e_a_t_1_1_geometry_1_1_two_refinement_data.html#af9fcc093d7f62b5962ee5af2d1acde1f",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#abc5f56ec104dea837fd782d257a7b920",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_system.html#a739412d365ec976c141d05fb3787c7b0",
"class_f_e_a_t_1_1_global_1_1_mean_filter.html#a5b04fceb520fb1abbc3be65d70a15726",
"class_f_e_a_t_1_1_global_1_1_p_m_d_c_d_s_c_matrix_3_01_global_1_1_matrix_3_01_l_a_f_e_m_1_1_spar5e6f49c15ccb43840cb894a3375f33c3.html#a92d328fdc9f305ed4a6ac7f3852a953c",
"class_f_e_a_t_1_1_global_1_1_synch_vector_ticket.html#ad20a302f39c916b5e9eea8846b8f114f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_container.html#ada2b8bd276f75ab9df2a22cb52f48e7c",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#a74c9279f1115605cd6ac24be25d4bdef",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#ab63eb9b24066b1e740e33644b27bc98b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_matrix_mirror_buffer.html#aa0f8a8c78b7e7b226dd476bd9fe924d3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_none_filter_blocked.html#a20f3a97e509dfd5920e9fab9e41375b2",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_filter.html#aab073914b60167bb30556976c78f479b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#a1118dd4122f83683d9c3c81d2db4ac39",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_slip_filter.html#a1ed744d4084df2eb3b4b51eec1408f46",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#ac096a8e63c8eac202322503a7c0ce868",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#ac0974c13476ac8da3f1996ed9fc9cdce",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded_1_1_image_iterator.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a3d2a2536488a06809e3b435c0d8435a4",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a44a72403da9cdefad5efc0eeccfc1912",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#aff83f81f0f25c35fd0781518293d2ac1",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix.html#a2975d97f66c2b4ee705213b2559e781f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_vector.html#a78fd8cd1ae56f9741becd64929ae3c2e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#a83fec5dca47080be3ee7d4e0fd9129bc",
"class_f_e_a_t_1_1_meshopt_1_1_concentration_function_default.html#afd87139f5c9e709b306eba1a5fd9b957",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#a91f916b1d0b8100d7c78a83793ee4484",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_quality_functional.html#a5f8a84b7c961a397be7869957d1b7304",
"class_f_e_a_t_1_1_runtime_1_1_scope_guard.html#afc7d9951b96e9ea3edc2bfd48afd685b",
"class_f_e_a_t_1_1_solver_1_1_a_d_p_solver_base_3_01_global_1_1_p_m_d_c_d_s_c_matrix_3_01_matrix_040b7e1caa6985ce9572fb7c7833a8b9.html#a4da3d5390745d8a5a1746d831e58ef70",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#acae5a4593ac598b83f710be5e29d3a7b",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#ad2de4b053e54def8d954544128edf243",
"class_f_e_a_t_1_1_solver_1_1_b_f_b_t.html#a97dc8d73b0423dab5b3c3782ce70eb11",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a0ab9bef4ea22df32065a43399a990561",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#a8d8f0c00e06fd7b9f1d939209f5ddca0",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#ad5a06be09b1b40ab2c5b7cceb8274ab8",
"class_f_e_a_t_1_1_solver_1_1_descent.html#aeebdc6fb0b2e5072cb5c0c52205ff58b",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#a4b0f71cc9382c03863f4d4c682e71e95",
"class_f_e_a_t_1_1_solver_1_1_euclid_precond.html#a9002bc7d51bcd73fa2322b7a5c9c5bbb",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#add5c03ded9d5568e14ec60cdfb1572b8",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#ae856346c0f498690e441ceefe6edd212",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#a0810de389dc450bea5cb6917b66ca051",
"class_f_e_a_t_1_1_solver_1_1_hypre_solver_base.html#a9002bc7d51bcd73fa2322b7a5c9c5bbb",
"class_f_e_a_t_1_1_solver_1_1_i_l_u_precond.html#a158c52877dea6c80d942a6aecc3735b2",
"class_f_e_a_t_1_1_solver_1_1_jacobi_precond.html#a3684c5805312d526006e0dcda44f6aaa",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#ae6fc9035510853c34e062754732c17ba",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#abb072260d9458116ae6aa240d5685fb3",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_hierarchy.html#adada1c798d6725f2dd0f75acd5822671",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a97f6e79dbe58d08330cec6ba5c0195c3",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#aa2b2f7492bf2c6b7d9f740b95e99969f",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#a9920812266b8f49921c9b133dd1451b4",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#aa5524e82d3d153e57af1787421c0ec90",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#aa8620b8ca6e6cce3137dfc29ecf2f08e",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#a6fd0f61a9458e3b6c32b4c502fcdc7f6",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#add5c03ded9d5568e14ec60cdfb1572b8",
"class_f_e_a_t_1_1_solver_1_1_para_sails_precond.html#a8c34d07a55cd09eaafef49269f140536",
"class_f_e_a_t_1_1_solver_1_1_polynomial_precond.html#a25e7d4394f26f6198f7b1a1711583bc9",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a0f794ef8c4c58303491aaac0d37f64ec",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#a6db06343bf3a13f014148a1e831c87cd",
"class_f_e_a_t_1_1_solver_1_1_r_g_c_r.html#ac6850204b9181f10a1221bebcf5a1b6b",
"class_f_e_a_t_1_1_solver_1_1_s_o_r_precond.html#add306b84e9c14ef2cb491f105e63df35",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a3216d57242be2e62ed958a84be4ff2e5",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#a3cd993260b72e64a99fae533594e0e70",
"class_f_e_a_t_1_1_solver_1_1_vanka_3_01_l_a_f_e_m_1_1_saddle_point_matrix_3_01_matrix_a___00_01_9ced3bd73ff865688c51cd8e1ad52f3a.html#a110b7d665022321e217ce5731f88ff57",
"class_f_e_a_t_1_1_space_1_1_argyris_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spa96d47ea9782e73b0166b29cca748b7d3.html#a039d795545bbd7867c97e252a4ed7113",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_64897ed7485eba824a67fb42fd5ab9c0.html#ac3f96b407162f030ff788d1a18c8e775",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_evaluator_3_01_space___00_01_trafo_evaluator__cefb60ca526d73ea9fbcb2c4bda04e64.html#a43908672b4097e72cfd3666bcd87f206",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___001b011599bed78529504db50cb49456d.html",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___0d5fd0fd65b9232ba00b88a0e0458fe94.html#ae4d17e132c111387fbe169819fc1c2aa",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_d256e15bbc7b3e4885146b0f6ddde539.html#af1791eca466aa31e4b2f1caead8ede7e",
"class_f_e_a_t_1_1_space_1_1_dof_mapping_uniform.html#ac52617aa7ce9ae5aefc98bd862360280",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spb4a6a2373ae94132a2de09a3720ca169.html#a05bec126e2b37bf1c13ca9c419cb4eeb",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s189823cc2b39e0d1a920e7698210faed.html#a6fa689aeb431d8551cf2977d30edd661",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_node_functional.html#ad091f3ff7d0c9c412833b3a8446a43ff",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s17de294e509ec5be8be361ba5ec7318e.html#a6bfdd412851516ebd9e42a508b51edb6",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_node_functional.html#a36c9f543c5dc1c565378914c77ffd903",
"class_f_e_a_t_1_1_space_1_1_parametric_evaluator.html#a6db7e5c3954c49646d2823108329675d",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_node_functional_3_01_space___00_01_shape_1_1_hypercub9798db6ea14cc0d2104186ccce2f66d8.html",
"class_f_e_a_t_1_1_string.html#a9dde92f2031160322d96b6e4ee67a5cd",
"class_f_e_a_t_1_1_tiny_1_1_matrix.html#a9f9957f9b9453b07ed920b9621365c2d",
"class_f_e_a_t_1_1_trafo_1_1_evaluator_base.html#a04385656ab1f3e606550f78430ae4b23",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_012_01_4_01_4.html#a766a02bc1af24aff8745b18eaf578007",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#acc3c782c86dc97d0c023ab1ad7e0f4e9",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree___00_01_shape_1_1_vertex_01_4.html#a97164a88685009b2dfcd657bc4e65679",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator.html#add9092926384fb5314de0140b987b0e9",
"class_f_e_a_t_1_1_voxel_assembly_1_1_voxel_burgers_velo_material_assembler_3_01_q2_standard_hype407418e5507f47d8c147d0c95b2976a9.html#a664d7da214324bf2f0a552c3d2e3e2e5",
"coding_conventions.html#codeconv_general_conventions",
"function_8hpp_source.html",
"likwid__marker_8hpp.html#a903341fdf602c35419e174726a19aa2f",
"namespace_f_e_a_t_1_1_assembly.html#a28186a335b7da11c0a4a2a3770462bf0a3c5097fb5c0edb92db17a1740f6ae4bf",
"namespace_f_e_a_t_1_1_pack.html#a422a077b3ced78b7b524168b89825943",
"namespace_f_e_a_t_1_1_space_1_1_argyris.html",
"preproc_macros.html#ppm_feat_compiler_intel",
"struct_f_e_a_t_1_1_assembly_1_1_f_e_interpolator_3_01_space_1_1_lagrange1_1_1_element_3_01_trafo88bd382a0f92c8baea41863ea8ff4cff.html#aa41ea98faf08f9869c459a3e791aff5b",
"struct_f_e_a_t_1_1_geometry_1_1_adaption_stats.html#aee9576bc34c1a08e954c1d123d7b06b7",
"struct_f_e_a_t_1_1_geometry_1_1_refinement_template_3_01_shape___00_01_shape_1_1_vertex_01_4.html#a99a9ba3eb6e2d8d07587b4b57462e9cb",
"struct_f_e_a_t_1_1_pack_1_1_type_helper_3_01_f_e_a_t_1_1_type_1_1_floating_class_00_01true_01_4.html",
"struct_f_e_a_t_1_1_space_1_1_lagrange2_1_1_dof_traits_3_01_shape___00_010_01_4.html",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___004d566ecc36dd8a41c0e23bcb70220e7b.html#a9500f4533e3e9c6c6a70208efce32894",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01signed_01short_01_4.html#a3b2f11d533d51e97c53ecc02365c57de",
"tools_meshtools.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';