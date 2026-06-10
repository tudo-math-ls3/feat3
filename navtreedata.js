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
      [ "Important Concepts used in FEAT", "index.html#main_sec_concepts", null ],
      [ "Mesh related Topics", "index.html#main_sec_mesh", null ],
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
    [ "LIKWID marker API with FEAT", "likwid_for_feat.html", [
      [ "LIKWID General Information", "likwid_for_feat.html#LIKWID_descr", null ],
      [ "Prerequisites for using LIWWID", "likwid_for_feat.html#LIKWID_prereques", null ],
      [ "Using the LIKWID marker API", "likwid_for_feat.html#LIKWID_use", null ]
    ] ],
    [ "The Memory Arbiter and its View classes", "memory_arbiter.html", [
      [ "Motivation", "memory_arbiter.html#memarb_sec_motivation", null ],
      [ "Overview", "memory_arbiter.html#memarb_sec_overview", null ],
      [ "The Memory Arbiter class", "memory_arbiter.html#memarb_sec_arbiter", [
        [ "Allocation and Initialization", "memory_arbiter.html#memarb_arbiter_alloc_init", null ],
        [ "Shallow and Deep Copies of Memory Arbiters", "memory_arbiter.html#memarb_arbiter_copies", null ],
        [ "Freeing Memory and Releasing Arbiters", "memory_arbiter.html#sec_memarb_arbiter_release", null ]
      ] ],
      [ "Memory Views, Locations and Access Rights", "memory_arbiter.html#memarb_views", [
        [ "Moving Views", "memory_arbiter.html#sec_memarb_view_move", null ],
        [ "Memory Locations", "memory_arbiter.html#sec_memarb_location", null ],
        [ "Memory Access Rights", "memory_arbiter.html#sec_memarb_access", null ]
      ] ],
      [ "Overlapping Memory Views", "memory_arbiter.html#sec_memarb_overlap", null ]
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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html#a5deb85b80b3dcd4565f468aa1a26fda9",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive_1_1_evaluator.html",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_frankes_function_1_1_evaluator.html#a1d733a386a49b31a267c71e0f98aa54d",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow.html#ad61f5cf78f6909cf082bd4c965b6e2cd",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_r_h_s2_d.html#a0616b2a6a58f7022f418ccf875d093f9",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_taylor_green_vortex_velo2_d_1_1_evaluator.html#ad7bd42502d2c18422db75e7998ab81a0",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_inverse_distance_function_1_1_evaluator.html#aca097573aea955a0ade7ef2d6088b577",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#ac79af4af8a677ec0176c6f4c24fc0682",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a76d73a4285209c7aef3365a12245da27",
"class_f_e_a_t_1_1_assembly_1_1_asm_traits1.html#a1ad4400bfc1d1cf9b3bb58e4049e36be",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembly_job_base.html#a6abf0a2d3fa7eeb7f230a8c89afc6a5d",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_diagonal_assembly_job.html#ae8f678f9ee0412453b48613066438644",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#af80337cfa1a32793a469d2656f292426",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_assembly_task_base.html#a9040fa67d8287413afc76c16ec5d1769",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_matrix_assembly_job_1_1_task.html#ab58683615104315bb0237c0561cecfc8",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_scalar_assembly_task_base.html#a41fe6ba2bab71139385c23ea9f051443",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_scalar_vector_assembly_job.html#a80d39ec9e681f0266744fe7d71904e42",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_matrix_assembly_job.html#a85a8c84b44a6d4685f865359fc2784c7",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#ad1ffa1d8dd3067c503578fdace46faed",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_force_functional_1_1_evaluator.html#a5dac3ee8e624548ef5d6092dc218359c",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_operator_blocked_1_1_evaluator.html#aafc3d06f8a01014085c35e2b3d3abedc",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a64adf5a8939e644040d355c251ac136f",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_basic_vector_task_c_r_t_p.html#a001fe27a4126b430ed66916ffc136c29",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_discrete_function_integral_job.html#a4ee7010c169e897ac9f3f41d0f7fe473",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_linear_functional_vector_job_1_1_task.html#ab8a7ea9d8a47fde64fb093e2ca29e8f5",
"class_f_e_a_t_1_1_assembly_1_1_normal_value_surface_integrator_job.html#a691a02fce8dcbd3aa532a1488c0665c1",
"class_f_e_a_t_1_1_assembly_1_1_surface_integrator.html#ac5c210dd5c4395cb82a3edf02203f7d4",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_matrix_task_c_r_t_p2.html#a887aed2bb8e6a28edd5fd1de5ef90d08",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_vector_task_c_r_t_p.html#abf14499736c69c330f710cb897c2284a",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_bilinear_operator_matrix_job2.html#aeee0d80710dcacdfb222677018288359",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_error_function_integral_job_1_1_task.html#a7f96794c68bd4d2965ebb020fd1e26f8",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_matrix_task_c_r_t_p1.html#ac5f8b4d8a394d5a281bac6d54e7077ca",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_linear_functional_vector_job_1_1_task.html#a86d14678b242e74bd5664879e245cc3a",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_stokes_vector_analysis_task_c_r_t_p.html#aea10a6b6ce76d147449e98f96a397ca1",
"class_f_e_a_t_1_1_control_1_1_blocked_unit_filter_system_level.html#a92ccf4288666a8bcda1a7a070f5e435b",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_adaptive_parti_domain_control.html#afa6a43e0bfb0b52b14ba9a19c17e94b3",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#a02d0cb614bcf85e53801eb43f16ec5b8",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control_base_1_1_ancestor.html#a2bc67b7da85af1e459ce527cff2287bf",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_voxel_domain_control.html#afddf5df0ac7908f4b270aa208ed76012",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_control_base.html#a68f9baf37018044e36b3f58560de9285",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_quadratic_system_level.html#a31d10a74e832b9fcac82c94b9586de4f",
"class_f_e_a_t_1_1_control_1_1_scalar_unit_filter_system_level.html#a26642cf398bbc5dbf74bda2e229d77c7",
"class_f_e_a_t_1_1_cubature_1_1_dunavant_driver_3_01_shape_1_1_simplex_3_012_01_4_01_4.html#a2eacfd10a6d98123f3d34c6684f424f7",
"class_f_e_a_t_1_1_cubature_1_1_scalar_1_1_rule.html#ab2b98e23aeb08f37b88bc4d181812027",
"class_f_e_a_t_1_1_dist_1_1_comm.html#add827514d1847da75268e4a851088c7b",
"class_f_e_a_t_1_1_dist_file_i_o.html#adbac17ec08b4fba595d2fadf12c47475",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_mesh_layer.html#a2f8afb482b60cab59f413fe79a6277e8",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_bezier_points_parser.html#a01ff5bb4414e2d4712064bf34bc0a8cc",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_circle_chart_parser.html#a593826ae6226710e09b6600baf2fb9ba",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_surface_mesh.html#a8846bbb31303df6a2cc109c8a89bb772",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#a96116e8891cf94396cab9e195c68404a",
"class_f_e_a_t_1_1_geometry_1_1_factory_3_01_mesh_part_3_01_mesh_type___01_4_01_4.html#a4573c9c8ca882515a025920e1dcf2d83",
"class_f_e_a_t_1_1_geometry_1_1_index_tree.html#a3b911ff616d3e76a4a7724208349f38a",
"class_f_e_a_t_1_1_geometry_1_1_mesh_atlas.html#ac4bff8d221d43f79a85a6c32bcd9a3c3",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node.html#ab21a3c549e2d2927a6f21c28621dabe6",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part_operations.html#a3d93a0eee9791bd0eaf68998890e59d5",
"class_f_e_a_t_1_1_geometry_1_1_parti_zoltan.html#a6854ca4ba927d1ba4e7425c6c880995d",
"class_f_e_a_t_1_1_geometry_1_1_patch_part_map_holder.html#af5228a718a6d2ac353db32dbbdd1bd0c",
"class_f_e_a_t_1_1_geometry_1_1_shape_convert_factory_3_01_mesh_part_3_01_conformal_mesh_3_01_shaef7d1b90fe5d8f8564ae577ae496f5bd.html#ad546af94da5d695f3337f4c81f4e2c5a",
"class_f_e_a_t_1_1_geometry_1_1_structured_mesh.html#ac92f314446aa141a8ccde12c058c288a",
"class_f_e_a_t_1_1_geometry_1_1_template_search_space_3_01_shape_1_1_vertex_00_01num__coords___01_4.html#a03a7b2cfe1b6a25ea1a79cdaf2b27a0e",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#a716fa19031b2457f0f5f81d052edcc94",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_system.html",
"class_f_e_a_t_1_1_global_1_1_inverse_mapping.html#a38999b8d8f12c25bb2a81a38a3c6d221",
"class_f_e_a_t_1_1_global_1_1_nonlinear_functional.html#a58438bc84f28a56bd333257a249b2a1a",
"class_f_e_a_t_1_1_global_1_1_symmetric_lumped_schur_matrix.html#a082c691f5fe3a4b3ec35bcb9b821d0d3",
"class_f_e_a_t_1_1_global_1_1_vector.html#a72d408e0b585dded66f234aa31b37684",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a834618d5a0078170e8eef67560c6ae45",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector_blocked.html#a0004af8000e4b65540748bd9a848d9ac",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_sequence.html#a4c9feb0d0e77b7afef0cd0e9b0ea88ac",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter_blocked.html#a9ef65ec93c9f83393439d32c099f5c78",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#ae32e9f4468717510196ff3f29a2d179a",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_diag_matrix.html#a2353664312d58b1227c83f3b11173973",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_mirror.html#a814486efd70a770aa7b2dbdf888f1283",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_saddle_point_matrix.html#a0bfc7d866e1de15d94bdfb3fa7d19271",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#a5aa4438c3792a45d0e485b579bb92618",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#a319f4512ada6844ad97fdd8dd89d07e7",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded.html#a7d9ecfa50a9478d839ffb897e6584796",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_c_r.html#abfdb3ce63ec1f76d27d4c1e7dbc12e7b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#ad338e98c1fd0978b9bfb5c4da26b5bbf",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#a2f2853fe053807e2d8250833032032a1",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_diag_matrix.html#a6a6c00e077edf435e3a172fb7f220251",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_matrix_row.html#ab1c12b1d919d0868ebc0b796be3df299",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter.html#abd4eb7bf12c373e76e1bc7f4a01116d8",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_vector_mirror.html#af5dc00e2137f8c37f5462b3423f7bd0d",
"class_f_e_a_t_1_1_memory_1_1_typed_view.html#ae6c9cdfe5f5b3f0f77b4d51e42198981",
"class_f_e_a_t_1_1_meshopt_1_1_du_dv_functional.html#a43dc6f51b94988bb999bc897a684df9a",
"class_f_e_a_t_1_1_meshopt_1_1_hyperelasticity_functional.html#ac794b4353251af905950edd42a4bc511",
"class_f_e_a_t_1_1_meshopt_1_1_rumpf_functional.html#a6decfabe26926dc3500fedbc77d16469",
"class_f_e_a_t_1_1_simple_arg_parser.html#a0975076dc12b65829226850bef1558a0",
"class_f_e_a_t_1_1_solver_1_1_a_d_p_solver_base_3_01_global_1_1_p_m_d_c_d_s_c_matrix_3_01_matrix_040b7e1caa6985ce9572fb7c7833a8b9.html#a4da3d5390745d8a5a1746d831e58ef70",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#acacbcaccbf7ec307778128928a901d59",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#ad29a0a21a3e91eaf8b46e12da46d0468",
"class_f_e_a_t_1_1_solver_1_1_b_f_b_t.html#a8b79414530c321b142fb1e6120dfc92d",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#a0810de389dc450bea5cb6917b66ca051",
"class_f_e_a_t_1_1_solver_1_1_boomer_a_m_g.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_chebyshev.html#ad355136ba5cd50b1f2ece6658ab1c0c8",
"class_f_e_a_t_1_1_solver_1_1_descent.html#aec76c75b2d2582ef2cf0ed928b583314",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#a3c2575b8dd54eeea43a0803c9e80be4f",
"class_f_e_a_t_1_1_solver_1_1_euclid_precond.html#a705fa6fa252ff552732351ee4f363c49",
"class_f_e_a_t_1_1_solver_1_1_f_g_m_r_e_s.html#ad0e0aa027e766929e7e741f6384f887a",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_g_m_r_e_s.html#a7d583ff48c7678f1cffa30fdd48315a7",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#abfeec9e7a424dee82dccb22be9f43fd4",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#a576f7cf3e1c847ac5611e477bd551b73",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#a44a696c771d3d261646f8de9a353cf38",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a576f7cf3e1c847ac5611e477bd551b73",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a41a9a9f80ce9e2355297cb3808d8e47f",
"class_f_e_a_t_1_1_solver_1_1_multi_grid.html#a50ff82108fdde3e797538ef0b7bc3d3a",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#a0f5f95eeb377ef784ac3224684fef393",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#afb3f9c6e873b3118aca41b5b3330552c",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_precond.html#af2ea6e8c96de45c8408f3dcc59b8f186",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#a0d4201e4324ec13233def993318e44f2",
"class_f_e_a_t_1_1_solver_1_1_nonlinear_operator_precond_wrapper.html#a89087e74542e4f0afdc7ac67ff132d15",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r.html#a3bd60beb7925eaeada45028627f85281",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#aa86287c6d142b78ac494420116688d66",
"class_f_e_a_t_1_1_solver_1_1_p_m_r.html#a0fd9bce681b27352bb6a37c2ef10e300",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a7e61ce0b51e7f97317042a184aba8613",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#a575442c5edf8721e18db2a49717f45a5",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html#a5079e20c7708d1409a0798b7a226e9de",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#aa7cf18c4f6e62acfd0c3cf3285c00b96",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#aeebdc6fb0b2e5072cb5c0c52205ff58b",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#a6063a9d31500d9ea5f299a97fa9f1e20",
"class_f_e_a_t_1_1_solver_1_1_saddle_umfpack_mean.html#ae0852a882b528b3c014402cb55ece2f5",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#aa95d79033c0912710a988a61e7312a15",
"class_f_e_a_t_1_1_solver_1_1_super_l_u.html#a07f5437d9e62afe0be9ead7bac2056ff",
"class_f_e_a_t_1_1_solver_1_1_umfpack.html#ad7e65dcf3170ce32b31d3302cde79317",
"class_f_e_a_t_1_1_solver_1_1_voxel_ama_vanka.html#a609024eacae0ad14e32f8086a5ea944c",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_element.html#abd80f2f43e19246bddf547654f60d852",
"class_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_element.html#a61ae0c9f6a924e658159f356d8378648",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_evaluator_3_01_space___00_01_trafo_evaluator_cd584e2ae53d2b4306ff680d7265a0ec.html#a3f2012e306e5fcd49316db981dcfeb13",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___03bb1a3d92225206453ee71c00db5aa7a.html#afb5ea011a2dd1d12c54dc0449111aa3f",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_964e26ad63d9177b3f747659af96be6e.html#ab70b19b6962179079f8a2d38ddddcd28",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_uniform.html#af5347b275c4767ab33631571fa08d2c5",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_element.html#a05c121f2092dbbfc574b8da1f745427b",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_node_functional.html#ac36f33ce44f2589fa5e008245c42dcb4",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_sbede9c86e9ef3e51ccbac5494ff49399.html#ae76abc18dc4af648989406abb35c9ae4",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_evaluator.html#afd49268c1770e25b0993037897e8acf2",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s5e9945e9db3d72bec06849dc7e957367.html#abeda000eff0270b2a43789f0148fb536",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_element_1_1_evaluator.html#a0d9853631c58ffb34ab456846cbf7987",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_2be5aeec58e5b7d19e38f7fb4e901a36.html#a8413b17c736e23b2476830b818dc9015",
"class_f_e_a_t_1_1_statistics.html#a8c5f522bb492231bc941e6237f9422d2",
"class_f_e_a_t_1_1_thread_fence.html",
"class_f_e_a_t_1_1_tiny_1_1_vector.html#a1eb9abbe04016410fe385690dd4c2368",
"class_f_e_a_t_1_1_trafo_1_1_inverse_mapping.html#a1ba9f7efd28af91b06f2f834e97daf2e",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_simplex_3_011_01_4_01_4.html#a3ab07eb0691a1b01baefab55b0a2d2b4",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_5b7e974ab9f9514411b05146b815324a.html#a25236e7b3dcb0d1a384f10aa46daaa9d",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree_fe815dfbe8417dc8e12ca8204196c316.html#a6c2c258b890f833c860b46bca281f3e7",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a9efefabb53d2efc80bd5c307c2daa4dd",
"class_f_e_a_t_1_1_xml_1_1_content_error.html#ac222f8bd0b0a8e6b9804ae43809b18d1",
"cuthill__mckee_8cpp_source.html",
"functions_q.html",
"matvecmult__dense__dense__mkl_8cpp_source.html",
"namespace_f_e_a_t_1_1_assembly.html#a105b50bbb2cca05a2fadb234d099560a",
"namespace_f_e_a_t_1_1_memory.html#a567578cf56f030248943beebab762ec6aefb2a684e4afb7d55e6147fbe5a332ee",
"namespace_f_e_a_t_1_1_solver.html#aa3e1a0d5d6876f6d4fe16f71003b257cadbb68c86754ccad2ee2908a92f34d3e4",
"parti__iterative_8hpp_source.html",
"standard__hexa_8cpp_source.html",
"struct_f_e_a_t_1_1_control_1_1_stokes_blocked_unit_velo_mean_pres_system_level.html#aeb334861993341bf2659c7f0db828291",
"struct_f_e_a_t_1_1_geometry_1_1_intern_1_1_element_key.html",
"struct_f_e_a_t_1_1_l_a_f_e_m_1_1_arch_1_1_mirror_dense_scatter.html",
"struct_f_e_a_t_1_1_space_1_1_bogner_fox_schmit_1_1_dof_traits.html",
"struct_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_0151719a61b7fd9347c54a261080665c24.html#a8056d2f66744f69f0069905890ade95b",
"struct_f_e_a_t_1_1_type_1_1_helper.html#a589aebd32dac2262fd82b43a61408f6c",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01unsigned_01short_01_4.html#a27c030bdf27eca375bc4b6a60bcd98f6",
"unit__filter__block__weak__bcsr__generic-eickt_8cpp_source.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';