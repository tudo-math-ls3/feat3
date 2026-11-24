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
"class_f_e_a_t_1_1_adjacency_1_1_composite_adjactor_1_1_image_iterator.html",
"class_f_e_a_t_1_1_analytic_1_1_auto_derive.html#a5b711c12c1465de1182d4e59acb12262",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_frankes_function.html#aec446a7b61f64af49501c7a033be5e84",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_poiseuille_pipe_flow.html#a4f8d060e44827f08b3403cbe8a69e679",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_sine_ring_vortex_pres2_d_1_1_evaluator.html#a8f688f2dc1dd4cf7291a819a751e4152",
"class_f_e_a_t_1_1_analytic_1_1_common_1_1_taylor_green_vortex_velo2_d_1_1_evaluator.html#a1b624e7d52058a5599660145723d4a45",
"class_f_e_a_t_1_1_analytic_1_1_distance_1_1_inverse_distance_function_1_1_evaluator.html",
"class_f_e_a_t_1_1_analytic_1_1_lambda_scalar_function3_d.html#a5759200287edd9018266ee196896bf79",
"class_f_e_a_t_1_1_analytic_1_1_parsed_vector_function.html#a2f395c33a4e56dfe1e3a490033642dd6",
"class_f_e_a_t_1_1_arg_parser.html#a6a70a019b035d166987bcafba555f9df",
"class_f_e_a_t_1_1_assembly_1_1_burgers_assembler_carreau.html#a2dfa25520ecd8bf0d0aea3a1f836b0b7",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_assembly_task_base.html#af80337cfa1a32793a469d2656f292426",
"class_f_e_a_t_1_1_assembly_1_1_burgers_blocked_matrix_assembly_job_1_1_task.html#a8858542347baf43b93ed4731d3047209",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_assembly_task_base.html#a120265fc700ee2f49b8f1f534739679f",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_matrix_assembly_job_1_1_task.html#a17939ba94669734c400c5442e57e8b1a",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_blocked_vector_assembly_job_1_1_task.html#af57e7535d9563c30127611f2b8b1857c",
"class_f_e_a_t_1_1_assembly_1_1_burgers_carreau_scalar_matrix_assembly_job_1_1_task.html#af2501e91e71a9575a6a3c783cb8c9f42",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_assembly_task_base.html#ad2052a0e5a80707c31351842cdb11821",
"class_f_e_a_t_1_1_assembly_1_1_burgers_scalar_vector_assembly_job_1_1_task.html#a56faffd7256b9d2b897a3ee00f631328",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_du_dv_operator_blocked_1_1_evaluator.html",
"class_f_e_a_t_1_1_assembly_1_1_common_1_1_laplace_operator_1_1_evaluator.html#a3445d9c7997262f2e2692c1c5aa769a8",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembler.html#a02a036ef4ee28397c32c9e6e2665018a",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_basic_matrix_task_c_r_t_p2.html#a2198c069fc9330d28cf68194288b034f",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_bilinear_operator_matrix_job2_1_1_task.html#a11e4760f3df5f79adfc59100578dd941",
"class_f_e_a_t_1_1_assembly_1_1_domain_assembly_linear_functional_vector_job.html#ab1990ed3d0a1515f13a07c53b2825b64",
"class_f_e_a_t_1_1_assembly_1_1_matrix_discrete_eval_data.html#aea3a60e60c5d3c162d106cec1e01b1c2",
"class_f_e_a_t_1_1_assembly_1_1_stokes_f_b_m_assembler.html#a192349b34b39806dc907e5ec992f692c",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_matrix_task_c_r_t_p1.html#a9c384a8d33215f36b13a291ed6f87c59",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_basic_task_base2.html#ac99df467f6e00bbe32224763726ce38c",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_bilinear_operator_matrix_job1_1_1_task.html#a4b6ad975b3ed3af60bc9c05dfcd5f4f3",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_discrete_function_integral_job_1_1_task.html#aa3d2a9ba34b6cb01567c801387d0f348",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_mass_matrix_job_1_1_task.html#ae6b12d0e8e77841c5d5ce496ba7d94b8",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_jump_task_base1.html#abe22ca72129ad05a0796a15f358ca0ea",
"class_f_e_a_t_1_1_assembly_1_1_trace_assembly_stokes_vector_analysis_task_c_r_t_p.html#a0ba9ad5a4a0efa91c7ca5cf4d6c9b8fa",
"class_f_e_a_t_1_1_control_1_1_blocked_basic_system_level.html#a4929d6fa47e6e4599fcd9c06f9f349c1",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_adaptive_parti_domain_control.html#a620d7fc91f55191c9e0d72c761598515",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_deslagged_domain_control.html#a844540e9168903501f5d266c2e35071d",
"class_f_e_a_t_1_1_control_1_1_domain_1_1_parti_domain_control.html#ac1c81fe91c0da812d53fa2188fc48ac6",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_du_dv_functional_control.html#a286a48a4eabd974accb892b1fb469239",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_meshopt_control_base.html#ac8c2e96cd77b1a85ee66f840c8aec7f7",
"class_f_e_a_t_1_1_control_1_1_meshopt_1_1_quadratic_system_level.html#a82e99f27e7f0c6dd64af9b753324df8f",
"class_f_e_a_t_1_1_control_1_1_scalar_unit_filter_system_level.html#ad81bdc8156d4b005ce5e4767598d6939",
"class_f_e_a_t_1_1_cubature_1_1_factory_explicit_wrapper_3_01_shape_1_1_simplex_3_013_01_4_01_4.html",
"class_f_e_a_t_1_1_cubature_1_1_shunn_ham_driver_3_01_shape_1_1_simplex_3_013_01_4_01_4.html#a80d329f6df20b8441cf3576ed21fbbde",
"class_f_e_a_t_1_1_dist_1_1_request.html#a18bb91893fb1e026e7daf502b4bcd078",
"class_f_e_a_t_1_1_file_error.html#ac3251651933f53a6f97b68c9348205c0",
"class_f_e_a_t_1_1_geometry_1_1_adaptive_mesh_layer.html#a9fe8c344fec9ff5ac06ab08c776fadab",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_c_g_a_l_surface_mesh.html#a7619901767a0842436e3b8c616b6df81",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_extrude.html#a3970c5cfb1da24ef41017d04afefd286",
"class_f_e_a_t_1_1_geometry_1_1_atlas_1_1_surface_mesh.html#ab5714d8e03017b8c10330ec21f64aeb6",
"class_f_e_a_t_1_1_geometry_1_1_conformal_mesh.html#ac2307143ce8b8b95e9b09ad6e6321e8f",
"class_f_e_a_t_1_1_geometry_1_1_factory_3_01_structured_mesh_3_01shape__dim___00_01num__coords___00_01_coord___01_4_01_4.html#a3efeb987a4b167b9a3e6fe422e1a8146",
"class_f_e_a_t_1_1_geometry_1_1_intern_1_1_adaptive_mesh_storage.html#a08252d9f556e24c68a9413df897e752c",
"class_f_e_a_t_1_1_geometry_1_1_mesh_distortion.html#ac2a6e8dc73cde941ce66035c33f37d5a",
"class_f_e_a_t_1_1_geometry_1_1_mesh_node_1_1_mesh_part_node_bin.html",
"class_f_e_a_t_1_1_geometry_1_1_mesh_part_operations.html#ae01a2eff04c94f0dacf30e929d3385b0",
"class_f_e_a_t_1_1_geometry_1_1_parti_zoltan.html#a74db2b8ee4a766fac4b45400da0d074e",
"class_f_e_a_t_1_1_geometry_1_1_patch_part_map_holder.html#af69b7e5d097548e21d72ca19921313e7",
"class_f_e_a_t_1_1_geometry_1_1_shape_convert_factory_3_01_mesh_part_3_01_conformal_mesh_3_01_shaef7d1b90fe5d8f8564ae577ae496f5bd.html#ae7a84a1a7c60ba74bddfaa762e99726f",
"class_f_e_a_t_1_1_geometry_1_1_structured_mesh.html#af0e0db5896fd215feb52ae36fe0f2a28",
"class_f_e_a_t_1_1_geometry_1_1_template_search_space_3_01_shape_1_1_vertex_00_01num__coords___01_4.html#a2aa172485412a541ef38b9c129a00f48",
"class_f_e_a_t_1_1_geometry_1_1_voxel_map.html#a76f77bec2918a480d8f19e323b0bf1fc",
"class_f_e_a_t_1_1_global_1_1_alg_dof_parti_system.html#a0030f06dbd18fed3915aa2bab4e4d5f9",
"class_f_e_a_t_1_1_global_1_1_global_inverse_mapping_data.html#a1f721e4050b2305f1de27dc3511a3170",
"class_f_e_a_t_1_1_global_1_1_nonlinear_functional.html",
"class_f_e_a_t_1_1_global_1_1_splitter.html#a81a9e5be48659ff08e5b1b770af7881a",
"class_f_e_a_t_1_1_global_1_1_vector.html#a227754d98d3bb90f267f511bf4ac7522",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_matrix.html#a484ec4cbdba1f480f18395bfc965648d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_dense_vector.html#ad486486693c22a3b456aad94630b2ada",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_filter_chain.html#a486e37ce9540f5d4d0d2defe45e05b4d",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_mean_filter.html#a7ff67497d511f19dbca41577bc6dd7ed",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_null_matrix.html#a5026b638e94169493a23ffcfb673d98f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_col_matrix.html#a8994516a57c288a94e88ce6a38dfd32f",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_full_matrix.html#a8a510e486b8ff9d7067264d852e34421",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_power_vector.html#a73b6472b2fd5f188fee60bc270b6b428",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_slip_filter.html#a00ef3ad57ba29ac341da441845adefb6",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_c_s_r.html#ac0e6a382d5819363a98d409c9c8b7053",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_b_wrapped_c_s_r.html#ab2814c85280cff88b2816a47bfa5adc0",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_banded_1_1_adjactor.html",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_matrix_c_s_r.html#a2aa8ee635d435a88f6af228ed48ee17e",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector.html#a0a94d17901d3eb78142e1d3d4993938b",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_sparse_vector_blocked.html#acbc7d933ca4ff8332198d6bbe9baf6fd",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_filter.html#a62da6456c47dcc1be8f339f9ee28e5d3",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_tuple_mirror.html#afac99d4fd61260624706878d48e77733",
"class_f_e_a_t_1_1_l_a_f_e_m_1_1_unit_filter_blocked.html#aabf2574f6026ff2cee668745fc6383ed",
"class_f_e_a_t_1_1_memory_1_1_arbiter.html#ac396469ef4110c2614f5aacb68e27f4f",
"class_f_e_a_t_1_1_mesh_opt_1_1_mesh_deformer.html#a1a921727049e09bcd06c75151b221f7a",
"class_f_e_a_t_1_1_meshopt_1_1_du_dv_functional.html#a87087a791b3e8c2e78f74b074d09bf55",
"class_f_e_a_t_1_1_meshopt_1_1_mesh_concentration_function.html#a0790d3e7264d3dce8906ed301937c3e4",
"class_f_e_a_t_1_1_meshopt_1_1_rumpf_functional_base.html#a5e6c5589ff118414fe38ba4f3b81f74a",
"class_f_e_a_t_1_1_simple_arg_parser.html#ab36fc3ab39bbf957c5c75fdfe9d9b0ac",
"class_f_e_a_t_1_1_solver_1_1_a_d_p_solver_base_3_01_global_1_1_p_m_d_c_d_s_c_matrix_3_01_matrix_040b7e1caa6985ce9572fb7c7833a8b9.html#abbe1392a88c229e6cd29cc1ab535067e",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_c_g.html#ae125ac188382e3e587aa9f1e13b08432",
"class_f_e_a_t_1_1_solver_1_1_a_l_g_l_i_b_min_l_b_f_g_s.html#add5c03ded9d5568e14ec60cdfb1572b8",
"class_f_e_a_t_1_1_solver_1_1_ama_vanka_g_m_r_e_s_3_01_global_1_1_matrix_3_01_local_matrix___00_0d762f8676132651836c5e02fb8066064.html#ab48c468782b6695d34bee1089a88b785",
"class_f_e_a_t_1_1_solver_1_1_analytic_function_operator.html#a0a4f2a2fcfef156c1a49ea01baa35b6f",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab.html#a7647edeb0b1945e4a74e8d46e95a9c4c",
"class_f_e_a_t_1_1_solver_1_1_bi_c_g_stab_l.html#aad1386af21c34fca491aca43cab04c58",
"class_f_e_a_t_1_1_solver_1_1_c_u_d_s_s.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_convert_precond_3_01_global_1_1_vector_3_01_local_vector_outer___00120b74f43424364c71ee3ecd9deab297.html#a34c36cfa7dbdf201a395c5711a4637f6",
"class_f_e_a_t_1_1_solver_1_1_direct_sparse_solver.html#a49bb14b2fdf83babbaf47fc72b2835d9",
"class_f_e_a_t_1_1_solver_1_1_direct_stokes_core_3_01_solver_d_t___00_01_solver_i_t___00_01_l_a_fba3df488bc06c314d472100c65203e90.html#abde855311e07004d13ba47b40985468d",
"class_f_e_a_t_1_1_solver_1_1_expression_base.html",
"class_f_e_a_t_1_1_solver_1_1_f_r_o_sch_preconditioner.html#a5cec7693df4ef68d8634d1073ace568d",
"class_f_e_a_t_1_1_solver_1_1_fixed_step_linesearch.html#aa9f9b03fe5fd3bc40f110c6321824ece",
"class_f_e_a_t_1_1_solver_1_1_gropp_p_c_g.html#ad29a0a21a3e91eaf8b46e12da46d0468",
"class_f_e_a_t_1_1_solver_1_1_i_d_r_s.html#a69fab2263b0c35bf03aa77d95e3a5917",
"class_f_e_a_t_1_1_solver_1_1_iterative_solver.html#a27d93b800c82dae6cf39f5ffc9a94015",
"class_f_e_a_t_1_1_solver_1_1_linesearch.html#a3216d57242be2e62ed958a84be4ff2e5",
"class_f_e_a_t_1_1_solver_1_1_m_q_c_linesearch.html#a1f17cb192fab20b3591eafdafd7c5f6e",
"class_f_e_a_t_1_1_solver_1_1_matrix_precond.html#ae162070c491d897ae6f91d263200af40",
"class_f_e_a_t_1_1_solver_1_1_multi_grid_level_std.html#a542ee96260b79c3ccc1df1b68b1c6b26",
"class_f_e_a_t_1_1_solver_1_1_n_l_c_g.html#ad5a06be09b1b40ab2c5b7cceb8274ab8",
"class_f_e_a_t_1_1_solver_1_1_n_l_opt_l_s.html#aec87d011de4859420af0f06559fbef2e",
"class_f_e_a_t_1_1_solver_1_1_n_l_s_d.html#adb4aee07c7b4c715bdf65b50a8beb4cc",
"class_f_e_a_t_1_1_solver_1_1_newton_raphson_linesearch.html#ae6fc9035510853c34e062754732c17ba",
"class_f_e_a_t_1_1_solver_1_1_p_c_g.html#af8e655f0cdc31c3b245587b58a14d5b3",
"class_f_e_a_t_1_1_solver_1_1_p_c_g_n_r_i_l_u.html#a6e0074b3a336238364169fb1a0fb0701",
"class_f_e_a_t_1_1_solver_1_1_p_c_r.html#ad5a06be09b1b40ab2c5b7cceb8274ab8",
"class_f_e_a_t_1_1_solver_1_1_p_s_d.html#a4fdec11243202ce55a49e67d259cab2a",
"class_f_e_a_t_1_1_solver_1_1_pipe_p_c_g.html#a0810de389dc450bea5cb6917b66ca051",
"class_f_e_a_t_1_1_solver_1_1_preconditioned_iterative_solver.html",
"class_f_e_a_t_1_1_solver_1_1_q_penalty.html#a6db06343bf3a13f014148a1e831c87cd",
"class_f_e_a_t_1_1_solver_1_1_r_bi_c_g_stab.html#abfeec9e7a424dee82dccb22be9f43fd4",
"class_f_e_a_t_1_1_solver_1_1_richardson.html#a25362c6f4bb0c5d1deac80d947442303",
"class_f_e_a_t_1_1_solver_1_1_s_s_o_r_precond_with_backend_3_01_preferred_backend_1_1generic_00_0586faebaaaa6acadbef21e802e831fb5.html#ac3efee3fdc3ac5ff4132ab7743dd7938",
"class_f_e_a_t_1_1_solver_1_1_secant_linesearch.html#a7647edeb0b1945e4a74e8d46e95a9c4c",
"class_f_e_a_t_1_1_solver_1_1_stokes_f_r_o_sch_preconditioner.html#aae0dbe1f8aa005dd40e12d0cb37f4cd6",
"class_f_e_a_t_1_1_solver_1_1_tpetra_solver_base.html#ae5ae58d38f8a98ab81805a03b51244a2",
"class_f_e_a_t_1_1_solver_1_1_vanka_factor_error.html#a76bcfc0b1c2ed8cc4d368d67c2c5c193",
"class_f_e_a_t_1_1_space_1_1_basis_data.html#aa394707a1671d9afe932bdf44115c76a",
"class_f_e_a_t_1_1_space_1_1_bernstein2_1_1_node_functional.html#aba5b55179b16f3517fee72d4ce4e0540",
"class_f_e_a_t_1_1_space_1_1_cai_dou_san_she_ye_1_1_element.html#af2b535319e8c52fa70246e6406d14320",
"class_f_e_a_t_1_1_space_1_1_cro_rav_ran_tur_1_1_evaluator_3_01_space___00_01_trafo_evaluator___03bb1a3d92225206453ee71c00db5aa7a.html#a15914083d87d19f34c58c4e6f9d5331e",
"class_f_e_a_t_1_1_space_1_1_discontinuous_1_1_element.html#af0ecd85c7feb185ca885caf9eaa48072",
"class_f_e_a_t_1_1_space_1_1_dof_assignment_null.html#a8b994c0cb6ba94e5475b21cea15768a0",
"class_f_e_a_t_1_1_space_1_1_eval_data_reduced.html#a850a2fb90cd4ed3368131ab9f7634f06",
"class_f_e_a_t_1_1_space_1_1_hermite3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_spd08b4bf16ff97c315f2f45da544757eb.html#a4b357a4785d45ec87feda720a8f8157a",
"class_f_e_a_t_1_1_space_1_1_lagrange1_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s76531dae8562723189ff92dab0342885.html#ae4d17e132c111387fbe169819fc1c2aa",
"class_f_e_a_t_1_1_space_1_1_lagrange2_1_1_evaluator.html#a24e648bd29f80680fe2028b3fe286673",
"class_f_e_a_t_1_1_space_1_1_lagrange3_1_1_evaluator_3_01_space___00_01_trafo_evaluator___00_01_s4e94f7d8b33a4b573a80bb38af238451.html#ae8fcaa4c9aae27813c462673db99a82a",
"class_f_e_a_t_1_1_space_1_1_p2_bubble_1_1_element.html",
"class_f_e_a_t_1_1_space_1_1_q1_t_b_n_p_1_1_evaluator.html",
"class_f_e_a_t_1_1_statistics.html#a1d9d72a66b62abcde8693acc6431562c",
"class_f_e_a_t_1_1_test_system_1_1_test_list.html#a329a2f833328ba0c2c884319c0bad255",
"class_f_e_a_t_1_1_tiny_1_1_tensor3.html#a857ae1645a90f95789a8013ccc2ab3cf",
"class_f_e_a_t_1_1_trafo_1_1_hierarchical_inverse_mapping.html#a179e492138d7615d19c4f5ad460db0c9",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_hypercube_3_012_01_4_01_4.html#a99aacf4b77b94a83e0558da0601a233f",
"class_f_e_a_t_1_1_trafo_1_1_iso_sphere_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#aeda24d9321bc45990ffed6598de74841",
"class_f_e_a_t_1_1_trafo_1_1_isoparam_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01degree___00_01_shape_1_1_vertex_01_4.html#ab5feec13b34277832ba7eb410667b1e4",
"class_f_e_a_t_1_1_trafo_1_1_standard_1_1_evaluator_3_01_trafo___00_01_eval_policy___00_01_shape_1_1_vertex_01_4.html#a04385656ab1f3e606550f78430ae4b23",
"class_f_e_a_t_1_1_voxel_assembly_1_1_voxel_burgers_velo_material_assembler_3_01_q2_standard_hype407418e5507f47d8c147d0c95b2976a9.html#a5682ead85a86fcb63e1ac2af72c5390e",
"coding_conventions.html#codeconv_nameconv_camelcase",
"feat_tpls.html#including_tpls_finding",
"likwid__marker_8hpp.html#a849dcba20310df11544762fd16d95d24",
"namespace_f_e_a_t.html#ad574925e947c823619c356d5d12b0815a2063c1608d6e0baf80249c42e2be5804",
"namespace_f_e_a_t_1_1_l_a_f_e_m.html#ada0e93baecfc2c4933d9c5eb9c96a646a48e8d7c8a7612515950a8bcbf57562a1",
"namespace_f_e_a_t_1_1_solver.html#a764bfc8e9007d466dea48527b04a8b4d",
"namespacemembers_func_f.html",
"scale__col__csr__generic-eickt_8cpp_source.html",
"struct_f_e_a_t_1_1_assembly_1_1_vector_error_info.html#a920b55120cc114bfc135eab7f22918e9",
"struct_f_e_a_t_1_1_geometry_1_1_entity_reference.html#a6f2b844847d73af17dc51c8e1248dc7d",
"struct_f_e_a_t_1_1_geometry_1_1_vertex_set.html",
"struct_f_e_a_t_1_1_pack_1_1_swap_helper.html",
"struct_f_e_a_t_1_1_space_1_1_hermite3_1_1_dof_traits_3_01_shape_1_1_simplex_3_01shape__dim___01_4_00_010_01_4.html#a8f624dbbb468623f1d4e18fffadef663",
"struct_f_e_a_t_1_1_trafo_1_1_standard_1_1_eval_helper_3_01_data_type___00_01_dom_point_type___002b260e771ecb4737ffdd677d6d46ebb5.html#a8582fd110567efa1a8f94c4434df3e39",
"struct_f_e_a_t_1_1_type_1_1_traits_3_01signed_01int_01_4.html#a0438c5f06fce51dcae6a0da21d717cc1",
"struct_f_e_a_t_1_1_voxel_assembly_1_1_space_helper_3_01_q2_standard_f_e_3_01_shape___01_4_00_01_d_t___00_01_i_t___01_4.html#ad7cba318aad480ed8edde03507e94655"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';