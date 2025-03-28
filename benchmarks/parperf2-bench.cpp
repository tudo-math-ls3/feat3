// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// ====================================================================================================================
// Linear Algebra Backend Parallel Performance Benchmark #2
// --------------------------------------------------------------------------------------------------------------------
// This benchmark measures the parallel performance of different linear algebra backends by performing a modified power
// iteration loop on the 2D FEM 9-point stencil matrix, which converges towards the eigenvector corresponding to the
// largest eigenvalue of the matrix.
// This benchmark application supports parallelization via MPI and OpenMP as well as through the MKL and CUDA backends
// of the LAFEM containers. This benchmark constructs the 9-point stencil matrix as well as the halo mirrors by hand,
// so no actual mesh/trafo/space/assembly functionality is used by this benchmark application to keep it simple.
//
// The (modified) power iteration loop
//
//   x_{k+1} := (1/|x_k|_2) * A * x_k
//
// is split up into five sub-steps, which are measured separately:
// 1. compute euclidean norm of x_k locally on each process
// 2. perform an MPI allreduce to compute the global euclidean norm of x_k
// 3. scale x_k by the inverse of its global norm and save the result in vector t_k
// 4. compute local matrix-vector product x_k := A*_k
// 5. synchronize vector x_k by performing a halo exchange
//
// This benchmark app measures the runtime of the five sub-steps separately and it computes both the memory bandwidth
// as well as the floating point operations per seconds of the three local sub-steps. Furthermore, this benchmark
// also includes LIKWID markers for the three sub-steps that enable profiling on a hardware-counter basis.
//
// --------------------------------------------------------------------------------------------------------------------
//
// USAGE: parperf-bench-1 <mx> <my> <nx> <ny> <min_time> [<backend>] [<dt-it>]
//
// Parameters:
// <mx> <my>    Dimensions of MPI process grid.
// <nx> <ny>    Dimensions of global mesh; must be multiples of process grid dimensions mx and my.
// <min_time>   Minimum wall-clock runtime of power iteration in seconds.
// <backend>    The desired linear algebra backend; must be one of:
//              'none'   : use raw OpenMP-parallelized loops
//              'generic': use generic (OpenMP-parallelized) backend.
//              'mkl':     use Intel MKL backend.
//              'cuda':    use CUDA backend.
// <dt-it>      The desired data/index-type pair; must be one of:
//              'hi':      half   + 32-bit unsigned int (requires CUDA)
//              'hl':      half   + 64-bit unsigned int (requires CUDA)
//              'fi':      float  + 32-bit unsigned int
//              'fl':      float  + 64-bit unsigned int
//              'di':      double + 32-bit unsigned int
//              'dl':      double + 64-bit unsigned int
//              'qi':      quad   + 32-bit unsigned int (requires quadmath)
//              'ql':      quad   + 64-bit unsigned int (requires quadmath)
//
// \author Peter Zajac
//

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/kahan_accumulator.hpp>
#include <kernel/util/likwid_marker.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/domain_assembler_helpers.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>

#ifdef FEAT_HAVE_OMP
#include <omp.h>
#endif

using namespace FEAT;

String stringify_flops(double flops, int precision = 3, int width = 7)
{
  if(flops >= 1e+18)
    return stringify_fp_fix(1e-18*flops, precision, width) + " EFlop";
  if(flops >= 1e+15)
    return stringify_fp_fix(1e-15*flops, precision, width) + " PFlop";
  if(flops >= 1e+12)
    return stringify_fp_fix(1e-12*flops, precision, width) + " TFlop";
  if(flops >= 1e+9)
    return stringify_fp_fix(1e-9*flops, precision, width) + " GFlop";
  if(flops >= 1e+6)
    return stringify_fp_fix(1e-6*flops, precision, width) + " MFlop";
  if(flops >= 1e+3)
    return stringify_fp_fix(1e-3*flops, precision, width) + " KFlop";
  return stringify_fp_fix(flops, precision, width) + " Flop";
}

Adjacency::Graph create_elems_at_rank(const Index mx, const Index my, const Index nx, Index ny)
{
  const Index m = mx*my;
  const Index n = nx*ny;
  const Index lx = nx / mx;
  const Index ly = ny / my;
  Adjacency::Graph graph(m, n, n);

  Index* dom_ptr = graph.get_domain_ptr();
  Index* img_idx = graph.get_image_idx();

  dom_ptr[0] = 0;

  // loop over process grid
  for(Index iy = 0; iy < my; ++iy)
  {
    for(Index ix = 0; ix < mx; ++ix)
    {
      Index p = dom_ptr[iy*mx + ix];// = (iy*mx + ix) * lx * ly;
      XASSERT(p == (iy*mx + ix) * lx * ly);
      for(Index jy = 0; jy < ly; ++jy)
      {
        for(Index jx = 0; jx < lx; ++jx, ++p)
        {
          img_idx[p] = (iy*ly + jy) * nx + ix*lx + jx;
        }
      }
      dom_ptr[iy*mx + ix + 1] = p;
    }
  }
  XASSERT(dom_ptr[mx*my] == nx*ny);

  return graph;
}

// forward declaration
template<typename DT_, typename IT_>
int run(const Dist::Comm& comm, const Index mx, const Index my, const Index nx, const Index ny, const Index level, double min_runtime, PreferredBackend preferred_backend);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  Dist::Comm comm(Dist::Comm::world());

  if(argc < 7)
  {
    comm.print("\nUSAGE: parperf-bench-2 <mx> <my> <nx> <ny> <min_time> [<backend>] [<dt-it>]\n");
    comm.print("This benchmark measures the parallel performance of the FEAT linear algebra");
    comm.print("backends by performing a power iteration to estimate the largest eigenvalue");
    comm.print("of the 2D 9-point stencil matrix stored in the standard CSR matrix format.");
    comm.print("\nCommand Line Parameters:");
    comm.print("<mx> <my>    Dimensions of MPI process grid.");
    comm.print("<nx> <ny>    Dimensions of global mesh; must be multiples of process grid dimensions.");
    comm.print("<min_time>   Minimum wall-clock runtime of power iteration in seconds.");
    comm.print("<backend>    The desired linear algebra backend; must be one of:");
    comm.print("             'generic': use generic (OpenMP-parallelized) backend.");
    comm.print("             'mkl':     use Intel MKL backend.");
    comm.print("             'cuda':    use CUDA backend.");
    comm.print("<dt-it>      The desired data/index-type pair; must be one of:");
    comm.print("             'hi':      half   + 32-bit unsigned int (requires CUDA)");
    comm.print("             'hl':      half   + 64-bit unsigned int (requires CUDA)");
    comm.print("             'fi':      float  + 32-bit unsigned int");
    comm.print("             'fl':      float  + 64-bit unsigned int");
    comm.print("             'di':      double + 32-bit unsigned int");
    comm.print("             'dl':      double + 64-bit unsigned int");
    comm.print("             'qi':      quad   + 32-bit unsigned int (requires quadmath)");
    comm.print("             'ql':      quad   + 64-bit unsigned int (requires quadmath)");
    return 0;
  }

  // by default, use generic backend
  Backend::set_preferred_backend(PreferredBackend::generic);

  // by default, we use the generic backend
  PreferredBackend preferred_backend = PreferredBackend::generic;

  // parse process and thread grid dimensions
  Index nx(0), ny(0), mx(0), my(0), level(0);
  double min_runtime(0.0);
  int dt_it = 8008; // 1000*sizeof(DT_) + sizeof(IT_)
  if(!String(argv[1]).parse(mx) || (mx < 1u))
  {
    comm.print("ERROR: Failed to parse '" + String(argv[1]) + "' as process grid dimension 'mx' >= 1");
    return 1;
  }
  if(!String(argv[2]).parse(my) || (my < 1u))
  {
    comm.print("ERROR: Failed to parse '" + String(argv[2]) + "' as process grid dimension 'my' >= 1");
    return 1;
  }
  if(!String(argv[3]).parse(nx) || (nx < 1u))
  {
    comm.print("ERROR: Failed to parse '" + String(argv[3]) + "' as global grid dimension 'nx' >= 1");
    return 1;
  }
  if(!String(argv[4]).parse(ny) || (ny < 1u))
  {
    comm.print("ERROR: Failed to parse '" + String(argv[4]) + "' as global grid dimension 'ny' >= 1");
    return 1;
  }
  if(!String(argv[5]).parse(level))
  {
    comm.print("ERROR: Failed to parse '" + String(argv[5]) + "' as refinement level");
    return 1;
  }
  if(!String(argv[6]).parse(min_runtime) || (min_runtime < 1e-3))
  {
    comm.print("ERROR: Failed to parse '" + String(argv[6]) + "' as minimum runtime > 0");
    return 1;
  }
  if(argc > 7)
  {
    String sback(argv[7]);
    if(sback.compare_no_case("generic") == 0)
      preferred_backend = PreferredBackend::generic;
    else if(sback.compare("cuda") == 0)
      preferred_backend = PreferredBackend::cuda;
    else if(sback.compare("mkl") == 0)
      preferred_backend = PreferredBackend::mkl;
    else
    {
      comm.print("ERROR: Failed to parse '" + sback + "' as backend; must be one of: 'generic', 'cuda', 'mkl'");
      return 1;
    }
  }
  if(argc > 8)
  {
    int sdt = 0;
    int sit = 0;
    switch(argv[8][0])
    {
    case 'h':
    case 'H':
      sdt = 2;
      break;

    case 'f':
    case 'F':
      sdt = 4;
      break;

    case 'd':
    case 'D':
      sdt = 8;
      break;

    case 'q':
    case 'Q':
      sdt = 16;
      break;

    default:
      comm.print("ERROR: Failed to parse '" + String(argv[8]) + "' as data/index-type pair");
      return 1;
    }
    switch(argv[8][1])
    {
    case 'i':
    case 'I':
      sit = 4;
      break;

    case 'l':
    case 'L':
      sit = 8;
      break;

    default:
      comm.print("ERROR: Failed to parse '" + String(argv[8]) + "' as data/index-type pair");
      return 1;
    }
#if !defined(FEAT_HAVE_HALFMATH)
    if(sdt == 2)
    {
      comm.print("ERROR: Data-Type 'half' selected, but half precision is not available!"
        "\nNote: Half precision is only available when compiling with CUDA enabled.");
      return 1;
    }
#endif
#if !defined(FEAT_HAVE_QUADMATH)
    if(sdt == 16)
    {
      comm.print("ERROR: Data-Type 'quad' selected, but quadruple precision is not available!"
        "\nNote: Quadruple precision is only available when compiling with libquadmath under gcc.");
      return 1;
    }
#endif
    dt_it = 1000*sdt + sit;
  }

  // ensure that the process grid dimensions are valid
  if(mx * my != Index(comm.size()))
  {
    comm.print("ERROR: invalid process grid dimensions: got " + stringify(mx) + "*" + stringify(my)
      + " = " + stringify(mx*my) + " but expected " + stringify(comm.size()));
    return 1;
  }

  // ensure that the global grid dimensions are multiples of the process grid dimensions
  if((nx < mx) || (nx % mx != 0u))
  {
    comm.print("ERROR: global grid dimension nx = " + stringify(nx) + " must be a multiple of process grid dimension mx = " + stringify(mx));
    return 1;
  }
  if((ny < my) || (ny % my != 0u))
  {
    comm.print("ERROR: global grid dimension ny = " + stringify(ny) + " must be a multiple of process grid dimension my = " + stringify(my));
    return 1;
  }

  // ensure that the selected backend is actually available
#if !defined(FEAT_HAVE_CUDA)
  if(preferred_backend == PreferredBackend::cuda)
  {
    comm.print("ERROR: Backend 'cuda' selected but CUDA is not available");
    return 1;
  }
#endif
#if !defined(FEAT_HAVE_MKL)
  if(preferred_backend == PreferredBackend::mkl)
  {
    comm.print("ERROR: Backend 'mkl' selected but MKL is not available");
    return 1;
  }
#endif

  // run benchmark
  switch(dt_it)
  {
#ifdef FEAT_HAVE_HALFMATH
  case 2004: // half, int32
    return run<Half, std::uint32_t>(comm, mx, my, nx, ny, level, min_runtime, preferred_backend);

  case 2008: // half, int64
    return run<Half, std::uint64_t>(comm, mx, my, nx, ny, level, min_runtime, preferred_backend);
#endif // FEAT_HAVE_HALFMATH

  case 4004: // float, int32
    return run<float, std::uint32_t>(comm, mx, my, nx, ny, level, min_runtime, preferred_backend);

  case 4008: // float, int64
    return run<float, std::uint64_t>(comm, mx, my, nx, ny, level, min_runtime, preferred_backend);

  case 8004: // double, int32
    return run<double, std::uint32_t>(comm, mx, my, nx, ny, level, min_runtime, preferred_backend);

  case 8008: // double, int64
    return run<double, std::uint64_t>(comm, mx, my, nx, ny, level, min_runtime, preferred_backend);

#ifdef FEAT_HAVE_QUADMATH
  case 16004: // quadruple, int32
    return run<__float128, std::uint32_t>(comm, mx, my, nx, ny, level, min_runtime, preferred_backend);

  case 16008: // quadruple, int64
    return run<__float128, std::uint64_t>(comm, mx, my, nx, ny, level, min_runtime, preferred_backend);
#endif // FEAT_HAVE_QUADMATH

  default:
    comm.print("ERROR: invalid data/index-type selected");
    return 1;
  }
}

template<typename DT_, typename IT_>
int run(const Dist::Comm& comm, const Index mx, const Index my, const Index nx, const Index ny, const Index level, double min_runtime, PreferredBackend preferred_backend)
{
  typedef DT_ DataType;
  typedef IT_ IndexType;

#ifdef FEAT_HAVE_MPI
  comm.print("Number of MPI Processes..:" + stringify(comm.size()).pad_front(18));
#else
  comm.print("Number of MPI Processes..:             -N/A-");
#endif
#ifdef FEAT_HAVE_OMP
  comm.print("OMP Threads per Process..:" + stringify(omp_get_max_threads()).pad_front(18));
#else
  comm.print("OMP Threads per Process..:             -N/A-");
#endif
#ifdef FEAT_HAVE_MKL
  comm.print("MKL Threads per Process..:" + stringify(mkl_get_max_threads()).pad_front(18));
#else
  comm.print("MKL Threads per Process..:             -N/A-");
#endif

  // local mesh dimensions
  const Index lx = nx / mx;
  const Index ly = ny / my;

  // level refinement factor per dimension = 2^level
  const Index lref = Index(1) << (level);

  // print problem setup and some dimensional statistics
  comm.print("Total   Grid Dimensions..:" + stringify(nx*lref).pad_front(8) + " x" + stringify(ny*lref).pad_front(8));
  comm.print("Process Grid Dimensions..:" + stringify(mx*lref).pad_front(8) + " x" + stringify(my*lref).pad_front(8));
  comm.print("Local   Grid Dimensions..:" + stringify(lx*lref).pad_front(8) + " x" + stringify(ly*lref).pad_front(8));

  // declare mesh, trafo and space types
  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceType;

  StopWatch watch_mesh_create, watch_setup;
  watch_setup.start();
  watch_mesh_create.start();

  // create a structured base mesh
  std::unique_ptr<MeshNodeType> base_mesh_node;
  {
    Geometry::StructUnitCubeFactory<MeshType> base_mesh_factory(nx, ny);
    base_mesh_node = MeshNodeType::make_unique(base_mesh_factory.make_unique());
  }

  // create patch mesh
  std::unique_ptr<MeshNodeType> mesh_node;
  std::vector<int> neighbor_ranks;
  {
    Adjacency::Graph elems_at_rank = create_elems_at_rank(mx, my, nx, ny);
    mesh_node = base_mesh_node->extract_patch(neighbor_ranks, elems_at_rank, comm.rank());
  }

  // refine up to desired level
  for(Index lvl(0); lvl < level; ++lvl)
    mesh_node = mesh_node->refine_unique();

  watch_mesh_create.stop();

  //Geometry::ExportVTK<MeshType> vtk(*mesh_node->get_mesh());
  //vtk.write("parperf2-bench", comm);

  // get mesh, create trafo and space
  MeshType& mesh = *mesh_node->get_mesh();
  TrafoType trafo(mesh);
  SpaceType space(trafo);

  XASSERT(mesh.get_num_elements() == (lx*ly*lref*lref));

  comm.print("Local Number of Elements.:" + stringify(lx*ly*lref*lref).pad_front(18));
  comm.print("Total Number of Elements.:" + stringify(nx*ny*lref*lref).pad_front(18));

  const Index num_dofs_l = space.get_num_dofs();

  comm.print("Local Number of DOFs.....:" + stringify(num_dofs_l).pad_front(18));

  // declare local types
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalMatrixType;
  typedef LAFEM::DenseVector<DataType, IndexType> LocalVectorType;
  typedef LAFEM::VectorMirror<DataType, IndexType> MirrorType;

  // create the gate
  StopWatch watch_asm_gate;
  watch_asm_gate.start();
  Global::Gate<LocalVectorType, MirrorType> gate(comm);
  for(const auto& it: mesh_node->get_halo_map())
  {
    MirrorType mirror;
    Assembly::MirrorAssembler::assemble_mirror(mirror, space, *it.second);
    XASSERT(!mirror.empty());
    gate.push(it.first, std::move(mirror));
  }
  gate.compile(LocalVectorType(num_dofs_l));
  const LocalVectorType& vec_f = gate.get_freqs();
  watch_asm_gate.stop();

  const Index num_dofs_g = gate.get_num_global_dofs();
  comm.print("Total Number of DOFs.....:" + stringify(num_dofs_g).pad_front(18));

  // create a domain assembler
  StopWatch watch_asm_domain_asm;
  watch_asm_domain_asm.start();
  Assembly::DomainAssembler<TrafoType> domain_asm(trafo);
  domain_asm.compile_all_elements();
  watch_asm_domain_asm.stop();

  // create local matrix and assemble structure
  StopWatch watch_asm_matrix_sym;
  watch_asm_matrix_sym.start();
  LocalMatrixType matrix;
  Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
  watch_asm_matrix_sym.stop();

  const Index num_nze_l = matrix.used_elements();
  const Index num_nze_g = mx*my*num_nze_l;
  comm.print("Local Number of NZEs.....:" + stringify(num_nze_l).pad_front(18));
  comm.print("Total Number of NZEs.....:" + stringify(num_nze_g).pad_front(18));

  // assemble laplace matrix
  StopWatch watch_asm_matrix_num;
  watch_asm_matrix_num.start();
  Assembly::Common::LaplaceOperator laplace_op;
  Assembly::assemble_bilinear_operator_matrix_1(domain_asm, matrix, laplace_op, space, "gauss-legendre:3");
  watch_asm_matrix_num.stop();

  // interpolate start vectors
  StopWatch watch_asm_vectors;
  watch_asm_vectors.start();
  LocalVectorType vec_x, vec_t = matrix.create_vector_l();
  Analytic::Common::Q2BubbleFunction<2> q2_bubble_function;
  Assembly::Interpolator::project(vec_x, q2_bubble_function, space);
  vec_t.format(Math::sqrt(DataType(3.14159) /  DataType(num_dofs_g)));
  watch_asm_vectors.stop();

  // compute bytes to be allocated (just for statistics)
  std::size_t bytes_mesh_l = mesh_node->bytes() + base_mesh_node->bytes();
  std::size_t bytes_gate_l = gate.bytes();
  std::size_t bytes_vec_l = vec_x.bytes();
  std::size_t bytes_mat_l = matrix.bytes();
  std::size_t bytes_mesh_g(0), bytes_gate_g(0), bytes_vec_g(0), bytes_mat_g(0);
  comm.allreduce(&bytes_mesh_l, &bytes_mesh_g, std::size_t(1), Dist::op_sum);
  comm.allreduce(&bytes_gate_l, &bytes_gate_g, std::size_t(1), Dist::op_sum);
  comm.allreduce(&bytes_vec_l, &bytes_vec_g, std::size_t(1), Dist::op_sum);
  comm.allreduce(&bytes_mat_l, &bytes_mat_g, std::size_t(1), Dist::op_sum);
  std::size_t bytes_tot_g = bytes_mat_g + 2*bytes_vec_g + bytes_gate_g;

  comm.print("Local Mesh Size..........:" + stringify_bytes(bytes_mesh_l, 3, 18));
  comm.print("Total Mesh Size..........:" + stringify_bytes(bytes_mesh_g, 3, 18));
  comm.print("Local Gate Size..........:" + stringify_bytes(bytes_gate_l, 3, 18));
  comm.print("Total Gate Size..........:" + stringify_bytes(bytes_gate_g, 3, 18));
  comm.print("Local Vector Size........:" + stringify_bytes(bytes_vec_l, 3, 18));
  comm.print("Total Vector Size........:" + stringify_bytes(bytes_vec_g, 3, 18));
  comm.print("Local Matrix Size........:" + stringify_bytes(bytes_mat_l, 3, 18));
  comm.print("Total Matrix Size........:" + stringify_bytes(bytes_mat_g, 3, 18));
  comm.print("Total Benchmark Size.....:" + stringify_bytes(bytes_tot_g, 3, 18));

  comm.print("Selected Backend.........: " + stringify(preferred_backend));
  comm.print("DataType.................: " + Type::Traits<DataType>::name());
  comm.print("DataType Size............: " + stringify(sizeof(DataType)) + " Bytes");
  comm.print("IndexType................: " + Type::Traits<IndexType>::name());
  comm.print("IndexType Size...........: " + stringify(sizeof(IndexType)) + " Bytes");

  // set backend
  Backend::set_preferred_backend(preferred_backend);

  // compute Euclid/Frobenius norms of matrices and vectors; this is primarily done to ensure that
  // all the memory is touched on the selected backend before the actual power iteration loop starts
  // note that these norms depend on the process grid size, so they vary for different numbers of
  // ranks even if the total grid dimensions are the same; we only print these norms to stdout to
  // ensure that the compiler does not optimize the statements away thus preventing the first touch
  comm.print("\nPerforming first touch on backend memory...");
  comm.print("|A| = " + stringify_fp_sci(gate.sum(matrix.norm_frobenius())));
  comm.print("|x| = " + stringify_fp_sci(gate.sum(vec_x.norm2sqr())));
  comm.print("|t| = " + stringify_fp_sci(gate.sum(vec_t.norm2sqr())));
  comm.print("|f| = " + stringify_fp_sci(gate.sum(vec_f.norm2sqr())));
  // ensure that the mirrors are also touched by firing a useless sync
  gate.sync_0(vec_t);

  watch_setup.stop();

  // define statistics variables
  TimeStamp stamp_begin, stamp_end, stamp_1, stamp_2, stamp_3, stamp_4, stamp_5;
  IndexType iter = 0;
  KahanAccumulator<double> accum_time_matvec;
  KahanAccumulator<double> accum_time_sync_0;
  KahanAccumulator<double> accum_time_tridot;
  KahanAccumulator<double> accum_time_reduce;
  KahanAccumulator<double> accum_time_scale;
  KahanAccumulator<double> accum_time_sanity;
  double time_total = 0.0;
  int is_finished = 0;

  // register likwid markers
  LIKWID_MARKER_REGISTER("MatVecMult");
  LIKWID_MARKER_REGISTER("VecTriDot");
  LIKWID_MARKER_REGISTER("VecScale");
  LIKWID_NVMARKER_REGISTER("NV_MatVecMult");
  LIKWID_NVMARKER_REGISTER("NV_VecTriDot");
  LIKWID_NVMARKER_REGISTER("NV_VecScale");

  comm.print("\nRunning power iteration for at least " + stringify(min_runtime) + " seconds, please be patient...");
  comm.print_flush();

  // power iteration loop
  do
  {
    // --------------------------------------------------------------------------------------------
    // STEP 1: compute local norm of x
    // --------------------------------------------------------------------------------------------
    stamp_1.stamp();

    // compute local norm of vector x by triple dot
    LIKWID_MARKER_START("VecTriDot");
    LIKWID_NVMARKER_START("NV_VecTriDot");

    DataType loc_norm = vec_f.triple_dot(vec_x, vec_x);

    LIKWID_MARKER_STOP("VecTriDot");
    LIKWID_NVMARKER_STOP("NV_VecTriDot");

    // --------------------------------------------------------------------------------------------
    // STEP 2: sum local norms to obtain global norm of x
    // --------------------------------------------------------------------------------------------

    stamp_2.stamp();

    // sum dot up via gate and compute euclid norm
    DataType glob_norm = Math::sqrt(gate.sum(loc_norm));

    // --------------------------------------------------------------------------------------------
    // STEP 3: normalize x by scaling it by its inverse norm: t := x/|x|
    // --------------------------------------------------------------------------------------------
    stamp_3.stamp();

    // normalize x as t
    LIKWID_MARKER_START("VecScale");
    LIKWID_NVMARKER_START("NV_VecScale");

    vec_t.scale(vec_x, DataType(1) / glob_norm);

    LIKWID_MARKER_STOP("VecScale");
    LIKWID_NVMARKER_STOP("NV_VecScale");

    // --------------------------------------------------------------------------------------------
    // STEP 4: compute local matrix-vector product x := A*t
    // --------------------------------------------------------------------------------------------
    stamp_4.stamp();

    // perform matrix-vector product x := A*t
    LIKWID_MARKER_START("MatVecMult");
    LIKWID_NVMARKER_START("NV_MatVecMult");

    matrix.apply(vec_x, vec_t);

    LIKWID_MARKER_STOP("MatVecMult");
    LIKWID_NVMARKER_STOP("NV_MatVecMult");

    // --------------------------------------------------------------------------------------------
    // STEP 5: synchronize vector x as type-0 vector to obtain type-1 vector
    // --------------------------------------------------------------------------------------------
    stamp_5.stamp();

    /// synchronize type-0 vector x (aka "halo exchange")
    gate.sync_0(vec_x);

    stamp_end.stamp();

    // --------------------------------------------------------------------------------------------
    // STEP 6: post-processing
    // --------------------------------------------------------------------------------------------

    ++iter;

    //comm.print(stringify(iter).pad_front(5) + ": " + stringify_fp_sci(glob_norm));

    // accumulate runtimes
    accum_time_tridot += stamp_2.elapsed(stamp_1);
    accum_time_reduce += stamp_3.elapsed(stamp_2);
    accum_time_scale  += stamp_4.elapsed(stamp_3);
    accum_time_matvec += stamp_5.elapsed(stamp_4);
    accum_time_sync_0 += stamp_end.elapsed(stamp_5);

    // update time stamp
    time_total = stamp_end.elapsed(stamp_begin);

    // are we there yet? (this may yield different results on different ranks!!!)
    is_finished = (time_total < min_runtime ? 0 : 1);

    // rank 0 decides whether we're done or not
    comm.bcast(&is_finished, std::size_t(1), 0);
  } while(is_finished == 0);

  comm.print("\nBenchmarking finished; iterations performed: "+ stringify(iter));

  // switch back to generic backend
  Backend::set_preferred_backend(PreferredBackend::generic);

  // compute final norm of vector x
  DataType norm_x = Math::sqrt(gate.sum(vec_f.triple_dot(vec_x, vec_x)));
  comm.print("\nFinal norm of x = " + stringify_fp_sci(norm_x));

  // compute sanity check time by summing up the accumulated part times
  accum_time_sanity += accum_time_tridot.value;
  accum_time_sanity += accum_time_reduce.value;
  accum_time_sanity += accum_time_scale.value;
  accum_time_sanity += accum_time_matvec.value;
  accum_time_sanity += accum_time_sync_0.value;

  // compute minimum/maximum/sum times
  static const std::size_t num_times = 12u;
  double times_min[num_times], times_max[num_times], times_sum[num_times], times_avg[num_times], times_tot[num_times];
  times_min[ 0] = times_max[ 0] = times_sum[ 0] = accum_time_sanity.value;
  times_min[ 1] = times_max[ 1] = times_sum[ 1] = accum_time_matvec.value;
  times_min[ 2] = times_max[ 2] = times_sum[ 2] = accum_time_tridot.value;
  times_min[ 3] = times_max[ 3] = times_sum[ 3] = accum_time_scale.value;
  times_min[ 4] = times_max[ 4] = times_sum[ 4] = accum_time_sync_0.value;
  times_min[ 5] = times_max[ 5] = times_sum[ 5] = accum_time_reduce.value;
  times_min[ 6] = times_max[ 6] = times_sum[ 6] = watch_mesh_create.elapsed();
  times_min[ 7] = times_max[ 7] = times_sum[ 7] = watch_asm_gate.elapsed();
  times_min[ 8] = times_max[ 8] = times_sum[ 8] = watch_asm_domain_asm.elapsed();
  times_min[ 9] = times_max[ 9] = times_sum[ 9] = watch_asm_matrix_sym.elapsed();
  times_min[10] = times_max[10] = times_sum[10] = watch_asm_matrix_num.elapsed();
  times_min[11] = times_max[11] = times_sum[11] = watch_asm_vectors.elapsed();

  comm.allreduce(times_min, times_min, num_times, Dist::op_min);
  comm.allreduce(times_max, times_max, num_times, Dist::op_max);
  comm.allreduce(times_sum, times_sum, num_times, Dist::op_sum);

  // compute average times
  for(std::size_t i = 0; i < num_times; ++i)
  {
    times_avg[i] = times_sum[i] / double(comm.size());
    times_tot[i] = times_avg[i] / double(comm.size()); // only required for total bandwidth/flops
  }

  double time_setup = watch_setup.elapsed();
  comm.allreduce(&time_setup, &time_setup, std::size_t(1), Dist::op_max);

  // print setup runtime summary
  comm.print("\nSetup Runtime Summary:" + String("Minimum").pad_front(18) + String("Average").pad_front(25) + String("Maximum").pad_front(25));
  comm.print("Total Setup Runtime......:" + stringify_fp_fix(time_setup, 6, 14));
  comm.print("Mesh Creation Time.......:" +
    stringify_fp_fix(times_min[6], 6, 14) + " [" + stringify_fp_fix(100.*times_min[6]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[6], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[6]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_max[6], 6, 14) + " [" + stringify_fp_fix(100.*times_max[6]/time_setup, 3, 7) + "%]");
  comm.print("Gate Assembly Time.......:" +
    stringify_fp_fix(times_min[7], 6, 14) + " [" + stringify_fp_fix(100.*times_min[7]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[7], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[7]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_max[7], 6, 14) + " [" + stringify_fp_fix(100.*times_max[7]/time_setup, 3, 7) + "%]");
  comm.print("DomainAsm Setup Time.....:" +
    stringify_fp_fix(times_min[8], 6, 14) + " [" + stringify_fp_fix(100.*times_min[8]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[8], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[8]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_max[8], 6, 14) + " [" + stringify_fp_fix(100.*times_max[8]/time_setup, 3, 7) + "%]");
  comm.print("Sym. Matrix Assembly Time:" +
    stringify_fp_fix(times_min[9], 6, 14) + " [" + stringify_fp_fix(100.*times_min[9]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[9], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[9]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_max[9], 6, 14) + " [" + stringify_fp_fix(100.*times_max[9]/time_setup, 3, 7) + "%]");
  comm.print("Num. Matrix Assembly Time:" +
    stringify_fp_fix(times_min[10], 6, 14) + " [" + stringify_fp_fix(100.*times_min[10]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[10], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[10]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_max[10], 6, 14) + " [" + stringify_fp_fix(100.*times_max[10]/time_setup, 3, 7) + "%]");
  comm.print("Vector Assembly Time.....:" +
    stringify_fp_fix(times_min[11], 6, 14) + " [" + stringify_fp_fix(100.*times_min[11]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[11], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[11]/time_setup, 3, 7) + "%]" +
    stringify_fp_fix(times_max[11], 6, 14) + " [" + stringify_fp_fix(100.*times_max[11]/time_setup, 3, 7) + "%]");

  // print runtime summary
  comm.print("\nRuntime Summary:" + String("Minimum").pad_front(24) + String("Average").pad_front(25) + String("Maximum").pad_front(25));
  comm.print("Loop Iteration Runtime...:" + stringify_fp_fix(time_total/double(iter), 6, 14));
  comm.print("Total Loop Runtime.......:" + stringify_fp_fix(time_total, 6, 14));
  comm.print("Time Sum Sanity Check....:" +
    stringify_fp_fix(times_min[0], 6, 14) + " [" + stringify_fp_fix(100.*times_min[0]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[0], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[0]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_max[0], 6, 14) + " [" + stringify_fp_fix(100.*times_max[0]/time_total, 3, 7) + "%]");
  comm.print("Local MatVec Mult Time...:" +
    stringify_fp_fix(times_min[1], 6, 14) + " [" + stringify_fp_fix(100.*times_min[1]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[1], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[1]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_max[1], 6, 14) + " [" + stringify_fp_fix(100.*times_max[1]/time_total, 3, 7) + "%]");
  comm.print("Local Triple-Dot Time....:" +
    stringify_fp_fix(times_min[2], 6, 14) + " [" + stringify_fp_fix(100.*times_min[2]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[2], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[2]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_max[2], 6, 14) + " [" + stringify_fp_fix(100.*times_max[2]/time_total, 3, 7) + "%]");
  comm.print("Local Vector Scale Time..:" +
    stringify_fp_fix(times_min[3], 6, 14) + " [" + stringify_fp_fix(100.*times_min[3]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[3], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[3]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_max[3], 6, 14) + " [" + stringify_fp_fix(100.*times_max[3]/time_total, 3, 7) + "%]");
  comm.print("Sync-0 Halo-Exchange Time:" +
    stringify_fp_fix(times_min[4], 6, 14) + " [" + stringify_fp_fix(100.*times_min[4]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[4], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[4]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_max[4], 6, 14) + " [" + stringify_fp_fix(100.*times_max[4]/time_total, 3, 7) + "%]");
  comm.print("Global All-Reduce Time...:" +
    stringify_fp_fix(times_min[5], 6, 14) + " [" + stringify_fp_fix(100.*times_min[5]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_avg[5], 6, 14) + " [" + stringify_fp_fix(100.*times_avg[5]/time_total, 3, 7) + "%]" +
    stringify_fp_fix(times_max[5], 6, 14) + " [" + stringify_fp_fix(100.*times_max[5]/time_total, 3, 7) + "%]");

  // Note: it is kind of impossible to compute the memory throughput accurately, because it is unclear
  // whether a pure write access, which appears in both the left side of a matrix vector multiplication
  //   x_i <- sum_j A_ij * t_j
  // as well as the left side of the scale operation
  //   t_i <- x_i * s,
  // has to be counted as a single memory transfer (write only) or double memory transfer (read+write).
  // In the computations below we count a write access as both read+write, because tests on several
  // modern systems show that this is apparently what our current compilers make of this.

  // this is the raw memory throughput ignoring cache effects of multiplicand vector,
  // i.e. we're assuming that every t_j is loaded separately from memory, which is overly pessimistic
  const double mbw_matvec = double(iter) * (
    sizeof(IndexType) * (double(num_dofs_l) + double(num_nze_l)) + // row_ptr + col_idx
    sizeof(DataType) * (2.0*double(num_dofs_l) + 2.0*double(num_nze_l))); // (vec_x [read+write]) + (a_val [read] + vec_t [read])

  // this is the cache-corrected throughput assuming that each value t_j of the multiplicand
  // vector is loaded only once from memory, which is pretty optimistic
  const double mbw_matvec_cc = double(iter) * (
    sizeof(IndexType) * (double(num_dofs_l) + double(num_nze_l)) + // row_ptr + col_idx
    sizeof(DataType) * (3.0*double(num_dofs_l) + double(num_nze_l))); // (vec_x [read+write] + vec_t [read]) + (a_val [read])

  // f*v*v; two of three vectors are the same
  const double mbw_tridot = 2.0 * double(iter) * sizeof(DataType) * double(num_dofs_l);

  // vec_t [read+write] + vec_x [read]
  const double mbw_vscale = 3.0 * double(iter) * sizeof(DataType) * double(num_dofs_l);

  // print memory bandwidth summary
  comm.print("\nMemory Bandwidth Summary:" + String("Minimum").pad_front(10) + String("Average").pad_front(21) +
    String("Maximum").pad_front(21) + String("Total").pad_front(21));
  comm.print("MatVec Mult Bandwidth...:" +
    stringify_bytes(std::size_t(mbw_matvec / times_max[1]), 3, 10) + "/sec" +
    stringify_bytes(std::size_t(mbw_matvec / times_avg[1]), 3, 13) + "/sec" +
    stringify_bytes(std::size_t(mbw_matvec / times_min[1]), 3, 13) + "/sec" +
    stringify_bytes(std::size_t(mbw_matvec / times_tot[1]), 3, 13) + "/sec");
  comm.print("MatVec Mult CC Bandwidth:" +
    stringify_bytes(std::size_t(mbw_matvec_cc / times_max[1]), 3, 10) + "/sec" +
    stringify_bytes(std::size_t(mbw_matvec_cc / times_avg[1]), 3, 13) + "/sec" +
    stringify_bytes(std::size_t(mbw_matvec_cc / times_min[1]), 3, 13) + "/sec" +
    stringify_bytes(std::size_t(mbw_matvec_cc / times_tot[1]), 3, 13) + "/sec");
  comm.print("Triple-Dot Bandwidth....:" +
    stringify_bytes(std::size_t(mbw_tridot / times_max[2]), 3, 10) + "/sec" +
    stringify_bytes(std::size_t(mbw_tridot / times_avg[2]), 3, 13) + "/sec" +
    stringify_bytes(std::size_t(mbw_tridot / times_min[2]), 3, 13) + "/sec" +
    stringify_bytes(std::size_t(mbw_tridot / times_tot[2]), 3, 13) + "/sec");
  comm.print("Vector Scale Bandwidth..:" +
    stringify_bytes(std::size_t(mbw_vscale / times_max[3]), 3, 10) + "/sec" +
    stringify_bytes(std::size_t(mbw_vscale / times_avg[3]), 3, 13) + "/sec" +
    stringify_bytes(std::size_t(mbw_vscale / times_min[3]), 3, 13) + "/sec" +
    stringify_bytes(std::size_t(mbw_vscale / times_tot[3]), 3, 13) + "/sec");

  // Note: it is kind of impossible to compute floating point operation counts accurately, because it is unclear
  // whether a multiply-then-add operation, which appears in both the matrix-vector multiplication
  //   x_i <- sum_j A_ij * t_j
  // as well as the tri-dot operation
  //   r   <- sum_i f_i*(v_i*v_i)
  // has to be counted as two operations (separate multiply + add) or a single operation (fused-multiply-add aka "FMA").
  // In the computations below we count a multiply-then-add operation as separate multiply and add operations, i.e. as 2 FLOPs.

  // print floating point operations summary
  comm.print("\nFloating Point Summary:" + String("Minimum").pad_front(12) + String("Average").pad_front(21) +
    String("Maximum").pad_front(21) + String("Total").pad_front(21));
  const double flop_matvec = 2.0 * double(iter) * double(num_nze_l);  // 1*mult + 1*add
  const double flop_tridot = 3.0 * double(iter) * double(num_dofs_l); // 2*mult + 1*add
  const double flop_vscale =       double(iter) * double(num_dofs_l); // 1*mult
  comm.print("MatVec Mult Flops.......:" +
    stringify_flops(flop_matvec / times_max[1], 3, 10) + "/sec" +
    stringify_flops(flop_matvec / times_avg[1], 3, 11) + "/sec" +
    stringify_flops(flop_matvec / times_min[1], 3, 11) + "/sec" +
    stringify_flops(flop_matvec / times_tot[1], 3, 11) + "/sec");
  comm.print("Triple-Dot Flops........:" +
    stringify_flops(flop_tridot / times_max[2], 3, 10) + "/sec" +
    stringify_flops(flop_tridot / times_avg[2], 3, 11) + "/sec" +
    stringify_flops(flop_tridot / times_min[2], 3, 11) + "/sec" +
    stringify_flops(flop_tridot / times_tot[2], 3, 11) + "/sec");
  comm.print("Vector Scale Flops......:" +
    stringify_flops(flop_vscale / times_max[3], 3, 10) + "/sec" +
    stringify_flops(flop_vscale / times_avg[3], 3, 11) + "/sec" +
    stringify_flops(flop_vscale / times_min[3], 3, 11) + "/sec" +
    stringify_flops(flop_vscale / times_tot[3], 3, 11) + "/sec");

  // print out peak physical memory usage
  MemoryUsage mem_use;
  std::size_t peak_p = mem_use.get_peak_physical();
  comm.allreduce(&peak_p, &peak_p, std::size_t(1), Dist::op_sum);
  comm.print("\nPeak Total Memory Usage.: " + stringify_bytes(peak_p));

  comm.print("\nThank you for choosing parperf2-bench, have a nice and productive day!");
  return 0;
}
