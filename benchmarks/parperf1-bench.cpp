// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// ====================================================================================================================
// Linear Algebra Backend Parallel Performance Benchmark #1
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
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>

#ifdef FEAT_HAVE_OMP
#include <omp.h>
#endif

using namespace FEAT;

/*void dump_matrix(Index lx, Index ly)
{
  typedef Geometry::StructuredMesh<2, 2> MeshType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  Index ne[] = {lx, ly};
  MeshType mesh(ne);
  auto& vtx = mesh.get_vertex_set();
  for(Index iy = 0; iy <= ly; ++iy)
  {
    for(Index ix = 0; ix <= lx; ++ix)
    {
      auto&v = vtx[iy*(lx+1)+ix];
      v[0] = double(ix);
      v[1] = double(iy);
    }
  }

  TrafoType trafo(mesh);
  SpaceType space(trafo);

  LAFEM::SparseMatrixCSR<double, Index> matrix;
  Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
  matrix.format();

  Cubature::DynamicFactory cubfac("gauss-legendre:3");
  Assembly::Common::LaplaceOperator lop;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, lop, space, cubfac, 6.0);

  matrix.write_out(LAFEM::FileMode::fm_mtx, "quxbar.mtx", false);
}*/

/// creates a mirror for one of the four boundary edges
template<typename DT_, typename IT_>
LAFEM::VectorMirror<DT_, IT_> create_mirror(Index n, Index k, Index i, Index o)
{
  LAFEM::VectorMirror<DT_, IT_> mir(n, k);
  IT_* idx = mir.indices();
  FEAT_PRAGMA_OMP(parallel for)
  for(Index x = 0; x < k; ++x)
    idx[x] = IT_(o + i*x);
  return mir;
}

/// creates a gate with all required mirrors for this rank
template<typename DT_, typename IT_>
Global::Gate<LAFEM::DenseVector<DT_, IT_>, LAFEM::VectorMirror<DT_, IT_>> create_gate(const Dist::Comm& comm, Index mx, Index my, Index lx, Index ly)
{
  const int rank = comm.rank();
  const Index num_dofs_l = (lx+1u)*(ly+1u);
  const Index imx = Index(rank) % mx;
  const Index imy = Index(rank) / mx;

  Global::Gate<LAFEM::DenseVector<DT_, IT_>, LAFEM::VectorMirror<DT_, IT_>> gate(comm);

  // add lower neighbor
  if(imy > 0u)
    gate.push(int((imy-1)*mx + imx), create_mirror<DT_, IT_>(num_dofs_l, lx+1u, 1u, 0u));
  // add upper neighbor
  if(imy+1 < my)
    gate.push(int((imy+1)*mx + imx), create_mirror<DT_, IT_>(num_dofs_l, lx+1u, 1u, ly*(lx+1u)));
  // add left neighbor
  if(imx > 0u)
    gate.push(int(imy*mx + imx - 1), create_mirror<DT_, IT_>(num_dofs_l, ly+1u, lx+1u, 0u));
  // add right neighbor
  if(imx+1 < mx)
    gate.push(int(imy*mx + imx + 1), create_mirror<DT_, IT_>(num_dofs_l, ly+1u, lx+1u, lx));

  // add lower left neighbor
  if((imy > 0u) && (imx > 0u))
    gate.push(int((imy-1)*mx + imx - 1u), create_mirror<DT_, IT_>(num_dofs_l, 1u, 1u, 0u));
  // add lower right neighbor
  if((imy > 0u) && (imx+1 < mx))
    gate.push(int((imy-1)*mx + imx + 1u), create_mirror<DT_, IT_>(num_dofs_l, 1u, 1u, lx));
  // add upper left neighbor
  if((imy+1 < my) && (imx > 0u))
    gate.push(int((imy+1)*mx + imx - 1u), create_mirror<DT_, IT_>(num_dofs_l, 1u, 1u, ly*(lx+1u)));
  // add upper right neighbor
  if((imy+1 < my) && (imx+1 < mx))
    gate.push(int((imy+1)*mx + imx + 1u), create_mirror<DT_, IT_>(num_dofs_l, 1u, 1u, (ly+1)*(lx+1u)-1u));

  // compile gate
  gate.compile(LAFEM::DenseVector<DT_, IT_>((lx+1u)*(ly+1u)));
  return gate;
}

/// creates a (type-0) 2D 9-point stencil matrix (with Neumann boundary)
template<typename DT_, typename IT_>
LAFEM::SparseMatrixCSR<DT_, IT_> create_matrix(Index lx, Index ly)
{
  const Index num_dofs_l = (lx+1u)*(ly+1u);
  const Index num_nze_l = (3u*(lx+1)-2u)*(3u*(ly+1)-2u);

  LAFEM::SparseMatrixCSR<DT_, IT_> matrix(num_dofs_l, num_dofs_l, num_nze_l);

  IT_* row_ptr = matrix.row_ptr();
  IT_* col_idx = matrix.col_ind();
  DT_* val = matrix.val();

  row_ptr[0u] = 0u;

  // loop over Y-rows
  FEAT_PRAGMA_OMP(parallel for)
  for(Index iy = 0; iy <= ly; ++iy)
  {
    // loop over elements in Y-row
    for(Index ix = 0; ix <= lx; ++ix)
    {
      // compute row number
      const IT_ irow = IT_(iy*(lx+1)+ix);

      // compute row offset (don't ask)
      IT_ inze = IT_((iy > 0u ? 3*iy-1u : 0u)*(3u*(lx+1)-2u) + (ix > 0u ? 3*ix-1u : 0u)*((iy > 0u) && (iy < ly) ? 3u : 2u));
#ifdef FEAT_HAVE_OMP
      row_ptr[irow] = inze;
#else
      // sanity check
      XASSERT(row_ptr[irow] == inze);
#endif

      // add lower row couplings
      if(iy > 0u)
      {
        Index ikl = (iy-1u)*(lx+1u);
        // add lower left coupling
        if(ix > 0u)
        {
          col_idx[inze] = IT_(ikl + ix-1u);
          val[inze] = -DT_(2);
          ++inze;
        }
        // add lower coupling
        {
          col_idx[inze] = IT_(ikl + ix);
          val[inze] = -DT_((ix > 0 ? 1 : 0) + (ix < lx ? 1 : 0));
          ++inze;
        }
        // add lower right coupling
        if(ix < lx)
        {
          col_idx[inze] = IT_(ikl + ix+1u);
          val[inze] = -DT_(2);
          ++inze;
        }
      }

      // add center row couplings
      {
        Index ikl = (iy)*(lx+1u);
        // add left coupling
        if(ix > 0u)
        {
          col_idx[inze] = IT_(ikl + ix-1u);
          val[inze] = -DT_((iy > 0 ? 1 : 0) + (iy < ly ? 1 : 0));
          ++inze;
        }
        // add center coupling
        {
          col_idx[inze] = IT_(ikl + ix);
          val[inze] = DT_(4 * ((iy > 0 ? 1 : 0) + (iy < ly ? 1 : 0)) * ((ix > 0 ? 1 : 0) + (ix < lx ? 1 : 0)));
          ++inze;
        }
        // add right coupling
        if(ix < lx)
        {
          col_idx[inze] = IT_(ikl + ix+1u);
          val[inze] = -DT_((iy > 0 ? 1 : 0) + (iy < ly ? 1 : 0));
          ++inze;
        }
      }

      // add upper row couplings
      if(iy < ly)
      {
        Index ikl = (iy+1u)*(lx+1u);
        // add lower left coupling
        if(ix > 0u)
        {
          col_idx[inze] = IT_(ikl + ix-1u);
          val[inze] = -DT_(2);
          ++inze;
        }
        // add lower coupling
        {
          col_idx[inze] = IT_(ikl + ix);
          val[inze] = -DT_((ix > 0 ? 1 : 0) + (ix < lx ? 1 : 0));
          ++inze;
        }
        // add lower right coupling
        if(ix < lx)
        {
          col_idx[inze] = IT_(ikl + ix+1u);
          val[inze] = -DT_(2);
          ++inze;
        }
      }

#ifndef FEAT_HAVE_OMP
      row_ptr[irow+1u] = inze;
#endif
    } // next ix
  } // next iy

#ifdef FEAT_HAVE_OMP
  row_ptr[num_dofs_l] = IT_(num_nze_l);
#else
  XASSERT(row_ptr[num_dofs_l] == IT_(num_nze_l));
#endif

  return matrix;
}
/// creates a 9-point stencil matrix
template<typename DT_, typename IT_>
LAFEM::DenseVector<DT_, IT_> create_vector_init(int rank, Index mx, Index my, Index lx, Index ly)
{
  LAFEM::DenseVector<DT_, IT_> vector((lx+1u)*(ly+1u));

  DT_* val = vector.elements();

  const Index jy = rank / mx;
  const Index jx = rank % mx;
  const DT_ imnx = DT_(1) / DT_(mx*lx+1);
  const DT_ imny = DT_(1) / DT_(my*ly+1);

  // loop over Y-rows
  FEAT_PRAGMA_OMP(parallel for)
  for(Index iy = 0; iy <= ly; ++iy)
  {
    // compute global Y-coordinate
    DT_ y = imny * DT_(jy*ly + iy);

    // loop over elements in Y-row
    for(Index ix = 0; ix <= lx; ++ix)
    {
      // compute global X-coordinate
      DT_ x = imnx * DT_(jx*lx + ix);

      // set value u(x,y) := x*(1-x)*y*(1-y)
      val[iy*(lx+1)+ix] = x*(DT_(1) - x)*y*(DT_(1) - y);
    } // next ix
  } // next iy

  return vector;
}

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

enum class SelectedBackend
{
  none,
  generic,
  cuda,
  mkl
};

// forward declaration
template<typename DT_, typename IT_>
int run(const Dist::Comm& comm, const Index mx, const Index my, const Index nx, const Index ny, double min_runtime, SelectedBackend selected_backend);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  Dist::Comm comm(Dist::Comm::world());

  if(argc < 6)
  {
    comm.print("\nUSAGE: parperf-bench-1 <mx> <my> <nx> <ny> <min_time> [<backend>] [<dt-it>]\n");
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

  // by default, we don't use any backend
  SelectedBackend selected_backend = SelectedBackend::none;

  // parse process and thread grid dimensions
  Index nx(0), ny(0), mx(0), my(0);
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
  if(!String(argv[5]).parse(min_runtime) || (min_runtime < 1e-3))
  {
    comm.print("ERROR: Failed to parse '" + String(argv[5]) + "' as minimum runtime > 0");
    return 1;
  }
  if(argc > 6)
  {
    String sback(argv[6]);
    if(sback.compare_no_case("none") == 0)
      selected_backend = SelectedBackend::none;
    else if(sback.compare_no_case("generic") == 0)
      selected_backend = SelectedBackend::generic;
    else if(sback.compare("cuda") == 0)
      selected_backend = SelectedBackend::cuda;
    else if(sback.compare("mkl") == 0)
      selected_backend = SelectedBackend::mkl;
    else
    {
      comm.print("ERROR: Failed to parse '" + sback + "' as backend; must be one of: 'none', 'generic', 'cuda', 'mkl'");
      return 1;
    }
  }
  if(argc > 7)
  {
    int sdt = 0;
    int sit = 0;
    switch(argv[7][0])
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
      comm.print("ERROR: Failed to parse '" + String(argv[7]) + "' as data/index-type pair");
      return 1;
    }
    switch(argv[7][1])
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
      comm.print("ERROR: Failed to parse '" + String(argv[7]) + "' as data/index-type pair");
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
  if(selected_backend == SelectedBackend::cuda)
  {
    comm.print("ERROR: Backend 'cuda' selected but CUDA is not available");
    return 1;
  }
#endif
#if !defined(FEAT_HAVE_MKL)
  if(selected_backend == SelectedBackend::mkl)
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
    return run<Half, std::uint32_t>(comm, mx, my, nx, ny, min_runtime, selected_backend);

  case 2008: // half, int64
    return run<Half, std::uint64_t>(comm, mx, my, nx, ny, min_runtime, selected_backend);
#endif // FEAT_HAVE_HALFMATH

  case 4004: // float, int32
    return run<float, std::uint32_t>(comm, mx, my, nx, ny, min_runtime, selected_backend);

  case 4008: // float, int64
    return run<float, std::uint64_t>(comm, mx, my, nx, ny, min_runtime, selected_backend);

  case 8004: // double, int32
    return run<double, std::uint32_t>(comm, mx, my, nx, ny, min_runtime, selected_backend);

  case 8008: // double, int64
    return run<double, std::uint64_t>(comm, mx, my, nx, ny, min_runtime, selected_backend);

#ifdef FEAT_HAVE_QUADMATH
  case 16004: // quadruple, int32
    return run<__float128, std::uint32_t>(comm, mx, my, nx, ny, min_runtime, selected_backend);

  case 16008: // quadruple, int64
    return run<__float128, std::uint64_t>(comm, mx, my, nx, ny, min_runtime, selected_backend);
#endif // FEAT_HAVE_QUADMATH

  default:
    comm.print("ERROR: invalid data/index-type selected");
    return 1;
  }
}

template<typename DT_, typename IT_>
int run(const Dist::Comm& comm, const Index mx, const Index my, const Index nx, const Index ny, double min_runtime, SelectedBackend selected_backend)
{
  typedef DT_ DataType;
  typedef IT_ IndexType;

  // local mesh dimensions
  const Index lx = nx / mx;
  const Index ly = ny / my;

  // matrix dimensions
  const Index num_dofs_l = (lx+1u)*(ly+1u);
  const Index num_nze_l = (3u*(lx+1)-2u)*(3u*(ly+1)-2u);
  const Index num_dofs_g = (nx+1)*(ny+1);
  const Index num_nze_g = num_nze_l*mx*my;

  // compute bytes to be allocated (just for statistics)
  const std::size_t bytes_vec_l = num_dofs_l*sizeof(DataType);
  const std::size_t bytes_vec_g = bytes_vec_l*mx*my;
  const std::size_t bytes_mat_l = sizeof(DataType)*num_nze_l + sizeof(IndexType)*(num_dofs_l+num_nze_l);
  const std::size_t bytes_mat_g = bytes_mat_l*mx*my;
  // total size: 1x matrix + 3x vector[x+t+f] + mirrors + mirror buffers
  const std::size_t bytes_tot_g = bytes_mat_g + 3*bytes_vec_g + (sizeof(DataType)+sizeof(IndexType))*((mx-1)*ly + (my-1)*lx);

  PreferredBackend preferred_backend = PreferredBackend::generic;
  switch(selected_backend)
  {
  case SelectedBackend::none:
    preferred_backend = PreferredBackend::generic;
    break;

  case SelectedBackend::generic:
    preferred_backend = PreferredBackend::generic;
    break;

  case SelectedBackend::cuda:
    preferred_backend = PreferredBackend::cuda;
    break;

  case SelectedBackend::mkl:
    preferred_backend = PreferredBackend::mkl;
    break;
  }

  // print problem setup and some dimensional statistics
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
  comm.print("Total   Grid Dimensions..:" + stringify(nx).pad_front(8) + " x" + stringify(ny).pad_front(8));
  comm.print("Process Grid Dimensions..:" + stringify(mx).pad_front(8) + " x" + stringify(my).pad_front(8));
  comm.print("Local   Grid Dimensions..:" + stringify(lx).pad_front(8) + " x" + stringify(ly).pad_front(8));
  comm.print("Local Number of Elements.:" + stringify(lx*ly).pad_front(18));
  comm.print("Total Number of Elements.:" + stringify(nx*ny).pad_front(18));
  comm.print("Local Number of DOFs.....:" + stringify(num_dofs_l).pad_front(18));
  comm.print("Total Number of DOFs.....:" + stringify(num_dofs_g).pad_front(18));
  comm.print("Local Number of NZEs.....:" + stringify(num_nze_l).pad_front(18));
  comm.print("Total Number of NZEs.....:" + stringify(num_nze_g).pad_front(18));
  comm.print("Local Vector Size........:" + stringify_bytes(bytes_vec_l, 3, 18));
  comm.print("Total Vector Size........:" + stringify_bytes(bytes_vec_g, 3, 18));
  comm.print("Local Matrix Size........:" + stringify_bytes(bytes_mat_l, 3, 18));
  comm.print("Total Matrix Size........:" + stringify_bytes(bytes_mat_g, 3, 18));
  comm.print("Total Benchmark Size.....:" + stringify_bytes(bytes_tot_g, 3, 18));
  if(selected_backend == SelectedBackend::none)
    comm.print("Selected Backend.........: none");
  else
    comm.print("Selected Backend.........: " + stringify(preferred_backend));
  comm.print("DataType.................: " + Type::Traits<DataType>::name());
  comm.print("DataType Size............: " + stringify(sizeof(DataType)) + " Bytes");
  comm.print("IndexType................: " + Type::Traits<IndexType>::name());
  comm.print("IndexType Size...........: " + stringify(sizeof(IndexType)) + " Bytes");

  // declare local types
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalMatrixType;
  typedef LAFEM::DenseVector<DataType, IndexType> LocalVectorType;
  typedef LAFEM::VectorMirror<DataType, IndexType> MirrorType;

  StopWatch watch_setup;
  watch_setup.start();

  // create the gate
  Global::Gate<LocalVectorType, MirrorType> gate = create_gate<DataType, IndexType>(comm, mx, my, lx, ly);
  const LocalVectorType& vec_f = gate.get_freqs();

  // create local matrix
  LocalMatrixType matrix = create_matrix<DataType, IndexType>(lx, ly);
  //matrix.write_out(LAFEM::FileMode::fm_mtx, "A.mtx");

  // create vectors (their content is more or less arbitrary)
  LocalVectorType vec_x = create_vector_init<DataType, IndexType>(comm.rank(), mx, my, lx, ly);
  LocalVectorType vec_t = matrix.create_vector_l();
  //vec_x.format(Math::sqrt(DataType(1) / DataType(num_dofs_g)));
  vec_t.format(Math::sqrt(DataType(3.14159) /  DataType(num_dofs_g)));

  // set backend
  Backend::set_preferred_backend(preferred_backend);

  // get vector and matrix arrays for raw implementation
  DataType* val_x = vec_x.elements();
  DataType* val_t = vec_t.elements();
  const DataType* val_f = vec_f.elements();
  const DataType* val_a = matrix.val();
  const IndexType* row_ptr = matrix.row_ptr();
  const IndexType* col_idx = matrix.col_ind();

  // compute Euclid/Frobenius norms of matrices and vectors; this is primarily done to ensure that
  // all the memory is touched on the selected backend before the actual power iteration loop starts
  // note that these norms depend on the process grid size, so they vary for different numbers of
  // ranks even if the total grid dimensions are the same; we only print these norms to stdout to
  // ensure that the compiler does not optimize the statements away thus preventing the first touch
  if(selected_backend != SelectedBackend::none)
  {
    comm.print("\nPerforming first touch on backend memory...");
    comm.print("|A| = " + stringify_fp_sci(gate.sum(matrix.norm_frobenius())));
    comm.print("|x| = " + stringify_fp_sci(gate.sum(vec_x.norm2sqr())));
    comm.print("|t| = " + stringify_fp_sci(gate.sum(vec_t.norm2sqr())));
    comm.print("|f| = " + stringify_fp_sci(gate.sum(vec_f.norm2sqr())));
    // ensure that the mirrors are also touched by firing a useless sync
    gate.sync_0(vec_t);
  }

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

  // register likwid markers; for none backend, we have to register the markers inside an OpenMP parallel region
  if(selected_backend == SelectedBackend::none)
  {
    FEAT_PRAGMA_OMP(parallel)
    {
      LIKWID_MARKER_THREADINIT;
      LIKWID_MARKER_REGISTER("MatVecMult");
      LIKWID_MARKER_REGISTER("VecTriDot");
      LIKWID_MARKER_REGISTER("VecScale");
      LIKWID_NVMARKER_REGISTER("NV_MatVecMult");
      LIKWID_NVMARKER_REGISTER("NV_VecTriDot");
      LIKWID_NVMARKER_REGISTER("NV_VecScale");
    }
  }
  else // some backend
  {
    LIKWID_MARKER_REGISTER("MatVecMult");
    LIKWID_MARKER_REGISTER("VecTriDot");
    LIKWID_MARKER_REGISTER("VecScale");
    LIKWID_NVMARKER_REGISTER("NV_MatVecMult");
    LIKWID_NVMARKER_REGISTER("NV_VecTriDot");
    LIKWID_NVMARKER_REGISTER("NV_VecScale");
  }

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
    DataType loc_norm = DataType(0);
    if(selected_backend == SelectedBackend::none)
    {
      FEAT_PRAGMA_OMP(parallel)
      {
        LIKWID_MARKER_START("VecTriDot");
        LIKWID_NVMARKER_START("NV_VecTriDot");

        FEAT_PRAGMA_OMP(for reduction(+:loc_norm))
        for(Index i = 0; i < num_dofs_l; ++i)
        {
          loc_norm += val_x[i]*val_x[i]*val_f[i];
        }

        LIKWID_MARKER_STOP("VecTriDot");
        LIKWID_NVMARKER_STOP("NV_VecTriDot");
      } // pragma omp parallel
    }
    else // use preferred backend
    {
      LIKWID_MARKER_START("VecTriDot");
      LIKWID_NVMARKER_START("NV_VecTriDot");

      loc_norm = vec_f.triple_dot(vec_x, vec_x);

      LIKWID_MARKER_STOP("VecTriDot");
      LIKWID_NVMARKER_STOP("NV_VecTriDot");
    }

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
    if(selected_backend == SelectedBackend::none)
    {
      const DataType s = DataType(1) / glob_norm;
      FEAT_PRAGMA_OMP(parallel)
      {
        LIKWID_MARKER_START("VecScale");
        LIKWID_NVMARKER_START("NV_VecScale");

        FEAT_PRAGMA_OMP(for)
        for(Index i = 0; i < num_dofs_l; ++i)
        {
          val_t[i] = s*val_x[i];
        }

        LIKWID_MARKER_STOP("VecScale");
        LIKWID_NVMARKER_STOP("NV_VecScale");
      } // pragma omp parallel
    }
    else // use preferred backend
    {
      LIKWID_MARKER_START("VecScale");
      LIKWID_NVMARKER_START("NV_VecScale");

      vec_t.scale(vec_x, DataType(1) / glob_norm);

      LIKWID_MARKER_STOP("VecScale");
      LIKWID_NVMARKER_STOP("NV_VecScale");
    }

    // --------------------------------------------------------------------------------------------
    // STEP 4: compute local matrix-vector product x := A*t
    // --------------------------------------------------------------------------------------------
    stamp_4.stamp();

    // perform matrix-vector product x := A*t
    if(selected_backend == SelectedBackend::none)
    {
      FEAT_PRAGMA_OMP(parallel)
      {
        LIKWID_MARKER_START("MatVecMult");
        LIKWID_NVMARKER_START("NV_MatVecMult");

        FEAT_PRAGMA_OMP(for)
        for(Index i = 0; i < num_dofs_l; ++i)
        {
          DataType r = DataType(0);
          for(IndexType j = row_ptr[i]; j < row_ptr[i+1]; ++j)
            r += val_a[j] * val_t[col_idx[j]];
          val_x[i] = r;
        }

        LIKWID_MARKER_STOP("MatVecMult");
        LIKWID_NVMARKER_STOP("NV_MatVecMult");
      } // pragma omp parallel
    }
    else // use preferred backend
    {
      LIKWID_MARKER_START("MatVecMult");
      LIKWID_NVMARKER_START("NV_MatVecMult");

      matrix.apply(vec_x, vec_t);

      LIKWID_MARKER_STOP("MatVecMult");
      LIKWID_NVMARKER_STOP("NV_MatVecMult");
    }

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
  times_min[0] = times_max[0] = times_sum[0] = accum_time_sanity.value;
  times_min[1] = times_max[1] = times_sum[1] = accum_time_matvec.value;
  times_min[2] = times_max[2] = times_sum[2] = accum_time_tridot.value;
  times_min[3] = times_max[3] = times_sum[3] = accum_time_scale.value;
  times_min[4] = times_max[4] = times_sum[4] = accum_time_sync_0.value;
  times_min[5] = times_max[5] = times_sum[5] = accum_time_reduce.value;

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

  // print runtime summary
  comm.print("\nRuntime Summary:" + String("Minimum").pad_front(24) + String("Average").pad_front(25) + String("Maximum").pad_front(25));
  comm.print("Total Setup Runtime......:" + stringify_fp_fix(time_setup, 6, 14));
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

  comm.print("\nThank you for choosing parperf1-bench, have a nice and productive day!");
  return 0;
}
