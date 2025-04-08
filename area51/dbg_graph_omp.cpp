
#include <kernel/runtime.hpp>                              // for Runtime
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/kahan_accumulator.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/util/hash.hpp>
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>            // for RefinedUnitCubeFactory
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")
#include <kernel/space/lagrange2/element.hpp>              // the Lagrange-2 Element (aka "Q2")
#include <kernel/space/lagrange3/element.hpp>              // the Lagrange-3 Element (aka "Q3")
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/analytic/common.hpp>                      // for ExpBubbleFunction
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/domain_assembler_helpers.hpp>    // for Assembly::assemble_***
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR

#ifdef FEAT_HAVE_OMP
#include <omp.h>
#endif

// We are using FEAT, so use the namespace here.
using namespace FEAT;

typedef Index Trp_;
typedef Index Tci_;
typedef Index Tcp_;
typedef Index Tri_;

void transpose_sequential(const Index m, const Index n, const Trp_ row_ptr[], const Tci_ col_idx[], Tcp_ col_ptr[], Tri_ row_idx[])
{
  if((m <= Index(0)) || (n <= Index(0)))
    return;

  XASSERT(row_ptr != nullptr);
  XASSERT(col_idx != nullptr);
  XASSERT(col_ptr != nullptr);
  XASSERT(row_idx != nullptr);

  // get number of non-zero entries
  const Index nnze = row_ptr[m];

  // format column pointer
  for(Index i = 0; i <= n; ++i)
    col_ptr[i] = Tcp_(0);

  // count column indices
  col_ptr[0] = 0u;
  for(Index i = 0; i < nnze; ++i)
    ++col_ptr[col_idx[i]+1u];

  // perform inclusive scan to obtain column pointer
  for(Index i = 0; i < n; ++i)
    col_ptr[i+1u] += col_ptr[i];

  // insert indices
  for(Index i = 0; i < m; ++i)
  {
    for(Trp_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
    {
      row_idx[col_ptr[col_idx[j]]++] = Tri_(i);
    }
  }

  // restore column pointer
  for(Index i = n; i > 0; --i)
    col_ptr[i] = col_ptr[i-1u];
  col_ptr[0] = Tcp_(0);
}

int omp_transpose(const Index m, const Index n, const Trp_ row_ptr[], const Tci_ col_idx[],
  Tcp_ col_ptr[], Tri_ row_idx[], int max_threads = 0, double max_memory = 0.0)
{
  if((m <= Index(0)) || (n <= Index(0)))
    return 0;

  XASSERT(row_ptr != nullptr);
  XASSERT(col_idx != nullptr);
  XASSERT(col_ptr != nullptr);
  XASSERT(row_idx != nullptr);
  XASSERT(m > Index(max_threads));

#if !defined(FEAT_HAVE_OMP)
  (void)max_threads;
  (void)max_memory;
  transpose_sequential(m, n, row_ptr, col_idx, col_ptr, row_idx);
  return 1;
#else // defined(FEAT_HAVE_OMP)

  // get maximum number of permitted threads
  if(max_threads <= 0)
    max_threads = omp_get_max_threads();

  // did the caller specify a memory usage limit?
  if(max_memory > 1e-3)
  {
    const double bytes_x = double(m)*double(sizeof(Trp_)); // bx
    const double bytes_y = double(n)*double(sizeof(Tci_)); // by
    // temporary memory usage for T threads is (T-1)*by and the memory limit L specifies that
    //      T*by <= L*(bx+by)
    // <==> T    <= L*(bx+by)/by
    // <==> T    <= L*(bx/by + 1)
    int limit_t = int(max_memory*(bytes_x/bytes_y + 1.0));
    if(limit_t < max_threads)
      max_threads = limit_t;
  }

  // don't parallelize unless we have at least 100 elements per thread
  if((max_threads <= 1) || (m < Index(max_threads)))
  {
    transpose_sequential(m, n, row_ptr, col_idx, col_ptr, row_idx);
    return 1;
  }

  // initialize first entry of column pointer
  col_ptr[0u] = Tcp_(0);

  // allocate temporary column pointer array vector
  std::vector<std::vector<Tcp_>> tmp_col_ptr(max_threads);

  // allocate offset for each individual thread
  std::vector<Tcp_> tmp_col_off(std::size_t(max_threads+1u), Tcp_(0));

  // parallel OpenMP region
  FEAT_PRAGMA_OMP(parallel num_threads(max_threads) shared(tmp_col_ptr, tmp_col_off))
  {
    const int num_threads(omp_get_num_threads());
    const int thread_id(omp_get_thread_num());

    // get column pointer array for this thread
    tmp_col_ptr[thread_id].resize(n+2u, Tcp_(0));
    Tcp_* my_col_ptr = tmp_col_ptr.at(thread_id).data();

    // compute first and last row for this thread
    const Index i0 = (m*Index(thread_id  )) / Index(num_threads);
    const Index i1 = (m*Index(thread_id+1)) / Index(num_threads);

    // compute first and last column for this thread
    const Index j0 = ((n+1u)*Index(thread_id  )) / Index(num_threads);
    const Index j1 = ((n+1u)*Index(thread_id+1)) / Index(num_threads);

    // loop over all rows in a technically thread-parallel fashion
    for(Index i = i0; i < i1; ++i)
    {
      for(Trp_ k = row_ptr[i]; k < row_ptr[i+1]; ++k)
      {
        ASSERT(col_idx[k] >= Tcp_(0));
        ASSERT(col_idx[k] < Tcp_(n));
        ++my_col_ptr[col_idx[k]+1u];
      }
    }
    // make sure each thread is finished
    FEAT_PRAGMA_OMP(barrier)

    // perform scan of each thread chunk size
    Tcp_ chunk_nnze = Tcp_(0);

    // loop over all columns in a technically thread-parallel fashion
    for(Index j = j0; j < j1; ++j)
    {
      // sum up the column pointers of all threads
      for(int k = 0; k < num_threads; ++k)
      {
        col_ptr[j] += tmp_col_ptr[k][j];
      }

      col_ptr[j] = (chunk_nnze += col_ptr[j]);
    }

    // save last scan as offset for next thread
    tmp_col_off[thread_id+1u] = chunk_nnze;

    // make sure each thread is finished
    FEAT_PRAGMA_OMP(barrier)

    // perform sequential scan of thread column offsets
    FEAT_PRAGMA_OMP(master)
    {
      for(int l = 0u; l < num_threads; ++l)
      {
        tmp_col_off[l+1u] += tmp_col_off[l];
      }

      // this one was omitted in the parallel loop above
      tmp_col_off[num_threads] += col_ptr[n];
    } // omp master

      // make sure the master is finished
    FEAT_PRAGMA_OMP(barrier)

    // add column offset to each thread chunk in parallel
    //Tcp_ off = tmp_col_off[thread_id];
    for(Index j = j0; j < j1; ++j)
    {
      Tcp_ cptr = (col_ptr[j] += tmp_col_off[thread_id]);
      for(int l = 0; l < num_threads; ++l)
      {
        Tcp_ tmp = tmp_col_ptr[l][j+1];
        tmp_col_ptr[l][j+1] = cptr;
        cptr += tmp;
      }
    }

    // make sure each thread is finished
    FEAT_PRAGMA_OMP(barrier)

    for(Index i = i0; i < i1; ++i) // parallel for
    {
      for(Trp_ k = row_ptr[i]; k < row_ptr[i+1]; ++k)
      {
        row_idx[my_col_ptr[col_idx[k]+1u]++] = i;
      }
    }
  } // omp parallel

  return max_threads;
#endif
}

  // Here's our main function
int main(int argc, char* argv[])
{
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceType;

  if(argc <= 2)
  {
    std::cout << "USAGE: " << argv[0] << " <level> <steps> <threads...>\n";
    return 0;
  }

  std::deque<int> viargs;
  for(int i = 1; i < argc; ++i)
  {
    int x = 0;
    if(!String(argv[i]).parse(x))
    {
      std::cout << "ERROR: Failed to parse '" << argv[i] << "'\n!";
      return 1;
    }
    viargs.push_back(x);
  }

  int level = viargs[0];
  int num_steps = viargs[1];

  viargs.pop_front();
  viargs.pop_front();

  Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);

  MeshType mesh(mesh_factory);

  TrafoType trafo(mesh);
  SpaceType space(trafo);

  std::uint32_t crc_row_ptr(0), crc_col_idx(0), crc_dom_ptr(0), crc_img_idx(0);

  //StopWatch watch_1, watch_2, watch_3, watch_4, watch_5, watch_6, watch_7;
  TimeStamp stamp_0, stamp_1, stamp_2, stamp_3, stamp_4, stamp_5, stamp_6;
  KahanAccumulator<double> time_1, time_2, time_3, time_4, time_5, time_6;
  std::cout << " NT:    DofRender   Transpose     Compose        Sort      Matrix      Format  ColPtrCRC  RowIdxCRC  RowPtrCRC  ColIdxCRC\n";
  for(; !viargs.empty(); viargs.pop_front())
  {
    int num_threads = viargs.front();
#ifdef FEAT_HAVE_OMP
    omp_set_num_threads(num_threads);
#endif
    for(int step = 0; step < num_steps; ++step)
    {
      stamp_0.stamp();

      Adjacency::Graph dof_graph(Space::DofMappingRenderer::render(space));

      stamp_1.stamp();

      //Adjacency::Graph dof_support(Adjacency::RenderType::transpose, dof_graph);
      Adjacency::Graph dof_support(dof_graph.get_num_nodes_image(), dof_graph.get_num_nodes_domain(), dof_graph.get_num_indices());

      omp_transpose(dof_graph.get_num_nodes_domain(), dof_graph.get_num_nodes_image(),
        dof_graph.get_domain_ptr(), dof_graph.get_image_idx(), dof_support.get_domain_ptr(), dof_support.get_image_idx(), num_threads);

      crc_dom_ptr = Hash::crc32(sizeof(Index)*(dof_support.get_num_nodes_domain()+1u), dof_support.get_domain_ptr());
      crc_img_idx = Hash::crc32(sizeof(Index)*dof_support.get_num_indices(), dof_support.get_image_idx());

      stamp_2.stamp();

      Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify, dof_support, dof_graph);

      stamp_3.stamp();

      dof_adjactor.sort_indices();

      stamp_4.stamp();

      LAFEM::SparseMatrixCSR<double, Index> matrix(dof_adjactor);

      stamp_5.stamp();

      matrix.format();

      stamp_6.stamp();

      crc_row_ptr = Hash::crc32(sizeof(Index)*(matrix.rows()+1u), matrix.row_ptr());
      crc_col_idx = Hash::crc32(sizeof(Index)*matrix.used_elements(), matrix.col_ind());

      time_1 += stamp_1.elapsed(stamp_0);
      time_2 += stamp_2.elapsed(stamp_1);
      time_3 += stamp_3.elapsed(stamp_2);
      time_4 += stamp_4.elapsed(stamp_3);
      time_5 += stamp_5.elapsed(stamp_4);
      time_6 += stamp_6.elapsed(stamp_5);
    }

    std::cout << stringify(num_threads).pad_front(3) << ": ";
    std::cout << stringify_fp_fix(time_1.value, 6, 12);
    std::cout << stringify_fp_fix(time_2.value, 6, 12);
    std::cout << stringify_fp_fix(time_3.value, 6, 12);
    std::cout << stringify_fp_fix(time_4.value, 6, 12);
    std::cout << stringify_fp_fix(time_5.value, 6, 12);
    std::cout << stringify_fp_fix(time_6.value, 6, 12);
    std::cout << "   " << std::hex << crc_dom_ptr;
    std::cout << "   " << std::hex << crc_img_idx;;
    std::cout << "   " << std::hex << crc_row_ptr;
    std::cout << "   " << std::hex << crc_col_idx;;
    std::cout << std::endl;

    time_1.clear();
    time_2.clear();
    time_3.clear();
    time_4.clear();
    time_5.clear();
    time_6.clear();
  }

  MemoryUsage mem_use;
  std::cout /*<< "Peak Physical Memory: "*/ << "\n" << mem_use.get_formatted_memory_usage() << std::endl;

  return 0;
}
