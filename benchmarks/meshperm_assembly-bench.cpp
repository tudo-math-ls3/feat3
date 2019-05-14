// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/assembly/burgers_assembly_job.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/basic_assembly_jobs.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace MeshPermAssemblyBench
{
  using namespace FEAT;

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> ScalarMatrixType;
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> ScalarVectorType;

  /// minimum numeric assembly count
  int min_num_asm_count = 3;

  /// minimum matrix-vector product count
  int min_mat_vec_count = 10;

  /// minimum numeric assembly time in seconds
  double min_num_asm_time = 1.0;

  /// minimum matrix-vector product time in seconds
  double min_mat_vec_time = 1.0;

  struct BenchResults
  {
    int num_worker_threads;
    /// permutation time
    double permute_time;
    /// total symbolic assembly time in seconds
    double sym_asm_time;
    /// total numeric assembly time in seconds
    double num_asm_time;
    /// total matrix-vector product time in seconds
    double mat_vec_time;
    /// number of numeric assemblies performed
    int num_asm_count;
    /// number of matrix-vector products performed
    int mat_vec_count;

    /// threaded total time
    double thread_time_total;
    double thread_time_assemble;
    double thread_time_wait;

    BenchResults() :
      num_worker_threads(1),
      permute_time(0.0),
      sym_asm_time(0.0), num_asm_time(0.0), mat_vec_time(0.0),
      num_asm_count(0), mat_vec_count(0),
      thread_time_total(0.0),
      thread_time_assemble(0.0),
      thread_time_wait(0.0)
    {
    }
  }; // struct BenchResults

  template<typename Shape_, int nc_, typename CT_>
  void bench_poisson(Geometry::ConformalMesh<Shape_, nc_, CT_>& mesh, BenchResults& bres)
  {
    typedef Geometry::ConformalMesh<Shape_, nc_, CT_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    TrafoType trafo(mesh);
    SpaceType space(trafo);

    ScalarMatrixType matrix;

    // assemble matrix structure
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    Cubature::DynamicFactory cubature("gauss-legendre:2");
    Assembly::Common::LaplaceOperator laplace_op;

    bres.num_asm_count = 0;
    bres.num_asm_time = 0.0;
    TimeStamp stamp_1, stamp_2;

    // numeric assembly loop
    do
    {
      matrix.format();
      Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, laplace_op, space, cubature);
      stamp_2.stamp();
      bres.num_asm_time = stamp_2.elapsed(stamp_1);
      bres.num_asm_count += 1;
    } while((bres.num_asm_count < min_num_asm_count) || (bres.num_asm_time < min_num_asm_time));
  }

  template<typename Shape_, int nc_, typename CT_>
  void bench_poisson_new(Geometry::ConformalMesh<Shape_, nc_, CT_>& mesh, BenchResults& bres, std::size_t nthreads)
  {
    typedef Geometry::ConformalMesh<Shape_, nc_, CT_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    TrafoType trafo(mesh);
    SpaceType space(trafo);

    ScalarMatrixType matrix;

    Assembly::ThreadingStrategy strategy = Assembly::ThreadingStrategy::layered;
    if(mesh.get_mesh_permutation().get_strategy() == Geometry::PermutationStrategy::colored)
      strategy = Assembly::ThreadingStrategy::colored;
    if(nthreads == std::size_t(0))
      strategy = Assembly::ThreadingStrategy::single;

    Assembly::DomainAssembler<TrafoType> dom_asm(trafo);
    dom_asm.set_threading_strategy(strategy);
    dom_asm.set_max_worker_threads(nthreads);
    dom_asm.compile_all_elements();
    bres.num_worker_threads = int(dom_asm.get_num_worker_threads());

    // assemble matrix structure
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
    Assembly::Common::LaplaceOperator laplace_op;
    Assembly::BilinearOperatorMatrixAssemblyJob1<Assembly::Common::LaplaceOperator, ScalarMatrixType, SpaceType>
      asm_job(laplace_op, matrix, space, "gauss-legendre:2");

    bres.num_asm_count = 0;
    bres.num_asm_time = 0.0;
    TimeStamp stamp_1, stamp_2;

    // numeric assembly loop
    do
    {
      matrix.format();
      //dom_asm.assemble_bilinear_operator_matrix_1(matrix, laplace_op, space, "gauss-legendre:2");
      dom_asm.assemble(asm_job);
      stamp_2.stamp();
      bres.num_asm_time = stamp_2.elapsed(stamp_1);
      bres.num_asm_count += 1;
    } while((bres.num_asm_count < min_num_asm_count) || (bres.num_asm_time < min_num_asm_time));

    const double sc = 1e-6 / double(bres.num_asm_count);

    // compute assembly timings
    auto ts = dom_asm.reduce_thread_stats();
    bres.thread_time_total    = sc * double(ts.micros_total);
    bres.thread_time_assemble = sc * double(ts.micros_assemble);
    bres.thread_time_wait     = sc * double(ts.micros_wait);
  }

  template<typename Shape_, int nc_, typename CT_>
  void bench_burgers(Geometry::ConformalMesh<Shape_, nc_, CT_>& mesh, BenchResults& bres)
  {
    typedef Geometry::ConformalMesh<Shape_, nc_, CT_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

    TrafoType trafo(mesh);
    SpaceType space(trafo);

    static constexpr int dim = Shape_::dimension;

    // assemble matrix structure
    LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, dim> matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    // interpolate convection vector
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> vector;
    Analytic::Common::CosineWaveFunction<dim> cosine_wave;
    Analytic::Gradient<decltype(cosine_wave)> conv_func(cosine_wave);
    Assembly::Interpolator::project(vector, conv_func, space);

    Cubature::DynamicFactory cubature("gauss-legendre:3");

    // setup burgers
    Assembly::BurgersAssembler<DataType, IndexType, dim> burgers;
    burgers.nu = burgers.beta = burgers.frechet_beta = burgers.theta = 1.0;

    bres.num_asm_count = 0;
    bres.num_asm_time = 0.0;
    TimeStamp stamp_1, stamp_2;

    // numeric assembly loop
    do
    {
      matrix.format();
      burgers.assemble_matrix(matrix, vector, space, cubature);
      stamp_2.stamp();
      bres.num_asm_time = stamp_2.elapsed(stamp_1);
      bres.num_asm_count += 1;
    } while((bres.num_asm_count < min_num_asm_count) || (bres.num_asm_time < min_num_asm_time));
  }

  template<typename Shape_, int nc_, typename CT_>
  void bench_burgers_new(Geometry::ConformalMesh<Shape_, nc_, CT_>& mesh, BenchResults& bres, std::size_t nthreads)
  {
    typedef Geometry::ConformalMesh<Shape_, nc_, CT_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

    TrafoType trafo(mesh);
    SpaceType space(trafo);

    static constexpr int dim = Shape_::dimension;

    Assembly::ThreadingStrategy strategy = Assembly::ThreadingStrategy::layered;
    if(mesh.get_mesh_permutation().get_strategy() == Geometry::PermutationStrategy::colored)
      strategy = Assembly::ThreadingStrategy::colored;
    if(nthreads == std::size_t(0))
      strategy = Assembly::ThreadingStrategy::single;

    Assembly::DomainAssembler<TrafoType> dom_asm(trafo);
    dom_asm.set_threading_strategy(strategy);
    dom_asm.set_max_worker_threads(nthreads);
    dom_asm.compile_all_elements();
    bres.num_worker_threads = int(dom_asm.get_num_worker_threads());

    // assemble matrix structure
    LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, dim> matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    // interpolate convection vector
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> vector;
    Analytic::Common::CosineWaveFunction<dim> cosine_wave;
    Analytic::Gradient<decltype(cosine_wave)> conv_func(cosine_wave);
    Assembly::Interpolator::project(vector, conv_func, space);

    // setup burgers
    Assembly::BurgersBlockedMatrixAssemblyJob<decltype(matrix), SpaceType>
      burgers(matrix, vector, space, "gauss-legendre:3");
    burgers.nu = burgers.beta = burgers.frechet_beta = burgers.theta = 1.0;

    bres.num_asm_count = 0;
    bres.num_asm_time = 0.0;
    TimeStamp stamp_1, stamp_2;

    // numeric assembly loop
    do
    {
      matrix.format();
      dom_asm.assemble(burgers);
      stamp_2.stamp();
      bres.num_asm_time = stamp_2.elapsed(stamp_1);
      bres.num_asm_count += 1;
    } while((bres.num_asm_count < min_num_asm_count) || (bres.num_asm_time < min_num_asm_time));

    const double sc = 1e-6 / double(bres.num_asm_count);

    // compute assembly timings
    auto ts = dom_asm.reduce_thread_stats();
    bres.thread_time_total    = sc * double(ts.micros_total);
    bres.thread_time_assemble = sc * double(ts.micros_assemble);
    bres.thread_time_wait     = sc * double(ts.micros_wait);
  }

  template<typename Mesh_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
  {
    // parse levels
    Index lvl_min(0);
    Index lvl_max(1);
    args.parse("level", lvl_max, lvl_min);

    // parse threads
    std::vector<std::size_t> num_threads;
    num_threads.push_back(0u);
    {
      auto* p = args.query("threads");
      if(p != nullptr)
      {
        for(std::size_t i(0); i < p->second.size(); ++i)
        {
          std::size_t t(0u);
          if(!p->second.at(i).parse(t))
          {
            std::cout << "ERROR: Failed to parse '" << p->second.at(i) << "' as thread count" << std::endl;
            return;
          }
          num_threads.push_back(t);
        }
      }
    }

    Random rng;

    // create an empty atlas and a root mesh node
    Geometry::MeshAtlas<Mesh_> atlas;
    std::deque<std::shared_ptr<Geometry::RootMeshNode<Mesh_>>> nodes;
    nodes.push_back(std::make_shared<Geometry::RootMeshNode<Mesh_>>(nullptr, &atlas));

    // try to parse the mesh file
#ifndef DEBUG
    try
#endif
    {
      std::cout << "Parsing mesh files..." << std::endl;
      // Now parse the mesh file
      mesh_reader.parse(*nodes.back(), atlas, nullptr);
    }
#ifndef DEBUG
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return;
    }
    catch(...)
    {
      std::cerr << "ERROR: unknown exception" << std::endl;
      return;
    }
#endif

    // refine
    std::cout << "Refining up to level " << lvl_max << "..." << std::endl;
    for(Index lvl(1); lvl <= lvl_max; ++lvl)
    {
      nodes.push_back(nodes.back()->refine_shared());
    }

    static constexpr std::size_t nperms(5u);

    std::deque<std::array<BenchResults, nperms>> res_poisson(nodes.size());
    std::deque<std::array<BenchResults, nperms>> res_burgers(nodes.size());

    std::deque<std::deque<std::array<BenchResults, nperms>>> res_poisson_smp(nodes.size());
    std::deque<std::deque<std::array<BenchResults, nperms>>> res_burgers_smp(nodes.size());

    const bool b_poisson = (args.check("no-poisson") < 0);
    const bool b_burgers = (args.check("no-burgers") < 0);

    std::cout << std::endl << "Performing benchmark..." << std::endl;
    for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
    {
      std::cout << "Level " << stringify(lvl).pad_front(2) << ": ";

      // allocate smp deques
      res_poisson_smp.at(lvl).resize(num_threads.size());
      res_burgers_smp.at(lvl).resize(num_threads.size());

      for(std::size_t perm(0); perm < nperms; ++perm)
      {
        std::cout << '.';
        std::cout.flush();

        // clone node
        auto node = nodes.at(lvl)->clone_shared();

        // permute
        TimeStamp stamp_1;
        switch(perm)
        {
        case 0:
          break;

        case 1:
          node->create_permutation(Geometry::PermutationStrategy::random);
          break;

        case 2:
          node->create_permutation(Geometry::PermutationStrategy::lexicographic);
          break;

        case 3:
          node->create_permutation(Geometry::PermutationStrategy::colored);
          break;

        case 4:
          /*continue; // skip this one
          node->create_permutation(Geometry::PermutationStrategy::cuthill_mckee);
          break;

        case 5:
          continue; // skip this one
          node->create_permutation(Geometry::PermutationStrategy::cuthill_mckee_reversed);
          break;

        case 6:*/
          node->create_permutation(Geometry::PermutationStrategy::geometric_cuthill_mckee);
          break;

        /*case 7:
          continue; // skip this one
          node->create_permutation(Geometry::PermutationStrategy::geometric_cuthill_mckee_reversed);
          break;*/
        }
        TimeStamp stamp_2;
        res_poisson.at(lvl).at(perm).permute_time =
        res_burgers.at(lvl).at(perm).permute_time = stamp_2.elapsed(stamp_1);

        if(b_poisson)
        {
          bench_poisson(*node->get_mesh(), res_poisson.at(lvl).at(perm));
          for(std::size_t i(0); i < num_threads.size(); ++i)
            bench_poisson_new(*node->get_mesh(), res_poisson_smp.at(lvl).at(i).at(perm), num_threads.at(i));
        }
        if(b_burgers)
        {
          bench_burgers(*node->get_mesh(), res_burgers.at(lvl).at(perm));
          for(std::size_t i(0); i < num_threads.size(); ++i)
            bench_burgers_new(*node->get_mesh(), res_burgers_smp.at(lvl).at(i).at(perm), num_threads.at(i));
        }
      }
      std::cout << " done!"<< std::endl;
    }

    for(std::size_t i(0); i < num_threads.size(); ++i)
    {
      std::cout << std::endl << "New Assembly chosen thread counts with " << num_threads[i] << " worker threads:" << std::endl;
      //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;

      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ": ";

        for(std::size_t perm(0); perm < nperms; ++perm)
        {
          std::cout << stringify(res_poisson_smp.at(lvl).at(i).at(perm).num_worker_threads).pad_front(12);
        }
        std::cout << std::endl;
      }
    }

    std::cout << std::endl << "Permutation Timings:" << std::endl;
    //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
    std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;

    for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
    {
      std::cout << stringify(lvl).pad_front(2) << ": ";
      for(std::size_t perm(0); perm < nperms; ++perm)
      {
        std::cout << stringify_fp_fix(res_poisson.at(lvl).at(perm).permute_time, 6, 12);
      }
      std::cout << std::endl;
    }

    if(b_poisson)
    {
      std::cout << std::endl << "Poisson Assembly Timings:" << std::endl;
      //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;

      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ": ";

        for(std::size_t perm(0); perm < nperms; ++perm)
        {
          const auto& b = res_poisson.at(lvl).at(perm);
          double t = b.num_asm_time / double(b.num_asm_count);
          std::cout << stringify_fp_fix(t, 6, 12);
        }
        std::cout << std::endl;
      }

      /*std::cout << std::endl << "Poisson Assembly Timing Relative to 2-Level:" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;

      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ": ";
        const double t2 = res_poisson.at(lvl).front().num_asm_time / double(res_poisson.at(lvl).front().num_asm_count);

        for(std::size_t perm(0); perm < nperms; ++perm)
        {
          const auto& b = res_poisson.at(lvl).at(perm);
          double t = b.num_asm_time / double(b.num_asm_count);
          std::cout << stringify_fp_fix(t/t2, 6, 12);
        }
        std::cout << std::endl;
      }*/


      for(std::size_t i(0); i < num_threads.size(); ++i)
      {
        std::cout << std::endl << "New Poisson Assembly Timings with " << num_threads[i] << " worker threads:" << std::endl;
        //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
        std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;

        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ": ";

          for(std::size_t perm(0); perm < nperms; ++perm)
          {
            const auto& b = res_poisson_smp.at(lvl).at(i).at(perm);
            double t = b.num_asm_time / double(b.num_asm_count);
            std::cout << stringify_fp_fix(t, 6, 12);
          }
          std::cout << std::endl;
        }
      }

      for(std::size_t i(0); i < num_threads.size(); ++i)
      {
        std::cout << std::endl << "New Poisson Assembly Timings with " << num_threads[i] << " worker threads: explicit assembly() Time" << std::endl;
      //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;
        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ": ";
          for(std::size_t perm(0); perm < nperms; ++perm)
            std::cout << stringify_fp_fix(res_poisson_smp.at(lvl).at(i).at(perm).thread_time_assemble, 6, 12);
          std::cout << std::endl;
        }
      }
      for(std::size_t i(0); i < num_threads.size(); ++i)
      {
        std::cout << std::endl << "New Poisson Assembly Timings with " << num_threads[i] << " worker threads: explicit wait() Time" << std::endl;
      //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;
        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ": ";
          for(std::size_t perm(0); perm < nperms; ++perm)
            std::cout << stringify_fp_fix(res_poisson_smp.at(lvl).at(i).at(perm).thread_time_wait, 6, 12);
          std::cout << std::endl;
        }
      }

      /*std::cout << std::endl << "New Poisson Assembly Timing Relative to 2-Level:" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;

      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ": ";
        const double t2 = res_poisson2.at(lvl).front().num_asm_time / double(res_poisson2.at(lvl).front().num_asm_count);

        for(std::size_t perm(0); perm < nperms; ++perm)
        {
          const auto& b = res_poisson2.at(lvl).at(perm);
          double t = b.num_asm_time / double(b.num_asm_count);
          std::cout << stringify_fp_fix(t/t2, 6, 12);
        }
        std::cout << std::endl;
      }*/
    }

    if(b_burgers)
    {
      std::cout << std::endl << "Burgers Assembly Timings:" << std::endl;
      //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;

      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ": ";

        for(std::size_t perm(0); perm < nperms; ++perm)
        {
          const auto& b = res_burgers.at(lvl).at(perm);
          double t = b.num_asm_time / double(b.num_asm_count);
          std::cout << stringify_fp_fix(t, 6, 12);
        }
        std::cout << std::endl;
      }

      /*std::cout << std::endl << "Burgers Assembly Timing Relative to 2-Level:" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;

      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ": ";
        const double t2 = res_burgers.at(lvl).front().num_asm_time / double(res_burgers.at(lvl).front().num_asm_count);

        for(std::size_t perm(0); perm < nperms; ++perm)
        {
          const auto& b = res_burgers.at(lvl).at(perm);
          double t = b.num_asm_time / double(b.num_asm_count);
          std::cout << stringify_fp_fix(t/t2, 6, 12);
        }
        std::cout << std::endl;
      }*/

      for(std::size_t i(0); i < num_threads.size(); ++i)
      {
        std::cout << std::endl << "New Burgers Assembly Timings with " << num_threads[i] << " worker threads:" << std::endl;
        //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
        std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;

        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ": ";

          for(std::size_t perm(0); perm < nperms; ++perm)
          {
            const auto& b = res_burgers_smp.at(lvl).at(i).at(perm);
            double t = b.num_asm_time / double(b.num_asm_count);
            std::cout << stringify_fp_fix(t, 6, 12);
          }
          std::cout << std::endl;
        }
      }

      for(std::size_t i(0); i < num_threads.size(); ++i)
      {
        std::cout << std::endl << "New Burgers Assembly Timings with " << num_threads[i] << " worker threads: explicit assembly() Time" << std::endl;
        //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
        std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;
        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ": ";
          for(std::size_t perm(0); perm < nperms; ++perm)
            std::cout << stringify_fp_fix(res_burgers_smp.at(lvl).at(i).at(perm).thread_time_assemble, 6, 12);
          std::cout << std::endl;
        }
      }
      for(std::size_t i(0); i < num_threads.size(); ++i)
      {
        std::cout << std::endl << "New Burgers Assembly Timings with " << num_threads[i] << " worker threads: explicit wait() Time" << std::endl;
        //std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;
        std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     GCMK" << std::endl;
        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ": ";
          for(std::size_t perm(0); perm < nperms; ++perm)
            std::cout << stringify_fp_fix(res_burgers_smp.at(lvl).at(i).at(perm).thread_time_wait, 6, 12);
          std::cout << std::endl;
        }
      }

      /*std::cout << std::endl << "New Burgers Assembly Timing Relative to 2-Level:" << std::endl;
      std::cout << "LVL     2LEVEL      RANDOM      LEXICO      COLORED     ACMK        ACMKR       GCMK        GCMKR" << std::endl;

      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ": ";
        const double t2 = res_burgers2.at(lvl).front().num_asm_time / double(res_burgers2.at(lvl).front().num_asm_count);

        for(std::size_t perm(0); perm < nperms; ++perm)
        {
          const auto& b = res_burgers2.at(lvl).at(perm);
          double t = b.num_asm_time / double(b.num_asm_count);
          std::cout << stringify_fp_fix(t/t2, 6, 12);
        }
        std::cout << std::endl;
      }*/
    }
  }

  void main(int argc, char** argv)
  {
    // This is the list of all supported meshes that could appear in the mesh file
    typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
    typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
    typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
    typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

    SimpleArgParser args(argc, argv);

    args.support("mesh");
    args.support("level");
    args.support("threads");
    args.support("no-poisson");
    args.support("no-burgers");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if( !unsupported.empty() )
    {
      // print all unsupported options to cerr
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
      return;
    }

    int num_mesh_files = args.check("mesh");
    if(num_mesh_files < 1)
    {
      std::cerr << "ERROR: You have to specify at least one meshfile with --mesh <files...>" << std::endl;
      return;
    }

    // get our filename deque
    auto mpars = args.query("mesh");
    XASSERT(mpars != nullptr);
    const std::deque<String>& filenames = mpars->second;

    // create an empty mesh file reader
    Geometry::MeshFileReader mesh_reader;
    mesh_reader.add_mesh_files(filenames);

    // read root markup
    mesh_reader.read_root_markup();

    // get mesh type
    const String mtype = mesh_reader.get_meshtype_string();

    std::cout << "Mesh Type: " << mtype << std::endl;

    if(mtype == "conformal:hypercube:2:2") run<H2M2D>(args, mesh_reader); else
    if(mtype == "conformal:hypercube:3:3") run<H3M3D>(args, mesh_reader); else
    if(mtype == "conformal:simplex:2:2") run<S2M2D>(args, mesh_reader); else
    if(mtype == "conformal:simplex:3:3") run<S3M3D>(args, mesh_reader); else
    std::cout << "ERROR: unsupported mesh type!" << std::endl;
  }
} // namespace MeshPermAssemblyBench

int main(int argc, char** argv)
{
  FEAT::Runtime::initialize(argc, argv);
  MeshPermAssemblyBench::main(argc, argv);
  return FEAT::Runtime::finalize();
}
