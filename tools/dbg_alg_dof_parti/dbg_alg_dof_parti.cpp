
// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime
#include <kernel/util/dist.hpp>                            // NEW: for Dist::Comm

// FEAT-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart
#include <kernel/geometry/mesh_node.hpp>                   // NEW: for RootMeshNode, MeshNode
#include <kernel/geometry/unit_cube_patch_generator.hpp>   // NEW: for UnitCubePatchGenerator

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>                      // for SineBubbleFunction
#include <kernel/adjacency/export_tga.hpp>


// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/error_computer.hpp>              // for L2/H1-error computation
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for LaplaceFunctional
#include <kernel/assembly/mirror_assembler.hpp>            // NEW: for MirrorAssembler

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/none_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/vector_mirror.hpp>                  // NEW: for VectorMirror

// FEAT-Global includes
#include <kernel/global/gate.hpp>                          // NEW: for Global::Gate
#include <kernel/global/filter.hpp>                        // NEW: for Global::Filter
#include <kernel/global/matrix.hpp>                        // NEW: for Global::Matrix
#include <kernel/global/vector.hpp>                        // NEW: for Global::Vector
#include <kernel/global/alg_dof_parti.hpp>

// FEAT-Solver includes
#include <kernel/solver/pcg.hpp>                           // for PCG
#include <kernel/solver/schwarz_precond.hpp>               // NEW: for SchwarzPrecond
#include <kernel/solver/ilu_precond.hpp>                   // NEW: for ILUPrecond

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial06
{
  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;
  typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> LocalMatrixType;
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> LocalVectorType;
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> LocalFilterType;
  typedef LAFEM::VectorMirror<MemType, DataType, IndexType> VectorMirrorType;
  typedef Global::Gate<LocalVectorType, VectorMirrorType> GateType;
  typedef Global::Matrix<LocalMatrixType, VectorMirrorType, VectorMirrorType> GlobalMatrixType;
  typedef Global::Vector<LocalVectorType, VectorMirrorType> GlobalVectorType;
  typedef Global::Filter<LocalFilterType, VectorMirrorType> GlobalFilterType;

  typedef Global::AlgDofPartiVector<LocalVectorType, VectorMirrorType> AlgDofPartiVectorType;
  typedef Global::AlgDofPartiMatrix<LocalMatrixType, VectorMirrorType> AlgDofPartiMatrixType;


  void write_tga(const String& filename, const AlgDofPartiMatrixType& adp_mat)
  {
    const Dist::Comm& comm = *adp_mat.get_comm();

    // lexicographic permutation for level 2
    static const Index p_25[25] =
    {
      0,4,1,11,9,
      6,8,7,14,13,
      2,5,3,12,10,
      18,20,19,24,23,
      15,17,16,22,21
    };
    // lexicographic permutation for level 3
    static const Index p_81[81] =
    {
      /*0,4,1,11,9,27,25,33,31,
      6,8,7,14,13,30,29,36,35,
      2,5,3,12,10,28,26,34,32,
      18,20,19,24,23,40,39,44,43,
      15,17,16,22,21,38,37,42,41,
      48,50,49,54,53,68,67,72,71,
      45,47,46,52,51,66,65,70,69,
      58,60,59,64,63,76,75,80,79,
      55,57,56,62,61,74,73,78,77*/
      0,4,1,11,9,17,15,23,21,
      6,8,7,14,13,20,19,26,25,
      2,5,3,12,10,18,16,24,22,
      30,32,31,36,35,40,39,44,43,
      27,29,28,34,33,38,37,42,41,
      48,50,49,54,53,58,57,62,61,
      45,47,46,52,51,56,55,60,59,
      66,68,67,72,71,76,75,80,79,
      63,65,64,70,69,74,73,78,77
    };

    static constexpr Index h_max = 800;
    static constexpr Index n_max = 100;
    Index q[n_max], v[n_max];
    unsigned int r[n_max], s[n_max], scan[h_max];
    for(Index i(0); i < n_max; ++i)
      v[i] = 0;

    const Index n = adp_mat.get_alg_dof_parti()->get_num_global_dofs();

    switch(n)
    {
    case 25:
      for(Index i(0); i < n; ++i)
        q[p_25[i]] = i;
      break;

    case 81:
      for(Index i(0); i < n; ++i)
        q[p_81[i]] = i;
      break;

    default:
      return;
    }

    const Index block_size = h_max / n;

    const LocalMatrixType& mat = adp_mat.owned();
    const IndexType* row_ptr = mat.row_ptr();
    const IndexType* col_idx = mat.col_ind();

    Adjacency::Graph graph(n, n, mat.used_elements());
    Index* dom_ptr = graph.get_domain_ptr();
    Index* img_idx = graph.get_image_idx();

    const Index glob_dof_off = adp_mat.get_alg_dof_parti()->get_global_dof_offset();

    for(Index i(0); i < mat.rows(); ++i)
    {
      v[q[glob_dof_off+i]] = (row_ptr[i+1]-row_ptr[i]);
    }

    dom_ptr[0] = 0;
    for(Index i(0); i < n; ++i)
    {
      dom_ptr[i+1] = dom_ptr[i] + v[i];
      v[i] = dom_ptr[i];
    }

    for(Index i(0); i < mat.rows(); ++i)
    {
      Index row = q[glob_dof_off + i];
      for(Index j(row_ptr[i]); j < row_ptr[i + 1]; ++j)
      {
        img_idx[v[row]++] = q[col_idx[j]];
      }
    }

    FILE* fo = nullptr;
    if(comm.rank() == 0)
    {
      fo = fopen(filename.c_str(), "wb");
      XASSERT(fo != nullptr);

      typedef unsigned char u8;
      // set up header
      unsigned char header[18];
      for(int i(0); i < 18; ++i)
        header[i] = 0;

      // get dimensions
      const Index w = block_size*n;
      const Index h = block_size*n;

      // set dimensions
      header[12] = u8( w       & 0xFF);
      header[13] = u8((w >> 8) & 0x7F);
      header[14] = u8( h       & 0xFF);
      header[15] = u8((h >> 8) & 0x7F);

      // set basic stuff
      header[ 2] = u8( 2); // datatype code
      header[16] = u8(32); // bits per pixel
      header[17] = u8( 8); // image descriptor

      fwrite(header, 1, 18, fo);
    }

    // this ranks color
    unsigned int mask = 0;
    switch(comm.rank() & 3)
    {
    case 0: mask = 0xFF7F7F7F; break;
    case 1: mask = 0xFFFF0000; break;
    case 2: mask = 0xFF00FF00; break;
    case 3: mask = 0xFF0000FF; break;
    }

    for(Index i(0); i < n; ++i)
    {
      // reset row
      for(Index j(0); j < n; ++j)
        s[j] = 0;

      // set indices
      for(Index j(dom_ptr[n-i-1]); j < dom_ptr[n-i]; ++j)
        s[img_idx[j]] = mask;

      // send to rank 0
      for(int k(1); k < comm.size(); ++k)
      {
        if(comm.rank() == 0)
        {
          comm.recv(r, n, k);
          for(Index j(0); j < n; ++j)
            s[j] |= r[j];
        }
        else if(comm.rank() == k)
          comm.send(s, n, 0);
      }

      if(comm.rank() == 0)
      {
        // write white line
        for(Index k = 0; k < block_size*n; ++k)
          scan[k] = 0xFF404040;
        fwrite(scan, 4, block_size*n, fo);

        // expand to block size
        for(Index j(0); j < n; ++j)
        {
          scan[j*block_size] = 0xFF404040;
          for(Index k(1); k < block_size; ++k)
            scan[j*block_size+k] = s[j];
        }
        for(Index k(1); k < block_size; ++k)
          fwrite(scan, 4, block_size*n, fo);
      }
    }

    // write footer
    if(fo != nullptr)
    {
      const char* footer = "\0\0\0\0\0\0\0\0TRUEVISION-XFILE.\0";
      fwrite(footer, 1, 26, fo);
      fclose(fo);
    }
  }

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function
  void main(Index level, int dbg_rank = -1)
  {
    Dist::Comm comm = Dist::Comm::world();
    comm.print("\nNumber of processes: " + stringify(comm.size()));
    if((comm.size() != 1) && (comm.size() != 4) && (comm.size() != 16) && (comm.size() != 64))
    {
      comm.print(std::cerr, "ERROR: You must run this tutorial with 1, 4, 16 or 64 processes!");
      Runtime::abort();
    }

    std::shared_ptr<RootMeshNodeType> root_mesh_node;

    std::vector<int> neighbour_ranks;

    Index lvl = Geometry::UnitCubePatchGenerator<MeshType>::create(
      comm.rank(), comm.size(), root_mesh_node, neighbour_ranks);

    String msg = "Neighbours of process " + stringify(comm.rank()) + ":";
    for(int i : neighbour_ranks)
      msg += " " + stringify(i);
    comm.allprint(msg);

    comm.print("\nBase Mesh Level: " + stringify(lvl));

    if(lvl < level)
    {
      comm.print("Refining Mesh to Level " + stringify(level) + "...");

      for(; lvl < level; ++lvl)
      {
        root_mesh_node = std::shared_ptr<RootMeshNodeType>(root_mesh_node->refine());
      }
    }

    MeshType& mesh = *root_mesh_node->get_mesh();
    TrafoType trafo(mesh);
    SpaceType space(trafo);


    GateType gate(comm);
    for(auto it = neighbour_ranks.begin(); it != neighbour_ranks.end(); ++it)
    {
      const int neighbour_rank = *it;
      const MeshPartType* neighbour_halo = root_mesh_node->find_halo_mesh_part(neighbour_rank);

      XASSERTM(neighbour_halo != nullptr, "Failed to retrieve neighbour halo!");
      VectorMirrorType neighbour_mirror;
      Assembly::MirrorAssembler::assemble_mirror(neighbour_mirror, space, *neighbour_halo);
      gate.push(neighbour_rank, std::move(neighbour_mirror));
    }

    gate.compile(LocalVectorType(space.get_num_dofs()));

    GlobalMatrixType matrix(&gate, &gate);

    LocalMatrixType& matrix_local = matrix.local();

    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_local, space);

    GlobalVectorType vec_sol = matrix.create_vector_r();
    GlobalVectorType vec_rhs = matrix.create_vector_l();
    GlobalVectorType vec_def = matrix.create_vector_l();
    GlobalVectorType vec_tmp = matrix.create_vector_l();
    GlobalVectorType vec_sol_hypre = matrix.create_vector_l();

    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    matrix_local.format();

    Assembly::Common::LaplaceOperator laplace_operator;
    //Assembly::Common::IdentityOperator laplace_operator;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_local, laplace_operator, space, cubature_factory);

    vec_rhs.local().format();

    // Again, we use the sine-bubble as a reference solution:
    Analytic::Common::ConstantFunction<ShapeType::dimension> rhs_function(1.0);
    //Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;

    //Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);
    Assembly::Common::ForceFunctional<decltype(rhs_function)> force_functional(rhs_function);

    Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs.local(), force_functional, space, cubature_factory);

    vec_rhs.sync_0();

    vec_sol.format();


    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    for(int ibnd(0); ibnd < 2*ShapeType::dimension; ++ibnd)
    {
      MeshPartType* bnd_mesh_part = root_mesh_node->find_mesh_part("bnd:" + stringify(ibnd));
      if(bnd_mesh_part != nullptr)
      {
        unit_asm.add_mesh_part(*bnd_mesh_part);
      }
    }

    GlobalFilterType filter;
    unit_asm.assemble(filter.local(), space);

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    comm.print(String(100, '#'));

    comm.print("Assembling AlgDofParti...");

    Global::AlgDofParti<LocalVectorType, VectorMirrorType> adp;

    // breakpoint for debugging
#ifdef FEAT_COMPILER_MICROSOFT
    if(comm.rank() == dbg_rank) __debugbreak();
#else
    (void)dbg_rank;
#endif

    adp.assemble_by_gate(gate);
    adp.assemble_allgather(true);

    /*
    {
      comm.print("Global AlgDofParti DOF Info:");
      String s;
      s += stringify(adp.get_global_dof_offset()).pad_front(4);
      s += " : ";
      s += stringify(adp.get_num_owned_dofs()).pad_front(4);
      s += " : ";
      s += stringify(adp.get_num_global_dofs()).pad_front(4);
      s += " : S[";
      for(const auto& x : adp._neighbours_owner)
      {
        s += " ";
        s += stringify(x.first);
        s += ":";
        s += stringify(x.second.num_indices());
      }
      s += " ] : R[";
      for(const auto& x : adp._neighbours_donee)
      {
        s += " ";
        s += stringify(x.first);
        s += ":";
        s += stringify(x.second.num_indices());
      }
      s += "] : D[";
      for(const auto& i : adp._glob_dof_idx)
      {
        s += " ";
        s += stringify(i);
      }
      s += "]";
      comm.allprint(s);
    }
    */

    /*
    comm.print("Testing...");

    // create a global vector
    Random rng(Random::get_seed());
    GlobalVectorType vector_g = vec_rhs.clone();
    GlobalVectorType vector_x = vec_rhs.clone();
    GlobalVectorType vector_e = vec_rhs.clone();
    AlgDofPartiVectorType adp_vector(&adp);

    for(int k(0); k < 5; ++k)
    {
      comm.print(String(100, '*'));
      vector_g.format(rng, -1.0, 1.0);
      vector_g.sync_1();
      DataType norm_g = vector_g.norm2();
      comm.print(String("Norm(g) = ") + stringify_fp_sci(norm_g, 12));

      // create formatted clone
      vector_x.format(std::nan(nullptr));

      // download g -> l
      adp_vector.download(vector_g);
      DataType norm_l, norm_l_l = adp_vector.local().norm2sqr();
      comm.allreduce(&norm_l_l, &norm_l, 1, Dist::op_sum);
      norm_l = Math::sqrt(norm_l);
      comm.print(String("Norm(l) = ") + stringify_fp_sci(norm_l, 12));

      // upload l -> g
      adp_vector.upload(vector_x);
      DataType norm_x = vector_x.norm2();
      comm.print(String("Norm(x) = ") + stringify_fp_sci(norm_x, 12));

      // compute error
      vector_e.axpy(vector_g, vector_x, -1.0);
      DataType norm_e = vector_e.norm2();
      comm.print(String("Norm(E) = ") + stringify_fp_sci(norm_e, 12));

      comm.allprint(String("Norm(l_l) = ") + stringify_fp_sci(norm_l_l, 12));
    }*/

    comm.print(String(100, '#'));
    comm.print("Creating AlgDofParti Matrix...");

    AlgDofPartiMatrixType adp_matrix(&adp);

    adp_matrix.upload_symbolic(matrix);

    {
      comm.print("AlgDofParti Matrix Info:");
      String s;
      s += stringify(adp_matrix.owned().rows()).pad_front(4);
      s += " x ";
      s += stringify(adp_matrix.owned().columns()).pad_front(4);
      s += " : ";
      s += stringify(adp_matrix.owned().used_elements()).pad_front(4);
      s += " vs ";
      s += stringify(matrix.local().rows()).pad_front(4);
      s += " x ";
      s += stringify(matrix.local().columns()).pad_front(4);
      s += " : ";
      s += stringify(matrix.local().used_elements()).pad_front(4);
      comm.allprint(s);
    }

    //write_tga("adp_matrix.tga", adp_matrix);


    comm.print("Initialising Matrix Values");

    adp_matrix.upload_numeric(matrix);

    adp_matrix.filter_matrix(filter);

    AlgDofPartiVectorType adp_vector_x(&adp);
    AlgDofPartiVectorType adp_vector_b(&adp);
    AlgDofPartiVectorType adp_vector_y(&adp);
    AlgDofPartiVectorType adp_vector_d(&adp);

    adp_vector_x.upload(vec_sol);
    adp_vector_b.upload(vec_rhs);

    // compute norm of rhs
    const double norm_b = vec_rhs.norm2();

    // ############################################################################################
    // ############################################################################################
    // ############################################################################################
/*
    //         "  0: 1.147408e-01 | 1.147408e-01 | 1.147408e-01"
    comm.print(String(100, '*') + "\n>>> Richardson\n");
    comm.print("ITS: Global       | AlgDofParti  | Hypre");
    comm.print("-----------------------------------------------");
    for(int step = 0; step < 20; ++step)
    {
      // compute defect vectors
      matrix.apply(vec_def, vec_sol, vec_rhs, -1.0);
      filter.filter_def(vec_def);

      //
      adp_matrix.apply(adp_vector_y, adp_vector_x);
      adp_vector_y.owned().scale(adp_vector_y.owned(), -1.0);
      adp_vector_y.owned().axpy(adp_vector_b.owned(), adp_vector_y.owned());

      // BLAS-1/2-style
      //HYPRE_Int HYPRE_ParVectorCopy( HYPRE_ParVector x , HYPRE_ParVector y );
      //HYPRE_Int HYPRE_ParCSRMatrixMatvec( HYPRE_Complex alpha , HYPRE_ParCSRMatrix A , HYPRE_ParVector x , HYPRE_Complex beta , HYPRE_ParVector y );
      HYPRE_ParVectorCopy(par_b, par_y);
      HYPRE_ParCSRMatrixMatvec(-1.0, parcsr_A, par_x, 1.0, par_y);

      // compute defect norms
      double norm_d = vec_def.norm2();
      double norm_alg_dof_parti_d = gate.norm2(adp_vector_y.owned().norm2());
      HYPRE_ParVectorInnerProd(par_y, par_y, &hy_tmp);
      double norm_hy_d = Math::sqrt(hy_tmp);

      comm.print(
        stringify(step).pad_front(3) + ": " + stringify_fp_sci(norm_d) + " | " +
        stringify_fp_sci(norm_alg_dof_parti_d) + " | " + stringify_fp_sci(norm_hy_d));

      double omega = 0.3;

      // update solution
      vec_sol.axpy(vec_def, vec_sol, omega);

      adp_vector_x.owned().axpy(adp_vector_y.owned(), adp_vector_x.owned(), omega);

      //HYPRE_Int HYPRE_ParVectorAxpy ( HYPRE_Complex alpha , HYPRE_ParVector x , HYPRE_ParVector y );
      HYPRE_ParVectorAxpy (omega, par_y, par_x);
    }
*/
    // ############################################################################################
    // ############################################################################################
    // ############################################################################################

    GlobalVectorType vec_dir = vec_sol.clone();
    AlgDofPartiVectorType adp_vector_q(&adp);

    // reset vectors
    vec_sol.format();
    adp_vector_x.owned().format();

    // compute initial defect
    vec_def.copy(vec_rhs);
    adp_vector_d.owned().copy(adp_vector_b.owned());

    // set initial direction
    vec_dir.copy(vec_def);
    adp_vector_y.owned().copy(adp_vector_d.owned());

    // compute initial gamma
    DataType gamma_g = vec_def.dot(vec_def);
    DataType gamma_l = gate.sum(adp_vector_d.owned().dot(adp_vector_d.owned()));

    //         "  0: 1.147408e-01 | 1.147408e-01 | 1.147408e-01"
    comm.print(String(100, '*') + "\n>>> PCG\n");
    comm.print("ITS: Global       | AlgDofParti  | Hypre");
    comm.print("-----------------------------------------------");
    for(int step = 1; step <= 50; ++step)
    {
      // compute matrix-vector product
      matrix.apply(vec_tmp, vec_dir);
      adp_matrix.apply(adp_vector_q, adp_vector_y);

      filter.filter_def(vec_tmp);

      // compute alpha
      DataType alpha_g = vec_tmp.dot(vec_dir);
      DataType alpha_l = gate.sum(adp_vector_q.owned().dot(adp_vector_y.owned()));

      alpha_g = gamma_g / alpha_g;
      alpha_l = gamma_l / alpha_l;

      // update solution
      vec_sol.axpy(vec_dir, vec_sol, alpha_g);
      adp_vector_x.owned().axpy(adp_vector_y.owned(), adp_vector_x.owned(), alpha_l);

      // update residual
      vec_def.axpy(vec_tmp, vec_def, -alpha_g);
      adp_vector_d.owned().axpy(adp_vector_q.owned(), adp_vector_d.owned(), -alpha_l);

      // update gamma
      DataType gamma_g2 = gamma_g;
      DataType gamma_l2 = gamma_l;

      gamma_g = vec_def.dot(vec_def);
      gamma_l = gate.sum(adp_vector_d.owned().dot(adp_vector_d.owned()));

      // compute defect norm
      DataType defnorm_g = Math::sqrt(gamma_g);
      DataType defnorm_l = Math::sqrt(gamma_l);
      DataType defnorm_h = 0.0;

      // print lines
      comm.print(
        stringify(step).pad_front(3) + ": " + stringify_fp_sci(defnorm_g) + " | " +
        stringify_fp_sci(defnorm_l) + " | " + stringify_fp_sci(defnorm_h));

      if(defnorm_g < 1E-8 * norm_b)
        break;

      // compute update direction
      vec_dir.axpy(vec_dir, vec_def, gamma_g / gamma_g2);
      adp_vector_y.owned().axpy(adp_vector_y.owned(), adp_vector_d.owned(), gamma_l / gamma_l2);
    }

    // ############################################################################################
    // ############################################################################################
    // ############################################################################################

    comm.print(String(100, '#'));

    comm.print("\nSolving using our own trusted PCG...");

    auto solver = Solver::new_pcg(matrix, filter);

    solver->set_tol_rel(1E-8);
    solver->set_plot_mode(Solver::PlotMode::iter);

    solver->init();

    vec_sol.format();
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    solver->done();


    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we want to export our solution to a (P)VTU file set.

    // In a parallel simulation, each process will write a separate VTU file, which contains the
    // data that is defined on the patch of the corresponding process. Moreover, one process
    // writes a single additional PVTU file, which can be read by ParaView to visualise
    // the whole domain that consists of all patches.

    // Build the VTK filename; we also append the number of processes to the filename:
/*
    String vtk_name = String("./dbg-adp-lvl") + stringify(level) + "-n" + stringify(comm.size());

    comm.print("Writing VTK file '" + vtk_name + ".pvtu'...");

    // Create a VTK exporter for our patch mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // Add the vertex-projection of our (local) solution and rhs vectors
    exporter.add_vertex_scalar("hypre", vec_sol_hypre.local().elements());
    exporter.add_vertex_scalar("sol", vec_sol_local.elements());
    exporter.add_vertex_scalar("rhs", vec_rhs_local.elements());

    // Finally, write the VTK files by calling the "write" function of the exporter and pass the
    // communicator as a second argument:
    exporter.write(vtk_name, comm);

    // Note: Do not forget the 'comm' argument in the call above as otherwise each process will
    // try to write to the same VTK file, resulting in garbage due to race conditions...
*/

    // That's all, folks.
    comm.print("Finished!");
  } // void main(...)
} // namespace Tutorial06

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialise the FEAT runtime environment:
  Runtime::initialise(argc, argv);

  // Specify the desired mesh refinement level, defaulted to 5.
  Index level(2);
  int dbg_rank = -1;

  // Now let's see if we have command line parameters: This tutorial supports passing
  // the refinement level as a command line parameter, to investigate the behaviour of the L2/H1
  // errors of the discrete solution.
  if(argc > 1)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[1]).parse(ilevel) || (ilevel < 1))
    {
      // Failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[1] << "' as refinement level." << std::endl;
      std::cerr << "Note: The last argument must be a positive integer." << std::endl;
      // Abort our runtime environment
      Runtime::abort();
    }
    // If parsing was successful, use the given information and notify the user
    level = Index(ilevel);
  }
  if(argc > 2)
  {
    String(argv[2]).parse(dbg_rank);
  }

  // call the tutorial's main function
  Tutorial06::main(level, dbg_rank);

  // And finally, finalise our runtime environment.
  return Runtime::finalise();
}
