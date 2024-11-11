// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/string.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/jacobi_precond.hpp>

/// \compilerhack ICC insists on warning #2259, even when enclosed in push/pop statements
/// The compiler does not actually issue the warning until the end of the file.
/// At that point the warning stack has been popped by the #pragma warning (pop).
/// So unfortunately this specific warning cannot be disabled by putting push/pop around the declaration.
#if defined(FEAT_COMPILER_INTEL) && defined(FEAT_HAVE_FLOATX)
_Pragma("warning(disable:2259)")
#endif

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;


#if defined(FEAT_HAVE_FLOATX) && !defined(FEAT_HAVE_CUDA)
  // IEEE 754-2008 binary16 aka half-precision
  // must use double as backend type because float is not supported
  typedef flx::floatx<5, 10, double> flx_f16;  // half-prec
  typedef flx::floatx<8, 7, double> flx_bf16;  // bfloat 16
  //typedef flx::floatx<8, 23, double> flx_f16;  // single-prec
#endif


// We're opening a new namespace for our tutorial.
namespace MultiPrecHierarchBench
{
  // Once again, we use quadrilaterals.
  //typedef Shape::Quadrilateral ShapeType;
  typedef Shape::Hypercube<2> ShapeType;
  // Use the unstructured conformal mesh class
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  // Define the corresponding mesh-part type
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  // Use the standard transformation mapping
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  // Use the Lagrange-1 element
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  // Our data arrays should be double precision.
  //typedef double DataType;
  // Use the default index type for indexing.
  typedef Index IndexType;

  double choose_tol(double)
  {
    return 1E-14;
  }

  float choose_tol(float)
  {
    return 1E-6F;
  }

#ifdef FEAT_HAVE_QUADMATH
  __float128 choose_tol(__float128)
  {
    return 1E-31Q;
  }
#endif

#if defined(FEAT_HAVE_FLOATX) && !defined(FEAT_HAVE_CUDA)
  flx_f16 choose_tol(flx_f16)
  {
    return flx_f16(1E-2);
  }
  flx_bf16 choose_tol(flx_bf16)
  {
    return flx_bf16(5E-2);
  }
#endif

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  class Level
  {
  public:
    // The mesh for this level
    MeshType mesh;
    // The trafo for this level
    TrafoType trafo;
    // The space for this level
    SpaceType space;

    explicit Level(Geometry::Factory<MeshType>& mesh_factory) :
      mesh(mesh_factory),
      trafo(mesh),
      space(trafo)
    {
    }
  }; // class Level

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  template<typename DT_, typename IT_>
  LAFEM::SparseMatrixCSR<DT_, IT_> expand_prol(const LAFEM::SparseMatrixCSR<DT_, IT_>& prol,
    const Index n, const DT_ alpha = DT_(1))
  {
    XASSERT(prol.rows() <= n);

    // expand prolongation matrix
    const Index pnr = prol.rows();
    const Index pnc = prol.columns();
    const Index pnz = prol.used_elements();
    const IT_* row_ptr_p = prol.row_ptr();
    const IT_* col_idx_p = prol.col_ind();
    const DT_* val_p = prol.val();

    // allocate expanded matrix
    LAFEM::SparseMatrixCSR<DT_, IT_> ex_prol(n, n, pnz + n - pnc);
    IT_* row_ptr_q = ex_prol.row_ptr();
    IT_* col_idx_q = ex_prol.col_ind();
    DT_* val_q = ex_prol.val();

    // copy coarse identity
    for(Index i(0); i < pnc; ++i)
    {
      row_ptr_q[i] = row_ptr_p[i];
      for(Index j(row_ptr_p[i]); j < row_ptr_p[i+1]; ++j)
      {
        col_idx_q[j] = col_idx_p[j];
        val_q[j] = alpha * val_p[j];
      }
    }
    row_ptr_q[pnc] = row_ptr_p[pnc];

    // copy interpolation + append identity
    for(Index i(pnc); i < pnr; ++i)
    {
      Index jq = row_ptr_q[i];
      for(Index j(row_ptr_p[i]); j < row_ptr_p[i+1]; ++j, ++jq)
      {
        col_idx_q[jq] = col_idx_p[j];
        val_q[jq] = alpha * val_p[j];
      }
      col_idx_q[jq] = i;
      val_q[jq] = DT_(1);
      row_ptr_q[i+1] = ++jq;
    }

    // append identity
    for(Index i(pnr); i < n; ++i)
    {
      Index jq = row_ptr_q[i];
      col_idx_q[jq] = i;
      val_q[jq] = DT_(1);
      row_ptr_q[i+1] = ++jq;
    }

    // export matrix
    return ex_prol;
  }

  template<typename DT_, typename IT_>
  void reduce_system(
    LAFEM::SparseMatrixCSR<DT_, IT_>& matrix,
    LAFEM::DenseVector<DT_, IT_>& vector,
    const LAFEM::SparseMatrixCSR<DT_, IT_>& prol,
    const DT_ alpha = DT_(1))
  {
    const Index n = matrix.rows();

    // expand prolongation matrix
    LAFEM::SparseMatrixCSR<DT_, IT_> mprol = expand_prol(prol, matrix.rows(), alpha);
    LAFEM::SparseMatrixCSR<DT_, IT_> mrest = mprol.transpose();

    // (R*A*P) * (P^-1x) = R*b

    // compute matrix structure
    Adjacency::Graph gap(Adjacency::RenderType::injectify_sorted, matrix, mprol);
    Adjacency::Graph grap(Adjacency::RenderType::injectify_sorted, mrest, gap);

    // compute matrix
    LAFEM::SparseMatrixCSR<DT_, IT_> matrix_r(grap);
    matrix_r.format();
    matrix_r.add_double_mat_product(mrest, matrix, mprol);

    // shrink matrix
    {
      const DT_* val = matrix_r.val();
      DT_ max_val = DT_(0);
      for(Index i(0); i < matrix_r.used_elements(); ++i)
        max_val = Math::max(max_val, Math::abs(val[i]));
      matrix_r.shrink(max_val * Math::pow(Math::eps<DT_>(), DT_(0.7)));
    }

    // compute rhs
    LAFEM::DenseVector<DT_, IT_> vector_r(n);
    mrest.apply(vector_r, vector);

    // overwrite input
    matrix = std::move(matrix_r);
    vector = std::move(vector_r);
  }

  template<typename DT_, typename IT_, typename Mesh_>
  void solve_by_thomas(
    LAFEM::DenseVector<DT_, IT_>& vec_sol,
    const LAFEM::DenseVector<DT_, IT_>& vec_rhs,
    const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix,
    const Mesh_& mesh)
  {
    const auto& vtx = mesh.get_vertex_set();
    const Index n = matrix.rows();
    XASSERT(n == vtx.get_num_vertices());

    std::map<double, Index> vtx_map;
    for(Index i(0); i < n; ++i)
      vtx_map.emplace(vtx[i][0], i);

    std::vector<Index> idx_map(vtx_map.size()), inv_map(vtx_map.size());
    {
      auto it = vtx_map.begin();
      for(Index i(0); i < n; ++i, ++it)
      {
        idx_map.at(i) = it->second;
        inv_map.at(it->second) = i;
      }
    }

    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    const DT_* mat_val = matrix.val();
    const DT_* rhs_val = vec_rhs.elements();
    DT_* sol_val = vec_sol.elements();

    std::vector<DT_> va(n), vb(n), vc(n), vd(n), vx(n);

    for(Index i(0); i < n; ++i)
    {
      const Index ii = inv_map[i];
      vd[i] = rhs_val[idx_map[i]];
      for(IT_ j(row_ptr[i]); j < row_ptr[i+1]; ++j)
      {
        const Index jj = inv_map[col_idx[j]];
        XASSERT((ii <= jj+1) && (jj <= ii+1));
        if(jj == ii)
          vb.at(ii) = mat_val[j];
        else if(jj < ii)
          va.at(ii) = mat_val[j];
        else
          vc.at(ii) = mat_val[j];
      }
    }

    // apply filter to matrix and vector
    vd.front() = vd.back() = DT_(0);
    vb.front() = DT_(1);
    vc.front() = DT_(0);
    va.back() = DT_(0);
    vb.back() = DT_(1);

    // forward sweep
    for(Index i(1); i < n; ++i)
    {
      DT_ wi = va[i] / vb[i-1];
      vb[i] -= wi*vc[i-1];
      vd[i] -= wi*vd[i-1];
    }

    // backward sweep
    vx.back() = vd.back() / vb.back();
    for(Index i(n-1); i > 0; )
    {
      --i;
      vx[i] = (vd[i] - vc[i]*vx[i+1]) / vb[i];
    }

    // unsort solution
    for(Index i(0); i < n; ++i)
      sol_val[idx_map[i]] = vx[i];
  }


  // Here's our tutorial's main function
  template<typename DataType, typename DTI_>
  void main(const Index level_max, const Index level_min, bool reduce = false)
  {
    typedef LAFEM::DenseVector<DataType, IndexType> VectorType;
    typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixType;
    typedef LAFEM::UnitFilter<DataType, IndexType> FilterType;

    typedef LAFEM::DenseVector<DTI_, IndexType> VectorTypeI;
    typedef LAFEM::SparseMatrixCSR<DTI_, IndexType> MatrixTypeI;
    typedef LAFEM::UnitFilter<DTI_, IndexType> FilterTypeI;

    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    //Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;
    //Analytic::Common::Q2BubbleFunction<ShapeType::dimension> sol_function;
    Analytic::Common::ExpBubbleFunction<ShapeType::dimension> sol_function;
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);
    Assembly::Common::LaplaceOperator laplace_operator;


    // First of all, let's create the coarse level
    std::shared_ptr<Level> level;
    {
      Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(1);
      level = std::make_shared<Level>(mesh_factory);
    }

    // vector for prolongation matrices
    std::vector<MatrixType> prol_mats;
    prol_mats.reserve(level_max);

    // print header line
    std::cout << "\n";
    std::cout << "Lvl      Norm2(RHS)";
    if(reduce)
      std::cout << "       #DOFs     #NZE(A)  #NZE(H[A])";
    if((MeshType::shape_dim > 1) || reduce)
      std::cout << "   Iter";
    std::cout << "   Abs-Def-Norm   Rel-Def-Norm       H0-Error       H1-Error     Time";
    std::cout << "\n";

    // Now, let's refine up to the maximum level
    for(Index ilevel(1); ilevel <= level_max; ++ilevel)
    {
      TimeStamp stamp;
      // refine mesh
      auto level_c = level;
      {
        Geometry::StandardRefinery<MeshType> mesh_factory(level_c->mesh);
        level = std::make_shared<Level>(mesh_factory);
      }

      // assemble prolongation matrix
      {
        MatrixType mat_prol;
        Assembly::SymbolicAssembler::assemble_matrix_2lvl(mat_prol, level->space, level_c->space);
        mat_prol.format();
        Assembly::GridTransfer::assemble_prolongation_direct(mat_prol, level->space, level_c->space, cubature_factory);
        mat_prol.shrink(Math::pow(Math::eps<DataType>(), DataType(0.7)));
        prol_mats.push_back(mat_prol.clone());
      }

      // skip if below level min
      if(ilevel < level_min)
        continue;

      std::cout << stringify(ilevel).pad_front(3) << ":";

      // Assemble the Laplace matrix:
      MatrixType matrix;
      {
        Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, level->space);
        matrix.format();
        Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, laplace_operator, level->space, cubature_factory);
      }

      // Assemble the unit filter
      FilterType filter;
      {
        Geometry::BoundaryFactory<MeshType> boundary_factory(level->mesh);
        MeshPartType boundary(boundary_factory);
        Assembly::UnitFilterAssembler<MeshType> unit_asm;
        unit_asm.add_mesh_part(boundary);
        unit_asm.assemble(filter, level->space);
      }

      // create vectorscl
      VectorType vec_sol = matrix.create_vector_r();
      VectorType vec_rhs = matrix.create_vector_r();
      vec_sol.format();
      vec_rhs.format();

      // And assemble the rhs vector:
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, force_functional, level->space, cubature_factory);

      // filter the vectors
      filter.filter_sol(vec_sol);
      filter.filter_rhs(vec_rhs);

      const DataType rhs_norm = vec_rhs.norm2();
      std::cout << stringify_fp_sci(rhs_norm).pad_front(15);

      // reduced system
      MatrixType matrix_r = matrix.clone(LAFEM::CloneMode::Deep);
      VectorType vec_sol_r = vec_sol.clone(LAFEM::CloneMode::Deep);
      VectorType vec_rhs_r = vec_rhs.clone(LAFEM::CloneMode::Deep);
      if(reduce)
      {
        // recude systems
        for(auto it = prol_mats.rbegin(); it != prol_mats.rend(); ++it)
        {
          reduce_system(matrix_r, vec_rhs_r, *it);
          filter.filter_rhs(vec_rhs_r);
        }
        std::cout << stringify(matrix.rows()).pad_front(12);
        std::cout << stringify(matrix.used_elements()).pad_front(12);
        std::cout << stringify(matrix_r.used_elements()).pad_front(12);
      }

      // SOLVER
      // convert system
      MatrixTypeI matrix_i;
      FilterTypeI filter_i;
      VectorTypeI vec_rhs_i, vec_sol_i;

      matrix_i.convert(matrix_r);
      filter_i.convert(filter);
      vec_rhs_i.convert(vec_rhs_r);
      vec_sol_i.convert(vec_sol_r);

      if((MeshType::shape_dim == 1) && !reduce)
      {
        // solve by thomas
        solve_by_thomas(vec_sol_i, vec_rhs_i, matrix_i, level->mesh);
      }
      else
      {
        auto jacobi = Solver::new_jacobi_precond(matrix_i, filter_i);
        auto solver = Solver::new_pcg(matrix_i, filter_i, jacobi);

        //solver->set_plot_mode(Solver::PlotMode::iter);
        solver->set_max_iter(100000);
        //solver->set_tol_rel(DTI_(0.001) * Math::eps<DTI_>());
        //solver->set_tol_rel(Math::pow(Math::eps<DTI_>(), DTI_(0.7)));
        //solver->set_tol_rel(DTI_(10.0) * Math::eps<DTI_>());
        solver->set_tol_rel(choose_tol(DTI_(0)));
        solver->init();
        solver->apply(vec_sol_i, vec_rhs_i);
        solver->done();
        std::cout << stringify(solver->get_num_iter()).pad_front(7);
      }
      vec_sol_r.convert(vec_sol_i);
      // END OF SOLVER

      if(reduce)
      {
        // unreduce solution
        for(auto it = prol_mats.begin(); it != prol_mats.end(); ++it)
        {
          auto mprol = expand_prol(*it, vec_sol_r.size());
          auto vtmp = vec_sol_r.clone();
          mprol.apply(vec_sol_r, vtmp);
          filter.filter_sol(vec_sol_r);
        }
      }

      vec_sol.copy(vec_sol_r);

      // compute final defect
      VectorType vec_def = matrix.create_vector_r();
      matrix.apply(vec_def, vec_sol, vec_rhs, -DataType(1));
      filter.filter_def(vec_def);
      const DataType def_norm = vec_def.norm2();

      // Compute and print errors
      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
        vec_sol, sol_function, level->space, cubature_factory);

      std::cout
        << stringify_fp_sci(def_norm).pad_front(15)
        << stringify_fp_sci(def_norm/rhs_norm).pad_front(15)
        << stringify_fp_sci(errors.norm_h0).pad_front(15)
        << stringify_fp_sci(errors.norm_h1).pad_front(15)
        << stamp.elapsed_string_now().pad_front(9)
        << "\n";
    } // end of level loop
  } // void main(...)
} // namespace MultiPrecHierarchBench

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialize our runtime environment
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  Index level_max(5), level_min(1);
  bool reduce(false);

  String prec_asm("dp"), prec_sol("dp");

  // First of all, let's see if we have command line parameters.
  if(argc > 1)
    prec_asm = argv[1];
  if(argc > 2)
    prec_sol = argv[2];
  if(argc > 3)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level->
    int ilevel(0);
    if(!String(argv[3]).parse(ilevel) || (ilevel < 1))
    {
      // failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[3] << "' as maximum refinement level->" << "\n";
      std::cerr << "Note: The first argument must be a positive integer." << "\n";
      Runtime::abort();
    }
    // parse successful
    level_max = Index(ilevel);
  }
  if(argc > 4)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level->
    int ilevel(0);
    if(!String(argv[4]).parse(ilevel) || (ilevel < 1))
    {
      // failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[4] << "' as minimum refinement level->" << "\n";
      std::cerr << "Note: The second argument must be a positive integer." << "\n";
      Runtime::abort();
    }
    // parse successful
    level_min = Index(ilevel);
  }
  if(argc > 5)
  {
    reduce = true;
  }

  // Perform level sanity checks
  if(level_max < level_min)
  {
    // minimum level greater than maximum level
    std::cerr << "ERROR: Minimum level " << level_min << " is greater than the maximum level " << level_max << "\n";
    Runtime::abort();
  }
  if(level_max == level_min)
  {
    // minimum and maximum levels are equal
    std::cout << "WARNING: Minimum and maximum levels are equal" << "\n";
  }

  String precs = prec_asm + ":" + prec_sol;

  // print selected levels
  std::cout << "Precision: " << precs << "\n";
  std::cout << "Level-Min: " << level_min << "\n";
  std::cout << "Level-Max: " << level_max << "\n";
  std::cout << "Reduce   : " << (reduce ? "yes" : "no") << "\n";

  // call the tutorial's main function
#ifdef FEAT_HAVE_QUADMATH
  if(precs == "qp:qp") MultiPrecHierarchBench::main<__float128, __float128>(level_max, level_min, reduce); else
  if(precs == "qp:dp") MultiPrecHierarchBench::main<__float128, double>(level_max, level_min, reduce); else
  if(precs == "qp:sp") MultiPrecHierarchBench::main<__float128, float>(level_max, level_min, reduce); else
#endif
  if(precs == "dp:dp") MultiPrecHierarchBench::main<double, double>(level_max, level_min, reduce); else
  if(precs == "dp:sp") MultiPrecHierarchBench::main<double, float>(level_max, level_min, reduce); else
  if(precs == "sp:sp") MultiPrecHierarchBench::main<float, float>(level_max, level_min, reduce); else
#if defined(FEAT_HAVE_FLOATX) && !defined(FEAT_HAVE_CUDA)
  if(precs == "dp:hp") MultiPrecHierarchBench::main<double, flx_f16>(level_max, level_min, reduce); else
  if(precs == "sp:hp") MultiPrecHierarchBench::main<float, flx_f16>(level_max, level_min, reduce); else
  if(precs == "dp:bp") MultiPrecHierarchBench::main<double, flx_bf16>(level_max, level_min, reduce); else
  if(precs == "sp:bp") MultiPrecHierarchBench::main<float, flx_bf16>(level_max, level_min, reduce); else
#endif
  {
    std::cout << "ERROR: unsupported precision combo " << precs << "\n";
    Runtime::abort();
  }

  MemoryUsage mem_use;
  std::cout << "\n" << mem_use.get_formatted_memory_usage() << "\n";

  // Finalize our runtime environment
  return 0;
}
