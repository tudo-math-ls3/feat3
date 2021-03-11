// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isosphere/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/richardson.hpp>

namespace DbgSphereFem
{
  using namespace FEAT;

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;
  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;
  typedef LAFEM::MeanFilter<MemType, DataType, IndexType> FilterType;
  typedef LAFEM::Transfer<MatrixType> TransferType;

  typedef std::shared_ptr<Solver::SolverBase<VectorType>> SolverTypePtr;

  typedef Shape::Triangle ShapeType;
  //typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType, 3> MeshType;

  //typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Trafo::IsoSphere::Mapping<MeshType> TrafoType;

  //typedef Space::Lagrange1::Element<TrafoType> SpaceType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceType;

  class Level
  {
  public:
    MeshType mesh;
    TrafoType trafo;
    SpaceType space;

    MatrixType matrix;
    FilterType filter;
    TransferType transfer;

    SolverTypePtr smoother;

  public:
    explicit Level(Geometry::Factory<MeshType>& factory) :
      mesh(factory),
      trafo(mesh),
      space(trafo)
    {
      auto& vtx = mesh.get_vertex_set();
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        vtx[i].normalize();
      }
    }

    std::shared_ptr<Level> refine() const
    {
      Geometry::StandardRefinery<MeshType> refinery(mesh);
      return std::make_shared<Level>(refinery);
    }
  };

  void main(int /*argc*/, char** /*argv*/)
  {
    Cubature::DynamicFactory cubature("auto-degree:" + stringify(2*SpaceType::local_degree+1));
    Assembly::Common::LaplaceBeltramiOperator oper;
    Analytic::Common::SphereSinBubbleFunction sol_func;
    Assembly::Common::LaplaceFunctional<decltype(sol_func)> force(sol_func);

    std::vector<DataType> vcr, vh0, vh1;
    std::vector<Index> vnit;

    const Index lvl_max = 5;

    // levels
    std::deque<std::shared_ptr<Level>> levels;

    // create levels
    for(Index ilvl(0); ilvl <= lvl_max; ++ilvl)
    {
      if(ilvl == 0)
      {
        std::cout << "Creating: 0";
        Geometry::RefinedUnitSphereFactory<MeshType> factory(0);
        levels.push_back(std::make_shared<Level>(factory));
      }
      else
      {
        std::cout << " " << stringify(ilvl);
        levels.push_front(levels.front()->refine());
      }

      // get the current level
      Level& lvl = *levels.front();

      // assemble matrix
      Assembly::SymbolicAssembler::assemble_matrix_std1(lvl.matrix, lvl.space);
      lvl.matrix.format();
      Assembly::BilinearOperatorAssembler::assemble_matrix1(lvl.matrix, oper, lvl.space, cubature);

      // assemble mean filter
      Assembly::MeanFilterAssembler::assemble(lvl.filter, lvl.space, cubature);

      if(ilvl == 0)
      {
        // create coarse grid solver
        //lvl.smoother = Solver::new_umfpack_mean(lvl.matrix, lvl.filter);
        lvl.smoother = Solver::new_pcg(lvl.matrix, lvl.filter);
        continue;
      }

      // assemble prolongation + restriction
      {
        MatrixType& prol = lvl.transfer.get_mat_prol();
        MatrixType& rest = lvl.transfer.get_mat_rest();
        prol.format();
        Assembly::SymbolicAssembler::assemble_matrix_2lvl(prol, lvl.space, levels.at(1u)->space);
        Assembly::GridTransfer::assemble_prolongation_direct(prol, lvl.space, levels.at(1u)->space, cubature);
        rest = prol.transpose();
      }

      // create smoother
      auto jacobi = Solver::new_jacobi_precond(lvl.matrix, lvl.filter);
      auto smoother = Solver::new_richardson(lvl.matrix, lvl.filter, DataType(0.5), jacobi);
      smoother->set_max_iter(4);
      smoother->set_min_iter(4);
      lvl.smoother = smoother;
    }
    std::cout << std::endl;

    // multigrid hierarchy
    auto mg_hierarchy = std::make_shared<Solver::MultiGridHierarchy<MatrixType, FilterType, TransferType>>(lvl_max+1);

    // create levels
    for(Index ilvl(0); ilvl <= lvl_max; ++ilvl)
    {
      // get the current level
      Level& lvl = *levels.at(ilvl);

      if(ilvl == lvl_max)
      {
        mg_hierarchy->push_level(lvl.matrix, lvl.filter, lvl.smoother);
      }
      else
      {
        mg_hierarchy->push_level(lvl.matrix, lvl.filter, lvl.transfer,
          lvl.smoother, lvl.smoother, lvl.smoother);
      }
    }

    mg_hierarchy->init();

    // Assemble and solve
    for(Index ilvl(1); ilvl <= lvl_max; ++ilvl)
    {
      std::cout << String(80, '*') << std::endl;
      std::cout << "Level: " << stringify(ilvl).pad_front(2) << std::endl;

      // get the current level
      Level& lvl = *levels.at(lvl_max-ilvl);
      std::cout << "Dofs.: " << stringify(lvl.matrix.rows()) << std::endl;

      // create vectors
      VectorType vec_sol = lvl.matrix.create_vector_l();
      VectorType vec_rhs = lvl.matrix.create_vector_l();

      vec_sol.format();
      vec_rhs.format();

      // assemble rHS
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, force, lvl.space, cubature);

      // create solver
      auto multigrid = Solver::new_multigrid(mg_hierarchy, Solver::MultiGridCycle::V, int(lvl_max-ilvl), -1);
      //multigrid->set_adapt_cgc(Solver::MultiGridAdaptCGC::MinEnergy);
      auto solver = Solver::new_richardson(lvl.matrix, lvl.filter, DataType(1), multigrid);
      solver->set_max_iter(1000);
      solver->set_plot_mode(Solver::PlotMode::summary);
      solver->init();
      solver->apply(vec_sol, vec_rhs);
      solver->done();

      vnit.push_back(solver->get_num_iter());
      vcr.push_back(solver->calc_convergence_rate());

      // try to compute error
      auto err = Assembly::ScalarErrorComputer<1, true>::compute(vec_sol, sol_func, lvl.space, cubature);
      std::cout << "H0-Error: " << stringify_fp_sci(err.norm_h0) << std::endl;
      std::cout << "H1-Error: " << stringify_fp_sci(err.norm_h1) << std::endl;

      vh0.push_back(err.norm_h0);
      vh1.push_back(err.norm_h1);

      for(std::size_t k(0); k < vh0.size(); ++k)
      {
        std::cout << "> " << stringify(k+1).pad_front(2);
        std::cout << ": " << stringify(vnit[k]).pad_front(3);
        std::cout << "  " << stringify_fp_fix(vcr[k], 3, 5);
        std::cout << " | " << stringify_fp_sci(vh0[k]) << "   " << stringify_fp_sci(vh1[k]);
        if(k > std::size_t(0))
        {
          std::cout << " | " << stringify_fp_fix(vh0[k-1]/vh0[k]) << "   " << stringify_fp_fix(vh1[k-1]/vh1[k]);
        }
        std::cout << std::endl;
      }
    }
    mg_hierarchy->done();
  }
} // namespace DbgSphereFem

int main(int argc, char** argv)
{
  FEAT::Runtime::initialize(argc, argv);
  DbgSphereFem::main(argc, argv);
  return FEAT::Runtime::finalize();
}
