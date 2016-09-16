#pragma once
#ifndef FEAT_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP
#define FEAT_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/util/mpi_cout.hpp>

#include <kernel/global/matrix.hpp>
#include <kernel/global/nonlinear_functional.hpp>
#include <kernel/global/vector.hpp>

#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>

#include <control/domain/domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_precond_factory.hpp>
#include <control/meshopt/meshopt_solver_factory.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Meshopt
    {
      /**
       * \brief Control class for HyperelasticityFunctionals
       *
       * \tparam Mem_
       * The memory architecture of the local HyperelasticityFunctional
       *
       * \tparam DT_
       * The floating point precision for the solver
       *
       * \tparam IT_
       * Index type
       *
       * \tparam DomainControl_
       * The domain control type this is based on
       *
       * \tparam Trafo_
       * The mesh's underlying transformation. At the time of writing, there is only Trafo::Standard, which means
       * P1/Q1 transformation.
       *
       * \tparam Hyperelasticity_
       * The (patch-) local HyperelasticityFunctional to use.
       *
       * \note Local Hyperelastiticy functionals are only implemented for Mem::Main,
       * \see FEAT::Meshopt::HyperElasticityFunctionalBase
       *
       * \author Jordi Paul
       *
       */
      template
      <
        typename Mem_, typename DT_, typename IT_, typename DomainControl_, typename Trafo_,
        template<typename, typename, typename, typename> class Hyperelasticity_
      >
      class HyperelasticityFunctionalControl
      : public MeshoptControlBase<DomainControl_, Trafo_>
      {
        public:
          /// Our memory architecture
          typedef Mem_ MemType;
          /// The floating point type
          typedef DT_ DataType;
          /// The index type
          typedef IT_ IndexType;
          /// The transformation we solve for
          typedef Trafo_ TrafoType;
          /// The type of the domain control
          typedef DomainControl_ DomainControlType;
          /// The underlying mesh type
          typedef typename DomainControl_::MeshType MeshType;
          /// The floating point type the mesh's coordinates use
          typedef typename MeshType::CoordType CoordType;
          /// The FE space the transformation lives in
          typedef typename FEAT::Meshopt::Intern::TrafoFE<Trafo_>::Space TrafoSpace;

          /// Our base class
          typedef MeshoptControlBase<DomainControl_, Trafo_> BaseClass;

          /// Type of the "system matrix" for the solver
          template<typename A, typename B, typename C>
          using LocalQualityFunctionalType = Hyperelasticity_<A, B, C, Trafo_>;

          /// Inter level transfer matrix
          typedef LAFEM::SparseMatrixBWrappedCSR<Mem_, DT_, IT_, MeshType::world_dim> TransferMatrixType;

          /// The system level type, holding all information about the nonlinear system of equations
          typedef MeshoptSystemLevel
          <
            Mem_, DT_, IT_,
            LocalQualityFunctionalType,
            Global::NonlinearFunctional
          > SystemLevelType;
          /// Type for holding inter level transfer information
          typedef MeshoptTransferLevel<SystemLevelType, TransferMatrixType> TransferLevelType;

          /// Domain layers
          typedef typename DomainControl_::LayerType DomainLayerType;
          /// Domain levels
          typedef typename DomainControl_::LevelType DomainLevelType;
          /// Type for assembling FE space based quantities like filters, gates etc.
          typedef MeshoptAssemblerLevel<TrafoSpace> AssemblerLevelType;

          /// Preconditioner type. Becaus the preconditioner is expensive to assemble symbolically and numerically,
          /// it is kept and recycled between nonlinear solver calls, so we have to keep it.
          typedef Solver::NLOptPrecond
          <
            typename SystemLevelType::GlobalSystemVectorR,
            typename SystemLevelType::GlobalSystemFilter
          > PrecondType;

          /// One assembler for every level
          std::deque<AssemblerLevelType*> _assembler_levels;
          /// These hold the system information for each level
          std::deque<SystemLevelType*> _system_levels;
          /// Inter level transfer information
          std::deque<TransferLevelType*> _transfer_levels;

        public:
          /// Number of refinement levels, same as in the domain control
          const Index num_levels;
          /// Solver configuration
          PropertyMap& solver_config;
          /// Name of the solver configuration from solver_config we want
          const String solver_name;
          /// The solver
          std::shared_ptr<Solver::IterativeSolver<typename SystemLevelType::GlobalSystemVectorR>> solver;
          /// The preconditioner. As this might involve a matrix to be assembled, we keep it between solves.
          std::shared_ptr<PrecondType> precond;

          /**
           * \brief Variadic template constructor
           *
           * \param[in] dom_ctrl
           * The domaincontrol holding all geometry information for all levels
           *
           * \param[in] dirichlet_list
           * List of meshpart identifiers for Dirichlet boundary conditions
           *
           * \param[in] slip_list
           * List of meshpart identifiers for slip boundary conditions
           *
           * \param[in] solver_name_
           * Name of the solver to select from the solver_config_
           *
           * \param[in] solver_config_
           * PropertyMap holding the solver configuration
           *
           */
          template<typename... Args_>
          explicit HyperelasticityFunctionalControl(
            DomainControl_& dom_ctrl,
            const std::deque<String>& dirichlet_list,
            const std::deque<String>& slip_list,
            const String& solver_name_,
            PropertyMap& solver_config_,
            Args_&&... args):
            _assembler_levels(),
            _system_levels(),
            _transfer_levels(),
            num_levels(dom_ctrl.get_levels().size()),
            solver_config(solver_config_),
            solver_name(solver_name_),
            precond(nullptr)
            {
              const DomainLayerType& layer = *dom_ctrl.get_layers().back();
              const std::deque<DomainLevelType*>& domain_levels = dom_ctrl.get_levels();

              for(Index i(0); i < num_levels; ++i)
              {
                // Push new assembler level first
                _assembler_levels.push_back(new AssemblerLevelType(*domain_levels.at(i), dirichlet_list, slip_list));
                // Push new system level, this needs some references to members of the assembler level
                _system_levels.push_back(new SystemLevelType(
                  dirichlet_list, slip_list,
                  domain_levels.at(i)->get_mesh_node(),
                  _assembler_levels.at(i)->trafo_space,
                  _assembler_levels.at(i)->dirichlet_asm,
                  _assembler_levels.at(i)->slip_asm,
                  std::forward<Args_>(args)...));

                // Call the operator's init() on the current level
                (*_system_levels.at(i)->op_sys).init();
                // If we are not on the coarsest level, create the new transfer level
                if(i > 0)
                {
                  _transfer_levels.push_back(new TransferLevelType(*_system_levels.at(i-1), *_system_levels.at(i)));
                }
              }

              for(Index i(0); i < num_levels; ++i)
              {
                _assembler_levels.at(i)->assemble_gates(layer, *_system_levels.at(i));
                // Assemble the system filter, all homogeneous
                _assembler_levels.at(i)->assemble_system_filter(*_system_levels.at(i));
              }

              for(Index i(0); (i+1) < num_levels; ++i)
              {
                _assembler_levels.at(i+1)->assemble_system_transfer(*_transfer_levels.at(i), *_assembler_levels.at(i));
              }

              //auto solver_section_p = solver_config.query_se(solver_name);
              //if(!solver_section_p.second)
              //  throw InternalError(__func__,__FILE__,__LINE__,"Could not find section for solver "+solver_name);

              auto* solver_section = solver_config.query_section(solver_name);
              if(solver_section == nullptr)
                throw InternalError(__func__,__FILE__,__LINE__,"Could not find section for solver "+solver_name);

              precond = MeshoptPrecondFactory::create_nlopt_precond(*this, dom_ctrl, solver_section);

              solver = Control::MeshoptSolverFactory::create_nonlinear_optimiser
                (_system_levels, _transfer_levels, &solver_config, solver_name, precond);
              solver->init();

            }

          /// Explicitly delete the default constructor
          HyperelasticityFunctionalControl(const HyperelasticityFunctionalControl&) = delete;

          /// \brief Virtual destructor
          virtual ~HyperelasticityFunctionalControl()
          {
            while(!_assembler_levels.empty())
            {
              delete _assembler_levels.back();
              _assembler_levels.pop_back();
            }

            while(!_system_levels.empty())
            {
              delete _system_levels.back();
              _system_levels.pop_back();
            }

            while(!_transfer_levels.empty())
            {
              delete _transfer_levels.back();
              _transfer_levels.pop_back();
            }

            solver->done();
            if(precond != nullptr)
              precond->done();
          }

          /// \copydoc BaseClass::compute_cell_size_defect()
          virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max,
              CoordType& vol_min, CoordType& vol_max) const override
          {
            return (*(_system_levels.back()->op_sys)).
              compute_cell_size_defect(lambda_min, lambda_max, vol_min, vol_max);
          }

          /// \copydoc BaseClass::name()
          virtual String name() const override
          {
            return "HyperelasticityFunctionalControl<>";
          }

          /// \copydoc BaseClass::print()
          virtual void print() const override
          {
            Util::mpi_cout(name()+" settings:\n");
            Util::mpi_cout_pad_line("Domain level min",_assembler_levels.front()->domain_level.get_level_index());
            Util::mpi_cout_pad_line("Domain level max",_assembler_levels.back()->domain_level.get_level_index());

            for(const auto& it : get_dirichlet_boundaries())
              Util::mpi_cout_pad_line("Displacement BC on",it);
            for(const auto& it : get_slip_boundaries())
              Util::mpi_cout_pad_line("Unilateral BC of place on",it);

            Util::mpi_cout_pad_line("DoF",_system_levels.back()->op_sys.columns());
            (*(_system_levels.back()->op_sys)).print();
            Util::mpi_cout("\n");

            FEAT::Statistics::expression_target = name();
            try
            {
              Util::mpi_cout_pad_line("Solver",FEAT::Statistics::get_formatted_solver_tree().trim() + "\n");
            }
            catch(std::exception& e)
            {
            }

            if(precond != nullptr)
            {
              Util::mpi_cout("Nonlinear preconditioner:\n");
              precond->print();
            }
          }

          /// \copydoc BaseClass::get_coords()
          virtual typename SystemLevelType::GlobalCoordsBuffer& get_coords() override
          {
            return _system_levels.back()->coords_buffer;
          }

          /// \copydoc BaseClass::buffer_to_mesh()
          virtual void buffer_to_mesh() override
          {
            // Write finest level
            (*(_system_levels.back()->op_sys)).buffer_to_mesh();

            // Get the coords buffer on the finest level
            const auto& coords_buffer_loc = *(_system_levels.back()->coords_buffer);

            // Transfer fine coords buffer to coarser levels and perform buffer_to_mesh
            for(size_t level(num_levels-1); level > 0; )
            {
              --level;
              Index ndofs(_assembler_levels.at(level)->trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::LocalCoordsBuffer
                vec_level(coords_buffer_loc, ndofs, Index(0));

              (*(_system_levels.at(level)->op_sys)).get_coords().copy(vec_level);
              (*(_system_levels.at(level)->op_sys)).buffer_to_mesh();
            }
          }

          /// \copydoc BaseClass::mesh_to_buffer()
          virtual void mesh_to_buffer() override
          {
            // Write finest level
            (*(_system_levels.back()->op_sys)).mesh_to_buffer();

            // Get the coords buffer on the finest level
            const auto& coords_buffer_loc = *(_system_levels.back()->coords_buffer);

            // Transfer fine coords buffer to coarser levels and perform buffer_to_mesh
            for(size_t level(num_levels-1); level > 0; )
            {
              --level;
              Index ndofs(_assembler_levels.at(level)->trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::LocalCoordsBuffer
                vec_level(coords_buffer_loc, ndofs, Index(0));

              (*(_system_levels.at(level)->op_sys)).get_coords().copy(vec_level);
            }

          }

          /// \copydoc BaseClass::get_dirichlet_boundaries()
          virtual std::deque<String> get_dirichlet_boundaries() const override
          {
            std::deque<String> dirichlet_boundaries;

            for(const auto& it:_assembler_levels.back()->dirichlet_asm)
              dirichlet_boundaries.push_back(it.first);

            return dirichlet_boundaries;
          }

          /// \copydoc BaseClass::get_slip_boundaries()
          virtual std::deque<String> get_slip_boundaries() const override
          {
            std::deque<String> slip_boundaries;

            for(const auto& it:_assembler_levels.back()->slip_asm)
              slip_boundaries.push_back(it.first);

            return slip_boundaries;
          }

          /// \copydoc BaseClass::add_to_vtk_exporter()
          virtual void add_to_vtk_exporter(
            Geometry::ExportVTK<MeshType>& exporter, const int deque_position) const override
          {
            const auto& sys_lvl = this->_system_levels.at(size_t(deque_position));
            (*(sys_lvl->op_sys)).add_to_vtk_exporter(exporter);

            auto grad = sys_lvl->op_sys.create_vector_r();
            sys_lvl->op_sys.compute_grad(grad);
            sys_lvl->filter_sys.filter_def(grad);
            exporter.add_vertex_vector("grad", *grad);

            for(auto& it:(*(sys_lvl->filter_sys)).template at<0>())
            {
              const String field_name("nu_"+it.first);
              // WARNING: This explicitly assumes that the filter vector belong to a P1/Q1 space and thus "lives"
              // in the mesh's vertices
              exporter.add_vertex_vector(field_name, it.second.get_filter_vector());
            }
          }

          /// \copydoc BaseClass::prepare()
          virtual void prepare(const typename SystemLevelType::GlobalSystemVectorR& vec_state) override
          {
            // Prepare the preconditioner before updating the status
            if(precond != nullptr)
              precond->init_numeric();

            typename SystemLevelType::LocalCoordsBuffer vec_buf;
            vec_buf.convert(*vec_state);

            for(size_t level(num_levels); level > 0; )
            {
              --level;
              Index ndofs(_assembler_levels.at(level)->trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::GlobalCoordsBuffer
                global_vec_level( &(_system_levels.at(level)->gate_sys), vec_buf, ndofs, Index(0));

              _system_levels.at(level)->op_sys.prepare
                (global_vec_level, _system_levels.at(level)->filter_sys);

              //(*(_system_levels.at(level)->op_sys)).init();
            }
          }

          /// \copydoc BaseClass::optimise()
          virtual void optimise() override
          {
            // fetch our finest levels
            //DomainLevelType& the_domain_level = *domain_levels.back();
            SystemLevelType& the_system_level = *_system_levels.back();
            AssemblerLevelType& the_asm_level = *_assembler_levels.back();

            // create our RHS and SOL vectors
            typename SystemLevelType::GlobalSystemVectorR vec_rhs(the_asm_level.assemble_rhs_vector(the_system_level));
            typename SystemLevelType::GlobalSystemVectorL vec_sol(the_asm_level.assemble_sol_vector(the_system_level));

            TimeStamp at;

            // solve
            //Solver::solve(*solver, vec_sol, vec_rhs, the_system_level.op_sys, the_system_level.filter_sys);

            // Let it be knownst to Statistics that it was Us who called the solver
            FEAT::Statistics::expression_target = name();

            the_system_level.op_sys.reset_num_evals();
            Solver::Status st = solver->correct(vec_sol, vec_rhs);
            TimeStamp bt;

            // Updates the mesh beyond the finest level etc.
            prepare(vec_sol);

            // Print solver summary
            if(Util::Comm::rank() == 0)
            {
              std::cout << solver->get_plot_name() << ": " << st << ", " << solver->get_num_iter();
              std::cout << " its, defect initial/final: " << stringify_fp_sci(solver->get_def_initial());
              std::cout << " / " << stringify_fp_sci(solver->get_def_final()) << std::endl;
              std::cout << "Needed evaluations: " << the_system_level.op_sys.get_num_func_evals() << " (func) / " << the_system_level.op_sys.get_num_grad_evals();
              std::cout <<  " (grad) / " << the_system_level.op_sys.get_num_hess_evals() << " (hess)" << std::endl;
            }

          }

      }; // class HyperelasticityFunctionalControl
    } // namespace Meshopt
  } // namespace Control
} // namespace FEAT

#endif // FEAT_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP
