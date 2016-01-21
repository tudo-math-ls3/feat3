#pragma once
#ifndef FEAST_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP
#define FEAST_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/nonlinear_functional.hpp>
#include <kernel/global/vector.hpp>
#include <control/domain/domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_precond_factory.hpp>
#include <control/meshopt/meshopt_solver_factory.hpp>

namespace FEAST
{
  namespace Control
  {
    namespace Meshopt
    {
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
          typedef typename FEAST::Meshopt::Intern::TrafoFE<Trafo_>::Space TrafoSpace;

          /// Our bas class
          typedef MeshoptControlBase<DomainControl_, Trafo_> BaseClass;

          /// Type of the "system matrix" for the solver
          template<typename A, typename B, typename C>
          using LocalQualityFunctionalType = Hyperelasticity_<A, B, C, Trafo_>;

          /// Inter level transfer matrix
          typedef TransferMatrixBlocked<Mem_, DT_, IT_, MeshType::world_dim> TransferMatrixType;

          typedef MeshoptSystemLevel<Mem_, DT_, IT_, LocalQualityFunctionalType, Global::NonlinearFunctional> SystemLevelType;
          typedef MeshoptTransferLevel<SystemLevelType, TransferMatrixType> TransferLevelType;

          typedef typename DomainControl_::LayerType DomainLayerType;
          typedef typename DomainControl_::LevelType DomainLevelType;
          typedef MeshoptAssemblerLevel<TrafoSpace> AssemblerLevelType;

          typedef Solver::NLOptPrecond
          <
            typename SystemLevelType::GlobalSystemVector,
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
          std::shared_ptr<Solver::IterativeSolver<typename SystemLevelType::GlobalSystemVector>> solver;
          /// The preconditioner. As this might involve a matrix to be assembled, we keep it between solves.
          std::shared_ptr<PrecondType> precond;

          /**
           * Constructor
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
                _assembler_levels.at(i)->assemble_coords_buffer(*_system_levels.at(i));
                // Assemble the system filter, all homogeneous
                _assembler_levels.at(i)->assemble_system_filter(*_system_levels.at(i));
              }

              for(Index i(0); (i+1) < num_levels; ++i)
              {
                _assembler_levels.at(i+1)->assemble_system_transfer(*_transfer_levels.at(i), *_assembler_levels.at(i));
              }

              precond = MeshoptPrecondFactory::create_nlopt_precond(*this, dom_ctrl);

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
          }

          /// \copydoc BaseClass::name()
          virtual String name() const override
          {
            return "HyperelasticityFunctionalControl<>";
          }

          /// \copydoc BaseClass::get_coords()
          virtual typename SystemLevelType::GlobalCoordsBuffer& get_coords() override
          {
            return _system_levels.back()->coords_buffer;
          }

          /// \copydoc BaseClass::buffer_to_mesh()
          virtual void buffer_to_mesh() override
          {
            const typename SystemLevelType::GlobalCoordsBuffer& coords_buffer(_system_levels.back()->coords_buffer);
            const auto& coords_buffer_loc = *coords_buffer;

            for(size_t level(0); level < _assembler_levels.size(); ++level)
            {
              auto& vertex_set = _assembler_levels.at(level)->mesh.get_vertex_set();
              for(Index i(0); i < vertex_set.get_num_vertices(); ++i)
                vertex_set[i] = coords_buffer_loc(i);
            }
          }

          /// \copydoc BaseClass::mesh_to_buffer()
          virtual void mesh_to_buffer() override
          {
            typename SystemLevelType::GlobalCoordsBuffer& coords_buffer(_system_levels.back()->coords_buffer);
            auto& coords_buffer_loc = *coords_buffer;

            const auto& vertex_set = _assembler_levels.back()->mesh.get_vertex_set();
            for(Index i(0); i < vertex_set.get_num_vertices(); ++i)
              coords_buffer_loc(i, vertex_set[i]);
          }

          /// \copydoc BaseClass::add_to_vtk_exporter()
          virtual void add_to_vtk_exporter(
            Geometry::ExportVTK<MeshType>& exporter, const int deque_position) const override
          {
            const auto& sys_lvl = this->_system_levels.at(size_t(deque_position));

            (*(sys_lvl->op_sys)).add_to_vtk_exporter(exporter);

            for(auto& it:(*(sys_lvl->filter_sys)).template at<0>())
            {
              String field_name("nu_"+it.first);
              exporter.add_vertex_vector(field_name, it.second.get_nu());
            }
          }

          /// \copydoc BaseClass::prepare()
          virtual void prepare(const typename SystemLevelType::GlobalSystemVector& vec_state) override
          {

            for(size_t level(num_levels); level > 0; )
            {
              --level;
              Index ndofs(_assembler_levels.at(level)->trafo_space.get_num_dofs());

              LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, MeshType::world_dim> vec_buf;
              vec_buf.convert(*vec_state);

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::GlobalCoordsBuffer
                global_vec_level( &(_system_levels.at(level)->gate_sys), vec_buf, ndofs, Index(0));

              (*(_system_levels.at(level)->op_sys)).prepare
                ((*global_vec_level), *(_system_levels.at(level)->filter_sys));

              (*(_system_levels.at(level)->op_sys)).init();
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
            typename SystemLevelType::GlobalSystemVector vec_rhs = the_asm_level.assemble_rhs_vector(the_system_level);
            typename SystemLevelType::GlobalSystemVector vec_sol = the_asm_level.assemble_sol_vector(the_system_level);

            Statistics::reset_flops();
            Statistics::reset_times();
            Statistics::reset_solver_statistics();

            TimeStamp at;

            // solve
            //Solver::solve(*solver, vec_sol, vec_rhs, the_system_level.op_sys, the_system_level.filter_sys);
            the_system_level.op_sys.reset_num_evals();
            //if(precond != nullptr)
            //  precond->init_numeric();

            Solver::Status st = solver->correct(vec_sol, vec_rhs);
            TimeStamp bt;

            auto filtered_grad = the_system_level.op_sys.create_vector_r();
            auto unfiltered_grad = the_system_level.op_sys.create_vector_r();
            the_system_level.op_sys.compute_grad(filtered_grad);
            the_system_level.op_sys.compute_grad(unfiltered_grad);

            the_system_level.filter_sys.filter_def(filtered_grad);

            std::cout << " |grad_f| = " << stringify_fp_sci(filtered_grad.norm2()) <<
              "|grad| = " << stringify_fp_sci(unfiltered_grad.norm2()) << std::endl;
            unfiltered_grad.axpy(filtered_grad, unfiltered_grad, DataType(-1));

            std::cout << " |grad_f - grad| = " << stringify_fp_sci(unfiltered_grad.norm2()) << std::endl;

            // Print solver summary
            if(Comm::rank() == 0)
            {
              std::cout << solver->get_plot_name() << ": " << st << ", " << solver->get_num_iter();
              std::cout << " its, defect initial/final: " << stringify_fp_sci(solver->get_def_initial());
              std::cout << " / " << stringify_fp_sci(solver->get_def_final()) << std::endl;
              std::cout << "Needed evaluations: " << the_system_level.op_sys.get_num_func_evals() << " (func) / " << the_system_level.op_sys.get_num_grad_evals();
              std::cout <<  " (grad) / " << the_system_level.op_sys.get_num_hess_evals() << " (hess)" << std::endl;
            }

            (*the_system_level.coords_buffer).convert(*vec_sol);

            buffer_to_mesh();
          }

      }; // class HyperelasticityFunctionalControl
    } // namespace Meshopt
  } // namespace Control
} // namespace FEAST

#endif // FEAST_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP
