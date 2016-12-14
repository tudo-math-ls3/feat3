#pragma once
#ifndef FEAT_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP
#define FEAT_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/global/matrix.hpp>
#include <kernel/global/nonlinear_functional.hpp>
#include <kernel/global/vector.hpp>

#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>

#include <control/domain/domain_control.hpp>
#include <control/solver_factory.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_precond_factory.hpp>

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
          typedef NonlinearSystemLevel
          <
            Mem_, DT_, IT_,
            LocalQualityFunctionalType
          > SystemLevelType;

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

        public:
          /// Solver configuration
          PropertyMap& solver_config;
          /// Name of the solver configuration from solver_config we want
          const String solver_name;
          /// The solver
          //std::shared_ptr<Solver::NLOptLS
          //<
          //  typename SystemLevelType::GlobalFunctional,
          //  typename SystemLevelType::GlobalSystemFilter>
          //> solver;
          std::shared_ptr<Solver::IterativeSolver<typename SystemLevelType::GlobalSystemVectorR>> solver;
          /// The preconditioner. As this might involve a matrix to be assembled, we keep it between solves.
          std::shared_ptr<PrecondType> precond;
          /// The level index (= number of refinements since the mesh file) to optimise the mesh on
          int meshopt_lvl;
          /// The position of this level in the deque of system levels
          size_t meshopt_lvl_pos;

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
            const int meshopt_lvl_,
            const std::deque<String>& dirichlet_list,
            const std::deque<String>& slip_list,
            const String& solver_name_,
            PropertyMap& solver_config_,
            Args_&&... args):
            BaseClass(dom_ctrl),
            _assembler_levels(),
            _system_levels(),
            solver_config(solver_config_),
            solver_name(solver_name_),
            precond(nullptr),
            meshopt_lvl(meshopt_lvl_),
            meshopt_lvl_pos(~size_t(1))
            {
              XASSERT(meshopt_lvl >= -1);

              // If the input level was set to -1, take the max level of the domain control
              if(meshopt_lvl == -1)
              {
                meshopt_lvl = dom_ctrl.get_levels().back()->get_level_index();
              }

              XASSERT(meshopt_lvl <= dom_ctrl.get_levels().back()->get_level_index());

              // Now find the position of the mesh optimisation level in the domain levels
              for(size_t i(0); i < dom_ctrl.get_levels().size(); ++i)
              {
                if(dom_ctrl.get_levels().at(i)->get_level_index() == meshopt_lvl)
                {
                  meshopt_lvl_pos  = i;
                }
              }

              XASSERT(meshopt_lvl_pos < dom_ctrl.get_levels().size());

              const DomainLayerType& layer = *dom_ctrl.get_layers().back();
              const std::deque<DomainLevelType*>& domain_levels = dom_ctrl.get_levels();

              for(Index i(0); i < domain_levels.size(); ++i)
              {
                // Push new assembler level first
                _assembler_levels.push_back(
                  new AssemblerLevelType(*domain_levels.at(i), dirichlet_list, slip_list));
                // Push new system level, this needs some references to members of the assembler level
                _system_levels.push_back(new SystemLevelType(
                  domain_levels.at(i)->get_level_index(),
                  dirichlet_list, slip_list,
                  domain_levels.at(i)->get_mesh_node(),
                  _assembler_levels.at(i)->trafo_space,
                  _assembler_levels.at(i)->dirichlet_asm,
                  _assembler_levels.at(i)->slip_asm,
                  std::forward<Args_>(args)...));
              }

              for(Index i(0); i < get_num_levels(); ++i)
              {
                _assembler_levels.at(i)->assemble_gates(layer, *_system_levels.at(i));
                // Assemble the system filter, all homogeneous
                _assembler_levels.at(i)->assemble_system_filter(*_system_levels.at(i));
                // Call the operator's init() on the current level. This needs the gates, so it cannot be called
                // earlier
                _system_levels.at(i)->global_functional.init();
              }

              for(Index i(0); i < get_num_levels(); ++i)
              {
                _assembler_levels.at(i)->assemble_system_transfer(*_system_levels.at(i), *_assembler_levels.at(i));
              }

              auto* solver_section = solver_config.query_section(solver_name);
              if(solver_section == nullptr)
              {
                throw InternalError(__func__,__FILE__,__LINE__,"Could not find section for solver "+solver_name);
              }

              precond = MeshoptPrecondFactory::create_nlopt_precond(*this, dom_ctrl, solver_section);

              solver = Control::SolverFactory::create_nonlinear_optimiser
                (_system_levels.at(meshopt_lvl_pos)->global_functional,
                _system_levels.at(meshopt_lvl_pos)->filter_sys, &solver_config, solver_name, precond);
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

            // Finish the solver
            solver->done();

            // Finish the preconditioner
            if(precond != nullptr)
            {
              precond->done();
            }
          }

          /// \copydoc BaseClass::compute_cell_size_defect()
          virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max,
              CoordType& vol_min, CoordType& vol_max, CoordType& vol) const override
          {

            _system_levels.back()->global_functional.local().compute_cell_size_defect_pre_sync(vol_min, vol_max, vol);

            vol_min = _system_levels.back()->gate_sys.min(vol_min);
            vol_max = _system_levels.back()->gate_sys.max(vol_max);
            vol = _system_levels.back()->gate_sys.sum(vol);

            CoordType cell_size_defect = _system_levels.back()->global_functional.local().compute_cell_size_defect_post_sync(
                lambda_min, lambda_max, vol_min, vol_max, vol);

            lambda_min = _system_levels.back()->gate_sys.min(lambda_min);
            lambda_max = _system_levels.back()->gate_sys.max(lambda_max);
            cell_size_defect = _system_levels.back()->gate_sys.sum(cell_size_defect);

            return cell_size_defect;
          }

          /// \copydoc BaseClass::name()
          virtual String name() const override
          {
            return "HyperelasticityFunctionalControl<>";
          }

          /// \copydoc BaseClass::get_num_levels()
          virtual size_t get_num_levels() const override
          {
            return _system_levels.size();
          }

          /// \copydoc BaseClass::print()
          virtual void print() const override
          {
            Index pad_width(30);
            Dist::Comm comm_world(Dist::Comm::world());

            String msg;

            msg = name().pad_back(pad_width, '.') + String(":");
            comm_world.print(msg);

            msg = String("level max / min").pad_back(pad_width, '.') + String(": ")
              + stringify(_assembler_levels.back()->domain_level.get_level_index()) + String(" / ")
              + stringify(_assembler_levels.front()->domain_level.get_level_index());
            comm_world.print(msg);

            msg = String("optimisation on level").pad_back(pad_width, '.') + String(": ")
              + stringify(meshopt_lvl);
            comm_world.print(msg);
            msg = String("optimisation on level").pad_back(pad_width, '.') + String(": ")
              + stringify(meshopt_lvl);
            comm_world.print(msg);

            for(const auto& it : get_dirichlet_boundaries())
            {
              msg = String("Displacement BC on").pad_back(pad_width, '.') + String(": ") + it;
              comm_world.print(msg);
            }

            for(const auto& it : get_slip_boundaries())
            {
              msg = String("Unilateral BC of place on").pad_back(pad_width, '.') + String(": ") + it;
              comm_world.print(msg);
            }

            msg = String("DoF").pad_back(pad_width, '.') + String(": ")
              + stringify(_system_levels.back()->global_functional.columns());
            comm_world.print(msg);

            FEAT::Statistics::expression_target = name();
            try
            {
              msg = String("Solver") .pad_back(pad_width, '.') + String(": ")
                + FEAT::Statistics::get_formatted_solver_tree().trim();
              comm_world.print(msg);
            }
            catch(std::exception& /*e*/)
            {
            }

            if(precond != nullptr)
            {
              msg = String("Nonlinear preconditioner:");
              comm_world.print(msg);
              precond->print();
            }

            (*(_system_levels.back()->global_functional)).print();
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
            _system_levels.back()->global_functional.local().buffer_to_mesh();

            // Get the coords buffer on the finest level
            const auto& coords_buffer_loc = _system_levels.back()->coords_buffer.local();

            // Transfer fine coords buffer to coarser levels and perform buffer_to_mesh
            for(size_t level(get_num_levels()-1); level > 0; )
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

              _system_levels.at(level)->global_functional.local().get_coords().copy(vec_level);
              _system_levels.at(level)->global_functional.local().buffer_to_mesh();
            }
          }

          /// \copydoc BaseClass::mesh_to_buffer()
          virtual void mesh_to_buffer() override
          {
            // Write finest level
            (*(_system_levels.back()->global_functional)).mesh_to_buffer();

            // Get the coords buffer on the finest level
            const auto& coords_buffer_loc = _system_levels.back()->coords_buffer.local();

            // Transfer fine coords buffer to coarser levels and perform buffer_to_mesh
            for(size_t level(get_num_levels()-1); level > 0; )
            {
              --level;
              Index ndofs(_assembler_levels.at(level)->trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::LocalCoordsBuffer vec_level(coords_buffer_loc, ndofs, Index(0));

              _system_levels.at(level)->global_functional.local().get_coords().copy(vec_level);
            }

          }

          /// \copydoc BaseClass::get_dirichlet_boundaries()
          virtual std::deque<String> get_dirichlet_boundaries() const override
          {
            std::deque<String> dirichlet_boundaries;

            for(const auto& it:_assembler_levels.back()->dirichlet_asm)
            {
              dirichlet_boundaries.push_back(it.first);
            }

            return dirichlet_boundaries;
          }

          /// \copydoc BaseClass::get_slip_boundaries()
          virtual std::deque<String> get_slip_boundaries() const override
          {
            std::deque<String> slip_boundaries;

            for(const auto& it:_assembler_levels.back()->slip_asm)
            {
              slip_boundaries.push_back(it.first);
            }

            return slip_boundaries;
          }

          /// \copydoc BaseClass::add_to_vtk_exporter()
          virtual void add_to_vtk_exporter(
            Geometry::ExportVTK<MeshType>& exporter, const int lvl_index) const override
          {

            for(size_t pos(0); pos < get_num_levels(); ++pos)
            {
              if(this->_system_levels.at(pos)->get_level_index() == lvl_index)
              {

                const auto& sys_lvl = this->_system_levels.at(pos);
                (*(sys_lvl->global_functional)).add_to_vtk_exporter(exporter);

                DataType fval(0);
                auto grad = sys_lvl->global_functional.create_vector_r();

                sys_lvl->global_functional.eval_fval_grad(fval, grad);
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

            for(size_t level(get_num_levels()); level > 0; )
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

              _system_levels.at(level)->global_functional.prepare
                (global_vec_level, _system_levels.at(level)->filter_sys);

              // Call init() here to recompute all scales with current_** ScaleComputation
              _system_levels.at(level)->global_functional.init();
            }

          }

          /// \copydoc BaseClass::optimise()
          virtual void optimise() override
          {
            // fetch our finest levels
            //DomainLevelType& the_domain_level = *domain_levels.back();
            SystemLevelType& the_system_level = *_system_levels.at(meshopt_lvl_pos);
            AssemblerLevelType& the_asm_level = *_assembler_levels.at(meshopt_lvl_pos);

            // create our RHS and SOL vectors
            typename SystemLevelType::GlobalSystemVectorR vec_rhs(the_asm_level.assemble_rhs_vector(the_system_level));
            typename SystemLevelType::GlobalSystemVectorL vec_sol(the_asm_level.assemble_sol_vector(the_system_level));

            TimeStamp at;

            // solve
            //Solver::solve(*solver, vec_sol, vec_rhs, the_system_level.global_functional, the_system_level.filter_sys);

            // Let it be knownst to Statistics that it was Us who called the solver
            FEAT::Statistics::expression_target = name();

            the_system_level.global_functional.reset_num_evals();
            Solver::Status st = solver->correct(vec_sol, vec_rhs);
            TimeStamp bt;

            solver->plot_summary(st);

            // If the mesh was not optimised on the finest domain level, we now need to prolongate the solution by:
            //  - refining the coarse vertex set using the StandardRefinery
            //  - copying the results to the CoordsBuffer
            //  - convert the CoordsBuffer to a vector type that we can filter with the system's Dirichlet filter
            //  - copy the filtered vector's contents back to the buffer and write the buffer to the mesh
            //  - apply the nonlinear filter representing unilateral BCs of place by calling adapt() for the
            //    corresponding meshparts
            for(size_t pos(meshopt_lvl_pos+1); pos < get_num_levels(); ++pos)
            {
              auto& coarse_mesh = this->_dom_ctrl.get_levels().at(pos-1)->get_mesh();
              auto& fine_mesh = this->_dom_ctrl.get_levels().at(pos)->get_mesh();
              auto& fine_vtx = fine_mesh.get_vertex_set();

              // Refine coarse vertex set and write the result to the CoordsBuffer
              Geometry::StandardRefinery<MeshType> refinery(coarse_mesh);
              refinery.fill_vertex_set(fine_vtx);
              _system_levels.at(pos)->global_functional.local().mesh_to_buffer();
              // Convert the buffer to a filterable vector
              typename SystemLevelType::GlobalSystemVectorL::LocalVectorType vec_sol_lvl;
              vec_sol_lvl.convert(*(_system_levels.at(pos)->coords_buffer));
              // Filter this vector, copy back the contents and write the changes to the mesh
              auto& dirichlet_filters_lvl = (*(_system_levels.at(pos)->filter_sys)).template at<1>();
              dirichlet_filters_lvl.filter_sol(vec_sol_lvl);
              _system_levels.at(pos)->coords_buffer.local().copy(vec_sol_lvl);
              _system_levels.at(pos)->global_functional.local().buffer_to_mesh();
              // Now call adapt() on the slip boundaries
              auto* fine_mesh_node = this->_dom_ctrl.get_levels().at(pos)->get_mesh_node();
              for(const auto& it:get_slip_boundaries())
              {
                fine_mesh_node->adapt_by_name(it);
              }
            }

            // Now we need to update all levels below the one we carried out the optimisation on
            for(size_t pos(meshopt_lvl_pos); pos > size_t(0); )
            {
              --pos;
              Index ndofs(_assembler_levels.at(pos)->trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::GlobalCoordsBuffer
                global_sol_level( &(_system_levels.at(pos)->gate_sys), vec_sol.local(), ndofs, Index(0));

              _system_levels.at(pos)->global_functional.prepare(global_sol_level, _system_levels.at(pos)->filter_sys);

              _assembler_levels.at(pos)->assemble_system_filter(*(_system_levels.at(pos)), global_sol_level);
            }
          }

      }; // class HyperelasticityFunctionalControl
    } // namespace Meshopt
  } // namespace Control
} // namespace FEAT

#endif // FEAT_CONTROL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_CONTROL_HPP
