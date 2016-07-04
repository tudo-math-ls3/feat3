#pragma once
#ifndef FEAT_CONTROL_MESHOPT_DUDV_FUNCTIONAL_CONTROL_HPP
#define FEAT_CONTROL_MESHOPT_DUDV_FUNCTIONAL_CONTROL_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/nonlinear_functional.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/meshopt/dudv_functional.hpp>

#include <control/domain/domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_solver_factory.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Meshopt
    {
      template<typename>
      class DuDvFunctionalAssemblerLevel;

      template<typename Mem_, typename DT_, typename IT_, typename DomainControl_, typename Trafo_>
      class DuDvFunctionalControl
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

          //typedef FEAT::Meshopt::DuDvFunctional<Mem_, DT_, IT_, Trafo_> MeshQualityFunctional;
          //typedef MeshoptSystemLevel<Mem_, DT_, IT_, MeshQualityFunctional::template MatrixTemplate, Global::Matrix>
          //  SystemLevelType;

          /// Template-alias away the Trafo so the SystemLevel can take it as a template template parameter
          template<typename A, typename B, typename C>
          using OperatorType =  FEAT::Meshopt::DuDvFunctional<A, B, C, Trafo_>;
          /// Linear system of equations on one refinement level
          typedef MeshoptSystemLevel<Mem_, DT_, IT_, OperatorType, Global::Matrix> SystemLevelType;

          /// Inter-level transfer matrix
          typedef TransferMatrixBlocked<Mem_, DT_, IT_, MeshType::world_dim> TransferMatrixType;
          /// Type to do all inter level information transfer
          typedef MeshoptTransferLevel<SystemLevelType, TransferMatrixType> TransferLevelType;

          typedef typename DomainControl_::LayerType DomainLayerType;
          typedef typename DomainControl_::LevelType DomainLevelType;
          typedef DuDvFunctionalAssemblerLevel<TrafoSpace> AssemblerLevelType;

          typedef typename SystemLevelType::GlobalSystemVectorL GlobalSystemVectorL;
          typedef typename SystemLevelType::GlobalSystemVectorR GlobalSystemVectorR;

          /// For every level of refinement we have one assembler level
          std::deque<AssemblerLevelType*> _assembler_levels;
          /// For every level of refinement, we have one system level
          std::deque<SystemLevelType*> _system_levels;
          /// Two subsequent levels can communicate through their transfer level
          std::deque<TransferLevelType*> _transfer_levels;

        public:
          /// Number of refinement levels
          const Index num_levels;
          /// Solver configuration
          PropertyMap& solver_config;
          /// The name of the section from solver_config we want to use
          String solver_name;

          std::shared_ptr<Solver::SolverBase<GlobalSystemVectorR>> solver;

          explicit DuDvFunctionalControl(
            DomainControl_& dom_ctrl, const std::deque<String>& dirichlet_list, const std::deque<String>& slip_list,
            const String& solver_name_, PropertyMap& solver_config_):
            _assembler_levels(),
            _system_levels(),
            _transfer_levels(),
            num_levels(dom_ctrl.get_levels().size()),
            solver_config(solver_config_),
            solver_name(solver_name_),
            solver(nullptr)
            {
              XASSERT(num_levels > Index(0));

              const DomainLayerType& layer = *dom_ctrl.get_layers().back();
              const std::deque<DomainLevelType*>& domain_levels = dom_ctrl.get_levels();

              for(Index i(0); i < num_levels; ++i)
              {
                _assembler_levels.push_back(new AssemblerLevelType(*domain_levels.at(i), dirichlet_list, slip_list));
                _system_levels.push_back(new SystemLevelType(dirichlet_list, slip_list,
                  domain_levels.at(i)->get_mesh_node(),
                  &(_assembler_levels.at(i)->trafo_space),
                  &(_assembler_levels.at(i)->dirichlet_asm),
                  &(_assembler_levels.at(i)->slip_asm)));

                // This assembles the system matrix symbolically
                (*(_system_levels.at(i)->op_sys)).init();

                if(i > 0)
                {
                  _transfer_levels.push_back(new TransferLevelType(*_system_levels.at(i-1), *_system_levels.at(i)));
                }
              }

              for(Index i(0); i < num_levels; ++i)
              {
                _assembler_levels.at(i)->assemble_gates(layer, *_system_levels.at(i));
                _assembler_levels.at(i)->assemble_coords_buffer(*_system_levels.at(i));
              }

              for(Index i(0); (i+1) < num_levels; ++i)
              {
                _assembler_levels.at(i+1)->assemble_system_transfer(*_transfer_levels.at(i), *_assembler_levels.at(i));
              }

              for(Index i(0); i < num_levels; ++i)
              {
                _assembler_levels.at(i)->assemble_system_matrix(*_system_levels.at(i));
              }

              // create our solver
              solver = Control::MeshoptSolverFactory::create_linear_solver
                (_system_levels, _transfer_levels, &solver_config, solver_name);
              // initialise
              solver->init();
            }

          /// Explicitly delete empty default constructor
          DuDvFunctionalControl() = delete;
          /// Explicitly delete move constructor
          DuDvFunctionalControl(DuDvFunctionalControl&&) = delete;

          /**
           * \brief Virtual destructor
           */
          virtual ~DuDvFunctionalControl()
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
            return "DuDvFunctionalControl<>";
          }

          /// \copydoc BaseClass::compute_cell_size_quality()
          virtual CoordType compute_cell_size_quality() const override
          {
            return (*(_system_levels.back()->op_sys)).compute_cell_size_quality();
          }

          /// \copydoc BaseClass::get_coords()
          virtual typename SystemLevelType::GlobalCoordsBuffer& get_coords() override
          {
            return _system_levels.back()->coords_buffer;
          }

          /// \copydoc BaseClass::buffer_to_mesh()
          virtual void buffer_to_mesh() override
          {
            // Write from control object to local mesh quality functional
            (*(_system_levels.back()->op_sys)).get_coords().clone(*(_system_levels.back()->coords_buffer));
            // Write finest level
            (*(_system_levels.back()->op_sys)).buffer_to_mesh();

            // Get the coords buffer on the finest level
            const typename SystemLevelType::GlobalCoordsBuffer& coords_buffer(_system_levels.back()->coords_buffer);
            const auto& coords_buffer_loc = *coords_buffer;

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

              (*(_system_levels.at(level)->op_sys)).get_coords().clone(vec_level, LAFEM::CloneMode::Deep);
              (*(_system_levels.at(level)->op_sys)).buffer_to_mesh();
            }

          }

          /// \copydoc BaseClass::mesh_to_buffer()
          virtual void mesh_to_buffer() override
          {
            // Write finest level
            (*(_system_levels.back()->op_sys)).mesh_to_buffer();
            // Write from local mesh quality functional to control object
            (*(_system_levels.back()->coords_buffer)).clone(
              (*(_system_levels.back()->op_sys)).get_coords(), LAFEM::CloneMode::Deep);

            // Get the coords buffer on the finest level
            const typename SystemLevelType::GlobalCoordsBuffer& coords_buffer(_system_levels.back()->coords_buffer);
            const auto& coords_buffer_loc = *coords_buffer;

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

              (*(_system_levels.at(level)->op_sys)).get_coords().clone(vec_level, LAFEM::CloneMode::Deep);
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

          virtual void init_numeric()
          {
            for(size_t lvl(0); lvl < size_t(num_levels); ++lvl)
            {
              _assembler_levels.at(lvl)->assemble_system_matrix(*(_system_levels.at(lvl)));
            }
          }

          /// \copydoc BaseClass::prepare()
          virtual void prepare(const GlobalSystemVectorR& vec_state) override
          {
            for(size_t level(num_levels); level > 0; )
            {
              --level;
              Index ndofs(_assembler_levels.at(level)->trafo_space.get_num_dofs());

              typename SystemLevelType::LocalCoordsBuffer vec_buf(ndofs, DT_(0));
              vec_buf.convert(*vec_state);

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::GlobalCoordsBuffer
                global_vec_level( &(_system_levels.at(level)->gate_sys), vec_buf, ndofs, Index(0));

              (*(_system_levels.at(level)->op_sys)).prepare(vec_buf, *(_system_levels.at(level)->filter_sys));

              _assembler_levels.at(level)->assemble_system_filter(*(_system_levels.at(level)), global_vec_level);
            }

          }

          virtual Solver::Status apply(GlobalSystemVectorR& vec_sol, const GlobalSystemVectorL& vec_rhs)
          {
            // Get our global system matrix and filter
            typename SystemLevelType::GlobalQualityFunctional& op_sys = (*_system_levels.back()).op_sys;
            typename SystemLevelType::GlobalSystemFilter& filter_sys = (*_system_levels.back()).filter_sys;

            return Solver::solve(*solver, vec_sol, vec_rhs, op_sys, filter_sys);
          }

          /// \copydoc BaseClass()::optimise()
          virtual void optimise() override
          {
            // fetch our finest levels
            //DomainLevelType& the_domain_level = *domain_levels.back();
            SystemLevelType& the_system_level = *_system_levels.back();
            AssemblerLevelType& the_asm_level = *_assembler_levels.back();

            // create our RHS and SOL vectors
            GlobalSystemVectorR vec_rhs = the_asm_level.assemble_rhs_vector(the_system_level);
            GlobalSystemVectorL vec_sol = the_asm_level.assemble_sol_vector(the_system_level);
            // solve
            this->apply(vec_sol, vec_rhs);

            (*the_system_level.coords_buffer).convert(*vec_sol);
            buffer_to_mesh();

            prepare(vec_sol);

          }

      };

      template<typename Space_>
      class DuDvFunctionalAssemblerLevel : public MeshoptAssemblerLevel<Space_>
      {
        public:
          typedef Control::Meshopt::MeshoptAssemblerLevel<Space_> BaseClass;
          typedef typename BaseClass::TrafoType TrafoType;

          static constexpr int shape_dim = BaseClass::MeshType::ShapeType::dimension;

          // Copy baseclass constructors
          using BaseClass::BaseClass;

          template<typename SystemLevel_>
          void assemble_system_matrix(SystemLevel_& sys_level)
          {
            (*(sys_level.op_sys)).assemble_system_matrix();
          }
      }; // class DuDvFunctionalAssemblerLevel

    } // namespace Meshopt
  } // namespace Control
} // namespace FEAT

#endif// FEAT_CONTROL_MESHOPT_DUDV_FUNCTIONAL_CONTROL_HPP
