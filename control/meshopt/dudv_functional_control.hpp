// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef FEAT_CONTROL_MESHOPT_DUDV_FUNCTIONAL_CONTROL_HPP
#define FEAT_CONTROL_MESHOPT_DUDV_FUNCTIONAL_CONTROL_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/nonlinear_functional.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/meshopt/dudv_functional.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/jacobi_precond.hpp>

#include <control/domain/domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Meshopt
    {
      template<typename>
      class DuDvFunctionalAssemblerLevel;

      /**
       * \brief Control class for DuDvFunctionals
       *
       * \tparam Mem_
       * The memory architecture of the local DuDvFunctional
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
       *
       * \author Jordi Paul
       *
       */
      template<typename Mem_, typename DT_, typename IT_, typename DomainControl_>
      class DuDvFunctionalControl
      : public MeshoptControlBase<DomainControl_>
      {
        public:
          /// Our memory architecture
          typedef Mem_ MemType;
          /// The floating point type
          typedef DT_ DataType;
          /// The index type
          typedef IT_ IndexType;

          /// Our base class
          typedef MeshoptControlBase<DomainControl_> BaseClass;

          /// The type of the domain control
          typedef DomainControl_ DomainControlType;
          /// Domain layers
          typedef typename DomainControl_::LayerType DomainLayerType;
          /// Domain levels
          typedef typename DomainControl_::LevelType DomainLevelType;

          /// The transformation we solve for
          typedef typename DomainLevelType::TrafoType TrafoType;

          /// Template-alias away the Trafo so the SystemLevel can take it as a template template parameter
          template<typename A, typename B, typename C>
          using LocalFunctionalType =  FEAT::Meshopt::DuDvFunctional<A, B, C, TrafoType>;

          /// The FE space the transformation lives in
          typedef typename LocalFunctionalType<Mem_, DT_, IT_>::SpaceType TrafoSpace;

          /// The underlying mesh type
          typedef typename DomainControl_::MeshType MeshType;
          /// The floating point type the mesh's coordinates use
          typedef typename MeshType::CoordType CoordType;

          /// Linear system of equations on one refinement level
          typedef QuadraticSystemLevel<Mem_, DT_, IT_, LocalFunctionalType> SystemLevelType;

          /// Inter-level transfer matrix
          typedef LAFEM::SparseMatrixBWrappedCSR<Mem_, DT_, IT_, MeshType::world_dim> TransferMatrixType;
          /// Global left vector type
          typedef typename SystemLevelType::GlobalSystemVectorL GlobalSystemVectorL;
          /// Global right vector type
          typedef typename SystemLevelType::GlobalSystemVectorR GlobalSystemVectorR;
          /// Global system matrix
          typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
          /// Global system filter
          typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;

          /// For every level of refinement, we have one system level
          std::deque<SystemLevelType*> _system_levels;

          /// The multigrid hierarchy for our solver
          std::shared_ptr<Solver::MultiGridHierarchy<
            GlobalSystemMatrix,
            GlobalSystemFilter,
            typename SystemLevelType::GlobalSystemTransfer
          >> _multigrid_hierarchy;

        public:
          /// Solver configuration
          PropertyMap& solver_config;
          /// The name of the section from solver_config we want to use
          String solver_name;
          /// The solver
          std::shared_ptr<Solver::SolverBase<GlobalSystemVectorR>> solver;
          /// Whether to reassemble the system matrix in every call of optimize
          const bool fixed_reference_domain;
          /// The level index (= number of refinements since the mesh file) to optimize the mesh on
          int meshopt_lvl;
          /// The position of this level in the deque of system levels
          size_t meshopt_lvl_pos;

          /**
           * \brief Constructor
           *
           * \param[in] dom_ctrl
           * The domaincontrol holding all geometry information for all levels
           *
           * \param[in] meshopt_lvl_
           * Index of the level to perform the mesh optimization on
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
           * \param[in] fixed_reference_domain_
           * If this is set to true, the system matrix is only assembled once, corresponding to a fixed reference
           * domain.
           *
           */
          explicit DuDvFunctionalControl(
            DomainControl_& dom_ctrl,
            const int meshopt_lvl_,
            const std::deque<String>& dirichlet_list, const std::deque<String>& slip_list,
            const String& solver_name_, PropertyMap& solver_config_,
            bool fixed_reference_domain_):
            BaseClass(dom_ctrl, dirichlet_list, slip_list),
            _system_levels(),
            solver_config(solver_config_),
            solver_name(solver_name_),
            solver(nullptr),
            fixed_reference_domain(fixed_reference_domain_),
            meshopt_lvl(meshopt_lvl_),
            meshopt_lvl_pos(~size_t(0))
          {
            XASSERT(meshopt_lvl >= -1);

            // If the input level was set to -1, take the max level of the domain control
            if(meshopt_lvl == -1)
            {
              meshopt_lvl = dom_ctrl.max_level_index();
            }

            XASSERT(meshopt_lvl<= dom_ctrl.max_level_index());

            // Now find the position of the coarsest mesh optimization level in the domain levels
            for(size_t i(0); i < dom_ctrl.size_physical(); ++i)
            {
              if(dom_ctrl.at(i)->get_level_index() == meshopt_lvl)
              {
                meshopt_lvl_pos  = i;
                break;
              }
            }

            XASSERT(meshopt_lvl_pos < dom_ctrl.size_physical());

            for(size_t i(0); i < dom_ctrl.size_physical(); ++i)
            {
              _system_levels.push_back( new SystemLevelType(
                dom_ctrl.at(i)->get_level_index(),
                dom_ctrl.at(i)->get_mesh_node(),
                dom_ctrl.at(i)->trafo,
                dirichlet_list, slip_list));

              // This assembles the system matrix numerically
              _system_levels.at(i)->local_functional.init();
            }

            // Now that _system_levels has the correct size, we can use get_num_levels()
            for(Index i(0); i < get_num_levels(); ++i)
            {
              _system_levels.at(i)->assemble_gates(dom_ctrl.at(i).layer());
            }

            // Assemble the transfer matrices on all levels except for the coarsest
            for(Index i(0); i+1 < get_num_levels(); ++i)
            {
              _system_levels.at(i)->assemble_system_transfer(*_system_levels.at(i+1));
            }

            // Now we need to assemble the filters on all levels before any shallow copies are made. This is done
            // by calling prepare()
            typename SystemLevelType::GlobalSystemVectorR vec_buf;
            vec_buf.local().convert(_system_levels.front()->coords_buffer.local());
            prepare(vec_buf);

            // create the solver
            this->_create_solver();
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
            while(!_system_levels.empty())
            {
              delete _system_levels.back();
              _system_levels.pop_back();
            }

            solver->done();
            _multigrid_hierarchy->done();
          }

          /// \copydoc BaseClass::name()
          virtual String name() const override
          {
            return "DuDvFunctionalControl<>";
          }

          /// \copydoc BaseClass::print()
          virtual String info() const override
          {
            const Index pad_width(30);

            String msg;

            msg += name().pad_back(pad_width, '.') + String(":") + String("\n");

            msg += String("level max/min").pad_back(pad_width, '.') + String(": ")
              + stringify(this->_dom_ctrl.max_level_index()) + String(" / ")
              + stringify(this->_dom_ctrl.min_level_index()) + String("\n");

            msg += String("optimization on level").pad_back(pad_width, '.') + String(": ")
              + stringify(meshopt_lvl) + String("\n");

            msg += String("Fixed reference domain").pad_back(pad_width, '.') + String(": ")
              + stringify(fixed_reference_domain) + String("\n");

            for(const auto& it : this->get_dirichlet_boundaries())
            {
              msg += String("Displacement BC on").pad_back(pad_width, '.') + String(": ") + it + String("\n");
            }

            for(const auto& it : this->get_slip_boundaries())
            {
              msg += String("Unilateral BC of place on").pad_back(pad_width, '.') + String(": ") + it + String("\n");
            }

            msg += String("DoF").pad_back(pad_width, '.') + String(": ")
              + stringify(_system_levels.front()->matrix_sys.columns()) + String("\n");

            try
            {
              msg += String("Solver") .pad_back(pad_width, '.') + String(": ")
                + FEAT::Statistics::get_formatted_solver_tree().trim() + String("\n");
            }
            catch(...)
            {
            }

            return msg;
          }

          /// \copydoc BaseClass::compute_cell_size_defect()
          virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max,
              CoordType& vol_min, CoordType& vol_max, CoordType& vol) const override
          {

            _system_levels.front()->local_functional.compute_cell_size_defect_pre_sync(vol_min, vol_max, vol);

            vol_min = _system_levels.front()->gate_sys.min(vol_min);
            vol_max = _system_levels.front()->gate_sys.max(vol_max);
            vol = _system_levels.front()->gate_sys.sum(vol);

            CoordType cell_size_defect =
              _system_levels.front()->local_functional.compute_cell_size_defect_post_sync(lambda_min, lambda_max, vol_min, vol_max, vol);

            lambda_min = _system_levels.front()->gate_sys.min(lambda_min);
            lambda_max = _system_levels.front()->gate_sys.max(lambda_max);
            cell_size_defect = _system_levels.front()->gate_sys.sum(cell_size_defect);

            return cell_size_defect;
          }

          /// \copydoc BaseClass::get_coords()
          virtual typename SystemLevelType::GlobalCoordsBuffer& get_coords() override
          {
            return _system_levels.front()->coords_buffer;
          }

          /// \copydoc BaseClass::get_num_levels()
          virtual size_t get_num_levels() const override
          {
            return _system_levels.size();
          }

          /// \copydoc BaseClass::buffer_to_mesh()
          virtual void buffer_to_mesh() override
          {
            // Write finest level
            _system_levels.front()->local_functional.buffer_to_mesh();

            // Get the coords buffer on the finest level
            const auto& coords_buffer_loc = _system_levels.front()->local_functional.get_coords();

            // Transfer fine coords buffer to coarser levels and perform buffer_to_mesh
            for(size_t level(0); level < get_num_levels(); ++level)
            {
              Index ndofs(_system_levels.at(level)->local_functional.trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::LocalCoordsBuffer
                vec_level(coords_buffer_loc, ndofs, Index(0));

              _system_levels.at(level)->local_functional.get_coords().copy(vec_level);
              _system_levels.at(level)->local_functional.buffer_to_mesh();
            }
          }

          /// \copydoc BaseClass::mesh_to_buffer()
          virtual void mesh_to_buffer() override
          {
            // Write finest level
            _system_levels.front()->local_functional.mesh_to_buffer();

            // Get the coords buffer on the finest level
            const auto& coords_buffer_loc = _system_levels.front()->coords_buffer.local();

            // Transfer fine coords buffer to coarser levels and perform buffer_to_mesh
            for(size_t level(0); level < get_num_levels(); ++level )
            {
              Index ndofs(_system_levels.at(level)->local_functional.trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::LocalCoordsBuffer
                vec_level(coords_buffer_loc, ndofs, Index(0));

              _system_levels.at(level)->local_functional.get_coords().copy(vec_level);
            }

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
                sys_lvl->local_functional.add_to_vtk_exporter(exporter);

                for(auto& it:sys_lvl->filter_sys.local().template at<0>())
                {
                  const String field_name("nu_"+it.first);
                  // WARNING: This explicitly assumes that the filter vector belong to a P1/Q1 space and thus "lives"
                  // in the mesh's vertices
                  exporter.add_vertex_vector(field_name, it.second.get_filter_vector());
                }
              }
            }
          }

          /**
           * \brief Numerically assembles the functional for evaluation
           *
           * This is mainly for using this functional as a preconditioner for the HyperelasticityFunctional so it
           * can be called separately.
           *
           */
          virtual void init_numeric()
          {
            if(!fixed_reference_domain)
            {
              for(size_t lvl(0); lvl < size_t(get_num_levels()); ++lvl)
              {
                _system_levels.at(lvl)->local_functional.assemble_system_matrix();
              }
            }
          }

          /// \copydoc BaseClass::prepare()
          virtual void prepare(const GlobalSystemVectorR& vec_state) override
          {
            typename SystemLevelType::LocalCoordsBuffer vec_buf;
            vec_buf.convert(vec_state.local());

            for(size_t level(0); level < get_num_levels(); ++level)
            {
              Index ndofs(_system_levels.at(level)->local_functional.trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::GlobalCoordsBuffer
                global_vec_level( &(_system_levels.at(level)->gate_sys), vec_buf, ndofs, Index(0));

              _system_levels.at(level)->local_functional.prepare(
                vec_buf, _system_levels.at(level)->filter_sys.local());

              _system_levels.at(level)->sync_system_filter();

            }

          }

          /**
           * \brief Applies the inverse of this functional's gradient to a right hand side
           *
           * \param[in,out] vec_sol
           * Initial guess, gets overwritten with the solution.
           *
           * \param[in] vec_rhs
           * Right hand side.
           *
           * This is for using this functional as a preconditioner i.e. for the HyperelasticityFunctional
           *
           * \returns
           * The status the solver finished with
           */
          virtual Solver::Status apply(GlobalSystemVectorR& vec_sol, const GlobalSystemVectorL& vec_rhs)
          {
            // Get our global system matrix and filter
            GlobalSystemMatrix& mat_sys = _system_levels.at(meshopt_lvl_pos)->matrix_sys;
            GlobalSystemFilter& filter_sys = _system_levels.at(meshopt_lvl_pos)->filter_sys;

            // Update the containers in the MatrixStock
            _multigrid_hierarchy->done_numeric();
            _multigrid_hierarchy->init_numeric();

            return Solver::solve(*solver, vec_sol, vec_rhs, mat_sys, filter_sys);
          }

          /// \copydoc BaseClass()::optimize()
          virtual void optimize() override
          {
            // Reassemble the system matrix
            init_numeric();

            // fetch our finest levels
            SystemLevelType& the_system_level = *_system_levels.at(meshopt_lvl_pos);

            // Create our RHS and SOL vectors
            GlobalSystemVectorR vec_rhs = the_system_level.assemble_rhs_vector();
            GlobalSystemVectorL vec_sol = the_system_level.assemble_sol_vector();

            // Let it be known to Statistics that it was Us who called the solver
            FEAT::Statistics::expression_target = name();

            // solve
            this->apply(vec_sol, vec_rhs);

            // Write the solution to the control object's buffer and the buffer to mesh
            typename SystemLevelType::LocalCoordsBuffer vec_buf;
            vec_buf.convert(vec_sol.local());
            the_system_level.coords_buffer.local().copy(vec_buf);

            the_system_level.local_functional.buffer_to_mesh();

            //If the mesh was not optimized on the finest domain level, we now need to prolongate the solution by:
            // - refine the coarse vertex set using the StandardRefinery
            // - copy the results to the CoordsBuffer
            // - convert the CoordsBuffer to a vector type that we can filter with the system's Dirichlet filter
            // - copy the filtered vector's contents back to the buffer and write the buffer to the mesh
            // - apply the nonlinear filter representing unilateral BCs of place by calling adapt() for the
            //   corresponding meshparts
            for(size_t pos(meshopt_lvl_pos); pos > size_t(0);)
            {
              --pos;

              auto& coarse_mesh = this->_dom_ctrl.at(pos+1)->get_mesh();
              auto& fine_mesh = this->_dom_ctrl.at(pos)->get_mesh();
              auto& fine_vtx = fine_mesh.get_vertex_set();

              // Refine coarse vertex set and write the result to the CoordsBuffer
              Geometry::StandardRefinery<MeshType> refinery(coarse_mesh);
              refinery.fill_vertex_set(fine_vtx);
              _system_levels.at(pos)->local_functional.mesh_to_buffer();
              // Convert the buffer to a filterable vector
              typename GlobalSystemVectorL::LocalVectorType vec_sol_lvl;
              vec_sol_lvl.convert(_system_levels.at(pos)->coords_buffer.local());
              // Filter this vector, copy back the contents and write the changes to the mesh
              auto& dirichlet_filters_lvl = _system_levels.at(pos)->filter_sys.local().template at<1>();
              dirichlet_filters_lvl.filter_sol(vec_sol_lvl);
              _system_levels.at(pos)->coords_buffer.local().copy(vec_sol_lvl);
              _system_levels.at(pos)->local_functional.buffer_to_mesh();
              // Now call adapt() on the slip boundaries
              auto* fine_mesh_node = this->_dom_ctrl.at(pos)->get_mesh_node();
              for(const auto& it:this->get_slip_boundaries())
              {
                fine_mesh_node->adapt_by_name(it);
              }
            }

            // Now we need to update all coarser levels
            for(size_t pos(meshopt_lvl_pos+1); pos < get_num_levels(); ++pos)
            {
              Index ndofs(_system_levels.at(pos)->local_functional.trafo_space.get_num_dofs());

              // At this point, what we really need is a primal restriction operator that restricts the FE function
              // representing the coordinate distribution to the coarser level. This is very simple for continuous
              // Lagrange elements (just discard the additional information from the fine level), but not clear in
              // the generic case. So we use an evil hack here:
              // Because of the underlying two level ordering, we just need to copy the first ndofs entries from
              // the fine level vector.
              typename SystemLevelType::GlobalCoordsBuffer
                global_sol_level( &(_system_levels.at(pos)->gate_sys), vec_sol.local(), ndofs, Index(0));

              // prepare() needs to be called of the local object, because there is a Global::Matrix around it
              // which knows nothing of prepare()
              _system_levels.at(pos)->local_functional.prepare(
                global_sol_level.local(), _system_levels.at(pos)->filter_sys.local());

              _system_levels.at(pos)->local_functional.get_coords().copy(global_sol_level.local());
              _system_levels.at(pos)->local_functional.buffer_to_mesh();
            }

          } // void optimize()

      protected:
        void _create_solver()
        {
          // first of all, let's fetch all the subsections for the DuDv solver
          auto* sec_linsol = solver_config.query_section(this->solver_name);
          auto* sec_mg_coarse = solver_config.query_section("DuDvMGCoarseSolver");
          auto* sec_mg_smooth = solver_config.query_section("DuDvMGSmoother");
          XASSERTM(sec_linsol    != nullptr, "mandatory DuDv solver config section is missing");
          XASSERTM(sec_mg_coarse != nullptr, "mandatory solver section [DuDvMGCoarseSolver] is missing");
          XASSERTM(sec_mg_smooth != nullptr, "mandatory solver section [DuDvMGSmoother] is missing");

          // our multigrid smoother is a simple Jacobi-smoother, so parse the damping parameter
          // as well as the number of smoothing steps
          Index smooth_steps = 4;
          auto smooth_steps_p = sec_mg_smooth->get_entry("steps");
          XASSERTM(smooth_steps_p.second, "DuDvMGSmoother.steps parameter is missing");
          smooth_steps_p.first.parse(smooth_steps);

          DT_ smooth_omega = DT_(0.5);
          auto smooth_damp_p = sec_mg_smooth->get_entry("omega");
          XASSERTM(smooth_damp_p.second, "DuDvMGSmoother.omega parameter is missing");
          smooth_damp_p.first.parse(smooth_omega);

          // now let's build the multigrid hierarchy
          _multigrid_hierarchy = std::make_shared<Solver::MultiGridHierarchy<GlobalSystemMatrix,
            GlobalSystemFilter, typename SystemLevelType::GlobalSystemTransfer>>(this->_dom_ctrl.size_virtual());

          // create smoothers for all levels except the coarsest one
          for(Index i(0); (i+1) < get_num_levels(); ++i)
          {
            const auto& lvl = *(_system_levels.at(i));
            auto jacobi = Solver::new_jacobi_precond(lvl.matrix_sys, lvl.filter_sys);
            auto smooth = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_omega, jacobi);
            smooth->set_min_iter(smooth_steps);
            smooth->set_max_iter(smooth_steps);
            _multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smooth, smooth, smooth);
          }

          // now let's create the coarse grid solver; this is always a PCG-Jacobi, which is configured
          // by the [DuDvMGCoarseSolver] section
          {
            const auto& lvl = *_system_levels.back();
            auto jacobi = Solver::new_jacobi_precond(lvl.matrix_sys, lvl.filter_sys);
            auto cg_pcg = Solver::new_pcg("DuDvMGCoarseSolver", sec_mg_coarse, lvl.matrix_sys, lvl.filter_sys, jacobi);
            _multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cg_pcg);
          }

          // now, let's create the actual multigrid solver object
          auto multigrid = Solver::new_multigrid(_multigrid_hierarchy, Solver::MultiGridCycle::V);

          // get the finest level
          const auto& lvl = *_system_levels.front();

          // finally, create the outer solver, which is always a PCG solver
          solver = Solver::new_pcg("DuDvLinearSolver", sec_linsol, lvl.matrix_sys, lvl.filter_sys, multigrid);

          // initialize the hierarchy and the solver
          _multigrid_hierarchy->init();
          solver->init();
        }
      };
    } // namespace Meshopt
  } // namespace Control
} // namespace FEAT

#endif// FEAT_CONTROL_MESHOPT_DUDV_FUNCTIONAL_CONTROL_HPP
