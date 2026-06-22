// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/adjacency/base.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/projection_filter_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/trace_assembler_basic_jobs.hpp>
#include <kernel/base_header.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/inverse_mapping.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/remote_lambda.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/projection_filter.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/memory_arbiter.hpp>

namespace FEAT::MeshOpt
{
  namespace
  {
    /// Utility type for treating a Memory::TypedView as a range
    template<typename T_>
    struct ViewRange
    {
      Memory::TypedView<T_> view;

      T_* begin()
      {
        return &view[0];
      }

      const T_* begin() const
      {
        return &view[0];
      }

      T_* end()
      {
        return &view[view.size()];
      }

      const T_* end() const
      {
        return &view[view.size()];
      }
    };
  } // namespace

  /// Quality measures of mesh deformation
  struct MeshDeformerQuality
  {
    /// Quality measure in L2 norm
    Real l2norm;
    /// Quality measure in point-wise maximum norm
    Real maxnorm;
  };

  /**
   * \brief Deforms mesh towards target cell sizes
   *
   * Based on "Mathematical and Numerical Analysis of a Robust and Efficient Grid
   * Deformation Method in the Finite Element Context" by Grajewski, Koester, and
   * Turek (https://doi.org/10.1137/050639387) the MeshDeformer can deform a mesh
   * towards a target cell size function. The mesh topology is left unchanged.
   *
   * \author Markus Muegge
   */
  template<typename MeshType_>
  class MeshDeformer
  {
  public:
    /// Mesh type
    using MeshType = MeshType_;
    /// MeshPart type
    using MeshPartType = Geometry::MeshPart<MeshType>;
    /// Number of coords of vertices
    static constexpr int dim = MeshType::world_dim;
    /// Vertex type
    using VertexType = typename MeshType::VertexType;
    /// Trafo type
    using TrafoType = Trafo::Standard::Mapping<MeshType>;
    /// Space type. This is the space that gets used to discretize the Poisson-Problem.
    using SpaceType = Space::Lagrange1::Element<TrafoType>;
    /// Data type
    using DataType = typename MeshType::CoordType;
    /// Local scalar vector type
    using LocalVectorType = LAFEM::DenseVector<DataType, Index>;
    /// Mirror type
    using MirrorType = LAFEM::VectorMirror<DataType, Index>;
    /// Gate type
    using GateType = Global::Gate<LocalVectorType, MirrorType>;
    /// Global scalar vector type
    using VectorType = Global::Vector<LocalVectorType, MirrorType>;
    /// Local blocked vector type
    using LocalBlockedVectorType = LAFEM::DenseVectorBlocked<DataType, Index, dim>;
    /// Global blocked vector type
    using BlockedVectorType = Global::Vector<LocalBlockedVectorType, MirrorType>;
    /// Local matrix type
    using LocalMatrixType = LAFEM::SparseMatrixCSR<DataType, Index>;
    /// Global matrix type
    using MatrixType = Global::Matrix<LocalMatrixType, MirrorType, MirrorType>;
    /// Local filter type
    using LocalFilterType = LAFEM::MeanFilter<DataType, Index>;
    /// Global filter type
    using FilterType = Global::Filter<LocalFilterType, MirrorType>;
    /// Local inverse mapping type
    using LocalInverseMappingType = Trafo::InverseMapping<TrafoType, DataType>;
    /// Global inverse mapping type. Used to determine rhs of ODEs.
    using InverseMappingType = Global::InverseMapping<LocalInverseMappingType>;
    /// Solver type
    using SolverType = Solver::SolverBase<VectorType>;
    /// Solver factory type for setting which solver to use for Poisson-Problem
    using SolverFactory = std::function<std::shared_ptr<SolverType>(const MatrixType&, const FilterType&)>;
    /// Projection filter type. For enforcing boundary conditions in a strong sense.
    using ProjectionFilterType = LAFEM::ProjectionFilter<MeshType, DataType, Index>;

  protected:
    /// Solver factory function
    SolverFactory _solver_factory = [](const MatrixType& m, const FilterType& f)
    { return Solver::new_gmres(m, f, 10); };

    /// Maximum number of iterations for each blend step
    Index _max_iters_per_blend_step = 1;

    /// Number of blend steps
    Index _blend_steps = 5;

    /// If true, boundary vertices are not allowed to move at all
    bool _keep_boundary_fixed{false};

    /// Gate for synchronization in MPI parallel build
    GateType* _gate = nullptr;
    /// MPI Communicator
    Dist::Comm _comm = Dist::Comm::world();

  public:
    /// Constructor
    explicit MeshDeformer(GateType* gate) : _gate(gate), _comm(_gate->get_comm()->comm_dup())
    {
    }

    /**
     * \brief Sets the maximum number of iterations during each blend step.
     *
     * \param[in] iters
     * Maximum number of iterations per blend step.
     *
     * The deformation depends on a good enough solution to the Poisson problem
     * created from the current cell sizes and the target cell sizes. We can
     * improve the solution by re-assembling the Poisson-problem (with the same
     * monitor function) on an already deformed mesh, solving that new problem,
     * and deforming the mesh further according to that improved solution.
     *
     * This sets the maximum number of such iterations. If the quality measure
     * is already sufficiently good before reaching the maximum number of steps
     * the deformer will continue to the next blend step (or finish, if it is on
     * the last blend step already).
     */
    void set_max_iters_per_blend_step(Index iters)
    {
      _max_iters_per_blend_step = iters;
    }

    /**
     * \brief Set number of blend steps
     *
     * \param[in] steps
     * Number of blend steps
     *
     * The deformation depends on a good-enough solution to the Poisson problem
     * created from the current cell sizes and the target cell sizes. If the
     * difference between those function is too large, the deformer might
     * produce widly degenerate meshes.
     *
     * The blend steps split the deformation into multiple steps, where each
     * steps monitor function is an interpolation between the current cell size
     * function and the target monitor function. That way each step can solve a
     * better conditioned Poisson problem and the mesh does not degenerate.
     */
    void set_blend_steps(Index steps)
    {
      _blend_steps = steps;
    }

    /**
     * \brief Set solver factory
     *
     * \param[in] factory
     * Solver factory
     *
     * The deformation needs to solve a Poisson problem in each iteration. Per
     * default the deformer uses a GMRES solver. You can set a factory function
     * to use a different solver if required.
     */
    void set_solver_factory(SolverFactory&& factory)
    {
      _solver_factory = std::move(factory);
    }

    /// If set to true, no vertex on the boundary will be moved
    void keep_boundary_fixed(bool fixed)
    {
      _keep_boundary_fixed = fixed;
    }

    /**
     * \brief Deform mesh to match monitor function
     *
     * \tparam MonitorFn_ Type of monitor function
     * Must have member DataType operator()(const VertexType&)
     *
     * \param[in,out] mesh_node
     * Mesh node whose mesh will be deformed
     *
     * \param[in] monitor_fn
     * Cell size target function
     */
    template<typename MonitorFn_>
    MeshDeformerQuality deform(Geometry::RootMeshNode<MeshType_>& mesh_node, const MonitorFn_& monitor_fn)
    {
      const Real blend_delta = Real(1) / Real(_blend_steps);

      // Prepare for deformation
      MeshType& mesh = *mesh_node.get_mesh();
      TrafoType trafo(mesh);
      SpaceType space(trafo);

      // Determine mesh boundary for handling boundary values
      Geometry::GlobalMaskedBoundaryFactory<MeshType> boundary_factory(mesh);
      for(auto& [rank, halo] : mesh_node.get_halo_map())
      {
        boundary_factory.add_halo(rank, *halo);
      }
      boundary_factory.compile(_comm);
      Geometry::MeshPart<MeshType> boundary(boundary_factory);

      // Determine surface features we want to maintain
      Geometry::MeshPart<MeshType> features = Assembly::FeatureMeshPartFactory<MeshType>::make(mesh, boundary);

      // Assemble projection filter to enforce strong boundary conditions
      ProjectionFilterType projection_filter;
      if(_keep_boundary_fixed)
      {
        const auto& tgs = boundary.template get_target_set<0>();
        const auto& vertices = mesh.get_vertex_set();
        for(Index i(0); i < boundary.get_num_entities(0); i++)
        {
          projection_filter.add_fixed_dof(tgs[i], vertices[tgs[i]]);
        }
      }
      else
      {
        // Normal boundary conditions
        // Project boundary back onto boundary or onto closest feature
        Assembly::ProjectionFilterAssembler<MeshType> projection_filter_assembler;
        projection_filter_assembler.add_mesh_part(boundary);
        projection_filter_assembler.add_feature_mesh_part(features);
        projection_filter_assembler.assemble(projection_filter, space, _gate);
      }

      Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
      domain_assembler.compile_all_elements();

      // Iteratively deform

      MeshDeformerQuality current_quality{};

      // Blend loop. Interpolate between current area function and target monitor function.
      // Loop i blends with factor (i + 1) * blend_delta
      for(Index i(0); i < _blend_steps; i++)
      {
        // Iterative loop
        // Re-assemble the Poisson problem on the already deformed mesh to get a better solution on the next solve.
        for(Index j(0); j < _max_iters_per_blend_step; j++)
        {
          // Deform the mesh
          _deform(
            mesh,
            trafo,
            space,
            domain_assembler,
            projection_filter,
            monitor_fn,
            DataType(i + 1) * blend_delta);

          // Determine quality of deformed mesh
          current_quality = _quality(mesh, space, domain_assembler, monitor_fn);

          // Exit early upon reaching sufficient quality
          // Value of 0.05 is taken from the paper
          if(current_quality.l2norm < 0.05)
          {
            break;
          }
        }
      }

      return current_quality;
    }

    /**
     * \brief Calculate the quality of the given mesh compared to the given monitor function
     *
     * \return The quality of the given mesh compared to the given monitor function
     */
    template<typename MonitorFn_>
    MeshDeformerQuality quality(Geometry::RootMeshNode<MeshType>& mesh_node, const MonitorFn_& monitor_fn)
    {
      // Set up trafo, space, and domain assembler
      MeshType& mesh = *mesh_node.get_mesh();

      TrafoType trafo(mesh);
      SpaceType space(trafo);

      Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
      domain_assembler.compile_all_elements();

      // Calculate quality
      return _quality(mesh, space, domain_assembler, monitor_fn);
    }

  protected:
    /// Deform the given mesh towards the cells sizes given by the monitor function.
    template<typename MonitorFn_>
    void _deform(
      MeshType& mesh,
      const TrafoType& trafo,
      const SpaceType& space,
      Assembly::DomainAssembler<TrafoType>& domain_assembler,
      const ProjectionFilterType& projection_filter,
      const MonitorFn_& monitor_lambda,
      Real blend_factor = 1.0)
    {
      // Compute area function g and monitor function f
      VectorType area_fn = _area_function(mesh);
      VectorType monitor_fn(_gate, space.get_num_dofs());
      {
        const auto& vertices = mesh.get_vertex_set();
        Memory::TypedView<DataType> view = monitor_fn.local().elements_view_w();
        for(Index i(0); i < space.get_num_dofs(); i++)
        {
          view[i] = monitor_lambda(vertices[i]);
        }
      }

      // We assume that the monitor_function will produce the same value for the same point on all ranks,
      // so no synchronization of monitor_fn is necessary.

      // Create inverse functions
      VectorType inv_area_fn(_gate, space.get_num_dofs());
      VectorType inv_monitor_fn(_gate, space.get_num_dofs());
      {
        if(blend_factor != DataType(1))
        {
          scale_to_same_integral_mean(domain_assembler, monitor_fn, area_fn, space);

          // Blend
          monitor_fn.scale(monitor_fn, blend_factor);
          monitor_fn.axpy(area_fn, DataType(1) - blend_factor);
        }

        // Invert
        inv_monitor_fn.component_invert(monitor_fn);
        inv_area_fn.component_invert(area_fn);

        scale_to_same_integral_mean(domain_assembler, inv_monitor_fn, inv_area_fn, space);
      }

      // Compute rhs (inv_monitor_fn - inv_area_fn)
      VectorType rhs(_gate, space.get_num_dofs());
      rhs.format();
      {
        VectorType tmp = inv_monitor_fn.clone();
        tmp.axpy(inv_area_fn, DataType(-1));

        Assembly::Common::IdentityOperator op;
        Assembly::assemble_bilinear_operator_apply_vector_1(
          domain_assembler,
          rhs.local(),
          tmp.local(),
          op,
          space,
          "auto-degree:5");
        rhs.sync_0();
      }

      // Solve -laplace w = rhs with integral mean 0
      VectorType solution(_gate, space.get_num_dofs());
      {
        solution.local().format();

        MatrixType matrix(_gate, _gate);
        Assembly::SymbolicAssembler::assemble_matrix_std1(matrix.local(), space);
        matrix.local().format();

        Assembly::Common::LaplaceOperator laplace_operator;
        Assembly::assemble_bilinear_operator_matrix_1(
          domain_assembler,
          matrix.local(),
          laplace_operator,
          space,
          "auto-degree:3");

        Assembly::MeanFilterAssembler mean_assembler;
        FilterType mean_filter;
        mean_assembler.assemble(mean_filter.local(), space, "auto-degree:3", DataType(0.0));

        mean_filter.local().filter_mat(matrix.local());
        mean_filter.filter_rhs(rhs);
        mean_filter.filter_sol(solution);

        auto solver = _solver_factory(matrix, mean_filter);

        solver->init();
        Solver::solve(*solver, solution, rhs, matrix, mean_filter);
        solver->done();
      }

      LocalBlockedVectorType velocity(space.get_num_dofs());
      velocity.format();
      Assembly::DiscreteVertexProjector::project_gradient(velocity, solution.local(), space);

      // Strictly enforce boundary conditions by projecting velocity onto the boundary
      projection_filter.filter_def(velocity);

      LocalInverseMappingType inv_mapping(trafo, 0.1, 1E-2);
      InverseMappingType global_unmapping(_comm, inv_mapping);

      // Solve ODEs
      {
        // Copy vertices into blocked vector so that we can do math on them
        LocalBlockedVectorType vertices(mesh.get_num_vertices());
        {
          const auto& mesh_vertices = mesh.get_vertex_set();
          Memory::TypedView<VertexType> view = vertices.elements_view_w();
          for(Index i(0); i < mesh.get_num_vertices(); i++)
          {
            view[i] = mesh_vertices[i];
          }
        }

        auto moved_vertices = runge_kutta_3(
          vertices,
          [&](DataType t, const LocalBlockedVectorType& vs)
          {
            LocalBlockedVectorType ode_rhs(mesh.get_num_vertices());
            Memory::TypedView<VertexType> view_out = ode_rhs.elements_view_w();

            // TODO: Abstract the following into a Global::DiscreteEvaluator class

            // Unmap all points
            std::vector<Global::GlobalInverseMappingData<DataType, dim>> unmapping =
              global_unmapping.unmap_points(ViewRange<VertexType>{vs.elements_view_r()});

            // Eval f, g, v
            struct EvalResult
            {
              DataType f;
              DataType g;
              VertexType v;

              static Dist::Typemap typemap()
              {
                return Dist::Typemap{}
                  .template add_entry<DataType>(offsetof(EvalResult, f))
                  .template add_entry<DataType>(offsetof(EvalResult, g))
                  .template add_entry<VertexType>(offsetof(EvalResult, v));
              }
            };

            std::vector<int> ranks;
            std::vector<Index> indices;
            std::vector<std::tuple<VertexType, Index, int>> args;

            for(Index i(0); i < unmapping.size(); i++)
            {
              if(unmapping[i].size() == 0)
              {
                // Unmapping for this point failed.
                // We call the lambda with the index of the vertex as a fallback
                ranks.push_back(_comm.rank());
                indices.push_back(i);
                args.emplace_back(VertexType(0.0), i, 0);
              }
              else
              {
                ranks.push_back(unmapping[i].ranks[0]);
                indices.push_back(i);
                args.emplace_back(unmapping[i].dom_points[0], unmapping[i].cells[0], 1);
              }
            }

            Global::ScalarRemoteLambda<EvalResult(VertexType, Index, int)> remote_eval(_comm);
            std::vector<EvalResult> evals = remote_eval.call(
              ranks,
              indices,
              args,
              [&](VertexType dom_point, Index cell, int success)
              {
                if(success == 0)
                {
                  const DataType f = inv_monitor_fn.local().elements_view_r()(cell);
                  const DataType g = inv_area_fn.local().elements_view_r()(cell);
                  const VertexType v = velocity.elements_view_r()(cell);

                  return EvalResult{f, g, v};
                }
                else
                {
                  Trafo::InverseMappingData<DataType, dim> inv_map_data;
                  inv_map_data.cells.push_back(cell);
                  inv_map_data.dom_points.push_back(dom_point);

                  const DataType f =
                    Assembly::DiscreteEvaluator::eval_fe_function(inv_map_data, inv_monitor_fn.local(), space)
                      .mean_value();
                  const DataType g =
                    Assembly::DiscreteEvaluator::eval_fe_function(inv_map_data, inv_area_fn.local(), space)
                      .mean_value();
                  const VertexType v =
                    Assembly::DiscreteEvaluator::eval_fe_function(inv_map_data, velocity, space).mean_value();

                  return EvalResult{f, g, v};
                }
              });

          FEAT_PRAGMA_OMP(parallel for)
          for(Index i = 0; i < mesh.get_num_vertices(); i++)
          {
            const DataType f = evals[i].f;
            const DataType g = evals[i].g;
            const VertexType v = evals[i].v;
            const DataType denom = t * f + (DataType(1.0) - t) * g;

            view_out[i] = (DataType(1.0) / denom) * v;
          }

          return std::move(ode_rhs);
          });

        // Copy moved vertices back into mesh
        {
          auto& mesh_vertices = mesh.get_vertex_set();
          Memory::TypedView<VertexType> view = moved_vertices.elements_view_r();
          for(Index i(0); i < mesh.get_num_vertices(); i++)
          {
            mesh_vertices[i] = view(i);
          }
        }
      }
    }

    template<typename MonitorFn_>
    MeshDeformerQuality _quality(
      MeshType& mesh,
      const SpaceType& space,
      Assembly::DomainAssembler<TrafoType>& assembler,
      const MonitorFn_& monitor_fn)
    {
      VectorType area = _area_function(mesh);
      VectorType monitor(_gate, mesh.get_num_vertices());
      {
        const auto& vertices = mesh.get_vertex_set();

        Memory::TypedView<DataType> monitor_view = monitor.local().elements_view_w();

        for(Index i(0); i < mesh.get_num_vertices(); i++)
        {
          monitor_view[i] = monitor_fn(vertices[i]);
        }
      }

      // Scale deformed_area_function
      scale_to_same_integral_mean(assembler, monitor, area, space);

      area.component_invert(area);
      VectorType quality_function(_gate, space.get_num_dofs());
      quality_function.component_product(area, monitor);

      double quality = Math::sqrt(
        Assembly::integrate_error_function<0>(
          assembler,
          Analytic::Common::ConstantFunction<dim>(1.0),
          quality_function.local(),
          space,
          "auto-degree:5")
          .norm_h0_sqr);
      double quality_max = Math::abs(quality_function.max_abs_element() - 1.0);

      _comm.allreduce(&quality, &quality, 1, Dist::op_sum);
      _comm.allreduce(&quality_max, &quality_max, 1, Dist::op_max);

      return {quality, quality_max};
    }

    /**
     * \brief Scale pair of vectors such that their integrals agree
     *
     * \tparam DomainAssembler_ DomainAssembler type
     * \tparam Vector_ Vector type
     * \tparam Space_ Space type
     *
     * \param[in] domain_assembler Domain assembler to use for calculating integrals
     * \param[in] a First vector
     * \param[inout] b Second Vector. Will be scaled to have same integral as \c a
     * \param[in] space Space of vectors \c a and \c b
     */
    void scale_to_same_integral_mean(
      Assembly::DomainAssembler<TrafoType>& domain_assembler,
      const VectorType& a,
      VectorType& b,
      const SpaceType& space)
    {
      auto int_a = Assembly::integrate_discrete_function<0>(domain_assembler, a.local(), space, "auto-degree:5").value;
      auto int_b = Assembly::integrate_discrete_function<0>(domain_assembler, b.local(), space, "auto-degree:5").value;

      // Add up integrals on all ranks
      _comm.allreduce(&int_a, &int_a, 1, Dist::op_sum);
      _comm.allreduce(&int_b, &int_b, 1, Dist::op_sum);

      // Choose alpha such that alpha * int_b = int_a
      b.scale(b, int_a / int_b);
    }

    /// Returns a vector of average cell areas of all cells surrounding each vertex
    VectorType _area_function(MeshType_& mesh)
    {
      TrafoType trafo(mesh);
      VectorType area_fn(_gate, mesh.get_num_vertices(), 0);

      // Precompute cell volumes
      std::vector<DataType> cell_volumes(mesh.get_num_entities(dim));
      for(Index i(0); i < mesh.get_num_entities(dim); i++)
      {
        cell_volumes[i] = trafo.compute_vol(i);
      }

      // Render cells@vertices
      const Adjacency::Graph c_at_v(Adjacency::RenderType::transpose, mesh.template get_index_set<dim, 0>());

      // Create frequency vector
      VectorType freqs(_gate, mesh.get_num_vertices());
      {
        Memory::TypedView<DataType> view = freqs.local().elements_view_w();
        for(Index i(0); i < mesh.get_num_vertices(); i++)
        {
          view[i] = DataType(c_at_v.degree(i));
        }
      }

      // Sync number of adjacent cells
      freqs.sync_0();

      {
        Memory::TypedView<DataType> view = area_fn.local().elements_view_w();
        Memory::TypedView<DataType> freq_view = freqs.local().elements_view_r();
        for(Index i(0); i < mesh.get_num_vertices(); i++)
        {
          // Compute average cell area
          DataType area(0.0);
          for(auto it = c_at_v.image_begin(i); it != c_at_v.image_end(i); ++it)
          {
            area += cell_volumes[*it];
          }

          view[i] = area / freq_view(i);
        }
      }

      area_fn.sync_0();

      return area_fn;
    }

    /// Solve ODE via RK3
    template<typename Rhs_>
    static LocalBlockedVectorType runge_kutta_3(const LocalBlockedVectorType& initial, const Rhs_& rhs)
    {
      const Real t_start(0);
      const Real t_end(1.0);
      const Real stepsize(0.05);

      LocalBlockedVectorType x = initial.clone();

      // Perform 100 steps of Runge-Kutta 3
      for(Real t(t_start); t < t_end; t += stepsize)
      {
        // std::cerr << "x = " << x << "\n";
        const LocalBlockedVectorType k1 = rhs(t, x);
        // std::cerr << "k1 = " << k1 << "\n";
        LocalBlockedVectorType x_tmp = x.clone();
        x_tmp.axpy(k1, stepsize / DataType(2));
        const LocalBlockedVectorType k2 = rhs(t + (stepsize / Real(2)), x_tmp);
        // std::cerr << "k2 = " << k2 << "\n";
        x_tmp = x.clone();
        x_tmp.axpy(k1, -stepsize);
        x_tmp.axpy(k2, DataType(2) * stepsize);
        const LocalBlockedVectorType k3 = rhs(t + stepsize, x_tmp);
        // std::cerr << "k3 = " << k3 << "\n";

        x.axpy(k1, stepsize / DataType(6));
        x.axpy(k2, stepsize / DataType(6) * 4);
        x.axpy(k3, stepsize / DataType(6));
        // std::cerr << "x = " << x << "\n";
      }

      return x;
    }
  };

} // namespace FEAT::MeshOpt
