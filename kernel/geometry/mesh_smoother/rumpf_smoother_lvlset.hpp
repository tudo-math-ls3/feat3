#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_SMOOTHER_LVLSET_HPP
#define KERNEL_GEOMETRY_RUMPF_SMOOTHER_LVLSET_HPP 1

#include <kernel/geometry/mesh_smoother/rumpf_smoother.hpp>
// For projecting the levelset function to a vertex vector
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional.hpp>
// DEBUG
#include <kernel/geometry/export_vtk.hpp>


namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Baseclass for a family of variational mesh optimisation algorithms using levelset functions.
     *
     * That is, for all mesh optimisation algorithms derived from Martin Rumpf's paper
     *
     * M. Rumpf: A variational approach to optimal meshes, Numerische Mathematik 72 (1996), pp. 523 - 540,
     *
     * extended by a levelset formulation for representing internal surfaces by Steffen Basting and Martin Weismann
     *
     * S. Basting and M. Weismann: A hybrid level set-front tracking finite element approach for fluid-structure
     * interaction and two-phase flow applications, Journal of Computational Physics 255(2013), 228 - 244.
     *
     * \tparam SpaceType_
     * Finite element space for the levelset function
     *
     * \tparam FunctionalType_
     * Functional used for defining mesh quality.
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \tparam H_EvalType
     * Local meshsize evaluator
     *
     * \author Jordi Paul
     *
     */
    template
    <
      typename SpaceType_,
      typename FunctionalType_,
      typename TrafoType_,
      typename DataType_,
      typename MemType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, DataType_>
    >
    class RumpfSmootherLevelset :
      public RumpfSmoother
      <
        FunctionalType_,
        TrafoType_,
        DataType_,
        MemType_,
        H_EvalType_
      >
    {
      public :
        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Type for the boundary mesh
        typedef typename Geometry::CellSubSet<ShapeType> BoundaryType;
        /// Factory for the boundary mesh
        typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;
        /// Who's my daddy?
        typedef RumpfSmoother<FunctionalType_,TrafoType_,DataType_, MemType_,H_EvalType_> BaseClass;
        /// Vector types for element sizes etc.
        typedef LAFEM::DenseVector<MemType_, DataType_> VectorType;
        /// Since the functional contains a ShapeType, these have to be the same
        //static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        //"ShapeTypes of the transformation / functional have to agree" );

        /// Vector with the original coordinates
        VectorType _coords_org[TrafoType_::MeshType::world_dim];

        /// FE space for the levelset function
        SpaceType_ _lvlset_space;
        /// DoF vector for the levelset function
        VectorType _lvlset_vec;
        /// DoF vector for the levelset function values in the mesh vertices
        VectorType _lvlset_vtx_vec;
        /// DoF vector the the gradient of the levelset function
        VectorType _lvlset_grad_vec[TrafoType_::MeshType::world_dim];
        /// DoF vector the the gradient of the levelset function
        VectorType _lvlset_grad_vtx_vec[TrafoType_::MeshType::world_dim];
        /// Align vertices/edges/faces with 0 levelset?
        bool _align_to_lvlset;
        /// Start mesh optimisation from original grid?
        bool _from_original;
        /// Use levelset-based r_adaptivity?
        bool _r_adaptivity;

      public:
        /// Tolerance up to which the levelset constraint has to be fulfilled
        DataType lvlset_constraint_tol;
        /// Last computed levelset constraint;
        DataType lvlset_constraint_last;

      public:
        /**
         * \copydoc FEAST::Geometry::RumpfSmoother::RumpfSmoother()
         *
         * \param[in] align_to_lvlset_
         * Align vertices/edges/faces with 0 levelset?
         *
         * \param[in] from_original_
         * Start mesh optimisation from original grid?
         *
         * \param[in] r_adaptivity_
         * Use levelset-based r_adaptivity?
         *
         **/
        explicit RumpfSmootherLevelset( TrafoType& trafo_, FunctionalType& functional_, bool align_to_lvlset_, bool from_original_, bool r_adaptivity_)
          : BaseClass(trafo_, functional_),
          _lvlset_space(trafo_),
          _lvlset_vec(_lvlset_space.get_num_dofs()),
          _lvlset_vtx_vec(this->_mesh.get_num_entities(0)),
          _align_to_lvlset(align_to_lvlset_),
          _from_original(from_original_),
          _r_adaptivity(r_adaptivity_),
          // This tolerance is pretty arbitrary and subject to change
          lvlset_constraint_tol( Math::pow(Math::eps<DataType_>(),DataType(0.5) ) ),
          lvlset_constraint_last(0)
          {
            // This sets the levelset penalty term to 0 if we do not align
            this->_functional.fac_lvlset = DataType(_align_to_lvlset);

            for(Index d = 0; d < this->_world_dim; ++d)
            {
              _coords_org[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
              _lvlset_grad_vec[d]= std::move(VectorType(_lvlset_space.get_num_dofs()));
              _lvlset_grad_vtx_vec[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
            }

            //for(Index d = 0; d < this->_world_dim; ++d)
            //  _grad[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
          }


        /// \copydoc ~RumpfSmoother()
        virtual ~RumpfSmootherLevelset()
        {
        };

        /// \brief Initialises member variables required for the first application of the smoother
        virtual void init() override
        {
          BaseClass::init();

          // Save original coordinates
          const typename MeshType::VertexSetType& vertex_set = this->_mesh.get_vertex_set();
          for(Index i(0); i < this->_nk; ++i)
          {
            for(Index d(0); d < this->_world_dim; ++d)
              _coords_org[d](i,DataType(vertex_set[i][d]));
          }

          this->prepare();
          this->compute_lambda();
          this->compute_h();
          this->initial_functional_value = compute_functional();
        }

        /// \copydoc RumpfSmoother::compute_functional()
        ///
        /// Variant containing the levelset penalty term
        virtual DataType compute_functional() override
        {
          DataType_ fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_, MeshType::world_dim> h;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <DataType_, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;

          // Reset last levelset constraint
          this->lvlset_constraint_last = DataType(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
            {
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));
              for(Index d = 0; d < this->_world_dim; d++)
              {
                // Get local coordinates
                x(d,j) = this->_coords[d](idx(cell,j));
                // Get local mesh size
                h(d) = this->_h[d](cell);
              }
            }

            fval += this->_lambda(cell)*this->_functional.compute_local_functional(x,h);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += this->_functional.compute_lvlset_penalty(lvlset_vals);
          }

          // Scale levelset constraint correctly and add it to the functional value
          fval += this->_functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last);
          return fval;
        } // compute_functional()

        /*
         * This makes the baseclass' compute_functional available, as one might want to compute the functional value
         * without the levelset contribution and now the compiler no longer complains about overloaded virtual
         */
        using BaseClass::compute_functional;
        /**
         * \copydoc RumpfSmootherLevelset::compute_functional()
         *
         * \param[in] func_norm
         * The contribution of the Frobenius norm for each cell
         *
         * \param[in] func_det
         * The contribution of the det term for each cell
         *
         * \param[in] func_det2
         * The contribution of the 1/det term for each cell
         *
         * \param[in] func_lvlset
         * The contribution of the levelset term for each cell
         *
         * Debug variant that saves the different contributions for each cell.
         **/
        virtual DataType compute_functional( DataType_* func_norm, DataType_* func_det, DataType_* func_det2, DataType* func_lvlset )
        {
          DataType_ fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_,MeshType::world_dim> h;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <DataType_, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;

          DataType_ norm_A(0), det_A(0), det2_A(0), lvlset_penalty(0);

          DataType_ func_norm_tot(0);
          DataType_ func_det_tot(0);
          DataType_ func_det2_tot(0);

          // Reset last levelset constraint
          this->lvlset_constraint_last = DataType(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
            {
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));
              for(Index d = 0; d < this->_world_dim; d++)
              {
                // Get local coordinates
                x(d,j) = this->_coords[d](idx(cell,j));
                // Get local mesh size
                h(d) = this->_h[d](cell);
              }
            }

            fval += this->_lambda(cell)*this->_functional.compute_local_functional(x,h, norm_A, det_A, det2_A);

            lvlset_penalty = this->_functional.compute_lvlset_penalty(lvlset_vals);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += lvlset_penalty;

            func_norm[cell] = this->_lambda(cell) * norm_A;
            func_det[cell] = this->_lambda(cell) * det_A;
            func_det2[cell] = this->_lambda(cell) * det2_A;
            func_lvlset[cell] =  lvlset_penalty;
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_det2_tot += func_det2[cell];
          }

          // Scale levelset constraint correctly and add it to the functional value
          fval += this->_functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last);

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) <<
            ", func_det2 = " << scientify(func_det2_tot) << ", func_lvlset = " <<
            scientify(this->_functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last)) << std::endl;

          return fval;
        } // compute_functional(func_norm, func_det, func_det2, func_lvlset)

        /// \copydoc RumpfSmoother::compute_gradient()
        virtual void compute_gradient()
        {
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_,MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <DataType_, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> lvlset_grad_vals;

          // Clear gradient vector
          for(Index i(0); i < this->_world_dim*this->_nk; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));
              for(Index d(0); d < this->_world_dim; ++d)
              {
                h(d) = this->_h[d](cell);
                // Get local coordinates
                x(d,j) = this->_coords[d](idx(cell,j));
                // Get levelset gradient
                lvlset_grad_vals(d,j) = this->_lvlset_grad_vtx_vec[d](idx(cell,j));
              }
            }

            // Compute local gradient and save it to grad_loc
            this->_functional.compute_local_grad(x, h, grad_loc);
            // This does not contain the weighting yet
            grad_loc *= this->_lambda(cell);
            // Add levelset penalty term, which is not weighted with lambda
            this->_functional.add_lvlset_penalty_grad(lvlset_vals, lvlset_grad_vals, grad_loc, lvlset_constraint_last);

            for(Index d(0); d < this->_world_dim; ++d)
            {
              // Add local contributions to global gradient vector
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                this->_grad[d*this->_nk + idx(cell,j)] += grad_loc(d,j);
            }
          }

          // Set gradient to 0 where bdry_id == -1
          // TODO: Convert this to proper filtering etc. and allow for more general BCs.
          for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
          {
            if(this->_bdry_id[i] == -1)
            {
              for(Index d(0); d < this->_world_dim; ++d)
                this->_grad[d*this->_nk + i] = DataType(0);
            }
          }

        } // compute_gradient

        /// \copydoc MeshSmoother::optimise()
        virtual void optimise() override
        {

          // Copy original coordinates if we want to start from the original grid
          if(this->_from_original)
          {
            for(Index d(0); d < this->_world_dim; ++d)
              this->_coords[d].copy(_coords_org[d]);

            this->set_coords();
          }

          // Quadratic penalty method for solving the constraint minimization problem as a sequence of unconstrained
          // minimisation problems if we want to align to the 0 levelset
          if(this->_align_to_lvlset)
          {
            // Value of the levelset constraint from the previous iteration
            DataType lvlset_constraint_previous(1);
            // The factor for the penalty term is increased by this factor in every iteration
            DataType big_factor_inc(1);

            int penalty_iterations(0);

            DataType fac_lvlset_org(this->_functional.fac_lvlset);
            while(lvlset_constraint_previous > this->lvlset_constraint_tol && this->_functional.fac_lvlset < DataType(0.5)*Math::Limits<DataType>::max())
            {
              penalty_iterations++;

              ALGLIBWrapper<RumpfSmootherLevelset<SpaceType_, FunctionalType_, TrafoType_, DataType_, MemType_, H_EvalType_>>::minimise_functional_cg(*this);
              // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
              this->set_coords();

              lvlset_constraint_previous = this->lvlset_constraint_last;
              // Increment for the penalty factor
              big_factor_inc *= Math::sqr(this->lvlset_constraint_last/lvlset_constraint_previous);

              // Increase penalty parameter by at least factor 5, very arbitrary
              this->_functional.fac_lvlset *= DataType(5)*FEAST::Math::max<DataType>(DataType(1),big_factor_inc);

              // DEBUG
              // std::cout << "Rumpf Smoother penalty iteration: " << penalty_iterations << ", last constraint: " << this->lvlset_constraint_last << ", penalty factor: " << _functional.fac_lvlset << ", fval = " << scientify(this->compute_functional()) << std::endl;
            }

            std::cout << "Rumpf Smoother penalty iterations: " << penalty_iterations <<
              ", last constraint: " << this->lvlset_constraint_last <<
              ", penalty factor: " << this->_functional.fac_lvlset <<
              ", fval = " << scientify(this->compute_functional())
              << std::endl;

            // Restore original values
            this->_functional.fac_lvlset = fac_lvlset_org;

          }
          else
          {
            // Standard method
            ALGLIBWrapper<RumpfSmootherLevelset<SpaceType_, FunctionalType_, TrafoType_, DataType_, MemType_, H_EvalType_>>::minimise_functional_cg(*this);
            // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
            this->set_coords();
          }

        }

        /**
         * \copydoc RumpfSmoother::prepare()
         *
         * In this case, it evaluates the levelset function on the current mesh and evaluates it in the current mesh's
         * vertices.
         *
         **/
        virtual void prepare()
        {
          // Project the levelset function to a grid vector on the new grid and...
          FEAST::Assembly::DiscreteVertexProjector::project(_lvlset_vtx_vec, _lvlset_vec, _lvlset_space);

          // ... do the same for its gradient
          FEAST::Assembly::DiscreteVertexProjector::project(this->_lvlset_grad_vtx_vec[0],this->_lvlset_grad_vec[0], this->_lvlset_space);
          FEAST::Assembly::DiscreteVertexProjector::project(this->_lvlset_grad_vtx_vec[1],this->_lvlset_grad_vec[1], this->_lvlset_space);

          // As the levelset might be used for r-adaptivity and it's vertex vector might have changed, re-compute
          if(this->_r_adaptivity)
          {
            this->compute_lambda();
            this->compute_h();
          }
        }

        /// \copydoc RumpfSmoother::compute_lambda()
        virtual void compute_lambda()
        {
          if(this->_r_adaptivity)
          {
            compute_lambda_lvlset();
          }
          else
          {
            BaseClass::compute_lambda_uniform();
          }

        } // compute_lambda

        /// \copydoc RumpfSmoother::compute_lambda()
        ///
        /// Includes condensation near the 0 levelset
        virtual void compute_lambda_lvlset()
        {
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));
          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_, MeshType::world_dim> h;

          DataType sum_lambda(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            DataType tmp(0);
            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              tmp += Math::sqr(this->_lvlset_vtx_vec(idx(cell,j)));

            this->_lambda(cell, DataType(1e-1)+Math::sqrt(tmp/DataType(1)));
            sum_lambda+=this->_lambda(cell);
          }

          // Scale so that sum(lambda) = 1
          sum_lambda = DataType(1)/sum_lambda;
          for(Index cell(0); cell < ncells; ++cell)
            this->_lambda(cell,sum_lambda*this->_lambda(cell));

        }

    }; // class RumpfSmootherLevelset

    /**
     * \brief Variant with the levelset function given explicitly
     *
     * \tparam AnalyticFunctionType_
     * Analytic function providing level set information
     *
     * \tparam AnalyticFunctionGrad0Type_
     * 1st component function of the AnalyticFunction's gradient
     *
     * \tparam AnalyticFunctionGrad0Type_
     * 2nd component function of the AnalyticFunction's gradient
     *
     * \tparam SpaceType_
     * FE space for the levelset function
     *
     * \tparam FunctionalType_
     * Functional used for determining mesh quality.
     *
     * \tparam TrafoType_
     * Our transformation.
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \tparam H_EvalType
     * Local meshsize evaluator
     *
     * \author Jordi Paul
     *
     **/
    template
    <
      typename AnalyticFunctionType_,
      typename AnalyticFunctionGrad0Type_,
      typename AnalyticFunctionGrad1Type_,
      typename SpaceType_,
      typename FunctionalType_,
      typename TrafoType_,
      typename DataType_,
      typename MemType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, DataType_>
    >
    class RumpfSmootherLevelsetAnalytic :
      public RumpfSmootherLevelset
      <
        SpaceType_,
        FunctionalType_,
        TrafoType_,
        DataType_,
        MemType_,
        H_EvalType_
      >
    {
      public:
        /// Transformation type
        typedef TrafoType_ TrafoType;
        /// Baseclass
        typedef RumpfSmootherLevelset
        <
          SpaceType_,
          FunctionalType_,
          TrafoType_,
          DataType_,
          MemType_,
          H_EvalType_
        > BaseClass;

        /// Analytic levelset function
        AnalyticFunctionType_& _analytic_lvlset;
        /// 1st component of its gradient
        AnalyticFunctionGrad0Type_& _analytic_lvlset_grad0;
        /// 2nd component of its gradient
        AnalyticFunctionGrad1Type_& _analytic_lvlset_grad1;
      public:
        /**
         * \copydoc RumpfSmootherLevelset()
         *
         * \param[in] align_to_lvlset_
         * Align vertices/edges/faces with 0 levelset?
         *
         * \param[in] from_original_
         * Start mesh optimisation from original grid?
         *
         * \param[in] r_adaptivity_
         * Use levelset-based r_adaptivity?
         *
         * \param[in] analytic_function_
         * The analytic function representing the levelset function
         *
         * \param[in] analytic_function_grad0_
         * The first component of the gradient
         *
         * \param[in] analytic_function_grad1_
         * The second component of the gradient
         *
         **/
        explicit RumpfSmootherLevelsetAnalytic(
          TrafoType& trafo_,
          FunctionalType_& functional_,
          bool align_to_lvlset_,
          bool from_original_,
          bool r_adaptivity_,
          AnalyticFunctionType_& analytic_function_,
          AnalyticFunctionGrad0Type_& analytic_function_grad0_,
          AnalyticFunctionGrad1Type_& analytic_function_grad1_)
          : BaseClass(trafo_, functional_, align_to_lvlset_, from_original_, r_adaptivity_),
          _analytic_lvlset(analytic_function_),
          _analytic_lvlset_grad0(analytic_function_grad0_),
          _analytic_lvlset_grad1(analytic_function_grad1_)
          {
          }

        /**
         * \copydoc RumpfSmootherLevelset::prepare()
         *
         * In this case, it evaluates the analytic levelset function on the current mesh and then evaluates it in the
         * current mesh's vertices.
         *
         **/
        virtual void prepare() override
        {
          this->set_coords();
          // Evaluate levelset function
          Assembly::Interpolator::project(this->_lvlset_vec, _analytic_lvlset, this->_lvlset_space);
          // Evaluate the gradient of the levelset function
          Assembly::Interpolator::project(this->_lvlset_grad_vec[0], _analytic_lvlset_grad0, this->_lvlset_space);
          Assembly::Interpolator::project(this->_lvlset_grad_vec[1], _analytic_lvlset_grad1, this->_lvlset_space);
          // BaseClass::prepare handles projection to grid vectors etc.
          BaseClass::prepare();

        }
    }; // class RumpfSmootherLevelsetAnalytic

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_RUMPF_SMOOTHER_LVLSET_HPP
