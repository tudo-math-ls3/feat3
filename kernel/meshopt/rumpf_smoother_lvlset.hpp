#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_SMOOTHER_LVLSET_HPP
#define KERNEL_MESHOPT_RUMPF_SMOOTHER_LVLSET_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/assembly/discrete_projector.hpp> // For projecting the levelset function to a vertex vector
#include <kernel/assembly/interpolator.hpp>
//#include <kernel/geometry/mesh_smoother/rumpf_functional.hpp>
#include <kernel/meshopt/rumpf_smoother.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /**
     * \brief Baseclass for a family of variational mesh optimisation algorithms using levelset functions.
     *
     * That is, for all mesh optimisation algorithms derived from Martin Rumpf's paper \cite Rum96
     *
     * extended by a levelset formulation for representing internal surfaces by Steffen Basting and Martin Weismann,
     * see \cite BW13
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \tparam FunctionalType_
     * Functional used for defining mesh quality.
     *
     * \tparam H_EvalType
     * Local meshsize evaluator
     *
     * \author Jordi Paul
     *
     */
    template
    <
      typename TrafoType_,
      typename FunctionalType_,
      typename LevelsetFunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, typename TrafoType_::MeshType::CoordType>
    >
    class RumpfSmootherLevelset :
      public RumpfSmootherBase
      <
        TrafoType_,
        FunctionalType_,
        H_EvalType_
      >
    {
      public :
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
        /// Type for the levelset part of the functional
        typedef LevelsetFunctionalType_ LevelsetFunctionalType;
        /// Own type for ALGLIBwrapper
        typedef RumpfSmootherLevelset<TrafoType, FunctionalType, LevelsetFunctionalType, H_EvalType_> MyType;

        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// We always use Index for now
        typedef Index IndexType;

        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<MemType, CoordType, IndexType> ScalarVectorType;
        /// Vector type for element scales etc.
        typedef LAFEM::DenseVectorBlocked<MemType, CoordType, IndexType, MeshType::world_dim> VectorType;
        /// Type for the filter enforcing boundary conditions
        typedef LAFEM::UnitFilterBlocked<MemType, CoordType, IndexType, MeshType::world_dim> FilterType;
        /// Finite Element space for the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space TrafoSpace;

        /// Who's my daddy?
        typedef RumpfSmootherBase<TrafoType_, FunctionalType_, H_EvalType_> BaseClass;
        /// Since the functional contains a ShapeType, these have to be the same
        static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value, "ShapeTypes of the transformation / functional have to agree" );

        /// Functional for levelset part
        LevelsetFunctionalType _lvlset_functional;
        /// FE space for the levelset function
        TrafoSpace _lvlset_space;
        /// DoF vector for the levelset function
        ScalarVectorType _lvlset_vec;
        /// DoF vector for the levelset function values in the mesh vertices
        ScalarVectorType _lvlset_vtx_vec;
        /// DoF vector the the gradient of the levelset function
        ScalarVectorType _lvlset_grad_vec[MeshType::world_dim];
        /// DoF vector the the gradient of the levelset function in the mesh vertices
        ScalarVectorType _lvlset_grad_vtx_vec[MeshType::world_dim];
        /// Vector with the original coordinates
        typename BaseClass::BaseClass::VertexVectorType _coords_org;

        /// Align vertices/edges/faces with 0 levelset?
        bool _align_to_lvlset;
        /// Use levelset-based r_adaptivity?
        bool _r_adaptivity;
        /// Recompute h in each iteration? By default, this is disabled due to poor numerical stability.
        bool _update_h;
        /// Regularisation parameter for levelset based optimal scales
        CoordType _r_adapt_reg;
        /// Exponent for levelset based optimal scales.
        /// The function is vol = vol(dist) = (_r_adapt_reg + dist)^_r_adapt_pow
        CoordType _r_adapt_pow;

      public:
        /// Tolerance up to which the levelset constraint has to be fulfilled
        CoordType lvlset_constraint_tol;
        /// Last computed levelset constraint;
        CoordType lvlset_constraint_last;

      public:
        /**
         * \copydoc FEAST::Geometry::RumpfSmoother::RumpfSmoother()
         *
         * \param[in] align_to_lvlset_
         * Align vertices/edges/faces with 0 levelset?
         *
         * \param[in] r_adaptivity_
         * Use levelset-based r_adaptivity?
         *
         **/
        explicit RumpfSmootherLevelset( TrafoType& trafo_, FunctionalType& functional_, LevelsetFunctionalType& lvlset_functional, bool align_to_lvlset_, bool r_adaptivity_)
          : BaseClass(trafo_, functional_),
          _lvlset_functional(lvlset_functional),
          _lvlset_space(trafo_),
          _lvlset_vec(_lvlset_space.get_num_dofs()),
          _lvlset_vtx_vec(trafo_.get_mesh().get_num_entities(0)),
          _coords_org(trafo_.get_mesh().get_num_entities(0)),
          _align_to_lvlset(align_to_lvlset_),
          _r_adaptivity(r_adaptivity_),
          _update_h(false),
          _r_adapt_reg(Math::pow( Math::eps<CoordType>(),CoordType(0.25)) ),
          _r_adapt_pow(CoordType(0.5)),
          // This tolerance is pretty arbitrary and subject to change
          lvlset_constraint_tol( Math::pow(Math::eps<CoordType>(),CoordType(0.75) ) ),
          lvlset_constraint_last(0)
          {
            // This sets the levelset penalty term to 0 if we do not align
            this->_lvlset_functional.fac_lvlset *= CoordType(_align_to_lvlset);

            for(int d(0); d < MeshType::world_dim; ++d)
            {
              _lvlset_grad_vec[d] = std::move(ScalarVectorType(_lvlset_space.get_num_dofs()));
              _lvlset_grad_vtx_vec[d] = std::move(ScalarVectorType(_lvlset_space.get_num_dofs()));
            }
          }

        /// \copydoc ~RumpfSmoother()
        virtual ~RumpfSmootherLevelset()
        {
        };

        /// \copydoc RumpfSmoother::print()
        virtual void print() override
        {
          std::cout << "RumpfSmootherLevelset characteristics: " << std::endl;

          std::cout << "Align to levelset:   " << _align_to_lvlset;
          if(_align_to_lvlset)
            std::cout << ", levelset constraint tolerance: " << scientify(lvlset_constraint_tol);
          std::cout << std::endl;

          std::cout << "Use r-adaptivity:    " << _r_adaptivity;
          if(_r_adaptivity)
            std::cout << ", update h in each iteration: " << _update_h;
            std::cout << ", r_adapt_reg = " << scientify(_r_adapt_reg)
            << ", r_adapt_pow = " << scientify(_r_adapt_pow);
          std::cout << std::endl;
          this->_functional.print();
        }

        /// \brief Initialises member variables required for the first application of the smoother
        virtual void init() override
        {
          BaseClass::init();

          // Save original coordinates
          const typename MeshType::VertexSetType& vertex_set = this->_mesh.get_vertex_set();

          for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
              _coords_org(i,vertex_set[i]);

          this->prepare();
          this->compute_lambda();
          this->compute_h();
          this->initial_functional_value = compute_functional();
        }

        /// \copydoc RumpfSmoother::compute_functional()
        ///
        /// Variant containing the levelset penalty term
        virtual CoordType compute_functional() override
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <CoordType, MeshType::world_dim> h;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;

          // Reset last levelset constraint
          this->lvlset_constraint_last = CoordType(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
            {
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,Index(j)));
              // Get local coordinates
              x[j] = this->_coords(idx(cell,Index(j)));
              // Get local mesh size
            }

            fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += this->_lvlset_functional.compute_lvlset_penalty(lvlset_vals);
          }

          // Scale levelset constraint correctly and add it to the functional value
          fval += this->_lvlset_functional.fac_lvlset/CoordType(2)*Math::sqr(lvlset_constraint_last);
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
         * \param[in] func_rec_det
         * The contribution of the 1/det term for each cell
         *
         * \param[in] func_lvlset
         * The contribution of the levelset term for each cell
         *
         * Debug variant that saves the different contributions for each cell.
         **/
        virtual CoordType compute_functional( CoordType* func_norm, CoordType* func_det, CoordType* func_rec_det, CoordType* func_lvlset )
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <CoordType,MeshType::world_dim> h;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;

          CoordType norm_A(0), det_A(0), rec_det_A(0), lvlset_penalty(0);

          CoordType func_norm_tot(0);
          CoordType func_det_tot(0);
          CoordType func_rec_det_tot(0);

          // Reset last levelset constraint
          this->lvlset_constraint_last = CoordType(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            // Get local mesh size
            h = this->_h(cell);
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,Index(j)));
              // Get local coordinates
              x[j] = this->_coords(idx(cell,Index(j)));
            }

            fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h, norm_A, det_A, rec_det_A);

            lvlset_penalty = this->_lvlset_functional.compute_lvlset_penalty(lvlset_vals);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += lvlset_penalty;

            func_norm[cell] = this->_mu(cell) * norm_A;
            func_det[cell] = this->_mu(cell) * det_A;
            func_rec_det[cell] = this->_mu(cell) * rec_det_A;
            func_lvlset[cell] =  lvlset_penalty;
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_rec_det_tot += func_rec_det[cell];
          }

          // Scale levelset constraint correctly and add it to the functional value
          fval += this->_lvlset_functional.fac_lvlset/CoordType(2)*Math::sqr(lvlset_constraint_last);

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) <<
            ", func_rec_det = " << scientify(func_rec_det_tot) << ", func_lvlset = " <<
            scientify(this->_lvlset_functional.fac_lvlset/CoordType(2)*Math::sqr(lvlset_constraint_last)) << std::endl;

          return fval;
        } // compute_functional(func_norm, func_det, func_rec_det, func_lvlset)

        /// \copydoc RumpfSmoother::compute_gradient()
        virtual void compute_gradient() override
        {
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <CoordType,MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> lvlset_grad_vals;

          // Clear gradient vector
          this->_grad.format();

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Global vertex/dof index
              Index i(idx(cell, Index(j)));
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(i);
              // Get local coordinates
              x[j] = this->_coords(i);
              // Get levelset gradient
              for(int d(0); d < MeshType::world_dim; ++d)
                lvlset_grad_vals(j,d) = this->_lvlset_grad_vtx_vec[d](i);
            }

            // Compute local gradient and save it to grad_loc
            this->_functional.compute_local_grad(x, h, grad_loc);
            // This does not contain the weighting yet
            grad_loc *= this->_mu(cell);
            // Add levelset penalty term, which is not weighted with lambda
            this->_lvlset_functional.add_lvlset_penalty_grad(lvlset_vals, lvlset_grad_vals, grad_loc, lvlset_constraint_last);

            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Global vertex/dof index
              Index i(idx(cell,Index(j)));
              Tiny::Vector<CoordType, MeshType::world_dim, MeshType::world_dim> tmp(this->_grad(i));
              tmp += grad_loc[j];

              this->_grad(i,tmp);
            }
          }

          this->_filter.filter_cor(this->_grad);

        } // compute_gradient

        /**
         * \brief Optimises the mesh using a fix point iteration wrt. h
         *
         * \param[in, out] total_grad_evals
         * Total number of gradient evaluations
         *
         * \param[in, out] total_iterations
         * Total number of optimiser iterations
         *
         * \param[out] termination_type
         * Termination type of the optimiser
         *
         */
        virtual void optimise_fixpoint_h(int& total_grad_evals, int& total_iterations, int& termination_type)
        {
          CoordType diff(1);
          CoordType tol(Math::pow(Math::eps<CoordType>(),CoordType(0.5) ));

          // coords_new = relaxation_parameter*coords + (1 - relaxation_parameter)*coords_old
          CoordType relaxation_parameter(0.9);

          for(Index iter(0); iter < 100; ++iter)
          {
            diff = CoordType(0);
            int iterations(0);
            int grad_evals(0);
            // Save old coordinates for explicit terms
            this->_coords_org.clone(this->_coords);

            ALGLIBWrapper<MyType>::minimise_functional_cg(grad_evals, iterations, termination_type,*this);

            total_grad_evals += grad_evals;
            total_iterations += iterations;

            for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
            {
              for(int d(0); d < MeshType::world_dim; ++d)
                diff+= Math::sqr(this->_coords_org(i)(d) - this->_coords(i)(d));
            }
            diff = Math::sqrt(diff/CoordType(this->_mesh.get_num_entities(0)));

            // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
            this->_coords_org.scale(this->_coords_org,CoordType(1)-relaxation_parameter);
            this->_coords.axpy(this->_coords, this->_coords_org, relaxation_parameter);

            this->set_coords();

            // DEBUG
            // std::cout << "Fixed point iteration " << iter <<", diff = " << scientify(diff) << ", mincg iterations: " << iterations << ", grad evals: " << grad_evals << std::endl;

            if(diff <= tol)
            {
              std::cout << iter+1 << " fixed point iterations, " << total_iterations << " mincg iterations, " << total_grad_evals << " grad evals, last terminationtype was " << termination_type << std::endl;
              return;
            }

            this->prepare();
            this->compute_lambda();
            this->compute_h();
          }
        }

        /// \copydoc MeshSmoother::optimise()
        virtual void optimise() override
        {
          int total_grad_evals(0);
          int total_iterations(0);
          int termination_type(0);

          // Copy original coordinates if we want to start from the original grid
          this->_coords_org.clone(this->_coords);

          this->prepare();
          this->compute_lambda();
          this->compute_h();

          if(this->_r_adaptivity && !this->_update_h)
            optimise_fixpoint_h(total_grad_evals, total_iterations, termination_type);
          // Standard method
          else
            if(!this->_align_to_lvlset)
          {
            ALGLIBWrapper<MyType>::minimise_functional_cg(total_grad_evals, total_iterations, termination_type,*this);
            // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
            this->set_coords();
            std::cout << total_iterations << " mincg iterations, " << total_grad_evals << " grad evals, terminationtype was " << termination_type << std::endl;
            return;
          }

          // Quadratic penalty method for solving the constraint minimization problem as a sequence of unconstrained
          // minimisation problems if we want to align to the 0 levelset
          if(this->_align_to_lvlset)
          {
            // The factor for the penalty term is increased by this factor in every iteration
            CoordType big_factor_inc(1);

            int penalty_iterations(0);
            this->compute_functional();
            // Value of the levelset constraint from the previous iteration
            CoordType lvlset_constraint_previous = this->lvlset_constraint_last;

            // Save old levelset penalty factor for restoring later
            CoordType fac_lvlset_org(this->_lvlset_functional.fac_lvlset);

            // The Penalty Iteration(TM)
            while(this->lvlset_constraint_last > this->lvlset_constraint_tol && this->_lvlset_functional.fac_lvlset < Math::sqrt(Math::Limits<CoordType>::max()))
            {
              int grad_evals(0);
              int iterations(0);

              penalty_iterations++;
              // Update penalty factor
              this->_lvlset_functional.fac_lvlset *= big_factor_inc;

              lvlset_constraint_previous = this->lvlset_constraint_last;

              ALGLIBWrapper<MyType>::minimise_functional_cg(grad_evals, iterations, termination_type,*this);

              total_grad_evals += grad_evals;
              total_iterations += iterations;

              // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
              this->set_coords();

              // DEBUG
              std::cout << "Rumpf Smoother penalty iteration: " << penalty_iterations << ", last constraint: "
              << this->lvlset_constraint_last << ", penalty factor: " << this->_lvlset_functional.fac_lvlset <<
              ", fval = " << scientify(this->compute_functional()) <<
              ", " << iterations << " mincg iterations, " << grad_evals <<
              " grad evals, terminationtype was " << termination_type << std::endl;

              // Increment for the penalty factor
              // Increase penalty parameter by at least factor 5, very arbitrary
              big_factor_inc = CoordType(5)*this->_lvlset_functional.fac_lvlset*
                Math::max<CoordType>(Math::sqr(lvlset_constraint_previous/this->lvlset_constraint_last), CoordType(1));
            }

            std::cout << "Rumpf Smoother penalty iterations: " << penalty_iterations <<
              ", last constraint: " << this->lvlset_constraint_last <<
              ", penalty factor: " << this->_lvlset_functional.fac_lvlset <<
              ", fval = " << scientify(this->compute_functional()) <<
              ", " << total_iterations << " mincg iterations, " << total_grad_evals <<
              " grad evals, terminationtype was " << termination_type << std::endl;

            // Restore original values
            this->_lvlset_functional.fac_lvlset = fac_lvlset_org;

          }

        }

        /**
         * \copydoc RumpfSmoother::prepare()
         *
         * In this case, it evaluates the levelset function on the current mesh and evaluates it in the current mesh's
         * vertices.
         *
         **/
        virtual void prepare() override
        {
          // Project the levelset function to a grid vector on the new grid and...
          FEAST::Assembly::DiscreteVertexProjector::project(_lvlset_vtx_vec, _lvlset_vec, _lvlset_space);

          // ... do the same for its gradient
          for(int d(0); d < MeshType::world_dim; ++d)
            FEAST::Assembly::DiscreteVertexProjector::project(this->_lvlset_grad_vtx_vec[d],this->_lvlset_grad_vec[d], this->_lvlset_space);

          // As the levelset might be used for r-adaptivity and it's vertex vector might have changed, re-compute
          if(this->_update_h && this->_r_adaptivity)
          {
            this->compute_lambda();

            // DEBUG
            //for(Index cell(0); cell < this->_mesh.get_num_entities(ShapeType::dimension); cell++ )
            //  this->_mu(cell,this->_lambda(cell));

            this->compute_h();
          }
        }

        /// \copydoc RumpfSmoother::compute_lambda()
        virtual void compute_lambda() override
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
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <CoordType, MeshType::world_dim> h;

          CoordType sum_lambda(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            CoordType tmp(0);

            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              tmp += this->_lvlset_vtx_vec(idx(cell,Index(j)));

            tmp = tmp/CoordType(Shape::FaceTraits<ShapeType,0>::count);
            tmp = Math::pow(this->_r_adapt_reg + Math::abs(tmp),this->_r_adapt_pow);

            this->_lambda(cell, tmp);
            sum_lambda+=this->_lambda(cell);
          }

          // Scale so that sum(lambda) = 1
          sum_lambda = CoordType(1)/sum_lambda;
          this->_lambda.scale(this->_lambda, sum_lambda);

        }

        /// \brief Sets parameters for compute_lambda_lvlset
        void set_r_adapt_params(CoordType reg_, CoordType pow_)
        {
          this->_r_adapt_reg = reg_;
          this->_r_adapt_pow = pow_;
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
     * \tparam CoordType
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
      typename TrafoType_,
      typename FunctionalType_,
      typename LevelsetFunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, typename TrafoType_::MeshType::CoordType>
    >
    class RumpfSmootherLevelsetAnalytic :
      public RumpfSmootherLevelset
      <
        TrafoType_,
        FunctionalType_,
        LevelsetFunctionalType_,
        H_EvalType_
      >
    {
      public:
        /// Type for the transformation
        typedef TrafoType_ TrafoType;

        /// Baseclass
        typedef RumpfSmootherLevelset
        <
          TrafoType_,
          FunctionalType_,
          LevelsetFunctionalType_,
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
          LevelsetFunctionalType_& lvlset_functional_,
          bool align_to_lvlset_,
          bool r_adaptivity_,
          AnalyticFunctionType_& analytic_function_,
          AnalyticFunctionGrad0Type_& analytic_function_grad0_,
          AnalyticFunctionGrad1Type_& analytic_function_grad1_)
          : BaseClass(trafo_, functional_, lvlset_functional_, align_to_lvlset_, r_adaptivity_),
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
          //TODO
          // Evaluate the gradient of the levelset function
          Assembly::Interpolator::project(this->_lvlset_grad_vec[0], _analytic_lvlset_grad0, this->_lvlset_space);
          Assembly::Interpolator::project(this->_lvlset_grad_vec[1], _analytic_lvlset_grad1, this->_lvlset_space);
          // BaseClass::prepare handles projection to grid vectors etc.
          BaseClass::prepare();

        }
    }; // class RumpfSmootherLevelsetAnalytic

  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_RUMPF_SMOOTHER_LVLSET_HPP
