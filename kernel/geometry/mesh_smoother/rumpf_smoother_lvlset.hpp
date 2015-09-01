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
     * That is, for all mesh optimisation algorithms derived from Martin Rumpf's paper \cite Rum96
     *
     * extended by a levelset formulation for representing internal surfaces by Steffen Basting and Martin Weismann,
     * see \cite BW13
     *
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \tparam FunctionalType_
     * Functional used for defining mesh quality.
     *
     * \tparam SpaceType_
     * Finite element space for the levelset function
     *
     * \tparam H_EvalType
     * Local meshsize evaluator
     *
     * \author Jordi Paul
     *
     */
    template
    <
      typename DataType_,
      typename MemType_,
      typename TrafoType_,
      typename FunctionalType_,
      typename LevelsetFunctionalType_,
      typename SpaceType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, DataType_>
    >
    class RumpfSmootherLevelset :
      public RumpfSmootherBase
      <
        DataType_,
        MemType_,
        TrafoType_,
        FunctionalType_,
        H_EvalType_
      >
    {
      public :
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
        /// Type for the levelset part of the functional
        typedef LevelsetFunctionalType_ LevelsetFunctionalType;
        /// Own type for ALGLIBwrapper
        typedef RumpfSmootherLevelset<DataType, MemType, TrafoType, FunctionalType, LevelsetFunctionalType, SpaceType_, H_EvalType_> MyType;

        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Type for the boundary mesh
        typedef typename Geometry::MeshPart<MeshType> BoundaryType;
        /// Factory for the boundary mesh
        typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;
        /// Who's my daddy?
        typedef RumpfSmootherBase<DataType_, MemType_, TrafoType_, FunctionalType_, H_EvalType_> BaseClass;
        /// Vector types for element sizes etc.
        typedef LAFEM::DenseVector<MemType_, DataType_> VectorType;
        /// Since the functional contains a ShapeType, these have to be the same
        static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value, "ShapeTypes of the transformation / functional have to agree" );

        /// Vector with the original coordinates
        VectorType _coords_org[TrafoType_::MeshType::world_dim];

        /// Functional for levelset part
        LevelsetFunctionalType _lvlset_functional;
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
        /// Use levelset-based r_adaptivity?
        bool _r_adaptivity;
        /// Recompute h in each iteration? By default, this is disabled due to poor numerical stability.
        bool _update_h;
        /// Regularisation parameter for levelset based optimal scales
        DataType _r_adapt_reg;
        /// Exponent for levelset based optimal scales.
        /// The function is vol = vol(dist) = (_r_adapt_reg + dist)^_r_adapt_pow
        DataType _r_adapt_pow;

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
          _align_to_lvlset(align_to_lvlset_),
          _r_adaptivity(r_adaptivity_),
          _update_h(false),
          _r_adapt_reg(Math::pow( Math::eps<DataType_>(),DataType(0.25)) ),
          _r_adapt_pow(DataType(0.5)),
          // This tolerance is pretty arbitrary and subject to change
          lvlset_constraint_tol( Math::pow(Math::eps<DataType_>(),DataType(0.75) ) ),
          lvlset_constraint_last(0)
          {
            // This sets the levelset penalty term to 0 if we do not align
            this->_lvlset_functional.fac_lvlset *= DataType(_align_to_lvlset);

            for(Index d = 0; d < MeshType::world_dim; ++d)
            {
              _coords_org[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
              _lvlset_grad_vec[d]= std::move(VectorType(_lvlset_space.get_num_dofs()));
              _lvlset_grad_vtx_vec[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
            }

            //for(Index d = 0; d < MeshType::world_dim; ++d)
            //  _grad[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
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
          {
            for(int d(0); d < MeshType::world_dim; ++d)
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
              for(Index d = 0; d < MeshType::world_dim; d++)
              {
                // Get local coordinates
                x(d,j) = this->_coords[d](idx(cell,j));
                // Get local mesh size
                h(d) = this->_h[d](cell);
              }
            }

            fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += this->_lvlset_functional.compute_lvlset_penalty(lvlset_vals);
          }

          // Scale levelset constraint correctly and add it to the functional value
          fval += this->_lvlset_functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last);
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
        virtual DataType compute_functional( DataType_* func_norm, DataType_* func_det, DataType_* func_rec_det, DataType* func_lvlset )
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

          DataType_ norm_A(0), det_A(0), rec_det_A(0), lvlset_penalty(0);

          DataType_ func_norm_tot(0);
          DataType_ func_det_tot(0);
          DataType_ func_rec_det_tot(0);

          // Reset last levelset constraint
          this->lvlset_constraint_last = DataType(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index j(0); j < Index(Shape::FaceTraits<ShapeType,0>::count); j++)
            {
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));
              for(Index d(0); d < Index(MeshType::world_dim); d++)
              {
                // Get local coordinates
                x(d,j) = this->_coords[d](idx(cell,j));
                // Get local mesh size
                h(d) = this->_h[d](cell);
              }
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
          fval += this->_lvlset_functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last);

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) <<
            ", func_rec_det = " << scientify(func_rec_det_tot) << ", func_lvlset = " <<
            scientify(this->_lvlset_functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last)) << std::endl;

          return fval;
        } // compute_functional(func_norm, func_det, func_rec_det, func_lvlset)

        /// \copydoc RumpfSmoother::compute_gradient()
        virtual void compute_gradient() override
        {
          // Total number of cells in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));
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
          for(Index i(0); i < MeshType::world_dim*nvertices; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index j(0); j < Index(Shape::FaceTraits<ShapeType,0>::count); ++j)
            {
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));
              for(Index d(0); d < Index(MeshType::world_dim); ++d)
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
            grad_loc *= this->_mu(cell);
            // Add levelset penalty term, which is not weighted with lambda
            this->_lvlset_functional.add_lvlset_penalty_grad(lvlset_vals, lvlset_grad_vals, grad_loc, lvlset_constraint_last);

            for(Index d(0); d < Index(MeshType::world_dim); ++d)
            {
              // Add local contributions to global gradient vector
              for(Index j(0); j < Index(Shape::FaceTraits<ShapeType,0>::count); ++j)
                this->_grad[d*nvertices + idx(cell,j)] += grad_loc(d,j);
            }
          }

          this->_filter_grad();

        } // compute_gradient

        virtual void optimise_fixpoint_h(int& total_grad_evals, int& total_iterations, int& termination_type)
        {
          DataType diff(1);
          DataType tol(Math::pow(Math::eps<DataType_>(),DataType(0.75) ));

          // coords_new = relaxation_parameter*coords + (1 - relaxation_parameter)*coords_old
          DataType relaxation_parameter(0.9);

          for(Index iter(0); iter < 1000; ++iter)
          {
            diff = DataType(0);
            int iterations(0);
            int grad_evals(0);
            // Save old coordinates for explicit terms
            for(Index d(0); d < MeshType::world_dim; ++d)
              this->_coords_org[d].clone(this->_coords[d]);

            ALGLIBWrapper<MyType>::minimise_functional_cg(grad_evals, iterations, termination_type,*this);

            total_grad_evals += grad_evals;
            total_iterations += iterations;

            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
                diff+= Math::sqr(this->_coords_org[d](i) - this->_coords[d](i));
            }
            diff = Math::sqrt(diff/DataType(this->_mesh.get_num_entities(0)));

            // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              this->_coords_org[d].scale(this->_coords_org[d],DataType(1)-relaxation_parameter);
              this->_coords[d].axpy(this->_coords[d], this->_coords_org[d], relaxation_parameter);
            }
            // DEBUG
            std::cout << "Fixed point iteration " << iter <<", diff = " << scientify(diff) << ", mincg iterations: " << iterations << ", grad evals: " << grad_evals << std::endl;

            this->set_coords();
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
          for(Index d(0); d < MeshType::world_dim; ++d)
            this->_coords_org[d].clone(this->_coords[d]);

          this->prepare();
          this->compute_lambda();
          this->compute_h();

          if(this->_r_adaptivity && !this->_update_h)
            optimise_fixpoint_h(total_grad_evals, total_iterations, termination_type);
          // Standard method
          else if(!this->_align_to_lvlset)
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
            DataType big_factor_inc(1);

            int penalty_iterations(0);
            this->compute_functional();
            // Value of the levelset constraint from the previous iteration
            DataType lvlset_constraint_previous = this->lvlset_constraint_last;

            // Save old levelset penalty factor for restoring later
            DataType fac_lvlset_org(this->_lvlset_functional.fac_lvlset);

            // The Penalty Iteration(TM)
            while(this->lvlset_constraint_last > this->lvlset_constraint_tol && this->_lvlset_functional.fac_lvlset < Math::sqrt(Math::Limits<DataType>::max()))
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
              big_factor_inc = DataType(5)*this->_lvlset_functional.fac_lvlset*
                Math::max<DataType>(Math::sqr(lvlset_constraint_previous/this->lvlset_constraint_last), DataType(1));
            }

            std::cout << "Rumpf Smoother penalty iterations: " << penalty_iterations <<
              ", last constraint: " << this->lvlset_constraint_last <<
              ", penalty factor: " << this->_lvlset_functional.fac_lvlset <<
              ", fval = " << scientify(this->compute_functional()) <<
              ", " << total_iterations << " mincg iterations, " << total_grad_evals <<
              " grad evals, terminationtype was " << termination_type << std::endl;

            // Restore original values
            this->_lvlset_functional.fac_lvlset = fac_lvlset_org;

            //bdry_id_org = new int[this->_mesh.get_num_entities(0)];
            //const DataType eps = DataType(1e-4);//Math::sqrt(Math::eps<DataType>());
            //for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
            //{
            //  bdry_id_org[i] = this->_bdry_id[i];
            //  if(Math::abs(this->_lvlset_vtx_vec(i)) < eps)
            //    this->_bdry_id[i] = -1;
            //}

          }

          //if(this->_align_to_lvlset)
          //{
          //  for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
          //    this->_bdry_id[i] = bdry_id_org[i];

          //  delete(bdry_id_org);
          //}

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
          FEAST::Assembly::DiscreteVertexProjector::project(this->_lvlset_grad_vtx_vec[0],this->_lvlset_grad_vec[0], this->_lvlset_space);
          FEAST::Assembly::DiscreteVertexProjector::project(this->_lvlset_grad_vtx_vec[1],this->_lvlset_grad_vec[1], this->_lvlset_space);

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
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_, MeshType::world_dim> h;

          DataType sum_lambda(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            DataType tmp(0);

            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              tmp += this->_lvlset_vtx_vec(idx(cell,j));

            tmp = tmp/DataType(Shape::FaceTraits<ShapeType,0>::count);
            tmp = Math::pow(this->_r_adapt_reg + Math::abs(tmp),this->_r_adapt_pow);

            this->_lambda(cell, tmp);
            sum_lambda+=this->_lambda(cell);
          }

          // Scale so that sum(lambda) = 1
          sum_lambda = DataType(1)/sum_lambda;
          for(Index cell(0); cell < ncells; ++cell)
          {
            this->_lambda(cell,sum_lambda*this->_lambda(cell));
            // DEBUG: Set mu to lambda
            //this->_mu(cell,this->_lambda(cell));
          }

          // DEBUG: Set mu to 1/lambda and scale such that sum(mu) = 1
          //DataType sum_mu(0);
          //for(Index cell(0); cell < ncells; ++cell)
          //{
          //  this->_mu(cell,DataType(1)/this->_lambda(cell));
          //  sum_mu+=this->_mu(cell);
          //}

          //// Scale so that sum(mu) = 1
          //sum_mu= DataType(1)/sum_mu;
          //for(Index cell(0); cell < ncells; ++cell)
          //  this->_mu(cell,sum_mu*this->_mu(cell));

        }

        /// \brief Sets parameters for compute_lambda_lvlset
        void set_r_adapt_params(DataType reg_, DataType pow_)
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
      typename DataType_,
      typename MemType_,
      typename TrafoType_,
      typename FunctionalType_,
      typename LevelsetFunctionalType_,
      typename SpaceType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, DataType_>
    >
    class RumpfSmootherLevelsetAnalytic :
      public RumpfSmootherLevelset
      <
        DataType_,
        MemType_,
        TrafoType_,
        FunctionalType_,
        LevelsetFunctionalType_,
        SpaceType_,
        H_EvalType_
      >
    {
      public:
        /// Transformation type
        typedef TrafoType_ TrafoType;
        /// Baseclass
        typedef RumpfSmootherLevelset
        <
          DataType_,
          MemType_,
          TrafoType_,
          FunctionalType_,
          LevelsetFunctionalType_,
          SpaceType_,
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
