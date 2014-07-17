#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_SMOOTHER_LEVELSET_HPP
#define KERNEL_GEOMETRY_RUMPF_SMOOTHER_LEVELSET_HPP 1

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
     * This is the baseclass for all mesh optimisation algorithms derived from Martin Rumpf's paper
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
     * \author Jordi Paul
     *
     */
    template
    <
      typename SpaceType_,
      typename FunctionalType_,
      typename TrafoType_,
      typename DataType_,
      typename MemType_
    >
    class RumpfSmootherLevelsetBase :
      public MeshSmoother<TrafoType_, DataType_, MemType_>
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
        typedef MeshSmoother<TrafoType_, DataType_, MemType_> BaseClass;
        /// Vector types for element sizes etc.
        typedef LAFEM::DenseVector<MemType_, DataType_> VectorType;

        /// Vector with the original coordinates
        VectorType _coords_org[TrafoType_::MeshType::world_dim];
        /// Global gradient of the functional
        //
        // Later this could be a
        // LAFEM::DenseVector<MemType_, DataType_> _grad[MeshType::world_dim];
        DataType_* _grad;
        /// The functional for determining mesh quality
        FunctionalType& _functional;
        /// Index vector for identifying boundary vertices and setting boundary conditions
        int* _bdry_id;

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
        /// Functional value before mesh optimisation.
        DataType initial_functional_value;
        /// Weights for computing the functional value.
        VectorType _lambda;
        /// Size parameters for the local reference element.
        VectorType _h[MeshType::world_dim];
        /// Tolerance up to which the levelset constraint has to be fulfilled
        DataType lvlset_constraint_tol;
        /// Last computed levelset constraint;
        DataType lvlset_constraint_last;

      public:
        /// \copydoc FEAST::Geometry::RumpfSmootherBase::RumpfSmootherBase()
        explicit RumpfSmootherLevelsetBase( TrafoType& trafo_, FunctionalType& functional_, bool align_to_lvlset_, bool from_original_, bool r_adaptivity_)
          : BaseClass(trafo_),
          _grad(new DataType_[this->_world_dim*this->_nk]),
          _functional(functional_),
          _lvlset_space(trafo_),
          _lvlset_vec(_lvlset_space.get_num_dofs()),
          _lvlset_vtx_vec(this->_mesh.get_num_entities(0)),
          _align_to_lvlset(align_to_lvlset_),
          _from_original(from_original_),
          _r_adaptivity(r_adaptivity_),
          _lambda(this->_mesh.get_num_entities(ShapeType::dimension)),
          _bdry_id(new int[this->_mesh.get_num_entities(0)]),
          lvlset_constraint_tol( Math::pow(Math::eps<DataType_>(),DataType(0.5) ) ),
          lvlset_constraint_last(0)
          {

            _functional.fac_lvlset = DataType(_align_to_lvlset);

            for(Index d = 0; d < this->_world_dim; ++d)
            {
              /// Type for the boundary mesh
              typedef typename Geometry::CellSubSet<ShapeType> BoundaryType;
              /// Factory for the boundary mesh
              typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;
              // Get the boundary set
              BoundaryFactoryType boundary_factory(this->_mesh);
              BoundaryType boundary(boundary_factory);
              Geometry::TargetSet boundary_set = boundary.template get_target_set<0>();

              // Zilch bdry_id
              for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
                _bdry_id[i] = 0;

              // Set id for boundary vertices
              for(Index i(0); i < boundary.get_num_entities(0); ++i)
                _bdry_id[boundary_set[i]] = -1;

              _coords_org[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
              _h[d]= std::move(VectorType(this->_mesh.get_num_entities(ShapeType::dimension)));
              _lvlset_grad_vec[d]= std::move(VectorType(_lvlset_space.get_num_dofs()));
              _lvlset_grad_vtx_vec[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
            }

            //for(Index d = 0; d < this->_world_dim; ++d)
            //  _grad[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
          }


        virtual ~RumpfSmootherLevelsetBase()
        {
          delete[] _grad;
        };

        /// \brief Initialises member variables required for the first application of the smoother
        virtual void init() override
        {
          BaseClass::init();

          const typename MeshType::VertexSetType& vertex_set = this->_mesh.get_vertex_set();
          for(Index i(0); i < this->_nk; ++i)
          {
            for(Index d(0); d < this->_world_dim; ++d)
              _coords_org[d](i,DataType(vertex_set[i][d]));
          }

          prepare();
          compute_lambda();
          compute_h();
          initial_functional_value = compute_functional();
        }

        /// \copydoc RumpfSmootherBase::compute_functional()
        ///
        /// Variant containing the levelset penalty term
        virtual DataType compute_functional()
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

            fval += this->_lambda(cell)*_functional.compute_local_functional(x,h,lvlset_vals,this->lvlset_constraint_last);
          }

          fval += _functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last);
          return fval;
        } // compute_functional()

        /**
         * \copydoc RumpfSmootherBase::compute_functional()
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

            fval += this->_lambda(cell)*_functional.compute_local_functional(x,h,lvlset_vals, norm_A, det_A, det2_A, lvlset_penalty, this->lvlset_constraint_last);

            func_norm[cell] = this->_lambda(cell) * norm_A;
            func_det[cell] = this->_lambda(cell) * det_A;
            func_det2[cell] = this->_lambda(cell) * det2_A;
            func_lvlset[cell] =  lvlset_penalty;
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_det2_tot += func_det2[cell];
          }

          fval += _functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last);

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) << ", func_det2 = " << scientify(func_det2_tot) << ", func_lvlset = " << scientify(_functional.fac_lvlset/DataType(2)*Math::sqr(lvlset_constraint_last)) << std::endl;

          return fval;
        } // compute_functional(func_norm, func_det, func_det2, func_lvlset)
        /*
        */

        /// \copydoc RumpfSmootherBase::compute_gradient()
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
            _functional.compute_local_grad(x, h, grad_loc);
            // This does not contain the weighting yet
            grad_loc *= this->_lambda(cell);
            // Add levelset penalty term, which is not weighted with lambda
            _functional.add_levelset_penalty_grad(lvlset_vals, lvlset_grad_vals, grad_loc, lvlset_constraint_last);

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

          /*
          for(Index i(0); i < this->_mesh.get_num_entities(0); i++)
          {
            for(Index d(0); d < this->_world_dim; ++d)
            {
              if( (Math::abs(this->_coords[d](i)) <= Math::eps<DataType_>() )
                  || (Math::abs(this->_coords[d](i) - DataType(1)) <= Math::eps<DataType_>() ) )
                this->_grad[d*this->_nk + i] = DataType(0);
            }
          }
          */


        }

        /// \copydoc MeshSmoother::optimise()
        virtual void optimise() override
        {

          if(this->_from_original)
          {
            for(Index d(0); d < this->_world_dim; ++d)
              this->_coords[d].copy(_coords_org[d]);

            this->set_coords();
          }

          // DEBUG
          //auto nt = this->_trafo.get_mesh().get_num_entities(2);
          //DataType* func_norm(new DataType[nt]);
          //DataType* func_det(new DataType[nt]);
          //DataType* func_det2(new DataType[nt]);
          //DataType* func_lvlset(new DataType[nt]);
          //
          if(this->_align_to_lvlset)
          {
            // Value of the levelset constraint from the previous iteration
            DataType lvlset_constraint_previous(1);
            // The factor for the penalty term is increased by this factor in every iteration
            DataType big_factor_inc(1);

            // Quadratic penalty method for solving the constraint minimization problem as a sequence of unconstrained
            // minimisation problems
            int penalty_iterations(0);

            DataType fac_lvlset_org(this->_functional.fac_lvlset);
            while(lvlset_constraint_previous > this->lvlset_constraint_tol && _functional.fac_lvlset < DataType(1e60))
            {
              penalty_iterations++;

              //DEBUG
              //filename = "iter_pre" + stringify(penalty_iterations) +".vtk";
              //Geometry::ExportVTK<MeshType> writer(this->_trafo.get_mesh());
              //writer.add_scalar_vertex("grad_0", &this->_grad[0]);
              //writer.add_scalar_vertex("grad_1", &this->_grad[this->_trafo.get_mesh().get_num_entities(0)]);
              //writer.add_scalar_vertex("levelset", this->_lvlset_vtx_vec.elements());
              //writer.write(filename);

              ALGLIBWrapper<RumpfSmootherLevelsetBase<SpaceType_, FunctionalType_, TrafoType_, DataType_, MemType_>>::minimise_functional_cg(*this);
              // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
              this->set_coords();

              lvlset_constraint_previous = this->lvlset_constraint_last;
              // Increment for the penalty factor
              big_factor_inc *= Math::sqr(this->lvlset_constraint_last/lvlset_constraint_previous);

              // Increase penalty parameter by at least factor 5, very arbitrary
              _functional.fac_lvlset *= DataType(5)*FEAST::Math::max<DataType>(DataType(1),big_factor_inc);

              // DEBUG
              // std::cout << "Rumpf Smoother penalty iteration: " << penalty_iterations << ", last constraint: " << this->lvlset_constraint_last << ", penalty factor: " << _functional.fac_lvlset << ", fval = " << scientify(this->compute_functional()) << std::endl;
            }

            std::cout << "Rumpf Smoother penalty iterations: " << penalty_iterations << ", last constraint: " << this->lvlset_constraint_last << ", penalty factor: " << _functional.fac_lvlset << ", fval = " << scientify(this->compute_functional()) << std::endl;

            // Restore original values
            _functional.fac_lvlset = fac_lvlset_org;

          }
          else
          {
            ALGLIBWrapper<RumpfSmootherLevelsetBase<SpaceType_, FunctionalType_, TrafoType_, DataType_, MemType_>>::minimise_functional_cg(*this);
            // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
            this->set_coords();
          }

          // DEBUG clean up
          //delete func_norm;
          //delete func_det;
          //delete func_det2;
          //delete func_lvlset;

        }

        /**
         * \copydoc RumpfSmootherBase::prepare()
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

          if(this->_r_adaptivity)
          {
            this->compute_lambda();
            this->compute_h();
          }
        }


      protected:
        /**
         * \brief Computes the volume of the optimal reference for each cell and saves it to _h.
         *
         */
        virtual void compute_h()
        {
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));
          DataType sum_det(0);
          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_, MeshType::world_dim> h;
          for(Index d(0); d < this->_world_dim; ++d)
            h(d) = DataType(1);

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // DEBUG
          //DataType_ target_vol(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d = 0; d < this->_world_dim; d++)
            {
              // Get local coordinates
              for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
                x(d,j) = this->_coords[d](idx(cell,j));
            }
            sum_det+=_functional.compute_det_A(x,h);
          }

          DataType_ exponent = DataType_(1)/DataType_(this->_world_dim);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < this->_world_dim; ++d)
              this->_h[d](cell,Math::pow(this->_lambda(cell)*sum_det,exponent));

            //target_vol += DataType_(0.25)*Math::sqrt(DataType(3))*Math::sqr(_h[0](cell));

          }
          //std::cout << "compute_h target vol = " << scientify(target_vol) << std::endl;

        } // compute_h

        /// \copydoc RumpfSmootherBase::compute_lambda()
        virtual void compute_lambda()
        {
          if(this->_r_adaptivity)
          {
            compute_lambda_lvlset();
            return;
          }

          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_, MeshType::world_dim> h;

          DataType sum_lambda(0);
          DataType fac(DataType(1)/DataType(ncells));
          for(Index cell(0); cell < ncells; ++cell)
          {
            _lambda(cell, fac);
            sum_lambda+=_lambda(cell);
          }

          // Scale so that sum(lambda) = 1
          sum_lambda = DataType(1)/sum_lambda;
          for(Index k(0); k < ncells; ++k)
            _lambda(k,sum_lambda*_lambda(k));

        } // compute_lambda

        /// \copydoc RumpfSmootherBase::compute_lambda()
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
              tmp += Math::abs(this->_lvlset_vtx_vec(idx(cell,j)));

            _lambda(cell, DataType(1e-1)+Math::sqrt(tmp/DataType(1)));
            sum_lambda+=_lambda(cell);
          }

          // Scale so that sum(lambda) = 1
          sum_lambda = DataType(1)/sum_lambda;
          for(Index k(0); k < ncells; ++k)
          {
            _lambda(k,sum_lambda*_lambda(k));
          //  std::cout << "lambda(" << k << ") = " << _lambda(k) << std::endl;
          }

        }

    }; // class RumpfSmootherLevelsetBase

    /**
     * \brief Generic template for specialisations of RumpfSmootherLevelsets
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
     * \author Jordi Paul
     *
     **/
    template<typename SpaceType_, typename FunctionalType_, typename TrafoType_, typename DataType_, typename MemType_ >
    class RumpfSmootherLevelset;

    /**
     * \brief Specialisation of RumpfSmootherLevelsetBase for Q1 transformations on Hypercube<2> meshes in 2d
     *
     * For each element, there is an optimal reference cell \f[ K^* = [-h_1, h_1] \times [-h_2, h_2] \f]
     *
     * \tparam SpaceType_
     * FE space for the levelset function
     *
     * \tparam FunctionalType_
     * Functional used for determining mesh quality.
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \author Jordi Paul
     *
     */
    template<typename SpaceType_,typename FunctionalType_, typename DataType_, typename MemType_ >
    class RumpfSmootherLevelset
    <
      SpaceType_,
      FunctionalType_,
      Trafo::Standard::Mapping
      <
        Geometry::ConformalMesh
        <
          Shape::Hypercube<2>,
          Shape::Hypercube<2>::dimension,
          Shape::Hypercube<2>::dimension,
          DataType_
        >
      >,
      DataType_,
      MemType_
    > : public RumpfSmootherLevelsetBase
    <
      SpaceType_,
      FunctionalType_,
      Trafo::Standard::Mapping<Geometry::ConformalMesh
      <
        Shape::Hypercube<2>,
        Shape::Hypercube<2>::dimension,
        Shape::Hypercube<2>::dimension,
        DataType_ >
      >,
      DataType_,
      MemType_
    >
    {
      public:
        /// Our functional type
        typedef FunctionalType_ FunctionalType;
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// ShapeType
        typedef Shape::Hypercube<2> ShapeType;
        /// Mesh of said ShapeType
        typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
        /// Type for the transformation
        typedef Trafo::Standard::Mapping<MeshType> TrafoType;
        /// Who's my daddy?
        typedef RumpfSmootherLevelsetBase<SpaceType_, FunctionalType_, TrafoType, DataType_, MemType_> BaseClass;

        /// \copydoc RumpfSmootherLevelsetBase()
        explicit RumpfSmootherLevelset( TrafoType& trafo_, FunctionalType& functional_, bool align_to_lvlset_, bool from_original_, bool r_adaptivity_)
          : BaseClass(trafo_, functional_, align_to_lvlset_, from_original_, r_adaptivity_)
        {
        }

      protected:
        /**
         * \copydoc RumpfSmootherLevelsetBase::compute_h
         *
         * TODO: The components of _h are supposed to reflect the aspect ratio of the real element
         */
        virtual void compute_h()
        {
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));
          DataType vol(0);
          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          for(Index cell(0); cell < ncells; ++cell)
            vol += this->_trafo.template compute_vol<ShapeType,DataType>(cell);

          DataType_ exponent = DataType_(1)/DataType_(this->_world_dim);
          // DEBUG
          //DataType_ target_vol(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < this->_world_dim; ++d)
              this->_h[d](cell,DataType(0.5)*Math::pow(this->_lambda(cell)*vol,exponent));

            //target_vol += DataType_(4) * this->_h[0](cell)*this->_h[1](cell);
          }
          //std::cout << "compute_h target_vol = " << scientify(target_vol) << std::endl;
        } // compute_h

        /**
         * \copydoc RumpfSmootherBase::compute_lambda()
         *
         * Includes condensation near the 0 levelset
         **/
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

            this->_lambda(cell,DataType(1e-5) + Math::sqrt(tmp/DataType(1)));
            sum_lambda+=this->_lambda(cell);
          }

          // Scale so that sum(lambda) = 1
          sum_lambda = DataType(1)/sum_lambda;
          for(Index k(0); k < ncells; ++k)
          {
            this->_lambda(k,sum_lambda*this->_lambda(k));
          //  std::cout << "lambda(" << k << ") = " << _lambda(k) << std::endl;
          }

        }
    }; // class RumpfSmootherLevelset<Hypercube<2>>

    /**
     * \brief Specialisation of RumpfSmootherLevelsetBase for P1 transformations on Simplex<2> meshes in 2d
     *
     * For each element, there is an optimal reference cell \f[ K^* = h \cdot \hat{K} \f], with \f$ \hat{K} \f$
     * the simplex defined by the vertices \f$ (0,0), (0,1) \left( \frac{1}{2}, \frac{\sqrt{3}}{2} \right) \f$
     *
     * \tparam SpaceType_
     * FE space for the levelset function
     *
     * \tparam FunctionalType_
     * Functional used for determining mesh quality.
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \author Jordi Paul
     *
     */
    template<typename SpaceType_, typename FunctionalType_, typename DataType_, typename MemType_>
    class RumpfSmootherLevelset
    <
      SpaceType_,
      FunctionalType_,
      Trafo::Standard::Mapping
      <
        Geometry::ConformalMesh
        <
          Shape::Simplex<2>,
          Shape::Simplex<2>::dimension,
          Shape::Simplex<2>::dimension,
          DataType_
        >
      >, DataType_, MemType_ >:
      public RumpfSmootherLevelsetBase
      <
        SpaceType_,
        FunctionalType_,
        Trafo::Standard::Mapping<Geometry::ConformalMesh
        <
          Shape::Simplex<2>,
          Shape::Simplex<2>::dimension,
          Shape::Simplex<2>::dimension,
          DataType_
        >
      >,
      DataType_,
      MemType_
    >
    {
      public:
        /// Our functional type
        typedef FunctionalType_ FunctionalType;
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// ShapeType
        typedef Shape::Simplex<2> ShapeType;
        /// Mesh of said ShapeType
        typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
        /// Type for the transformation
        typedef Trafo::Standard::Mapping<MeshType> TrafoType;
        /// Who's my daddy?
        typedef RumpfSmootherLevelsetBase< SpaceType_, FunctionalType_, TrafoType, DataType_, MemType_ > BaseClass;

        /// \copydoc RumpfSmootherLevelsetBase()
        explicit RumpfSmootherLevelset( TrafoType& trafo_, FunctionalType& functional_, bool align_to_lvlset_, bool from_original_, bool r_adaptivity_)
          : BaseClass(trafo_, functional_, align_to_lvlset_, from_original_, r_adaptivity_)
          {
          }
    };

    /**
     * \brief Generic template for specialisations of RumpfSmootherLevelsets with analytic function
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
      typename MemType_
    >
    class RumpfSmootherLevelsetAnalytic;

    /**
     * \brief Specialisation of RumpfSmootherLevelsetAnalytic for P1 transformations in 2d
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
      typename DataType_,
      typename MemType_
    >
    class RumpfSmootherLevelsetAnalytic
    <
      AnalyticFunctionType_,
      AnalyticFunctionGrad0Type_,
      AnalyticFunctionGrad1Type_,
      SpaceType_,
      FunctionalType_,
      Trafo::Standard::Mapping
      <
        Geometry::ConformalMesh
        <
          Shape::Simplex<2>,
          Shape::Simplex<2>::dimension,
          Shape::Simplex<2>::dimension,
          DataType_
        >
      >,
      DataType_,
      MemType_
    > :
    public RumpfSmootherLevelset
    <
      SpaceType_,
      FunctionalType_,
      Trafo::Standard::Mapping
      <
        Geometry::ConformalMesh
        <
          Shape::Simplex<2>,
          Shape::Simplex<2>::dimension,
          Shape::Simplex<2>::dimension,
          DataType_
        >
      >,
      DataType_,
      MemType_
    >
    {
      public:
        /// Transformation type
        typedef typename Trafo::Standard::Mapping
        <
          Geometry::ConformalMesh
          <
            Shape::Simplex<2>,
            Shape::Simplex<2>::dimension,
            Shape::Simplex<2>::dimension,
            DataType_
          >
        > TrafoType;
        /// Baseclass
        typedef RumpfSmootherLevelset
        <
          SpaceType_,
          FunctionalType_,
          Trafo::Standard::Mapping
          <
            Geometry::ConformalMesh
            <
              Shape::Simplex<2>,
              Shape::Simplex<2>::dimension,
              Shape::Simplex<2>::dimension,
              DataType_
            >
          >,
          DataType_,
          MemType_
        > BaseClass;

        /// Analytic levelset function
        AnalyticFunctionType_& _analytic_lvlset;
        /// 1st component of its gradient
        AnalyticFunctionGrad0Type_& _analytic_lvlset_grad0;
        /// 2nd component of its gradient
        AnalyticFunctionGrad1Type_& _analytic_lvlset_grad1;
      public:
        /**
         * \copydoc RumpfSmootherLevelsetBase()
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
         * \copydoc RumpfSmootherLevelsetBase::prepare()
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
    };

    /**
     * \brief Specialisation of RumpfSmootherLevelsetAnalytic for Q1 transformations in 2d
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
      typename DataType_,
      typename MemType_
    >
    class RumpfSmootherLevelsetAnalytic
    <
      AnalyticFunctionType_,
      AnalyticFunctionGrad0Type_,
      AnalyticFunctionGrad1Type_,
      SpaceType_,
      FunctionalType_,
      Trafo::Standard::Mapping
      <
        Geometry::ConformalMesh
        <
          Shape::Hypercube<2>,
          Shape::Hypercube<2>::dimension,
          Shape::Hypercube<2>::dimension,
          DataType_
        >
      >,
      DataType_,
      MemType_
    > :
    public RumpfSmootherLevelset
    <
      SpaceType_,
      FunctionalType_,
      Trafo::Standard::Mapping
      <
        Geometry::ConformalMesh
        <
          Shape::Hypercube<2>,
          Shape::Hypercube<2>::dimension,
          Shape::Hypercube<2>::dimension,
          DataType_
        >
      >,
      DataType_,
      MemType_
    >
    {
      public:
        /// Our functional type
        typedef FunctionalType_ FunctionalType;
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// ShapeType
        typedef Shape::Hypercube<2> ShapeType;
        /// Mesh of said ShapeType
        typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
        /// Type for the transformation
        typedef Trafo::Standard::Mapping<MeshType> TrafoType;
        /// Who's my daddy?
        typedef RumpfSmootherLevelset< SpaceType_, FunctionalType_, TrafoType, DataType_, MemType_ > BaseClass;

        /// Analytic levelset function
        AnalyticFunctionType_& _analytic_lvlset;
        /// 1st component of its gradient
        AnalyticFunctionGrad0Type_& _analytic_lvlset_grad0;
        /// 2nd component of its gradient
        AnalyticFunctionGrad1Type_& _analytic_lvlset_grad1;

      public:
        /**
         * \copydoc RumpfSmootherLevelsetBase()
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
         * \copydoc RumpfSmootherBase::prepare()
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
          // Evaluate the gradient of the levelset function...
          Assembly::Interpolator::project(this->_lvlset_grad_vec[0], _analytic_lvlset_grad0, this->_lvlset_space);
          Assembly::Interpolator::project(this->_lvlset_grad_vec[1], _analytic_lvlset_grad1, this->_lvlset_space);
          // BaseClass::prepare handles projection to grid vectors etc.
          BaseClass::prepare();
        }

    }; // class RumpfSmootherLevelset
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_RUMPF_SMOOTHER_LEVELSET_HPP
