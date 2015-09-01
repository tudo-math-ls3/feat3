#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP
#define KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP 1

#include <kernel/geometry/mesh_smoother/mesh_smoother.hpp>
// ALGLIB includes
FEAST_DISABLE_WARNINGS
#include <thirdparty/alglib/cpp/src/optimization.h>
FEAST_RESTORE_WARNINGS
#include <kernel/geometry/mesh_smoother/h_evaluator.hpp>
// DEBUG
#include <iostream>
#include <kernel/geometry/export_vtk.hpp>

namespace FEAST
{
  namespace Geometry
  {

    template<typename RumpfSmootherType_>
    struct ALGLIBWrapper;

    /**
     * \brief Baseclass for a family of variational mesh optimisation algorithms.
     *
     * That is, mesh optimisation algorithms derived from Martin Rumpf's paper \cite Rum96
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
      typename DataType_,
      typename MemType_,
      typename TrafoType_,
      typename FunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, DataType_>
    >
    class RumpfSmootherBase:
      public MeshSmoother<DataType_, MemType_, TrafoType_>
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
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Who's my daddy?
        typedef MeshSmoother<DataType_, MemType_, TrafoType_> BaseClass;
        /// Vector types for element sizes etc.
        typedef LAFEM::DenseVector<MemType_, DataType_> VectorType;
        /// Since the functional contains a ShapeType, these have to be the same
        static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        "ShapeTypes of the transformation / functional have to agree" );

        /// Global gradient of the functional
        // Later this could be a
        // VectorType _grad[MeshType::world_dim];
        DataType_* _grad;
        /// The functional for determining mesh quality
        FunctionalType& _functional;
        /// Index vector for identifying boundary vertices and setting boundary conditions
        int* _bdry_id;

        // This is public for debugging purposes
      public:
        /// Functional value before mesh optimisation.
        DataType initial_functional_value;
        /// Weights for the local contributions to the global functional value.
        VectorType _mu;
        /// Weights for local mesh size
        VectorType _lambda;
        /// Size parameters for the local reference element.
        VectorType _h[MeshType::world_dim];
        // DEBUG
        //int iteration_count;

      public:
        /**
         * \brief Prints some characteristics of the RumpfSmoother object
         **/
        virtual void print()
        {
          _functional.print();
        }
        /**
         * \brief Constructor
         *
         * \param[in] trafo_
         * Reference to the underlying transformation
         *
         * \param[in] functional_
         * Reference to the functional used
         *
         */
        explicit RumpfSmootherBase(const TrafoType& trafo_, FunctionalType& functional_)
          : BaseClass(trafo_),
          _grad(new DataType_[MeshType::world_dim*trafo_.get_mesh().get_num_entities(0)]),
          _functional(functional_),
          _bdry_id(new int[trafo_.get_mesh().get_num_entities(0)]),
          _mu(trafo_.get_mesh().get_num_entities(ShapeType::dimension),DataType(1)/DataType(trafo_.get_mesh().get_num_entities(ShapeType::dimension))),
          _lambda(trafo_.get_mesh().get_num_entities(ShapeType::dimension))
          {
            /// Type for the boundary mesh
            typedef typename Geometry::MeshPart<MeshType> BoundaryType;
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
              _bdry_id[boundary_set[i]] = -MeshType::world_dim;

            for(Index d(0); d < MeshType::world_dim; ++d)
              _h[d]= std::move(VectorType(this->_mesh.get_num_entities(ShapeType::dimension)));

            //for(Index d = 0; d < MeshType::world_dim; ++d)
            //  _grad[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
          }

        /**
         * \brief Computes a quality indicator concerning the cell sizes
         *
         * In a truly optimal mesh (consisting ONLY of Rumpf reference cells of the right size), every cell's volume is
         * exaclty lambda(cell). This is especially the goal for r-adaptivity.
         * So in an optimal mesh, size(cell)/lambda(cell) = 1, so we compute the Euclidean norm of this vector, scaled
         * by the number of cells so it is independant of the refinement level. Not sure if the scaling part is
         * sensible, though.
         *
         **/
        DataType cell_size_quality()
        {
          VectorType tmp(this->_mesh.get_num_entities(ShapeType::dimension));
          DataType my_vol(0);

          for(Index cell(0); cell < this->_mesh.get_num_entities(ShapeType::dimension); ++cell)
          {
            my_vol = this->_trafo.template compute_vol<ShapeType, DataType>(cell);
            tmp(cell, Math::abs(DataType(1) - my_vol/this->_lambda(cell)));
          }

          return tmp.norm2()/Math::sqrt(DataType(this->_mesh.get_num_entities(ShapeType::dimension)));
        }

        /// \brief Destructor
        virtual ~RumpfSmootherBase()
        {
          delete[] _grad;
          delete[] _bdry_id;
        };

        /// \brief Initialises member variables required for the first application of the smoother
        virtual void init()
        {
          BaseClass::init();
          compute_lambda();
          compute_h();
          compute_mu();
          initial_functional_value = compute_functional();
        }

        /**
         * \brief Prepares the functional for evaluation.
         *
         * Needs to be called whenever any data like the mesh, the levelset function etc. changed.
         *
         **/
        virtual void prepare()
        {
        }

        virtual DataType compute_functional() = 0;
        virtual void compute_gradient() = 0;

      protected:

        /**
         * \brief Computes the volume of the optimal reference for each cell and saves it to _h.
         *
         **/
        virtual void compute_h()
        {
          H_EvalType_::compute_h(_h, this->_coords, _lambda, this->_trafo);
        }

        /// \brief Computes the weights _lambda.
        virtual void compute_lambda()
        {
          compute_lambda_current();
          //compute_lambda_uniform();
        }

        /// \brief Computes the uniformly distributed weights _lambda.
        virtual void compute_lambda_current()
        {
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_, MeshType::world_dim> h;

          DataType sum_lambda(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            _lambda(cell, this->_trafo.template compute_vol<ShapeType, DataType>(cell));
            sum_lambda+=_lambda(cell);
          }

          // Scale so that sum(lambda) = 1
          sum_lambda = DataType(1)/sum_lambda;
          for(Index k(0); k < ncells; ++k)
            _lambda(k,sum_lambda*_lambda(k));
        }

        /// \brief Computes the uniformly distributed weights _lambda.
        virtual void compute_lambda_uniform()
        {
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

        }

        /// \brief Computes the weights mu
        virtual void compute_mu()
        {
          DataType fac = DataType(1)/DataType(this->_mesh.get_num_entities(ShapeType::dimension));
          for(Index cell(0); cell < this->_mesh.get_num_entities(ShapeType::dimension); cell++)
            this->_mu(cell,fac);
        }

        /// \brief Filters the functional gradient
        virtual void _filter_grad()
        {
          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));

          // Set gradient to 0 where bdry_id == -(world_dim+1)
          // TODO: Convert this to proper filtering etc. and allow for more general BCs.
          for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
          {
            // Only filter if the id is negative
            if(this->_bdry_id[i] < 0)
            {
              // Dirichlet in all directions
              if(this->_bdry_id[i] == -(MeshType::world_dim))
              {
                for(Index d(0); d < MeshType::world_dim; ++d)
                  this->_grad[d*nvertices + i] = DataType(0);
              }
              // Check for other boundary conditions
              else
              {
                // Quick'n'dirty hack for just fixing one component
                for(Index d(0); d < MeshType::world_dim; ++d)
                {
                  if(this->_bdry_id[i] == -int(d+1))
                    this->_grad[d*nvertices + i] = DataType(0);
                }
              }
            }
          }
        } // _filter_grad

    }; // class RumpfSmootherBase

    template
    <
      typename DataType_,
      typename MemType_,
      typename TrafoType_,
      typename FunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, DataType_>
    >
    class RumpfSmoother:
      public RumpfSmootherBase<DataType_, MemType_, TrafoType_, FunctionalType_, H_EvalType_>
    {
      public :
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Who's my daddy?
        typedef RumpfSmootherBase<DataType, MemType, TrafoType, FunctionalType> BaseClass;
        /// Vector types for element sizes etc.
        typedef LAFEM::DenseVector<MemType_, DataType_> VectorType;
        /// Since the functional contains a ShapeType, these have to be the same
        //static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        //"ShapeTypes of the transformation / functional have to agree" );

      public:
        /**
         * \brief Constructor
         *
         * \param[in] trafo_
         * Reference to the underlying transformation
         *
         * \param[in] functional_
         * Reference to the functional used
         *
         */
        explicit RumpfSmoother(const TrafoType& trafo_, FunctionalType& functional_)
          : BaseClass(trafo_, functional_)
          {
          }

        /// \brief Destructor
        virtual ~RumpfSmoother()
        {
        };

        /**
         * \brief Computes the functional value on the current mesh.
         *
         * \returns
         * The functional value
         *
         **/
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
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d = 0; d < MeshType::world_dim; d++)
            {
              h(d) = this->_h[d](cell);
              for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
              {
                x(d,j) = this->_coords[d](idx(cell,j));
              }
            }
            fval += this->_mu(cell) * this->_functional.compute_local_functional(x,h);
          }

          return fval;
        } // compute_functional

        /**
         * \brief Computes the functional value on the current mesh.
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
         * \returns
         * The functional value
         *
         * Debug variant that saves the different contributions for each cell.
         **/
        virtual DataType compute_functional( DataType_* func_norm, DataType_* func_det, DataType_* func_rec_det )
        {
          DataType_ fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;

          DataType_ norm_A(0), det_A(0), rec_det_A(0);

          DataType_ func_norm_tot(0);
          DataType_ func_det_tot(0);
          DataType_ func_rec_det_tot(0);
          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d = 0; d < MeshType::world_dim; d++)
            {
              h(d) = this->_h[d](cell);
              // Get local coordinates
              for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
                x(d,j) = this->_coords[d](idx(cell,j));
            }

            // Scale local functional value with lambda
            fval += this->_mu(cell) * this->_functional.compute_local_functional(x,h, norm_A, det_A, rec_det_A);

            func_norm[cell] = this->_mu(cell) * norm_A;
            func_det[cell] = this->_mu(cell) * det_A;
            func_rec_det[cell] = this->_mu(cell) * rec_det_A;
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_rec_det_tot += func_rec_det[cell];
          }

          std::cout << "fval = " << scientify(fval) << " func_norm = " << scientify(func_norm_tot)
          << ", func_det = " << scientify(func_det_tot) << ", func_rec_det = " << scientify(func_rec_det_tot) << std::endl;

          return fval;
        } // compute_functional

        /// \brief Computes the gradient of the functional with regard to the nodal coordinates.
        virtual void compute_gradient()
        {
          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix<DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> grad_loc;

          // Clear gradient vector
          for(Index i(0); i < MeshType::world_dim*nvertices; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              h(d) = this->_h[d](cell);
              // Get local coordinates
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                x(d,j) = this->_coords[d](idx(cell,j));
            }

            this->_functional.compute_local_grad(x, h, grad_loc);

            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              // Add local contributions to global gradient vector
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                this->_grad[d*nvertices + idx(cell,j)] += this->_mu(cell)*grad_loc(d,j);
            }
          }

          this->_filter_grad();

        } // compute_gradient


        /// \copydoc MeshSmoother::optimise()
        virtual void optimise()
        {
          int total_grad_evals(0);
          int total_iterations(0);
          int termination_type(0);

          ALGLIBWrapper<RumpfSmoother<DataType_, MemType_, TrafoType_, FunctionalType_, H_EvalType_>>::minimise_functional_cg(total_grad_evals, total_iterations, termination_type,*this);
          // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
          this->set_coords();
          std::cout << total_iterations << " mincg iterations, " << total_grad_evals << " grad evals, terminationtype was " << termination_type << std::endl;
        }

        /**
         * \brief Prepares the functional for evaluation.
         *
         * Needs to be called whenever any data like the mesh, the levelset function etc. changed.
         *
         **/
        virtual void prepare()
        {
        }

    }; // class RumpfSmoother

    template
    <
      typename DataType_,
      typename MemType_,
      typename TrafoType_,
      typename FunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, DataType_>
    >
    class RumpfSmootherSplit:
      public RumpfSmoother<DataType_, MemType_, TrafoType_, FunctionalType_, H_EvalType_>
    {
      public :
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Who's my daddy?
        typedef RumpfSmoother<DataType, MemType, TrafoType, FunctionalType> BaseClass;
        /// Vector types for element sizes etc.
        typedef LAFEM::DenseVector<MemType_, DataType_> VectorType;
        /// Since the functional contains a ShapeType, these have to be the same
        //static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        //"ShapeTypes of the transformation / functional have to agree" );

        /// Coordinates of the vertices, as they get changed in the optimisation process
        VectorType _coords_old[MeshType::world_dim];
        DataType relaxation_parameter;
      public:
        /**
         * \brief Constructor
         *
         * \param[in] trafo_
         * Reference to the underlying transformation
         *
         * \param[in] functional_
         * Reference to the functional used
         *
         */
        explicit RumpfSmootherSplit(const TrafoType& trafo_, FunctionalType& functional_)
          : BaseClass(trafo_, functional_),
          relaxation_parameter(0)
          {
            for(Index d = 0; d < MeshType::world_dim; ++d)
              _coords_old[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
          }

        /// \brief Destructor
        virtual ~RumpfSmootherSplit()
        {
        };

        /**
         * \brief Computes the functional value on the current mesh.
         *
         * \returns
         * The functional value
         *
         **/
        virtual DataType compute_functional()
        {
          DataType_ fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x_old;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d = 0; d < MeshType::world_dim; d++)
            {
              h(d) = this->_h[d](cell);
              for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
              {
                x(d,j) = this->_coords[d](idx(cell,j));
                x_old(d,j) = this->_coords_old[d](idx(cell,j));
              }
            }
            fval += this->_mu(cell) * (this->_functional._fac_norm*this->_functional.compute_norm_A(x,h)
                + this->_functional._fac_det*this->_functional.compute_det_A(x_old,h)
                +this->_functional._fac_rec_det*this->_functional.compute_rec_det_A(x_old,h));
          }

          return fval;
        } // compute_functional

        /**
         * \brief Computes the functional value on the current mesh.
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
         * \returns
         * The functional value
         *
         * Debug variant that saves the different contributions for each cell.
         **/
        virtual DataType compute_functional( DataType_* func_norm, DataType_* func_det, DataType_* func_rec_det )
        {
          DataType_ fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;

          DataType_ norm_A(0), det_A(0), rec_det_A(0);

          DataType_ func_norm_tot(0);
          DataType_ func_det_tot(0);
          DataType_ func_rec_det_tot(0);
          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d = 0; d < MeshType::world_dim; d++)
            {
              h(d) = this->_h[d](cell);
              // Get local coordinates
              for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
                x(d,j) = this->_coords[d](idx(cell,j));
            }

            // Scale local functional value with lambda
            fval += this->_mu(cell) * this->_functional.compute_local_functional(x,h, norm_A, det_A, rec_det_A);

            func_norm[cell] = this->_mu(cell) * norm_A;
            func_det[cell] = this->_mu(cell) * det_A;
            func_rec_det[cell] = this->_mu(cell) * rec_det_A;
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_rec_det_tot += func_rec_det[cell];
          }

          std::cout << "fval = " << scientify(fval) << " func_norm = " << scientify(func_norm_tot)
          << ", func_det = " << scientify(func_det_tot) << ", func_rec_det = " << scientify(func_rec_det_tot) << std::endl;

          return fval;
        } // compute_functional

        /// \brief Computes the gradient of the functional with regard to the nodal coordinates.
        virtual void compute_gradient()
        {
          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix<DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> grad_loc;

          // Clear gradient vector
          for(Index i(0); i < MeshType::world_dim*nvertices; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              h(d) = this->_h[d](cell);
              // Get local coordinates
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                x(d,j) = this->_coords[d](idx(cell,j));
            }

            this->_functional.compute_grad_norm(x, h, grad_loc);

            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              // Add local contributions to global gradient vector
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                this->_grad[d*nvertices + idx(cell,j)] += this->_mu(cell)*grad_loc(d,j);
            }
          }

          this->_filter_grad();

        } // compute_gradient

        /// \copydoc MeshSmoother::optimise()
        virtual void optimise()
        {
          DataType diff(1);
          DataType tol(Math::pow(Math::eps<DataType_>(),DataType(0.75) ));

          int total_grad_evals(0);
          int total_iterations(0);
          int termination_type(0);

          // coords_new = relaxation_parameter*coords + (1 - relaxation_parameter)*coords_old
          relaxation_parameter = DataType(0.1);
          for(Index iter(0); iter < 1000; ++iter)
          {
            diff = DataType(0);
            int iterations(0);
            int grad_evals(0);
            // Save old coordinates for explicit terms
            for(Index d(0); d < MeshType::world_dim; ++d)
              this->_coords_old[d].clone(this->_coords[d]);

            ALGLIBWrapper<RumpfSmoother<DataType_, MemType_, TrafoType_, FunctionalType_, H_EvalType_>>::minimise_functional_cg(grad_evals, iterations, termination_type,*this);

            total_grad_evals += grad_evals;
            total_iterations += iterations;

            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
                diff+= Math::sqr(this->_coords_old[d](i) - this->_coords[d](i));
            }
            diff = Math::sqrt(diff/DataType(this->_mesh.get_num_entities(0)));

            // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              this->_coords_old[d].template scale(this->_coords_old[d],DataType(1)-relaxation_parameter);
              this->_coords[d].template axpy(this->_coords[d], this->_coords_old[d], relaxation_parameter);
            }

            this->set_coords();
            if(diff < DataType(tol))
            {
              std::cout << iter << " fixed point iterations, " << total_iterations << " mincg iterations, " << total_grad_evals << " grad evals, last terminationtype was " << termination_type << std::endl;
              return;
            }
          }
        }

        /**
         * \brief Prepares the functional for evaluation.
         *
         * Needs to be called whenever any data like the mesh, the levelset function etc. changed.
         *
         **/
        virtual void prepare()
        {
        }
    }; // class RumpfSmoother

    /**
     * \brief Wrapper struct for calling and passing data to ALGLIB.
     *
     * \tparam TrafoType_
     * Type of the transformation that the mesh smoother involved is defined for
     *
     * \tparam DataType_
     * Our data type
     *
     * \tparam MemType_
     * Memory architecture
     *
     * \author Jordi Paul
     *
     */
    template<typename RumpfSmootherType_>
    struct ALGLIBWrapper
    {
      /// Datatype from the RumpfSmoother
      typedef typename RumpfSmootherType_::DataType DataType;
      /**
       * \brief Minimises the functional using ALGLIB's nonlinear CG
       *
       * \param[in] my_smoother
       * Mesh smoother providing data and routines for evaluating the functional and its gradient.
       *
       *
       */
      static void minimise_functional_cg(int& grad_eval_count, int& iteration_count, int& termination_type,
      RumpfSmootherType_& my_smoother)
      {

        const Index world_dim(RumpfSmootherType_::MeshType::world_dim);
        const Index num_vert(my_smoother._mesh.get_num_entities(0));

        // Array for the coordinates for passing to ALGLIB's optimiser
        alglib::real_1d_array x;
        // Warning: ae_int_t/Index conversion
        x.setlength(alglib::ae_int_t(world_dim*num_vert));
        for(Index d(0); d < world_dim; ++d)
        {
          for(Index j(0); j < num_vert; ++j)
            // Warning: Index/ae_int_t conversion
            x[alglib::ae_int_t(d*num_vert + j)] = double(my_smoother._coords[d](j));
        }

        double epsg = 1.e-10;
        double epsf = 0.;
        double epsx = 0.;
        alglib::ae_int_t maxits = 10000; // 1000
        alglib::mincgstate state;
        alglib::mincgreport rep;

        alglib::mincgcreate(x, state);
        alglib::mincgsetcgtype(state, -1);
        alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
        // mincgsuggeststep(state, 1.e-4);
        alglib::mincgoptimize(state, ALGLIBWrapper::functional_grad, nullptr, &my_smoother);
        alglib::mincgresults(state, x, rep);

        //alglib::mincgcreatef(x, 1e-6, state);
        //alglib::mincgsetcgtype(state, -1);
        //alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
        //// mincgsuggeststep(state, 1.e-4);
        //alglib::mincgoptimize(state, ALGLIBWrapper::functional, nullptr, &my_smoother);
        //alglib::mincgresults(state, x, rep);

        //std::cout << "mincg: terminationtype " << rep.terminationtype << ", " << rep.iterationscount << " its, " << rep.nfev << " grad evals" << std::endl;

        // Warning: ae_int_t to int conversions
        iteration_count = int(rep.iterationscount);
        grad_eval_count = int(rep.nfev);
        termination_type = int(rep.terminationtype);
        if(rep.terminationtype == -8)
          throw InternalError(__func__, __FILE__, __LINE__, "Optimizer stopped with status -8");

      }

      /**
       * \brief Computes the functional value
       *
       * The actual work is done by the mesh smoother's routines. The mesh smoother is passed as ptr and then casted
       *
       * \param[in] x
       * Vertex coordinates for the mesh in its current state in the optimisation process.
       *
       * \param[out] func
       * Functional value.
       *
       * \param[out] grad
       * Gradient of the functional.
       *
       * \param[in] ptr
       * The mesh smoother.
       *
       */
      static void functional(const alglib::real_1d_array& x, double& func, void* ptr)
      {
        // Evil downcast, but we know what we are doing, right?
        RumpfSmootherType_* my_smoother = reinterpret_cast<RumpfSmootherType_*> (ptr);

        const Index world_dim(my_smoother->get_world_dim());
        const Index num_vert(my_smoother->get_num_vert());

        // Copy back the vertex coordinates, needed for computing the gradient on the modified mesh
        for(Index d(0); d < world_dim; ++d)
        {
          for(Index j(0); j < num_vert; ++j)
          {
            // Skip the bounday vertices to preserve their value
            if(my_smoother->_bdry_id[j] >=0)
              // Warning: Index/ae_int_t conversion
              my_smoother->_coords[d](j,DataType(x[alglib::ae_int_t(d*num_vert + j)]));
          }
        }

        // Let the functional initialise stuff if it needs to, like evaluating the levelset function on the current
        // mesh
        my_smoother->prepare();
        // Compute functional value
        func = double(my_smoother->compute_functional());

      }
      /**
       * \brief Computes the functional value and its gradient
       *
       * The actual work is done by the mesh smoother's routines. The mesh smoother is passed as ptr and then casted
       *
       * \param[in] x
       * Vertex coordinates for the mesh in its current state in the optimisation process.
       *
       * \param[out] func
       * Functional value.
       *
       * \param[out] grad
       * Gradient of the functional.
       *
       * \param[in] ptr
       * The mesh smoother.
       *
       */
      static void functional_grad(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)
      {
        // Evil downcast, but we know what we are doing, right?
        RumpfSmootherType_* my_smoother = reinterpret_cast<RumpfSmootherType_*> (ptr);

        const Index world_dim(RumpfSmootherType_::MeshType::world_dim);
        const Index num_vert(my_smoother->_mesh.get_num_entities(0));

        // Copy back the vertex coordinates, needed for computing the gradient on the modified mesh
        for(Index d(0); d < world_dim; ++d)
        {
          for(Index j(0); j < num_vert; ++j)
            // Warning: Index/ae_int_t conversion
            my_smoother->_coords[d](j,DataType(x[alglib::ae_int_t(d*num_vert + j)]));
        }

        // Let the functional initialise stuff if it needs to, like evaluating the levelset function on the current
        // mesh
        my_smoother->prepare();
        // Compute functional value
        func = double(my_smoother->compute_functional());
        // Compute functional gradient
        my_smoother->compute_gradient();

        // Copy to array for the gradient for passing to ALGLIB's optimiser
        //DataType norm_grad(0);
        for(Index j(0); j < world_dim*num_vert; ++j)
        {
          grad[alglib::ae_int_t(j)] = my_smoother->_grad[j];
          //norm_grad+=Math::sqr(my_smoother->_grad[j]);
        }
        //norm_grad = Math::sqrt(norm_grad/DataType(world_dim*num_vert));
        //std::cout << "|| grad ||_l2 = " << scientify(norm_grad) << std::endl;

        //std::string filename = "iter_" + stringify(my_smoother->iteration_count) + ".vtk";
        //Geometry::ExportVTK<typename RumpfSmootherType_::MeshType> writer_pre(my_smoother->_mesh);

        //writer_pre.add_scalar_cell("lambda", my_smoother->_lambda.elements() );
        //writer_pre.add_scalar_cell("h_0", my_smoother->_h[0].elements() );
        //writer_pre.add_scalar_cell("h_1", my_smoother->_h[1].elements() );
        //writer_pre.add_scalar_vertex("grad_0", &my_smoother->_grad[0]);
        //writer_pre.add_scalar_vertex("grad_1", &my_smoother->_grad[my_smoother->_mesh.get_num_entities(0)]);
        ////writer_pre.add_scalar_cell("norm", func_norm);
        ////writer_pre.add_scalar_cell("det", func_det);
        ////writer_pre.add_scalar_cell("rec_det", func_rec_det);
        ////writer_pre.add_scalar_vertex("levelset", my_smoother->_lvlset_vtx_vec.elements());
        ////writer_pre.add_scalar_cell("levelset_constraint", func_lvlset );
        ////std::cout << "Writing " << filename << std::endl;
        //writer_pre.write(filename);
        //my_smoother->iteration_count++;

      }
    }; // struct ALGLIBWrapper

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP
