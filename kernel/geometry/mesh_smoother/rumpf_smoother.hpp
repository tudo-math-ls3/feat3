#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP
#define KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP 1

#include <kernel/geometry/mesh_smoother/mesh_smoother.hpp>
// ALGLIB includes
FEAST_DISABLE_WARNINGS
#include <thirdparty/alglib/cpp/src/optimization.h>
FEAST_RESTORE_WARNINGS

namespace FEAST
{
  namespace Geometry
  {

    template<typename RumpfSmootherType_>
    struct ALGLIBWrapper;

    /**
     * \brief Baseclass for a family of variational mesh optimisation algorithms.
     *
     * This abstract class is the baseclass for all mesh optimisation algorithms derived from Martin Rumpf's paper
     *
     * M. Rumpf: A variational approach to optimal meshes, Numerische Mathematik 72 (1996), pp. 523 - 540.
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
    template<typename FunctionalType_, typename TrafoType_, typename DataType_, typename MemType_>
    class RumpfSmootherBase :
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
        /// Who's my daddy?
        typedef MeshSmoother<TrafoType_, DataType_, MemType_> BaseClass;
        /// Vector types for element sizes etc.
        typedef LAFEM::DenseVector<MemType_, DataType_> VectorType;

        /// Global gradient of the functional
        //
        // Later this could be a
        // VectorType _grad[MeshType::world_dim];
        DataType_* _grad;
        /// The functional for determining mesh quality
        FunctionalType& _functional;
        /// Index vector for identifying boundary vertices and setting boundary conditions
        int* _bdry_id;

      protected :
        /// Functional value before mesh optimisation.
        DataType initial_functional_value;
        /// Weights for computing the functional value.
        VectorType _lambda;

        // This is public for debugging purposes
      public:
        /// Size parameters for the local reference element.
        VectorType _h[MeshType::world_dim];

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
        explicit RumpfSmootherBase(const TrafoType& trafo_, FunctionalType& functional_)
          : BaseClass(trafo_),
          _grad(new DataType_[this->_world_dim*this->_nk]),
          _functional(functional_),
          _lambda(this->_mesh.get_num_entities(ShapeType::dimension)),
          _bdry_id(new int[this->_mesh.get_num_entities(0)])
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

        for(Index d = 0; d < this->_world_dim; ++d)
          _h[d]= std::move(VectorType(this->_mesh.get_num_entities(ShapeType::dimension)));

        //for(Index d = 0; d < this->_world_dim; ++d)
        //  _grad[d]= std::move(VectorType(this->_mesh.get_num_entities(0)));
      }


        virtual ~RumpfSmootherBase()
        {
          delete[] _grad;
          delete _bdry_id;
        };

        /// \brief Initialises member variables required for the first application of the smoother
        virtual void init()
        {
          BaseClass::init();
          compute_lambda();
          compute_h();
          initial_functional_value = compute_functional();
        }

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
            for(Index d = 0; d < this->_world_dim; d++)
            {
              h(d) = this->_h[d](cell);
              for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
              {
                x(d,j) = this->_coords[d](idx(cell,j));
              }
            }
            fval += this->_lambda(cell)*_functional.compute_local_functional(x,h);
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
         * \param[in] func_det2
         * The contribution of the 1/det term for each cell
         *
         * \returns
         * The functional value
         *
         * Debug variant that saves the different contributions for each cell.
         **/
        virtual DataType compute_functional( DataType_* func_norm, DataType_* func_det, DataType_* func_det2 )
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

          DataType_ norm_A(0), det_A(0), det2_A(0);

          DataType_ func_norm_tot(0);
          DataType_ func_det_tot(0);
          DataType_ func_det2_tot(0);
          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d = 0; d < this->_world_dim; d++)
            {
              h(d) = this->_h[d](cell);
              // Get local coordinates
              for(Index j = 0; j < Shape::FaceTraits<ShapeType,0>::count; j++)
                x(d,j) = this->_coords[d](idx(cell,j));
            }

            fval += this->_lambda(cell)*_functional.compute_local_functional(x,h, norm_A, det_A, det2_A);

            func_norm[cell] = this->_lambda(cell) * norm_A;
            func_det[cell] = this->_lambda(cell) * det_A;
            func_det2[cell] = this->_lambda(cell) * det2_A;
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_det2_tot += func_det2[cell];
          }

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) << ", func_det2 = " << scientify(func_det2_tot) << std::endl;

          return fval;
        } // compute_functional

        /// \brief Computes the gradient of the functional with regard to the nodal coordinates.
        virtual void compute_gradient()
        {
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
          for(Index i(0); i < this->_world_dim*this->_nk; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < this->_world_dim; ++d)
            {
              h(d) = this->_h[d](cell);
              // Get local coordinates
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                x(d,j) = this->_coords[d](idx(cell,j));
            }

            _functional.compute_local_grad(x, h, grad_loc);

            for(Index d(0); d < this->_world_dim; ++d)
            {
              // Add local contributions to global gradient vector
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                this->_grad[d*this->_nk + idx(cell,j)] += this->_lambda(cell)*grad_loc(d,j);
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
        virtual void optimise()
        {
          ALGLIBWrapper<RumpfSmootherBase<FunctionalType_, TrafoType_, DataType_, MemType_>>::minimise_functional_cg(*this);
          // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
          this->set_coords();
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

          DataType_ target_vol(0);
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

            target_vol += DataType_(0.25)*Math::sqrt(DataType(3))*Math::sqr(_h[0](cell));

          }
          std::cout << "compute_h target vol = " << scientify(target_vol) << std::endl;
        }

        /// Computes the weights _lambda.
        virtual void compute_lambda()
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

    }; // class RumpfSmootherBase

    /**
     * \brief Generic template for specialisations of RumpfSmoothers
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
    template< typename FunctionalType_, typename TrafoType_, typename DataType_, typename MemType_ >
    class RumpfSmoother;

    /**
     * \brief Specialisation of RumpfSmoother for Q1 transformations on Hypercube<2> meshes in 2d
     *
     * For each element, there is an optimal reference cell \f[ K^* = [-h_1, h_1] \times [-h_2, h_2] \f]
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
    template<typename FunctionalType_, typename DataType_, typename MemType_ >
    class RumpfSmoother
    <
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
    > : public RumpfSmootherBase
    <
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
        typedef RumpfSmootherBase< FunctionalType_, TrafoType, DataType_, MemType_ > BaseClass;

        /// \copydoc RumpfSmootherBase()
        explicit RumpfSmoother( const TrafoType& trafo_, FunctionalType& functional_)
          : BaseClass(trafo_, functional_)
        {
        }

      protected:
        /**
         * \copydoc RumpfSmootherBase::compute_h
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
          DataType_ target_vol(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < this->_world_dim; ++d)
              this->_h[d](cell,DataType(0.5)*Math::pow(this->_lambda(cell)*vol,exponent));

            target_vol += DataType_(4) * this->_h[0](cell)*this->_h[1](cell);
          }
          std::cout << "compute_h target_vol = " << scientify(target_vol) << std::endl;
        } // compute_h

    }; // class RumpfSmoother<Hypercube<2>>

    /**
     * \brief Specialisation of RumpfSmoother for P1 transformations on Simplex<2> meshes in 2d
     *
     * For each element, there is an optimal reference cell \f[ K^* = h \cdot \hat{K} \f], with \f$ \hat{K} \f$
     * the simplex defined by the vertices \f$ (0,0), (0,1) \left( \frac{1}{2}, \frac{\sqrt{3}}{2} \right) \f$
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
    template<typename FunctionalType_, typename DataType_, typename MemType_>
    class RumpfSmoother
    <
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
      public RumpfSmootherBase
      <
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
          typedef RumpfSmootherBase< FunctionalType_, TrafoType, DataType_, MemType_ > BaseClass;

          /// \copydoc RumpfSmootherBase()
          explicit RumpfSmoother( const TrafoType& trafo_, FunctionalType& functional_)
            : BaseClass(trafo_, functional_)
          {
          }
      }; // class RumpfSmoother<Simplex<2>>

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
      static void minimise_functional_cg(RumpfSmootherType_& my_smoother)
      {

        const Index world_dim(my_smoother.get_world_dim());
        const Index num_vert(my_smoother.get_num_vert());

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
        double epsx = 0;
        alglib::ae_int_t maxits = 10000; // 1000
        alglib::mincgstate state;
        alglib::mincgreport rep;

        alglib::mincgcreate(x, state);
        alglib::mincgsetcgtype(state, -1);
        alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
        // mincgsuggeststep(state, 1.e-4);
        alglib::mincgoptimize(state, ALGLIBWrapper::functional_grad, nullptr, &my_smoother);
        alglib::mincgresults(state, x, rep);

        // Copy back the vertex coordinates
        for(Index d(0); d < world_dim; ++d)
        {
          for(Index j(0); j < num_vert; ++j)
            // Warning: Index/ae_int_t conversion
            my_smoother._coords[d](j,DataType(x[alglib::ae_int_t(d*num_vert + j)]));
        }

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

        const Index world_dim(my_smoother->get_world_dim());
        const Index num_vert(my_smoother->get_num_vert());

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
        DataType norm_grad(0);
        for(Index j(0); j < world_dim*num_vert; ++j)
        {
          grad[alglib::ae_int_t(j)] = my_smoother->_grad[j];
          norm_grad+=Math::sqr(my_smoother->_grad[j]);
        }
        norm_grad = Math::sqrt(norm_grad/DataType(world_dim*num_vert));
        //std::cout << "|| grad ||_l2 = " << scientify(norm_grad) << std::endl;

      }
    }; // struct ALGLIBWrapper

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP
