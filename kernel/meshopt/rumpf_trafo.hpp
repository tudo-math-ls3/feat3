#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_TRAFO_HPP
#define KERNEL_MESHOPT_RUMPF_TRAFO_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/comm_base.hpp>

namespace FEAT
{
  namespace Meshopt
  {
    /**
     * \brief Computes quantities associated with the transformation to Rumpf reference cells
     *
     * Generic class template for RumpfTrafos.
     *
     * \tparam TrafoType_
     * Transformation type, the ShapeType is essential
     *
     * \tparam DataType_
     * Our data type
     *
     * For details on the reference cells, refer to \cite Rum96.
     *
     * \author Jordi Paul
     *
     **/
#ifndef DOXYGEN
    template<typename TrafoType_, typename DataType_>
    struct RumpfTrafo;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation
    // of specialisations in TrafoType_.
#else
    template<typename TrafoType_, typename DataType_>
    struct RumpfTrafo
    {
        /// Our data type
        typedef DataType_ DataType;
        /// Our transformation type
        typedef TrafoType_ TrafoType;
        /// Our shape type, as the Rumpf smoother needs to know what we are working with
        typedef typename TrafoType::ShapeType ShapeType;
        /// Type of the transformation's underlying mesh
        typedef typename TrafoType::MeshType MeshType;
        /// Type for a pack of local vertex coordinates
        typedef Tiny::Matrix<DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> Tx;
        /// Type for the material tensor
        typedef Tiny::Matrix<DataType_, MeshType::world_dim, MeshType::world_dim> MatTensorType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<Mem::Main, DataType, Index> ScalarVectorType;
        /// Vector type for coordinate vectors etc.
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType, Index, MeshType::world_dim> VectorType;

      /**
       * \brief Compute the transformation's determinant on a cell
       *
       * Meaning the determinant of the transformation from the Rumpf reference element onto the current element.
       *
       * \param[in] x
       * Cell coordinates
       *
       * \returns
       * The transformation's determinant
       *
       **/
      static DataType compute_det(Tx& x);
      /**
       * \brief Computes the sum of the local transformation's determinants
       *
       * This is needed to correctly calculate the optimal local scales h.
       *
       * \param[in] coords_
       * Vector of vertex coordinates
       *
       * \param[in] mesh_
       * The underlying mesh.
       *
       * \returns The sum of all local transformation's determinants
       *
       **/
      static DataType_ compute_sum_det(const VectorType& coords_, const MeshType& mesh_);

      /**
       * \brief Computes the optimal local mesh size
       *
       * The optimal local mesh size is calculated from the weights lambda_ such that the volumes of all cells
       * add up to the domain's volume, if each local cell was transformed to a Rumpf reference element with the
       * appropriate size.
       *
       * \param[in] h_
       * Vector of local mesh sizes
       *
       * \param[in] coords_
       * Vector of vertex coordinates
       *
       * \param[in] lambda_
       * Vector of element weights
       *
       * \param[in] mesh_
       * The underlying mesh.
       *
       **/
      static void compute_h(ScalarVectorType& h_, const VectorType& coords_, const ScalarVectorType& lambda_,
      const MeshType& mesh_);

      /**
       * \brief Computes the transformation's determinant's gradient on a cell
       *
       * Meaning the determinant of the transformation from the Rumpf reference element onto the current element.
       *
       * \param[in] x
       * Cell coordinates
       *
       * \param[in,out] grad
       * The gradient to compute
       *
       * \returns
       * The gradient of the transformation's determinant
       *
       **/
      static void compute_grad_det(Tx& grad, const Tx& x);

      /**
       * \brief Computes the gradient of sum of all cell's local transformation determinants
       *
       * This is the sum of all local gradients and needed to correctly calculate the gradient of the optimal local
       * scales h.
       *
       * \param[out] grad_
       * The global gradient wrt. the local vertices
       *
       * \param[in] coords_
       * Global vector of vertex coordinates
       *
       * \param[in] mesh_
       * The underlying mesh_.
       *
       */
      static void compute_grad_sum_det(VectorType& grad_, const VectorType& coords_, const MeshType& mesh_);

      /**
       * \brief Computes the local matrix describing the material behaviour
       *
       * \param[in] x
       * The coordinates of the local cell.
       *
       * \param[in] h
       * The local optimal edge length.
       *
       * \returns
       * The local material matrix.
       */
      static MatTensorType compute_mat_tensor(const Tx& x, const DataType_& h);
    }; // struct RumpfTrafo;
#endif

    /// \cond internal
    /**
     * \brief Computes local mesh size h for Rumpf smoothers
     *
     * Specialisation for Hypercube<2> meshes
     *
     * \tparam DataType_
     * Our data type
     *
     **/
    template<typename DataType_>
    struct RumpfTrafo<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>,2,2,DataType_>>,DataType_>
    {
      public:
        /// Our transformation type
        typedef Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>,2,2,DataType_>> TrafoType;
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type, as the Rumpf smoother needs to know what we are working with
        typedef typename TrafoType::ShapeType ShapeType;
        /// Type of the transformation's underlying mesh
        typedef typename TrafoType::MeshType MeshType;
        /// Type for a pack of local vertex coordinates
        typedef Tiny::Matrix<DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> Tx;
        /// Type for the material tensor
        typedef Tiny::Matrix<DataType_, MeshType::world_dim, MeshType::world_dim> MatTensorType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<Mem::Main, DataType, Index> ScalarVectorType;
        /// Vector type for coordinate vectors etc.
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType, Index, MeshType::world_dim> VectorType;

      private:
        /**
         * \brief Compute the transformation's determinant on a cell
         */
        static DataType compute_det(Tx& x)
        {
          return (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(0,1) * x(1,0) + x(0,1) * x(2,0)
              + x(1,0) * x(3,1) - x(1,1) * x(3,0) - x(2,0) * x(3,1) + x(2,1) * x(3,0)) / DataType(2);
        }

        /**
         * \brief Computes the transformation's determinant's gradient on a cell
         */
        static void compute_grad_det(Tx& grad, const Tx& x)
        {
          grad(0,0) = ( x(1,1) - x(2,1)) / DataType(2);
          grad(0,1) = (-x(1,0) + x(2,0)) / DataType(2);
          grad(1,0) = (-x(0,1) + x(3,1)) / DataType(2);
          grad(1,1) = ( x(0,0) - x(3,0)) / DataType(2);
          grad(2,0) = ( x(0,1) - x(3,1)) / DataType(2);
          grad(2,1) = (-x(0,0) + x(3,0)) / DataType(2);
          grad(3,0) = (-x(1,1) + x(2,1)) / DataType(2);
          grad(3,1) = ( x(1,0) - x(2,0)) / DataType(2);
        }

      public:
        /**
         * \brief Computes the optimal local mesh size
         */
        static void compute_h(ScalarVectorType& h_, const VectorType& coords_, const ScalarVectorType& lambda_,
        const MeshType& mesh_)
        {
          DataType_ sum_det = compute_sum_det(coords_, mesh_);
          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);

          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // For hypercubes, h is half the refence cell's edge length
            h_(cell, DataType(0.5)*Math::pow(lambda_(cell)*sum_det,exponent));
          }
        }

        /**
         * \brief Computes the sum of the local transformation's determinants
         **/
        static DataType_ compute_sum_det(const VectorType& coords_, const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          Tiny::Matrix <DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;

          // Index set for local/global numbering
          const auto& idx = mesh_.template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = coords_(idx(cell,Index(j)));

            sum_det += compute_det(x);
          }

#ifdef FEAT_HAVE_MPI
          Util::Comm::allreduce(&sum_det, &sum_det, 1, Util::CommOperationSum());
#endif
          return sum_det;
        }

        /**
         * \brief Computes the gradient of sum of the local transformation's determinants
         **/
        static void compute_grad_sum_det(VectorType& grad_, const VectorType& coords_, const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          Tx x(DataType_(0));

          // Index set for local/global numbering
          const auto& idx = mesh_.template get_index_set<ShapeType::dimension,0>();

          Tx local_grad(DataType_(0));

          grad_.format();

          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = coords_(idx(cell,Index(j)));

            compute_grad_det(local_grad, x);

            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              grad_(i, grad_(i) + local_grad[j]);
            }

          }
        } // compute_grad_sum_det

      /**
       * \brief Computes the local matrix describing the material behaviour
       */
        static MatTensorType compute_mat_tensor(const Tx& DOXY(x), const DataType_& h)
        {
          MatTensorType mat_tensor(DataType_(0));
          DataType scal(DataType(1)/h);

          for(int d(0); d < MeshType::world_dim; ++d)
          {
            mat_tensor(d,d) = scal;
          }

          return mat_tensor;
        }

    }; // RumpfTrafo<Hypercube<2>>

    /**
     * \brief Computes local mesh size h for Rumpf smoothers
     *
     * Specialisation for Hypercube<3> meshes
     *
     * \tparam DataType_
     * Our data type
     *
     **/
    template<typename DataType_>
    struct RumpfTrafo<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3>,3,3,DataType_>>,DataType_>
    {
      public:
        /// Our transformation type
        typedef Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3>,3,3,DataType_>> TrafoType;
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type, as the Rumpf smoother needs to know what we are working with
        typedef typename TrafoType::ShapeType ShapeType;
        /// Type of the transformation's underlying mesh
        typedef typename TrafoType::MeshType MeshType;
        /// Type for a pack of local vertex coordinates
        typedef Tiny::Matrix<DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> Tx;
        /// Type for the material tensor
        typedef Tiny::Matrix<DataType_, MeshType::world_dim, MeshType::world_dim> MatTensorType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<Mem::Main, DataType, Index> ScalarVectorType;
        /// Vector type for coordinate vectors etc.
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType, Index, MeshType::world_dim> VectorType;

      private:
        /**
         * \brief Compute the transformation's determinant on a cell
         */
        static DataType compute_det(Tx& x)
        {
          return DataType_(12)*(-x(0,0) * x(1,1) * x(2,2) - x(0,0) * x(1,1) * x(3,2) + x(0,0) * x(1,1) * x(4,2) + x(0,0) * x(1,1) * x(5,2) + x(0,0) * x(1,2) * x(2,1) + x(0,0) * x(1,2) * x(3,1) - x(0,0) * x(1,2) * x(4,1) - x(0,0) * x(1,2) * x(5,1) + x(0,0) * x(2,1) * x(3,2) - x(0,0) * x(2,1) * x(4,2) - x(0,0) * x(2,1) * x(6,2) - x(0,0) * x(2,2) * x(3,1) + x(0,0) * x(2,2) * x(4,1) + x(0,0) * x(2,2) * x(6,1) - x(0,0) * x(4,1) * x(5,2) + x(0,0) * x(4,1) * x(6,2) + x(0,0) * x(4,2) * x(5,1) - x(0,0) * x(4,2) * x(6,1) + x(0,1) * x(1,0) * x(2,2) + x(0,1) * x(1,0) * x(3,2) - x(0,1) * x(1,0) * x(4,2) - x(0,1) * x(1,0) * x(5,2) - x(0,1) * x(1,2) * x(2,0) - x(0,1) * x(1,2) * x(3,0) + x(0,1) * x(1,2) * x(4,0) + x(0,1) * x(1,2) * x(5,0) - x(0,1) * x(2,0) * x(3,2) + x(0,1) * x(2,0) * x(4,2) + x(0,1) * x(2,0) * x(6,2) + x(0,1) * x(2,2) * x(3,0) - x(0,1) * x(2,2) * x(4,0) - x(0,1) * x(2,2) * x(6,0) + x(0,1) * x(4,0) * x(5,2) - x(0,1) * x(4,0) * x(6,2) - x(0,1) * x(4,2) * x(5,0) + x(0,1) * x(4,2) * x(6,0) - x(0,2) * x(1,0) * x(2,1) - x(0,2) * x(1,0) * x(3,1) + x(0,2) * x(1,0) * x(4,1) + x(0,2) * x(1,0) * x(5,1) + x(0,2) * x(1,1) * x(2,0) + x(0,2) * x(1,1) * x(3,0) - x(0,2) * x(1,1) * x(4,0) - x(0,2) * x(1,1) * x(5,0) + x(0,2) * x(2,0) * x(3,1) - x(0,2) * x(2,0) * x(4,1) - x(0,2) * x(2,0) * x(6,1) - x(0,2) * x(2,1) * x(3,0) + x(0,2) * x(2,1) * x(4,0) + x(0,2) * x(2,1) * x(6,0) - x(0,2) * x(4,0) * x(5,1) + x(0,2) * x(4,0) * x(6,1) + x(0,2) * x(4,1) * x(5,0) - x(0,2) * x(4,1) * x(6,0) + x(1,0) * x(2,1) * x(3,2) - x(1,0) * x(2,2) * x(3,1) + x(1,0) * x(3,1) * x(5,2) + x(1,0) * x(3,1) * x(7,2) - x(1,0) * x(3,2) * x(5,1) - x(1,0) * x(3,2) * x(7,1) - x(1,0) * x(4,1) * x(5,2) + x(1,0) * x(4,2) * x(5,1) - x(1,0) * x(5,1) * x(7,2) + x(1,0) * x(5,2) * x(7,1) - x(1,1) * x(2,0) * x(3,2) + x(1,1) * x(2,2) * x(3,0) - x(1,1) * x(3,0) * x(5,2) - x(1,1) * x(3,0) * x(7,2) + x(1,1) * x(3,2) * x(5,0) + x(1,1) * x(3,2) * x(7,0) + x(1,1) * x(4,0) * x(5,2) - x(1,1) * x(4,2) * x(5,0) + x(1,1) * x(5,0) * x(7,2) - x(1,1) * x(5,2) * x(7,0) + x(1,2) * x(2,0) * x(3,1) - x(1,2) * x(2,1) * x(3,0) + x(1,2) * x(3,0) * x(5,1) + x(1,2) * x(3,0) * x(7,1) - x(1,2) * x(3,1) * x(5,0) - x(1,2) * x(3,1) * x(7,0) - x(1,2) * x(4,0) * x(5,1) + x(1,2) * x(4,1) * x(5,0) - x(1,2) * x(5,0) * x(7,1) + x(1,2) * x(5,1) * x(7,0) - x(2,0) * x(3,1) * x(6,2) - x(2,0) * x(3,1) * x(7,2) + x(2,0) * x(3,2) * x(6,1) + x(2,0) * x(3,2) * x(7,1) + x(2,0) * x(4,1) * x(6,2) - x(2,0) * x(4,2) * x(6,1) + x(2,0) * x(6,1) * x(7,2) - x(2,0) * x(6,2) * x(7,1) + x(2,1) * x(3,0) * x(6,2) + x(2,1) * x(3,0) * x(7,2) - x(2,1) * x(3,2) * x(6,0) - x(2,1) * x(3,2) * x(7,0) - x(2,1) * x(4,0) * x(6,2) + x(2,1) * x(4,2) * x(6,0) - x(2,1) * x(6,0) * x(7,2) + x(2,1) * x(6,2) * x(7,0) - x(2,2) * x(3,0) * x(6,1) - x(2,2) * x(3,0) * x(7,1) + x(2,2) * x(3,1) * x(6,0) + x(2,2) * x(3,1) * x(7,0) + x(2,2) * x(4,0) * x(6,1) - x(2,2) * x(4,1) * x(6,0) + x(2,2) * x(6,0) * x(7,1) - x(2,2) * x(6,1) * x(7,0) - x(3,0) * x(5,1) * x(7,2) + x(3,0) * x(5,2) * x(7,1) + x(3,0) * x(6,1) * x(7,2) - x(3,0) * x(6,2) * x(7,1) + x(3,1) * x(5,0) * x(7,2) - x(3,1) * x(5,2) * x(7,0) - x(3,1) * x(6,0) * x(7,2) + x(3,1) * x(6,2) * x(7,0) - x(3,2) * x(5,0) * x(7,1) + x(3,2) * x(5,1) * x(7,0) + x(3,2) * x(6,0) * x(7,1) - x(3,2) * x(6,1) * x(7,0) + x(4,0) * x(5,1) * x(6,2) + x(4,0) * x(5,1) * x(7,2) - x(4,0) * x(5,2) * x(6,1) - x(4,0) * x(5,2) * x(7,1) - x(4,0) * x(6,1) * x(7,2) + x(4,0) * x(6,2) * x(7,1) - x(4,1) * x(5,0) * x(6,2) - x(4,1) * x(5,0) * x(7,2) + x(4,1) * x(5,2) * x(6,0) + x(4,1) * x(5,2) * x(7,0) + x(4,1) * x(6,0) * x(7,2) - x(4,1) * x(6,2) * x(7,0) + x(4,2) * x(5,0) * x(6,1) + x(4,2) * x(5,0) * x(7,1) - x(4,2) * x(5,1) * x(6,0) - x(4,2) * x(5,1) * x(7,0) - x(4,2) * x(6,0) * x(7,1) + x(4,2) * x(6,1) * x(7,0) - x(5,0) * x(6,1) * x(7,2) + x(5,0) * x(6,2) * x(7,1) + x(5,1) * x(6,0) * x(7,2) - x(5,1) * x(6,2) * x(7,0) - x(5,2) * x(6,0) * x(7,1) + x(5,2) * x(6,1) * x(7,0));
        }

        /**
         * \brief Computes the transformation's determinant's gradient on a cell
         */
        static void compute_grad_det(Tx& grad, const Tx& x)
        {
          grad(0,0) = -x(1,1) * x(2,2) - x(1,1) * x(3,2) + x(1,1) * x(4,2) + x(1,1) * x(5,2) + x(1,2) * x(2,1) + x(1,2) * x(3,1) - x(1,2) * x(4,1) - x(1,2) * x(5,1) + x(2,1) * x(3,2) - x(2,1) * x(4,2) - x(2,1) * x(6,2) - x(2,2) * x(3,1) + x(2,2) * x(4,1) + x(2,2) * x(6,1) - x(4,1) * x(5,2) + x(4,1) * x(6,2) + x(4,2) * x(5,1) - x(4,2) * x(6,1);
          grad(0,1) = x(1,0) * x(2,2) + x(1,0) * x(3,2) - x(1,0) * x(4,2) - x(1,0) * x(5,2) - x(1,2) * x(2,0) - x(1,2) * x(3,0) + x(1,2) * x(4,0) + x(1,2) * x(5,0) - x(2,0) * x(3,2) + x(2,0) * x(4,2) + x(2,0) * x(6,2) + x(2,2) * x(3,0) - x(2,2) * x(4,0) - x(2,2) * x(6,0) + x(4,0) * x(5,2) - x(4,0) * x(6,2) - x(4,2) * x(5,0) + x(4,2) * x(6,0);
          grad(0,2) = -x(1,0) * x(2,1) - x(1,0) * x(3,1) + x(1,0) * x(4,1) + x(1,0) * x(5,1) + x(1,1) * x(2,0) + x(1,1) * x(3,0) - x(1,1) * x(4,0) - x(1,1) * x(5,0) + x(2,0) * x(3,1) - x(2,0) * x(4,1) - x(2,0) * x(6,1) - x(2,1) * x(3,0) + x(2,1) * x(4,0) + x(2,1) * x(6,0) - x(4,0) * x(5,1) + x(4,0) * x(6,1) + x(4,1) * x(5,0) - x(4,1) * x(6,0);

          grad(1,0) = x(0,1) * x(2,2) + x(0,1) * x(3,2) - x(0,1) * x(4,2) - x(0,1) * x(5,2) - x(0,2) * x(2,1) - x(0,2) * x(3,1) + x(0,2) * x(4,1) + x(0,2) * x(5,1) + x(2,1) * x(3,2) - x(2,2) * x(3,1) + x(3,1) * x(5,2) + x(3,1) * x(7,2) - x(3,2) * x(5,1) - x(3,2) * x(7,1) - x(4,1) * x(5,2) + x(4,2) * x(5,1) - x(5,1) * x(7,2) + x(5,2) * x(7,1);
          grad(1,1) = -x(0,0) * x(2,2) - x(0,0) * x(3,2) + x(0,0) * x(4,2) + x(0,0) * x(5,2) + x(0,2) * x(2,0) + x(0,2) * x(3,0) - x(0,2) * x(4,0) - x(0,2) * x(5,0) - x(2,0) * x(3,2) + x(2,2) * x(3,0) - x(3,0) * x(5,2) - x(3,0) * x(7,2) + x(3,2) * x(5,0) + x(3,2) * x(7,0) + x(4,0) * x(5,2) - x(4,2) * x(5,0) + x(5,0) * x(7,2) - x(5,2) * x(7,0);
          grad(1,2) = x(0,0) * x(2,1) + x(0,0) * x(3,1) - x(0,0) * x(4,1) - x(0,0) * x(5,1) - x(0,1) * x(2,0) - x(0,1) * x(3,0) + x(0,1) * x(4,0) + x(0,1) * x(5,0) + x(2,0) * x(3,1) - x(2,1) * x(3,0) + x(3,0) * x(5,1) + x(3,0) * x(7,1) - x(3,1) * x(5,0) - x(3,1) * x(7,0) - x(4,0) * x(5,1) + x(4,1) * x(5,0) - x(5,0) * x(7,1) + x(5,1) * x(7,0);

          grad(2,0) = -x(0,1) * x(1,2) - x(0,1) * x(3,2) + x(0,1) * x(4,2) + x(0,1) * x(6,2) + x(0,2) * x(1,1) + x(0,2) * x(3,1) - x(0,2) * x(4,1) - x(0,2) * x(6,1) - x(1,1) * x(3,2) + x(1,2) * x(3,1) - x(3,1) * x(6,2) - x(3,1) * x(7,2) + x(3,2) * x(6,1) + x(3,2) * x(7,1) + x(4,1) * x(6,2) - x(4,2) * x(6,1) + x(6,1) * x(7,2) - x(6,2) * x(7,1);
          grad(2,1) = x(0,0) * x(1,2) + x(0,0) * x(3,2) - x(0,0) * x(4,2) - x(0,0) * x(6,2) - x(0,2) * x(1,0) - x(0,2) * x(3,0) + x(0,2) * x(4,0) + x(0,2) * x(6,0) + x(1,0) * x(3,2) - x(1,2) * x(3,0) + x(3,0) * x(6,2) + x(3,0) * x(7,2) - x(3,2) * x(6,0) - x(3,2) * x(7,0) - x(4,0) * x(6,2) + x(4,2) * x(6,0) - x(6,0) * x(7,2) + x(6,2) * x(7,0);
          grad(2,2) = -x(0,0) * x(1,1) - x(0,0) * x(3,1) + x(0,0) * x(4,1) + x(0,0) * x(6,1) + x(0,1) * x(1,0) + x(0,1) * x(3,0) - x(0,1) * x(4,0) - x(0,1) * x(6,0) - x(1,0) * x(3,1) + x(1,1) * x(3,0) - x(3,0) * x(6,1) - x(3,0) * x(7,1) + x(3,1) * x(6,0) + x(3,1) * x(7,0) + x(4,0) * x(6,1) - x(4,1) * x(6,0) + x(6,0) * x(7,1) - x(6,1) * x(7,0);

          grad(3,0) = -x(0,1) * x(1,2) + x(0,1) * x(2,2) + x(0,2) * x(1,1) - x(0,2) * x(2,1) + x(1,1) * x(2,2) - x(1,1) * x(5,2) - x(1,1) * x(7,2) - x(1,2) * x(2,1) + x(1,2) * x(5,1) + x(1,2) * x(7,1) + x(2,1) * x(6,2) + x(2,1) * x(7,2) - x(2,2) * x(6,1) - x(2,2) * x(7,1) - x(5,1) * x(7,2) + x(5,2) * x(7,1) + x(6,1) * x(7,2) - x(6,2) * x(7,1);
          grad(3,1) = x(0,0) * x(1,2) - x(0,0) * x(2,2) - x(0,2) * x(1,0) + x(0,2) * x(2,0) - x(1,0) * x(2,2) + x(1,0) * x(5,2) + x(1,0) * x(7,2) + x(1,2) * x(2,0) - x(1,2) * x(5,0) - x(1,2) * x(7,0) - x(2,0) * x(6,2) - x(2,0) * x(7,2) + x(2,2) * x(6,0) + x(2,2) * x(7,0) + x(5,0) * x(7,2) - x(5,2) * x(7,0) - x(6,0) * x(7,2) + x(6,2) * x(7,0);
          grad(3,2) = -x(0,0) * x(1,1) + x(0,0) * x(2,1) + x(0,1) * x(1,0) - x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,0) * x(5,1) - x(1,0) * x(7,1) - x(1,1) * x(2,0) + x(1,1) * x(5,0) + x(1,1) * x(7,0) + x(2,0) * x(6,1) + x(2,0) * x(7,1) - x(2,1) * x(6,0) - x(2,1) * x(7,0) - x(5,0) * x(7,1) + x(5,1) * x(7,0) + x(6,0) * x(7,1) - x(6,1) * x(7,0);

          grad(4,0) = x(0,1) * x(1,2) - x(0,1) * x(2,2) + x(0,1) * x(5,2) - x(0,1) * x(6,2) - x(0,2) * x(1,1) + x(0,2) * x(2,1) - x(0,2) * x(5,1) + x(0,2) * x(6,1) + x(1,1) * x(5,2) - x(1,2) * x(5,1) - x(2,1) * x(6,2) + x(2,2) * x(6,1) + x(5,1) * x(6,2) + x(5,1) * x(7,2) - x(5,2) * x(6,1) - x(5,2) * x(7,1) - x(6,1) * x(7,2) + x(6,2) * x(7,1);
          grad(4,1) = -x(0,0) * x(1,2) + x(0,0) * x(2,2) - x(0,0) * x(5,2) + x(0,0) * x(6,2) + x(0,2) * x(1,0) - x(0,2) * x(2,0) + x(0,2) * x(5,0) - x(0,2) * x(6,0) - x(1,0) * x(5,2) + x(1,2) * x(5,0) + x(2,0) * x(6,2) - x(2,2) * x(6,0) - x(5,0) * x(6,2) - x(5,0) * x(7,2) + x(5,2) * x(6,0) + x(5,2) * x(7,0) + x(6,0) * x(7,2) - x(6,2) * x(7,0);
          grad(4,2) = x(0,0) * x(1,1) - x(0,0) * x(2,1) + x(0,0) * x(5,1) - x(0,0) * x(6,1) - x(0,1) * x(1,0) + x(0,1) * x(2,0) - x(0,1) * x(5,0) + x(0,1) * x(6,0) + x(1,0) * x(5,1) - x(1,1) * x(5,0) - x(2,0) * x(6,1) + x(2,1) * x(6,0) + x(5,0) * x(6,1) + x(5,0) * x(7,1) - x(5,1) * x(6,0) - x(5,1) * x(7,0) - x(6,0) * x(7,1) + x(6,1) * x(7,0);

          grad(5,0) = x(0,1) * x(1,2) - x(0,1) * x(4,2) - x(0,2) * x(1,1) + x(0,2) * x(4,1) + x(1,1) * x(3,2) - x(1,1) * x(4,2) + x(1,1) * x(7,2) - x(1,2) * x(3,1) + x(1,2) * x(4,1) - x(1,2) * x(7,1) + x(3,1) * x(7,2) - x(3,2) * x(7,1) - x(4,1) * x(6,2) - x(4,1) * x(7,2) + x(4,2) * x(6,1) + x(4,2) * x(7,1) - x(6,1) * x(7,2) + x(6,2) * x(7,1);
          grad(5,1) = -x(0,0) * x(1,2) + x(0,0) * x(4,2) + x(0,2) * x(1,0) - x(0,2) * x(4,0) - x(1,0) * x(3,2) + x(1,0) * x(4,2) - x(1,0) * x(7,2) + x(1,2) * x(3,0) - x(1,2) * x(4,0) + x(1,2) * x(7,0) - x(3,0) * x(7,2) + x(3,2) * x(7,0) + x(4,0) * x(6,2) + x(4,0) * x(7,2) - x(4,2) * x(6,0) - x(4,2) * x(7,0) + x(6,0) * x(7,2) - x(6,2) * x(7,0);
          grad(5,2) = x(0,0) * x(1,1) - x(0,0) * x(4,1) - x(0,1) * x(1,0) + x(0,1) * x(4,0) + x(1,0) * x(3,1) - x(1,0) * x(4,1) + x(1,0) * x(7,1) - x(1,1) * x(3,0) + x(1,1) * x(4,0) - x(1,1) * x(7,0) + x(3,0) * x(7,1) - x(3,1) * x(7,0) - x(4,0) * x(6,1) - x(4,0) * x(7,1) + x(4,1) * x(6,0) + x(4,1) * x(7,0) - x(6,0) * x(7,1) + x(6,1) * x(7,0);

          grad(6,0) = -x(0,1) * x(2,2) + x(0,1) * x(4,2) + x(0,2) * x(2,1) - x(0,2) * x(4,1) - x(2,1) * x(3,2) + x(2,1) * x(4,2) - x(2,1) * x(7,2) + x(2,2) * x(3,1) - x(2,2) * x(4,1) + x(2,2) * x(7,1) - x(3,1) * x(7,2) + x(3,2) * x(7,1) + x(4,1) * x(5,2) + x(4,1) * x(7,2) - x(4,2) * x(5,1) - x(4,2) * x(7,1) + x(5,1) * x(7,2) - x(5,2) * x(7,1);
          grad(6,1) = x(0,0) * x(2,2) - x(0,0) * x(4,2) - x(0,2) * x(2,0) + x(0,2) * x(4,0) + x(2,0) * x(3,2) - x(2,0) * x(4,2) + x(2,0) * x(7,2) - x(2,2) * x(3,0) + x(2,2) * x(4,0) - x(2,2) * x(7,0) + x(3,0) * x(7,2) - x(3,2) * x(7,0) - x(4,0) * x(5,2) - x(4,0) * x(7,2) + x(4,2) * x(5,0) + x(4,2) * x(7,0) - x(5,0) * x(7,2) + x(5,2) * x(7,0);
          grad(6,2) = -x(0,0) * x(2,1) + x(0,0) * x(4,1) + x(0,1) * x(2,0) - x(0,1) * x(4,0) - x(2,0) * x(3,1) + x(2,0) * x(4,1) - x(2,0) * x(7,1) + x(2,1) * x(3,0) - x(2,1) * x(4,0) + x(2,1) * x(7,0) - x(3,0) * x(7,1) + x(3,1) * x(7,0) + x(4,0) * x(5,1) + x(4,0) * x(7,1) - x(4,1) * x(5,0) - x(4,1) * x(7,0) + x(5,0) * x(7,1) - x(5,1) * x(7,0);

          grad(7,0) = x(1,1) * x(3,2) - x(1,1) * x(5,2) - x(1,2) * x(3,1) + x(1,2) * x(5,1) - x(2,1) * x(3,2) + x(2,1) * x(6,2) + x(2,2) * x(3,1) - x(2,2) * x(6,1) - x(3,1) * x(5,2) + x(3,1) * x(6,2) + x(3,2) * x(5,1) - x(3,2) * x(6,1) + x(4,1) * x(5,2) - x(4,1) * x(6,2) - x(4,2) * x(5,1) + x(4,2) * x(6,1) - x(5,1) * x(6,2) + x(5,2) * x(6,1);
          grad(7,1) = -x(1,0) * x(3,2) + x(1,0) * x(5,2) + x(1,2) * x(3,0) - x(1,2) * x(5,0) + x(2,0) * x(3,2) - x(2,0) * x(6,2) - x(2,2) * x(3,0) + x(2,2) * x(6,0) + x(3,0) * x(5,2) - x(3,0) * x(6,2) - x(3,2) * x(5,0) + x(3,2) * x(6,0) - x(4,0) * x(5,2) + x(4,0) * x(6,2) + x(4,2) * x(5,0) - x(4,2) * x(6,0) + x(5,0) * x(6,2) - x(5,2) * x(6,0);
          grad(7,2) = x(1,0) * x(3,1) - x(1,0) * x(5,1) - x(1,1) * x(3,0) + x(1,1) * x(5,0) - x(2,0) * x(3,1) + x(2,0) * x(6,1) + x(2,1) * x(3,0) - x(2,1) * x(6,0) - x(3,0) * x(5,1) + x(3,0) * x(6,1) + x(3,1) * x(5,0) - x(3,1) * x(6,0) + x(4,0) * x(5,1) - x(4,0) * x(6,1) - x(4,1) * x(5,0) + x(4,1) * x(6,0) - x(5,0) * x(6,1) + x(5,1) * x(6,0);

          // Important: Scale the matrix with the missing factor 1/12

          grad *= (DataType(1)/DataType(12));

        }

      public:
        /**
         * \brief Computes the optimal local mesh size
         */
        static void compute_h(ScalarVectorType& h_, const VectorType& coords_, const ScalarVectorType& lambda_,
        const MeshType& mesh_)
        {
          DataType_ sum_det = compute_sum_det(coords_, mesh_);
          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);

          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // For hypercubes, h is half the refence cell's edge length
            h_(cell, DataType(0.5)*Math::pow(lambda_(cell)*sum_det,exponent));
          }
        }

        /**
         * \brief Computes the sum of the local transformation's determinants
         **/
        static DataType_ compute_sum_det(const VectorType& coords_, const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          Tiny::Matrix <DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;

          // Index set for local/global numbering
          const auto& idx = mesh_.template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = coords_(idx(cell,Index(j)));

            sum_det += compute_det(x);
          }

#ifdef FEAT_HAVE_MPI
          Util::Comm::allreduce(&sum_det, &sum_det, 1, Util::CommOperationSum());
#endif
          return sum_det;
        }

        /**
         * \brief Computes the gradient of sum of the local transformation's determinants
         **/
        static void compute_grad_sum_det(VectorType& grad_, const VectorType& coords_, const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          Tx x(DataType_(0));

          // Index set for local/global numbering
          const auto& idx = mesh_.template get_index_set<ShapeType::dimension,0>();

          Tx local_grad(DataType_(0));

          grad_.format();

          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = coords_(idx(cell,Index(j)));

            compute_grad_det(local_grad, x);

            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              grad_(i, grad_(i) + local_grad[j]);
            }

          }
        } // compute_grad_sum_det

      /**
       * \brief Computes the local matrix describing the material behaviour
       */
        static MatTensorType compute_mat_tensor(const Tx& DOXY(x), const DataType_& h)
        {
          MatTensorType mat_tensor(DataType_(0));
          DataType scal(DataType(1)/h);

          for(int d(0); d < MeshType::world_dim; ++d)
          {
            mat_tensor(d,d) = scal;
          }

          return mat_tensor;
        }

    }; // RumpfTrafo<Hypercube<3>>

    /**
     * \brief Computes local mesh size h for Rumpf smoothers
     *
     * Specialisation for Simplex<2> meshes
     *
     * \tparam DataType_
     * Our data type
     *
     **/
    template<typename DataType_>
    struct RumpfTrafo<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>,2,2,DataType_>>,DataType_>
    {
      public:
        /// Our transformation type
        typedef Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>,2,2,DataType_>> TrafoType;
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type, as the Rumpf smoother needs to know what we are working with
        typedef typename TrafoType::ShapeType ShapeType;
        /// Type of the transformation's underlying mesh
        typedef typename TrafoType::MeshType MeshType;
        /// Type for a pack of local vertex coordinates
        typedef Tiny::Matrix <DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> Tx;
        /// Type for the material tensor
        typedef Tiny::Matrix<DataType_, MeshType::world_dim, MeshType::world_dim> MatTensorType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<Mem::Main, DataType, Index> ScalarVectorType;
        /// Vector type for coordinate vectors etc.
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType, Index, MeshType::world_dim> VectorType;

      private:
        /**
         * \brief Compute the transformation's determinant on a cell
         **/
        static DataType compute_det(Tx& x)
        {
          return DataType(2) / Math::sqrt(DataType(3)) *
            (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(2,0) * x(0,1) + x(1,0) * x(2,1) - x(2,0) * x(1,1));
        }

        /**
         * \brief Computes the transformation's determinant's gradient on a cell
         */
        static inline void compute_grad_det(Tx& grad, const Tx& x)
        {
          const DataType fac(DataType(2) / Math::sqrt(DataType(3)));

          grad(0,0) = fac * ( x(1,1) - x(2,1));
          grad(0,1) = fac * (-x(1,0) + x(2,0));
          grad(1,0) = fac * (-x(0,1) + x(2,1));
          grad(1,1) = fac * ( x(0,0) - x(2,0));
          grad(2,0) = fac * ( x(0,1) - x(1,1));
          grad(2,1) = fac * (-x(0,0) + x(1,0));
        }

      public:
        /**
         * \brief Computes the optimal local mesh size
         */
        static void compute_h(ScalarVectorType& h_, const VectorType& coords_, const ScalarVectorType& lambda_,
        const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          DataType_ sum_det = compute_sum_det(coords_, mesh_);
          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);

          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            h_(cell, Math::pow(lambda_(cell)*sum_det,exponent));
          }
        }

        /**
         * \brief Computes the sum of the local transformation's determinants
         **/
        static DataType_ compute_sum_det(const VectorType& coords_, const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          Tiny::Matrix <DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;

          // Index set for local/global numbering
          const auto& idx = mesh_.template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = coords_(idx(cell,Index(j)));

            sum_det += compute_det(x);
          }

#ifdef FEAT_HAVE_MPI
          Util::Comm::allreduce(&sum_det, &sum_det, 1, Util::CommOperationSum());
#endif
          return sum_det;
        }

        /**
         * \brief Computes the gradient of sum of the local transformation's determinants
         */
        static void compute_grad_sum_det(VectorType& grad_, const VectorType& coords_, const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          Tx x(DataType_(0));

          // Index set for local/global numbering
          const auto& idx = mesh_.template get_index_set<ShapeType::dimension,0>();

          Tx local_grad(DataType_(0));

          grad_.format();

          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = coords_(idx(cell,Index(j)));

            compute_grad_det(local_grad, x);

            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              grad_(i, grad_(i) + local_grad[j]);
            }

          }
        } // compute_grad_sum_det

      /**
       * \brief Computes the local matrix describing the material behaviour
       */
        static MatTensorType compute_mat_tensor(const Tx& DOXY(x), const DataType_& h)
        {
          MatTensorType mat_tensor(DataType_(0));
          DataType scal(DataType(1)/h);

          //mat_tensor(0,0) = scal;

          //mat_tensor(1,0) = -scal/Math::sqrt(DataType(3));
          //mat_tensor(1,1) = scal*DataType_(2)/Math::sqrt(DataType(3));

          mat_tensor(0,0) = scal;
          mat_tensor(0,1) = -scal/Math::sqrt(DataType(3));

          mat_tensor(1,1) = scal*DataType_(2)/Math::sqrt(DataType(3));

          return mat_tensor;
        }

    }; // RumpfTrafo<Simplex<2>>

    /**
     * \brief Computes local mesh size h for Rumpf smoothers
     *
     * Specialisation for Simplex<3> meshes
     *
     * \tparam DataType_
     * Our data type
     *
     **/
    template<typename DataType_>
    struct RumpfTrafo<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<3>,3,3,DataType_>>,DataType_>
    {
      public:
        /// Our transformation type
        typedef Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3,DataType_>> TrafoType;
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type, as the Rumpf smoother needs to know what we are working with
        typedef typename TrafoType::ShapeType ShapeType;
        /// Type of the transformation's underlying mesh
        typedef typename TrafoType::MeshType MeshType;
        /// Type for a pack of local vertex coordinates
        typedef Tiny::Matrix <DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> Tx;
        /// Type for the material tensor
        typedef Tiny::Matrix<DataType_, MeshType::world_dim, MeshType::world_dim> MatTensorType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<Mem::Main, DataType, Index> ScalarVectorType;
        /// Vector type for coordinate vectors etc.
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType, Index, MeshType::world_dim> VectorType;

      private:
        /**
         * \brief Compute the transformation's determinant on a cell
         **/
        static DataType compute_det(Tx& x)
        {
          return -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (x(0,0) * x(1,1) * x(2,2) - x(0,0) * x(1,1) * x(3,2) - x(0,0) * x(1,2) * x(2,1) + x(0,0) * x(1,2) * x(3,1) + x(0,0) * x(2,1) * x(3,2) - x(0,0) * x(2,2) * x(3,1) - x(0,1) * x(1,0) * x(2,2) + x(0,1) * x(1,0) * x(3,2) + x(0,1) * x(1,2) * x(2,0) - x(0,1) * x(1,2) * x(3,0) - x(0,1) * x(2,0) * x(3,2) + x(0,1) * x(2,2) * x(3,0) + x(0,2) * x(1,0) * x(2,1) - x(0,2) * x(1,0) * x(3,1) - x(0,2) * x(1,1) * x(2,0) + x(0,2) * x(1,1) * x(3,0) + x(0,2) * x(2,0) * x(3,1) - x(0,2) * x(2,1) * x(3,0) - x(1,0) * x(2,1) * x(3,2) + x(1,0) * x(2,2) * x(3,1) + x(1,1) * x(2,0) * x(3,2) - x(1,1) * x(2,2) * x(3,0) - x(1,2) * x(2,0) * x(3,1) + x(1,2) * x(2,1) * x(3,0)) / DataType(3);
        }

        /**
         * \brief Computes the transformation's determinant's gradient on a cell
         */
        static inline void compute_grad_det(Tx& grad, const Tx& x)
        {
          grad(0,0) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (x(1,1) * x(2,2) - x(1,1) * x(3,2) - x(2,1) * x(1,2) + x(3,1) * x(1,2) + x(2,1) * x(3,2) - x(3,1) * x(2,2)) / DataType(3);
          grad(0,1) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (-x(1,0) * x(2,2) + x(1,0) * x(3,2) + x(2,0) * x(1,2) - x(3,0) * x(1,2) - x(2,0) * x(3,2) + x(3,0) * x(2,2)) / DataType(3);
          grad(0,2) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (x(1,0) * x(2,1) - x(1,0) * x(3,1) - x(2,0) * x(1,1) + x(3,0) * x(1,1) + x(2,0) * x(3,1) - x(3,0) * x(2,1)) / DataType(3);

          grad(1,0) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (-x(0,1) * x(2,2) + x(0,1) * x(3,2) + x(2,1) * x(0,2) - x(3,1) * x(0,2) - x(2,1) * x(3,2) + x(3,1) * x(2,2)) / DataType(3);
          grad(1,1) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (x(0,0) * x(2,2) - x(0,0) * x(3,2) - x(2,0) * x(0,2) + x(3,0) * x(0,2) + x(2,0) * x(3,2) - x(3,0) * x(2,2)) / DataType(3);
          grad(1,2) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (-x(0,0) * x(2,1) + x(0,0) * x(3,1) + x(2,0) * x(0,1) - x(3,0) * x(0,1) - x(2,0) * x(3,1) + x(3,0) * x(2,1)) / DataType(3);

          grad(2,0) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (x(0,1) * x(1,2) - x(0,1) * x(3,2) - x(1,1) * x(0,2) + x(3,1) * x(0,2) + x(1,1) * x(3,2) - x(3,1) * x(1,2)) / DataType(3);
          grad(2,1) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (-x(0,0) * x(1,2) + x(0,0) * x(3,2) + x(1,0) * x(0,2) - x(3,0) * x(0,2) - x(1,0) * x(3,2) + x(3,0) * x(1,2)) / DataType(3);
          grad(2,2) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (x(0,0) * x(1,1) - x(0,0) * x(3,1) - x(1,0) * x(0,1) + x(3,0) * x(0,1) + x(1,0) * x(3,1) - x(3,0) * x(1,1)) / DataType(3);

          grad(3,0) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (-x(0,1) * x(1,2) + x(0,1) * x(2,2) + x(1,1) * x(0,2) - x(2,1) * x(0,2) - x(1,1) * x(2,2) + x(2,1) * x(1,2)) / DataType(3);
          grad(3,1) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (x(0,0) * x(1,2) - x(0,0) * x(2,2) - x(1,0) * x(0,2) + x(2,0) * x(0,2) + x(1,0) * x(2,2) - x(2,0) * x(1,2)) / DataType(3);
          grad(3,2) = -Math::sqrt(DataType(3)) * Math::sqrt(DataType(6)) * (-x(0,0) * x(1,1) + x(0,0) * x(2,1) + x(1,0) * x(0,1) - x(2,0) * x(0,1) - x(1,0) * x(2,1) + x(2,0) * x(1,1)) / DataType(3);


        }

      public:
        /**
         * \brief Computes the optimal local mesh size
         */
        static void compute_h(ScalarVectorType& h_, const VectorType& coords_, const ScalarVectorType& lambda_,
        const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          DataType_ sum_det = compute_sum_det(coords_, mesh_);
          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);

          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            h_(cell, Math::pow(lambda_(cell)*sum_det,exponent));
          }
        }

        /**
         * \brief Computes the sum of the local transformation's determinants
         **/
        static DataType_ compute_sum_det(const VectorType& coords_, const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          Tiny::Matrix <DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;

          // Index set for local/global numbering
          const auto& idx = mesh_.template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              x[j] = coords_(idx(cell,Index(j)));
            }

            sum_det += compute_det(x);
          }

#ifdef FEAT_HAVE_MPI
          Util::Comm::allreduce(&sum_det, &sum_det, 1, Util::CommOperationSum());
#endif
          return sum_det;
        }

        /**
         * \brief Computes the gradient of sum of the local transformation's determinants
         */
        static void compute_grad_sum_det(VectorType& grad_, const VectorType& coords_, const MeshType& mesh_)
        {
          // This will hold the coordinates for one element for passing to other routines
          Tx x(DataType_(0));

          // Index set for local/global numbering
          const auto& idx = mesh_.template get_index_set<ShapeType::dimension,0>();

          Tx local_grad(DataType_(0));

          grad_.format();

          for(Index cell(0); cell < mesh_.get_num_entities(ShapeType::dimension); ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              x[j] = coords_(idx(cell,Index(j)));
            }

            compute_grad_det(local_grad, x);

            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              grad_(i, grad_(i) + local_grad[j]);
            }

          }
        } // compute_grad_sum_det

      /**
       * \brief Computes the local matrix describing the material behaviour
       */
        static MatTensorType compute_mat_tensor(const Tx& DOXY(x), const DataType_& h)
        {
          MatTensorType mat_tensor(DataType_(0));
          DataType scal(DataType(1)/h);

          mat_tensor(0,0) = scal;
          mat_tensor(0,1) = -scal/Math::sqrt(DataType(3));
          mat_tensor(0,2) = -scal/Math::sqrt(DataType(6));

          mat_tensor(1,1) = scal*DataType_(2)/Math::sqrt(DataType(3));
          mat_tensor(1,2) = -scal/Math::sqrt(DataType(6));

          mat_tensor(2,2) = scal*DataType(0.5)*Math::sqrt(DataType(6));

          return mat_tensor;
        }


    }; // RumpfTrafo<Simplex<3>>
#ifdef FEAT_EICKT
    extern template struct RumpfTrafo<
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>, double >;
    extern template struct RumpfTrafo<
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, double>>, double >;
    extern template struct RumpfTrafo<
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, double>>, double >;
    extern template struct RumpfTrafo<
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, double>>, double >;
#endif // FEAT_EICKT

    /// \endcond
  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_RUMPF_TRAFO_HPP
