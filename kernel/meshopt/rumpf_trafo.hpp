#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_TRAFO_HPP
#define KERNEL_MESHOPT_RUMPF_TRAFO_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/comm_base.hpp>

namespace FEAST
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
      /**
       * \brief Compute the transformation's determinant on a cell
       *
       * Meaning the determinant of the transformation from the Rumpf reference element onto the current element.
       *
       * \tparam Tx_
       * Type for the cell coordinates x
       *
       * \param[in] x
       * Cell coordinates
       *
       * \returns
       * The transformation's determinant
       *
       **/
      template<typename Tx_>
      static DataType compute_det(Tx_& x);
      /**
       * \brief Computes the sum of the local transformation's determinants
       *
       * This is needed to correctly calculate the optimal local scales h.
       *
       * \tparam Tcoords_
       * Type of the vector of coordinates of the mesh vertices
       *
       * \param[in] coords_
       * Vector of vertex coordinates
       *
       * \param[in] trafo_
       * The underlying transformation for accessing mesh information
       *
       * \returns The sum of all local transformation's determinants
       *
       **/
      template<typename Tcoords_>
      static DataType_ compute_sum_det(const Tcoords_& coords_, const TrafoType& trafo_);

      /**
       * \brief Computes the optimal local mesh size
       *
       * The optimal local mesh size is calculated from the weights lambda_ such that the volumes of all cells
       * add up to the domain's volume, if each local cell was transformed to a Rumpf reference element with the
       * appropriate size.
       *
       * \tparam Th_
       * Type of the vector of local mesh sizes
       *
       * \tparam Tcoords_
       * Type of the vector of coordinates of the mesh vertices
       *
       * \tparam Tlambda_
       * Type of the vector of (element wise) weights
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
       * \param[in] trafo_
       * The underlying transformation for accessing mesh information
       *
       **/
      template<typename Th_, typename Tcoords_, typename Tlambda_>
      static void compute_h(Th_& h_, const Tcoords_& coords_, const Tlambda_& lambda_, const TrafoType& trafo_);

      /**
       * \brief Computes the transformation's determinant's gradient on a cell
       *
       * Meaning the determinant of the transformation from the Rumpf reference element onto the current element.
       *
       * \tparam Tgrad_
       * Type for the gradient grad
       *
       * \tparam Tx_
       * Type for the cell coordinates x
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
      template<typename Tgrad_, typename Tx_>
      static void compute_grad_det(Tgrad_& grad, const Tx_& x);

      /**
       * \brief Computes the gradient of sum of all cell's local transformation determinants
       *
       * This is the sum of all local gradients and needed to correctly calculate the gradient of the optimal local
       * scales h.
       *
       * \tparam Tcoords_
       * Type of the global vector of coordinates of the mesh vertices
       *
       * \param[out] grad_
       * The global gradient wrt. the local vertices
       *
       * \param[in] coords_
       * Global vector of vertex coordinates
       *
       * \param[in] trafo_
       * The underlying transformation for accessing mesh information
       *
       **/
      template<typename Tcoords_>
      static void compute_grad_sum_det(Tcoords_& grad_, const Tcoords_& coords_, const TrafoType& trafo_);
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
    struct RumpfTrafo
    <
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
      DataType_
    >
    {
      public:
        /// Our transformation type
        typedef typename Trafo::Standard::Mapping
        <
          Geometry::ConformalMesh
          <
            Shape::Hypercube<2>,
            Shape::Hypercube<2>::dimension,
            Shape::Hypercube<2>::dimension,
            DataType_
          >
        > TrafoType;
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type, as the Rumpf smoother needs to know what we are working with
        typedef typename TrafoType::ShapeType ShapeType;
        /// Type of the transformation's underlying mesh
        typedef typename TrafoType::MeshType MeshType;

      private:
        /**
         * \brief Compute the transformation's determinant on a cell
         */
        template<typename Tx_>
        static DataType compute_det(Tx_& x)
        {
          return (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(0,1) * x(1,0) + x(0,1) * x(2,0)
              + x(1,0) * x(3,1) - x(1,1) * x(3,0) - x(2,0) * x(3,1) + x(2,1) * x(3,0)) / DataType(2);
        }

        /**
         * \brief Computes the transformation's determinant's gradient on a cell
         */
        template<typename Tgrad_, typename Tx_>
        static void compute_grad_det(Tgrad_& grad, const Tx_& x)
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
        template<typename Th_, typename Tcoords_, typename Tlambda_>
        static void compute_h(Th_& h_, const Tcoords_& coords_, const Tlambda_& lambda_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          DataType_ sum_det = compute_sum_det(coords_, trafo_);
          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);

          for(Index cell(0); cell < ncells; ++cell)
          {
            Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp;
            for(int d(0); d < MeshType::world_dim; ++d)
              // For hypercubes, h is half the refence cell's edge length
              tmp(d) = DataType(0.5)*Math::pow(lambda_(cell)*sum_det,exponent);

            h_(cell,tmp);
          }
        }

        /**
         * \brief Computes the sum of the local transformation's determinants
         **/
        template<typename Tcoords_>
        static DataType_ compute_sum_det(const Tcoords_& coords_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = coords_(idx(cell,Index(j)));

            sum_det += compute_det(x);
          }
          Comm::allreduce(&sum_det, 1, &sum_det);
          return sum_det;
        }

        /**
         * \brief Computes the gradient of sum of the local transformation's determinants
         **/
        template<typename Tcoords_>
        static void compute_grad_sum_det(Tcoords_& grad_, const Tcoords_& coords_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix<DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          FEAST::Tiny::Matrix<DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim>
            local_grad(DataType_(0));

          grad_.format();

          for(Index cell(0); cell < ncells; ++cell)
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

    }; // RumpfTrafo<Hypercube<2>>

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
    struct RumpfTrafo
    <
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
      DataType_
    >
    {
      public:
        /// Our transformation type
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
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type, as the Rumpf smoother needs to know what we are working with
        typedef typename TrafoType::ShapeType ShapeType;
        /// Type of the transformation's underlying mesh
        typedef typename TrafoType::MeshType MeshType;

      private:
        /**
         * \brief Compute the transformation's determinant on a cell
         **/
        template<typename Tx_>
        static DataType compute_det(Tx_& x)
        {
          return DataType(2) / Math::sqrt(DataType(3)) *
            (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(2,0) * x(0,1) + x(1,0) * x(2,1) - x(2,0) * x(1,1));
        }

        /**
         * \brief Computes the transformation's determinant's gradient on a cell
         **/
        template<typename Tgrad_, typename Tx_>
        static void compute_grad_det(Tgrad_& grad, const Tx_& x)
        {
          DataType fac(DataType(2) / Math::sqrt(DataType(3)));

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
         **/
        template<typename Th_, typename Tcoords_, typename Tlambda_>
        static void compute_h(Th_& h_, const Tcoords_& coords_, const Tlambda_& lambda_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          DataType_ sum_det = compute_sum_det(coords_, trafo_);
          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);

          for(Index cell(0); cell < ncells; ++cell)
          {
            Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim>
              tmp(Math::pow(lambda_(cell)*sum_det,exponent));
            h_(cell, tmp);
          }
        }

        /**
         * \brief Computes the sum of the local transformation's determinants
         **/
        template<typename Tcoords_>
        static DataType_ compute_sum_det(const Tcoords_& coords_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = coords_(idx(cell,Index(j)));

            sum_det += compute_det(x);
          }
          Comm::allreduce(&sum_det, 1, &sum_det);
          return sum_det;
        }

        /**
         * \brief Computes the gradient of sum of the local transformation's determinants
         **/
        template<typename Tcoords_>
        static void compute_grad_sum_det(Tcoords_& grad_, const Tcoords_& coords_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix<DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          FEAST::Tiny::Matrix<DataType_, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim>
            local_grad(DataType_(0));

          grad_.format();

          for(Index cell(0); cell < ncells; ++cell)
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

    }; // RumpfTrafo<Simplex<2>>
    /// \endcond
  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_RUMPF_TRAFO_HPP
