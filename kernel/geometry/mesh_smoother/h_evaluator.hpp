#pragma once
#ifndef KERNEL_GEOMETRY_H_EVALUATOR_HPP
#define KERNEL_GEOMETRY_H_EVALUATOR_HPP 1

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Computes local mesh size h
     *
     * Generic class template for H_Evaluators.
     *
     * \tparam TrafoType_
     * Transformation type, the ShapeType is essential
     *
     * \tparam DataType_
     * Our data type
     *
     *
     * \autho Jordi Paul
     *
     **/
    template<typename TrafoType_, typename DataType_>
    struct H_Evaluator;

    template<typename TrafoType_, typename DataType_>
    struct H_EvaluatorQ1Hack;

    /**
     * \brief Computes local mesh size h for the Q1Hack Rumpf smoothers
     *
     * Specialisation for Hypercube<d> meshes
     *
     * \tparam DataType_
     * Our data type
     *
     **/
    template<typename DataType_>
    struct H_EvaluatorQ1Hack
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
        static DataType compute_det(Tx_& x)
        {
          return DataType(0.5000000000e0 * x(0,3) * x(1,2) + 0.5000000000e0 * x(0,1) * x(1,3) - 0.5000000000e0 * x(0,0) * x(1,2) - 0.5000000000e0 * x(0,1) * x(1,0) - 0.5000000000e0 * x(0,2) * x(1,3) + 0.5000000000e0 * x(0,2) * x(1,0) - 0.5000000000e0 * x(0,3) * x(1,1) + 0.5000000000e0 * x(0,0) * x(1,1));

        }
      public:
        /**
         * \brief Computes the optimal local mesh size
         *
         * The optimal local mesh size is calculated from the weights lambda_ such that the volumes of all cells
         * add up to the domain's volume, if each local cell was transformed to a Rumpf reference element with the
         * appropriate size.
         *
         * The functionals use simplical meshes, but because of the different scaling of the terms norm_A and
         * det_A / 1/det_A, two different optimal mesh sizes have to be computed: h[0] is for the norm_A term
         * and h[1] for everything else.
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
        static void compute_h(Th_& h_, const Tcoords_& coords_, const Tlambda_& lambda_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              // Get local coordinates
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                x(d,j) = coords_[d](idx(cell,j));
            }
            sum_det += compute_det(x);
          }

          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);
          DataType_ const1 = DataType(2)/(Math::sqrt(DataType(3)));
          DataType_ const2 = DataType(2)/(Math::sqrt(DataType(2)*Math::sqrt(DataType(3))));
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              h_[0](cell,const1*Math::pow(lambda_(cell)*sum_det,exponent));
              h_[1](cell,const2*Math::pow(lambda_(cell)*sum_det,exponent));
            }
          }
        } // compute_h

    }; // struct H_EvaluatorQ1Hack


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
    struct H_Evaluator
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
        static DataType compute_det(Tx_& x)
        {
          return DataType(0.5000000000e0 * x(0,3) * x(1,2) + 0.5000000000e0 * x(0,1) * x(1,3) - 0.5000000000e0 * x(0,0) * x(1,2) - 0.5000000000e0 * x(0,1) * x(1,0) - 0.5000000000e0 * x(0,2) * x(1,3) + 0.5000000000e0 * x(0,2) * x(1,0) - 0.5000000000e0 * x(0,3) * x(1,1) + 0.5000000000e0 * x(0,0) * x(1,1));

        }

      public:
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
        static void compute_h(Th_& h_, const Tcoords_& coords_, const Tlambda_& lambda_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          //DataType_ target_vol(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              // Get local coordinates
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                x(d,j) = coords_[d](idx(cell,j));
            }
            sum_det += compute_det(x);
          }

          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
              // For hypercubes, h is half the refence cell's edge length
              h_[d](cell,DataType(0.5)*Math::pow(lambda_(cell)*sum_det,exponent));

          }
        }

    };

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
    struct H_Evaluator
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
        static DataType compute_det(Tx_& x)
        {
          return DataType(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) );
        }

      public:
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
        static void compute_h(Th_& h_, const Tcoords_& coords_, const Tlambda_& lambda_, const TrafoType& trafo_)
        {
          Index ncells( trafo_.get_mesh().get_num_entities(ShapeType::dimension) );

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          DataType sum_det(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              // Get local coordinates
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                x(d,j) = coords_[d](idx(cell,j));
            }
            sum_det += compute_det(x);
          }

          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
              h_[d](cell,Math::pow(lambda_(cell)*sum_det,exponent));

          }
        }

    };
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_H_EVALUATOR_HPP
