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
     * \author Jordi Paul
     *
     **/
#ifndef DOXYGEN
    template<typename TrafoType_, typename DataType_>
    struct H_Evaluator;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation
    // of specialisations in TrafoType_.
#else
    template<typename TrafoType_, typename DataType_>
    struct H_Evaluator
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
       * \brief Computes the gradient of sum of the local transformation's determinants
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
    }; // struct H_Evaluator;
#endif

    /**
     * \brief Computes local mesh size h for Q1 hack funtionals
     *
     * As the functionals use hypercubes split into simplices, the optimal scales have to be computed differently.
     *
     * \tparam TrafoType_
     * Transformation type, the ShapeType is essential
     *
     * \tparam DataType_
     * Our data type
     *
     *
     * \author Jordi Paul
     *
     **/
#ifndef DOXYGEN
    template<typename TrafoType_, typename DataType_>
    struct H_EvaluatorQ1Hack;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation
    // of specialisations in TrafoType_.
#else
    template<typename TrafoType_, typename DataType_>
    struct H_EvaluatorQ1Hack
    {
      /**
       * \copydoc H_Evaluator::compute_det()
       **/
      template<typename Tx_>
      static DataType compute_det(Tx_& x);
      /**
       * \copydoc H_Evaluator::compute_sum_det()
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
      static void compute_h(Th_& h_, const Tcoords_& coords_, const Tlambda_& lambda_, const TrafoType& trafo_);

      /**
       * \copydoc H_Evaluator::compute_grad_det()
       **/
      template<typename Tgrad_, typename Tx_>
      static void compute_grad_det(Tgrad_& grad, const Tx_& x);

      /**
       * \copydoc H_Evaluator::compute_grad_sum_det()
       **/
      template<typename Tcoords_>
      static void compute_grad_sum_det(Tcoords_& grad_, const Tcoords_& coords_, const TrafoType& trafo_);
    }; // struct H_Evaluator;
#endif

    /// \cond internal
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
         **/
        template<typename Tx_>
        static DataType compute_det(Tx_& x)
        {
          return DataType(0.5000000000e0 * x(0,3) * x(1,2) + 0.5000000000e0 * x(0,1) * x(1,3) - 0.5000000000e0 * x(0,0) * x(1,2) - 0.5000000000e0 * x(0,1) * x(1,0) - 0.5000000000e0 * x(0,2) * x(1,3) + 0.5000000000e0 * x(0,2) * x(1,0) - 0.5000000000e0 * x(0,3) * x(1,1) + 0.5000000000e0 * x(0,0) * x(1,1));

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

          DataType_ exponent = DataType_(1)/DataType_(MeshType::world_dim);
          DataType_ sum_det = compute_sum_det(coords_, trafo_);

          // For given lambda = volume(cell), how do we need to chose h such that two Rumpf reference simplices of
          // scale h have the same volume?
          // They have volume = 2 * sqrt(3)/4 * h, so h = 2/sqrt(3) * volume =: const1 * volume
          DataType_ const1 = DataType(2)/(Math::sqrt(DataType(3)));
          // For given lambda = volume(cell), how do we need to chose h such that the det + 1/det part of the
          // functional has a local minimum for two Rumpf reference simplices of scale h?
          // Lengthy computation with maple says that this is the right constant.
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

        /**
         * \brief Computes the sum of the local transformation's determinants
         **/
        template<typename Tcoords_>
        static DataType_ compute_sum_det(const Tcoords_& coords_, const TrafoType& trafo_)
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
          return sum_det;
        }

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
         **/
        template<typename Tx_>
        static DataType compute_det(Tx_& x)
        {
          return DataType(0.5000000000e0 * x(0,3) * x(1,2) + 0.5000000000e0 * x(0,1) * x(1,3) - 0.5000000000e0 * x(0,0) * x(1,2) - 0.5000000000e0 * x(0,1) * x(1,0) - 0.5000000000e0 * x(0,2) * x(1,3) + 0.5000000000e0 * x(0,2) * x(1,0) - 0.5000000000e0 * x(0,3) * x(1,1) + 0.5000000000e0 * x(0,0) * x(1,1));
        }

        /**
         * \brief Computes the transformation's determinant's gradient on a cell
         **/
        template<typename Tgrad_, typename Tx_>
        static void compute_grad_det(Tgrad_& grad, const Tx_& x)
        {
          grad(0,0) = x(1,1) / DataType(0.2e1) - x(1,2) / DataType(0.2e1);
          grad(0,1) = x(1,3) / DataType(0.2e1) - x(1,0) / DataType(0.2e1);
          grad(0,2) = -x(1,3) / DataType(0.2e1) + x(1,0) / DataType(0.2e1);
          grad(0,3) = -x(1,1) / DataType(0.2e1) + x(1,2) / DataType(0.2e1);
          grad(1,0) = -x(0,1) / DataType(0.2e1) + x(0,2) / DataType(0.2e1);
          grad(1,1) = x(0,0) / DataType(0.2e1) - x(0,3) / DataType(0.2e1);
          grad(1,2) = -x(0,0) / DataType(0.2e1) + x(0,3) / DataType(0.2e1);
          grad(1,3) = -x(0,2) / DataType(0.2e1) + x(0,1) / DataType(0.2e1);
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
            for(Index d(0); d < MeshType::world_dim; ++d)
              // For hypercubes, h is half the refence cell's edge length
              h_[d](cell,DataType(0.5)*Math::pow(lambda_(cell)*sum_det,exponent));
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
          FEAST::Tiny::Matrix<DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          FEAST::Tiny::Matrix<DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count>
            local_grad(DataType_(0));

          for(Index d(0); d < MeshType::world_dim; ++d)
            grad_[d].format();

          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              // Get local coordinates
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                x(d,j) = coords_[d](idx(cell,j));
            }

            compute_grad_det(local_grad, x);

            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              {
                DataType tmp = grad_[d](idx(cell,j)) + local_grad(d,j);
                grad_[d](idx(cell,j), tmp);
              }
            }

          }
        } // compute_grad_sum_det

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
         **/
        template<typename Tx_>
        static DataType compute_det(Tx_& x)
        {
          return DataType(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) );
        }

        /**
         * \brief Computes the transformation's determinant's gradient on a cell
         **/
        template<typename Tgrad_, typename Tx_>
        static void compute_grad_det(Tgrad_& grad, const Tx_& x)
        {
          grad(0,0) = DataType_(0.2e1) / DataType_(0.3e1) * Math::sqrt(DataType_(0.3e1)) * (x(1,1) - x(1,2));
          grad(0,1) = DataType_(0.2e1) / DataType_(0.3e1) * (-x(1,0) + x(1,2)) * Math::sqrt(DataType_(0.3e1));
          grad(1,0) = DataType_(0.2e1) / DataType_(0.3e1) * Math::sqrt(DataType_(0.3e1)) * (-x(0,1) + x(0,2));
          grad(1,1) = DataType_(0.2e1) / DataType_(0.3e1) * Math::sqrt(DataType_(0.3e1)) * (x(0,0) - x(0,2));
          grad(0,2) = DataType_(0.2e1) / DataType_(0.3e1) * Math::sqrt(DataType_(0.3e1)) * (x(1,0) - x(1,1));
          grad(1,2) = DataType_(0.2e1) / DataType_(0.3e1) * (-x(0,0) + x(0,1)) * Math::sqrt(DataType_(0.3e1));
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
            for(Index d(0); d < MeshType::world_dim; ++d)
              h_[d](cell,Math::pow(lambda_(cell)*sum_det,exponent));
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
          FEAST::Tiny::Matrix<DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;

          // Index set for local/global numbering
          auto& idx = trafo_.get_mesh().template get_index_set<ShapeType::dimension,0>();

          FEAST::Tiny::Matrix<DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count>
            local_grad(DataType_(0));

          for(Index d(0); d < MeshType::world_dim; ++d)
            grad_[d].format();

          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              // Get local coordinates
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                x(d,j) = coords_[d](idx(cell,j));
            }

            compute_grad_det(local_grad, x);

            for(Index d(0); d < MeshType::world_dim; ++d)
            {
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              {
                DataType tmp = grad_[d](idx(cell,j)) + local_grad(d,j);
                grad_[d](idx(cell,j), tmp);
              }
            }

          }
        } // compute_grad_sum_det

    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_H_EVALUATOR_HPP
