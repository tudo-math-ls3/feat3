#pragma once
#ifndef KERNEL_TRAFO_STANDARD_VOLUME_HPP
#define KERNEL_TRAFO_STANDARD_VOLUME_HPP 1

// includes, FEAST
#include <kernel/util/math.hpp>
#include <kernel/util/meta_math.hpp>

namespace FEAST
{
  namespace Trafo
  {
    namespace Standard
    {

      /**
       * \brief Cell volume computation
       *
       * This class that wraps routines for computing cell volumes for hypercubes and simplices based on the
       * standard transformation in 1D/2D/3D. Computation of lower dimensional volumes of objects immersed in higher
       * dimensional spaces (i.e. triangles immersed in 3D) is supported.
       *
       * \tparam ShapeType_
       * The underlying mesh's shape type
       *
       * \author Jordi Paul
       */
      template<typename ShapeType_>
      struct CellVolumeEvaluator
      {
        template
        <
          /// Data type
          typename DataType_,
          /// Transformation type
          typename TrafoType_,
          /// Type for index sets
          typename IndexSetType_,
          /// Type for vertex sets
          typename VertexSetType_
        >
        /**
         * \brief The routine that does the actual work
         *
         * \param[in] trafo
         * The underlying transformation, needed if the cell volume can only be computed by integration.
         *
         * \param[in] idx
         * The index set of the underlying mesh of trafo
         *
         * \param[in] vtx
         * The vertex set of the underlying mesh of trafo
         *
         * \param[in] cell
         * Cell index
         *
         * \returns The volume of cell.
         **/
        static DataType_ compute_vol(
          const TrafoType_ trafo,
          const IndexSetType_ idx,
          const VertexSetType_ vtx,
          const Index cell);
      };

      /// \copydoc CellVolumeEvaluator
      /// Specialisation for Hypercube<1> shape
      template<>
        struct CellVolumeEvaluator<Shape::Hypercube<1>>
        {
          template
          <
            typename DataType_,
            typename TrafoType_,
            typename IndexSetType_,
            typename VertexSetType_
          >
          /// \copydoc CellVolumeEvaluator::compute_vol()
          static DataType_ compute_vol(
            const TrafoType_& DOXY(trafo),
            const IndexSetType_& idx,
            const VertexSetType_& vtx,
            const Index cell)
            {
              DataType_ v = DataType_(0.);
              for(Index d = 0; d < Index(TrafoType_::MeshType::world_dim); d++)
              {
                v += Math::sqr(DataType_(vtx[idx(cell,1)][d]) - DataType_(vtx[idx(cell,0)][d]));
              }
              return Math::sqrt(v);
            }
        };

      /// Specialisation for Hypercubes of dimension 2 and 3
      ///
      /// As the Hypercube<2> could be immersed in 3D, the Gaussian formula for general polyhedra in 2D is not
      /// applicable, so integration is used.
      ///
      /// For Hypercube<3> shape based on the standard Q1 transformation, faces no longer need to be planar, so
      /// integration is necessary anyway.
      template<int shape_dim_>
        struct CellVolumeEvaluator<Shape::Hypercube<shape_dim_>>
        {
          template
          <
            typename DataType_,
            typename TrafoType_,
            typename IndexSetType_,
            typename VertexSetType_
          >
          /// \copydoc CellVolumeEvaluator::compute_vol()
          static DataType_ compute_vol(
            const TrafoType_& trafo,
            const IndexSetType_& DOXY(idx),
            const VertexSetType_& DOXY(vtx),
            const Index cell)
            {
              typedef typename TrafoType_::template Evaluator<Shape::Hypercube<shape_dim_>, DataType_>::Type Evaluator;
              typedef typename Evaluator::DomainPointType DomainPointType;
              typedef typename Evaluator::JacobianMatrixType JacobianMatrixType;

              /* The integration is done with a hard coded tensor product of the two point 1d Gauss-Legendre
               * quadrature, which acts on [-1, 1] with quadrature points xq \in {-sqrt(1/3), sqrt(1/3)} and weights
               * wq = 1
               */
              const DataType_ x = DataType_(1.)/Math::sqrt(DataType_(3.));

              // Evaluator for the trafo
              Evaluator trafo_eval(trafo);
              trafo_eval.prepare(cell);
              // Jacobian of the trafo
              JacobianMatrixType jac_mat;
              // This is where the integration points go
              DomainPointType xq;

              DataType_ vol = DataType_(0.);
              // Integrate
              for(Index i(0); i < Index(1 << shape_dim_); ++i)
              {
                // Set point coords using bloody bitshifts
                for(Index j(0); j < Index(shape_dim_); ++j)
                {
                  xq(j) = (DataType_(((i >> j) & 1) << 1) - DataType_(1)) * x;
                }

                // Get Jacobian
                trafo_eval.calc_jac_mat(jac_mat, xq);
                // vol += sqrt( det( J^T J ) )
                vol += jac_mat.vol();
              }

              trafo_eval.finish();

              return vol;
            }
        };

      /// Specialisation for simplices, dimension independent.
      template<int shape_dim_>
        struct CellVolumeEvaluator<Shape::Simplex<shape_dim_>>
        {
          /// \copydoc CellVolumeEvaluator::compute_vol()
          template<
            typename DataType_,
            typename TrafoType_,
            typename IndexSetType_,
            typename VertexSetType_>
            static DataType_ compute_vol(
                const TrafoType_& DOXY(trafo),
                const IndexSetType_& idx,
                const VertexSetType_& vtx,
                const Index cell)
            {
              // A consist of all vectors of type [v_(d+1) - v_0, ..., v_1 - v_0], where v_0, ... , v_(d+1) are
              // the world_dim_ vectors of the vertices comprising the shape_dim_ simplex.
              FEAST::Tiny::Matrix< DataType_,
              TrafoType_::MeshType::world_dim, shape_dim_,
              TrafoType_::MeshType::world_dim, shape_dim_ > A(0.);

              // Fill the matrix
              for(Index i = 0; i < Index(TrafoType_::MeshType::world_dim); i++ )
              {
                for(Index j = 0; j < shape_dim_; j++)
                {
                  A(i,j) = DataType_(vtx[idx(cell,j+1)][i]) - DataType_(vtx[idx(cell,0)][i]);
                }
              }

              // Never forget to scale with the volume of the standard simplex.
              return A.vol()/DataType_(MetaMath::Factorial<shape_dim_>::value);
            }
        };

    } // namespace Standard
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_STANDARD_VOLUME_HPP
