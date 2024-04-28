#pragma once
#ifndef KERNEL_TRAFO_STANDARD_DETAILS_HPP
#define KERNEL_TRAFO_STANDARD_DETAILS_HPP 1

// This file defines detail implementations as templated functions.
// This is mainly done due to avoid code duplication for the voxel assembly.
// \author Maximilian Esser, Peter Zajac

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>
#ifndef __CUDA_ARCH__
#include <kernel/util/math.hpp>
#else
#include <kernel/util/cuda_math.cuh>
#endif

namespace FEAT
{
  namespace Trafo
  {
    namespace Standard
    {
      /**
       * \brief Evalautator helper class for standard trafo
       *
       * This wrapper around a number of static functions provides easy access to the implementation details of the standard trafo evaluator
       * for different shape types. All provided static function are defined as device and host functions in the cuda context.
       *
       * \tparam DataType_ The datatype of the trafo coefficients.
       * \tparam DomPointType_ The datatype, i.e. a coordinate Vector, of the domain space.
       * \tparam ImgPointType_ The datatype of the image space.
       * \tparam JacMatType_ The datatype of the jacobian of the trafo.
       * \tparam JacMatInvType_ The type of the inverse jacobian.
       * \tparam JacDetType_ The datatype of the determinant of the jacobian.
       * \tparam HessType_ The type of the hessian tensor3.
       * \tparam HessInvType_ The type of the inverse of the hessian tensor.
       * \tparam Shape_ The shapetype of the trafo.
       * \tparam img_dim The dimension of the target finite element space (i.e. the dimension of the coefficients).
       *
       * \author Maximilian Esser
       */
      template<typename DataType_, typename DomPointType_, typename ImgPointType_, typename JacMatType_, typename JacMatInvType_, typename JacDetType_,
               typename HessType_, typename HessInvType_, typename Shape_, int img_dim>
      struct EvalHelper;


      template<typename DataType_, typename DomPointType_, typename ImgPointType_, typename JacMatType_, typename JacMatInvType_, typename JacDetType_,
               typename HessType_, typename HessInvType_, int img_dim>
      struct EvalHelper<DataType_, DomPointType_, ImgPointType_, JacMatType_, JacMatInvType_, JacDetType_, HessType_, HessInvType_, Shape::Simplex<1>, img_dim>
      {
        /// evaluation data type
        typedef DataType_ DataType;
        /// domain point type
        typedef DomPointType_ DomainPointType;
        /// image point type
        typedef ImgPointType_ ImagePointType;
        /// jacobian matrix type
        typedef JacMatType_ JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef JacMatInvType_ JacobianInverseType;
        /// jacobian determinant type
        typedef JacDetType_ JacobianDeterminantType;
        /// hessian tensor type
        typedef HessType_ HessianTensorType;
        /// hessian inverse tensor type
        typedef HessInvType_ HessianInverseType;
        /// the shape type
        typedef Shape::Simplex<1> ShapeType;

        static constexpr int num_verts = Shape::FaceTraits<ShapeType, 0>::count;
        static constexpr int image_dim = img_dim;

        /**
         * \brief Sets the coefficient vector based on the underlying standard transformation.
         *
         * \tparam VertexType_ The type of the vertex array.
         * \tparam IndexTupleType_ A type mapping the local index to the dofs. Has to support operator[](Index i).
         *
         * \param[out] coeffs
         * A two dimensional array holding the local coefficients of the trafo after the call.
         *
         * \param[in] index_tuple
         * A mapping from local indices to global.
         *
         * \param[in] vertex_set
         * A pointer to the beginnnig of the vertex set.
         */
        template<typename VertexType_, typename IndexTupleType_>
        CUDA_HOST_DEVICE static void inline set_coefficients(DataType (&coeffs)[image_dim][num_verts], const IndexTupleType_& index_tuple, const VertexType_* vertex_set)
        {
          typedef VertexType_ VertexType;

          const VertexType& v0 = vertex_set[index_tuple[0]];
          const VertexType& v1 = vertex_set[index_tuple[1]];


          for(int i(0); i < image_dim; ++i)
          {
            coeffs[i][0] = DataType(v0[i]);
            coeffs[i][1] = DataType(v1[i] - v0[i]);
          }

        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void map_point(ImagePointType& img_point, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = coeffs[i][0] + coeffs[i][1] * dom_point[0];
          }
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& DOXY(dom_point), const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < img_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& DOXY(dom_point), const DataType DOXY(coeffs)[image_dim][num_verts])
        {
          hess_ten.format();
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The volume of the current cell.
         */
        CUDA_HOST_DEVICE static inline DataType volume(const DataType (&coeffs)[image_dim][num_verts])
        {
          DataType v = DataType(0);
          #ifndef __CUDA_ARCH__
          for(int i(0); i < image_dim; ++i)
            v += Math::sqr(coeffs[i][1]);
          return Math::sqrt(v);
          #else
          for(int i(0); i < image_dim; ++i)
            v += CudaMath::cuda_sqr(coeffs[i][1]);
          return CudaMath::cuda_sqrt(v);
          #endif
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        CUDA_HOST_DEVICE static inline DataType width_directed(const ImagePointType& DOXY(ray), const DataType (&coeffs)[image_dim][num_verts])
        {
          // in 1D, the width is always equal to the volume
          return volume(coeffs);
        }
      }; // EvalHelper<Simplex<1>,...>

      template<typename DataType_, typename DomPointType_, typename ImgPointType_, typename JacMatType_, typename JacMatInvType_, typename JacDetType_,
               typename HessType_, typename HessInvType_, int img_dim>
      struct EvalHelper<DataType_, DomPointType_, ImgPointType_, JacMatType_, JacMatInvType_, JacDetType_, HessType_, HessInvType_, Shape::Simplex<2>, img_dim>
      {
        /// evaluation data type
        typedef DataType_ DataType;
        /// domain point type
        typedef DomPointType_ DomainPointType;
        /// image point type
        typedef ImgPointType_ ImagePointType;
        /// jacobian matrix type
        typedef JacMatType_ JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef JacMatInvType_ JacobianInverseType;
        /// jacobian determinant type
        typedef JacDetType_ JacobianDeterminantType;
        /// hessian tensor type
        typedef HessType_ HessianTensorType;
        /// hessian inverse tensor type
        typedef HessInvType_ HessianInverseType;
        /// the shape type
        typedef Shape::Simplex<2> ShapeType;

        static constexpr int num_verts = Shape::FaceTraits<ShapeType, 0>::count;
        static constexpr int image_dim = img_dim;

        /**
         * \brief Sets the coefficient vector based on the underlying standard transformation.
         *
         * \tparam VertexType_ The type of the vertex array.
         * \tparam IndexTupleType_ A type mapping the local index to the dofs. Has to support operator[](Index i).
         *
         * \param[out] coeffs
         * A two dimensional array holding the local coefficients of the trafo after the call.
         *
         * \param[in] index_tuple
         * A mapping from local indices to global.
         *
         * \param[in] vertex_set
         * A pointer to the beginnnig of the vertex set.
         *
         */
        template<typename VertexType_, typename IndexTupleType_>
        CUDA_HOST_DEVICE static void inline set_coefficients(DataType (&coeffs)[img_dim][num_verts], const IndexTupleType_& index_tuple, const VertexType_* vertex_set)
        {
          typedef VertexType_ VertexType;

          const VertexType& v0 = vertex_set[index_tuple[0]];
          const VertexType& v1 = vertex_set[index_tuple[1]];
          const VertexType& v2 = vertex_set[index_tuple[2]];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
          {
            coeffs[i][0] = DataType(v0[i]);
            coeffs[i][1] = DataType(v1[i] - v0[i]);
            coeffs[i][2] = DataType(v2[i] - v0[i]);
          }

        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void map_point(ImagePointType& img_point, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = coeffs[i][0] + coeffs[i][1] * dom_point[0] + coeffs[i][2] * dom_point[1];
          }
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& DOXY(dom_point), const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& DOXY(dom_point), const DataType DOXY(coeffs)[image_dim][num_verts])
        {
          hess_ten.format();
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The volume of the current cell.
         */
        CUDA_HOST_DEVICE static inline DataType volume(const DataType (&coeffs)[image_dim][num_verts])
        {
          JacobianMatrixType jac_mat;
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
          }
          return jac_mat.vol() / DataType(2);
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        CUDA_HOST_DEVICE static inline DataType width_directed(const ImagePointType& ray, const DataType (&coeffs)[image_dim][num_verts])
        {
          // This one is a little tricky:
          // We follow a similar approach as on quadrilaterals here, i.e. we first map the ray
          // vector onto the reference cells by multiplying it to the inverse jacobian matrix.
          // However, our reference cell is the "right triangle" spanned by the three points
          // (0,0) (1,0) (0,1) and therefore its edges have different lengths, which would
          // result in a "skew" width. To circumvent this, we will afterwards map the ray
          // vector onto a "equilateral unit triangle" (i.e. the triangle where all edges
          // have length = 1) and we compute the norm of that mapped vector then.
          JacobianMatrixType jac_mat;
          JacobianInverseType jac_inv;
          DomainPointType ref_ray;

          // compute jacobian matrix
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
          }

          // invert jacobian matrix and multiply by ray vector
          jac_inv.set_inverse(jac_mat);
          ref_ray.set_mat_vec_mult(jac_inv, ray);

          // We have mapped the ray onto our reference triangle. Now let's map that one
          // onto a equilateral unit triangle, which can be performed by multiplying
          // the following matrix by the reference ray vector:
          //
          //  eut_ray := [ -1       +1      ] * ref_ray
          //             [ sqrt(3)  sqrt(3) ]
          //
          // Finally, we have to compute the norm of that vector, which is explicitly
          // written down in the denominator of the following expression:
          #ifndef __CUDA_ARCH__
          return DataType(2) / Math::sqrt(Math::sqr(ref_ray[1]-ref_ray[0]) + DataType(3) * Math::sqr(ref_ray[0]+ref_ray[1]));
          #else
          return DataType(2) * CudaMath::cuda_rsqrt(Math::sqr(ref_ray[1]-ref_ray[0]) + DataType(3) * Math::sqr(ref_ray[0]+ref_ray[1]));
          #endif
        }
      }; // EvalHelper<Simplex<2>,...>

      template<typename DataType_, typename DomPointType_, typename ImgPointType_, typename JacMatType_, typename JacMatInvType_, typename JacDetType_,
               typename HessType_, typename HessInvType_, int img_dim>
      struct EvalHelper<DataType_, DomPointType_, ImgPointType_, JacMatType_, JacMatInvType_, JacDetType_, HessType_, HessInvType_, Shape::Simplex<3>, img_dim>
      {
        /// evaluation data type
        typedef DataType_ DataType;
        /// domain point type
        typedef DomPointType_ DomainPointType;
        /// image point type
        typedef ImgPointType_ ImagePointType;
        /// jacobian matrix type
        typedef JacMatType_ JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef JacMatInvType_ JacobianInverseType;
        /// jacobian determinant type
        typedef JacDetType_ JacobianDeterminantType;
        /// hessian tensor type
        typedef HessType_ HessianTensorType;
        /// hessian inverse tensor type
        typedef HessInvType_ HessianInverseType;
        /// the shape type
        typedef Shape::Simplex<3> ShapeType;

        static constexpr int num_verts = Shape::FaceTraits<ShapeType, 0>::count;
        static constexpr int image_dim = img_dim;

        /**
         * \brief Sets the coefficient vector based on the underlying standard transformation.
         *
         * \tparam VertexType_ The type of the vertex array.
         * \tparam IndexTupleType_ A type mapping the local index to the dofs. Has to support operator[](Index i).
         *
         * \param[out] coeffs
         * A two dimensional array holding the local coefficients of the trafo after the call.
         *
         * \param[in] index_tuple
         * A mapping from local indices to global.
         *
         * \param[in] vertex_set
         * A pointer to the beginnnig of the vertex set.
         */
        template<typename VertexType_, typename IndexTupleType_>
        CUDA_HOST_DEVICE static void inline set_coefficients(DataType (&coeffs)[img_dim][num_verts], const IndexTupleType_& index_tuple, const VertexType_* vertex_set)
        {
          // fetch the vertices of the edge
          typedef VertexType_ VertexType;
          const VertexType& v0 = vertex_set[index_tuple[0]];
          const VertexType& v1 = vertex_set[index_tuple[1]];
          const VertexType& v2 = vertex_set[index_tuple[2]];
          const VertexType& v3 = vertex_set[index_tuple[3]];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
          {
            coeffs[i][0] = DataType(v0[i]);
            coeffs[i][1] = DataType(v1[i] - v0[i]);
            coeffs[i][2] = DataType(v2[i] - v0[i]);
            coeffs[i][3] = DataType(v3[i] - v0[i]);
          }

        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void map_point(ImagePointType& img_point, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = coeffs[i][0] + coeffs[i][1] * dom_point[0]
                                        + coeffs[i][2] * dom_point[1]
                                        + coeffs[i][3] * dom_point[2];
          }
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& DOXY(dom_point), const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
            jac_mat(i,2) = coeffs[i][3];
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& DOXY(dom_point), const DataType DOXY(coeffs)[image_dim][num_verts])
        {
          hess_ten.format();
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The volume of the current cell.
         */
        CUDA_HOST_DEVICE static inline DataType volume(const DataType (&coeffs)[image_dim][num_verts])
        {
          JacobianMatrixType jac_mat;
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
            jac_mat(i,2) = coeffs[i][3];
          }
          return jac_mat.vol() / DataType(6);
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        CUDA_HOST_DEVICE static inline DataType width_directed(const ImagePointType& ray, const DataType (&coeffs)[image_dim][num_verts])
        {
          // We follow the same approach as for triangles here:
          // First, map the ray onto the reference element and then
          // map it to a regular tetrahedron to compute its norm.
          JacobianMatrixType jac_mat;
          JacobianInverseType jac_inv;
          DomainPointType ref_ray;

          // compute jacobian matrix
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
            jac_mat(i,2) = coeffs[i][3];
          }

          // invert jacobian matrix and multiply by ray vector
          jac_inv.set_inverse(jac_mat);
          ref_ray.set_mat_vec_mult(jac_inv, ray);

          // The transformation from the reference tetrahedron to the regular unit tetrahedron
          // can be performed by applying the following matrix-vector product:
          //
          //             [ 2  1         1      ]
          //  eut_ray := [ 0  1        -1      ] * ref_ray
          //             [ 0  sqrt(2)  sqrt(2) ]
          //
          // Finally, we have to compute the norm of that vector, which is explicitly
          // written down in the denominator of the following expression:
          #ifndef __CUDA_ARCH__
          return DataType(2) / Math::sqrt(Math::sqr(DataType(2) * ref_ray[0] + ref_ray[1] + ref_ray[2]) +
            Math::sqr(ref_ray[1] - ref_ray[2]) + DataType(2) * Math::sqr(ref_ray[1] + ref_ray[2]));
          #else
          return DataType(2) * CudaMath::cuda_rsqrt(CudaMath::sqr(DataType(2) * ref_ray[0] + ref_ray[1] + ref_ray[2]) +
            CudaMath::cuda_sqr(ref_ray[1] - ref_ray[2]) + DataType(2) * CudaMath::cuda_sqr(ref_ray[1] + ref_ray[2]));
          #endif
        }
      }; // EvalHelper<Simplex<3>,...>

      template<typename DataType_, typename DomPointType_, typename ImgPointType_, typename JacMatType_, typename JacMatInvType_, typename JacDetType_,
               typename HessType_, typename HessInvType_, int img_dim>
      struct EvalHelper<DataType_, DomPointType_, ImgPointType_, JacMatType_, JacMatInvType_, JacDetType_, HessType_, HessInvType_, Shape::Hypercube<1>, img_dim>
      {
        /// evaluation data type
        typedef DataType_ DataType;
        /// domain point type
        typedef DomPointType_ DomainPointType;
        /// image point type
        typedef ImgPointType_ ImagePointType;
        /// jacobian matrix type
        typedef JacMatType_ JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef JacMatInvType_ JacobianInverseType;
        /// jacobian determinant type
        typedef JacDetType_ JacobianDeterminantType;
        /// hessian tensor type
        typedef HessType_ HessianTensorType;
        /// hessian inverse tensor type
        typedef HessInvType_ HessianInverseType;
        /// the shape type
        typedef Shape::Hypercube<1> ShapeType;

        static constexpr int num_verts = Shape::FaceTraits<ShapeType, 0>::count;
        static constexpr int image_dim = img_dim;

        /**
         * \brief Sets the coefficient vector based on the underlying standard transformation.
         *
         * \tparam VertexType_ The type of the vertex array.
         * \tparam IndexTupleType_ A type mapping the local index to the dofs. Has to support operator[](Index i).
         *
         * \param[out] coeffs
         * A two dimensional array holding the local coefficients of the trafo after the call.
         *
         * \param[in] index_tuple
         * A mapping from local indices to global.
         *
         * \param[in] vertex_set
         * A pointer to the beginnnig of the vertex set.
         */
        template<typename VertexType_, typename IndexTupleType_>
        CUDA_HOST_DEVICE static void inline set_coefficients(DataType (&coeffs)[img_dim][num_verts], const IndexTupleType_& index_tuple, const VertexType_* vertex_set)
        {
          typedef VertexType_ VertexType;

          const VertexType& v0 = vertex_set[index_tuple[0]];
          const VertexType& v1 = vertex_set[index_tuple[1]];


          for(int i(0); i < img_dim; ++i)
          {
            coeffs[i][0] = DataType(0.5) * DataType( v0[i] + v1[i]);
            coeffs[i][1] = DataType(0.5) * DataType(-v0[i] + v1[i]);
          }

        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void map_point(ImagePointType& img_point, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = coeffs[i][0] + coeffs[i][1] * dom_point[0];
          }
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& DOXY(dom_point), const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < img_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& DOXY(dom_point), const DataType DOXY(coeffs)[image_dim][num_verts])
        {
          hess_ten.format();
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The volume of the current cell.
         */
        CUDA_HOST_DEVICE static inline DataType volume(const DataType (&coeffs)[image_dim][num_verts])
        {
          DataType v = DataType(0);
          #ifndef __CUDA_ARCH__
          for(int i(0); i < image_dim; ++i)
            v += Math::sqr(coeffs[i][1]);
          return DataType(2) * Math::sqrt(v);
          #else
          for(int i(0); i < image_dim; ++i)
            v += CudaMath::cuda_sqr(coeffs[i][1]);
          return DataType(2) * CudaMath::cuda_sqrt(v);
          #endif
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        CUDA_HOST_DEVICE static inline DataType width_directed(const ImagePointType& DOXY(ray), const DataType (&coeffs)[image_dim][num_verts])
        {
          // in 1D, the width is always equal to the volume
          return volume(coeffs);
        }
      }; // EvalHelper<Hypercube<1>,...>

      template<typename DataType_, typename DomPointType_, typename ImgPointType_, typename JacMatType_, typename JacMatInvType_, typename JacDetType_,
               typename HessType_, typename HessInvType_, int img_dim>
      struct EvalHelper<DataType_, DomPointType_, ImgPointType_, JacMatType_, JacMatInvType_, JacDetType_, HessType_, HessInvType_, Shape::Hypercube<2>, img_dim>
      {
        /// evaluation data type
        typedef DataType_ DataType;
        /// domain point type
        typedef DomPointType_ DomainPointType;
        /// image point type
        typedef ImgPointType_ ImagePointType;
        /// jacobian matrix type
        typedef JacMatType_ JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef JacMatInvType_ JacobianInverseType;
        /// jacobian determinant type
        typedef JacDetType_ JacobianDeterminantType;
        /// hessian tensor type
        typedef HessType_ HessianTensorType;
        /// hessian inverse tensor type
        typedef HessInvType_ HessianInverseType;
        /// the shape type
        typedef Shape::Hypercube<2> ShapeType;

        static constexpr int num_verts = Shape::FaceTraits<ShapeType, 0>::count;
        static constexpr int image_dim = img_dim;

        /**
         * \brief Sets the coefficient vector based on the underlying standard transformation.
         *
         * \tparam VertexType_ The type of the vertex array.
         * \tparam IndexTupleType_ A type mapping the local index to the dofs. Has to support operator[](Index i).
         *
         * \param[out] coeffs
         * A two dimensional array holding the local coefficients of the trafo after the call.
         *
         * \param[in] index_tuple
         * A mapping from local indices to global.
         *
         * \param[in] vertex_set
         * A pointer to the beginnnig of the vertex set.
         */
        template<typename VertexType_, typename IndexTupleType_>
        CUDA_HOST_DEVICE static void inline set_coefficients(DataType (&coeffs)[img_dim][num_verts], const IndexTupleType_& index_tuple, const VertexType_* vertex_set)
        {
          typedef VertexType_ VertexType;

          const VertexType& v0 = vertex_set[index_tuple[0]];
          const VertexType& v1 = vertex_set[index_tuple[1]];
          const VertexType& v2 = vertex_set[index_tuple[2]];
          const VertexType& v3 = vertex_set[index_tuple[3]];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
          {
            coeffs[i][0] = DataType(0.25) * DataType( v0[i] + v1[i] + v2[i] + v3[i]);
            coeffs[i][1] = DataType(0.25) * DataType(-v0[i] + v1[i] - v2[i] + v3[i]);
            coeffs[i][2] = DataType(0.25) * DataType(-v0[i] - v1[i] + v2[i] + v3[i]);
            coeffs[i][3] = DataType(0.25) * DataType( v0[i] - v1[i] - v2[i] + v3[i]);
          }

        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void map_point(ImagePointType& img_point, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = coeffs[i][0] + coeffs[i][1] * dom_point[0] +
              (coeffs[i][2] + coeffs[i][3] * dom_point[0]) * dom_point[1];
          }
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1] + dom_point[1] * coeffs[i][3];
            jac_mat(i,1) = coeffs[i][2] + dom_point[0] * coeffs[i][3];
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& DOXY(dom_point), const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            hess_ten(i,0,0) = hess_ten(i,1,1) = DataType(0);
            hess_ten(i,0,1) = hess_ten(i,1,0) = coeffs[i][3];
          }
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The volume of the current cell.
         */
        CUDA_HOST_DEVICE static inline DataType volume(const DataType (&coeffs)[image_dim][num_verts])
        {
          // According to Varignon's theorem, the area/volume of a quadrilateral is
          // equal to twice the area of the dual parallelogram of the quadrilateral,
          // which is spanned by the four edge midpoints of the quadrilateral.
          // The Jacobian matrix of this transformation evaluated at the barycentre
          // of the reference element spans a parallelogram, which intersects with
          // our original quadrilateral in the edge midpoints and therefore (again
          // using Varignon's theorem) has the same area as the original quadrilateral.
          // Now the area of the "Jacobian parallelogram" is equal to four times
          // the determinant of its Jacobian determinant, which finally gives us a
          // formula for our quadrilateral area: 4*det(Jacobian(0,0)).
          // Note that this is not a lousy approximation, but a true identity :)

          // compute jacobian matrix at barycentre
          JacobianMatrixType jac_mat;
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
          }

          // return scaled volume
          return DataType(4) * jac_mat.vol();
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        CUDA_HOST_DEVICE static inline DataType width_directed(const ImagePointType& ray, const DataType (&coeffs)[image_dim][num_verts])
        {
          JacobianMatrixType jac_mat;
          JacobianInverseType jac_inv;
          DomainPointType ref_ray;

          // compute jacobian matrix at barycentre
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
          }

          // invert jacobian matrix and multiply by ray vector
          jac_inv.set_inverse(jac_mat);
          ref_ray.set_mat_vec_mult(jac_inv, ray);

          // return scaled inverse ray norm
          return DataType(2) / ref_ray.norm_euclid();
        }
      }; // EvalHelper<Hypercube<2>,...>

      template<typename DataType_, typename DomPointType_, typename ImgPointType_, typename JacMatType_, typename JacMatInvType_, typename JacDetType_,
               typename HessType_, typename HessInvType_, int img_dim>
      struct EvalHelper<DataType_, DomPointType_, ImgPointType_, JacMatType_, JacMatInvType_, JacDetType_, HessType_, HessInvType_, Shape::Hypercube<3>, img_dim>
      {
        /// evaluation data type
        typedef DataType_ DataType;
        /// domain point type
        typedef DomPointType_ DomainPointType;
        /// image point type
        typedef ImgPointType_ ImagePointType;
        /// jacobian matrix type
        typedef JacMatType_ JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef JacMatInvType_ JacobianInverseType;
        /// jacobian determinant type
        typedef JacDetType_ JacobianDeterminantType;
        /// hessian tensor type
        typedef HessType_ HessianTensorType;
        /// hessian inverse tensor type
        typedef HessInvType_ HessianInverseType;
        /// the shape type
        typedef Shape::Hypercube<3> ShapeType;

        static constexpr int num_verts = Shape::FaceTraits<ShapeType, 0>::count;
        static constexpr int image_dim = img_dim;

        /**
         * \brief Sets the coefficient vector based on the underlying standard transformation.
         *
         * \tparam VertexType_ The type of the vertex array.
         * \tparam IndexTupleType_ A type mapping the local index to the dofs. Has to support operator[](Index i).
         *
         * \param[out] coeffs
         * A two dimensional array holding the local coefficients of the trafo after the call.
         *
         * \param[in] index_tuple
         * A mapping from local indices to global.
         *
         * \param[in] vertex_set
         * A pointer to the beginnnig of the vertex set.
         */
        template<typename VertexType_, typename IndexTupleType_>
        CUDA_HOST_DEVICE static void inline set_coefficients(DataType (&coeffs)[img_dim][num_verts], const IndexTupleType_& index_tuple, const VertexType_* vertex_set)
        {
          typedef VertexType_ VertexType;

          const VertexType& v0 = vertex_set[index_tuple[0]];
          const VertexType& v1 = vertex_set[index_tuple[1]];
          const VertexType& v2 = vertex_set[index_tuple[2]];
          const VertexType& v3 = vertex_set[index_tuple[3]];
          const VertexType& v4 = vertex_set[index_tuple[4]];
          const VertexType& v5 = vertex_set[index_tuple[5]];
          const VertexType& v6 = vertex_set[index_tuple[6]];
          const VertexType& v7 = vertex_set[index_tuple[7]];

          // calculate transformation coefficients
          // j = _coeff[i][j] for all j = 0....7
          // v = 0 + 1*x + 2*y + 3*z + x*y*4 + x*z*5 + y*z*6 + x*y*z*7
          for(int i(0); i < image_dim; ++i)
          {
            coeffs[i][0] = DataType(0.125) * DataType( + v0[i] + v1[i] + v2[i] + v3[i] + v4[i] + v5[i] + v6[i] + v7[i]);
            coeffs[i][1] = DataType(0.125) * DataType( - v0[i] + v1[i] - v2[i] + v3[i] - v4[i] + v5[i] - v6[i] + v7[i]);
            coeffs[i][2] = DataType(0.125) * DataType( - v0[i] - v1[i] + v2[i] + v3[i] - v4[i] - v5[i] + v6[i] + v7[i]);
            coeffs[i][3] = DataType(0.125) * DataType( - v0[i] - v1[i] - v2[i] - v3[i] + v4[i] + v5[i] + v6[i] + v7[i]);
            coeffs[i][4] = DataType(0.125) * DataType( + v0[i] - v1[i] - v2[i] + v3[i] + v4[i] - v5[i] - v6[i] + v7[i]);
            coeffs[i][5] = DataType(0.125) * DataType( + v0[i] - v1[i] + v2[i] - v3[i] - v4[i] + v5[i] - v6[i] + v7[i]);
            coeffs[i][6] = DataType(0.125) * DataType( + v0[i] + v1[i] - v2[i] - v3[i] - v4[i] - v5[i] + v6[i] + v7[i]);
            coeffs[i][7] = DataType(0.125) * DataType( - v0[i] + v1[i] + v2[i] - v3[i] + v4[i] - v5[i] - v6[i] + v7[i]);
          }

        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void map_point(ImagePointType& img_point, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = coeffs[i][0]
              + dom_point[0] * (coeffs[i][1])
              + dom_point[1] * (coeffs[i][2] + dom_point[0]*coeffs[i][4])
              + dom_point[2] * (coeffs[i][3] + dom_point[0]*coeffs[i][5]
              + dom_point[1] * (coeffs[i][6] + dom_point[0]*coeffs[i][7]));
          }
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1] + dom_point[1] * coeffs[i][4] + dom_point[2] * (coeffs[i][5] + dom_point[1] * coeffs[i][7]);
            jac_mat(i,1) = coeffs[i][2] + dom_point[0] * coeffs[i][4] + dom_point[2] * (coeffs[i][6] + dom_point[0] * coeffs[i][7]);
            jac_mat(i,2) = coeffs[i][3] + dom_point[0] * coeffs[i][5] + dom_point[1] * (coeffs[i][6] + dom_point[0] * coeffs[i][7]);
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         */
        CUDA_HOST_DEVICE static inline void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& dom_point, const DataType (&coeffs)[image_dim][num_verts])
        {
          for(int i(0); i < image_dim; ++i)
          {
            hess_ten(i,0,0) = hess_ten(i,1,1) = hess_ten(i,2,2) = DataType(0);
            hess_ten(i,0,1) = hess_ten(i,1,0) = coeffs[i][4] + coeffs[i][7] * dom_point[2];
            hess_ten(i,0,2) = hess_ten(i,2,0) = coeffs[i][5] + coeffs[i][7] * dom_point[1];
            hess_ten(i,1,2) = hess_ten(i,2,1) = coeffs[i][6] + coeffs[i][7] * dom_point[0];
          }
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The volume of the current cell.
         */
        CUDA_HOST_DEVICE static inline DataType volume(const DataType (&coeffs)[image_dim][num_verts])
        {
          // In contrast to 2D, it is not sufficient to evaluate the Jacobian determinant
          // in the barycentre of the cell to compute the cell's volume, as this will give
          // very inaccurate results for non-parallelepiped cells. Instead, we have to use
          // the 2x2x2 Gauss-Legendre cubature rule to integrate the volume of the cell,
          // which is hard-coded in the code below.

          // compute 1D Gauss-Legendre root
          const DataType cx = DataType(FEAT_F128C(0.57735026918962576450914878050195745564760175127));

          JacobianMatrixType jac_mat;
          DomainPointType cub_pt;
          DataType vol = DataType(0);

          // loop over all 8 cubature points
          for(int i(0); i < 8; ++i)
          {
            // set cubature point coords by magic bitshifts
            for(int j(0); j < 3; ++j)
              cub_pt[j] = DataType((((i >> j) & 1) << 1) - 1) * cx;

            // compute jacobian matrix and add its volume
            calc_jac_mat(jac_mat, cub_pt, coeffs);
            vol += jac_mat.vol();
          }

          return vol;
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \param[in] coeffs
         * A two dimensional array holding the local coefficients of the trafo.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        CUDA_HOST_DEVICE static inline DataType width_directed(const ImagePointType& ray, const DataType (&coeffs)[image_dim][num_verts])
        {
          JacobianMatrixType jac_mat;
          JacobianInverseType jac_inv;
          DomainPointType ref_ray;

          // compute jacobian matrix at barycentre
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = coeffs[i][1];
            jac_mat(i,1) = coeffs[i][2];
            jac_mat(i,2) = coeffs[i][3];
          }

          // invert jacobian matrix and multiply by ray vector
          jac_inv.set_inverse(jac_mat);
          ref_ray.set_mat_vec_mult(jac_inv, ray);

          // return scaled inverse ray norm
          return DataType(2) / ref_ray.norm_euclid();
        }
      }; // EvalHelper<Hypercube<3>,...>

      // template<typename DataType_, typename DomPointType_, typename ImgPointType_, typename JacMatType_, typename JacMatInvType_, typename JacDetType_,
      //          typename HessType_, typename HessInvType_, int img_dim>
      // struct EvalHelper<DataType_, DomPointType_, ImgPointType_, JacMatType_, JacMatInvType_, JacDetType_, HessType_, HessInvType_, Shape::Vertex, img_dim>
      // {
      //   /// evaluation data type
      //   typedef DataType_ DataType;
      //   /// domain point type
      //   typedef DomPointType_ DomainPointType;
      //   /// image point type
      //   typedef ImgPointType_ ImagePointType;
      //   /// jacobian matrix type
      //   typedef JacMatType_ JacobianMatrixType;
      //   /// jacobian inverse matrix type
      //   typedef JacMatInvType_ JacobianInverseType;
      //   /// jacobian determinant type
      //   typedef JacDetType_ JacobianDeterminantType;
      //   /// hessian tensor type
      //   typedef HessType_ HessianTensorType;
      //   /// hessian inverse tensor type
      //   typedef HessInvType_ HessianInverseType;
      //   /// the shape type
      //   typedef Shape::Vertex ShapeType;

      //   static constexpr int num_verts = Shape::FaceTraits<ShapeType, 0>::count;
      //   static constexpr int image_dim = img_dim;

      //   template<typename VertexSetType_, typename IndexTupleType_>
      //   CUDA_HOST_DEVICE static void inline set_coefficients(DataType (&coeffs)[image_dim][num_verts], const IndexTupleType_& DOXY(index_tuple), const VertexSetType_& vertex_set)
      //   {
      //     const auto& vtx = vertex_set;

      //     for(int i(0); i < image_dim; ++i)
      //     {
      //       coeffs[i][0] = DataType(vtx[i]);
      //     }

      //   }

      //   CUDA_HOST_DEVICE static inline void map_point(ImagePointType& img_point, const DomainPointType& DOXY(dom_point), const DataType (&coeffs)[image_dim][num_verts])
      //   {
      //     for(int i(0); i < image_dim; ++i)
      //     {
      //       img_point[i] = coeffs[i][0];
      //     }
      //   }

      //   CUDA_HOST_DEVICE static inline void calc_jac_mat(JacobianMatrixType& DOXY(jac_mat), const DomainPointType& DOXY(dom_point), const DataType DOXY(coeffs)[image_dim][num_verts])
      //   {
      //   }

      //   CUDA_HOST_DEVICE static inline void calc_hess_ten(HessianTensorType& DOXY(hess_ten), const DomainPointType& DOXY(dom_point), const DataType DOXY(coeffs)[image_dim][num_verts])
      //   {
      //   }

      //   CUDA_HOST_DEVICE static inline DataType volume(const DataType DOXY(coeffs)[image_dim][num_verts])
      //   {
      //   }

      //   CUDA_HOST_DEVICE static inline DataType width_directed(const ImagePointType& DOXY(ray), const DataType DOXY(coeffs)[image_dim][num_verts])
      //   {
      //   }
      // };
    }
  }
}


#endif
