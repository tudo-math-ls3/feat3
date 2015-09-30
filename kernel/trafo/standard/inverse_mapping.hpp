#pragma once
#ifndef KERNEL_TRAFO_STANDARD_INVERSE_MAPPING
#define KERNEL_TRAFO_STANDARD_INVERSE_MAPPING 1
#include <kernel/base_header.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAST
{
  namespace Trafo
  {
    namespace Standard
    {
      /**
       * \brief Computes barycentric coordinates wrt. a Simplex
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam world_dim
       * Dimension of the points to compute the barycentric coordinates for
       *
       * \tparam sb_
       * Stride for the vector of barycentric coordinates
       *
       * \tparam sp_
       * Stride for the point vector
       *
       * \tparam smx_, snx_
       * Row and column strides for the matrix holding the vertices
       *
       * \param[out] bary
       * Vector of barycentric coordinates
       *
       * \param[in] point
       * The point to compute the barycentric coordinates for
       *
       * \param[in] coords
       * Coordinates of the simplex on which we compute \c bary
       *
       * Assume we have a non-degenerate \c Simplex<s> called \f$ S \f$ in \f$ \mathbb{R}^d \f$ where either
       * \f$ s = d-1 \f$ or \f$ s = d \f$. Then \f$ S \f$ is defined by vertices
       * \f$ x^j \in \mathbb{R}^d, j = 0,\dots,s \f$.
       *
       * Assume first that \f$ s = d \f$. Then
       * \f[
       *   \forall x \in \mathbb{R}^d: x = x^0 + \sum_{j=1}^s \lambda_j (x^j - x^0)
       * \f]
       * and \f$ \lambda_j, j=0,\dots,s \f$ with \f$ \lambda_0 := 1 - sum_{j=1}^s \f$ are called the <b> barycentric
       * coordinates of \f$ x \f$ wrt. S </b>. The above equation can be rewritten as
       * \f[
       *   \forall x \in \mathbb{R}^d: x = \sum_{j=0}^s \lambda_j x^j.
       * \f]
       *
       * It is easy to see that
       * \f[
       *   x \in S \Leftrightarrow  \forall j = 1, \dots, s: \lambda_j \in [0, 1] .
       * \f]
       * If \f$ \exists j \in \{ 0, \dots, s \}: \lambda_j < 0 \f$, then \f$ x \notin S \f$ and \f$ x \f$ lies on the
       * far side of the plane defined by the facet opposite vertex \f$ j \f$. This makes the barycentric coordinates
       * very handy for finding out in which direction of a given simplex a point lies.
       *
       * This can even be done if the simplex is of co-dimension 1, i.e. a Simplex<2> in 3d. Assume now that
       * \f$ s = d-1 \f$. Then define the \f$d\f$-dimensional simplex \$ \hat{S} \$ by the vertices
       * \f$ x^j, j=0, \dots, s+1 \f$, where \f$ x^{s+1} \f$ is defined by
       * \f[
       *    \forall j=0, \dots, s: x^{s+1}-x^0 \perp x^j - x^0, \| x^{s+1} - x^0\|_2 = \left(\Pi_{j=0}^s \| x^j - x^0 \|_2
       *    \right)^{\frac{1}{d}-1}
       * \f]
       *
       * The scaling of \f$ x^{s+1} - x^0 \f$ makes sure that the length of the additional edge in the "ficticious"
       * simplex \f$ \widehat{S} \f$ is similar to the lengths of the edges
       *
       */
      template<typename DT_, int world_dim, int sb_, int sp_, int smx_, int snx_>
      void inverse_mapping(
        Tiny::Vector<DT_, world_dim, sb_>& coeffs,
        const Tiny::Vector<DT_, world_dim, sp_>& point,
        const Tiny::Matrix<DT_, world_dim+1, world_dim, smx_, snx_>& x)
        {
          // A will contain the transformation matrix
          Tiny::Matrix<DT_, world_dim, world_dim> A(DT_(0));
          for(int i(0); i < world_dim; ++i)
          {
            for(int j(0); j < world_dim; ++j)
            {
              A(i,j) = x(j+1,i) - x(0,i);
            }
          }

          // This will be the inverse of the transformation matrix
          Tiny::Matrix<DT_, world_dim, world_dim> Ainv(DT_(0));
          Ainv.set_inverse(A);

          // This is the solution of A u = point - x[0]
          coeffs = Ainv*(point - x[0]);

        }

      template<typename DT_, int world_dim, int sb_, int sp_, int smx_, int snx_>
      void inverse_mapping(
        Tiny::Vector<DT_, 2, sb_>& coeffs,
        const Tiny::Vector<DT_, world_dim, sp_>& point,
        const Tiny::Matrix<DT_, 2, world_dim, smx_, snx_>& x)
      {
        static_assert( (world_dim == 2 || world_dim == 3),
        "world dim has to be 2 or 3 for complementary barycentric coordinates");

        auto tmp = x[1]-x[0];
        DT_ sp(Tiny::dot(point - x[0],tmp));
        DT_ nsqr(Math::sqr(tmp.norm_euclid()));
        coeffs[0] = sp/nsqr;
        tmp = point - (x[0] + coeffs[0]*(x[1]-x[0]));
        coeffs[1] = tmp.norm_euclid();
      }

      template<typename DT_, int sb_, int sp_, int smx_, int snx_>
      void inverse_mapping(
        Tiny::Vector<DT_, 3, sb_>& coeffs,
        const Tiny::Vector<DT_, 3, sp_>& point,
        const Tiny::Matrix<DT_, 3, 3, smx_, snx_>& x)
        {
          static constexpr int world_dim = 3;
          static constexpr int shape_dim = world_dim-1;

          // A will contain the transformation matrix
          Tiny::Matrix<DT_, world_dim, world_dim> A(DT_(0));
          // Fill the all rows in the first shape_dim = world_dim-1 columns
          for(int i(0); i < world_dim; ++i)
          {
            for(int j(0); j < shape_dim; ++j)
              A(i,j) = x(j+1,i) - x(0,i);
          }

          // The last column is the additional direction for our augmented simplex and it is orthogonal to the rest
          Tiny::Vector<DT_, world_dim> ortho = Tiny::orthogonal(A.template size_cast<world_dim, shape_dim>());

          // In 3d, the 2-norm of ortho is the 2d volume of the parallelogram defined by A[0], A[1]. That makes this
          // axis badly scaled if the parallelogram is either very small or very large. So we rescale ortho to unity
          // so the resulting coefficient is just the distance
          ortho.normalise();

          // Set the last column in A
          for(int i(0); i < world_dim; ++i)
            A(i,world_dim-1) = ortho(i);

          // This will be the inverse of the transformation matrix
          Tiny::Matrix<DT_, world_dim, world_dim> Ainv(DT_(0));
          Ainv.set_inverse(A);

          // This is the solution of A u = point - x[0]
          coeffs = Ainv*(point - x[0]);
        }

    } // namespace Standard
  } // namespace Trafo
} // namespace FEAST
#endif //KERNEL_TRAFO_STANDARD_INVERSE_MAPPING
