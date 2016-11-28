#pragma once
#ifndef KERNEL_TRAFO_STANDARD_INVERSE_MAPPING
#define KERNEL_TRAFO_STANDARD_INVERSE_MAPPING 1
#include <kernel/base_header.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  namespace Trafo
  {
    namespace Standard
    {
      /**
       * \brief Computes the inverse coordinate mapping wrt. a full-dimensional Simplex
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam world_dim
       * Dimension of the points to compute the inverse mapping for
       *
       * \tparam sc_
       * Stride for the vector of coefficients
       *
       * \tparam sp_
       * Stride for the point vector
       *
       * \tparam smx_, snx_
       * Row and column strides for the matrix holding the vertices
       *
       * \param[out] coeffs
       * The coefficients for the mapping from the reference cell to the real cell
       *
       * \param[in] point
       * The point to compute the coefficients for
       *
       * \param[in] x
       * Coordinates of the simplex on which we compute \c coeffs.
       *
       * Assume we have a non-degenerate \c Simplex<d> called \f$ S \f$ in \f$ \mathbb{R}^d \f$ defined by vertices
       * \f$ x^j \in \mathbb{R}^d, j = 0, \dots, d \f$.
       *
       * Then
       * \f[
       *   \forall x \in \mathbb{R}^d: x = x^0 + \sum_{j=1}^s \lambda_j (x^j - x^0)
       * \f]
       * and \f$ \lambda_j, j=1,\dots,s \f$ with \f$ \lambda_0 := 1 - sum_{j=1}^s \f$ are called the <b> barycentric
       * coordinates of \f$ x \f$ wrt. S </b>. Since \f$ \lambda_0 \f$ is redundant, we just compute the
       * \c coefficients \f$ \lambda_1, \dots, \lambda_s \f$. Note that the above equation can be rewritten as
       * \f[
       *   \forall x \in \mathbb{R}^d: x = \sum_{j=0}^s \lambda_j x^j.
       * \f]
       *
       * It is easy to see that
       * \f[
       *   x \in S \Leftrightarrow  \forall j = 0, \dots, d: \lambda_j \in [0, 1] .
       * \f]
       * If \f$ \exists j \in \{ 0, \dots, s \}: \lambda_j < 0 \f$, then \f$ x \notin S \f$ and \f$ x \f$ lies on the
       * far side of the plane defined by the facet opposite vertex \f$ j \f$. This makes the barycentric coordinates
       * very handy for finding out in which direction of a given simplex a point lies.
       *
       * \author Jordi Paul
       */
      template<typename DT_, int world_dim, int sc_, int sp_, int smx_, int snx_>
      void inverse_mapping(
        Tiny::Vector<DT_, world_dim, sc_>& coeffs,
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

      /**
       * \brief Computes the inverse coordinate mapping wrt. a Simplex<1> embedded in 2d or 3d
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam world_dim
       * Dimension of the points to compute the inverse mapping for
       *
       * \tparam sc_
       * Stride for the vector of coefficients
       *
       * \tparam sp_
       * Stride for the point vector
       *
       * \tparam smx_, snx_
       * Row and column strides for the matrix holding the vertices
       *
       * \param[out] coeffs
       * The coefficients for the mapping from the reference cell to the real cell
       *
       * \param[in] point
       * The point to compute the coefficients for
       *
       * \param[in] x
       * Coordinates of the simplex on which we compute \c coeffs.
       *
       * Assume we have a non-degenerate \c Simplex<1> called \f$ S \f$ in \f$ \mathbb{R}^d \f$ where \f$ d=2, 3 \f$
       * Then \f$ S \f$ is defined by vertices \f$ x^0, x^1 \in \mathbb{R}^d \f$.
       *
       * For a point \f$ x \in \mathbb{R}^d \f$ we are looking for its projection \f$ p \f$ onto the straight line
       * defined by \f$ x^0, x^1 \f$, which means that \f$ (x^1 - x^0, x - p) = 0 \f$. If this holds, \f$ p \f$ has
       * the form \f$ p = x^0 + \omega (x^1 - x^0), \omega \in \mathbb{R} \f$.
       *
       * It is easy to see that
       * \f[
       *   x \in S \Leftrightarrow \omega \in [0, 1]
       * \f]
       * and that
       * \f[
       *   \omega = \frac{(x^1 - x^0, x - x^0)}{\| x^1 - x^0 \|_2^2}.
       * \f]
       *
       * This routine computes the coefficient \f$ \omega \f$ and saves it to \f$ \mathrm{coeffs[0]} \f$ and the distance
       * \f$ \mathrm{coeffs[1]} = \| x - p \|_2 \f$.
       *
       * \author Jordi Paul
       */
      template<typename DT_, int world_dim, int sc_, int sp_, int smx_, int snx_>
      void inverse_mapping(
        Tiny::Vector<DT_, 2, sc_>& coeffs,
        const Tiny::Vector<DT_, world_dim, sp_>& point,
        const Tiny::Matrix<DT_, 2, world_dim, smx_, snx_>& x)
      {
        static_assert( (world_dim == 2 || world_dim == 3),
        "world dim has to be 2 or 3 for complementary barycentric coordinates");

        auto tmp = x[1]-x[0];
        DT_ sp(Tiny::dot(point - x[0],tmp));
        DT_ nsqr(Math::sqr(tmp.norm_euclid()));
        // This is omega
        coeffs[0] = sp/nsqr;
        tmp = point - (x[0] + coeffs[0]*(x[1]-x[0]));
        // This is the distance of point to the straight line defined by x[0], x[1]
        coeffs[1] = tmp.norm_euclid();
      }

      /**
       * \brief Computes the inverse coordinate mapping wrt. a Simplex<2> in 3d
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam world_dim
       * Dimension of the points to compute the inverse mapping for
       *
       * \tparam sc_
       * Stride for the vector of coefficients
       *
       * \tparam sp_
       * Stride for the point vector
       *
       * \tparam smx_, snx_
       * Row and column strides for the matrix holding the vertices
       *
       * \param[out] coeffs
       * The coefficients for the mapping from the reference cell to the real cell
       *
       * \param[in] point
       * The point to compute the coefficients for
       *
       * \param[in] x
       * Coordinates of the simplex on which we compute \c coeffs.
       *
       *
       * Assume we have a non-degenerate \c Simplex<2> called \f$ S \f$ in \f$ \mathbb{R}^d \f$ defined by vertices
       * \f$ x^j \in \mathbb{R}^d, j = 0, \dots, 2 \f$. Then
       * \f[
       *   v^3 := (x^1 - x^0) \times (x^2 - x^0) \Rightarrow v^3 \perp S.
       * \f]
       *
       * Using \f$ x^3 := \| v^3 \|_2^{-1} v^3 \f$, the vertices \f$ x^0, \dots, x^3 \f$ define a fictitious
       * \c Simplex<3> \f$ S' \f$. Then we can proceed with computing the coefficients as in the \c Simplex<d> in
       * \f$ \mathbb{R}^d\f$ variant of this function, and the last coefficient is just the distance of the point
       * \f$ x \f$ to the plane in which \f$ S \f$ lies.
       *
       * If \f$ \exists j \in \{ 0, \dots, s \}: \lambda_j < 0 \f$, then \f$ x \notin S \f$ and \f$ x \f$ lies on the
       * far side of the plane defined by the facet opposite vertex \f$ j \f$. This makes the barycentric coordinates
       * very handy for finding out in which direction of a given simplex a point lies.
       *
       * \author Jordi Paul
       */
      template<typename DT_, int sc_, int sp_, int smx_, int snx_>
      void inverse_mapping(
        Tiny::Vector<DT_, 3, sc_>& coeffs,
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
} // namespace FEAT
#endif //KERNEL_TRAFO_STANDARD_INVERSE_MAPPING
