#pragma once
#ifndef FIBER_TENSOR_OPERATIONS_HPP
#define FIBER_TENSOR_OPERATIONS_HPP 1

//includes, FEAT
#include <kernel/util/tiny_algebra.hpp>



namespace FEAT
{
  namespace CCND_FIBER
  {
#ifdef DOXYGEN
    /**
     * \brief Computes the outer contraction of two m_ sized vectors over a hyper-symmetric (m_)^4 fourth order tensor resulting in a m_ x m_ matrix.
     *        Due to its symmetry, the number of operations can be reduced considerably, which in turn needs specialization for each dimension.
     *
     * \tparam DT_
     * Data type.
     *
     * \tparam m_
     * Number of rows of the matrix.
     *
     * \tparam snv_
     * Stride of the first vector.
     *
     * \tparam snw_
     * Stride of the second vector.
     *
     * \tparam slt_
     * Stride of the fourth order moment vector.
     *
     * \param[in] x
     * A m_ sized vector.
     *
     * \param[in] y
     * A m_ sized vector.
     *
     * \param[in] t
     * A m_ * (m_ + 1) * (m_ + 2) * (m_ + 3) / 24 sized vector representing the tensor entries of the fourth moment tensor.
     *
     * \param[in] alpha
     * A scaling parameter.
     *
     * \returns
     * A symmetric m_ times m_ matrix.
     *
     * \note The order of x and y does not matter due to the symmetry of the tensor.
     *
     * \note So far, this is only implemented for m_ = 2,3 due to the specific ordering needed for our tensors Aijkl, which are defined as follows:
     *             2D                                           |                            3D
     *        t(0) = A1111                                      |                       t(0)  = A1111
     *        t(1) = A2222                                      |                       t(1)  = A2222
     *        t(2) = A1112                                      |                       t(2)  = A3333
     *        t(3) = A2221                                      |                       t(3)  = A1112
     *        t(4) = A1122                                      |                       t(4)  = A1113
     *                                                          |                       t(5)  = A2221
     *                                                          |                       t(6)  = A2223
     *                                                          |                       t(7)  = A3331
     *                                                          |                       t(8)  = A3332
     *                                                          |                       t(9)  = A1122
     *                                                          |                       t(10) = A1133
     *                                                          |                       t(11) = A2233
     *                                                          |                       t(12) = A1123
     *                                                          |                       t(13) = A2213
     *                                                          |                       t(14) = A3312
     */
    template<typename DT_, int m_, int snv_, int snw_, int slt_>
    Tiny::Matrix<DT_, m_, m_> fourth_order_contraction(const Tiny::Vector<DT_, m_, snv_>& x, const Tiny::Vector<DT_, m_, snw_>& y, const Tiny::Vector<DT_, m_ * (m_ + 1) * (m_ + 2) * (m_ + 3) / 24, slt_>& t, typename Tiny::Matrix<DT_, m_, m_>::DataType alpha);
#endif

    /// \cond internal
    template<typename DT_, int snv_, int snw_, int slt_>
    Tiny::Matrix<DT_, 2, 2> fourth_order_contraction(const Tiny::Vector<DT_, 2, snv_>& x, const Tiny::Vector<DT_, 2, snw_>& y, const Tiny::Vector<DT_, 5, slt_>& t, typename Tiny::Matrix<DT_, 2, 2>::DataType alpha)
    {
      Tiny::Matrix<DT_, 2, 2> nu(DT_(0));

      typedef typename Tiny::Matrix<DT_, 2, 2>::DataType DataType;
      //do some calulations beforehand
      const DataType tmp00 = x[0] * y[0];
      const DataType tmp11 = x[1] * y[1];
      const DataType tmp01 = x[1] * y[0] + x[0] * y[1];

      nu[0][0] = alpha * (tmp00 * t[0] + tmp11 * t[4] + t[2] * tmp01);
      nu[0][1] = alpha * (tmp00 * t[2] + tmp11 * t[3] + t[4] * tmp01);
      nu[1][0] = alpha * (tmp00 * t[2] + tmp11 * t[3] + t[4] * tmp01);
      nu[1][1] = alpha * (tmp00 * t[4] + tmp11 * t[1] + t[3] * tmp01);

      return nu;
    }

    template<typename DT_, int snv_, int snw_, int slt_>
    Tiny::Matrix<DT_, 3, 3> fourth_order_contraction(const Tiny::Vector<DT_, 3, snv_>& x, const Tiny::Vector<DT_, 3, snw_>& y, const Tiny::Vector<DT_, 15, slt_>& t, typename Tiny::Matrix<DT_, 3, 3>::DataType alpha)
    {
      Tiny::Matrix<DT_, 3, 3> nu(DT_(0));

      typedef typename Tiny::Matrix<DT_, 3, 3>::DataType DataType;
      //do some calulations beforehand
      const DataType tmp00 = x[0] * y[0];
      const DataType tmp11 = x[1] * y[1];
      const DataType tmp22 = x[2] * y[2];
      const DataType tmp01 = x[0] * y[1] + x[1] * y[0];
      const DataType tmp02 = x[0] * y[2] + x[2] * y[0];
      const DataType tmp12 = x[1] * y[2] + x[2] * y[1];

      nu[0][0] = alpha * (t[0] * tmp00 + t[3] * tmp01 + t[4] * tmp02 + t[9] * tmp11 + t[12] * tmp12 + t[10] * tmp22);
      nu[0][1] = alpha * (t[3] * tmp00 + t[9] * tmp01 + t[12] * tmp02 + t[5] * tmp11 + t[13] * tmp12 + t[14] * tmp22);
      nu[1][0] = alpha * (t[3] * tmp00 + t[9] * tmp01 + t[12] * tmp02 + t[5] * tmp11 + t[13] * tmp12 + t[14] * tmp22);
      nu[0][2] = alpha * (t[4] * tmp00 + t[12] * tmp01 + t[10] * tmp02 + t[13] * tmp11 + t[14] * tmp12 + t[7] * tmp22);
      nu[2][0] = alpha * (t[4] * tmp00 + t[12] * tmp01 + t[10] * tmp02 + t[13] * tmp11 + t[14] * tmp12 + t[7] * tmp22);
      nu[1][1] = alpha * (t[9] * tmp00 + t[5] * tmp01 + t[13] * tmp02 + t[1] * tmp11 + t[6] * tmp12 + t[11] * tmp22);
      nu[1][2] = alpha * (t[12] * tmp00 + t[13] * tmp01 + t[14] * tmp02 + t[6] * tmp11 + t[11] * tmp12 + t[8] * tmp22);
      nu[2][1] = alpha * (t[12] * tmp00 + t[13] * tmp01 + t[14] * tmp02 + t[6] * tmp11 + t[11] * tmp12 + t[8] * tmp22);
      nu[2][2] = alpha * (t[10] * tmp00 + t[14] * tmp01 + t[7] * tmp02 + t[11] * tmp11 + t[8] * tmp12 + t[2] * tmp22);

      return nu;
    }
    /// \endcond


    /**
     * \brief Adds the given tiny matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
     *
     * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a m-by-n-by-k-by-l input tensor,
     * respresented by a m * n * k * l sized input vector \e t, where T_{i,j,p,q} = t(i*(n*k*l) + j*(k*l) + p*l + q), with i = 0,...,m-1 ; j = 0,...,n-1 ;
     * p = 0,...,k-1 ; q = 0,...,l-1.
     * This operation computes:
     * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{i,j, p,q} \cdot w_q \f]
     *
     * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
     *
     * \param[out] A
     * The m- by n-sized matrix, on which to be added.
     *
     * \param[in] x
     * The k-size vector, which is the left multiplicand.
     *
     * \param[in] y
     * The l-sized vector, which is the right muliplicand.
     *
     * \param[in] t
     * The m*n*k*l sized vector representing the tensor T.
     *
     * \param[in] alpha
     * A scaling factor for the product.
     *
     */
    template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
    void add_tensor4_outer_product_contraction_34(
      Tiny::Matrix<DT_, m_, n_>& A,
      const Tiny::Vector<DT_, k_, snv_>& x,
      const Tiny::Vector<DT_, l_, snw_>& y,
      const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
      DT_ alpha = DT_(1))
    {
      const int dim_size_3 = n_ * k_ * l_;
      const int dim_size_2 = k_ * l_;
      const int dim_size_1 = l_;
      for(int i(0); i < m_; ++i)
      {
        for(int j(0); j < n_; ++j)
        {
          DT_ r(0);
          for(int p(0); p < k_; ++p)
          {
            for(int q(0); q < l_; ++q)
            {
              const int mapping = i * dim_size_3 + j * dim_size_2 + p * dim_size_1 + q;
              r += x(p) * t(mapping) * y(q);
            }
          }
          A[i][j] += alpha * r;
        }
      }

      return;
    }

    /**
       * \brief Adds this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a m-by-k-by-n-by-l input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{i,p,j,q} = t(i*(k*n*l) + p*(n*l) + j*l + q), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_k \cdot T_{i,p,j,q} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand.
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand.
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      void add_tensor4_outer_product_contraction_24(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        const int dim_size_3 = k_ * n_ * l_;
        const int dim_size_2 = n_ * l_;
        const int dim_size_1 = l_;
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DT_ r(0);
            for(int p(0); p < k_; ++p)
            {
              for(int q(0); q < l_; ++q)
              {
                const int mapping = i * dim_size_3 + p * dim_size_2 + j * dim_size_1 + q;
                r += x(p) * t(mapping) * y(q);
              }
            }
            A[i][j] += alpha * r;
          }
        }

        return;
      }

      /**
       * \brief Adds this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a m-by-k-by-l-by-n input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{i,p,q,j} = t(i*(k*l*n) + p*(l*n) + q*n + j), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{i,p,q,j} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      void add_tensor4_outer_product_contraction_23(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        const int dim_size_3 = k_ * l_ * n_;
        const int dim_size_2 = l_ * n_;
        const int dim_size_1 = n_;
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DT_ r(0);
            for(int p(0); p < k_; ++p)
            {
              for(int q(0); q < l_; ++q)
              {
                const int mapping = i * dim_size_3 + p * dim_size_2 + q * dim_size_1 + j;
                r += x(p) * t(mapping) * y(q);
              }
            }
            A[i][j] += alpha * r;
          }
        }

        return;
      }

      /**
       * \brief Adds this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a k-by-m-by-n-by-l input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{p,i,j,q} = t(p*(m*n*l) + i*(n*l) + j*l + q), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{p,i,j,q} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand.
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand.
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      void add_tensor4_outer_product_contraction_14(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        const int dim_size_3 = m_ * n_ * l_;
        const int dim_size_2 = n_ * l_;
        const int dim_size_1 = l_;
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DT_ r(0);
            for(int p(0); p < k_; ++p)
            {
              for(int q(0); q < l_; ++q)
              {
                const int mapping = p * dim_size_3 + i * dim_size_2 + j * dim_size_1 + q;
                r += x(p) * t(mapping) * y(q);
              }
            }
            A[i][j] += alpha * r;
          }
        }

        return;
      }

      /**
       * \brief Adds this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a k-by-m-by-l-by-n input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{p,i,q,j} = t(p*(m*l*n) + i*(l*n) + q*n + j), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{p,i,q,j} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      void add_tensor4_outer_product_contraction_13(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        const int dim_size_3 = m_ * l_ * n_;
        const int dim_size_2 = l_ * n_;
        const int dim_size_1 = n_;
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DT_ r(0);
            for(int p(0); p < k_; ++p)
            {
              for(int q(0); q < l_; ++q)
              {
                const int mapping = p * dim_size_3 + i * dim_size_2 + q * dim_size_1 + j;
                r += x(p) * t(mapping) * y(q);
              }
            }
            A[i][j] += alpha * r;
          }
        }

        return;
      }

      /**
       * \brief Adds this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a k-by-l-by-m-by-n input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{p,l,i,j} = t(p*(l*m*n) + q*(m*n) + i*n + j), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{p,q,i,j} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand.
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand.
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      void add_tensor4_outer_product_contraction_12(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        const int dim_size_3 = l_ * m_ * n_;
        const int dim_size_2 = m_ * n_;
        const int dim_size_1 = n_;
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DT_ r(0);
            for(int p(0); p < k_; ++p)
            {
              for(int q(0); q < l_; ++q)
              {
                const int mapping = p * dim_size_3 + q * dim_size_2 + i * dim_size_1 + j;
                r += x(p) * t(mapping) * y(q);
              }
            }
            A[i][j] += alpha * r;
          }
        }

        return;
      }

      /**
       * \brief Adds this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors and the tensor is symmetric
       *
       * Let \e A denote a m-by-m matrix, \e v the m-sized and \e w the m-sized input vectors an \e T a m-by-m-by-m-by-m input tensor,
       * respresented by a m*(m+1)*(m+2)*(m+3)/24 sized input vector \e t, due to how the tensor is mapped, we will only specialize this for m = 3 dimensions for now..
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{p,q,i,j} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. Since the tensor is symmetric, there is no difference over which entries we reduce exactly
       *
       * \note The tensor is mapped to a 15 dimension vector, where the entries which differ from one another are represented in the following way:
       *                              t(0)  = A1111
       *                              t(1)  = A2222
       *                              t(2)  = A3333
       *                              t(3)  = A1112
       *                              t(4)  = A1113
       *                              t(5)  = A2221
       *                              t(6)  = A2223
       *                              t(7)  = A3331
       *                              t(8)  = A3332
       *                              t(9)  = A1122
       *                              t(10) = A1133
       *                              t(11) = A2233
       *                              t(12) = A1123
       *                              t(13) = A2213
       *                              t(14) = A3312
       *
       * \param[out] A
       * The m- by m-sized matrix, on which to be added.
       *
       * \param[in] x
       * The m-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The m-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*(m+1)*(m+2)*(m+3)/24 sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int snv_, int snw_, int slt_>
      void add_tensor4_outer_product_contraction_symmetric(
        Tiny::Matrix<DT_, m_, m_>& A,
        const Tiny::Vector<DT_, m_, snv_>& x,
        const Tiny::Vector<DT_, m_, snw_>& y,
        const Tiny::Vector<DT_, m_ * (m_ + 1) * (m_ + 2) * (m_ + 3) / 24, slt_>& t,
        DT_ alpha = DT_(1))
      {
        static_assert(((m_ == 3) || (m_ == 2)), "this function works only for 2x2 and 3x3 matrices");

        A += fourth_order_contraction(x, y, t, alpha);

        return;
      }

      /**
       * \brief Sets this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a m-by-n-by-k-by-l input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{i,j,p,q} = t(i*(n*k*l) + j*(k*l) + p*l + q), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{i,j, p,q} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand.
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand.
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \return \p *this
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      inline void set_tensor4_outer_product_contraction_34(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        A.format();
        add_tensor4_outer_product_contraction_34(A, x, y, t, alpha);
      }

      /**
       * \brief Sets this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a m-by-k-by-n-by-l input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{i,p,j,q} = t(i*(k*n*l) + p*(n*l) + j*l + q), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_k \cdot T_{i,p,j,q} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      inline void set_tensor4_outer_product_contraction_24(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        A.format();
        add_tensor4_outer_product_contraction_24(A, x, y, t, alpha);
      }

      /**
       * \brief Sets this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a m-by-k-by-l-by-n input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{i,p,q,j} = t(i*(k*l*n) + p*(l*n) + q*n + j), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{i,p,q,j} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      inline void set_tensor4_outer_product_contraction_23(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        A.format();
        add_tensor4_outer_product_contraction_23(A, x, y, t, alpha);
      }

      /**
       * \brief Sets this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a k-by-m-by-n-by-l input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{p,i,j,q} = t(p*(m*n*l) + i*(n*l) + j*l + q), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{p,i,j,q} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      inline void set_tensor4_outer_product_contraction_14(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        A.format();
        add_tensor4_outer_product_contraction_14(A, x, y, t, alpha);
      }

      /**
       * \brief Sets this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a k-by-m-by-l-by-n input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{p,i,q,j} = t(p*(m*l*n) + i*(l*n) + q*n + j), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{p,i,q,j} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      inline void set_tensor4_outer_product_contraction_13(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        A.format();
        add_tensor4_outer_product_contraction_13(A, x, y, t, alpha);
      }

      /**
       * \brief Sets this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors
       *
       * Let \e A denote a m-by-n matrix, \e v the k-sized and \e w the l-sized input vectors an \e T a k-by-l-by-m-by-n input tensor,
       * respresented by a m * n * k * l sized input vector \e t, where T_{p,l,i,j} = t(p*(l*m*n) + q*(m*n) + i*n + j), with i = 0,...,m-1 ; j = 0,...,n-1 ;
       * p = 0,...,k-1 ; q = 0,...,l-1.
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{p,q,i,j} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \param[out] A
       * The m- by n-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int n_, int k_, int l_, int snv_, int snw_, int slt_>
      inline void set_tensor4_outer_product_contraction_12(
        Tiny::Matrix<DT_, m_, n_>& A,
        const Tiny::Vector<DT_, k_, snv_>& x,
        const Tiny::Vector<DT_, l_, snw_>& y,
        const Tiny::Vector<DT_, m_ * n_ * k_ * l_, slt_>& t,
        DT_ alpha = DT_(1))
      {
        A.format();
        add_tensor4_outer_product_contraction_12(A, x, y, t, alpha);
      }

      /**
       * \brief Sets this matrix to the result of a tensor4-matrix contraction, where the matrix is given as an outer product of two vectors and the tensor is symmetric
       *
       * Let \e A denote a m-by-m matrix, \e v the m-sized and \e w the m-sized input vectors an \e T a m-by-m-by-m-by-m input tensor,
       * respresented by a m*(m+1)*(m+2)*(m+3)/24 sized input vector \e t, due to how the tensor is mapped, we will only specialize this for m = 3 dimensions for now...
       * This operation computes:
       * \f[ \forall i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \sum_{p=0}^{k-1}\sum_{q=0}^{l-1} v_p \cdot T_{p,q,i,j} \cdot w_q \f]
       *
       * \note this function is used in the assembly of tensor-diffusion formulations. There are other variants, where we sum regarding other indizes.
       *
       * \note see the add method for how the tensor is mapped
       *
       * \param[out] A
       * The m- by m-sized matrix, on which to be added.
       *
       * \param[in] x
       * The k-size vector, which is the left multiplicand
       *
       * \param[in] y
       * The l-sized vector, which is the right muliplicand
       *
       * \param[in] t
       * The m*n*k*l sized vector representing the tensor T.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       */
      template<typename DT_, int m_, int snv_, int snw_, int slt_>
      inline void set_tensor4_outer_product_contraction_symmetric(
        Tiny::Matrix<DT_, m_, m_>& A,
        const Tiny::Vector<DT_, m_, snv_>& x,
        const Tiny::Vector<DT_, m_, snw_>& y,
        const Tiny::Vector<DT_, m_ * (m_ + 1) * (m_ + 2) * (m_ + 3) / 24, slt_>& t,
        DT_ alpha = DT_(1))
      {
        static_assert(((m_ == 3) || (m_ == 2)), "this function works only for 2x2 and 3x3 matrices");
        A.format();
        add_tensor4_outer_product_contraction_symmetric(A, x, y, t, alpha);
      }


  }
}
#endif