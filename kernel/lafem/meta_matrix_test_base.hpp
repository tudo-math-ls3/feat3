#pragma once
#ifndef KERNEL_LAFEM_META_MATRIX_TEST_BASE_HPP
#define KERNEL_LAFEM_META_MATRIX_TEST_BASE_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/power_diag_matrix.hpp>
#include <kernel/lafem/power_col_matrix.hpp>
#include <kernel/lafem/power_row_matrix.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /// \cond internal
    // Helper class: selects the matrix and vector types
    template<typename AlgoType_, typename DataType_>
    struct MetaMatrixTestHelper
    {
      /// scalar vector type
      typedef DenseVector<typename AlgoType_::MemType, DataType_> ScalarVector;
      /// scalar matrix type A
      typedef SparseMatrixCSR<typename AlgoType_::MemType, DataType_> ScalarMatrixA;
      /// scalar matrix type B
      typedef SparseMatrixELL<typename AlgoType_::MemType, DataType_> ScalarMatrixB;
      /// scalar matrix type D
      typedef SparseMatrixCOO<typename AlgoType_::MemType, DataType_> ScalarMatrixD;
    };

    // MKL specialisation: There is no ELL implementation, so we choose COO for D matrices
    template<typename DataType_>
    struct MetaMatrixTestHelper<Algo::MKL, DataType_>
    {
      /// scalar vector type
      typedef DenseVector<Mem::Main, DataType_> ScalarVector;
      /// scalar matrix type A
      typedef SparseMatrixCSR<Mem::Main, DataType_> ScalarMatrixA;
      /// scalar matrix type B
      typedef SparseMatrixCSR<Mem::Main, DataType_> ScalarMatrixB;
      /// scalar matrix type D
      typedef SparseMatrixCOO<Mem::Main, DataType_> ScalarMatrixD;
    };

    // CUDA specialisation: There is no COO implementation, so we choose ELL for D matrices
    template<typename DataType_>
    struct MetaMatrixTestHelper<Algo::CUDA, DataType_>
    {
      /// scalar vector type
      typedef DenseVector<Mem::CUDA, DataType_> ScalarVector;
      /// scalar matrix type A
      typedef SparseMatrixCSR<Mem::CUDA, DataType_> ScalarMatrixA;
      /// scalar matrix type B
      typedef SparseMatrixELL<Mem::CUDA, DataType_> ScalarMatrixB;
      /// scalar matrix type D
      typedef SparseMatrixELL<Mem::CUDA, DataType_> ScalarMatrixD;
    };
    /// \endcond

    /**
     * \brief Abstract base class for meta-matrix tests.
     *
     * This class is used by various meta-matrix tests. It offers the functionality of
     * generating a saddle-point system, i.e. a matrix, an rhs and a corresponding solution vector.
     *
     * The matrix, which is generated, has the following structure:
     *
     * \verbatim
       S = / A  B \
           \ D  0 /

       A = / A1  0 \  B = / B1 \
           \ 0  A2 /      \ B2 /

       D = ( D1 D2 )
       \endverbatim
     * where:
     *  - A1, A2 are of type SparseMatrixCSR
     *  - B1, B2 are of type SparseMatrixELL
     *  - D1, D2 are of type SparseMatrixCOO
     *  - A is of type PowerDiagMatrix
     *  - B is of type PowerColMatrix
     *  - D is of type PowerRowMatrix
     *  - S is of type SaddlePointMatrix
     *
     * Putting all this together leads to
     * \verbatim
       SaddlePointMatrix<
         PowerDiagMatrix< SparseMatrixCSR<...> >,
         PowerColMatrix< SparseMatrixELL<...> >,
         PowerRowMatrix< SpaceMatrixCOO<...> >
         >
     * \endverbatim
     *
     * \author Peter Zajac
     */
    template<typename Algo_, typename DataType_>
    class MetaMatrixTestBase
      : public FEAST::TestSystem::TaggedTest<typename Algo_::MemType, DataType_, Algo_>
    {
    public:
      typedef Algo_ AlgoType;
      typedef typename AlgoType::MemType MemType;
      typedef DataType_ DataType;
      typedef MetaMatrixTestHelper<AlgoType, DataType> Helper;

      /// scalar vector type
      typedef typename Helper::ScalarVector ScalarVector;
      /// velocity vector type
      typedef PowerVector<ScalarVector, 2> VeloVector;
      /// pressure vector type
      typedef ScalarVector PresVector;
      /// system vector type
      typedef TupleVector<VeloVector, PresVector> SystemVector;

      /// scalar matrix type A
      typedef typename Helper::ScalarMatrixA ScalarMatrixA;
      /// scalar matrix type B
      typedef typename Helper::ScalarMatrixB ScalarMatrixB;
      /// scalar matrix type D
      typedef typename Helper::ScalarMatrixD ScalarMatrixD;

      /// velocity matrix type
      typedef PowerDiagMatrix<ScalarMatrixA, 2> VeloMatrix;
      /// gradient matrix type
      typedef PowerColMatrix<ScalarMatrixB, 2> GradMatrix;
      /// divergence matrix type
      typedef PowerRowMatrix<ScalarMatrixD, 2> DiveMatrix;

      /// system matrix type
      typedef SaddlePointMatrix<VeloMatrix, GradMatrix, DiveMatrix> SystemMatrix;

      explicit MetaMatrixTestBase(const String & name) :
        FEAST::TestSystem::TaggedTest<typename Algo_::MemType, DataType_, Algo_>(name)
      {
      }

      /// generate test matrix
      static void gen_system(Index m, SystemMatrix& mat_sys, SystemVector& vec_sol, SystemVector& vec_rhs)
      {
        /// create two pointstars
        PointstarFactoryFD<DataType_> ps_fd(m, Index(2));
        PointstarFactoryFE<DataType_> ps_fe(m);

        /// generate the corresponding CSR matrices
        const SparseMatrixCSR<Mem::Main, DataType_> mat_fd(ps_fd.matrix_csr());
        const SparseMatrixCSR<Mem::Main, DataType_> mat_fe(ps_fe.matrix_csr());

        // generate Q2-bubble and eigenvector
        const DenseVector<Mem::Main, DataType_> vec_eigen(ps_fd.eigenvector_min());
        const DenseVector<Mem::Main, DataType_> vec_bubble(ps_fd.vector_q2_bubble());

        // set system matrix
        mat_sys.template at<Index(0),Index(0)>().template at<Index(0),Index(0)>().convert(mat_fe);
        mat_sys.template at<Index(0),Index(0)>().template at<Index(1),Index(1)>().convert(mat_fe);
        mat_sys.template at<Index(0),Index(1)>().template at<Index(0),Index(0)>().convert(mat_fd);
        mat_sys.template at<Index(0),Index(1)>().template at<Index(1),Index(0)>().convert(mat_fd);
        mat_sys.template at<Index(1),Index(0)>().template at<Index(0),Index(0)>().convert(mat_fd);
        mat_sys.template at<Index(1),Index(0)>().template at<Index(0),Index(1)>().convert(mat_fd);

        // set solution vector
        vec_sol.template at<Index(0)>().template at<Index(0)>().convert(vec_bubble);
        vec_sol.template at<Index(0)>().template at<Index(1)>().convert(vec_eigen);
        vec_sol.template at<Index(1)>().convert(vec_eigen);

        // create vectors for rhs computation
        DenseVector<Mem::Main, DataType_> vec_rhs1(vec_bubble.size());
        DenseVector<Mem::Main, DataType_> vec_rhs2(vec_bubble.size());
        DenseVector<Mem::Main, DataType_> vec_rhs3(vec_bubble.size());

        // compute rhs vector (by exploiting the eigenvector property)
        mat_fe.template apply<Algo::Generic>(vec_rhs1, vec_bubble);
        vec_rhs1.template axpy<Algo::Generic>(vec_eigen, vec_rhs1, ps_fd.lambda_min());
        vec_rhs2.template scale<Algo::Generic>(vec_eigen, ps_fd.lambda_min() + ps_fe.lambda_min());
        mat_fd.template apply<Algo::Generic>(vec_rhs3, vec_bubble);
        vec_rhs3.template axpy<Algo::Generic>(vec_eigen, vec_rhs3, ps_fd.lambda_min());

        // set rhs vector
        vec_rhs.template at<Index(0)>().template at<Index(0)>().convert(vec_rhs1);
        vec_rhs.template at<Index(0)>().template at<Index(1)>().convert(vec_rhs2);
        vec_rhs.template at<Index(1)>().convert(vec_rhs3);
      }
    }; // MetaVectorTestBase
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_META_MATRIX_TEST_BASE_HPP
