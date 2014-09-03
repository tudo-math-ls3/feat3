#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/vector_mirror.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for VectorMirror class template
 *
 * \author Peter Zajac
 */
template<
  typename Algo_,
  typename DT_,
  typename IT_>
class VectorMirrorTest
  : public FullTaggedTest<typename Algo_::MemType, Algo_, DT_, IT_>
{
  typedef DenseVector<typename Algo_::MemType, DT_, IT_> VectorType;
  typedef DenseVector<typename Algo_::MemType, IT_, IT_> IVectorType;
  typedef SparseMatrixCSR<typename Algo_::MemType, DT_, IT_> MatrixType;
  typedef VectorMirror<typename Algo_::MemType, DT_, IT_> MirrorType;

public:
  VectorMirrorTest()
    : FullTaggedTest<typename Algo_::MemType, Algo_, DT_, IT_>("VectorMirrorTest")
  {
  }

  virtual void run() const override
  {
    test_1(); // tests gather_prim/scatter_prim
    test_2(); // tests gather_axpy_prim/scatter_axpy_prim
  }

  void test_1() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

    VectorType a1(Index(3), DT_(0));
    VectorType a2(Index(3), DT_(0));
    VectorType b1(Index(3), DT_(0));
    VectorType b2(Index(3), DT_(0));
    VectorType c1(Index(3), DT_(0));
    VectorType c2(Index(3), DT_(0));

    // initialise local vectors
    a1(Index(0), DT_(4));
    a1(Index(1), DT_(7));
    a1(Index(2), DT_(2));
    b1(Index(0), DT_(2));
    b1(Index(1), DT_(1));
    b1(Index(2), DT_(3));
    c1(Index(0), DT_(1));
    c1(Index(1), DT_(5));
    c1(Index(2), DT_(2));

    // initialise global vectors
    a2(Index(0), DT_(7));
    a2(Index(1), DT_(9));
    a2(Index(2), DT_(3));
    b2(Index(0), DT_(7));
    b2(Index(1), DT_(3));
    b2(Index(2), DT_(8));
    c2(Index(0), DT_(7));
    c2(Index(1), DT_(8));
    c2(Index(2), DT_(9));

    // create mirror vectors
    IVectorType col_idx0(Index(1), IT_(0));
    IVectorType col_idx1(Index(1), IT_(1));
    IVectorType col_idx2(Index(1), IT_(2));
    VectorType mir_val(Index(1), DT_(1));
    IVectorType row_ptr(Index(2));
    row_ptr(Index(0), IT_(0));
    row_ptr(Index(1), IT_(1));
    MatrixType mat_gather0(Index(1), Index(3), col_idx0, mir_val, row_ptr);
    MatrixType mat_gather1(Index(1), Index(3), col_idx1, mir_val, row_ptr);
    MatrixType mat_gather2(Index(1), Index(3), col_idx2, mir_val, row_ptr);
    MatrixType mat_scatter0(mat_gather0.transpose());
    MatrixType mat_scatter1(mat_gather1.transpose());
    MatrixType mat_scatter2(mat_gather2.transpose());

    // create mirrors
    MirrorType mirror0(std::move(mat_gather0), std::move(mat_scatter0));
    MirrorType mirror1(std::move(mat_gather1), std::move(mat_scatter1));
    MirrorType mirror2(std::move(mat_gather2), std::move(mat_scatter2));

    // create four buffer and one temporary vectors
    VectorType vec_buf_ab(Index(1), DT_(0));
    VectorType vec_buf_ac(Index(1), DT_(0));
    VectorType vec_buf_bc(Index(1), DT_(0));
    VectorType vec_buf_abc(Index(1), DT_(0));
    VectorType vec_buf_tmp(Index(1), DT_(0));

    // gather from a
    mirror0.template gather_prim<Algo_>(vec_buf_tmp, a1);
    vec_buf_abc.template axpy<Algo_>(vec_buf_tmp, vec_buf_abc);
    mirror1.template gather_prim<Algo_>(vec_buf_tmp, a1);
    vec_buf_ac.template axpy<Algo_>(vec_buf_tmp, vec_buf_ac);
    mirror2.template gather_prim<Algo_>(vec_buf_tmp, a1);
    vec_buf_ab.template axpy<Algo_>(vec_buf_tmp, vec_buf_ab);

    // gather from b
    mirror0.template gather_prim<Algo_>(vec_buf_tmp, b1);
    vec_buf_abc.template axpy<Algo_>(vec_buf_tmp, vec_buf_abc);
    mirror1.template gather_prim<Algo_>(vec_buf_tmp, b1);
    vec_buf_ab.template axpy<Algo_>(vec_buf_tmp, vec_buf_ab);
    mirror2.template gather_prim<Algo_>(vec_buf_tmp, b1);
    vec_buf_bc.template axpy<Algo_>(vec_buf_tmp, vec_buf_bc);

    // gather from c
    mirror0.template gather_prim<Algo_>(vec_buf_tmp, c1);
    vec_buf_abc.template axpy<Algo_>(vec_buf_tmp, vec_buf_abc);
    mirror1.template gather_prim<Algo_>(vec_buf_tmp, c1);
    vec_buf_bc.template axpy<Algo_>(vec_buf_tmp, vec_buf_bc);
    mirror2.template gather_prim<Algo_>(vec_buf_tmp, c1);
    vec_buf_ac.template axpy<Algo_>(vec_buf_tmp, vec_buf_ac);

    // scatter to a
    mirror0.template scatter_prim<Algo_>(a1, vec_buf_abc);
    mirror1.template scatter_prim<Algo_>(a1, vec_buf_ac);
    mirror2.template scatter_prim<Algo_>(a1, vec_buf_ab);

    // scatter to b
    mirror0.template scatter_prim<Algo_>(b1, vec_buf_abc);
    mirror1.template scatter_prim<Algo_>(b1, vec_buf_ab);
    mirror2.template scatter_prim<Algo_>(b1, vec_buf_bc);

    // scatter to c
    mirror0.template scatter_prim<Algo_>(c1, vec_buf_abc);
    mirror1.template scatter_prim<Algo_>(c1, vec_buf_bc);
    mirror2.template scatter_prim<Algo_>(c1, vec_buf_ac);

    // subtract reference results
    a1.template axpy<Algo_>(a2, a1, -DT_(1));
    b1.template axpy<Algo_>(b2, b1, -DT_(1));
    c1.template axpy<Algo_>(c2, c1, -DT_(1));

    // check norms
    TEST_CHECK_EQUAL_WITHIN_EPS(a1.template norm2<Algo_>(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b1.template norm2<Algo_>(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(c1.template norm2<Algo_>(), DT_(0), tol);
  }

  void test_2() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

    VectorType a1(Index(3), DT_(0));
    VectorType a2(Index(3), DT_(0));
    VectorType b1(Index(3), DT_(0));
    VectorType b2(Index(3), DT_(0));
    VectorType c1(Index(3), DT_(0));
    VectorType c2(Index(3), DT_(0));

    // initialise local vectors
    a1(Index(0), DT_(4));
    a1(Index(1), DT_(7));
    a1(Index(2), DT_(2));
    b1(Index(0), DT_(2));
    b1(Index(1), DT_(1));
    b1(Index(2), DT_(3));
    c1(Index(0), DT_(1));
    c1(Index(1), DT_(5));
    c1(Index(2), DT_(2));

    // initialise global vectors
    a2(Index(0), DT_(7));
    a2(Index(1), DT_(9));
    a2(Index(2), DT_(3));
    b2(Index(0), DT_(7));
    b2(Index(1), DT_(3));
    b2(Index(2), DT_(8));
    c2(Index(0), DT_(7));
    c2(Index(1), DT_(8));
    c2(Index(2), DT_(9));

    // create mirror vectors
    IVectorType col_idx0(Index(1), IT_(0));
    IVectorType col_idx1(Index(1), IT_(1));
    IVectorType col_idx2(Index(1), IT_(2));
    VectorType mir_val(Index(1), DT_(1));
    IVectorType row_ptr(Index(2));
    row_ptr(Index(0), IT_(0));
    row_ptr(Index(1), IT_(1));
    MatrixType mat_gather0(Index(1), Index(3), col_idx0, mir_val, row_ptr);
    MatrixType mat_gather1(Index(1), Index(3), col_idx1, mir_val, row_ptr);
    MatrixType mat_gather2(Index(1), Index(3), col_idx2, mir_val, row_ptr);
    MatrixType mat_scatter0(mat_gather0.transpose());
    MatrixType mat_scatter1(mat_gather1.transpose());
    MatrixType mat_scatter2(mat_gather2.transpose());

    // create mirrors
    MirrorType mirror0(std::move(mat_gather0), std::move(mat_scatter0));
    MirrorType mirror1(std::move(mat_gather1), std::move(mat_scatter1));
    MirrorType mirror2(std::move(mat_gather2), std::move(mat_scatter2));

    // create four buffer and one temporary vectors
    VectorType vec_buf_ab(Index(1), DT_(0));
    VectorType vec_buf_ac(Index(1), DT_(0));
    VectorType vec_buf_bc(Index(1), DT_(0));
    VectorType vec_buf_abc(Index(1), DT_(0));
    VectorType vec_buf_tmp(Index(1), DT_(0));

    // gather from a
    mirror0.template gather_axpy_prim<Algo_>(vec_buf_abc, a1);
    mirror1.template gather_axpy_prim<Algo_>(vec_buf_ac, a1);
    mirror2.template gather_axpy_prim<Algo_>(vec_buf_ab, a1);

    // gather from b
    mirror0.template gather_axpy_prim<Algo_>(vec_buf_abc, b1);
    mirror1.template gather_axpy_prim<Algo_>(vec_buf_ab, b1);
    mirror2.template gather_axpy_prim<Algo_>(vec_buf_bc, b1);

    // gather from c
    mirror0.template gather_axpy_prim<Algo_>(vec_buf_abc, c1);
    mirror1.template gather_axpy_prim<Algo_>(vec_buf_bc, c1);
    mirror2.template gather_axpy_prim<Algo_>(vec_buf_ac, c1);

    // scatter to a
    a1.copy(a2);
    mirror0.template scatter_axpy_prim<Algo_>(a1, vec_buf_abc, -DT_(1));
    mirror1.template scatter_axpy_prim<Algo_>(a1, vec_buf_ac, -DT_(1));
    mirror2.template scatter_axpy_prim<Algo_>(a1, vec_buf_ab, -DT_(1));

    // scatter to b
    b1.copy(b2);
    mirror0.template scatter_axpy_prim<Algo_>(b1, vec_buf_abc, -DT_(1));
    mirror1.template scatter_axpy_prim<Algo_>(b1, vec_buf_ab, -DT_(1));
    mirror2.template scatter_axpy_prim<Algo_>(b1, vec_buf_bc, -DT_(1));

    // scatter to c
    c1.copy(c2);
    mirror0.template scatter_axpy_prim<Algo_>(c1, vec_buf_abc, -DT_(1));
    mirror1.template scatter_axpy_prim<Algo_>(c1, vec_buf_bc, -DT_(1));
    mirror2.template scatter_axpy_prim<Algo_>(c1, vec_buf_ac, -DT_(1));

    // check norms
    TEST_CHECK_EQUAL_WITHIN_EPS(a1.template norm2<Algo_>(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b1.template norm2<Algo_>(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(c1.template norm2<Algo_>(), DT_(0), tol);
  }
};

VectorMirrorTest<Algo::Generic, double, Index> vector_mirror_test_generic_d;
