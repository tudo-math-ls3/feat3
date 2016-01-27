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
  typename MemType_,
  typename DT_,
  typename IT_>
class VectorMirrorTest
  : public FullTaggedTest<MemType_, DT_, IT_>
{
  typedef DenseVector<MemType_, DT_, IT_> VectorType;
  typedef DenseVector<MemType_, IT_, IT_> IVectorType;
  typedef SparseMatrixCSR<MemType_, DT_, IT_> MatrixType;
  typedef VectorMirror<MemType_, DT_, IT_> MirrorType;

public:
  VectorMirrorTest()
    : FullTaggedTest<MemType_, DT_, IT_>("VectorMirrorTest")
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
    auto vec_buf_abc(mirror0.create_buffer_vector());
    VectorType vec_buf_tmp(Index(1), DT_(0));

    // gather from a
    mirror0.gather_prim(vec_buf_tmp, a1);
    vec_buf_abc.axpy(vec_buf_tmp, vec_buf_abc);
    mirror1.gather_prim(vec_buf_tmp, a1);
    vec_buf_ac.axpy(vec_buf_tmp, vec_buf_ac);
    mirror2.gather_prim(vec_buf_tmp, a1);
    vec_buf_ab.axpy(vec_buf_tmp, vec_buf_ab);

    // gather from b
    mirror0.gather_prim(vec_buf_tmp, b1);
    vec_buf_abc.axpy(vec_buf_tmp, vec_buf_abc);
    mirror1.gather_prim(vec_buf_tmp, b1);
    vec_buf_ab.axpy(vec_buf_tmp, vec_buf_ab);
    mirror2.gather_prim(vec_buf_tmp, b1);
    vec_buf_bc.axpy(vec_buf_tmp, vec_buf_bc);

    // gather from c
    mirror0.gather_prim(vec_buf_tmp, c1);
    vec_buf_abc.axpy(vec_buf_tmp, vec_buf_abc);
    mirror1.gather_prim(vec_buf_tmp, c1);
    vec_buf_bc.axpy(vec_buf_tmp, vec_buf_bc);
    mirror2.gather_prim(vec_buf_tmp, c1);
    vec_buf_ac.axpy(vec_buf_tmp, vec_buf_ac);

    // scatter to a
    mirror0.scatter_prim(a1, vec_buf_abc);
    mirror1.scatter_prim(a1, vec_buf_ac);
    mirror2.scatter_prim(a1, vec_buf_ab);

    // scatter to b
    mirror0.scatter_prim(b1, vec_buf_abc);
    mirror1.scatter_prim(b1, vec_buf_ab);
    mirror2.scatter_prim(b1, vec_buf_bc);

    // scatter to c
    mirror0.scatter_prim(c1, vec_buf_abc);
    mirror1.scatter_prim(c1, vec_buf_bc);
    mirror2.scatter_prim(c1, vec_buf_ac);

    // subtract reference results
    a1.axpy(a2, a1, -DT_(1));
    b1.axpy(b2, b1, -DT_(1));
    c1.axpy(c2, c1, -DT_(1));

    // check norms
    TEST_CHECK_EQUAL_WITHIN_EPS(a1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(c1.norm2(), DT_(0), tol);
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
    auto vec_buf_ab(mirror0.create_buffer_vector());
    VectorType vec_buf_ac(Index(1), DT_(0));
    VectorType vec_buf_bc(Index(1), DT_(0));
    VectorType vec_buf_abc(Index(1), DT_(0));
    VectorType vec_buf_tmp(Index(1), DT_(0));

    // gather from a
    mirror0.gather_axpy_prim(vec_buf_abc, a1);
    mirror1.gather_axpy_prim(vec_buf_ac, a1);
    mirror2.gather_axpy_prim(vec_buf_ab, a1);

    // gather from b
    mirror0.gather_axpy_prim(vec_buf_abc, b1);
    mirror1.gather_axpy_prim(vec_buf_ab, b1);
    mirror2.gather_axpy_prim(vec_buf_bc, b1);

    // gather from c
    mirror0.gather_axpy_prim(vec_buf_abc, c1);
    mirror1.gather_axpy_prim(vec_buf_bc, c1);
    mirror2.gather_axpy_prim(vec_buf_ac, c1);

    // scatter to a
    a1.copy(a2);
    mirror0.scatter_axpy_prim(a1, vec_buf_abc, -DT_(1));
    mirror1.scatter_axpy_prim(a1, vec_buf_ac, -DT_(1));
    mirror2.scatter_axpy_prim(a1, vec_buf_ab, -DT_(1));

    // scatter to b
    b1.copy(b2);
    mirror0.scatter_axpy_prim(b1, vec_buf_abc, -DT_(1));
    mirror1.scatter_axpy_prim(b1, vec_buf_ab, -DT_(1));
    mirror2.scatter_axpy_prim(b1, vec_buf_bc, -DT_(1));

    // scatter to c
    c1.copy(c2);
    mirror0.scatter_axpy_prim(c1, vec_buf_abc, -DT_(1));
    mirror1.scatter_axpy_prim(c1, vec_buf_bc, -DT_(1));
    mirror2.scatter_axpy_prim(c1, vec_buf_ac, -DT_(1));

    // check norms
    TEST_CHECK_EQUAL_WITHIN_EPS(a1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(c1.norm2(), DT_(0), tol);
  }
};

VectorMirrorTest<Mem::Main, double, Index> vector_mirror_test_generic_d;

/**
 * \brief Test class for VectorMirrorBlocked class template
 *
 * \author Jordi Paul
 */
template < typename Mem_, typename DT_, typename IT_, int BS_>
class VectorMirrorBlockedTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
  typedef DenseVectorBlocked<Mem_, DT_, IT_, BS_> VectorType;
  typedef DenseVector<Mem_, IT_, IT_> IVectorType;
  typedef SparseMatrixCSR<Mem_, DT_, IT_> MatrixType;
  typedef DenseVector<Mem_, DT_, IT_> BufferVectorType;
  typedef VectorMirrorBlocked<Mem_, DT_, IT_, BS_> MirrorType;

  static constexpr int BlockSize = BS_;

public:
  VectorMirrorBlockedTest()
    : FullTaggedTest<Mem_, DT_, IT_>("VectorMirrorBlockedTest")
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

    Tiny::Vector<DT_, BS_> tiny_tmp(DT_(0));

    // Vectors to gather/scatter
    tiny_tmp(0) = DT_(-3), tiny_tmp(1) = DT_(0.5), tiny_tmp(BlockSize-1) = -DT_(1);
    a1(Index(0), tiny_tmp);
    tiny_tmp(0) = DT_(-7), tiny_tmp(1) = DT_(5), tiny_tmp(BlockSize-1) = -DT_(1000);
    a1(Index(1), tiny_tmp);
    tiny_tmp(0) = DT_(-0.3), tiny_tmp(1) = DT_(0.1), tiny_tmp(BlockSize-1) = DT_(2);
    a1(Index(2), tiny_tmp);

    tiny_tmp(0) = DT_(0.3), tiny_tmp(1) = DT_(0.001), tiny_tmp(BlockSize-1) = DT_(2.2);
    b1(Index(0), tiny_tmp);
    tiny_tmp(0) = DT_(-3.7), tiny_tmp(1) = DT_(8.001), tiny_tmp(BlockSize-1) = -DT_(9.002);
    b1(Index(1), tiny_tmp);
    tiny_tmp(0) = DT_(19), tiny_tmp(1) = DT_(-111), tiny_tmp(BlockSize-1) = -DT_(111);
    b1(Index(2), tiny_tmp);

    tiny_tmp(0) = DT_(4711), tiny_tmp(1) = DT_(0), tiny_tmp(BlockSize-1) = DT_(0);
    c1(Index(0), tiny_tmp);
    tiny_tmp(0) = DT_(0), tiny_tmp(1) = DT_(1), tiny_tmp(BlockSize-1) = DT_(7);
    c1(Index(1), tiny_tmp);
    tiny_tmp(0) = DT_(0.001), tiny_tmp(1) = DT_(1), tiny_tmp(BlockSize-1) = -DT_(7.7);
    c1(Index(2), tiny_tmp);

    // Supposed results. See the gather operations below for this to make sense.
    a2(Index(0), a1(0)+b1(0)+c1(0));
    a2(Index(1), a1(1)+c1(2));
    a2(Index(2), a1(2)+b1(1));

    b2(Index(0), a1(0)+b1(0)+c1(0));
    b2(Index(1), a1(2)+b1(1));
    b2(Index(2), b1(2)+c1(1));

    c2(Index(0), a1(0)+b1(0)+c1(0));
    c2(Index(1), b1(2)+c1(1));
    c2(Index(2), a1(1)+c1(2));

    // Create mirror vectors. Every mirror gathers one entry from / scatters one entry to a DenseVectorBlocked of
    // length 3
    IVectorType col_idx0(Index(1), IT_(0));
    IVectorType col_idx1(Index(1), IT_(1));
    IVectorType col_idx2(Index(1), IT_(2));
    BufferVectorType mir_val(Index(1), DT_(1));
    IVectorType row_ptr(Index(2));
    row_ptr(Index(0), IT_(0));
    row_ptr(Index(1), IT_(1));
    MatrixType mat_gather0(Index(1), Index(3), col_idx0, mir_val, row_ptr);
    MatrixType mat_gather1(Index(1), Index(3), col_idx1, mir_val, row_ptr);
    MatrixType mat_gather2(Index(1), Index(3), col_idx2, mir_val, row_ptr);
    MatrixType mat_scatter0(mat_gather0.transpose());
    MatrixType mat_scatter1(mat_gather1.transpose());
    MatrixType mat_scatter2(mat_gather2.transpose());

    // Create mirrors. mirror{k} gathers/scatters to/from entry k
    MirrorType mirror0(std::move(mat_gather0), std::move(mat_scatter0));
    MirrorType mirror1(std::move(mat_gather1), std::move(mat_scatter1));
    MirrorType mirror2(std::move(mat_gather2), std::move(mat_scatter2));

    // Create four buffer and one temporary vectors
    // a1(2)+b1(1)
    BufferVectorType vec_buf_ab(Index(BlockSize)*Index(1), DT_(0));
    // a1(1)+c1(2)
    BufferVectorType vec_buf_ac(Index(BlockSize)*Index(1), DT_(0));
    // b1(2)+c(1)
    BufferVectorType vec_buf_bc(Index(BlockSize)*Index(1), DT_(0));
    // a1(0) + b1(0) + c1(0)
    auto vec_buf_abc(mirror0.create_buffer_vector());
    // Generic buffer for adding
    BufferVectorType vec_buf_tmp(Index(BlockSize)*Index(1), DT_(0));

    // Gather from a
    mirror0.gather_prim(vec_buf_tmp, a1);
    vec_buf_abc.axpy(vec_buf_tmp, vec_buf_abc);
    mirror1.gather_prim(vec_buf_tmp, a1);
    vec_buf_ac.axpy(vec_buf_tmp, vec_buf_ac);
    mirror2.gather_prim(vec_buf_tmp, a1);
    vec_buf_ab.axpy(vec_buf_tmp, vec_buf_ab);

    // Gather from b
    mirror0.gather_prim(vec_buf_tmp, b1);
    vec_buf_abc.axpy(vec_buf_tmp, vec_buf_abc);
    mirror1.gather_prim(vec_buf_tmp, b1);
    vec_buf_ab.axpy(vec_buf_tmp, vec_buf_ab);
    mirror2.gather_prim(vec_buf_tmp, b1);
    vec_buf_bc.axpy(vec_buf_tmp, vec_buf_bc);

    // Gather from c
    mirror0.gather_prim(vec_buf_tmp, c1);
    vec_buf_abc.axpy(vec_buf_tmp, vec_buf_abc);
    mirror1.gather_prim(vec_buf_tmp, c1);
    vec_buf_bc.axpy(vec_buf_tmp, vec_buf_bc);
    mirror2.gather_prim(vec_buf_tmp, c1);
    vec_buf_ac.axpy(vec_buf_tmp, vec_buf_ac);

    // Now we no longer need a1, b1, c1 and overwrite them by scattering the buffers back into them
    mirror0.scatter_prim(a1, vec_buf_abc);
    mirror1.scatter_prim(a1, vec_buf_ac);
    mirror2.scatter_prim(a1, vec_buf_ab);

    mirror0.scatter_prim(b1, vec_buf_abc);
    mirror1.scatter_prim(b1, vec_buf_ab);
    mirror2.scatter_prim(b1, vec_buf_bc);

    mirror0.scatter_prim(c1, vec_buf_abc);
    mirror1.scatter_prim(c1, vec_buf_bc);
    mirror2.scatter_prim(c1, vec_buf_ac);

    // Subtract reference results
    a1.axpy(a2, a1, -DT_(1));
    b1.axpy(b2, b1, -DT_(1));
    c1.axpy(c2, c1, -DT_(1));

    // Check norms
    TEST_CHECK_EQUAL_WITHIN_EPS(a1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(c1.norm2(), DT_(0), tol);
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

    Tiny::Vector<DT_, BS_> tiny_tmp(DT_(0));

    // Vectors to gather/scatter
    tiny_tmp(0) = DT_(-0.3), tiny_tmp(1) = DT_(0.1), tiny_tmp(BlockSize-1) = DT_(1.999);
    a1(Index(0), tiny_tmp);
    tiny_tmp(0) = DT_(3.01), tiny_tmp(1) = -DT_(0.5), tiny_tmp(BlockSize-1) = -DT_(219);
    a1(Index(1), tiny_tmp);
    tiny_tmp(0) = DT_(7), tiny_tmp(1) = -DT_(5), tiny_tmp(BlockSize-1) = -DT_(1000);
    a1(Index(2), tiny_tmp);

    tiny_tmp(0) = DT_(0.3), tiny_tmp(1) = DT_(1.00), tiny_tmp(BlockSize-1) = -DT_(2.9);
    b1(Index(0), tiny_tmp);
    tiny_tmp(0) = -DT_(9), tiny_tmp(1) = DT_(11), tiny_tmp(BlockSize-1) = -DT_(1011);
    b1(Index(1), tiny_tmp);
    tiny_tmp(0) = DT_(-302.7), tiny_tmp(1) = DT_(8.1), tiny_tmp(BlockSize-1) = DT_(2.009);
    b1(Index(2), tiny_tmp);

    tiny_tmp(0) = -DT_(0.1), tiny_tmp(1) = DT_(1), tiny_tmp(BlockSize-1) = -DT_(8.7);
    c1(Index(0), tiny_tmp);
    tiny_tmp(0) = DT_(0), tiny_tmp(1) = DT_(1), tiny_tmp(BlockSize-1) = DT_(79);
    c1(Index(1), tiny_tmp);
    tiny_tmp(0) = DT_(4711), tiny_tmp(1) = DT_(9.6), tiny_tmp(BlockSize-1) = DT_(0);
    c1(Index(2), tiny_tmp);

    // Supposed results. See the gather operations below for this to make sense.
    a2(Index(0), a1(0)+b1(0)+c1(0));
    a2(Index(1), a1(1)+c1(2));
    a2(Index(2), a1(2)+b1(1));

    b2(Index(0), a1(0)+b1(0)+c1(0));
    b2(Index(1), a1(2)+b1(1));
    b2(Index(2), b1(2)+c1(1));

    c2(Index(0), a1(0)+b1(0)+c1(0));
    c2(Index(1), b1(2)+c1(1));
    c2(Index(2), a1(1)+c1(2));

    // Create mirror vectors. Every mirror gathers one entry from / scatters one entry to a DenseVectorBlocked of
    // length 3
    IVectorType col_idx0(Index(1), IT_(0));
    IVectorType col_idx1(Index(1), IT_(1));
    IVectorType col_idx2(Index(1), IT_(2));
    BufferVectorType mir_val(Index(1), DT_(1));
    IVectorType row_ptr(Index(2));
    row_ptr(Index(0), IT_(0));
    row_ptr(Index(1), IT_(1));
    MatrixType mat_gather0(Index(1), Index(3), col_idx0, mir_val, row_ptr);
    MatrixType mat_gather1(Index(1), Index(3), col_idx1, mir_val, row_ptr);
    MatrixType mat_gather2(Index(1), Index(3), col_idx2, mir_val, row_ptr);
    MatrixType mat_scatter0(mat_gather0.transpose());
    MatrixType mat_scatter1(mat_gather1.transpose());
    MatrixType mat_scatter2(mat_gather2.transpose());

    // Create mirrors. mirror{k} gathers/scatters to/from entry k
    MirrorType mirror0(std::move(mat_gather0), std::move(mat_scatter0));
    MirrorType mirror1(std::move(mat_gather1), std::move(mat_scatter1));
    MirrorType mirror2(std::move(mat_gather2), std::move(mat_scatter2));

    // Create four buffer and one temporary vectors
    // a1(2)+b1(1)
    BufferVectorType vec_buf_ab(Index(BlockSize)*Index(1), DT_(0));
    // a1(1)+c1(2)
    auto vec_buf_ac(mirror1.create_buffer_vector());
    // b1(2)+c(1)
    BufferVectorType vec_buf_bc(Index(BlockSize)*Index(1), DT_(0));
    // a1(0) + b1(0) + c1(0)
    BufferVectorType vec_buf_abc(Index(BlockSize)*Index(1), DT_(0));
    BufferVectorType vec_buf_tmp(Index(BlockSize)*Index(1), DT_(0));

    // Gather from a
    mirror0.gather_axpy_prim(vec_buf_abc, a1);
    mirror1.gather_axpy_prim(vec_buf_ac, a1);
    mirror2.gather_axpy_prim(vec_buf_ab, a1);

    // Gather from b
    mirror0.gather_axpy_prim(vec_buf_abc, b1);
    mirror1.gather_axpy_prim(vec_buf_ab, b1);
    mirror2.gather_axpy_prim(vec_buf_bc, b1);

    // Gather from c
    mirror0.gather_axpy_prim(vec_buf_abc, c1);
    mirror1.gather_axpy_prim(vec_buf_bc, c1);
    mirror2.gather_axpy_prim(vec_buf_ac, c1);

    // Now we no longer need a1, b1, c1 and overwrite them by scattering the buffers back into them
    a1.copy(a2);
    mirror0.scatter_axpy_prim(a1, vec_buf_abc, -DT_(1));
    mirror1.scatter_axpy_prim(a1, vec_buf_ac, -DT_(1));
    mirror2.scatter_axpy_prim(a1, vec_buf_ab, -DT_(1));

    b1.copy(b2);
    mirror0.scatter_axpy_prim(b1, vec_buf_abc, -DT_(1));
    mirror1.scatter_axpy_prim(b1, vec_buf_ab, -DT_(1));
    mirror2.scatter_axpy_prim(b1, vec_buf_bc, -DT_(1));

    c1.copy(c2);
    mirror0.scatter_axpy_prim(c1, vec_buf_abc, -DT_(1));
    mirror1.scatter_axpy_prim(c1, vec_buf_bc, -DT_(1));
    mirror2.scatter_axpy_prim(c1, vec_buf_ac, -DT_(1));

    // Check norms
    TEST_CHECK_EQUAL_WITHIN_EPS(a1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(b1.norm2(), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(c1.norm2(), DT_(0), tol);
  }
};

VectorMirrorBlockedTest<Mem::Main, float, Index, 2> vector_mirror_blocked_test_generic_f;
VectorMirrorBlockedTest<Mem::Main, double, unsigned int, 3> vector_mirror_blocked_test_generic_d;
