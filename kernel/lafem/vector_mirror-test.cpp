// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_vector.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for VectorMirror class template
 *
 * \author Peter Zajac
 */
template<
  typename DT_,
  typename IT_>
class VectorMirrorDenseVectorTest
  : public UnitTest
{
  typedef DenseVector<DT_, IT_> BufferVectorType;
  typedef DenseVector<DT_, IT_> VectorType;
  typedef VectorMirror<DT_, IT_> MirrorType;

public:
  explicit VectorMirrorDenseVectorTest(PreferredBackend backend)
    : UnitTest("VectorMirrorDenseVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    // create mirrors
    MirrorType mirror0(Index(3), Index(1));
    MirrorType mirror1(Index(3), Index(1));
    MirrorType mirror2(Index(3), Index(1));
    mirror0.indices_view_w()[0] = Index(0);
    mirror1.indices_view_w()[0] = Index(1);
    mirror2.indices_view_w()[0] = Index(2);

    // Test for DenseVector
    {
      BufferVectorType vec_buf_ab_a(Index(1), DT_(0));
      BufferVectorType vec_buf_ab_b(Index(1), DT_(0));
      BufferVectorType vec_buf_ac_a(Index(1), DT_(0));
      BufferVectorType vec_buf_ac_c(Index(1), DT_(0));
      BufferVectorType vec_buf_bc_b(Index(1), DT_(0));
      BufferVectorType vec_buf_bc_c(Index(1), DT_(0));
      BufferVectorType vec_buf_abc_a(Index(1), DT_(0));
      BufferVectorType vec_buf_abc_b(Index(1), DT_(0));
      BufferVectorType vec_buf_abc_c(Index(1), DT_(0));

      VectorType a1(Index(3), DT_(0));
      VectorType a2(Index(3), DT_(0));
      VectorType b1(Index(3), DT_(0));
      VectorType b2(Index(3), DT_(0));
      VectorType c1(Index(3), DT_(0));
      VectorType c2(Index(3), DT_(0));

      // initialize local vectors
      {
        Memory::TypedView<DT_> va1 = a1.elements_view_w();
        Memory::TypedView<DT_> vb1 = b1.elements_view_w();
        Memory::TypedView<DT_> vc1 = c1.elements_view_w();
        va1[Index(0)] = DT_(4);
        va1[Index(1)] = DT_(7);
        va1[Index(2)] = DT_(2);
        vb1[Index(0)] = DT_(2);
        vb1[Index(1)] = DT_(1);
        vb1[Index(2)] = DT_(3);
        vc1[Index(0)] = DT_(1);
        vc1[Index(1)] = DT_(5);
        vc1[Index(2)] = DT_(2);
      }

      // initialize global vectors
      {
        Memory::TypedView<DT_> va2 = a2.elements_view_w();
        Memory::TypedView<DT_> vb2 = b2.elements_view_w();
        Memory::TypedView<DT_> vc2 = c2.elements_view_w();
        va2[Index(0)] = DT_(7);
        va2[Index(1)] = DT_(9);
        va2[Index(2)] = DT_(3);
        vb2[Index(0)] = DT_(7);
        vb2[Index(1)] = DT_(3);
        vb2[Index(2)] = DT_(8);
        vc2[Index(0)] = DT_(7);
        vc2[Index(1)] = DT_(8);
        vc2[Index(2)] = DT_(9);
      }

      // gather from a
      mirror0.gather(vec_buf_abc_a, a1);
      mirror1.gather(vec_buf_ac_a, a1);
      mirror2.gather(vec_buf_ab_a, a1);

      // gather from b
      mirror0.gather(vec_buf_abc_b, b1);
      mirror1.gather(vec_buf_ab_b, b1);
      mirror2.gather(vec_buf_bc_b, b1);

      // gather from c
      mirror0.gather(vec_buf_abc_c, c1);
      mirror1.gather(vec_buf_bc_c, c1);
      mirror2.gather(vec_buf_ac_c, c1);

      // scatter to a
      mirror0.scatter_axpy(a1, vec_buf_abc_b);
      mirror0.scatter_axpy(a1, vec_buf_abc_c);
      mirror1.scatter_axpy(a1, vec_buf_ac_c);
      mirror2.scatter_axpy(a1, vec_buf_ab_b);

      // scatter to b
      mirror0.scatter_axpy(b1, vec_buf_abc_a);
      mirror0.scatter_axpy(b1, vec_buf_abc_c);
      mirror1.scatter_axpy(b1, vec_buf_ab_a);
      mirror2.scatter_axpy(b1, vec_buf_bc_c);

      // scatter to c
      mirror0.scatter_axpy(c1, vec_buf_abc_a);
      mirror0.scatter_axpy(c1, vec_buf_abc_b);
      mirror1.scatter_axpy(c1, vec_buf_bc_b);
      mirror2.scatter_axpy(c1, vec_buf_ac_a);

      // check against reference results
      TEST_CHECK_LESS_THAN(a1.max_rel_diff(a2), tol);
      TEST_CHECK_LESS_THAN(b1.max_rel_diff(b2), tol);
      TEST_CHECK_LESS_THAN(c1.max_rel_diff(c2), tol);
    } // End of test for DenseVector
  }
};

SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, float, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, double, std::uint64_t, PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif

/**
 * \brief Test class for VectorMirror class template
 *
 * \author Peter Zajac
 */
template<
  typename DT_,
  typename IT_>
class VectorMirrorSparseVectorTest
  : public UnitTest
{
  typedef DenseVector<DT_, IT_> BufferVectorType;
  typedef SparseVector<DT_, IT_> SparseVectorType;
  typedef VectorMirror<DT_, IT_> MirrorType;

public:
  explicit VectorMirrorSparseVectorTest(PreferredBackend backend)
    : UnitTest("VectorMirrorSparseVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    // create mirrors
    MirrorType mirror0(Index(3), Index(1));
    MirrorType mirror1(Index(3), Index(1));
    MirrorType mirror2(Index(3), Index(1));
    mirror0.indices_view_w()[0] = Index(0);
    mirror1.indices_view_w()[0] = Index(1);
    mirror2.indices_view_w()[0] = Index(2);

    // Test for SparseVectors
    {
      BufferVectorType vec_buf_ab_a(Index(1), DT_(0));
      BufferVectorType vec_buf_ab_b(Index(1), DT_(0));
      BufferVectorType vec_buf_ac_a(Index(1), DT_(0));
      BufferVectorType vec_buf_ac_c(Index(1), DT_(0));
      BufferVectorType vec_buf_bc_b(Index(1), DT_(0));
      BufferVectorType vec_buf_bc_c(Index(1), DT_(0));
      BufferVectorType vec_buf_abc_a(Index(1), DT_(0));
      BufferVectorType vec_buf_abc_b(Index(1), DT_(0));
      BufferVectorType vec_buf_abc_c(Index(1), DT_(0));

      SparseVectorType a1(Index(3), Index(2));
      SparseVectorType a2(Index(3), Index(2));
      SparseVectorType b1(Index(3), Index(2));
      SparseVectorType b2(Index(3), Index(2));
      SparseVectorType c1(Index(3), Index(2));
      SparseVectorType c2(Index(3), Index(2));

      // Vectors to gather/scatter. Note that we use the emplacement operator which also sets the sparsity pattern.
      {
        Memory::TypedView<IT_> a1i = a1.indices_view_w();
        Memory::TypedView<DT_> a1v = a1.elements_view_w();
        a1i[0] = IT_(0);
        a1i[1] = IT_(1);
        a1v[0] = DT_(4);
        a1v[1] = DT_(7);
      }

      {
        Memory::TypedView<IT_> b1i = b1.indices_view_w();
        Memory::TypedView<DT_> b1v = b1.elements_view_w();
        b1i[0] = IT_(1);
        b1i[1] = IT_(2);
        b1v[0] = DT_(1);
        b1v[1] = DT_(3);
      }

      {
        Memory::TypedView<IT_> c1i = c1.indices_view_w();
        Memory::TypedView<DT_> c1v = c1.elements_view_w();
        c1i[0] = IT_(0);
        c1i[1] = IT_(2);
        c1v[0] = DT_(1);
        c1v[1] = DT_(2);
      }

      // Supposed results. See the gather operations below for this to make sense.
      // Note that each sparse vector is missing one entry and this defines the "sparsity" pattern
      {
        Memory::TypedView<IT_> a2i = a2.indices_view_w();
        Memory::TypedView<DT_> a2v = a2.elements_view_w();
        a2i[0] = IT_(0);
        a2i[1] = IT_(1);
        a2v[0] = DT_(5);
        a2v[1] = DT_(9);
      }

      {
        Memory::TypedView<IT_> b2i = b2.indices_view_w();
        Memory::TypedView<DT_> b2v = b2.elements_view_w();
        b2i[0] = IT_(1);
        b2i[1] = IT_(2);
        b2v[0] = DT_(1);
        b2v[1] = DT_(3);
      }

      {
        Memory::TypedView<IT_> c2i = c2.indices_view_w();
        Memory::TypedView<DT_> c2v = c2.elements_view_w();
        c2i[0] = IT_(0);
        c2i[1] = IT_(2);
        c2v[0] = DT_(5);
        c2v[1] = DT_(9);
      }

      // gather from a
      mirror0.gather(vec_buf_abc_a, a1);
      mirror1.gather(vec_buf_ac_a, a1);
      mirror2.gather(vec_buf_ab_a, a1);

      // gather from b
      mirror0.gather(vec_buf_abc_b, b1);
      mirror1.gather(vec_buf_ab_b, b1);
      mirror2.gather(vec_buf_bc_b, b1);

      // gather from c
      mirror0.gather(vec_buf_abc_c, c1);
      mirror1.gather(vec_buf_bc_c, c1);
      mirror2.gather(vec_buf_ac_c, c1);

      // scatter to a
      mirror0.scatter_axpy(a1, vec_buf_abc_b);
      mirror0.scatter_axpy(a1, vec_buf_abc_c);
      mirror1.scatter_axpy(a1, vec_buf_ac_c);
      mirror2.scatter_axpy(a1, vec_buf_ab_b);

      // scatter to b
      mirror0.scatter_axpy(b1, vec_buf_abc_a);
      mirror0.scatter_axpy(b1, vec_buf_abc_c);
      mirror1.scatter_axpy(b1, vec_buf_ab_a);
      mirror2.scatter_axpy(b1, vec_buf_bc_c);

      // scatter to c
      mirror0.scatter_axpy(c1, vec_buf_abc_a);
      mirror0.scatter_axpy(c1, vec_buf_abc_b);
      mirror1.scatter_axpy(c1, vec_buf_bc_b);
      mirror2.scatter_axpy(c1, vec_buf_ac_a);

      TEST_CHECK_EQUAL(a1.num_nzes(), a2.num_nzes());
      TEST_CHECK(a1.same_layout(a2));
      TEST_CHECK_LESS_THAN(a1.max_rel_diff(a2), tol);

      TEST_CHECK_EQUAL(b1.num_nzes(), b2.num_nzes());
      TEST_CHECK(b1.same_layout(b2));
      TEST_CHECK_LESS_THAN(b1.max_rel_diff(b2), tol);

      TEST_CHECK_EQUAL(c1.num_nzes(), c2.num_nzes());
      TEST_CHECK(c1.same_layout(c2));
      TEST_CHECK_LESS_THAN(c1.max_rel_diff(c2), tol);
    } // End of test for SparseVector
  }
};

SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, float, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, double, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, double, std::uint64_t, PreferredBackend::cuda);
//#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, Half, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorTest, Half, std::uint64_t, PreferredBackend::cuda);
//#endif
//#endif

/**
 * \brief Test class for VectorMirror class template
 *
 * \author Peter Zajac
 */
template<
  typename DT_,
  typename IT_>
class VectorMirrorDenseVectorBlockedTest
  : public UnitTest
{
  static constexpr int block_size = 3;

  typedef DenseVector<DT_, IT_> BufferVectorType;
  typedef DenseVectorBlocked<DT_, IT_, block_size> BlockedVectorType;
  typedef VectorMirror<DT_, IT_> MirrorType;

public:
  explicit VectorMirrorDenseVectorBlockedTest(PreferredBackend backend)
    : UnitTest("VectorMirrorDenseVectorBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    // create mirrors
    MirrorType mirror0(Index(3), Index(1));
    MirrorType mirror1(Index(3), Index(1));
    MirrorType mirror2(Index(3), Index(1));
    mirror0.indices_view_w()[0] = Index(0);
    mirror1.indices_view_w()[0] = Index(1);
    mirror2.indices_view_w()[0] = Index(2);

    // Tests for DenseVectorBlocked
    {
      // a1(2)+b1(1)
      BufferVectorType vec_buf_ab_a(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_ab_b(Index(block_size)*Index(1), DT_(0));
      // a1(1)+c1(2)
      BufferVectorType vec_buf_ac_a(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_ac_c(Index(block_size)*Index(1), DT_(0));
      // b1(2)+c(1)
      BufferVectorType vec_buf_bc_b(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_bc_c(Index(block_size)*Index(1), DT_(0));
      // a1(0) + b1(0) + c1(0)
      BufferVectorType vec_buf_abc_a(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_abc_b(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_abc_c(Index(block_size)*Index(1), DT_(0));

      BlockedVectorType a1(Index(3), DT_(0));
      BlockedVectorType a2(Index(3), DT_(0));
      BlockedVectorType b1(Index(3), DT_(0));
      BlockedVectorType b2(Index(3), DT_(0));
      BlockedVectorType c1(Index(3), DT_(0));
      BlockedVectorType c2(Index(3), DT_(0));

      typedef typename BlockedVectorType::ValueType ValueType;

      // Vectors to gather/scatter
      {
        Memory::TypedView<ValueType> va1 = a1.elements_view_rw();
        va1[0][0] = DT_(-3);
        va1[0][1] = DT_(0.5);
        va1[0][2] = -DT_(1);
        va1[1][0] = DT_(-7);
        va1[1][1] = DT_(5);
        va1[1][2] = -DT_(1000);
        va1[2][0] = DT_(-0.3);
        va1[2][1] = DT_(0.1);
        va1[2][2] = DT_(2);
      }

      {
        Memory::TypedView<ValueType> vb1 = b1.elements_view_rw();
        vb1[0][0] = DT_(0.3);
        vb1[0][1] = DT_(0.001);
        vb1[0][2] = DT_(2.2);
        vb1[1][0] = DT_(-3.7);
        vb1[1][1] = DT_(8.001);
        vb1[1][2] = -DT_(9.002);
        vb1[2][0] = DT_(19);
        vb1[2][1] = DT_(-111);
        vb1[2][2] = -DT_(111);
      }

      {
        Memory::TypedView<ValueType> vc1 = c1.elements_view_rw();
        vc1[0][0] = DT_(2303);
        vc1[0][1] = DT_(0);
        vc1[0][2] = DT_(0);
        vc1[1][0] = DT_(0);
        vc1[1][1] = DT_(1);
        vc1[1][2] = DT_(7);
        vc1[2][0] = DT_(0.001);
        vc1[2][1] = DT_(1);
        vc1[2][2] = -DT_(7.7);
      }

      // Supposed results. See the gather operations below for this to make sense.
      {
        Memory::TypedView<ValueType> va2 = a2.elements_view_rw();
        va2[0][0] = DT_(-3)  + DT_(0.3)   + DT_(2303);
        va2[0][1] = DT_(0.5) + DT_(0.001) + DT_(0);
        va2[0][2] = -DT_(1)  + DT_(2.2)   + DT_(0);
        va2[1][0] = DT_(-7)    + DT_(0.001);
        va2[1][1] = DT_(5)     + DT_(1);
        va2[1][2] = -DT_(1000) - DT_(7.7);
        va2[2][0] = DT_(-0.3) + DT_(-3.7);
        va2[2][1] = DT_(0.1)  + DT_(8.001);
        va2[2][2] = DT_(2)    - DT_(9.002);
      }

      {
        Memory::TypedView<ValueType> vb2 = b2.elements_view_rw();
        vb2[0][0] = DT_(-3)  + DT_(0.3)   + DT_(2303);
        vb2[0][1] = DT_(0.5) + DT_(0.001) + DT_(0);
        vb2[0][2] = -DT_(1)  + DT_(2.2)   + DT_(0);
        vb2[1][0] = DT_(-0.3) + DT_(-3.7);
        vb2[1][1] = DT_(0.1)  + DT_(8.001);
        vb2[1][2] = DT_(2)    - DT_(9.002);
        vb2[2][0] = DT_(19)   + DT_(0);
        vb2[2][1] = DT_(-111) + DT_(1);
        vb2[2][2] = -DT_(111) + DT_(7);
      }

      {
        Memory::TypedView<ValueType> vc2 = c2.elements_view_rw();
        vc2[0][0] = DT_(-3)  + DT_(0.3)   + DT_(2303);
        vc2[0][1] = DT_(0.5) + DT_(0.001) + DT_(0);
        vc2[0][2] = -DT_(1)  + DT_(2.2)   + DT_(0);
        vc2[1][0] = DT_(0) + DT_(19);
        vc2[1][1] = DT_(1) + DT_(-111);
        vc2[1][2] = DT_(7) - DT_(111);
        vc2[2][0] = DT_(-7)    + DT_(0.001);
        vc2[2][1] = DT_(5)     + DT_(1);
        vc2[2][2] = -DT_(1000) - DT_(7.7);
      }

      // gather from a
      mirror0.gather(vec_buf_abc_a, a1);
      mirror1.gather(vec_buf_ac_a, a1);
      mirror2.gather(vec_buf_ab_a, a1);

      // gather from b
      mirror0.gather(vec_buf_abc_b, b1);
      mirror1.gather(vec_buf_ab_b, b1);
      mirror2.gather(vec_buf_bc_b, b1);

      // gather from c
      mirror0.gather(vec_buf_abc_c, c1);
      mirror1.gather(vec_buf_bc_c, c1);
      mirror2.gather(vec_buf_ac_c, c1);

      // scatter to a
      mirror0.scatter_axpy(a1, vec_buf_abc_b);
      mirror0.scatter_axpy(a1, vec_buf_abc_c);
      mirror1.scatter_axpy(a1, vec_buf_ac_c);
      mirror2.scatter_axpy(a1, vec_buf_ab_b);

      // scatter to b
      mirror0.scatter_axpy(b1, vec_buf_abc_a);
      mirror0.scatter_axpy(b1, vec_buf_abc_c);
      mirror1.scatter_axpy(b1, vec_buf_ab_a);
      mirror2.scatter_axpy(b1, vec_buf_bc_c);

      // scatter to c
      mirror0.scatter_axpy(c1, vec_buf_abc_a);
      mirror0.scatter_axpy(c1, vec_buf_abc_b);
      mirror1.scatter_axpy(c1, vec_buf_bc_b);
      mirror2.scatter_axpy(c1, vec_buf_ac_a);

      // check against reference results
      TEST_CHECK_LESS_THAN(a1.max_rel_diff(a2), tol);
      TEST_CHECK_LESS_THAN(b1.max_rel_diff(b2), tol);
      TEST_CHECK_LESS_THAN(c1.max_rel_diff(c2), tol);
    } // End of test for DenseVectorBlocked
  }
};

SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, float, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, double, std::uint64_t, PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(VectorMirrorDenseVectorBlockedTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif

/**
 * \brief Test class for VectorMirror class template
 *
 * \author Peter Zajac
 */
template<
  typename DT_,
  typename IT_>
class VectorMirrorSparseVectorBlockedTest
  : public UnitTest
{
  static constexpr int block_size = 3;

  typedef DenseVector<DT_, IT_> BufferVectorType;
  typedef SparseVectorBlocked<DT_, IT_, block_size> BlockedSparseVectorType;
  typedef VectorMirror<DT_, IT_> MirrorType;

public:
  explicit VectorMirrorSparseVectorBlockedTest(PreferredBackend backend)
    : UnitTest("VectorMirrorSparseVectorBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    // create mirrors
    MirrorType mirror0(Index(3), Index(1));
    MirrorType mirror1(Index(3), Index(1));
    MirrorType mirror2(Index(3), Index(1));
    mirror0.indices_view_w()[0] = Index(0);
    mirror1.indices_view_w()[0] = Index(1);
    mirror2.indices_view_w()[0] = Index(2);

    // Tests for SparseVectorBlocked
    {
      // a1(2)+b1(1)
      BufferVectorType vec_buf_ab_a(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_ab_b(Index(block_size)*Index(1), DT_(0));
      // a1(1)+c1(2)
      BufferVectorType vec_buf_ac_a(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_ac_c(Index(block_size)*Index(1), DT_(0));
      // b1(2)+c(1)
      BufferVectorType vec_buf_bc_b(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_bc_c(Index(block_size)*Index(1), DT_(0));
      // a1(0) + b1(0) + c1(0)
      BufferVectorType vec_buf_abc_a(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_abc_b(Index(block_size)*Index(1), DT_(0));
      BufferVectorType vec_buf_abc_c(Index(block_size)*Index(1), DT_(0));

      // clear buffers
      vec_buf_abc_a.format(DT_(0));
      vec_buf_abc_b.format(DT_(0));
      vec_buf_abc_c.format(DT_(0));
      vec_buf_ab_a.format(DT_(0));
      vec_buf_ab_b.format(DT_(0));
      vec_buf_bc_b.format(DT_(0));
      vec_buf_bc_c.format(DT_(0));
      vec_buf_ac_a.format(DT_(0));
      vec_buf_ac_c.format(DT_(0));

      BlockedSparseVectorType a1(Index(3), Index(2));
      BlockedSparseVectorType a2(Index(3), Index(2));
      BlockedSparseVectorType b1(Index(3), Index(2));
      BlockedSparseVectorType b2(Index(3), Index(2));
      BlockedSparseVectorType c1(Index(3), Index(2));
      BlockedSparseVectorType c2(Index(3), Index(2));

      typedef typename BlockedSparseVectorType::ValueType ValueType;

      // Vectors to gather/scatter
      {
        Memory::TypedView<IT_> a1i = a1.indices_view_w();
        Memory::TypedView<ValueType> a1v = a1.elements_view_w();
        a1i[0] = IT_(0);
        a1i[1] = IT_(1);
        a1v[0][0] = DT_(-3);
        a1v[0][1] = DT_(0.5);
        a1v[0][2] = -DT_(1);
        a1v[1][0] = DT_(-7);
        a1v[1][1] = DT_(5);
        a1v[1][2] = -DT_(1000);
      }

      {
        Memory::TypedView<IT_> b1i = b1.indices_view_w();
        Memory::TypedView<ValueType> b1v = b1.elements_view_w();
        b1i[0] = IT_(1);
        b1i[1] = IT_(2);
        b1v[0][0] = DT_(-3.7);
        b1v[0][1] = DT_(8.001);
        b1v[0][2] = -DT_(9.002);
        b1v[1][0] = DT_(19);
        b1v[1][1] = DT_(-111);
        b1v[1][2] = -DT_(111);
      }

      {
        Memory::TypedView<IT_> c1i = c1.indices_view_w();
        Memory::TypedView<ValueType> c1v = c1.elements_view_w();
        c1i[0] = IT_(0);
        c1i[1] = IT_(2);
        c1v[0][0] = DT_(2303);
        c1v[0][1] = DT_(0);
        c1v[0][2] = DT_(0);
        c1v[1][0] = DT_(0.001);
        c1v[1][1] = DT_(1);
        c1v[1][2] = -DT_(7.7);
      }

      // Supposed results. See the gather operations below for this to make sense.
      {
        Memory::TypedView<IT_> a2i = a2.indices_view_w();
        Memory::TypedView<ValueType> a2v = a2.elements_view_w();
        a2i[0] = IT_(0);
        a2i[1] = IT_(1);
        a2v[0][0] = DT_(-3)  + DT_(2303);
        a2v[0][1] = DT_(0.5) + DT_(0);
        a2v[0][2] = -DT_(1)  + DT_(0);
        a2v[1][0] = DT_(-7)    + DT_(0.001);
        a2v[1][1] = DT_(5)     + DT_(1);
        a2v[1][2] = -DT_(1000) - DT_(7.7);
      }
      //a2(Index(0), a1(0)+b1(0)+c1(0));
      //a2(Index(1), a1(1)+c1(2));

      {
        Memory::TypedView<IT_> b2i = b2.indices_view_w();
        Memory::TypedView<ValueType> b2v = b2.elements_view_w();
        b2i[0] = IT_(1);
        b2i[1] = IT_(2);
        b2v[0][0] = DT_(-3.7);
        b2v[0][1] = DT_(8.001);
        b2v[0][2] = -DT_(9.002);
        b2v[1][0] = DT_(19);
        b2v[1][1] = DT_(-111);
        b2v[1][2] = -DT_(111);
      }
      //b2(Index(1), a1(2)+b1(1));
      //b2(Index(2), b1(2)+c1(1));

      {
        Memory::TypedView<IT_> c2i = c2.indices_view_w();
        Memory::TypedView<ValueType> c2v = c2.elements_view_w();
        c2i[0] = IT_(0);
        c2i[1] = IT_(2);
        c2v[0][0] = DT_(-3)  + DT_(2303);
        c2v[0][1] = DT_(0.5) + DT_(0);
        c2v[0][2] = -DT_(1)  + DT_(0);
        c2v[1][0] = DT_(0.001) + DT_(-7);
        c2v[1][1] = DT_(1)     + DT_(5);
        c2v[1][2] = -DT_(7.7)  - DT_(1000);
      }
      //c2(Index(0), a1(0)+b1(0)+c1(0));
      //c2(Index(2), a1(1)+c1(2));

      // gather from a
      mirror0.gather(vec_buf_abc_a, a1);
      mirror1.gather(vec_buf_ac_a, a1);
      mirror2.gather(vec_buf_ab_a, a1);

      // gather from b
      mirror0.gather(vec_buf_abc_b, b1);
      mirror1.gather(vec_buf_ab_b, b1);
      mirror2.gather(vec_buf_bc_b, b1);

      // gather from c
      mirror0.gather(vec_buf_abc_c, c1);
      mirror1.gather(vec_buf_bc_c, c1);
      mirror2.gather(vec_buf_ac_c, c1);

      // scatter to a
      mirror0.scatter_axpy(a1, vec_buf_abc_b);
      mirror0.scatter_axpy(a1, vec_buf_abc_c);
      mirror1.scatter_axpy(a1, vec_buf_ac_c);
      mirror2.scatter_axpy(a1, vec_buf_ab_b);

      // scatter to b
      mirror0.scatter_axpy(b1, vec_buf_abc_a);
      mirror0.scatter_axpy(b1, vec_buf_abc_c);
      mirror1.scatter_axpy(b1, vec_buf_ab_a);
      mirror2.scatter_axpy(b1, vec_buf_bc_c);

      // scatter to c
      mirror0.scatter_axpy(c1, vec_buf_abc_a);
      mirror0.scatter_axpy(c1, vec_buf_abc_b);
      mirror1.scatter_axpy(c1, vec_buf_bc_b);
      mirror2.scatter_axpy(c1, vec_buf_ac_a);

      TEST_CHECK_EQUAL(a1.num_nzes(), a2.num_nzes());
      TEST_CHECK(a1.same_layout(a2));
      TEST_CHECK_LESS_THAN(a1.max_rel_diff(a2), tol);

      TEST_CHECK_EQUAL(b1.num_nzes(), b2.num_nzes());
      TEST_CHECK(b1.same_layout(b2));
      TEST_CHECK_LESS_THAN(b1.max_rel_diff(b2), tol);

      TEST_CHECK_EQUAL(c1.num_nzes(), c2.num_nzes());
      TEST_CHECK(c1.same_layout(c2));
      TEST_CHECK_LESS_THAN(c1.max_rel_diff(c2), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, float, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, double, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, double, std::uint64_t, PreferredBackend::cuda);
//#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, Half, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(VectorMirrorSparseVectorBlockedTest, Half, std::uint64_t, PreferredBackend::cuda);
//#endif
//#endif
