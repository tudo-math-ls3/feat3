// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
class VectorMirrorTest
  : public UnitTest
{
  static constexpr int block_size = 3;

  typedef DenseVector<DT_, IT_> BufferVectorType;

  typedef DenseVector<DT_, IT_> VectorType;
  typedef SparseVector<DT_, IT_> SparseVectorType;
  typedef DenseVectorBlocked<DT_, IT_, block_size> BlockedVectorType;
  typedef SparseVectorBlocked<DT_, IT_, block_size> BlockedSparseVectorType;


  //typedef DenseVector<IT_, IT_> IVectorType;
  typedef SparseMatrixCSR<DT_, IT_> MatrixType;
  typedef VectorMirror<DT_, IT_> MirrorType;

public:
  VectorMirrorTest(PreferredBackend backend)
    : UnitTest("VectorMirrorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~VectorMirrorTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

    // create mirrors
    MirrorType mirror0(Index(3), Index(1));
    MirrorType mirror1(Index(3), Index(1));
    MirrorType mirror2(Index(3), Index(1));
    mirror0.indices()[0] = Index(0);
    mirror1.indices()[0] = Index(1);
    mirror2.indices()[0] = Index(2);

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
      a1(Index(0), DT_(4));
      a1(Index(1), DT_(7));
      a1(Index(2), DT_(2));
      b1(Index(0), DT_(2));
      b1(Index(1), DT_(1));
      b1(Index(2), DT_(3));
      c1(Index(0), DT_(1));
      c1(Index(1), DT_(5));
      c1(Index(2), DT_(2));

      // initialize global vectors
      a2(Index(0), DT_(7));
      a2(Index(1), DT_(9));
      a2(Index(2), DT_(3));
      b2(Index(0), DT_(7));
      b2(Index(1), DT_(3));
      b2(Index(2), DT_(8));
      c2(Index(0), DT_(7));
      c2(Index(1), DT_(8));
      c2(Index(2), DT_(9));

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

      // subtract reference results
      a1.axpy(a2, -DT_(1));
      b1.axpy(b2, -DT_(1));
      c1.axpy(c2, -DT_(1));

      // check norms
      TEST_CHECK_EQUAL_WITHIN_EPS(a1.norm2(), DT_(0), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(b1.norm2(), DT_(0), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(c1.norm2(), DT_(0), tol);
    } // End of test for DenseVector

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

      SparseVectorType a1(Index(3));
      SparseVectorType a2(Index(3));
      SparseVectorType b1(Index(3));
      SparseVectorType b2(Index(3));
      SparseVectorType c1(Index(3));
      SparseVectorType c2(Index(3));

      // Vectors to gather/scatter. Note that we use the emplacement operator which also sets the sparsity pattern.
      a1(Index(0), DT_(4));
      a1(Index(1), DT_(7));

      b1(Index(1), DT_(1));
      b1(Index(2), DT_(3));

      c1(Index(0), DT_(1));
      c1(Index(2), DT_(2));

      // Supposed results. See the gather operations below for this to make sense.
      // Note that each sparse vector is missing one entry and this defines the "sparsity" pattern
      a2(Index(0), a1(0)+b1(0)+c1(0));
      a2(Index(1), a1(1)+c1(2));

      b2(Index(1), a1(2)+b1(1));
      b2(Index(2), b1(2)+c1(1));

      c2(Index(0), a1(0)+b1(0)+c1(0));
      c2(Index(2), a1(1)+c1(2));

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

      TEST_CHECK_MSG(a1.used_elements() == a2.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < a1.used_elements(); ++i)
      {
        Index i1(a1.indices()[i]);
        Index i2(a2.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS(a1.elements()[i], a2.elements()[i], tol);
      }

      TEST_CHECK_MSG(b1.used_elements() == b2.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < b1.used_elements(); ++i)
      {
        Index i1(b1.indices()[i]);
        Index i2(b2.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS(b1.elements()[i], b2.elements()[i], tol);
      }

      TEST_CHECK_MSG(c1.used_elements() == c2.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < c1.used_elements(); ++i)
      {
        Index i1(c1.indices()[i]);
        Index i2(c2.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS(c1.elements()[i], c2.elements()[i], tol);
      }

    } // End of test for SparseVector

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

      Tiny::Vector<DT_, block_size> tiny_tmp(DT_(0));

      // Vectors to gather/scatter
      tiny_tmp(0) = DT_(-3), tiny_tmp(1) = DT_(0.5), tiny_tmp(block_size-1) = -DT_(1);
      a1(Index(0), tiny_tmp);
      tiny_tmp(0) = DT_(-7), tiny_tmp(1) = DT_(5), tiny_tmp(block_size-1) = -DT_(1000);
      a1(Index(1), tiny_tmp);
      tiny_tmp(0) = DT_(-0.3), tiny_tmp(1) = DT_(0.1), tiny_tmp(block_size-1) = DT_(2);
      a1(Index(2), tiny_tmp);

      tiny_tmp(0) = DT_(0.3), tiny_tmp(1) = DT_(0.001), tiny_tmp(block_size-1) = DT_(2.2);
      b1(Index(0), tiny_tmp);
      tiny_tmp(0) = DT_(-3.7), tiny_tmp(1) = DT_(8.001), tiny_tmp(block_size-1) = -DT_(9.002);
      b1(Index(1), tiny_tmp);
      tiny_tmp(0) = DT_(19), tiny_tmp(1) = DT_(-111), tiny_tmp(block_size-1) = -DT_(111);
      b1(Index(2), tiny_tmp);

      tiny_tmp(0) = DT_(2303), tiny_tmp(1) = DT_(0), tiny_tmp(block_size-1) = DT_(0);
      c1(Index(0), tiny_tmp);
      tiny_tmp(0) = DT_(0), tiny_tmp(1) = DT_(1), tiny_tmp(block_size-1) = DT_(7);
      c1(Index(1), tiny_tmp);
      tiny_tmp(0) = DT_(0.001), tiny_tmp(1) = DT_(1), tiny_tmp(block_size-1) = -DT_(7.7);
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

      // Subtract reference results
      a1.axpy(a2, -DT_(1));
      b1.axpy(b2, -DT_(1));
      c1.axpy(c2, -DT_(1));

      // Check norms
      TEST_CHECK_EQUAL_WITHIN_EPS(a1.norm2(), DT_(0), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(b1.norm2(), DT_(0), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(c1.norm2(), DT_(0), tol);
    } // End of test for DenseVectorBlocked

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

      BlockedSparseVectorType a1(Index(3));
      BlockedSparseVectorType a2(Index(3));
      BlockedSparseVectorType b1(Index(3));
      BlockedSparseVectorType b2(Index(3));
      BlockedSparseVectorType c1(Index(3));
      BlockedSparseVectorType c2(Index(3));

      Tiny::Vector<DT_, block_size> tiny_tmp(DT_(0));

      // Vectors to gather/scatter
      tiny_tmp(0) = DT_(-3), tiny_tmp(1) = DT_(0.5), tiny_tmp(block_size-1) = -DT_(1);
      a1(Index(0), tiny_tmp);
      tiny_tmp(0) = DT_(-7), tiny_tmp(1) = DT_(5), tiny_tmp(block_size-1) = -DT_(1000);
      a1(Index(1), tiny_tmp);

      tiny_tmp(0) = DT_(-3.7), tiny_tmp(1) = DT_(8.001), tiny_tmp(block_size-1) = -DT_(9.002);
      b1(Index(1), tiny_tmp);
      tiny_tmp(0) = DT_(19), tiny_tmp(1) = DT_(-111), tiny_tmp(block_size-1) = -DT_(111);
      b1(Index(2), tiny_tmp);

      tiny_tmp(0) = DT_(2303), tiny_tmp(1) = DT_(0), tiny_tmp(block_size-1) = DT_(0);
      c1(Index(0), tiny_tmp);
      tiny_tmp(0) = DT_(0.001), tiny_tmp(1) = DT_(1), tiny_tmp(block_size-1) = -DT_(7.7);
      c1(Index(2), tiny_tmp);

      // Supposed results. See the gather operations below for this to make sense.
      a2(Index(0), a1(0)+b1(0)+c1(0));
      a2(Index(1), a1(1)+c1(2));

      b2(Index(1), a1(2)+b1(1));
      b2(Index(2), b1(2)+c1(1));

      c2(Index(0), a1(0)+b1(0)+c1(0));
      c2(Index(2), a1(1)+c1(2));

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

      TEST_CHECK_MSG(a1.used_elements() == a2.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < a1.used_elements(); ++i)
      {
        Index i1(a1.indices()[i]);
        Index i2(a2.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS((a1.elements()[i]-a2.elements()[i]).norm_euclid(), DT_(0), tol);
      }

      TEST_CHECK_MSG(b1.used_elements() == b2.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < b1.used_elements(); ++i)
      {
        Index i1(b1.indices()[i]);
        Index i2(b2.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS((b1.elements()[i]-b2.elements()[i]).norm_euclid(), DT_(0), tol);
      }

      TEST_CHECK_MSG(c1.used_elements() == c2.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < c1.used_elements(); ++i)
      {
        Index i1(c1.indices()[i]);
        Index i2(c2.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS((c1.elements()[i]-c2.elements()[i]).norm_euclid(), DT_(0), tol);
      }
    }
  }
};

VectorMirrorTest <double, std::uint32_t> vector_mirror_test_main_double_uint32(PreferredBackend::generic);
VectorMirrorTest <float, std::uint32_t> vector_mirror_test_main_float_uint32(PreferredBackend::generic);
VectorMirrorTest <double, std::uint64_t> vector_mirror_test_main_double_uint64(PreferredBackend::generic);
VectorMirrorTest <float, std::uint64_t> vector_mirror_test_main_float_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
VectorMirrorTest <float, std::uint64_t> mkl_vector_mirror_test_float_uint64(PreferredBackend::mkl);
VectorMirrorTest <double, std::uint64_t> mkl_vector_mirror_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
VectorMirrorTest <__float128, std::uint64_t> vector_mirror_test_float128_uint64(PreferredBackend::generic);
VectorMirrorTest <__float128, std::uint32_t> vector_mirror_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
VectorMirrorTest <Half, std::uint32_t> vector_mirror_test_half_uint32(PreferredBackend::generic);
VectorMirrorTest <Half, std::uint64_t> vector_mirror_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
VectorMirrorTest <float, std::uint32_t> cuda_vector_mirror_test_float_uint32(PreferredBackend::cuda);
VectorMirrorTest <double, std::uint32_t> cuda_vector_mirror_test_double_uint32(PreferredBackend::cuda);
VectorMirrorTest <float, std::uint64_t> cuda_vector_mirror_test_float_uint64(PreferredBackend::cuda);
VectorMirrorTest <double, std::uint64_t> cuda_vector_mirror_test_double_uint64(PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
VectorMirrorTest <Half, std::uint32_t> cuda_vector_mirror_test_half_uint32(PreferredBackend::cuda);
VectorMirrorTest <Half, std::uint64_t> cuda_vector_mirror_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
