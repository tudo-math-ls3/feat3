#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
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
  typename MemType_,
  typename DT_,
  typename IT_>
class VectorMirrorTest
  : public FullTaggedTest<MemType_, DT_, IT_>
{
  static constexpr int block_size = 3;

  typedef DenseVector<Mem::Main, DT_, IT_> BufferVectorType;

  typedef DenseVector<MemType_, DT_, IT_> VectorType;
  typedef SparseVector<MemType_, DT_, IT_> SparseVectorType;
  typedef DenseVectorBlocked<MemType_, DT_, IT_, block_size> BlockedVectorType;
  typedef SparseVectorBlocked<MemType_, DT_, IT_, block_size> BlockedSparseVectorType;


  //typedef DenseVector<MemType_, IT_, IT_> IVectorType;
  typedef SparseMatrixCSR<MemType_, DT_, IT_> MatrixType;
  typedef VectorMirror<MemType_, DT_, IT_> MirrorType;

public:
  VectorMirrorTest()
    : FullTaggedTest<MemType_, DT_, IT_>("VectorMirrorTest")
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
      a1.axpy(a2, a1, -DT_(1));
      b1.axpy(b2, b1, -DT_(1));
      c1.axpy(c2, c1, -DT_(1));

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

      // There is no axpy for SparseVector yet, so for now download the vectors (if necessary) and do it by hand.
      LAFEM::SparseVector<Mem::Main, DT_, IT_> a1_main; a1_main.convert(a1);
      LAFEM::SparseVector<Mem::Main, DT_, IT_> a2_main; a2_main.convert(a2);

      TEST_CHECK_MSG(a1_main.used_elements() == a2_main.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < a1_main.used_elements(); ++i)
      {
        Index i1(a1_main.indices()[i]);
        Index i2(a2_main.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS(a1_main.elements()[i], a2_main.elements()[i], tol);
      }

      LAFEM::SparseVector<Mem::Main, DT_, IT_> b1_main; b1_main.convert(b1);
      LAFEM::SparseVector<Mem::Main, DT_, IT_> b2_main; b2_main.convert(b2);

      TEST_CHECK_MSG(b1_main.used_elements() == b2_main.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < b1_main.used_elements(); ++i)
      {
        Index i1(b1_main.indices()[i]);
        Index i2(b2_main.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS(b1_main.elements()[i], b2_main.elements()[i], tol);
      }

      LAFEM::SparseVector<Mem::Main, DT_, IT_> c1_main; c1_main.convert(c1);
      LAFEM::SparseVector<Mem::Main, DT_, IT_> c2_main; c2_main.convert(c2);

      TEST_CHECK_MSG(c1_main.used_elements() == c2_main.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < c1_main.used_elements(); ++i)
      {
        Index i1(c1_main.indices()[i]);
        Index i2(c2_main.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS(c1_main.elements()[i], c2_main.elements()[i], tol);
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
      a1.axpy(a2, a1, -DT_(1));
      b1.axpy(b2, b1, -DT_(1));
      c1.axpy(c2, c1, -DT_(1));

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

      // There is no axpy for SparseVector yet, so for now download the vectors (if necessary) and do it by hand.
      LAFEM::SparseVectorBlocked<Mem::Main, DT_, IT_, block_size> a1_main; a1_main.convert(a1);
      LAFEM::SparseVectorBlocked<Mem::Main, DT_, IT_, block_size> a2_main; a2_main.convert(a2);

      TEST_CHECK_MSG(a1_main.used_elements() == a2_main.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < a1_main.used_elements(); ++i)
      {
        Index i1(a1_main.indices()[i]);
        Index i2(a2_main.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS((a1_main.elements()[i]-a2_main.elements()[i]).norm_euclid(), DT_(0), tol);
      }

      LAFEM::SparseVectorBlocked<Mem::Main, DT_, IT_, block_size> b1_main; b1_main.convert(b1);
      LAFEM::SparseVectorBlocked<Mem::Main, DT_, IT_, block_size> b2_main; b2_main.convert(b2);

      TEST_CHECK_MSG(b1_main.used_elements() == b2_main.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < b1_main.used_elements(); ++i)
      {
        Index i1(b1_main.indices()[i]);
        Index i2(b2_main.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS((b1_main.elements()[i]-b2_main.elements()[i]).norm_euclid(), DT_(0), tol);
      }

      LAFEM::SparseVectorBlocked<Mem::Main, DT_, IT_, block_size> c1_main; c1_main.convert(c1);
      LAFEM::SparseVectorBlocked<Mem::Main, DT_, IT_, block_size> c2_main; c2_main.convert(c2);

      TEST_CHECK_MSG(c1_main.used_elements() == c2_main.used_elements(),"Wrong number of nonzeros.");
      for(Index i(0); i < c1_main.used_elements(); ++i)
      {
        Index i1(c1_main.indices()[i]);
        Index i2(c2_main.indices()[i]);
        TEST_CHECK_MSG(i1 == i2,"Error in sparsity pattern.");
        TEST_CHECK_EQUAL_WITHIN_EPS((c1_main.elements()[i]-c2_main.elements()[i]).norm_euclid(), DT_(0), tol);
      }
    }
  }
};

VectorMirrorTest<Mem::Main, double, unsigned long> vector_mirror_test_main_d_ul;
VectorMirrorTest<Mem::Main, float, unsigned int> vector_mirror_test_main_f_ui;
/// \todo Add cuda vector mirror tests
//VectorMirrorTest<Mem::CUDA, double, unsigned long> vector_mirror_test_cuda_d_ul;
//VectorMirrorTest<Mem::CUDA, float, unsigned int> vector_mirror_test_cuda_f_ui;
