// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse matrix csr blocked class.
 *
 * \test test description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  void test_vector_types() const
  {
    // define a hand-full of sparse matrix BCSR
    SparseMatrixBCSR<DT_, IT_, 3, 3> bcsr_3x3;
    SparseMatrixBCSR<DT_, IT_, 3, 1> bcsr_3x1;
    SparseMatrixBCSR<DT_, IT_, 1, 3> bcsr_1x3;
    SparseMatrixBCSR<DT_, IT_, 1, 1> bcsr_1x1;

    // now create the left/right vectors
    TEST_CHECK_EQUAL(bcsr_3x3.create_vector_l().name(), "DenseVectorBlocked");
    TEST_CHECK_EQUAL(bcsr_3x3.create_vector_r().name(), "DenseVectorBlocked");
    TEST_CHECK_EQUAL(bcsr_3x1.create_vector_l().name(), "DenseVectorBlocked");
    TEST_CHECK_EQUAL(bcsr_3x1.create_vector_r().name(), "DenseVector");
    TEST_CHECK_EQUAL(bcsr_1x3.create_vector_l().name(), "DenseVector");
    TEST_CHECK_EQUAL(bcsr_1x3.create_vector_r().name(), "DenseVectorBlocked");
    TEST_CHECK_EQUAL(bcsr_1x1.create_vector_l().name(), "DenseVector");
    TEST_CHECK_EQUAL(bcsr_1x1.create_vector_r().name(), "DenseVector");
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();
    test_vector_types();

    SparseMatrixBCSR<DT_, IT_, 2, 3> zero;
    TEST_CHECK(zero.empty());

    DenseVector<DT_, IT_> dv1(12);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1)/DT_(7*i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(2, 2, dv3, dv2, dv1);

    TEST_CHECK(!c.empty());
    TEST_CHECK_EQUAL(c.val_view_r()(1)(1,1), DT_(10+1)/DT_(7*10+1));

    SparseMatrixBCSR<DT_, IT_, 2, 3> d;
    d.convert(c);
    TEST_CHECK_EQUAL(d.num_rows(), c.num_rows());
    TEST_CHECK_EQUAL(d.num_cols(), c.num_cols());
    TEST_CHECK_EQUAL(d.num_nzes(), c.num_nzes());
    TEST_CHECK_LESS_THAN(d.max_rel_diff(c), eps);
    TEST_CHECK_EQUAL(d.val_arbiter(), c.val_arbiter());
    TEST_CHECK_EQUAL(d.row_ptr_arbiter(), c.row_ptr_arbiter());
    SparseMatrixBCSR<DT_, IT_, 2, 3> e;
    e.clone(c);
    TEST_CHECK_LESS_THAN(e.max_rel_diff(c), eps);
    TEST_CHECK_NOT_EQUAL(e.val_arbiter(), c.val_arbiter());
    TEST_CHECK_EQUAL(e.row_ptr_arbiter(), c.row_ptr_arbiter());
    e = c.clone(CloneMode::Deep);
    TEST_CHECK_LESS_THAN(e.max_rel_diff(c), eps);
    TEST_CHECK_NOT_EQUAL(e.val_arbiter(), c.val_arbiter());
    TEST_CHECK_NOT_EQUAL(e.row_ptr_arbiter(), c.row_ptr_arbiter());

    SparseMatrixBCSR<DT_, IT_, 2, 3> f(c.layout());
    TEST_CHECK_EQUAL(f.num_rows(), c.num_rows());
    TEST_CHECK_EQUAL(f.num_cols(), c.num_cols());
    TEST_CHECK_EQUAL(f.num_nzes(), c.num_nzes());
    TEST_CHECK_NOT_EQUAL(f.val_arbiter(), c.val_arbiter());
    TEST_CHECK_EQUAL(f.row_ptr_arbiter(), c.row_ptr_arbiter());
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRSerializeTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRSerializeTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    DenseVector<DT_, IT_> dv1(12);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1)/DT_(7*i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(2, 2, dv3, dv2, dv1);

    BinaryStream bs;
    c.write_out(FileMode::fm_bcsr, bs);
    bs.seekg(0);
    SparseMatrixBCSR<DT_, IT_, 2, 3> g(FileMode::fm_bcsr, bs);
    TEST_CHECK_LESS_THAN(g.max_rel_diff(c), eps);

    //std::stringstream ts;
    //f.write_out(FileMode::fm_mtx, ts);
    //SparseMatrixCSR<DT_, IT_> j(FileMode::fm_mtx, ts);
    //TEST_CHECK_EQUAL(j, f);

    //std::stringstream ts2;
    //f.write_out_mtx(ts2, true);
    //SparseMatrixCSR<DT_, IT_> j2(FileMode::fm_mtx, ts2);
    //TEST_CHECK_EQUAL(j2, f);

    auto kp = c.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixBCSR<DT_, IT_, 2, 3> k(kp);
    TEST_CHECK_LESS_THAN(k.max_rel_diff(c), eps);

#ifdef FEAT_HAVE_ZLIB
    auto zl = c.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixBCSR<DT_, IT_, 2, 3> zlib(zl);
    TEST_CHECK_LESS_THAN(k.max_rel_diff(c), eps);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = c.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixBCSR<DT_, IT_, 2, 3> zfp(zf);
    TEST_CHECK_LESS_THAN(zfp.max_rel_diff(c), DT_(1e-5));
#endif
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

/**
 * \brief Test class for the sparse matrix csr blocked apply method.
 *
 * \test test description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRApplyTest
  : public UnitTest
{
public:
   explicit SparseMatrixBCSRApplyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRApplyTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    DenseVector<DT_, IT_> dv1(12);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1)/DT_(7*i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(2, 2, dv3, dv2, dv1);

    DenseVector<DT_, IT_> x(c.num_cols_raw());
    DenseVector<DT_, IT_> y(c.num_rows_raw());
    DenseVector<DT_, IT_> r(c.num_rows_raw());
    DenseVector<DT_, IT_> ref(c.num_rows_raw());
    {
      Memory::TypedView<DT_> vx = x.elements_view_w();
      for (Index i(0) ; i < x.size() ; ++i)
      {
        vx[i] = DT_(i);
      }
    }
    {
      Memory::TypedView<DT_> vy = y.elements_view_w();
      Memory::TypedView<DT_> vr = r.elements_view_w();
      for (Index i(0) ; i < r.size() ; ++i)
      {
        vr[i] = DT_(4711);
        vy[i] = DT_(i % 100);
      }
    }
    DenseVectorBlocked<DT_, IT_, 3> xb(x);
    DenseVectorBlocked<DT_, IT_, 2> yb(y);
    DenseVectorBlocked<DT_, IT_, 2> rb(r);

    SparseMatrixCSR<DT_, IT_> csr;
    csr.convert(c);
    csr.apply(ref, x);

    c.apply(r, x);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, x);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(r, xb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    //defect
    DT_ alpha(-1);
    csr.apply(ref, x, y, alpha);
    c.apply(r, x, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, x, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(r, xb, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb, y, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    //axpy
    alpha = DT_(1.234);
    csr.apply(ref, x, y, alpha);
    c.apply(r, x, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    // &r == &y
    r.copy(y);
    c.apply(r, x, r, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, x, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(r, xb, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb, y, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

/**
 * \brief Test class for the sparse matrix csr blocked apply transposed method.
 *
 * \test test description missing
 *
 * \author Pia Ritter
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRApplyTransposedTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRApplyTransposedTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRApplyTransposedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();
    DenseVector<DT_, IT_> dv1(12);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1)/DT_(7*i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 3, 2> c(2, 2, dv3, dv2, dv1);

    DenseVector<DT_, IT_> x(c.num_rows_raw());
    DenseVector<DT_, IT_> y(c.num_cols_raw());
    DenseVector<DT_, IT_> r(c.num_cols_raw());
    DenseVector<DT_, IT_> ref(c.num_cols_raw());
    {
      Memory::TypedView<DT_> vx = x.elements_view_w();
      for (Index i(0) ; i < x.size() ; ++i)
      {
        vx[i] = DT_(i);
      }
    }
    {
      Memory::TypedView<DT_> vy = y.elements_view_w();
      Memory::TypedView<DT_> vr = r.elements_view_w();
      Memory::TypedView<DT_> vref = ref.elements_view_w();
      for (Index i(0) ; i < r.size() ; ++i)
      {
        vr[i] = DT_(4711);
        vref[i] = DT_(4711);
        vy[i] = DT_(i % 100);
      }
    }

    DenseVectorBlocked<DT_, IT_, 3> xb(x);
    DenseVectorBlocked<DT_, IT_, 2> yb(y);
    DenseVectorBlocked<DT_, IT_, 2> rb(r);

    SparseMatrixCSR<DT_, IT_> csr;
    csr.convert(c);
    csr.apply_transposed(ref, x);

    c.apply_transposed(r, x);

    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(rb, x);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(r, xb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(rb, xb);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    //defect
    DT_ alpha(-1);
    csr.apply_transposed(ref, x, y, alpha);
    c.apply_transposed(r, x, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(rb, x, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(r, xb, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(rb, xb, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(rb, xb, y, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    //axpy
    alpha = DT_(1.234);
    csr.apply_transposed(ref, x, y, alpha);
    c.apply_transposed(r, x, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    // &r == &y
    r.copy(y);
    c.apply_transposed(r, x, r, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(rb, x, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(r, xb, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(rb, xb, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply_transposed(rb, xb, y, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, double, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplyTransposedTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

/**
 * \brief Test class for the sparse matrix csr blocked apply method.
 *
 * \test test description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRApplySquareTest
  : public UnitTest
{
public:
   explicit SparseMatrixBCSRApplySquareTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRApplySquareTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRApplySquareTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    DenseVector<DT_, IT_> dv1(18);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1)/DT_(7*i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 3, 3> c(2, 2, dv3, dv2, dv1);

    DenseVector<DT_, IT_> x(c.num_cols_raw());
    DenseVector<DT_, IT_> y(c.num_rows_raw());
    DenseVector<DT_, IT_> r(c.num_rows_raw());
    DenseVector<DT_, IT_> ref(c.num_rows_raw());
    {
      Memory::TypedView<DT_> vx = x.elements_view_w();
      for (Index i(0) ; i < x.size() ; ++i)
      {
        vx[i] = DT_(i);
      }
    }
    {
      Memory::TypedView<DT_> vy = y.elements_view_w();
      Memory::TypedView<DT_> vr = r.elements_view_w();
      Memory::TypedView<DT_> vref = ref.elements_view_w();
      for (Index i(0) ; i < r.size() ; ++i)
      {
        vr[i] = DT_(4711);
        vref[i] = DT_(4711);
        vy[i] = DT_(i % 100);
      }
    }
    DenseVectorBlocked<DT_, IT_, 3> xb(x);
    DenseVectorBlocked<DT_, IT_, 3> yb(y);
    DenseVectorBlocked<DT_, IT_, 3> rb(r);

    SparseMatrixCSR<DT_, IT_> csr;
    csr.convert(c);
    csr.apply(ref, x);

    c.apply(r, x);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, x);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(r, xb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    // defect
    DT_ alpha(-1);
    csr.apply(ref, x, y, alpha);
    c.apply(r, x, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    // &r == &y
    r.copy(y);
    c.apply(r, x, r, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, x, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(r, xb, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb, y, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    // axpy
    alpha = DT_(1.234);
    csr.apply(ref, x, y, alpha);
    c.apply(r, x, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    // &r == &y
    r.copy(y);
    c.apply(r, x, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, x, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(r, xb, y, alpha);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb, yb, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);

    c.apply(rb, xb, y, alpha);
    r.convert(rb);
    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref), eps);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRApplySquareTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRDiagTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRDiagTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRDiagTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    DenseVector<DT_, IT_> dv1(18);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 3, 3> c(2, 2, dv3, dv2, dv1);

    DenseVectorBlocked<DT_, IT_, 3> dref(2);
    {
      auto vd(dref.elements_view_w());
      vd[0][0] = DT_(1);
      vd[0][1] = DT_(5);
      vd[0][2] = DT_(9);
      vd[1][0] = DT_(10);
      vd[1][1] = DT_(14);
      vd[1][2] = DT_(18);
    }
    auto diag = c.extract_diag();
    TEST_CHECK_LESS_THAN(diag.max_rel_diff(dref), tol);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRDiagTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

/**
 * \brief Test class for the sparse matrix csr blocked scale method.
 *
 * \test test description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRScaleTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRScaleTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();
    DenseVector<DT_, IT_> dv1(12);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1)/DT_(7*i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 2, 3> a(2, 2, dv3, dv2, dv1);
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(a.layout());

    DT_ scal = DT_(4711);

    SparseMatrixCSR<DT_, IT_> a_csr;
    a_csr.convert(a);
    SparseMatrixCSR<DT_, IT_> ref(a_csr.layout());
    SparseMatrixCSR<DT_, IT_> result_csr;
    ref.scale(a_csr, scal);

    c.scale(a, scal);
    result_csr.convert(c);
    TEST_CHECK_LESS_THAN(result_csr.max_rel_diff(ref), eps);

    a.scale(a, scal);
    result_csr.convert(a);
    TEST_CHECK_LESS_THAN(result_csr.max_rel_diff(ref), eps);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRScaleRowColTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRScaleRowColTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRScaleRowColTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ pi(Math::pi<DT_>());
    const DT_ tol = TestSystem::tol<DT_>();

    DenseVector<DT_, IT_> dv1(30);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1);
    }
    DenseVector<IT_, IT_> dv2(5);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
      v2[2] = IT_(0);
      v2[3] = IT_(2);
      v2[4] = IT_(3);
    }
    DenseVector<IT_, IT_> dv3(4);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(2);
      v3[2] = IT_(5);
      v3[3] = IT_(5);
    }
    SparseMatrixBCSR<DT_, IT_, 2, 3> a(3, 4, dv3, dv2, dv1);

    DenseVectorBlocked<DT_, IT_, 2> s(a.num_rows());
    {
      auto vs = s.elements_view_w();
      for (Index i(0); i < s.size(); ++i)
      {
        vs[i][0] = pi * DT_(i % 3 + 1) - DT_(5.21) + DT_(i);
        vs[i][1] = DT_(i) * DT_(12.31);
      }
    }
    DenseVectorBlocked<DT_, IT_, 3> t(a.num_cols());
    {
      auto vt = t.elements_view_w();
      for (Index i(0); i < t.size(); ++i)
      {
        vt[i][0] = pi * DT_(i % 3 + 1) - DT_(5.21) + DT_(i);
        vt[i][1] = DT_(i) * DT_(12.31);
        vt[i][2] = DT_(i) *DT_(-324.21) - DT_(13.37);
      }
    }

    SparseMatrixBCSR<DT_, IT_, 2, 3> ar = a.clone(LAFEM::CloneMode::Weak);
    SparseMatrixBCSR<DT_, IT_, 2, 3> ac = a.clone(LAFEM::CloneMode::Weak);

    {
      const Memory::TypedView<IT_> row_ptr = a.row_ptr_view_r();
      const Memory::TypedView<IT_> col_idx = a.col_idx_view_r();
      auto var = ar.val_view_rw();
      auto vac = ac.val_view_rw();

      for(Index i = 0; i < a.num_rows(); ++i)
      {
        const DT_ s0 = pi * DT_(i % 3 + 1) - DT_(5.21) + DT_(i);
        const DT_ s1 = DT_(i) * DT_(12.31);

        for(IT_ k = row_ptr[i]; k < row_ptr[i+1]; ++k)
        {
          IT_ j = col_idx[k];
          const DT_ t0 = pi * DT_(j % 3 + 1) - DT_(5.21) + DT_(j);
          const DT_ t1 = DT_(j) * DT_(12.31);
          const DT_ t2 = DT_(j) *DT_(-324.21) - DT_(13.37);

          var[k][0][0] *= s0;
          var[k][0][1] *= s0;
          var[k][0][2] *= s0;
          var[k][1][0] *= s1;
          var[k][1][1] *= s1;
          var[k][1][2] *= s1;

          vac[k][0][0] *= t0;
          vac[k][1][0] *= t0;
          vac[k][0][1] *= t1;
          vac[k][1][1] *= t1;
          vac[k][0][2] *= t2;
          vac[k][1][2] *= t2;
        }
      }
    }

    SparseMatrixBCSR<DT_, IT_, 2, 3> b(a.clone(LAFEM::CloneMode::Layout));

    b.scale_rows(a, s);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(ar), tol);

    b.scale_cols(a, t);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(ac), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRScaleRowColTest, double, std::uint32_t, PreferredBackend::cuda);
#endif


/**
 * \brief Test class for the sparse matrix csr blocked norm method.
 *
 * \test test description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRNormTest
  : UnitTest
{
public:
   explicit SparseMatrixBCSRNormTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRNormTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    DenseVector<DT_, IT_> dv1(12);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1)/DT_(7*i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 2, 3> a(2, 2, dv3, dv2, dv1);

    SparseMatrixCSR<DT_, IT_> a_csr;
    a_csr.convert(a);
    DT_ ref = a_csr.norm_frobenius();

    DT_ result = a.norm_frobenius();
    TEST_CHECK_RELATIVE(result, ref, tol);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRNormTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

/**
 * \brief Test class for the sparse matrix csr blocked axpy method.
 *
 * \test test description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRAxpyTest
  : public UnitTest
{
public:
   explicit SparseMatrixBCSRAxpyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRAxpyTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    DenseVector<DT_, IT_> dv1(12);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    DenseVector<DT_, IT_> dv4(12);
    {
      Memory::TypedView<DT_> v4(dv4.elements_view_w());
      for (Index i(0) ; i < dv4.size() ; ++i)
        v4[i] = DT_(i) - DT_(1);
    }

    SparseMatrixBCSR<DT_, IT_, 2, 3> a(2, 2, dv3, dv2, dv1);
    SparseMatrixBCSR<DT_, IT_, 2, 3> b(2, 2, dv3, dv2, dv4);
    SparseMatrixBCSR<DT_, IT_, 2, 3> ref(a.layout());
    SparseMatrixBCSR<DT_, IT_, 2, 3> ref2(a.layout());

    DT_ scal = DT_(1.234);

    {
      Memory::TypedView<DT_> vr1 = ref.val_view_raw_w();
      Memory::TypedView<DT_> vr2 = ref2.val_view_raw_w();
      Memory::TypedView<DT_> v1(dv1.elements_view_r());
      Memory::TypedView<DT_> v4(dv4.elements_view_r());
      for(Index i(0); i < ref.num_nzes_raw(); ++i)
      {
        vr1[i] = scal*v1(i) + v4(i);
        vr2[i] = scal*v4(i) + v4(i);
      }
    }

    // r != x
    a.scale(a, scal);
    a.axpy(b); /// \todo use axpby here
    TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), eps);

    // r == x
    b.axpy(b, scal);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(ref2), eps);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRAxpyTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

/**
 * \brief Test class for the sparse matrix csr blocked permute method.
 *
 * \test test description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRPermuteTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRPermuteTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRPermuteTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    DenseVector<DT_, IT_> dv1(18);
    {
      Memory::TypedView<DT_> v1(dv1.elements_view_w());
      for (Index i(0) ; i < dv1.size() ; ++i)
        v1[i] = DT_(i+1)/DT_(7*i+1);
    }
    DenseVector<IT_, IT_> dv2(2);
    {
      Memory::TypedView<IT_> v2(dv2.elements_view_w());
      v2[0] = IT_(0);
      v2[1] = IT_(1);
    }
    DenseVector<IT_, IT_> dv3(3);
    {
      Memory::TypedView<IT_> v3(dv3.elements_view_w());
      v3[0] = IT_(0);
      v3[1] = IT_(1);
      v3[2] = IT_(2);
    }
    SparseMatrixBCSR<DT_, IT_, 3, 3> a(2, 2, dv3, dv2, dv1);

    DenseVector<DT_, IT_> x(a.num_cols_raw());
    //DenseVector<DT_, IT_> r(a.num_rows_raw());
    {
      Memory::TypedView<DT_> vx = x.elements_view_w();
      for (Index i(0) ; i < x.size() ; ++i)
      {
        vx[i] = DT_(i);
      }
    }
    DenseVectorBlocked<DT_, IT_, 3> xb(x);
    DenseVectorBlocked<DT_, IT_, 3> rb(xb.size());

    /////////////////////////

    std::cout<<xb<<"\n";
    a.apply(rb, xb);
    DT_ ref_norm = rb.norm2();

    auto a_backup = a.clone(CloneMode::Deep);

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    Adjacency::Permutation perm(a.num_rows(), rng);

    a.permute(perm, perm);
    xb.permute(perm);

    std::cout<<xb<<"\n";
    a.apply(rb, xb);
    DT_ norm = rb.norm2();
    TEST_CHECK_RELATIVE(norm, ref_norm, eps);

    a = a_backup.clone(CloneMode::Deep);
    auto perm_inv = perm.inverse();
    a.permute(perm_inv, perm);
    a.permute(perm, perm_inv);
    TEST_CHECK_LESS_THAN(a.max_rel_diff(a_backup), eps);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRPermuteTest, double, std::uint32_t, PreferredBackend::cuda);
#endif


template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRMinMaxAbsElementTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRMinMaxAbsElementTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRMinMaxAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    SparseMatrixBCSR<DT_, IT_, 2, 3> a(2, 2, 12);
    {
      Memory::TypedView<DT_> val_raw(a.val_view_raw_w());
      Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
      Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
      row_ptr[0] = IT_(0);
      row_ptr[1] = IT_(1);
      col_idx[0] = IT_(0);
      col_idx[1] = IT_(1);
      col_idx[2] = IT_(2);
      for (Index i(0) ; i < a.num_nzes_raw() ; ++i)
        val_raw[i] = DT_(i+1) * (i%2 == 0 ? DT_(1) : DT_(-1));
    }

    // min/max_element is not available for MKL and CUDA
    if(Backend::get_preferred_backend() == PreferredBackend::generic)
    {
      TEST_CHECK_EQUAL(a.min_element(), DT_(-72));
      TEST_CHECK_EQUAL(a.max_element(), DT_(71));
    }

    TEST_CHECK_EQUAL(a.min_abs_element(), DT_(1));
    TEST_CHECK_EQUAL(a.max_abs_element(), DT_(72));
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMinMaxAbsElementTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRConvertTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRConvertTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRConvertTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    SparseMatrixFactory<DT_, IT_> b_fac(IT_(6), IT_(4));

    b_fac.add(IT_(0), IT_(0), DT_(2));
    b_fac.add(IT_(0), IT_(1), DT_(3));
    b_fac.add(IT_(1), IT_(0), DT_(4));
    b_fac.add(IT_(1), IT_(2), DT_(5));
    b_fac.add(IT_(2), IT_(1), DT_(6));
    b_fac.add(IT_(2), IT_(2), DT_(7));
    b_fac.add(IT_(5), IT_(2), DT_(8));
    b_fac.add(IT_(5), IT_(3), DT_(9));

    SparseMatrixCSR<DT_, IT_> csr(b_fac.make_csr());

    constexpr int height = 3;
    constexpr int width = 2;

    SparseMatrixBCSR<DT_, IT_, height, width> bcsr;
    bcsr.convert(csr);

    TEST_CHECK_EQUAL(bcsr.num_rows(), Index(2));
    TEST_CHECK_EQUAL(bcsr.num_cols(), Index(2));
    TEST_CHECK_EQUAL(bcsr.num_nzes(), Index(3));

    {
      Memory::TypedView<IT_> row_ptr = bcsr.row_ptr_view_r();
      TEST_CHECK_EQUAL(row_ptr(0), IT_(0));
      TEST_CHECK_EQUAL(row_ptr(1), IT_(2));
      TEST_CHECK_EQUAL(row_ptr(2), IT_(3));
    }

    {
      Memory::TypedView<IT_> col_idx = bcsr.col_idx_view_r();
      TEST_CHECK_EQUAL(col_idx(0), IT_(0));
      TEST_CHECK_EQUAL(col_idx(1), IT_(1));
      TEST_CHECK_EQUAL(col_idx(2), IT_(1));
    }

    {
      Memory::TypedView<DT_> val_raw = bcsr.val_view_raw_r();
      TEST_CHECK_RELATIVE(val_raw(0), DT_(2), tol);
      TEST_CHECK_RELATIVE(val_raw(1), DT_(3), tol);
      TEST_CHECK_RELATIVE(val_raw(2), DT_(4), tol);
      TEST_CHECK_EQUAL(val_raw(3), DT_(0));
      TEST_CHECK_EQUAL(val_raw(4), DT_(0));
      TEST_CHECK_RELATIVE(val_raw(5), DT_(6), tol);
      TEST_CHECK_EQUAL(val_raw(6), DT_(0));
      TEST_CHECK_EQUAL(val_raw(7), DT_(0));
      TEST_CHECK_RELATIVE(val_raw(8), DT_(5), tol);
      TEST_CHECK_EQUAL(val_raw(9), DT_(0));
      TEST_CHECK_RELATIVE(val_raw(10), DT_(7), tol);
      TEST_CHECK_EQUAL(val_raw(11), DT_(0));
      TEST_CHECK_EQUAL(val_raw(12), DT_(0));
      TEST_CHECK_EQUAL(val_raw(13), DT_(0));
      TEST_CHECK_EQUAL(val_raw(14), DT_(0));
      TEST_CHECK_EQUAL(val_raw(15), DT_(0));
      TEST_CHECK_RELATIVE(val_raw(16), DT_(8), tol);
      TEST_CHECK_RELATIVE(val_raw(17), DT_(9), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRConvertTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRConvertTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRConvertTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRConvertTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRConvertTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRConvertTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRMaxRelDiffTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(12); size < Index(50); size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size + 2);

      for (Index row(0); row < a_fac.num_rows(); ++row)
      {
        a_fac.add(row, row, DT_(0));
        a_fac.add(row, row+1, DT_(0));
        a_fac.add(row, row+2, DT_(0));
      }

      SparseMatrixCSR<DT_, IT_> a_csr(a_fac.make_csr());

      SparseMatrixBCSR<DT_, IT_, 2, 3> a, b;
      a.convert(a_csr);
      b.clone(a, LAFEM::CloneMode::Layout);

      const Index nnzes = a.num_nzes_raw();

      const Index off0 = (3*nnzes) / 8;
      const Index off1 = (1*nnzes) / 8;
      const Index off2 = (6*nnzes) / 8;

      // a = i, b = i
      {
        Memory::TypedView<DT_> va = a.val_view_raw_w();
        Memory::TypedView<DT_> vb = b.val_view_raw_w();
        for (Index i(0); i < nnzes; ++i)
        {
          va[i] = vb[i] = DT_(int(i) - int(off0)) * DT_(0.123);
        }
      }

      // identical vectors, result should be zero
      TEST_CHECK_LESS_THAN(a.max_rel_diff(b), tol*tol);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(a), tol*tol);

      // two values close to zero
      const DT_ delta_a0(Math::sqrt(Math::eps<DT_>()));
      const DT_ delta_b0(Math::sqr(Math::eps<DT_>()));
      const DT_ ref0 = (delta_a0 + delta_b0) / (DT_(1) + delta_a0 + delta_b0);
      a.val_view_raw_rw()[off0] += delta_a0;
      b.val_view_raw_rw()[off0] -= delta_b0;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref0, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref0, tol);

      const DT_ delta1 = DT_(0.17);
      const DT_ ref1 = delta1 / (DT_(off0 - off1)*DT_(0.246) + delta1 + DT_(1));
      a.val_view_raw_rw()[off1] -= delta1;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref1, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref1, tol);

      const DT_ delta2 = DT_(0.73);
      const DT_ ref2 = delta2 / (DT_(off2 - off0)*DT_(0.246) + delta2 + DT_(1));
      b.val_view_raw_rw()[off2] += delta2;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref2, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref2, tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRSameLayoutTest
  : public UnitTest
{
public:
  explicit SparseMatrixBCSRSameLayoutTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    for (Index size(4); size < Index(128); size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> fac_a(size, size+2);
      for (IT_ i = 0; i < IT_(size); ++i)
        fac_a.add(i, i, DT_(i));

      SparseMatrixCSR<DT_, IT_> a_csr(fac_a.make_csr());
      SparseMatrixBCSR<DT_, IT_, 2, 2> a;
      a.convert(a_csr);

      // weak copy
      auto b = a.clone(CloneMode::Weak);
      TEST_CHECK(a.same_layout(b));

      // shallow copy
      auto c = a.clone(CloneMode::Shallow);
      TEST_CHECK(a.same_layout(c));

      // different values at same position
      SparseMatrixFactory<DT_, IT_> fac_d(size, size+2);
      for (IT_ i = 0; i < IT_(size); ++i)
        fac_d.add(i, i , DT_(i + 1));
      SparseMatrixCSR<DT_, IT_> d_csr(fac_d.make_csr());
      SparseMatrixBCSR<DT_, IT_, 2, 2> d;
      d.convert(d_csr);
      TEST_CHECK(a.same_layout(d));

      // values at different position
      SparseMatrixFactory<DT_, IT_> fac_e(size, size+2);
      for (IT_ i = 0; i < IT_(size - 1); ++i)
        fac_e.add(i, i + 2, DT_(i));
      SparseMatrixCSR<DT_, IT_> e_csr(fac_e.make_csr());
      SparseMatrixBCSR<DT_, IT_, 2, 2> e;
      e.convert(e_csr);
      TEST_CHECK(!a.same_layout(e));

      // different sizes
      SparseMatrixFactory<DT_, IT_> fac_f(size + 2, size + 2);
      for (IT_ i = 0; i < IT_(size + 2); ++i)
        fac_f.add(i, i, DT_(i));
      SparseMatrixCSR<DT_, IT_> f_csr(fac_f.make_csr());
      SparseMatrixBCSR<DT_, IT_, 2, 2> f;
      f.convert(f_csr);
      TEST_CHECK(!a.same_layout(f));
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBCSRSameLayoutTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
