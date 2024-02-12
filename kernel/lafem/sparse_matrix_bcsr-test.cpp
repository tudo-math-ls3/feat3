// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>

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
   SparseMatrixBCSRTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRTest()
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
    test_vector_types();

    SparseMatrixBCSR<DT_, IT_, 2, 3> zero1;
    SparseMatrixBCSR<DT_, IT_, 2, 3> zero2;
    TEST_CHECK_EQUAL(zero2, zero1);
    zero2.convert(zero1);

    DenseVector<DT_, IT_> dv1(12);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i+1)/DT_(7*i+1));
    DenseVector<IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(2, 2, dv2, dv1, dv3);

    TEST_CHECK_EQUAL(c(1,0)(0,0), DT_(0));
    TEST_CHECK_EQUAL(c(1,1)(1,1), DT_(10+1)/DT_(7*10+1));

    SparseMatrixBCSR<DT_, IT_, 2, 3> d;
    d.convert(c);
    TEST_CHECK_EQUAL(d.rows(), c.rows());
    TEST_CHECK_EQUAL(d.columns(), c.columns());
    TEST_CHECK_EQUAL(d.used_elements(), c.used_elements());
    TEST_CHECK_EQUAL(d, c);
    TEST_CHECK_EQUAL((void*)d.template val<Perspective::pod>(), (void*)c.template val<Perspective::pod>());
    TEST_CHECK_EQUAL((void*)d.row_ptr(), (void*)c.row_ptr());
    SparseMatrixBCSR<DT_, IT_, 2, 3> e;
    e.clone(c);
    TEST_CHECK_EQUAL(e, c);
    TEST_CHECK_NOT_EQUAL((void*)e.template val<Perspective::pod>(), (void*)c.template val<Perspective::pod>());
    TEST_CHECK_EQUAL((void*)e.row_ptr(), (void*)c.row_ptr());
    e = c.clone(CloneMode::Deep);
    TEST_CHECK_EQUAL(e, c);
    TEST_CHECK_NOT_EQUAL((void*)e.template val<Perspective::pod>(), (void*)c.template val<Perspective::pod>());
    TEST_CHECK_NOT_EQUAL((void*)e.row_ptr(), (void*)c.row_ptr());

    SparseMatrixBCSR<DT_, IT_, 2, 3> f(c.layout());
    TEST_CHECK_EQUAL(f.rows(), c.rows());
    TEST_CHECK_EQUAL(f.columns(), c.columns());
    TEST_CHECK_EQUAL(f.used_elements(), c.used_elements());
    TEST_CHECK_NOT_EQUAL((void*)f.template val<Perspective::pod>(), (void*)c.template val<Perspective::pod>());
    TEST_CHECK_EQUAL((void*)f.row_ptr(), (void*)c.row_ptr());

  }
};
SparseMatrixBCSRTest <float, std::uint64_t> cpu_sparse_matrix_bcsr_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRTest <double, std::uint64_t> cpu_sparse_matrix_bcsr_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRTest <float, std::uint32_t> cpu_sparse_matrix_bcsr_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRTest <double, std::uint32_t> cpu_sparse_matrix_bcsr_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRTest <float, std::uint64_t> mkl_cpu_sparse_matrix_bcsr_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRTest <double, std::uint64_t> mkl_cpu_sparse_matrix_bcsr_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRTest <__float128, std::uint64_t> cpu_sparse_matrix_bcsr_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRTest <__float128, std::uint32_t> cpu_sparse_matrix_bcsr_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRTest <Half, std::uint32_t> cpu_sparse_matrix_bcsr_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRTest <Half, std::uint64_t> cpu_sparse_matrix_bcsr_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRTest <float, std::uint64_t> cuda_sparse_matrix_bcsr_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBCSRTest <double, std::uint64_t> cuda_sparse_matrix_bcsr_test_double_uint64(PreferredBackend::cuda);
SparseMatrixBCSRTest <float, std::uint32_t> cuda_sparse_matrix_bcsr_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBCSRTest <double, std::uint32_t> cuda_sparse_matrix_bcsr_test_double_uint32(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRSerializeTest
  : public UnitTest
{
public:
   SparseMatrixBCSRSerializeTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRSerializeTest()
  {
  }
  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv1(12);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i+1)/DT_(7*i+1));
    DenseVector<IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(2, 2, dv2, dv1, dv3);

    BinaryStream bs;
    c.write_out(FileMode::fm_bcsr, bs);
    bs.seekg(0);
    SparseMatrixBCSR<DT_, IT_, 2, 3> g(FileMode::fm_bcsr, bs);
    TEST_CHECK_EQUAL(g, c);

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
    TEST_CHECK_EQUAL(k, c);

#ifdef FEAT_HAVE_ZLIB
    auto zl = c.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixBCSR<DT_, IT_, 2, 3> zlib(zl);
    TEST_CHECK_EQUAL(k, c);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = c.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixBCSR<DT_, IT_, 2, 3> zfp(zf);
    for(Index i(0) ; i < c.rows() ; ++i)
    {
      for(Index j(0) ; j < c.columns() ; ++j)
      {
        for(int q(0) ; q < c(i,j).m ; ++q)
        {
          for(int w(0) ; w < c(i,j).n ; ++w)
          {
            TEST_CHECK_EQUAL_WITHIN_EPS(c(i,j)(q,w), zfp(i,j)(q,w), DT_(1e-5));
          }
        }
      }
    }
#endif
  }
};

SparseMatrixBCSRSerializeTest <float, std::uint64_t> cpu_sparse_matrix_bcsr_serialize_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRSerializeTest <double, std::uint64_t> cpu_sparse_matrix_bcsr_serialize_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRSerializeTest <float, std::uint32_t> cpu_sparse_matrix_bcsr_serialize_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRSerializeTest <double, std::uint32_t> cpu_sparse_matrix_bcsr_serialize_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRSerializeTest <float, std::uint64_t> mkl_cpu_sparse_matrix_bcsr_serialize_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRSerializeTest <double, std::uint64_t> mkl_cpu_sparse_matrix_bcsr_serialize_test_double_uint64(PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SparseMatrixBCSRSerializeTest <__float128, std::uint64_t> cpu_sparse_matrix_bcsr_serialize_test_float128_uint64(PreferredBackend::generic);
//SparseMatrixBCSRSerializeTest <__float128, std::uint32_t> cpu_sparse_matrix_bcsr_serialize_test_float128_uint32(PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRSerializeTest <Half, std::uint32_t> cpu_sparse_matrix_bcsr_serialize_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRSerializeTest <Half, std::uint64_t> cpu_sparse_matrix_bcsr_serialize_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRSerializeTest <float, std::uint64_t> cuda_sparse_matrix_bcsr_serialize_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBCSRSerializeTest <double, std::uint64_t> cuda_sparse_matrix_bcsr_serialize_test_double_uint64(PreferredBackend::cuda);
SparseMatrixBCSRSerializeTest <float, std::uint32_t> cuda_sparse_matrix_bcsr_serialize_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBCSRSerializeTest <double, std::uint32_t> cuda_sparse_matrix_bcsr_serialize_test_double_uint32(PreferredBackend::cuda);
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
   SparseMatrixBCSRApplyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRApplyTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv1(12);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i+1));
    DenseVector<IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(2, 2, dv2, dv1, dv3);

    DenseVector<DT_, IT_> x(c.template columns<Perspective::pod>());
    DenseVector<DT_, IT_> y(c.template rows<Perspective::pod>());
    DenseVector<DT_, IT_> r(c.template rows<Perspective::pod>());
    DenseVector<DT_, IT_> ref(c.template rows<Perspective::pod>());
    for (Index i(0) ; i < x.size() ; ++i)
    {
      x(i, DT_(i));
    }
    for (Index i(0) ; i < r.size() ; ++i)
    {
      r(i, DT_(4711));
      ref(i, DT_(4711));
      y(i, DT_(i % 100));
    }
    DenseVectorBlocked<DT_, IT_, 3> xb(x);
    DenseVectorBlocked<DT_, IT_, 2> yb(y);
    DenseVectorBlocked<DT_, IT_, 2> rb(r);

    SparseMatrixCSR<DT_, IT_> csr;
    csr.convert(c);
    csr.apply(ref, x);

    c.apply(r, x);

    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, x);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(r, xb);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, xb);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    //defect
    DT_ alpha(-1);
    csr.apply(ref, x, y, alpha);
    c.apply(r, x, y, alpha);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, x, yb, alpha);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(r, xb, y, alpha);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, xb, yb, alpha);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, xb, y, alpha);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    //axpy
    alpha = DT_(1.234);
    csr.apply(ref, x, y, alpha);
    c.apply(r, x, y, alpha);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    // &r == &y
    r.copy(y);
    c.apply(r, x, r, alpha);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    c.apply(rb, x, yb, alpha);
    r.convert(rb);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    c.apply(r, xb, y, alpha);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    c.apply(rb, xb, yb, alpha);
    r.convert(rb);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    c.apply(rb, xb, y, alpha);
    r.convert(rb);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));
  }
};
SparseMatrixBCSRApplyTest <float, std::uint64_t> cpu_sm_bcsr_apply_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRApplyTest <double, std::uint64_t> cpu_sm_bcsr_apply_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRApplyTest <float, std::uint32_t> cpu_sm_bcsr_apply_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRApplyTest <double, std::uint32_t> cpu_sm_bcsr_apply_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRApplyTest <float, std::uint64_t> mkl_cpu_sm_bcsr_apply_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRApplyTest <double, std::uint64_t> mkl_cpu_sm_bcsr_apply_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRApplyTest <__float128, std::uint64_t> cpu_sm_bcsr_apply_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRApplyTest <__float128, std::uint32_t> cpu_sm_bcsr_apply_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRApplyTest <Half, std::uint32_t> cpu_sm_bcsr_apply_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRApplyTest <Half, std::uint64_t> cpu_sm_bcsr_apply_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRApplyTest <float, std::uint64_t> cuda_sm_bcsr_apply_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBCSRApplyTest <double, std::uint64_t> cuda_sm_bcsr_apply_test_double_uint64(PreferredBackend::cuda);
SparseMatrixBCSRApplyTest <float, std::uint32_t> cuda_sm_bcsr_apply_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBCSRApplyTest <double, std::uint32_t> cuda_sm_bcsr_apply_test_double_uint32(PreferredBackend::cuda);
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
   SparseMatrixBCSRApplySquareTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRApplySquareTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRApplySquareTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv1(18);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i+1));
    DenseVector<IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 3, 3> c(2, 2, dv2, dv1, dv3);

    DenseVector<DT_, IT_> x(c.template columns<Perspective::pod>());
    DenseVector<DT_, IT_> y(c.template rows<Perspective::pod>());
    DenseVector<DT_, IT_> r(c.template rows<Perspective::pod>());
    DenseVector<DT_, IT_> ref(c.template rows<Perspective::pod>());
    for (Index i(0) ; i < x.size() ; ++i)
    {
      x(i, DT_(i));
    }
    for (Index i(0) ; i < r.size() ; ++i)
    {
      r(i, DT_(4711));
      ref(i, DT_(4711));
      y(i, DT_(i % 100));
    }
    DenseVectorBlocked<DT_, IT_, 3> xb(x);
    DenseVectorBlocked<DT_, IT_, 3> yb(y);
    DenseVectorBlocked<DT_, IT_, 3> rb(r);

    SparseMatrixCSR<DT_, IT_> csr;
    csr.convert(c);
    csr.apply(ref, x);

    c.apply(r, x);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, x);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(r, xb);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, xb);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    // defect
    DT_ alpha(-1);
    csr.apply(ref, x, y, alpha);
    c.apply(r, x, y, alpha);
    TEST_CHECK_EQUAL(r, ref);

    // &r == &y
    r.copy(y);
    c.apply(r, x, r, alpha);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, x, yb, alpha);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(r, xb, y, alpha);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, xb, yb, alpha);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    c.apply(rb, xb, y, alpha);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    // axpy
    alpha = DT_(1.234);
    csr.apply(ref, x, y, alpha);
    c.apply(r, x, y, alpha);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    // &r == &y
    r.copy(y);
    c.apply(r, x, y, alpha);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    c.apply(rb, x, yb, alpha);
    r.convert(rb);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    c.apply(r, xb, y, alpha);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    c.apply(rb, xb, yb, alpha);
    r.convert(rb);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));

    c.apply(rb, xb, y, alpha);
    r.convert(rb);
    for (Index i(0) ; i < r.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-3));
  }
};
SparseMatrixBCSRApplySquareTest <float, std::uint64_t> cpu_sm_bcsr_apply_square_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRApplySquareTest <double, std::uint64_t> cpu_sm_bcsr_apply_square_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRApplySquareTest <float, std::uint32_t> cpu_sm_bcsr_apply_square_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRApplySquareTest <double, std::uint32_t> cpu_sm_bcsr_apply_square_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRApplySquareTest <float, std::uint64_t> mkl_cpu_sm_bcsr_apply_square_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRApplySquareTest <double, std::uint64_t> mkl_cpu_sm_bcsr_apply_square_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRApplySquareTest <__float128, std::uint64_t> cpu_sm_bcsr_apply_square_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRApplySquareTest <__float128, std::uint32_t> cpu_sm_bcsr_apply_square_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRApplySquareTest <Half, std::uint32_t> cpu_sm_bcsr_apply_square_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRApplySquareTest <Half, std::uint64_t> cpu_sm_bcsr_apply_square_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRApplySquareTest <float, std::uint32_t> gpu_sm_bcsr_apply_square_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBCSRApplySquareTest <double, std::uint32_t> gpu_sm_bcsr_apply_square_test_double_uint32(PreferredBackend::cuda);
SparseMatrixBCSRApplySquareTest <float, std::uint64_t> gpu_sm_bcsr_apply_square_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBCSRApplySquareTest <double, std::uint64_t> gpu_sm_bcsr_apply_square_test_double_uint64(PreferredBackend::cuda);
#endif


template<
  typename DT_,
  typename IT_>
class SparseMatrixBCSRDiagTest
  : public UnitTest
{
public:
   SparseMatrixBCSRDiagTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRDiagTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRDiagTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv1(18);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i+1));
    DenseVector<IT_, IT_> dv2(3);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    dv2(2, IT_(2));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 3, 3> smb(2, 2, dv2, dv1, dv3);

    auto diag = smb.extract_diag();
    TEST_CHECK_EQUAL(diag.template elements<Perspective::pod>()[0], 1);
    TEST_CHECK_EQUAL(diag.template elements<Perspective::pod>()[1], 5);
    TEST_CHECK_EQUAL(diag.template elements<Perspective::pod>()[2], 9);
    TEST_CHECK_EQUAL(diag.template elements<Perspective::pod>()[3], 10);
    TEST_CHECK_EQUAL(diag.template elements<Perspective::pod>()[4], 14);
    TEST_CHECK_EQUAL(diag.template elements<Perspective::pod>()[5], 18);
  }
};
SparseMatrixBCSRDiagTest <float, std::uint64_t> cpu_sm_bcsr_diag_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRDiagTest <double, std::uint64_t> cpu_sm_bcsr_diag_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRDiagTest <float, std::uint32_t> cpu_sm_bcsr_diag_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRDiagTest <double, std::uint32_t> cpu_sm_bcsr_diag_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRDiagTest <float, std::uint64_t> mkl_cpu_sm_bcsr_diag_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRDiagTest <double, std::uint64_t> mkl_cpu_sm_bcsr_diag_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRDiagTest <__float128, std::uint64_t> cpu_sm_bcsr_diag_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRDiagTest <__float128, std::uint32_t> cpu_sm_bcsr_diag_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRDiagTest <Half, std::uint32_t> cpu_sm_bcsr_diag_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRDiagTest <Half, std::uint64_t> cpu_sm_bcsr_diag_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRDiagTest <float, std::uint64_t> cuda_sm_bcsr_diag_test_float_uint64_cuda(PreferredBackend::cuda);
SparseMatrixBCSRDiagTest <double, std::uint64_t> cuda_sm_bcsr_diag_test_double_uint64_cuda(PreferredBackend::cuda);
SparseMatrixBCSRDiagTest <float, std::uint32_t> cuda_sm_bcsr_diag_test_float_uint32_cuda(PreferredBackend::cuda);
SparseMatrixBCSRDiagTest <double, std::uint32_t> cuda_sm_bcsr_diag_test_double_uint32_cuda(PreferredBackend::cuda);
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
   SparseMatrixBCSRScaleTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRScaleTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv1(12);
    for (Index i(0) ; i < dv1.size() ; ++i)
    {
      dv1(i, DT_(i+1));
    }
    DenseVector<IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 2, 3> a(2, 2, dv2, dv1, dv3);
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(a.layout());

    DT_ scal = DT_(4711);

    SparseMatrixCSR<DT_, IT_> a_csr;
    a_csr.convert(a);
    SparseMatrixCSR<DT_, IT_> ref(a_csr.layout());
    SparseMatrixCSR<DT_, IT_> result_csr;
    ref.scale(a_csr, scal);

    c.scale(a, scal);
    result_csr.convert(c);
    TEST_CHECK_EQUAL(result_csr, ref);

    a.scale(a, scal);
    result_csr.convert(a);
    TEST_CHECK_EQUAL(result_csr, ref);
  }
};
SparseMatrixBCSRScaleTest <float, std::uint64_t> cpu_sm_bcsr_scale_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRScaleTest <double, std::uint64_t> cpu_sm_bcsr_scale_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRScaleTest <float, std::uint32_t> cpu_sm_bcsr_scale_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRScaleTest <double, std::uint32_t> cpu_sm_bcsr_scale_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRScaleTest <float, std::uint64_t> mkl_cpu_sm_bcsr_scale_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRScaleTest <double, std::uint64_t> mkl_cpu_sm_bcsr_scale_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRScaleTest <__float128, std::uint64_t> cpu_sparse_matrix_bcsr_scale_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRScaleTest <__float128, std::uint32_t> cpu_sparse_matrix_bcsr_scale_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRScaleTest <Half, std::uint32_t> cpu_sm_bcsr_scale_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRScaleTest <Half, std::uint64_t> cpu_sm_bcsr_scale_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRScaleTest <float, std::uint64_t> cuda_sm_bcsr_scale_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBCSRScaleTest <double, std::uint64_t> cuda_sm_bcsr_scale_test_double_uint64(PreferredBackend::cuda);
SparseMatrixBCSRScaleTest <float, std::uint32_t> cuda_sm_bcsr_scale_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBCSRScaleTest <double, std::uint32_t> cuda_sm_bcsr_scale_test_double_uint32(PreferredBackend::cuda);
#endif


template<
  typename DT_,
  typename IT_>
  class SparseMatrixBCSRScaleRowColTest
  : public UnitTest
{
public:
  SparseMatrixBCSRScaleRowColTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRScaleRowColTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRScaleRowColTest()
  {
  }

  virtual void run() const override
  {

    const DT_ pi(Math::pi<DT_>());
    const DT_ eps(Math::pow(Math::eps<DT_>(), DT_(0.8)));

    DenseVector<DT_, IT_> dv1(30);
    for (Index i(0) ; i < dv1.size() ; ++i)
    {
      dv1(i, DT_(i+1));
    }
    DenseVector<IT_, IT_> dv2(5);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    dv2(2, IT_(0));
    dv2(3, IT_(2));
    dv2(4, IT_(3));

    DenseVector<IT_, IT_> dv3(4);
    dv3(0, IT_(0));
    dv3(1, IT_(2));
    dv3(2, IT_(5));
    dv3(3, IT_(5));

    SparseMatrixBCSR<DT_, IT_, 2, 3> a(3, 4, dv2, dv1, dv3);
    SparseMatrixBCSR<DT_, IT_, 2, 3> b;
    b.clone(a);

    // Scale rows

    DenseVectorBlocked<DT_, IT_, 2> s(a.rows());
    for (Index i(0); i < s.size(); ++i)
    {
      s(i, Tiny::Vector<DT_, 2>{pi * (i % 3 + 1) - DT_(5.21) + DT_(i),  DT_(i) * DT_(12.31) });
    }
    DenseVectorBlocked<DT_, IT_, 3> t(a.columns());
    for (Index i(0); i < t.size(); ++i)
    {
      t(i, Tiny::Vector<DT_, 3>{pi * (i % 3 + 1) - DT_(5.21) + DT_(i),  DT_(i) * DT_(12.31), DT_(i) *DT_(-324.21) - DT_(13.37) });
    }

    b.scale_rows(b, s);
    for (Index row(0); row < a.rows(); ++row)
    {
      for (Index col(0); col < a.columns(); ++col)
      {
        auto ab = a(row, col);
        auto bb = b(row, col);
        auto sb = s(row);
        for (int irow(0); irow < 2; ++irow)
        {
          for (int icol(0); icol < 3; ++icol)
          {
            auto bv = bb(irow, icol);
            auto av = ab(irow, icol);
            TEST_CHECK_EQUAL_WITHIN_EPS(bv, av * sb(irow), eps);
          }
        }
      }
    }

    a.clone(b);
    b.scale_cols(b, t);
    for (Index row(0); row < a.rows(); ++row)
    {
      for (Index col(0); col < a.columns(); ++col)
      {
        auto ab = a(row, col);
        auto bb = b(row, col);
        auto tb = t(col);
        for (int irow(0); irow < 2; ++irow)
        {
          for (int icol(0); icol < 3; ++icol)
          {
            auto bv = bb(irow, icol);
            auto av = ab(irow, icol);
            TEST_CHECK_EQUAL_WITHIN_EPS(bv, av * tb(icol), eps);
          }
        }
      }
    }
  }
};

SparseMatrixBCSRScaleRowColTest <float, std::uint64_t>  cpu_sm_bcsr_row_col_scale_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRScaleRowColTest <double, std::uint64_t> cpu_sm_bcsr_row_col_scale_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRScaleRowColTest <float, std::uint32_t>  cpu_sm_bcsr_row_col_scale_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRScaleRowColTest <double, std::uint32_t> cpu_sm_bcsr_row_col_scale_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRScaleRowColTest <float, std::uint64_t>  mkl_cpu_sm_bcsr_row_col_scale_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRScaleRowColTest <double, std::uint64_t> mkl_cpu_sm_bcsr_row_col_scale_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRScaleRowColTest <__float128, std::uint64_t> cpu_sparse_matrix_bcsr_row_col_scale_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRScaleRowColTest <__float128, std::uint32_t> cpu_sparse_matrix_bcsr_row_col_scale_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRScaleRowColTest <Half, std::uint32_t> cpu_sm_bcsr_row_col_scale_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRScaleRowColTest <Half, std::uint64_t> cpu_sm_bcsr_row_col_scale_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRScaleRowColTest <float, std::uint64_t>  cuda_sm_bcsr_row_col_scale_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBCSRScaleRowColTest <double, std::uint64_t> cuda_sm_bcsr_row_col_scale_test_double_uint64(PreferredBackend::cuda);
SparseMatrixBCSRScaleRowColTest <float, std::uint32_t>  cuda_sm_bcsr_row_col_scale_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBCSRScaleRowColTest <double, std::uint32_t> cuda_sm_bcsr_row_col_scale_test_double_uint32(PreferredBackend::cuda);
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
   SparseMatrixBCSRNormTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRNormTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRNormTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv1(12);
    for (Index i(0) ; i < dv1.size() ; ++i)
    {
      dv1(i, DT_(i+1));
    }
    DenseVector<IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 2, 3> a(2, 2, dv2, dv1, dv3);

    SparseMatrixCSR<DT_, IT_> a_csr;
    a_csr.convert(a);
    DT_ ref = a_csr.norm_frobenius();

    DT_ result = a.norm_frobenius();
    TEST_CHECK_EQUAL(result, ref);
  }
};
SparseMatrixBCSRNormTest <float, std::uint64_t> cpu_sm_bcsr_norm_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRNormTest <double, std::uint64_t> cpu_sm_bcsr_norm_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRNormTest <float, std::uint32_t> cpu_sm_bcsr_norm_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRNormTest <double, std::uint32_t> cpu_sm_bcsr_norm_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRNormTest <float, std::uint64_t> mkl_cpu_sm_bcsr_norm_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRNormTest <double, std::uint64_t> mkl_cpu_sm_bcsr_norm_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRNormTest <__float128, std::uint64_t> cpu_sm_bcsr_norm_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRNormTest <__float128, std::uint32_t> cpu_sm_bcsr_norm_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRNormTest <Half, std::uint32_t> cpu_sm_bcsr_norm_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRNormTest <Half, std::uint64_t> cpu_sm_bcsr_norm_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRNormTest <float, std::uint64_t> cuda_sm_bcsr_norm_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBCSRNormTest <double, std::uint64_t> cuda_sm_bcsr_norm_test_double_uint64(PreferredBackend::cuda);
SparseMatrixBCSRNormTest <float, std::uint32_t> cuda_sm_bcsr_norm_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBCSRNormTest <double, std::uint32_t> cuda_sm_bcsr_norm_test_double_uint32(PreferredBackend::cuda);
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
   SparseMatrixBCSRAxpyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRAxpyTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv1(12);
    for (Index i(0) ; i < dv1.size() ; ++i)
    {
      dv1(i, DT_(i+1));
    }
    DenseVector<IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 2, 3> a(2, 2, dv2, dv1, dv3);
    DenseVector<DT_, IT_> dv4(12);
    for (Index i(0) ; i < dv4.size() ; ++i)
    {
      dv4(i, DT_(i-1));
    }
    SparseMatrixBCSR<DT_, IT_, 2, 3> b(2, 2, dv2, dv4, dv3);
    SparseMatrixBCSR<DT_, IT_, 2, 3> c(a.layout());
    SparseMatrixBCSR<DT_, IT_, 2, 3> ref(a.layout());

    DT_ scal = DT_(1.234);

    for(Index i(0) ; i < ref.template used_elements<Perspective::pod>() ; ++i)
    {
      ref.template val<Perspective::pod>()[i] = scal * dv1(i) + dv4(i);
    }

    c.axpy(a, b, scal);
    TEST_CHECK_EQUAL(c, ref);

    c.copy(b);
    c.axpy(a, c, scal);
    TEST_CHECK_EQUAL(c, ref);

    c.copy(a);
    c.axpy(c, b, scal);
    TEST_CHECK_EQUAL(c, ref);

    scal = DT_(0);
    ref.clone(b);
    c.axpy(a, b, scal);
    TEST_CHECK_EQUAL(c, ref);

    scal = DT_(1);
    for(Index i(0) ; i < ref.template used_elements<Perspective::pod>() ; ++i)
    {
      ref.template val<Perspective::pod>()[i] = dv1(i) + dv4(i);
    }
    c.axpy(a, b, scal);
    TEST_CHECK_EQUAL(c, ref);

    scal = DT_(-1);
    for(Index i(0) ; i < ref.template used_elements<Perspective::pod>() ; ++i)
    {
      ref.template val<Perspective::pod>()[i] = dv4(i) - dv1(i);
    }
    c.axpy(a, b, scal);
    TEST_CHECK_EQUAL(c, ref);
  }
};
SparseMatrixBCSRAxpyTest <float, std::uint64_t> cpu_sm_bcsr_axpy_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRAxpyTest <double, std::uint64_t> cpu_sm_bcsr_axpy_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRAxpyTest <float, std::uint32_t> cpu_sm_bcsr_axpy_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRAxpyTest <double, std::uint32_t> cpu_sm_bcsr_axpy_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRAxpyTest <float, std::uint64_t> mkl_cpu_sm_bcsr_axpy_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRAxpyTest <double, std::uint64_t> mkl_cpu_sm_bcsr_axpy_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRAxpyTest <__float128, std::uint64_t> cpu_sm_bcsr_axpy_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRAxpyTest <__float128, std::uint32_t> cpu_sm_bcsr_axpy_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRAxpyTest <Half, std::uint32_t> cpu_sm_bcsr_axpy_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRAxpyTest <Half, std::uint64_t> cpu_sm_bcsr_axpy_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRAxpyTest <float, std::uint64_t> cuda_sm_bcsr_axpy_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBCSRAxpyTest <double, std::uint64_t> cuda_sm_bcsr_axpy_test_double_uint64(PreferredBackend::cuda);
SparseMatrixBCSRAxpyTest <float, std::uint32_t> cuda_sm_bcsr_axpy_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBCSRAxpyTest <double, std::uint32_t> cuda_sm_bcsr_axpy_test_double_uint32(PreferredBackend::cuda);
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
   SparseMatrixBCSRPermuteTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBCSRPermuteTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBCSRPermuteTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv1(18);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i+1));
    DenseVector<IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixBCSR<DT_, IT_, 3, 3> a(2, 2, dv2, dv1, dv3);

    DenseVector<DT_, IT_> x(a.template columns<Perspective::pod>());
    DenseVector<DT_, IT_> r(a.template rows<Perspective::pod>());
    for (Index i(0) ; i < x.size() ; ++i)
    {
      x(i, DT_(i));
    }
    DenseVectorBlocked<DT_, IT_, 3> xb(x);
    DenseVectorBlocked<DT_, IT_, 3> rb(r);

    /////////////////////////

    std::cout<<xb<<std::endl;
    a.apply(rb, xb);
    DT_ ref_norm = rb.norm2();

    auto a_backup = a.clone(CloneMode::Deep);

    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rng(seed);
    Adjacency::Permutation perm(a.rows(), rng);

    a.permute(perm, perm);
    xb.permute(perm);

    std::cout<<xb<<std::endl;
    a.apply(rb, xb);
    DT_ norm = r.norm2();
    TEST_CHECK_EQUAL_WITHIN_EPS(norm, ref_norm, DT_(1e-4));

    a = a_backup.clone(CloneMode::Deep);
    auto perm_inv = perm.inverse();
    a.permute(perm_inv, perm);
    a.permute(perm, perm_inv);
    TEST_CHECK_EQUAL(a, a_backup);
  }
};
SparseMatrixBCSRPermuteTest <float, std::uint64_t> cpu_sm_bcsr_permute_test_float_uint64(PreferredBackend::generic);
SparseMatrixBCSRPermuteTest <double, std::uint64_t> cpu_sm_bcsr_permute_test_double_uint64(PreferredBackend::generic);
SparseMatrixBCSRPermuteTest <float, std::uint32_t> cpu_sm_bcsr_permute_test_float_uint32(PreferredBackend::generic);
SparseMatrixBCSRPermuteTest <double, std::uint32_t> cpu_sm_bcsr_permute_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBCSRPermuteTest <float, std::uint64_t> mkl_cpu_sm_bcsr_permute_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBCSRPermuteTest <double, std::uint64_t> mkl_cpu_sm_bcsr_permute_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBCSRPermuteTest <__float128, std::uint64_t> cpu_sm_bcsr_permute_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBCSRPermuteTest <__float128, std::uint32_t> cpu_sm_bcsr_permute_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBCSRPermuteTest <Half, std::uint32_t> cpu_sm_bcsr_permute_test_half_uint32(PreferredBackend::generic);
SparseMatrixBCSRPermuteTest <Half, std::uint64_t> cpu_sm_bcsr_permute_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBCSRPermuteTest <float, std::uint64_t> cuda_sm_bcsr_permute_test_float_uint64_cuda(PreferredBackend::cuda);
SparseMatrixBCSRPermuteTest <double, std::uint64_t> cuda_sm_bcsr_permute_test_double_uint64_cuda(PreferredBackend::cuda);
SparseMatrixBCSRPermuteTest <float, std::uint32_t> cuda_sm_bcsr_permute_test_float_uint32_cuda(PreferredBackend::cuda);
SparseMatrixBCSRPermuteTest <double, std::uint32_t> cuda_sm_bcsr_permute_test_double_uint32_cuda(PreferredBackend::cuda);
#endif
