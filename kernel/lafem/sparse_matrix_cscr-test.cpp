// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_cscr.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>
//#include <kernel/adjacency/cuthill_mckee.hpp>

#include <sstream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse matrix csr class.
 *
 * \test test description missing
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRTest
  : public UnitTest
{
public:
   explicit SparseMatrixCSCRTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    SparseMatrixCSCR<DT_, IT_> zero;
    TEST_CHECK(zero.empty());

    SparseMatrixCSCR<DT_, IT_> zero3(10, 11, 3, 12);
    TEST_CHECK_EQUAL(zero3.num_nzes(), 12);
    TEST_CHECK_EQUAL(zero3.num_rows(), 10);
    TEST_CHECK_EQUAL(zero3.num_cols(), 11);
    TEST_CHECK_EQUAL(zero3.num_nonzero_rows(), 3);


    SparseMatrixCSCR<DT_, IT_> empty1(11, 12, 3, 111);
    SparseMatrixCSCR<DT_, IT_> empty2;
    SparseMatrixCSCR<DT_, IT_> empty3;
    empty2.convert(empty1);
    empty3.convert(empty2);
    TEST_CHECK_EQUAL(empty1.num_rows(), empty3.num_rows());
    TEST_CHECK_EQUAL(empty1.num_cols(), empty3.num_cols());
    TEST_CHECK_EQUAL(empty1.num_nzes(), empty3.num_nzes());
    TEST_CHECK_EQUAL(empty1.num_nonzero_rows(), empty3.num_nonzero_rows());

    SparseMatrixCSCR<DT_, IT_> empty4(empty2.layout());
    SparseMatrixCSCR<DT_, IT_> empty5(empty3.layout());
    empty4.convert(empty1);
    empty5.convert(empty4);
    TEST_CHECK_EQUAL(empty5.num_rows(), empty5.num_rows());
    TEST_CHECK_EQUAL(empty5.num_cols(), empty5.num_cols());
    TEST_CHECK_EQUAL(empty5.num_nzes(), empty5.num_nzes());
    TEST_CHECK_EQUAL(empty5.num_nonzero_rows(), empty5.num_nonzero_rows());
    empty5.convert(zero);
    TEST_CHECK_EQUAL(empty5.num_rows(), 0);
    TEST_CHECK_EQUAL(empty5.num_cols(), 0);
    TEST_CHECK_EQUAL(empty5.num_nzes(), 0);
    TEST_CHECK_EQUAL(empty5.num_nonzero_rows(), 0);

    DenseVector<DT_, IT_> val(6);
    DenseVector<IT_, IT_> col_idx(6);
    DenseVector<IT_, IT_> row_ptr(4);
    DenseVector<IT_, IT_> row_idx(3);
    {
      Memory::TypedView<DT_> vval = val.elements_view_w();
      Memory::TypedView<IT_> vci = col_idx.elements_view_w();
      Memory::TypedView<IT_> vrp = row_ptr.elements_view_w();
      Memory::TypedView<IT_> vri = row_idx.elements_view_w();
      for (Index i(0) ; i < val.size() ; ++i)
      {
        vval[i] = DT_(i+1);
        vci[i] = IT_(i + 3);
      }
      vrp[0] = IT_(0);
      vrp[1] = IT_(2);
      vrp[2] = IT_(5);
      vrp[3] = IT_(6);
      vri[0] = IT_(0);
      vri[1] = IT_(4);
      vri[2] = IT_(8);
    }

    SparseMatrixCSCR<DT_, IT_> a(10, 10, row_ptr, row_idx, col_idx, val);
    TEST_CHECK(!a.empty());
    TEST_CHECK_EQUAL(a.num_nzes(), 6);
    TEST_CHECK_EQUAL(a.num_rows(), 10);
    TEST_CHECK_EQUAL(a.num_cols(), 10);
    TEST_CHECK_EQUAL(a.num_nonzero_rows(), 3);
    //TEST_CHECK_EQUAL(a(0,0), DT_(0));
    //TEST_CHECK_EQUAL(a(0,3), DT_(1));
    //TEST_CHECK_EQUAL(a(4,4), DT_(0));
    //TEST_CHECK_EQUAL(a(4,6), DT_(4));
    //TEST_CHECK_EQUAL(a(9,9), DT_(0));

    SparseMatrixCSCR<DT_, IT_> b;
    b.convert(a);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(a), tol);
    SparseMatrixCSCR<DT_, IT_> c;
    c.convert(b);
    TEST_CHECK_LESS_THAN(c.max_rel_diff(b), tol);
    SparseMatrixCSCR<DT_, IT_> d(a.clone());
    TEST_CHECK_LESS_THAN(d.max_rel_diff(c), tol);


    SparseMatrixFactory<DT_, IT_> fac(IT_(10), IT_(10));
    fac.add(0, 3, DT_(1));
    fac.add(0, 4, DT_(2));
    fac.add(1, 5, DT_(7));
    fac.add(3, 3, DT_(8));
    fac.add(4, 5, DT_(3));
    fac.add(4, 6, DT_(4));
    fac.add(4, 7, DT_(5));
    fac.add(6, 6, DT_(9));
    fac.add(8, 8, DT_(6));
    SparseMatrixCSR<DT_, IT_> csr(10, 10);
    csr.convert(fac.make_csr());
    VectorMirror<DT_, IT_> mirror_main(10, 3);
    {
      Memory::TypedView<IT_> vmi = mirror_main.indices_view_w();
      vmi[0] = IT_(0);
      vmi[1] = IT_(4);
      vmi[2] = IT_(8);
    }
    VectorMirror<DT_, IT_> mirror;
    mirror.convert(mirror_main);
    SparseMatrixCSCR<DT_, IT_> f(csr, mirror);
    TEST_CHECK(f.same_layout(a));
    TEST_CHECK_LESS_THAN(f.max_rel_diff(a), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRSerializeTest
  : public UnitTest
{
public:
   explicit SparseMatrixCSCRSerializeTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    SparseMatrixCSCR<DT_, IT_> a(10, 10, 3, 6);
    {
      Memory::TypedView<DT_> av = a.val_view_w();
      Memory::TypedView<IT_> ci = a.col_idx_view_w();
      Memory::TypedView<IT_> rp = a.row_ptr_view_w();
      Memory::TypedView<IT_> ri = a.row_idx_view_w();
      for (Index i(0) ; i < a.num_nzes() ; ++i)
      {
        av[i] = DT_(i+1);
        ci[i] = IT_(i + 3);
      }
      rp[0] = IT_(0);
      rp[1] = IT_(2);
      rp[2] = IT_(5);
      rp[3] = IT_(6);
      ri[0] = IT_(0);
      ri[1] = IT_(4);
      ri[2] = IT_(8);
    }

    BinaryStream bs;
    a.write_out(FileMode::fm_cscr, bs);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));
    bs.seekg(0);
    SparseMatrixCSCR<DT_, IT_> e(FileMode::fm_cscr, bs);
    TEST_CHECK_EQUAL(e.num_rows(), 10);
    TEST_CHECK_EQUAL(e.num_cols(), 10);
    TEST_CHECK_EQUAL(e.num_nonzero_rows(), 3);
    TEST_CHECK(e.same_layout(a));
    TEST_CHECK_LESS_THAN(e.max_rel_diff(a), tol);

    auto l = a.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixCSCR<DT_, IT_> tst(l);
    TEST_CHECK_EQUAL(tst.num_rows(), a.num_rows());
    TEST_CHECK_EQUAL(tst.num_cols(), a.num_cols());
    TEST_CHECK_EQUAL(tst.num_nonzero_rows(), a.num_nonzero_rows());
    TEST_CHECK(tst.same_layout(a));
    TEST_CHECK_LESS_THAN(tst.max_rel_diff(a), tol);

#ifdef FEAT_HAVE_ZLIB
    auto zl = a.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixCSCR<DT_, IT_> zlib(zl);
    TEST_CHECK_EQUAL(zlib.num_rows(), a.num_rows());
    TEST_CHECK_EQUAL(zlib.num_cols(), a.num_cols());
    TEST_CHECK_EQUAL(zlib.num_nonzero_rows(), a.num_nonzero_rows());
    TEST_CHECK(zlib.same_layout(a));
    TEST_CHECK_LESS_THAN(zlib.max_rel_diff(a), tol);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = a.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixCSCR<DT_, IT_> zfp(zf);
    TEST_CHECK_EQUAL(zfp.num_rows(), a.num_rows());
    TEST_CHECK_EQUAL(zfp.num_cols(), a.num_cols());
    TEST_CHECK_EQUAL(zfp.num_nonzero_rows(), a.num_nonzero_rows());
    TEST_CHECK(zfp.same_layout(a));
    TEST_CHECK_LESS_THAN(zfp.max_rel_diff(a), DT_(1e-7));
#endif

    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRApplyTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSCRApplyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ alpha(DT_(0.693));

    //for (Index size(3); size < Index(100); size *= 3)
    for(Index nzrows(1); nzrows < Index(5); nzrows *= 3)
    {
      const Index size = nzrows*3;
      SparseMatrixCSCR<DT_, IT_> a(size, size+2, nzrows, 3*nzrows);
      DenseVector<DT_, IT_> x(size+2);
      DenseVector<DT_, IT_> y(size);
      DenseVector<DT_, IT_> r(size);
      DenseVector<DT_, IT_> ref1(size);
      DenseVector<DT_, IT_> ref2(size);

      ref1.format();
      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> row_idx(a.row_idx_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());
        Memory::TypedView<DT_> vx(x.elements_view_w());
        Memory::TypedView<DT_> vy(y.elements_view_w());
        Memory::TypedView<DT_> vr1(ref1.elements_view_w());
        Memory::TypedView<DT_> vr2(ref2.elements_view_w());

        for (Index i(0); i < size+2; ++i)
          vx[i] = DT_((i*524287ll) % 97) * DT_(2.303) - DT_(0.7);

        for (Index i(0); i < size; ++i)
          vy[i] = DT_((i*138937ll) % 73) * DT_(1.618) - DT_(0.6);

        for (Index i(0); i < nzrows; ++i)
        {
          row_ptr[i] = IT_(3*i);
          row_idx[i] = IT_(3*i+1);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);

          vr1[row_idx[i]] = a_val[3*i+0] * vx[i+0] + a_val[3*i+1] * vx[i+1] + a_val[3*i+2] * vx[i+2];
        }
        row_ptr[nzrows] = IT_(3*nzrows);

        for (Index i(0); i < size; ++i)
          vr2[i] = vy[i] + alpha*vr1[i];
      }

      a.apply(r, x);
      TEST_CHECK_LESS_THAN(r.max_rel_diff(ref1), tol);

      a.apply(r, x, y, DT_(0.0));
      TEST_CHECK_LESS_THAN(r.max_rel_diff(y), tol);

      a.apply(r, x, y, alpha);
      TEST_CHECK_LESS_THAN(r.max_rel_diff(ref2), tol);

      r.copy(y);
      a.apply(r, x, r, alpha);
      TEST_CHECK_LESS_THAN(r.max_rel_diff(ref2), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, double, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, double, std::uint32_t, PreferredBackend::cuda);
//#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, Half, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTest, Half, std::uint64_t, PreferredBackend::cuda);
//#endif
//#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRApplyTransposedTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSCRApplyTransposedTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRApplyTransposedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ alpha(DT_(0.693));

    for(Index nzrows(1); nzrows < Index(5); nzrows *= 3)
    {
      const Index size = nzrows*3;
      SparseMatrixCSCR<DT_, IT_> a(size, size+2, nzrows, 3*nzrows);
      DenseVector<DT_, IT_> x(size);
      DenseVector<DT_, IT_> y(size+2);
      DenseVector<DT_, IT_> r(size+2);
      DenseVector<DT_, IT_> ref1(size+2);
      DenseVector<DT_, IT_> ref2(size+2);
      ref1.format();

      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> row_idx(a.row_idx_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());
        Memory::TypedView<DT_> vx(x.elements_view_w());
        Memory::TypedView<DT_> vy(y.elements_view_w());
        Memory::TypedView<DT_> vr1(ref1.elements_view_w());
        Memory::TypedView<DT_> vr2(ref2.elements_view_w());

        for (Index i(0); i < size; ++i)
          vx[i] = DT_((i*524287ll) % 97) * DT_(2.303) - DT_(0.7);

        for (Index i(0); i < size+2; ++i)
          vy[i] = DT_((i*138937ll) % 73) * DT_(1.618) - DT_(0.6);

        for (Index i(0); i < nzrows; ++i)
        {
          row_ptr[i] = IT_(3*i);
          row_idx[i] = IT_(3*i+1);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);
          //vx[i] = DT_((i*524287ll) % 97) * DT_(2.303) - DT_(0.7);

          vr1[i+0] += a_val[3*i+0] * vx[row_idx[i]];
          vr1[i+1] += a_val[3*i+1] * vx[row_idx[i]];
          vr1[i+2] += a_val[3*i+2] * vx[row_idx[i]];
        }

        for (Index i(0); i < size+2; ++i)
          vr2[i] = vy[i] + alpha*vr1[i];

        row_ptr[nzrows] = IT_(3*nzrows);
      }

      a.apply_transposed(r, x);
      TEST_CHECK_LESS_THAN(r.max_rel_diff(ref1), tol);

      a.apply_transposed(r, x, y, DT_(0.0));

      TEST_CHECK_LESS_THAN(r.max_rel_diff(y), tol);

      a.apply_transposed(r, x, y, alpha);
      TEST_CHECK_LESS_THAN(r.max_rel_diff(ref2), tol);

      r.copy(y);
      a.apply_transposed(r, x, r, alpha);
      TEST_CHECK_LESS_THAN(r.max_rel_diff(ref2), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, double, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, double, std::uint32_t, PreferredBackend::cuda);
//#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, Half, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRApplyTransposedTest, Half, std::uint64_t, PreferredBackend::cuda);
//#endif
//#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRMaxRelDiffTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSCRMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(24); size < Index(100); size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size + 2);

      // fill only every 3rd row
      for (Index row(2); row < a_fac.num_rows(); row += 3)
      {
        a_fac.add(row, row, DT_(0));
        a_fac.add(row, row+1, DT_(0));
        a_fac.add(row, row+2, DT_(0));
      }

      SparseMatrixCSCR<DT_, IT_> a, b;
      a.convert(a_fac.make_csr());
      b.clone(a, LAFEM::CloneMode::Layout);

      const Index nnzes = a.num_nzes();

      const Index off0 = (3*nnzes) / 8;
      const Index off1 = (1*nnzes) / 8;
      const Index off2 = (6*nnzes) / 8;

      // a = i, b = i
      {
        Memory::TypedView<DT_> va = a.val_view_w();
        Memory::TypedView<DT_> vb = b.val_view_w();
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
      a.val_view_rw()[off0] += delta_a0;
      b.val_view_rw()[off0] -= delta_b0;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref0, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref0, tol);

      const DT_ delta1 = DT_(0.17);
      const DT_ ref1 = delta1 / (DT_(off0 - off1)*DT_(0.246) + delta1 + DT_(1));
      a.val_view_rw()[off1] -= delta1;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref1, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref1, tol);

      const DT_ delta2 = DT_(0.73);
      const DT_ ref2 = delta2 / (DT_(off2 - off0)*DT_(0.246) + delta2 + DT_(1));
      b.val_view_rw()[off2] += delta2;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref2, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref2, tol);
    }
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::generic);
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::mkl);
//#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRSameLayoutTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSCRSameLayoutTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const Index size = 10;
    const Index diff_row = 4;
    const Index diff_col = 6;
    const DT_ initial_value = DT_(10.0);

    // ref matrix a
    SparseMatrixFactory<DT_, IT_> fac_a(size, size);
    fac_a.add(diff_row+1, diff_col, DT_(1));
    fac_a.add(diff_row, diff_col-1, DT_(2));
    fac_a.add(diff_row, diff_col, initial_value);

    // convert to SparseMatrixCSCR
    SparseMatrixCSCR<DT_, IT_> a;
    a.convert(fac_a.make_csr());

    // weak copy
    auto b = a.clone(CloneMode::Weak);
    TEST_CHECK(a.same_layout(b));

    // shallow copy
    auto c = a.clone(CloneMode::Shallow);
    TEST_CHECK(a.same_layout(c));

    // different values at same position
    SparseMatrixFactory<DT_, IT_> fac_d(size, size);
    fac_d.add(diff_row+1, diff_col, DT_(1));
    fac_d.add(diff_row, diff_col-1, DT_(2));
    fac_d.add(diff_row, diff_col, DT_(0.5));
    SparseMatrixCSCR<DT_, IT_> d;
    d.convert(fac_d.make_csr());
    TEST_CHECK(a.same_layout(d));

    // value at different position
    SparseMatrixFactory<DT_, IT_> fac_e(size, size);
    fac_e.add(diff_row + 1, diff_col, DT_(1));
    fac_e.add(diff_row, diff_col-1, DT_(2));
    fac_e.add(diff_row + 2, diff_col, initial_value);
    SparseMatrixCSCR<DT_, IT_> e;
    e.convert(fac_e.make_csr());
    TEST_CHECK(!a.same_layout(e));

    // different sizes
    SparseMatrixFactory<DT_, IT_> fac_f(size + 2, size);
    fac_f.add(diff_row + 1, diff_col, DT_(1));
    fac_f.add(diff_row, diff_col-1, DT_(2));
    fac_f.add(diff_row, diff_col, initial_value);
    SparseMatrixCSR<DT_, IT_> f_csr(size + 2, size);
    SparseMatrixCSCR<DT_, IT_> f;
    f.convert(fac_f.make_csr());
    TEST_CHECK(!a.same_layout(f));
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSCRSameLayoutTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
