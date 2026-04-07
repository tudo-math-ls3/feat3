// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/random.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

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
class SparseMatrixCSRTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();
    SparseMatrixCSR<DT_, IT_> zero;
    TEST_CHECK(zero.empty());

    SparseMatrixCSR<DT_, IT_> zero3(10, 11, 12);
    TEST_CHECK_EQUAL(zero3.num_nzes(), 12);
    TEST_CHECK_EQUAL(zero3.num_rows(), 10);
    TEST_CHECK_EQUAL(zero3.num_cols(), 11);


    SparseMatrixCSR<DT_, IT_> empty1(11, 12, 111);
    SparseMatrixCSR<DT_, IT_> empty2;
    SparseMatrixCSR<DT_, IT_> empty3;
    empty2.convert(empty1);
    empty3.convert(empty2);
    TEST_CHECK_EQUAL(empty1.num_rows(), empty3.num_rows());
    TEST_CHECK_EQUAL(empty1.num_cols(), empty3.num_cols());
    TEST_CHECK_EQUAL(empty1.num_nzes(), empty3.num_nzes());

    SparseMatrixCSR<DT_, IT_> empty4(empty2.layout());
    SparseMatrixCSR<DT_, IT_> empty5(empty3.layout());
    empty4.convert(empty1);
    empty5.convert(empty4);
    TEST_CHECK_EQUAL(empty5.num_rows(), empty5.num_rows());
    TEST_CHECK_EQUAL(empty5.num_cols(), empty5.num_cols());
    TEST_CHECK_EQUAL(empty5.num_nzes(), empty5.num_nzes());
    empty5.convert(zero);
    TEST_CHECK_EQUAL(empty5.num_rows(), 0);
    TEST_CHECK_EQUAL(empty5.num_cols(), 0);
    TEST_CHECK_EQUAL(empty5.num_nzes(), 0);

    SparseMatrixFactory<DT_, IT_> a01(IT_(10), IT_(10));
    SparseMatrixCSR<DT_, IT_> b01(a01.make_csr());
    TEST_CHECK_EQUAL(b01.num_rows(), 10ul);
    TEST_CHECK_EQUAL(b01.num_cols(), 10ul);
    TEST_CHECK_EQUAL(b01.num_nzes(), 0ul);

    SparseMatrixFactory<DT_, IT_> a(IT_(10), IT_(10));
    a.add(IT_(1), IT_(2), DT_(7));
    a.add(IT_(5), IT_(5), DT_(2));
    a.add(IT_(5), IT_(7), DT_(3));
    a.add(IT_(5), IT_(2), DT_(4));
    SparseMatrixCSR<DT_, IT_> b(a.make_csr());
    TEST_CHECK(!b.empty());
    TEST_CHECK_EQUAL(b.num_nzes(), a.num_nzes());
    {
      const Memory::TypedView<DT_> b_val = b.val_view_r();
      const Memory::TypedView<IT_> b_rp = b.row_ptr_view_r();
      const Memory::TypedView<IT_> b_ci = b.col_idx_view_r();

      TEST_CHECK_EQUAL(b_val(0), DT_(7));
      TEST_CHECK_EQUAL(b_val(1), DT_(4));
      TEST_CHECK_EQUAL(b_val(2), DT_(2));
      TEST_CHECK_EQUAL(b_val(3), DT_(3));
    }

    Index bandw, bandw_idx;
    b.bandwidth_row(bandw, bandw_idx);
    TEST_CHECK_EQUAL(bandw, Index(6));
    TEST_CHECK_EQUAL(bandw_idx, Index(5));
    //b.bandwidth_column(bandw, bandw_idx);
    //TEST_CHECK_EQUAL(bandw, Index(5));
    //TEST_CHECK_EQUAL(bandw_idx, Index(2));

    Index radius, radius_idx;
    b.radius_row(radius, radius_idx);
    TEST_CHECK_EQUAL(radius, Index(3));
    TEST_CHECK_EQUAL(radius_idx, Index(5));
    //b.radius_column(radius, radius_idx);
    //TEST_CHECK_EQUAL(radius, Index(3));
    //TEST_CHECK_EQUAL(radius_idx, Index(2));

    SparseMatrixCSR<DT_, IT_> bl(b.layout());
    TEST_CHECK_EQUAL(bl.num_nzes(), b.num_nzes());
    TEST_CHECK_EQUAL(bl.num_rows(), b.num_rows());
    TEST_CHECK_EQUAL(bl.num_cols(), b.num_cols());

    bl = b.layout();
    TEST_CHECK_EQUAL(bl.num_nzes(), b.num_nzes());
    TEST_CHECK_EQUAL(bl.num_rows(), b.num_rows());
    TEST_CHECK_EQUAL(bl.num_cols(), b.num_cols());

    typename SparseLayout<IT_, SparseLayoutId::lt_csr>::template MatrixType<DT_> x(b.layout());
    TEST_CHECK_EQUAL(x.row_ptr_arbiter(), b.row_ptr_arbiter());
    TEST_CHECK_NOT_EQUAL(x.val_arbiter(), b.val_arbiter());
    /// \compilerhack icc 14.x and msvc do not understand the following single line, so we need a typedef detour here
#if defined(FEAT_COMPILER_MICROSOFT) || (defined(FEAT_COMPILER_INTEL) && __INTEL_COMPILER < 1500)
    typedef decltype(b.layout()) LayoutId;
    typename LayoutId::template MatrixType<DT_> y(b.layout());
#else
    typename decltype(b.layout())::template MatrixType<DT_> y(b.layout());
#endif
    TEST_CHECK_EQUAL(y.row_ptr_arbiter(), b.row_ptr_arbiter());
    TEST_CHECK_NOT_EQUAL(y.val_arbiter(), b.val_arbiter());


    SparseMatrixCSR<DT_, IT_> z;
    z.convert(b);
    TEST_CHECK_EQUAL(z.num_nzes(), 4ul);
    TEST_CHECK_EQUAL(z.num_rows(), a.num_rows());
    TEST_CHECK_EQUAL(z.num_cols(), a.num_cols());
    {
      Memory::TypedView<DT_> z_view = z.val_view_r();
      TEST_CHECK_EQUAL(z_view(0), DT_(7));
      TEST_CHECK_EQUAL(z_view(1), DT_(4));
    }

    SparseMatrixCSR<DT_, IT_> c;
    c.clone(b);
    TEST_CHECK_NOT_EQUAL(c.val_arbiter(), b.val_arbiter());
    TEST_CHECK_EQUAL(c.col_idx_arbiter(), b.col_idx_arbiter());
    c = b.clone(CloneMode::Deep);
    TEST_CHECK_NOT_EQUAL(c.val_arbiter(), b.val_arbiter());
    TEST_CHECK_NOT_EQUAL(c.col_idx_arbiter(), b.col_idx_arbiter());

    DenseVector<IT_, IT_> col_idx(c.num_nzes(), c.col_idx_arbiter().attach());
    DenseVector<DT_, IT_> val(c.num_nzes(), c.val_arbiter().attach());
    DenseVector<IT_, IT_> row_ptr(c.num_rows() + 1, c.row_ptr_arbiter().attach());
    SparseMatrixCSR<DT_, IT_> d(c.num_rows(), c.num_cols(), row_ptr, col_idx, val);
    TEST_CHECK_LESS_THAN(d.max_rel_diff(c), eps);

    SparseMatrixCSR<DT_, IT_> e;
    e.convert(c);
    TEST_CHECK_LESS_THAN(e.max_rel_diff(c), eps);
    e.copy(c);
    TEST_CHECK_LESS_THAN(e.max_rel_diff(c), eps);
    e.clone(c);
    b.clone(e);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(c), eps);

    // new clone testing
    auto clone1 = b.clone(CloneMode::Deep);
    TEST_CHECK_LESS_THAN(clone1.max_rel_diff(b), eps);
    clone1.val_view_rw()[1] = DT_(132);
    TEST_CHECK_LESS_THAN(eps, clone1.max_rel_diff(b));
    TEST_CHECK_NOT_EQUAL(clone1.val_arbiter(), b.val_arbiter());
    TEST_CHECK_NOT_EQUAL(clone1.row_ptr_arbiter(), b.row_ptr_arbiter());
    auto clone2 = clone1.clone(CloneMode::Layout);
    clone2.format(DT_(4713));
    //TEST_CHECK_NOT_EQUAL(clone2(5, 5), clone1(5, 5));
    TEST_CHECK_NOT_EQUAL(clone2.val_view_r()(1), clone1.val_view_r()(1));
    TEST_CHECK_NOT_EQUAL(clone2.val_arbiter(), clone1.val_arbiter());
    TEST_CHECK_EQUAL(clone2.row_ptr_arbiter(), clone1.row_ptr_arbiter());
    auto clone3 = clone1.clone(CloneMode::Weak);
    TEST_CHECK_LESS_THAN(clone3.max_rel_diff(clone1), eps);
    clone3.val_view_rw()[1] = DT_(133);
    TEST_CHECK_LESS_THAN(eps, clone3.max_rel_diff(clone1));
    TEST_CHECK_NOT_EQUAL(clone3.val_arbiter(), clone1.val_arbiter());
    TEST_CHECK_EQUAL(clone3.row_ptr_arbiter(), clone1.row_ptr_arbiter());
    auto clone4 = clone1.clone(CloneMode::Shallow);
    TEST_CHECK_LESS_THAN(clone4.max_rel_diff(clone1), eps);
    clone4.val_view_rw()[1] = DT_(134);
    TEST_CHECK_LESS_THAN(clone4.max_rel_diff(clone1), eps);
    TEST_CHECK_EQUAL(clone4.val_arbiter(), clone1.val_arbiter());
    TEST_CHECK_EQUAL(clone4.row_ptr_arbiter(), clone1.row_ptr_arbiter());

    // shrink test
    SparseMatrixFactory<DT_, IT_> ffac(IT_(10), IT_(12));
    for (IT_ row(0); row < ffac.num_rows(); ++row)
    {
      ffac.add(row, row, DT_(2));
      ffac.add(row, row+1, DT_(1));
    }
    SparseMatrixCSR<DT_, IT_> f(ffac.make_csr());
    f.shrink(DT_(1.5));
    TEST_CHECK_EQUAL(f.num_nzes(), 10ul);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRSerializeTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRSerializeTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol_b = TestSystem::tol<DT_>();
    const DT_ tol_t = TestSystem::tol<DT_>();

    Index size = 37;
    SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
    {
      Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
      Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
      Memory::TypedView<DT_> a_val(a.val_view_w());

      for (Index i(0); i < size; ++i)
      {
        row_ptr[i] = IT_(3*i);
        col_idx[3*i+0] = IT_(i);
        col_idx[3*i+1] = IT_(i+1);
        col_idx[3*i+2] = IT_(i+2);
        a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
        a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
        a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);
      }
      row_ptr[size] = IT_(3*size);
    }

    SparseMatrixFactory<DT_, IT_> a_sym_fac(size, size);
    for (Index i(0); i < size; ++i)
    {
      a_sym_fac.add(i, i, DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3));
      if(i > Index(0))
        a_sym_fac.add(i, i-1, DT_(((i-1)*65537ll) % 37) * DT_(2.713) - DT_(0.2));
      if(i+1 < size)
        a_sym_fac.add(i, i+1, DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2));
    }
    SparseMatrixCSR<DT_, IT_> a_sym(a_sym_fac.make_csr());

    BinaryStream bs;
    a.write_out(FileMode::fm_csr, bs);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));
    bs.seekg(0);
    SparseMatrixCSR<DT_, IT_> g(FileMode::fm_csr, bs);
    TEST_CHECK_LESS_THAN(g.max_rel_diff(a), tol_b);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));

    // unsymmetric MTX format
    std::stringstream ts;
    a.write_out(FileMode::fm_mtx, ts);
    SparseMatrixCSR<DT_, IT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_LESS_THAN(j.max_rel_diff(a), tol_t);

    // symmetric MTX format
    std::stringstream ts2;
    a_sym.write_out(FileMode::fm_mtx, ts2, true);
    SparseMatrixCSR<DT_, IT_> j2(FileMode::fm_mtx, ts2);
    TEST_CHECK_LESS_THAN(j2.max_rel_diff(a_sym), tol_t);

    auto kp = a.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixCSR<DT_, IT_> k(kp);
    TEST_CHECK_LESS_THAN(k.max_rel_diff(a), tol_b);
#ifdef FEAT_HAVE_ZLIB
    auto zl = a.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixCSR<DT_, IT_> zlib(zl);
    TEST_CHECK_LESS_THAN(zlib.max_rel_diff(a), tol_b);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = a.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixCSR<DT_, IT_> zfp(zf);
    TEST_CHECK_LESS_THAN(zlib.max_rel_diff(a), DT_(1e-7));
#endif
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRApplyTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRApplyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ alpha(DT_(0.693));

    for (Index size(1); size < Index(100); size *= 3)
    {
      SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
      DenseVector<DT_, IT_> x(size+2);
      DenseVector<DT_, IT_> y(size);
      DenseVector<DT_, IT_> r(size);
      DenseVector<DT_, IT_> ref1(size);
      DenseVector<DT_, IT_> ref2(size);

      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());
        Memory::TypedView<DT_> vx(x.elements_view_w());
        Memory::TypedView<DT_> vy(y.elements_view_w());
        Memory::TypedView<DT_> vr1(ref1.elements_view_w());
        Memory::TypedView<DT_> vr2(ref2.elements_view_w());

        for (Index i(0); i < size+2; ++i)
          vx[i] = DT_((i*524287ll) % 97) * DT_(2.303) - DT_(0.7);

        for (Index i(0); i < size; ++i)
        {
          row_ptr[i] = IT_(3*i);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);
          vy[i] = DT_((i*138937ll) % 73) * DT_(1.618) - DT_(0.6);
          vr1[i] = a_val[3*i+0] * vx[i+0] + a_val[3*i+1] * vx[i+1] + a_val[3*i+2] * vx[i+2];
          vr2[i] = vy[i] + alpha*vr1[i];
        }
        row_ptr[size] = IT_(3*size);
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

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, double, std::uint32_t, PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRApplyTransposedTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRApplyTransposedTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRApplyTransposedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ alpha(DT_(0.693));

    for (Index size(1); size < Index(100); size *= 3)
    {
      SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
      DenseVector<DT_, IT_> x(size);
      DenseVector<DT_, IT_> y(size+2);
      DenseVector<DT_, IT_> r(size+2);
      DenseVector<DT_, IT_> ref1(size+2);
      DenseVector<DT_, IT_> ref2(size+2);

      ref1.format();

      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());
        Memory::TypedView<DT_> vx(x.elements_view_w());
        Memory::TypedView<DT_> vy(y.elements_view_w());
        Memory::TypedView<DT_> vr1(ref1.elements_view_w());
        Memory::TypedView<DT_> vr2(ref2.elements_view_w());

        for (Index i(0); i < size+2; ++i)
        {
          vy[i] = DT_((i*138937ll) % 73) * DT_(1.618) - DT_(0.6);
        }

        for (Index i(0); i < size; ++i)
        {
          row_ptr[i] = IT_(3*i);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);
          vx[i] = DT_((i*524287ll) % 97) * DT_(2.303) - DT_(0.7);
          vr1[i+0] += a_val[3*i+0] * vx[i];
          vr1[i+1] += a_val[3*i+1] * vx[i];
          vr1[i+2] += a_val[3*i+2] * vx[i];
        }
        row_ptr[size] = IT_(3*size);

        for (Index i(0); i < size+2; ++i)
          vr2[i] = vy[i] + alpha*vr1[i];
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

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, double, std::uint32_t, PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRApplyTransposedTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRBApplyTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRBApplyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRBApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ alpha(DT_(0.693));

    for (IT_ size(1); size < Index(100); size *= IT_(3))
    {
      SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
      DenseVectorBlocked<DT_, IT_, 3> x(size+2);
      DenseVectorBlocked<DT_, IT_, 3> y(size);
      DenseVectorBlocked<DT_, IT_, 3> r(size);
      DenseVectorBlocked<DT_, IT_, 3> ref1(size);
      DenseVectorBlocked<DT_, IT_, 3> ref2(size);

      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());
        Memory::TypedView<Tiny::Vector<DT_, 3>> vx(x.elements_view_w());
        Memory::TypedView<Tiny::Vector<DT_, 3>> vy(y.elements_view_w());
        Memory::TypedView<Tiny::Vector<DT_, 3>> vr1(ref1.elements_view_w());
        Memory::TypedView<Tiny::Vector<DT_, 3>> vr2(ref2.elements_view_w());

        for (Index i(0); i < size+2; ++i)
        {
          vx[i][0] = DT_((i*524287ll) % 97) * DT_(2.303) - DT_(0.7);
          vx[i][1] = DT_(0.7) * vx[i][0];
          vx[i][2] = DT_(1.3) * vx[i][0];
        }

        for (Index i(0); i < size; ++i)
        {
          row_ptr[i] = IT_(3*i);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);

          vy[i][0] = DT_((i*138937ll) % 73) * DT_(1.618) - DT_(0.6);
          vy[i][1] = vy[i][0] * DT_(0.3);
          vy[i][2] = vy[i][0] * DT_(1.7);

          vr1[i][0] = a_val[3*i+0] * vx[i+0][0] + a_val[3*i+1] * vx[i+1][0] + a_val[3*i+2] * vx[i+2][0];
          vr1[i][1] = DT_(0.7) * vr1[i][0];
          vr1[i][2] = DT_(1.3) * vr1[i][0];

          vr2[i][0] = vy[i][0] + alpha * vr1[i][0];
          vr2[i][1] = vy[i][1] + alpha * vr1[i][1];
          vr2[i][2] = vy[i][2] + alpha * vr1[i][2];
        }
        row_ptr[size] = IT_(3*size);
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

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRBApplyTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRScaleTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRScaleTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ alpha(DT_(0.693));

    for (Index size(1); size < Index(300); size *= 2)
    {
      SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());

        for (Index i(0); i < size; ++i)
        {
          row_ptr[i] = IT_(3*i);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);
        }
        row_ptr[size] = IT_(3*size);
      }

      SparseMatrixCSR<DT_, IT_> ref(a.clone(LAFEM::CloneMode::Weak));
      {
        Memory::TypedView<DT_> r_val(ref.val_view_rw());
        for(Index i(0); i < ref.num_nzes(); ++i)
          r_val[i] *= alpha;
      }

      SparseMatrixCSR<DT_, IT_> b(a.clone(LAFEM::CloneMode::Weak));

      b.scale(a, alpha);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref), tol);

      a.scale(a, alpha);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRScaleRowColTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRScaleRowColTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRScaleRowColTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }


  virtual void run() const override
  {
    const DT_ pi(Math::pi<DT_>());
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(2); size < Index(300); size *= 3)
    {
      DenseVector<DT_, IT_> sr(size);
      DenseVector<DT_, IT_> sc(size + 2);

      {
        Memory::TypedView<DT_> vs(sr.elements_view_w());
        for (Index i(0); i < sr.size(); ++i)
        {
          vs[i] = pi * DT_(i % 3 + 1) - DT_(5.21) + DT_(i % 17);
        }
      }

      {
        Memory::TypedView<DT_> vs(sc.elements_view_w());
        for (Index i(0); i < sc.size(); ++i)
        {
          vs[i] = pi * DT_(i % 5 + 2) - DT_(7.32) + DT_(i % 13);
        }
      }

      SparseMatrixFactory<DT_, IT_> a_fac(size, size + 2);
      SparseMatrixFactory<DT_, IT_> ar_fac(size, size + 2);
      SparseMatrixFactory<DT_, IT_> ac_fac(size, size + 2);
      for (IT_ row(0); row < a_fac.num_rows(); ++row)
      {
        {
          Index col = row;
          a_fac.add(row, col, DT_(2));
          ar_fac.add(row, col, DT_(2) * (pi * DT_(row % 3 + 1) - DT_(5.21) + DT_(row % 17)));
          ac_fac.add(row, col, DT_(2) * (pi * DT_(col % 5 + 2) - DT_(7.32) + DT_(col % 13)));
        }
        {
          Index col = row+1;
          a_fac.add(row, col,  DT_(col % 7));
          ar_fac.add(row, col, DT_(col % 7) * (pi * DT_(row % 3 + 1) - DT_(5.21) + DT_(row % 17)));
          ac_fac.add(row, col, DT_(col % 7) * (pi * DT_(col % 5 + 2) - DT_(7.32) + DT_(col % 13)));
        }
        {
          Index col = row+2;
          a_fac.add(row, col,  DT_(col % 11));
          ar_fac.add(row, col, DT_(col % 11) * (pi * DT_(row % 3 + 1) - DT_(5.21) + DT_(row % 17)));
          ac_fac.add(row, col, DT_(col % 11) * (pi * DT_(col % 5 + 2) - DT_(7.32) + DT_(col % 13)));
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> ar(ar_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> ac(ac_fac.make_csr());

      SparseMatrixCSR<DT_, IT_> b(a.clone(LAFEM::CloneMode::Layout));

      // Scale rows
      b.scale_rows(a, sr);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ar), tol);

      // Scale cols
      b.scale_cols(a, sc);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ac), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRScaleRowColTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRTranspositionTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRTranspositionTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRTranspositionTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(2); size < Index(200); size *= 3)
    {
      const Index m = size;
      const Index n = size + 2;
      SparseMatrixFactory<DT_, IT_> a_fac(m, n), b_fac(n ,m);

      for(Index i = 0; i < m; ++i)
      {
        const Index j1 = (i * Index(131071)) % n;
        const Index j2 = (i * Index(524287)) % n;
        const Index j3 = (i * Index(2147483647)) % n;
        a_fac.add(i, j1, DT_(3*i + 1));
        a_fac.add(i, j2, DT_(3*i + 2));
        a_fac.add(i, j3, DT_(3*i + 3));
        b_fac.add(j1, i, DT_(3*i + 1));
        b_fac.add(j2, i, DT_(3*i + 2));
        b_fac.add(j3, i, DT_(3*i + 3));
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());

      SparseMatrixCSR<DT_, IT_> c = a.transpose();

      TEST_CHECK_EQUAL(c.num_rows(), b.num_rows());
      TEST_CHECK_EQUAL(c.num_cols(), b.num_cols());
      TEST_CHECK_EQUAL(c.num_nzes(), b.num_nzes());

      TEST_CHECK(c.same_layout(b));

      TEST_CHECK_LESS_THAN(c.max_rel_diff(b), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRTranspositionTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRPermuteTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRPermuteTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRPermuteTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (IT_ size(25); size < Index(1000); size *= IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size);

      DenseVector<DT_, IT_> x(size);
      {
        Memory::TypedView<DT_> vx(x.elements_view_w());
        for (Index i(0); i < size; ++i)
        {
          vx[i] = DT_(i % 100) * DT_(1.234);
        }
      }

      for (IT_ row(0); row < a_fac.num_rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.num_cols(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
          }
          else if ((row == col + 7) || (row + 7 == col))
          {
            a_fac.add(row, col, DT_(-1));
          }
          else if ((row == col + 15) || (row + 15 == col))
          {
            a_fac.add(row, col, DT_(-2));
          }
          else if ((row == col + a_fac.num_cols() / 2) || (row + a_fac.num_rows() / 2 == col))
          {
            a_fac.add(row, col, DT_(1));
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());

      DenseVector<DT_, IT_> r(size);
      a.apply(r, x);
      DT_ ref_norm = r.norm2();

      auto a_backup = a.clone(CloneMode::Deep);

      Random rng;
      std::cout << "RNG Seed: " << rng.get_seed() << "\n";
      Adjacency::Permutation perm(a.num_rows(), rng);

      a.permute(perm, perm);
      x.permute(perm);

      a.apply(r, x);
      DT_ norm = r.norm2();
      DT_ deviation = norm / ref_norm;
      TEST_CHECK(deviation > DT_(0.99));
      TEST_CHECK(deviation < DT_(1.01));

      a = a_backup.clone(CloneMode::Deep);
      auto perm_inv = perm.inverse();
      a.permute(perm_inv, perm);
      a.permute(perm, perm_inv);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(a_backup), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRPermuteTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRDiagTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRDiagTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRDiagTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(2); size < Index(200); size *= 3)
    {
      SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
      DenseVector<DT_, IT_> ref(size);
      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());
        Memory::TypedView<DT_> vr(ref.elements_view_w());

        for (Index i(0); i < size; ++i)
        {
          row_ptr[i] = IT_(3*i);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);
          vr[i] = a_val[3*i+0];
        }
        row_ptr[size] = IT_(3*size);
      }

      auto diag = a.extract_diag();
      TEST_CHECK_LESS_THAN(diag.max_rel_diff(ref), tol);
    }

    SparseMatrixFactory<DT_, IT_> b_fac(16, 16);
    b_fac.add(0, 0, DT_(1.0));
    b_fac.add(2, 2, DT_(2.0));
    b_fac.add(3, 3, DT_(3.0));
    b_fac.add(7, 7, DT_(7.0));

    SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());

    auto ref = b.create_vector_l();
    ref.format(DT_(0.0));
    {
      Memory::TypedView<DT_> vr(ref.elements_view_rw());
      vr[0] = DT_(1.0);
      vr[2] = DT_(2.0);
      vr[3] = DT_(3.0);
      vr[7] = DT_(7.0);
    }
    auto diag = b.extract_diag();
    TEST_CHECK_LESS_THAN(diag.max_rel_diff(ref), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRDiagTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRAxpyTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRAxpyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ alpha(DT_(0.693));

    for (Index size(2); size < Index(200); size *= 3)
    {

      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> b_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_fac1(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_fac2(IT_(size), IT_(size + 2));

      for (IT_ row(0); row < a_fac.num_rows(); ++row)
      {
        {
          Index col = row;
          a_fac.add(row, col, DT_(2));
          b_fac.add(row, col, DT_((row + col) % 15));
          ref_fac1.add(row, col, DT_(2) * alpha + DT_((row + col) % 15));
          ref_fac2.add(row, col, DT_((row + col) % 15) * (DT_(1) + alpha));
        }

        {
          Index col = row+1;
          a_fac.add(row, col, DT_(-1));
          b_fac.add(row, col, DT_((row + col + 1) % 15));
          ref_fac1.add(row, col, DT_(-1) * alpha + DT_((row + col + 1) % 15));
          ref_fac2.add(row, col, DT_((row + col + 1) % 15) * (DT_(1) + alpha));
        }

        {
          Index col = row+2;
          a_fac.add(row, col, DT_(1));
          b_fac.add(row, col, DT_((row + col + 3) % 11));
          ref_fac1.add(row, col, DT_(1) * alpha + DT_((row + col + 3) % 11));
          ref_fac2.add(row, col, DT_((row + col + 3) % 11) * (DT_(1) + alpha));
        }
      }

      SparseMatrixCSR<DT_, IT_> ref1(ref_fac1.make_csr());
      SparseMatrixCSR<DT_, IT_> ref2(ref_fac2.make_csr());
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());

      // r != x
      a.scale(a, alpha);
      a.axpy(b);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref1), tol);

      // r == x
      b.axpy(b, alpha);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref2), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAxpyTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRFrobeniusTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRFrobeniusTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRFrobeniusTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = DT_(10) * TestSystem::tol<DT_>();
    for (Index size(2); size < Index(200); size *= 3)
    {
      DT_ ref = DT_(0);
      SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());

        for (Index i(0); i < size; ++i)
        {
          row_ptr[i] = IT_(3*i);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);
          ref += Math::sqr(a_val[3*i+0]) + Math::sqr(a_val[3*i+1]) + Math::sqr(a_val[3*i+2]);
        }
        row_ptr[size] = IT_(3*size);
      }
      ref = Math::sqrt(ref);

      TEST_CHECK_EQUAL_WITHIN_EPS(a.norm_frobenius(), ref, tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRFrobeniusTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRLumpTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRLumpTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRLumpTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(1); size < Index(100); size *= 3)
    {
      SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
      DenseVector<DT_, IT_> ref(size);
      {
        Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
        Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
        Memory::TypedView<DT_> a_val(a.val_view_w());
        Memory::TypedView<DT_> vr(ref.elements_view_w());

        for (Index i(0); i < size; ++i)
        {
          row_ptr[i] = IT_(3*i);
          col_idx[3*i+0] = IT_(i);
          col_idx[3*i+1] = IT_(i+1);
          col_idx[3*i+2] = IT_(i+2);
          a_val[3*i+0] = DT_((i*131071ll) % 31) * DT_(3.141) - DT_(0.3);
          a_val[3*i+1] = DT_((i*65537ll) % 37) * DT_(2.713) - DT_(0.2);
          a_val[3*i+2] = DT_((i*40487ll) % 53) * DT_(1.414) - DT_(0.1);
          vr[i] = a_val[3*i+0] + a_val[3*i+1] + a_val[3*i+2];
        }
        row_ptr[size] = IT_(3*size);
      }

      DenseVector<DT_, IT_> lump = a.lump_rows();
      TEST_CHECK_LESS_THAN(lump.max_rel_diff(ref), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRLumpTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRCompressionTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRCompressionTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRCompressionTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    Index mat_rows = 56;
    Index mat_cols = 56;
    Index mat_nzes = 400;

    SparseMatrixCSR<DT_, IT_> a(mat_rows, mat_cols, mat_nzes);

    {
      Memory::TypedView<IT_> row_ptr = a.row_ptr_view_w();
      Memory::TypedView<IT_> col_idx = a.col_idx_view_w();
      Memory::TypedView<DT_> val = a.val_view_w();
      for (Index i(0); i < 400; ++i)
      {
        col_idx[i] = IT_(i / 40 + (3 * i) % 16);
        val[i] = DT_(7) / DT_(3 * (i + 1)) - DT_(13) / DT_((i + 1) * (i + 1));
      }
      row_ptr[0] = 0;
      for (Index i(1); i < mat_rows; ++i)
      {
        row_ptr[i] = IT_(i * Index(400 / mat_rows - 1));
      }
      row_ptr[mat_rows] = 400;
    }

#ifdef FEAT_HAVE_ZLIB
#ifdef FEAT_HAVE_ZFP
    LAFEM::SerialConfig config(false, false);
    config.set_tolerance(FEAT::Real(1e-2));
    std::vector<char> uncompressed = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_zlib);
    std::vector<char> elements_compressed_zlib = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_zfp);
    std::vector<char> elements_compressed_zfp = a.serialize(config);
    config.set_indices_compression(LAFEM::CompressionModes::indices_zlib);
    std::vector<char> elements_indices_compressed_zfp = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_zlib);
    std::vector<char> elements_indices_compressed_zlib = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_off);
    std::vector<char> indices_compressed_zlib = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_zfp);
    config.set_indices_compression(LAFEM::CompressionModes::indices_off);
    config.set_tolerance(FEAT::Real(1e-4));
    std::vector<char> elements_compressed_zfp_e4 = a.serialize(config);
    config.set_tolerance(FEAT::Real(1e-6));
    std::vector<char> elements_compressed_zfp_e6 = a.serialize(config);

    XASSERTM(uncompressed.size() > elements_compressed_zlib.size(), "ele zlib is not smaller than uncomp");
    XASSERTM(elements_compressed_zlib.size() > elements_compressed_zfp.size(), "ele zfp is not smaleer than ele zlib");
    XASSERTM(elements_compressed_zfp.size() > elements_indices_compressed_zfp.size(), "ele + ind zfp is not smaller than ele zfp");
    XASSERTM(elements_compressed_zlib.size() > elements_indices_compressed_zlib.size(), "ele + ind zlib is not smaller than ele zlib");
    XASSERTM(uncompressed.size() > indices_compressed_zlib.size(), "ind zlib is not smaller than uncomp");
    XASSERTM(indices_compressed_zlib.size() > elements_indices_compressed_zlib.size(), "ele + ind zlib is not smaller than ind zlib");
    XASSERTM(elements_compressed_zfp_e4.size() > elements_compressed_zfp.size(), "e4 smaller than e2");
    XASSERTM(elements_compressed_zfp_e6.size() > elements_compressed_zfp_e4.size(), "e6 smaller than e4");
#endif
#endif
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, __float128, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, __float128, std::uint64_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRCompressionTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRMinMaxElementTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRMinMaxElementTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRMinMaxElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const Index size = 19;

    SparseMatrixCSR<DT_, IT_> a(size, size+2, 3*size);
    {
      Memory::TypedView<IT_> row_ptr(a.row_ptr_view_w());
      Memory::TypedView<IT_> col_idx(a.col_idx_view_w());
      Memory::TypedView<DT_> a_val(a.val_view_w());

      for (Index i(0); i < size; ++i)
      {
        row_ptr[i] = IT_(3*i);
        col_idx[3*i+0] = IT_(i);
        col_idx[3*i+1] = IT_(i+1);
        col_idx[3*i+2] = IT_(i+2);
        a_val[3*i+0] = DT_((i*131071ll) % 31) - DT_(3.5);
        a_val[3*i+1] = (DT_((i*65537ll) % 37) - DT_(1.875));
        a_val[3*i+2] = DT_((i*40487ll) % 53) - DT_(1.25);
      }
      row_ptr[size] = IT_(3*size);
    }

    if(Backend::get_preferred_backend() == PreferredBackend::generic)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(a.min_element(), DT_(-3.5), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(a.max_element(), DT_(49.75), tol);
    }

    TEST_CHECK_EQUAL_WITHIN_EPS(a.min_abs_element(), DT_(0.125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(a.max_abs_element(), DT_(49.75), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMinMaxElementTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<typename DT_, typename IT_>
class SparseMatrixCSRAddDoubleMatMultTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRAddDoubleMatMultTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRAddDoubleMatMultTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    SparseMatrixFactory<DT_, IT_> b_fac(IT_(4), IT_(3));
    SparseMatrixFactory<DT_, IT_> d_fac(IT_(3), IT_(4));
    SparseMatrixFactory<DT_, IT_> s_fac(IT_(3), IT_(3));
    SparseMatrixFactory<DT_, IT_> s_fac2(IT_(3), IT_(3));

    b_fac.add(IT_(0), IT_(0), DT_(2));
    b_fac.add(IT_(0), IT_(1), DT_(3));
    b_fac.add(IT_(1), IT_(0), DT_(4));
    b_fac.add(IT_(1), IT_(2), DT_(5));
    b_fac.add(IT_(2), IT_(1), DT_(6));
    b_fac.add(IT_(2), IT_(2), DT_(7));
    b_fac.add(IT_(3), IT_(1), DT_(8));
    d_fac.add(IT_(0), IT_(0), DT_(5));
    d_fac.add(IT_(0), IT_(1), DT_(7));
    d_fac.add(IT_(1), IT_(0), DT_(3));
    d_fac.add(IT_(1), IT_(2), DT_(2));
    d_fac.add(IT_(2), IT_(1), DT_(8));
    d_fac.add(IT_(2), IT_(2), DT_(6));
    d_fac.add(IT_(2), IT_(3), DT_(4));
    s_fac.add(IT_(0), IT_(0), DT_(104));
    s_fac.add(IT_(0), IT_(1), DT_(30));
    s_fac.add(IT_(0), IT_(2), DT_(105));
    s_fac.add(IT_(1), IT_(0), DT_(12));
    s_fac.add(IT_(1), IT_(1), DT_(66));
    s_fac.add(IT_(1), IT_(2), DT_(56));
    s_fac.add(IT_(2), IT_(0), DT_(96));
    s_fac.add(IT_(2), IT_(1), DT_(304));
    s_fac.add(IT_(2), IT_(2), DT_(288));
    s_fac2.add(IT_(0), IT_(0), DT_(696));
    s_fac2.add(IT_(0), IT_(1), DT_(120));
    s_fac2.add(IT_(0), IT_(2), DT_(770));
    s_fac2.add(IT_(1), IT_(0), DT_(48));
    s_fac2.add(IT_(1), IT_(1), DT_(504));
    s_fac2.add(IT_(1), IT_(2), DT_(504));
    s_fac2.add(IT_(2), IT_(0), DT_(704));
    s_fac2.add(IT_(2), IT_(1), DT_(2896));
    s_fac2.add(IT_(2), IT_(2), DT_(2392));

    SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());
    SparseMatrixCSR<DT_, IT_> d(d_fac.make_csr());
    SparseMatrixCSR<DT_, IT_> s(s_fac.make_csr());

    DenseVector<DT_, IT_> a(4u);
    {
      Memory::TypedView<DT_> va = a.elements_view_w();
      va[0] = DT_(2);
      va[1] = DT_(3);
      va[2] = DT_(4);
      va[3] = DT_(5);
    }

    // perform double matrix product and check result
    s.add_double_mat_product(d, a, b, -DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(s.norm_frobenius(), DT_(0), tol);

    // convert matrices to blocked variants
    SparseMatrixBCSR<DT_, IT_, 2, 1> b2(b.layout());
    SparseMatrixBCSR<DT_, IT_, 1, 2> d2(d.layout());
    SparseMatrixCSR<DT_, IT_> s2(s_fac2.make_csr());
    {
      auto vb2 = b2.val_view_w();
      Memory::TypedView<DT_> vb = b.val_view_r();
      for(Index i(0); i < b.num_nzes();  ++i)
        vb2[i][0][0] = DT_(2) * (vb2[i][1][0] = vb(i));
    }
    {
      auto vd2 = d2.val_view_w();
      Memory::TypedView<DT_> vd = d.val_view_r();
      for(Index i(0); i < d.num_nzes();  ++i)
        vd2[i][0][0] = DT_(3) * (vd2[i][0][1] = vd(i));
    }

    DenseVectorBlocked<DT_, IT_, 2> a2(4u);
    {
      auto va2 = a2.elements_view_w();
      va2[0][0] = DT_(1);
      va2[0][1] = DT_(2);
      va2[1][0] = DT_(3);
      va2[1][1] = DT_(4);
      va2[2][0] = DT_(5);
      va2[2][1] = DT_(6);
      va2[3][0] = DT_(7);
      va2[3][1] = DT_(8);
    }

    s2.add_double_mat_product(d2, a2, b2, -DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(s2.norm_frobenius(), DT_(0), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAddDoubleMatMultTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAddDoubleMatMultTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAddDoubleMatMultTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRAddDoubleMatMultTest, double, std::uint64_t, PreferredBackend::generic);

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRMaxRelDiffTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(12); size < Index(100); size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size + 2);

      for (Index row(0); row < a_fac.num_rows(); ++row)
      {
        a_fac.add(row, row, DT_(0));
        a_fac.add(row, row+1, DT_(0));
        a_fac.add(row, row+2, DT_(0));
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> b(a.clone(LAFEM::CloneMode::Layout));

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
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::generic);
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::mkl);
//#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, double, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRMaxRelDiffTest, double, std::uint64_t, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRSameLayoutTest
  : public UnitTest
{
public:
  explicit SparseMatrixCSRSameLayoutTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    for (Index size(4); size < Index(128); size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> fac_a(size, size);
      for (IT_ i = 0; i < IT_(size); ++i)
        fac_a.add(i, i, DT_(i));
      SparseMatrixCSR<DT_, IT_> a(fac_a.make_csr());

      // weak copy
      auto b = a.clone(CloneMode::Weak);
      TEST_CHECK(a.same_layout(b));

      // shallow copy
      auto c = a.clone(CloneMode::Shallow);
      TEST_CHECK(a.same_layout(c));

      // different values at same position
      SparseMatrixFactory<DT_, IT_> fac_d(size, size);
        for (IT_ i = 0; i < IT_(size); ++i)
          fac_d.add(i, i , DT_(i + 1));
      SparseMatrixCSR<DT_, IT_> d(fac_d.make_csr());
      TEST_CHECK(a.same_layout(d));

      // values at different position
      SparseMatrixFactory<DT_, IT_> fac_e(size, size);
        for (IT_ i = 0; i < IT_(size - 1); ++i)
          fac_e.add(i + 1, i + 1, DT_(i));
      SparseMatrixCSR<DT_, IT_> e(fac_e.make_csr());
      TEST_CHECK(!a.same_layout(e));

      // different sizes
      SparseMatrixFactory<DT_, IT_> fac_f(size + 2, size + 2);
        for (IT_ i = 0; i < IT_(size + 2); ++i)
          fac_f.add(i, i, DT_(i));
      SparseMatrixCSR<DT_, IT_> f(fac_f.make_csr());
      TEST_CHECK(!a.same_layout(f));
    }
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixCSRSameLayoutTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
