// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/util/binary_stream.hpp>

#include <kernel/util/random.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <sstream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse matrix banded class.
 *
 * \test test description missing
 *
 * \tparam DT_
 * description missing
 *
 * \author Christoph Lohmann
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixBandedTest
  : public UnitTest
{
  typedef SparseMatrixBanded<DT_, IT_> MatrixType;

public:
  explicit SparseMatrixBandedTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBandedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    // create random matrix
    const Index tsize(100);
    const Index rows(tsize + rng(Index(0), Index(20)));
    const Index columns(tsize + rng(Index(0), Index(20)));

    const Index num_bands(5 + rng(Index(0), Index(10)));

    DenseVector<IT_, IT_> vec_offsets(num_bands);
    DenseVector<DT_, IT_> vec_val(rng, num_bands * rows, DT_(0.1), DT_(1.0));

    // create random vector of offsets
    FEAT::Adjacency::Permutation permutation(rows + columns - 1, rng);
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      for (Index i(0); i < num_bands; ++i)
      {
        voff[i] = IT_(permutation.get_perm_pos()[i]);
      }
      std::sort(voff.get_w(), voff.get_w() + num_bands);
    }

    SparseMatrixBanded<DT_, IT_> a(rows, columns, vec_val, vec_offsets);

    auto b = a.clone();
    SparseMatrixBanded<DT_, IT_> bl(b.layout());
    TEST_CHECK_EQUAL(bl.num_nzes(), b.num_nzes());
    TEST_CHECK_EQUAL(bl.num_rows(), b.num_rows());
    TEST_CHECK_EQUAL(bl.num_cols(), b.num_cols());

    bl = b.layout();
    TEST_CHECK_EQUAL(bl.num_nzes(), b.num_nzes());
    TEST_CHECK_EQUAL(bl.num_rows(), b.num_rows());
    TEST_CHECK_EQUAL(bl.num_cols(), b.num_cols());

    typename SparseLayout<IT_, SparseLayoutId::lt_banded>::template MatrixType<DT_> x(b.layout());
    // icc 14.0.2 does not understand the following line, so we need a typedef here
    //typename decltype(b.layout())::template MatrixType<DT_> y(b.layout());
    typedef decltype(b.layout()) LayoutId;
    typename LayoutId::template MatrixType<DT_> y(b.layout());
    TEST_CHECK_EQUAL(x.offsets_arbiter(), b.offsets_arbiter());
    TEST_CHECK_NOT_EQUAL(x.val_arbiter(), b.val_arbiter());

    SparseMatrixBanded<DT_, IT_> c;
    c.clone(a, LAFEM::CloneMode::Weak);

    SparseMatrixBanded<DT_, IT_> z;
    z.convert(a);
    TEST_CHECK_LESS_THAN(a.max_rel_diff(z), tol);


    TEST_CHECK_NOT_EQUAL(c.val_arbiter(), a.val_arbiter());
    TEST_CHECK_EQUAL(c.offsets_arbiter(), a.offsets_arbiter());
    c = z.clone(CloneMode::Deep);
    TEST_CHECK_NOT_EQUAL(c.val_arbiter(), z.val_arbiter());
    TEST_CHECK_NOT_EQUAL(c.offsets_arbiter(), z.offsets_arbiter());

    DenseVector<IT_, IT_> offsets_d(c.num_bands(), c.offsets_arbiter().attach());
    DenseVector<DT_, IT_> val_d(c.num_bands() * c.num_rows(), c.val_arbiter().attach());
    SparseMatrixBanded<DT_, IT_> d(c.num_rows(), c.num_cols(), val_d, offsets_d);
    TEST_CHECK_LESS_THAN(d.max_rel_diff(c), tol);

    SparseMatrixBanded<DT_, IT_> e;
    e.convert(c);
    TEST_CHECK_LESS_THAN(e.max_rel_diff(c), tol);
    e.copy(c);
    TEST_CHECK_LESS_THAN(e.max_rel_diff(c), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBandedConvertTest
  : public UnitTest
{
public:
  explicit SparseMatrixBandedConvertTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBandedConvertTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    // create random matrix
    const Index tsize(10);
    const Index rows(tsize + rng(Index(0), Index(5)));
    const Index columns(tsize + rng(Index(0), Index(5)));

    const Index num_bands(5 + rng(Index(0), Index(3)));

    DenseVector<IT_, IT_> vec_offsets(num_bands);
    DenseVector<DT_, IT_> vec_val(rng, num_bands * rows, DT_(0.1), DT_(1.0));

    // create random vector of offsets
    FEAT::Adjacency::Permutation permutation(rows + columns - 1, rng);
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      for (Index i(0); i < num_bands; ++i)
      {
        voff[i] = IT_(permutation.get_perm_pos()[i]);
      }
      std::sort(voff.get_w(), voff.get_w() + num_bands);
    }

    SparseMatrixBanded<DT_, IT_> a(rows, columns, vec_val, vec_offsets);

    SparseMatrixCSR<DT_, IT_> csr;
    csr.convert(a);
    TEST_CHECK_EQUAL(csr.num_rows(), a.num_rows());
    TEST_CHECK_EQUAL(csr.num_cols(), a.num_cols());
    TEST_CHECK_EQUAL(csr.num_nzes(), a.num_nzes());

    SparseMatrixBanded<DT_, IT_> b;
    b.convert(csr);
    TEST_CHECK_EQUAL(csr.num_rows(), b.num_rows());
    TEST_CHECK_EQUAL(csr.num_cols(), b.num_cols());
    TEST_CHECK_EQUAL(csr.num_nzes(), b.num_nzes());

    TEST_CHECK(a.same_layout(b));

    // max_rel_diff will fail because the non-zero entries outside of the actual matrix window
    // are undefined, so we test mat-vec-mult with some random vectors instead
    //TEST_CHECK_LESS_THAN(a.max_rel_diff(b), tol);
    for(int k = 0; k < 10; ++k)
    {
      DenseVector<DT_, IT_> y1(a.num_rows()), y2(a.num_rows()), y3(a.num_rows());
      DenseVector<DT_, IT_> x(rng, a.num_cols(), DT_(0.1), DT_(5.0));

      a.apply(y1, x);
      b.apply(y2, x);
      csr.apply(y3, x);
      TEST_CHECK_LESS_THAN(y1.max_rel_diff(y2), tol);
      TEST_CHECK_LESS_THAN(y2.max_rel_diff(y3), tol);
      TEST_CHECK_LESS_THAN(y3.max_rel_diff(y1), tol);
    }
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedConvertTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBandedSerializeTest
  : public UnitTest
{
  typedef SparseMatrixBanded<DT_, IT_> MatrixType;

public:
  explicit SparseMatrixBandedSerializeTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBandedSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    // create random matrix
    const Index tsize(100);
    const Index rows(tsize + rng(Index(0), Index(20)));
    const Index columns(tsize + rng(Index(0), Index(20)));

    const Index num_bands(5 + rng(Index(0), Index(10)));

    DenseVector<IT_, IT_> vec_offsets(num_bands);
    DenseVector<DT_, IT_> vec_val(rng, num_bands * rows, DT_(0.1), DT_(1.0));

    // create random vector of offsets
    FEAT::Adjacency::Permutation permutation(rows + columns - 1, rng);
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      for (Index i(0); i < num_bands; ++i)
      {
        voff[i] = IT_(permutation.get_perm_pos()[i]);
      }
      std::sort(voff.get_w(), voff.get_w() + num_bands);
    }

    SparseMatrixBanded<DT_, IT_> c(rows, columns, vec_val, vec_offsets);

    BinaryStream bs;
    c.write_out(FileMode::fm_bm, bs);
    bs.seekg(0);
    SparseMatrixBanded<DT_, IT_> f(FileMode::fm_bm, bs);
    TEST_CHECK_LESS_THAN(f.max_rel_diff(c), tol);

    auto kp = c.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixBanded<DT_, IT_> k(kp);
    TEST_CHECK_LESS_THAN(k.max_rel_diff(c), tol);

#ifdef FEAT_HAVE_ZLIB
    auto zl = c.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixBanded<DT_, IT_> zlib(zl);
    TEST_CHECK_LESS_THAN(zlib.max_rel_diff(c), tol);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = c.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixBanded<DT_, IT_> zfp(zf);
    TEST_CHECK_LESS_THAN(zfp.max_rel_diff(c), DT_(1E-7));
#endif
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBandedApplyTest
  : public UnitTest
{
  typedef SparseMatrixBanded<DT_, IT_> MatrixType;

private:
  const Index _opt;

public:
  explicit SparseMatrixBandedApplyTest(PreferredBackend backend) :
    UnitTest("SparseMatrixBandedApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend), _opt(0)
  {
  }

  explicit SparseMatrixBandedApplyTest(PreferredBackend backend, const Index opt) :
    UnitTest("SparseMatrixBandedApplyTest: " + stringify(opt) + " offsets", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend), _opt(opt)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    DT_ alpha(DT_(4.321));

    // create random matrix
    const Index tsize(100);
    const Index rows(tsize + rng(Index(0), Index(20)));
    const Index columns(tsize + rng(Index(0), Index(20)));

    Index num_bands;
    if (_opt == 0)
      num_bands = 5 + rng(Index(0), Index(10));
    else
      num_bands = _opt;

    DenseVector<IT_, IT_> vec_offsets(num_bands);
    DenseVector<DT_, IT_> vec_val(rng, num_bands * rows, DT_(0), DT_(10));

    // create random vector of offsets
    FEAT::Adjacency::Permutation permutation(rows + columns - 1, rng);
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      for (Index i(0); i < num_bands; ++i)
      {
        voff[i] = IT_(permutation.get_perm_pos()[i]);
      }
      std::sort(voff.get_w(), voff.get_w() + num_bands);
    }

    // create test-matrix
    SparseMatrixBanded<DT_, IT_> sys(rows, columns, vec_val, vec_offsets);

    DenseMatrix<DT_, IT_> dense_mat(rows, columns);
    dense_mat.format();
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_r();
      Memory::TypedView<DT_> vval = vec_val.elements_view_r();
      Memory::TypedView<DT_> vdns = dense_mat.elements_view_rw();
      for (Index off_idx(0); off_idx < num_bands; ++off_idx)
      {
        Index offset = Index(voff(off_idx));
        Index start_col = offset >= rows ? offset + 1 - rows : 0;
        Index end_col = offset < columns ? offset + 1 : columns;
        for (Index col(start_col); col  < end_col; ++col)
        {
          Index row = col + rows - offset  - 1;
          vdns[row * columns +  col] = vval(off_idx * rows + row);
        }
      }
    }

    auto x(sys.create_vector_r());
    DenseVector<DT_, IT_> r(rows, DT_(0));
    auto ref1(sys.create_vector_l());
    auto ref2(sys.create_vector_l());

    x.format(rng, DT_(0.1), DT_(2.1));
    ref1.format();

    sys.apply(r, x);
    dense_mat.apply(ref1, x);

    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref1), tol);

    {
      Memory::TypedView<DT_> vref = ref1.elements_view_rw();
      for (Index i(0); i < r.size(); ++i)
        vref[i] += Math::cos(DT_(i)) + DT_(1);
    }

    sys.apply(r, x, ref1, alpha);
    dense_mat.apply(ref2, x, ref1, alpha);

    TEST_CHECK_LESS_THAN(r.max_rel_diff(ref2), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, double, std::uint32_t, PreferredBackend::generic);

SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, float, std::uint64_t, PreferredBackend::generic, 9);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, double, std::uint64_t, PreferredBackend::generic, 9);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, float, std::uint32_t, PreferredBackend::generic, 9);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, double, std::uint32_t, PreferredBackend::generic, 9);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedApplyTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBandedScaleTest
  : public UnitTest
{
  typedef SparseMatrixBanded<DT_, IT_> MatrixType;

public:
  explicit SparseMatrixBandedScaleTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBandedScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    DT_ alpha(DT_(4.321));

    // create random matrix
    const Index tsize(100);
    const Index rows(tsize + rng(Index(0), Index(20)));
    const Index columns(tsize + rng(Index(0), Index(20)));

    const Index num_bands(5 + rng(Index(0), Index(10)));

    DenseVector<IT_, IT_> vec_offsets(num_bands);
    DenseVector<DT_, IT_> vec_val_a(num_bands * rows);
    DenseVector<DT_, IT_> vec_val_ref(num_bands * rows, DT_(1));

    // create random vector of offsets
    FEAT::Adjacency::Permutation permutation(rows + columns - 1, rng);
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      for (Index i(0); i < num_bands; ++i)
      {
        voff[i] = IT_(permutation.get_perm_pos()[i]);
      }
      std::sort(voff.get_w(), voff.get_w() + num_bands);
    }

    // fill data-array
    {
      Memory::TypedView<DT_> va = vec_val_a.elements_view_w();
      Memory::TypedView<DT_> vr = vec_val_ref.elements_view_w();
      for (Index i(0); i < vec_val_a.size(); ++i)
      {
        va[i] = rng(DT_(0), DT_(10));
        vr[i] = va[i] * alpha;
      }
    }

    // create test-matrix
    SparseMatrixBanded<DT_, IT_> ref(rows, columns, vec_val_ref, vec_offsets);
    SparseMatrixBanded<DT_, IT_> a(rows, columns, vec_val_a, vec_offsets);
    SparseMatrixBanded<DT_, IT_> b;
    b.clone(a);

    b.scale(a, alpha);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(ref), tol);

    a.scale(a, alpha);
    TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedScaleTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBandedAxpyTest
  : public UnitTest
{
  typedef SparseMatrixBanded<DT_, IT_> MatrixType;

public:
  explicit SparseMatrixBandedAxpyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBandedAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    DT_ alpha(DT_(4.321));

    // create random matrix
    const Index tsize(100);
    const Index rows(tsize + rng(Index(0), Index(20)));
    const Index columns(tsize + rng(Index(0), Index(20)));

    const Index num_bands(5 + rng(Index(0), Index(10)));

    DenseVector<IT_, IT_> vec_offsets(num_bands);
    DenseVector<DT_, IT_> vec_val_a(num_bands * rows);
    DenseVector<DT_, IT_> vec_val_b(num_bands * rows);
    DenseVector<DT_, IT_> vec_val_ref(num_bands * rows, DT_(1));
    DenseVector<DT_, IT_> vec_val_ref2(num_bands * rows, DT_(1));

    // create random vector of offsets
    FEAT::Adjacency::Permutation permutation(rows + columns - 1, rng);
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      for (Index i(0); i < num_bands; ++i)
      {
        voff[i] = IT_(permutation.get_perm_pos()[i]);
      }
      std::sort(voff.get_w(), voff.get_w() + num_bands);
    }

    // fill data-array
    {
      Memory::TypedView<DT_> va = vec_val_a.elements_view_w();
      Memory::TypedView<DT_> vb = vec_val_b.elements_view_w();
      Memory::TypedView<DT_> vr1 = vec_val_ref.elements_view_w();
      Memory::TypedView<DT_> vr2 = vec_val_ref2.elements_view_w();
      for (Index i(0); i < vec_val_a.size(); ++i)
      {
        va[i] = rng(DT_(0), DT_(10));
        vb[i] = rng(DT_(0), DT_(10));
        vr1[i] = va(i) + vb(i) * alpha;
        vr2[i] = vb(i) + vb(i) * alpha;
      }
    }

    // create test-matrix
    SparseMatrixBanded<DT_, IT_> ref(rows, columns, vec_val_ref, vec_offsets);
    SparseMatrixBanded<DT_, IT_> ref2(rows, columns, vec_val_ref2, vec_offsets);
    SparseMatrixBanded<DT_, IT_> a(rows, columns, vec_val_a, vec_offsets);
    SparseMatrixBanded<DT_, IT_> b(rows, columns, vec_val_b, vec_offsets);

    // r != x
    a.axpy(b, alpha);
    TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), tol);

    // r == x
    b.axpy(b, alpha);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(ref2), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedAxpyTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBandedMaxRelDiffTest
  : public UnitTest
{
  typedef SparseMatrixBanded<DT_, IT_> MatrixType;

public:
  explicit SparseMatrixBandedMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBandedMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    const DT_ delta = DT_(123.5);
    const DT_ initial_value = DT_(10.0);

    const Index rows = 50;
    const Index columns = 50;
    const Index num_bands = 3;

    DenseVector<IT_, IT_> vec_offsets(num_bands);
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      voff[0] = IT_(rows - 1);
      voff[1] = IT_(rows);
      voff[2] = IT_(rows + 1);
    }

    const Index vec_val_size = num_bands * rows;
    DenseVector<DT_, IT_> vec_val_b(vec_val_size, initial_value);

    // create reference matrix b
    MatrixType b(rows, columns, vec_val_b, vec_offsets);

    // copy b into a
    MatrixType a = b.clone();

    // create delta matrix with only one value
    const Index diff_offset_idx = 1;
    const Index diff_row_idx = 20;
    const Index diff_val_index = (diff_offset_idx * rows) + diff_row_idx;
    DenseVector<DT_, IT_> vec_val_delta(vec_val_size, DT_(0.));
    vec_val_delta.elements_view_rw()[diff_val_index] = delta;
    MatrixType delta_mat(rows, columns, vec_val_delta, vec_offsets);

    // a = a + 1.0 * delta_mat
    a.axpy(delta_mat, DT_(1.0));

    // reference value
    const DT_ ref = delta / (DT_(2) * initial_value + delta + DT_(1));

    // test ||a-b||_infty
    const DT_ diff_1 = a.max_rel_diff(b);
    TEST_CHECK_RELATIVE(diff_1, ref, tol);

    // test ||b-a||_infty
    const DT_ diff_2 = b.max_rel_diff(a);
    TEST_CHECK_RELATIVE(diff_2, ref, tol);
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedMaxRelDiffTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixBandedSameLayoutTest
  : public UnitTest
{
public:
  explicit SparseMatrixBandedSameLayoutTest(PreferredBackend backend)
    : UnitTest("SparseMatrixBandedSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  typedef SparseMatrixBanded<DT_, IT_> MatrixType;

  virtual void run() const override
  {
    const DT_ initial_value = DT_(10.0);
    const Index rows = 50;
    const Index columns = 50;
    const Index num_bands = 3;

    DenseVector<IT_, IT_> vec_offsets(num_bands);
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      voff[0] = IT_(rows - 1);
      voff[1] = IT_(rows);
      voff[2] = IT_(rows + 1);
    }

    const Index vec_val_size = num_bands * rows;
    DenseVector<DT_, IT_> vec_val_a(vec_val_size, initial_value);

    // create reference matrix a
    MatrixType a(rows, columns, vec_val_a, vec_offsets);

    // weak copy
    auto b = a.clone(CloneMode::Weak);
    TEST_CHECK(a.same_layout(b));

    // shallow copy
    auto c = a.clone(CloneMode::Shallow);
    TEST_CHECK(a.same_layout(c));

    // different values at same position
    DenseVector<DT_, IT_> vec_val_d(vec_val_size, DT_(0.5));
    MatrixType d(rows, columns, vec_val_d, vec_offsets);
    TEST_CHECK(a.same_layout(d));

    // values at different position
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      voff[0] = IT_(rows + 1);
      voff[1] = IT_(rows + 2);
      voff[2] = IT_(rows + 3);
    }
    MatrixType e(rows, columns, vec_val_a, vec_offsets);
    TEST_CHECK(!a.same_layout(e));

    // different sizes
    {
      Memory::TypedView<IT_> voff = vec_offsets.elements_view_w();
      voff[0] = IT_(rows - 1);
      voff[1] = IT_(rows);
      voff[2] = IT_(rows + 1);
    }
    DenseVector<DT_, IT_> vec_val_f(num_bands * (rows - 2), initial_value);
    MatrixType f(rows - 2, columns - 2, vec_val_f, vec_offsets);
    TEST_CHECK(!a.same_layout(f));
  }
};
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixBandedSameLayoutTest, double, std::uint32_t, PreferredBackend::cuda);
#endif
