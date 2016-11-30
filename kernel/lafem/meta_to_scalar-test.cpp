#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/power_diag_matrix.hpp>
#include <kernel/lafem/power_col_matrix.hpp>
#include <kernel/lafem/power_row_matrix.hpp>
#include <kernel/lafem/power_full_matrix.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>
#include <kernel/util/random.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;


/**
 * \brief MetaToScalarTest
 *
 * This class defines the MetaToScalarTest which transforms a meta-matrix to a scalar matrix
 * It depends on the class "MetaMatrixTestBase"
 *
 * \author Christoph Lohmann
 */
template<typename MT_>
class MetaToScalarTest
  : public MetaMatrixTestBase<typename MT_::MemType, typename MT_::DataType, typename MT_::IndexType>
{
public:
  typedef typename MT_::DataType DataType;
  typedef typename MT_::IndexType IndexType;
  typedef MetaMatrixTestBase<typename MT_::MemType, DataType, IndexType> BaseClass;
  typedef typename BaseClass::SystemDiagMatrix SystemMatrix;
  typedef typename BaseClass::SystemVector SystemVector;

  typedef typename SystemVector::MemType Mem_;
  typedef DataType DT_;
  typedef IndexType IT_;

   MetaToScalarTest()
    : BaseClass("meta_to_scalar_test: " + MT_::name())
  {
  }

  virtual ~MetaToScalarTest()
  {
  }

  virtual void run() const override
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.7)));

    // generate a test system: A,x,b
    SystemMatrix mat_sys;
    SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_sys, vec_sol, vec_rhs);

    // test t <- A*x; t <- t - b
    vec_sol.scale(vec_sol, DT_(0));

    MT_ mat_sys_scalar;
    mat_sys_scalar.convert(mat_sys);

    DenseVector<Mem_, DT_, IT_> vec_rhs_scalar(mat_sys_scalar.rows());
    DenseVector<Mem_, DT_, IT_> vec_rhs_scalar2(mat_sys_scalar.rows());
    DenseVector<Mem_, DT_, IT_> vec_sol_scalar(mat_sys_scalar.columns());
    vec_sol_scalar.format();



    for (Index i(0); i < vec_sol.size(); ++i)
    {
      vec_sol_scalar(i, DT_(1));
      vec_sol_scalar.copy_inv(vec_sol);
      mat_sys.apply(vec_rhs, vec_sol);
      vec_sol_scalar.convert(vec_sol);
      vec_rhs_scalar2.copy(vec_rhs);
      mat_sys_scalar.apply(vec_rhs_scalar, vec_sol_scalar);

      for (Index j(0); j < vec_rhs_scalar.size(); ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j), vec_rhs_scalar2(j), tol);
      }

      vec_sol_scalar(i, DT_(0));
    }

    /**
     * check VecMetaToScalar and "VecScalarToMeta"
     */
    vec_rhs.scale(vec_rhs, DT_(0));

    for (Index j(0); j < vec_rhs_scalar.size(); ++j)
    {
      vec_rhs_scalar(j, DT_(j));
    }

    vec_rhs_scalar.copy_inv(vec_rhs);

    for (Index j(0); j < vec_rhs.first().first().size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j), vec_rhs.first().first()(j), tol);
    }
    for (Index j(0); j < vec_rhs.first().rest().size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j + vec_rhs.first().first().size()), vec_rhs.first().rest()(j), tol);
    }
    for (Index j(0); j < vec_rhs.rest().size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j + vec_rhs.first().size()), vec_rhs.rest()(j), tol);
    }

    vec_rhs_scalar.scale(vec_rhs_scalar, DT_(0));

    vec_rhs_scalar.copy(vec_rhs);

    for (Index j(0); j < vec_rhs.first().first().size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j), vec_rhs.first().first()(j), tol);
    }
    for (Index j(0); j < vec_rhs.first().rest().size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j + vec_rhs.first().first().size()), vec_rhs.first().rest()(j), tol);
    }
    for (Index j(0); j < vec_rhs.rest().size(); ++j)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_rhs_scalar(j + vec_rhs.first().size()), vec_rhs.rest()(j), tol);
    }
  }
};

MetaToScalarTest<SparseMatrixCOO<Mem::Main, float, unsigned int> > cpu_meta_matrix_to_coo_test_generic_float_uint;
MetaToScalarTest<SparseMatrixCOO<Mem::Main, double, unsigned int> > cpu_meta_matrix_to_coo_test_generic_double_uint;
MetaToScalarTest<SparseMatrixCOO<Mem::Main, float, unsigned long> > cpu_meta_matrix_to_coo_test_generic_float_ulong;
MetaToScalarTest<SparseMatrixCOO<Mem::Main, double, unsigned long> > cpu_meta_matrix_to_coo_test_generic_double_ulong;

MetaToScalarTest<SparseMatrixCSR<Mem::Main, float, unsigned int> > cpu_meta_matrix_to_csr_test_generic_float_uint;
MetaToScalarTest<SparseMatrixCSR<Mem::Main, double, unsigned int> > cpu_meta_matrix_to_csr_test_generic_double_uint;
MetaToScalarTest<SparseMatrixCSR<Mem::Main, float, unsigned long> > cpu_meta_matrix_to_csr_test_generic_float_ulong;
MetaToScalarTest<SparseMatrixCSR<Mem::Main, double, unsigned long> > cpu_meta_matrix_to_csr_test_generic_double_ulong;

MetaToScalarTest<SparseMatrixELL<Mem::Main, float, unsigned int> > cpu_meta_matrix_to_ell_test_generic_float_uint;
MetaToScalarTest<SparseMatrixELL<Mem::Main, double, unsigned int> > cpu_meta_matrix_to_ell_test_generic_double_uint;
MetaToScalarTest<SparseMatrixELL<Mem::Main, float, unsigned long> > cpu_meta_matrix_to_ell_test_generic_float_ulong;
MetaToScalarTest<SparseMatrixELL<Mem::Main, double, unsigned long> > cpu_meta_matrix_to_ell_test_generic_double_ulong;

#ifdef FEAT_HAVE_QUADMATH
MetaToScalarTest<SparseMatrixCOO<Mem::Main, __float128, unsigned int> > cpu_meta_matrix_to_coo_test_generic_float128_uint;
MetaToScalarTest<SparseMatrixCOO<Mem::Main, __float128, unsigned long> > cpu_meta_matrix_to_coo_test_generic_float128_ulong;

MetaToScalarTest<SparseMatrixCSR<Mem::Main, __float128, unsigned int> > cpu_meta_matrix_to_csr_test_generic_float128_uint;
MetaToScalarTest<SparseMatrixCSR<Mem::Main, __float128, unsigned long> > cpu_meta_matrix_to_csr_test_generic_float128_ulong;

MetaToScalarTest<SparseMatrixELL<Mem::Main, __float128, unsigned int> > cpu_meta_matrix_to_ell_test_generic_float128_uint;
MetaToScalarTest<SparseMatrixELL<Mem::Main, __float128, unsigned long> > cpu_meta_matrix_to_ell_test_generic_float128_ulong;
#endif

#ifdef FEAT_HAVE_CUDA
// MetaToScalarTest<SparseMatrixCOO<Mem::CUDA, float, unsigned int> > cuda_meta_matrix_to_coo_test_generic_float_uint;
// MetaToScalarTest<SparseMatrixCOO<Mem::CUDA, double, unsigned int> > cuda_meta_matrix_to_coo_test_generic_double_uint;
// MetaToScalarTest<SparseMatrixCOO<Mem::CUDA, float, unsigned long> > cuda_meta_matrix_to_coo_test_generic_float_ulong;
// MetaToScalarTest<SparseMatrixCOO<Mem::CUDA, double, unsigned long> > cuda_meta_matrix_to_coo_test_generic_double_ulong;

MetaToScalarTest<SparseMatrixCSR<Mem::CUDA, float, unsigned int> > cuda_meta_matrix_to_csr_test_generic_float_uint;
MetaToScalarTest<SparseMatrixCSR<Mem::CUDA, double, unsigned int> > cuda_meta_matrix_to_csr_test_generic_double_uint;
MetaToScalarTest<SparseMatrixCSR<Mem::CUDA, float, unsigned long> > cuda_meta_matrix_to_csr_test_generic_float_ulong;
MetaToScalarTest<SparseMatrixCSR<Mem::CUDA, double, unsigned long> > cuda_meta_matrix_to_csr_test_generic_double_ulong;

MetaToScalarTest<SparseMatrixELL<Mem::CUDA, float, unsigned int> > cuda_meta_matrix_to_ell_test_generic_float_uint;
MetaToScalarTest<SparseMatrixELL<Mem::CUDA, double, unsigned int> > cuda_meta_matrix_to_ell_test_generic_double_uint;
MetaToScalarTest<SparseMatrixELL<Mem::CUDA, float, unsigned long> > cuda_meta_matrix_to_ell_test_generic_float_ulong;
MetaToScalarTest<SparseMatrixELL<Mem::CUDA, double, unsigned long> > cuda_meta_matrix_to_ell_test_generic_double_ulong;
#endif


/**
 * \brief MetaBCSRToScalarTest
 *
 * This class defines the MetaBCSRToScalarTest which transforms a meta-matrix with a block-csr-matrix to a scalar matrix
 *
 * \author Christoph Lohmann
 */
template<typename MT_>
class MetaBCSRToScalarTest
  : public FEAT::TestSystem::FullTaggedTest<typename MT_::MemType, typename MT_::DataType, typename MT_::IndexType>
{
public:
  typedef typename MT_::DataType DataType;
  typedef typename MT_::IndexType IndexType;
  typedef typename MT_::MemType Mem_;
  typedef DataType DT_;
  typedef IndexType IT_;

   MetaBCSRToScalarTest()
    : FEAT::TestSystem::FullTaggedTest<typename MT_::MemType, typename MT_::DataType, typename MT_::IndexType>
      ("meta_block_csr_to_scalar_test: " + MT_::name())
  {
  }

  virtual ~MetaBCSRToScalarTest()
  {
  }

  virtual void run() const override
  {
    Random::SeedType seed(Random::SeedType(time(nullptr)));
    Random random(seed);
    std::cout << "seed: " << seed << std::endl;

    DenseVector<Mem_, DT_, IT_> dv11(18);
    for (Index i(0) ; i < dv11.size() ; ++i)
      dv11(i, random(DT_(0), DT_(10)));
    DenseVector<Mem_, IT_, IT_> dv12(3);
    dv12(0, IT_(0));
    dv12(1, IT_(1));
    dv12(2, IT_(2));
    DenseVector<Mem_, IT_, IT_> dv13(3);
    dv13(0, IT_(0));
    dv13(1, IT_(1));
    dv13(2, IT_(3));
    SparseMatrixBCSR<Mem_, DT_, IT_, 2, 3> c1(2, 3, dv12, dv11, dv13);

    DenseVector<Mem_, DT_, IT_> dv21(24);
    for (Index i(0) ; i < dv21.size() ; ++i)
      dv21(i, random(DT_(0), DT_(10)));
    DenseVector<Mem_, IT_, IT_> dv22(4);
    dv22(0, IT_(0));
    dv22(1, IT_(1));
    dv22(2, IT_(1));
    dv22(3, IT_(2));
    DenseVector<Mem_, IT_, IT_> dv23(4);
    dv23(0, IT_(0));
    dv23(1, IT_(1));
    dv23(2, IT_(2));
    dv23(3, IT_(4));
    SparseMatrixBCSR<Mem_, DT_, IT_, 2, 3> c2(3, 3, dv22, dv21, dv23);

    typedef SparseMatrixBCSR<Mem_, DT_, IT_, 2, 3> BCSRMatrix;
    typedef PowerColMatrix<BCSRMatrix, 2> ColMatrix;
    typedef PowerRowMatrix<BCSRMatrix, 2> RowMatrix;
    typedef SaddlePointMatrix<BCSRMatrix, RowMatrix, ColMatrix> SaddleMatrix;

    typedef SparseMatrixCSR<Mem_, DT_, IT_> CSRMatrix;
    typedef PowerColMatrix<CSRMatrix, 2> ColMatrix2;
    typedef PowerRowMatrix<CSRMatrix, 2> RowMatrix2;
    typedef SaddlePointMatrix<CSRMatrix, RowMatrix2, ColMatrix2> SaddleMatrix2;

    SaddleMatrix c4;
    c4.template at<0,0>().convert(c2);
    c4.template at<0,1>().template at<0,0>().convert(c2);
    c4.template at<0,1>().template at<0,1>().convert(c2);
    c4.template at<1,0>().template at<0,0>().convert(c1);
    c4.template at<1,0>().template at<1,0>().convert(c1);

    SaddleMatrix2 c5;
    c5.convert(c4);

    MT_ c6;
    c6.convert(c5);

    MT_ c7;
    c7.convert(c4);

    TEST_CHECK_EQUAL(c6, c7);
  }
};

MetaBCSRToScalarTest<SparseMatrixCOO<Mem::Main, float, unsigned int> > cpu_meta_bcsr_matrix_to_coo_test_generic_float_uint;
MetaBCSRToScalarTest<SparseMatrixCOO<Mem::Main, double, unsigned int> > cpu_meta_bcsr_matrix_to_coo_test_generic_double_uint;
MetaBCSRToScalarTest<SparseMatrixCOO<Mem::Main, float, unsigned long> > cpu_meta_bcsr_matrix_to_coo_test_generic_float_ulong;
MetaBCSRToScalarTest<SparseMatrixCOO<Mem::Main, double, unsigned long> > cpu_meta_bcsr_matrix_to_coo_test_generic_double_ulong;

MetaBCSRToScalarTest<SparseMatrixCSR<Mem::Main, float, unsigned int> > cpu_meta_bcsr_matrix_to_csr_test_generic_float_uint;
MetaBCSRToScalarTest<SparseMatrixCSR<Mem::Main, double, unsigned int> > cpu_meta_bcsr_matrix_to_csr_test_generic_double_uint;
MetaBCSRToScalarTest<SparseMatrixCSR<Mem::Main, float, unsigned long> > cpu_meta_bcsr_matrix_to_csr_test_generic_float_ulong;
MetaBCSRToScalarTest<SparseMatrixCSR<Mem::Main, double, unsigned long> > cpu_meta_bcsr_matrix_to_csr_test_generic_double_ulong;

MetaBCSRToScalarTest<SparseMatrixELL<Mem::Main, float, unsigned int> > cpu_meta_bcsr_matrix_to_ell_test_generic_float_uint;
MetaBCSRToScalarTest<SparseMatrixELL<Mem::Main, double, unsigned int> > cpu_meta_bcsr_matrix_to_ell_test_generic_double_uint;
MetaBCSRToScalarTest<SparseMatrixELL<Mem::Main, float, unsigned long> > cpu_meta_bcsr_matrix_to_ell_test_generic_float_ulong;
MetaBCSRToScalarTest<SparseMatrixELL<Mem::Main, double, unsigned long> > cpu_meta_bcsr_matrix_to_ell_test_generic_double_ulong;

#ifdef FEAT_HAVE_QUADMATH
MetaBCSRToScalarTest<SparseMatrixCOO<Mem::Main, __float128, unsigned int> > cpu_meta_bcsr_matrix_to_coo_test_generic_float128_uint;
MetaBCSRToScalarTest<SparseMatrixCOO<Mem::Main, __float128, unsigned long> > cpu_meta_bcsr_matrix_to_coo_test_generic_float128_ulong;

MetaBCSRToScalarTest<SparseMatrixCSR<Mem::Main, __float128, unsigned int> > cpu_meta_bcsr_matrix_to_csr_test_generic_float128_uint;
MetaBCSRToScalarTest<SparseMatrixCSR<Mem::Main, __float128, unsigned long> > cpu_meta_bcsr_matrix_to_csr_test_generic_float128_ulong;

MetaBCSRToScalarTest<SparseMatrixELL<Mem::Main, __float128, unsigned int> > cpu_meta_bcsr_matrix_to_ell_test_generic_float128_uint;
MetaBCSRToScalarTest<SparseMatrixELL<Mem::Main, __float128, unsigned long> > cpu_meta_bcsr_matrix_to_ell_test_generic_float128_ulong;
#endif

#ifdef FEAT_HAVE_CUDA
// MetaBCSRToScalarTest<SparseMatrixCOO<Mem::CUDA, float, unsigned int> > cuda_meta_bcsr_matrix_to_coo_test_generic_float_uint;
// MetaBCSRToScalarTest<SparseMatrixCOO<Mem::CUDA, double, unsigned int> > cuda_meta_bcsr_matrix_to_coo_test_generic_double_uint;
// MetaBCSRToScalarTest<SparseMatrixCOO<Mem::CUDA, float, unsigned long> > cuda_meta_bcsr_matrix_to_coo_test_generic_float_ulong;
// MetaBCSRToScalarTest<SparseMatrixCOO<Mem::CUDA, double, unsigned long> > cuda_meta_bcsr_matrix_to_coo_test_generic_double_ulong;

MetaBCSRToScalarTest<SparseMatrixCSR<Mem::CUDA, float, unsigned int> > cuda_meta_bcsr_matrix_to_csr_test_generic_float_uint;
MetaBCSRToScalarTest<SparseMatrixCSR<Mem::CUDA, double, unsigned int> > cuda_meta_bcsr_matrix_to_csr_test_generic_double_uint;
MetaBCSRToScalarTest<SparseMatrixCSR<Mem::CUDA, float, unsigned long> > cuda_meta_bcsr_matrix_to_csr_test_generic_float_ulong;
MetaBCSRToScalarTest<SparseMatrixCSR<Mem::CUDA, double, unsigned long> > cuda_meta_bcsr_matrix_to_csr_test_generic_double_ulong;

MetaBCSRToScalarTest<SparseMatrixELL<Mem::CUDA, float, unsigned int> > cuda_meta_bcsr_matrix_to_ell_test_generic_float_uint;
MetaBCSRToScalarTest<SparseMatrixELL<Mem::CUDA, double, unsigned int> > cuda_meta_bcsr_matrix_to_ell_test_generic_double_uint;
MetaBCSRToScalarTest<SparseMatrixELL<Mem::CUDA, float, unsigned long> > cuda_meta_bcsr_matrix_to_ell_test_generic_float_ulong;
MetaBCSRToScalarTest<SparseMatrixELL<Mem::CUDA, double, unsigned long> > cuda_meta_bcsr_matrix_to_ell_test_generic_double_ulong;
#endif



/**
 * \brief VecMetaToScalarTest
 *
 * This class defines the VecMetaToScalarTest which tests the conversion of
 * meta-vectors with a densevector with different template-parameters
 *
 * \author Christoph Lohmann
 */
template<
  typename Mem_,
  typename DT_,
  typename IT_>
class VecMetaToScalarTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  typedef DenseVector<Mem_, DT_, IT_> DenseVec;
  typedef DenseVectorBlocked<Mem_, DT_, IT_, 3> DenseVecBlocked;
  typedef PowerVector<DenseVec, 1> PowerVec1;
  typedef PowerVector<DenseVec, 2> PowerVec2;
  typedef PowerVector<DenseVec, 3> PowerVec3;
  typedef TupleVector<PowerVec3, PowerVec2, PowerVec1, DenseVec, DenseVecBlocked> TupleVec;

   VecMetaToScalarTest()
    : FullTaggedTest<Mem_, DT_, IT_>("vec_meta_to_scalar_test")
  {
  }

  virtual ~VecMetaToScalarTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps(Math::pow(Math::eps<DT_>(), DT_(0.8)));

    Random::SeedType seed(Random::SeedType(time(nullptr)));
    Random random(seed);
    std::cout << "seed: " << seed << std::endl;

    DenseVec dv[7];
    for (Index j(0); j < 7; ++j)
    {
      DenseVec tdv(random(Index(10), Index(15)));
      dv[j].convert(tdv);
      for (Index i(0); i < dv[j].size(); ++i)
      {
        dv[j](i, random(DT_(3.14), DT_(3456.2)));
      }
    }

    DenseVecBlocked dvb(2);
    Tiny::Vector<DT_, 3> tv1, tv2(0);
    tv1(0) = DT_(0);
    tv1(1) = DT_(3);
    tv1(2) = DT_(5);
    tv2(0) = DT_(0.1);
    tv2(1) = DT_(3);
    tv2(2) = DT_(-3.4);
    dvb(0, tv1);
    dvb(1, tv2);

    Index glob_size(0);
    for (Index j(0); j < 7; ++j)
    {
      glob_size += dv[j].size();
    }
    glob_size += dvb.template size<Perspective::native>();

    Index glob_size_pod(0);
    for (Index j(0); j < 7; ++j)
    {
      glob_size_pod += dv[j].size();
    }
    glob_size_pod += dvb.template size<Perspective::pod>();

    PowerVec3 pv3;
    pv3.template at<0>().convert(dv[0]);
    pv3.template at<1>().convert(dv[1]);
    pv3.template at<2>().convert(dv[2]);

    PowerVec2 pv2;
    pv2.template at<0>().convert(dv[3]);
    pv2.template at<1>().convert(dv[4]);

    PowerVec1 pv1;
    pv1.template at<0>().convert(dv[5]);

    TupleVec tv;
    tv.template at<0>().convert(pv3);
    tv.template at<1>().convert(pv2);
    tv.template at<2>().convert(pv1);
    tv.template at<3>().convert(dv[6]);
    tv.template at<4>().convert(dvb);

    typename TupleVec::template ContainerType<Mem::Main, DT_, IT_> tv_main;
    typename TupleVec::template ContainerType<Mem_, DT_, IT_> tv_convert;
    tv_main.convert(tv);

    TEST_CHECK_EQUAL(tv.name(), tv_convert.name());
    TEST_CHECK_EQUAL(tv_main.size(), glob_size);
    TEST_CHECK_EQUAL(tv_main.template size<Perspective::pod>(), glob_size_pod);

    for (Index i(0); i < tv.template at<2>().template at<0>().size(); ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(tv_main.template at<2>().template at<0>()(i), tv.template at<2>().template at<0>()(i), eps);
    }
    for (Index i(0); i < tv.template at<4>().size(); ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(tv_main.template at<4>()(i)(0), tv.template at<4>()(i)(0), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(tv_main.template at<4>()(i)(1), tv.template at<4>()(i)(1), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(tv_main.template at<4>()(i)(2), tv.template at<4>()(i)(2), eps);
    }

    DenseVector<Mem::Main, DT_, IT_> tv_scalar_main;
    tv_scalar_main.convert(tv);

    std::cout << tv.template at<1>().size() << std::endl;

    for (Index i(0); i < tv.template at<2>().template at<0>().size(); ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(tv_scalar_main(i + tv.template at<0>().template size<Perspective::pod>()
                                                 + tv.template at<1>().template size<Perspective::pod>()),
                                  tv.template at<2>().template at<0>()(i), eps);
    }
    for (Index i(0); i < tv.template at<4>().size(); ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(tv_scalar_main(Index(3) * i + Index(0) + tv.template at<0>().template size<Perspective::pod>()
                                                 + tv.template at<1>().template size<Perspective::pod>()
                                                 + tv.template at<2>().template size<Perspective::pod>()
                                                 + tv.template at<3>().template size<Perspective::pod>()),
                                  tv.template at<4>()(i)(0), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(tv_scalar_main(Index(3) * i + Index(1) + tv.template at<0>().template size<Perspective::pod>()
                                                 + tv.template at<1>().template size<Perspective::pod>()
                                                 + tv.template at<2>().template size<Perspective::pod>()
                                                 + tv.template at<3>().template size<Perspective::pod>()),
                                  tv.template at<4>()(i)(1), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(tv_scalar_main(Index(3) * i + Index(2) + tv.template at<0>().template size<Perspective::pod>()
                                                 + tv.template at<1>().template size<Perspective::pod>()
                                                 + tv.template at<2>().template size<Perspective::pod>()
                                                 + tv.template at<3>().template size<Perspective::pod>()),
                                  tv.template at<4>()(i)(2), eps);
    }
  }
};
VecMetaToScalarTest<Mem::Main, float, unsigned int> cpu_vec_meta_to_scalar_test_generic_float_uint;
VecMetaToScalarTest<Mem::Main, float, unsigned long> cpu_vec_meta_to_scalar_test_generic_float_ulong;
VecMetaToScalarTest<Mem::Main, double, unsigned int> cpu_vec_meta_to_scalar_test_generic_double_uint;
VecMetaToScalarTest<Mem::Main, double, unsigned long> cpu_vec_meta_to_scalar_test_generic_double_ulong;

#ifdef FEAT_HAVE_QUADMATH
VecMetaToScalarTest<Mem::Main, __float128, unsigned int> cpu_vec_meta_to_scalar_test_generic_float128_uint;
VecMetaToScalarTest<Mem::Main, __float128, unsigned long> cpu_vec_meta_to_scalar_test_generic_float128_ulong;
#endif

#ifdef FEAT_HAVE_CUDA
VecMetaToScalarTest<Mem::CUDA, float, unsigned int> cuda_vec_meta_to_scalar_test_generic_float_uint;
VecMetaToScalarTest<Mem::CUDA, float, unsigned long> cuda_vec_meta_to_scalar_test_generic_float_ulong;
VecMetaToScalarTest<Mem::CUDA, double, unsigned int> cuda_vec_meta_to_scalar_test_generic_double_uint;
VecMetaToScalarTest<Mem::CUDA, double, unsigned long> cuda_vec_meta_to_scalar_test_generic_double_ulong;
#endif
