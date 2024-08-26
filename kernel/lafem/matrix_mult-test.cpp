// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/util/random.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

template<
  typename DT_,
  typename IT_>
class MatrixMultTest
  : public UnitTest
{
  typedef SparseMatrixCSR<DT_, IT_> MatrixType;
  typedef DenseVector<DT_, IT_> VectorType;

public:
  MatrixMultTest(PreferredBackend backend)
    : UnitTest("MatrixMultTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~MatrixMultTest()
  {
  }

  static void mmult(MatrixType& x, const MatrixType& a, const MatrixType& y, const MatrixType& bt)
  {
    // This implementation is taken from an old (but well tested)
    // version of the MatrixMirror class template.

    const MatrixType row_mir_mat(a.clone());
    const MatrixType col_mir_mat(bt.transpose());

    VectorType vec_work = y.create_vector_l();
    vec_work.format();
    DT_* _work = vec_work.elements();

    // fetch row mirror arrays
    const IT_* row_ptr_a(row_mir_mat.row_ptr());
    const IT_* col_idx_a(row_mir_mat.col_ind());
    const DT_* av(row_mir_mat.val());

    // fetch col mirror arrays
    const IT_* row_ptr_b(col_mir_mat.row_ptr());
    const IT_* col_idx_b(col_mir_mat.col_ind());
    const DT_* bv(col_mir_mat.val());

    // fetch system matrix arrays
    const IT_* row_ptr_y(y.row_ptr());
    const IT_* col_idx_y(y.col_ind());
    const DT_* yv(y.val());

    // fetch buffer arrays
    const IT_* row_ptr_x(x.row_ptr());
    const IT_* col_idx_x(x.col_ind());
    DT_* xv(x.val());

    // In the following, we have to compute:
    //    X := A * Z := A * Y * B^T,
    // where:
    // X is the buffer matrix
    // Y is the system matrix
    // A is the row-mirror gather matrix
    // B is the col-mirror gather matrix

    // loop over all buffer rows (X)
    for(Index irow_x(0); irow_x < x.rows(); ++irow_x)
    {
      Index irow_a(irow_x); // row of a := row of x

      // loop over all non-zeroes in the buffer row (X_i.)
      for(IT_ ix(row_ptr_x[irow_x]); ix < row_ptr_x[irow_x + 1]; ++ix)
      {
        // init result
        DT_ x_ij(DT_(0));

        // fetch the column index
        IT_ irow_b(col_idx_x[ix]); // row of b := col of x

        // loop over all non-zeroes of the col-mirror (B_j.)
        for(IT_ ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
        {
          // and densify the sparse row B_j.
          _work[col_idx_b[ib]] = bv[ib];
        }

        // loop over all non-zeroes of the row-mirror (A_i.)
        for(IT_ ia(row_ptr_a[irow_a]); ia < row_ptr_a[irow_a + 1]; ++ia)
        {
          // fetch the column index
          IT_ irow_y(col_idx_a[ia]); // row of y := col of a

          // temporary entry: Z_kj := Y_k. * B_j.
          DT_ z_kj(DT_(0));

          // loop over all non-zeroes of the system matrix (Y_k.)
          for(IT_ iy(row_ptr_y[irow_y]); iy < row_ptr_y[irow_y + 1]; ++iy)
          {
            z_kj += DT_(yv[iy]) * DT_(_work[col_idx_y[iy]]);
          }

          // update x_ij += a_ik * z_kj
          x_ij += DT_(av[ia]) * z_kj;
        }

        // store X_ij
        xv[ix] = x_ij;

        // reset temporary data
        for(IT_ ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
        {
          _work[col_idx_b[ib]] = DT_(0);
        }
      }
    }
  }


  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));

    // create 3 matrices a,b and d
    LAFEM::PointstarFactoryFD<DT_, IT_> psf(17, 2);
    MatrixType a = psf.matrix_csr();
    MatrixType b = a.clone();
    MatrixType d = a.clone();

    // fill a, b and d with random values
    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rng(seed);
    for(Index i(0); i < a.used_elements(); ++i)
    {
      a.val()[i] = DT_(rng(0.0, 1.0));
      b.val()[i] = DT_(rng(0.0, 1.0));
      d.val()[i] = DT_(rng(0.0, 1.0));
    }

    // create matrix structure for X = D*A*B
    Adjacency::Graph graph_da(Adjacency::RenderType::injectify, d, a);
    Adjacency::Graph graph_dab(Adjacency::RenderType::injectify_sorted, graph_da, b);
    MatrixType x(graph_dab);

    // compute reference: X = D*A*B
    mmult(x, d, a, b);

    // subtract: X -= D*A*B
    x.add_double_mat_product(d, a, b, -DT_(1));

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(x.norm_frobenius(), DT_(0), tol);
  }
};

MatrixMultTest <float, std::uint32_t> matrix_mult_test_float_uint32(PreferredBackend::generic);
MatrixMultTest <double, std::uint32_t> matrix_mult_test_double_uint32(PreferredBackend::generic);
MatrixMultTest <float, std::uint64_t> matrix_mult_test_float_uint64(PreferredBackend::generic);
MatrixMultTest <double, std::uint64_t> matrix_mult_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MatrixMultTest <float, std::uint64_t> mkl_matrix_mult_test_float_uint64(PreferredBackend::mkl);
MatrixMultTest <double, std::uint64_t> mkl_matrix_mult_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MatrixMultTest <__float128, std::uint64_t> matrix_mult_test_float128_uint64(PreferredBackend::generic);
MatrixMultTest <__float128, std::uint32_t> matrix_mult_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MatrixMultTest <Half, std::uint32_t> matrix_mult_test_half_uint32(PreferredBackend::generic);
MatrixMultTest <Half, std::uint64_t> matrix_mult_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MatrixMultTest <float, std::uint64_t> cuda_matrix_mult_test_float_uint64(PreferredBackend::cuda);
MatrixMultTest <double, std::uint64_t> cuda_matrix_mult_test_double_uint64(PreferredBackend::cuda);
MatrixMultTest <float, std::uint32_t> cuda_matrix_mult_test_float_uint32(PreferredBackend::cuda);
MatrixMultTest <double, std::uint32_t> cuda_matrix_mult_test_double_uint32(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class MatrixMultTest2
  : public UnitTest
{
  typedef SparseMatrixCSR<DT_, IT_> MatrixType;
  typedef DenseVector<DT_, IT_> VectorType;

public:
  MatrixMultTest2(PreferredBackend backend)
    : UnitTest("MatrixMultTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~MatrixMultTest2()
  {
  }
  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create CSR matrix D of dimension (7x8)
    SparseMatrixFactory<DT_, IT_> factory_D(IT_(7), IT_(8));
    factory_D.add(IT_(0), IT_(5), DT_(2));
    factory_D.add(IT_(1), IT_(0), DT_(7));
    factory_D.add(IT_(2), IT_(3), DT_(5));
    factory_D.add(IT_(3), IT_(1), DT_(13));
    factory_D.add(IT_(4), IT_(0), DT_(4));
    factory_D.add(IT_(5), IT_(2), DT_(10));
    factory_D.add(IT_(5), IT_(6), DT_(3));
    factory_D.add(IT_(6), IT_(7), DT_(3));

    MatrixType d(factory_D.make_csr());


    //create CSR matrix B of dimension (8x3)
    SparseMatrixFactory<DT_, IT_> factory_B(IT_(8), IT_(3));
    factory_B.add(IT_(0), IT_(0), DT_(1));
    factory_B.add(IT_(1), IT_(2), DT_(1));
    factory_B.add(IT_(2), IT_(0), DT_(7));
    factory_B.add(IT_(3), IT_(1), DT_(2));
    factory_B.add(IT_(4), IT_(0), DT_(9));
    factory_B.add(IT_(6), IT_(1), DT_(5));
    factory_B.add(IT_(6), IT_(2), DT_(7));
    factory_B.add(IT_(7), IT_(0), DT_(3));
    MatrixType b(factory_B.make_csr());


   //create matrix adjacency structure for X=D*B
    Adjacency::Graph graph_db(Adjacency::RenderType::injectify_sorted, d, b);
    MatrixType x(graph_db);

    // calculate X=D*B
    x.add_mat_mat_product(d, b);

    //initialize reference Vectors
    const DT_ val_ref[8] = { DT_(7),DT_(10),DT_(13),DT_(4), DT_(70), DT_(15), DT_(21), DT_(9) };
    const IT_ col_ind_ref[8] = { IT_(0),IT_(1),IT_(2),IT_(0),IT_(0),IT_(1),IT_(2),IT_(0) };
    const IT_ row_ptr_ref[8] = { IT_(0),IT_(0),IT_(1),IT_(2),IT_(3),IT_(4),IT_(7),IT_(8) };


    // pointer on the arrays of the CSR matrix x
    IT_* row_ptr = x.row_ptr();
    IT_* col_ind = x.col_ind();
    DT_* val = x.val();

    // check if x is implemented correct
    for (IT_ i(0); i < 8; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(val_ref[i], val[i], tol);
      TEST_CHECK_EQUAL(col_ind_ref[i],col_ind[i]);
      TEST_CHECK_EQUAL(row_ptr_ref[i], row_ptr[i]);
    }


  }
};

MatrixMultTest2 <float, std::uint32_t> matrix_mult_test2_float_uint32(PreferredBackend::generic);
MatrixMultTest2 <double, std::uint32_t> matrix_mult_test2_double_uint32(PreferredBackend::generic);
MatrixMultTest2 <float, std::uint64_t> matrix_mult_test2_float_uint64(PreferredBackend::generic);
MatrixMultTest2 <double, std::uint64_t> matrix_mult_test2_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MatrixMultTest2 <float, std::uint64_t> mkl_matrix_mult_test2_float_uint64(PreferredBackend::mkl);
MatrixMultTest2 <double, std::uint64_t> mkl_matrix_mult_test2_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MatrixMultTest2 <__float128, std::uint64_t> matrix_mult_test2_float128_uint64(PreferredBackend::generic);
MatrixMultTest2 <__float128, std::uint32_t> matrix_mult_test2_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MatrixMultTest2 <Half, std::uint32_t> matrix_mult_test2_half_uint32(PreferredBackend::generic);
MatrixMultTest2 <Half, std::uint64_t> matrix_mult_test2_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MatrixMultTest2 <float, std::uint64_t> cuda_matrix_mult_test2_float_uint64(PreferredBackend::cuda);
MatrixMultTest2 <double, std::uint64_t> cuda_matrix_mult_test2_double_uint64(PreferredBackend::cuda);
MatrixMultTest2 <float, std::uint32_t> cuda_matrix_mult_test2_float_uint32(PreferredBackend::cuda);
MatrixMultTest2 <double, std::uint32_t> cuda_matrix_mult_test2_double_uint32(PreferredBackend::cuda);
#endif
