#include <kernel/util/linear_algebra.hpp>
#include <kernel/util/factorial.hpp>
#include <kernel/util/random.hpp>
#include <test_system/test_system.hpp>
#include <limits>
#include <stdint.h>

using namespace FEAST;
using namespace LinAlg;
using namespace TestSystem;

using std::abs;
using std::sin;
using std::cos;
using std::min;
using std::max;
using Intern::sqr;

// helper function: calculates the cube of a value
template<typename T_>
static inline T_ cub(T_ x)
{
  return x*x*x;
}

// helper function: calculates the binomial coefficient n over k
/*static inline u64 binomial(u64 n, u64 k)
{
  u64 a = 1;
  u64 b = 1;

  for(u64 i = k + 1 ; i <= n ; ++i)
    a *= i;
  for(u64 i = n - k ; i >= 2 ; --i)
    b *= i;

  return a / b;
}*/

/*
  This list keeps track of which linear algebra functions are tested by which test drivers.

  vec_clear                    test_vec_clear_copy
  vec_copy                     test_vec_clear_copy
  vec_swap                     --
  vec_shiftl                   --
  vec_shiftr                   --
  vec_scale                    --
  vec_dot                      test_vec_dot
  vec_axpy                     test_vec_axpy_lcomb
  vec_lcomb                    test_vec_axpy_lcomb
  vec_comp_mult                --
  vec_norm_asum                test_vec_norm_asum
  vec_norm_euclid              test_vec_norm_euclid
  vec_norm_max                 test_vec_norm_max, various others
  mat_clear                    --
  mat_identity                 test_mat_mat_mult
  mat_copy<false>              --
  mat_copy<true>               --
  mat_swap<false>              --
  mat_swap<true>               --
  mat_scale                    --
  mat_axpy<false>              test_mat_axpy, various others
  mat_axpy<true>               test_mat_plu, various others
  mat_mat_mult<false,false>    test_mat_mat_mult, test_mat_mat_solve
  mat_mat_mult<false,true>     test_mat_mat_mult, test_mat_mat_solve
  mat_mat_mult<true,false>     test_mat_mat_mult
  mat_mat_mult<true,true>      test_mat_mat_mult
  mat_vec_mult<false>          test_mat_vec_solve
  mat_vec_mult<true>           test_mat_vec_solve
  mat_norm_row_sum             --
  mat_norm_frobenius           --
  mat_norm_max                 test_mat_axpy, various others
  mat_factorise                test_mat_vec_solve
  mat_vec_solve<false>         test_mat_vec_solve
  mat_vec_solve<true>          test_mat_vec_solve
  mat_mat_solve<false>         test_mat_mat_solve
  mat_mat_solve<true>          test_mat_mat_solve
  mat_invert<n>                test_mat_inv_hilbert
*/

/**
 * \brief Tagged-Test for low-level Linear Algebra.
 *
 * \author Peter Zajac
 */
template<typename DataType_>
class LinAlgTest :
  public TaggedTest<Archs::None, DataType_>
{
#define ZERO DataType_(0)
#define ONE DataType_(1)
#define TWO DataType_(2)

  private:
    // machine exactness
    const DataType_ _eps;
    // random number generator
    mutable Random _random;

    // initialise Hilbert-Matrix; see http://en.wikipedia.org/wiki/Hilbert_matrix
    static void _init_hilbert_mat(size_t n, DataType_ a[])
    {
      for(size_t i(0) ; i < n ; ++i)
      {
        for(size_t j(0) ; j < n ; ++j)
        {
          a[i * n + j] = ONE / DataType_(i + j + 1);
        }
      }
    }

    // initialise Hilbert-Matrix inverse
    static void _init_hilbert_inv(size_t n, DataType_ a[])
    {
      for(size_t i(0) ; i < n ; ++i)
      {
        for(size_t j(0) ; j < n ; ++j)
        {
          a[i * n + j] = ((i+j) & 0x1 ? -ONE : ONE) * DataType_((i + j + 1)
            * binomial(n + i, n - j - 1) * binomial(n + j, n - i - 1)
            * sqr(binomial(i + j, i)));
        }
      }
    }

    // initialise Lehmer-Matrix; see http://en.wikipedia.org/wiki/Lehmer_matrix
    static void _init_lehmer_mat(size_t n, DataType_ a[])
    {
      for(size_t i(0) ; i < n ; ++i)
      {
        for(size_t j(0) ; j < n ; ++j)
        {
          a[i * n + j] = DataType_(min(i, j) + 1) / DataType_(max(i, j) + 1);
        }
      }
    }

    // initialise Lehmer-Matrix inverse
    static void _init_lehmer_inv(size_t n, DataType_ a[])
    {
      DataType_ b;
      mat_clear(n, n, n, a, ZERO);
      a[0] = DataType_(4) / DataType_(3);
      a[1] = b = -DataType_(2) / DataType_(3);
      for(size_t i(1) ; i < n - 1 ; ++i)
      {
        a[i * (n + 1) - 1] = b;
        a[i * (n + 1)    ] = DataType_(4*cub(i+1)) / DataType_(4*sqr(i+1) - 1);
        a[i * (n + 1) + 1] = b = -DataType_((i+2)*(i+1)) / DataType_(2*i + 3);
      }
      a[(n - 1) * (n + 1) - 1] = b;
      a[(n - 1) * (n + 1)    ] = DataType_(sqr(n)) / DataType_(2*n - 1);
    }

    // helper function: generate a random number in range [0,1]
    inline DataType_ _rnd() const
    {
      DataType_ t;
      _random >> t;
      return t;
    }

  public:
    // constructor
    LinAlgTest(const std::string& id) :
      TaggedTest<Archs::None, DataType_>(id),
      _eps(std::numeric_limits<DataType_>::epsilon())
    {
    }

    bool test_vec_clear_copy() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.9)));

#define N 16
      DataType_ x[N], y[N];

      // x := 1
      vec_clear(N, x, ONE);
      // y := x
      vec_copy(N, y, x);

      for(size_t i = 0 ; i < N ; ++i)
      {
        if(abs(y[i] - ONE) >= tol)
        {
          return false;
        }
      }

      return true;
#undef N
    }

    bool test_vec_axpy_lcomb() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.9)));

#define N 16
      DataType_ x[N], y[N], z[N];

      // initialise x[i] = cos(i)^2; y[i] = sin(i)^2
      for(size_t i = 0 ; i < N ; ++i)
      {
        x[i] = sqr(cos(DataType_(i)));
        y[i] = sqr(sin(DataType_(i)));
      }

      // calculate z := x + y
      vec_lin_comb(N, z, x, y, ONE, ONE);
      // calculate y += x
      vec_axpy(N, y, x, ONE);

      // check y[i] nad z[i] against one
      for(size_t i = 0 ; i < N ; ++i)
      {
        // test y
        if(abs(y[i] - ONE) >= tol)
        {
          return false;
        }
        // test z
        if(abs(z[i] - ONE) >= tol)
        {
          return false;
        }
      }

      return true;
#undef N
    }

    void test_vec_dot() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.9)));

#define N 16
      DataType_ x[N], y[N];

      // x[i] = i+1; y[i] = 1/(i+1)
      for(size_t i = 0 ; i < N ; ++i)
      {
        x[i] = DataType_(i+1);
        y[i] = ONE / x[i];
      }

      TEST_CHECK_EQUAL_WITHIN_EPS(vec_dot<DataType_>(N, x, y), DataType_(N), tol);
#undef N
    }

    void test_vec_norm_asum() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.9)));

#define N 16
      DataType_ x[N];

      // x[i] := 1/2^i
      for(size_t i = 0 ; i < N ; ++i)
      {
        x[i] = ONE / DataType_(1 << i);
      }

      static DataType_ r = (TWO - ONE/DataType_(1 << (N-1)));
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_norm_asum(N, x), r, tol);
#undef N
    }

    void test_vec_norm_euclid() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.9)));

#define N 16
      DataType_ x[N];

      // x[i] := 1/sqrt(2^i)
      for(size_t i = 0 ; i < N; ++i)
      {
        x[i] = ONE / sqrt(DataType_(1 << i));
      }

      static DataType_ r = sqrt(TWO - ONE/DataType_(1 << (N-1)));
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_norm_euclid(N, x), r, tol);
#undef N
    }

    void test_vec_norm_max() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.9)));

#define N 16
      static DataType_ r = (DataType_(N - 1) / DataType_(N));
      DataType_ x[N];

      // x[i] := 1 - 1/(i+1)
      for(size_t i = 0 ; i < N; ++i)
      {
        x[i] = ONE - (ONE / DataType_(i + 1));
      }

      TEST_CHECK_EQUAL_WITHIN_EPS(vec_norm_max(N, x), r, tol);
#undef N
    }

    void test_mat_axpy() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.9)));

#define M 8
#define N 16
      DataType_ a[N*M], b[N*M], c[M*N];

      // initialise a random matrix and its transpose
      for(size_t i = 0 ; i < M ; ++i)
      {
        for(size_t j = 0 ; j < N ; ++j)
        {
          a[j * M + i] = b[j * M + i] = c[i * N + j] = _rnd();
        }
      }

      // B := B - A
      mat_axpy<false>(N, M, M, b, M, a, -ONE);
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(N, M, M, b), ZERO, tol);

      // C := C - A^T
      mat_axpy<true>(M, N, N, c, M, a, -ONE);
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(M, N, N, c), ZERO, tol);
#undef M
#undef N
    }

    void test_mat_mat_mult() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.8)));

#define N 16
      // initialise a Lehmer matrix and its inverse
      DataType_ a[N*N], b[N*N], c[N*N];
      _init_lehmer_mat(N, a);
      _init_lehmer_inv(N, b);

      // C := C - A * B
      mat_identity(N, N, c);
      mat_mat_mult<false, false>(N, N, N, N, c, N, a, N, b, -ONE);
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(N, N, N, c), ZERO, tol);

      // C := C - A^T * B
      mat_identity(N, N, c);
      mat_mat_mult<true, false>(N, N, N, N, c, N, a, N, b, -ONE);
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(N, N, N, c), ZERO, tol);

      // C := C - A * B^T
      mat_identity(N, N, c);
      mat_mat_mult<false, true>(N, N, N, N, c, N, a, N, b, -ONE);
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(N, N, N, c), ZERO, tol);

      // C := C - A^T * B^T
      mat_identity(N, N, c);
      mat_mat_mult<true, true>(N, N, N, N, c, N, a, N, b, -ONE);
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(N, N, N, c), ZERO, tol);

#undef N
    }

    void test_mat_vec_solve() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.6)));

#define N 16
      DataType_ a[N*N], b[N*N], x[N], y[N], z[N];
      size_t p[N];

      // initialise a Lehmer-Matrix and its inverse
      _init_lehmer_mat(N, a);
      _init_lehmer_inv(N, b);

      // initialise a random vector
      for(size_t i = 0 ; i < N ; ++i)
      {
        x[i] = y[i] = z[i] = _rnd();
      }

      // factorise the matrix
      mat_factorise(N, N, N, a, p);

      // solve
      mat_solve_vec<false>(N, x, N, a, p);
      mat_solve_vec<true> (N, y, N, a, p);

      // subtract analytic solution
      mat_vec_mult<false>(N, N, N, x, b, z, -ONE);
      mat_vec_mult<true> (N, N, N, y, b, z, -ONE);

      TEST_CHECK_EQUAL_WITHIN_EPS(vec_norm_max(N, x), ZERO, tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_norm_max(N, x), ZERO, tol);
#undef N
    }

    void test_mat_mat_solve() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.6)));

#define M 8
#define N 16
      DataType_ a[N*N], b[N*N], x[N*M], y[N*M], z[N*M];
      size_t p[N];

      // initialise a Lehmer-Matrix and its inverse
      _init_lehmer_mat(N, a);
      _init_lehmer_inv(N, b);

      // initialise a random matrix
      for(size_t i = 0 ; i < N*M ; ++i)
      {
        x[i] = y[i] = z[i] = _rnd();
      }

      // factorise the matrix
      mat_factorise(N, N, N, a, p);

      // solve
      mat_solve_mat<false>(N, M, M, x, N, a, p);
      mat_solve_mat<true> (N, M, M, y, N, a, p);

      // subtract analytic solution
      mat_mat_mult<false,false>(N, M, N, M, x, N, b, M, z, -ONE);
      mat_mat_mult<true ,false>(N, M, N, M, y, N, b, M, z, -ONE);

      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(N, M, M, x), ZERO, tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(N, M, M, y), ZERO, tol);
#undef N
#undef M
    }

    /*
    template<size_t n_>
    void test_mat_inv_hilbert() const
    {
      // check whether to skip this test
      if(MatInvHilbert<DataType_, n_>::skip)
        return;

      // initialise a Hilbert-Matrix: a_ij := 1 / (i + j + 1); for i,j = 0,...
      DataType_ a[n_ * n_], b[n_ * n_];
      _init_hilbert_mat(n_, a);
      _init_hilbert_inv(n_, b);

      // invert matrix
      DataType_ c[n_ * n_];
      mat_invert<n_>(c, a);

      // subtract analytic inverse
      mat_axpy<false>(n_, n_, n_, c, n_, b, -ONE);

      // check norm of error
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(n_, n_, n_, c), ZERO, (MatInvHilbert<DataType_, n_>::tol()));
    }*/

    template<size_t n_>
    void test_mat_inv_lehmer() const
    {
      // set tolerance
      const DataType_ tol(std::pow(_eps, DataType_(0.8)));

      // initialise a Lehmer-Matrix
      DataType_ a[n_ * n_], b[n_ * n_];
      _init_lehmer_mat(n_, a);
      _init_lehmer_inv(n_, b);

      // invert matrix
      DataType_ c[n_ * n_];
      mat_invert<n_>(c, a);

      // subtract analytic inverse
      mat_axpy<false>(n_, n_, n_, c, n_, b, -ONE);

      // check norm of error
      TEST_CHECK_EQUAL_WITHIN_EPS(mat_norm_max(n_, n_, n_, c), ZERO, tol);
    }

    virtual void run() const
    {
      // test vector operations
      TEST_CHECK(test_vec_clear_copy());
      TEST_CHECK(test_vec_axpy_lcomb());
      test_vec_dot();
      test_vec_norm_asum();
      test_vec_norm_euclid();
      test_vec_norm_max();

      // test matrix operations
      test_mat_axpy();
      test_mat_mat_mult();
      test_mat_vec_solve();
      test_mat_mat_solve();

      // test matrix inversion
      test_mat_inv_lehmer<2>();
      test_mat_inv_lehmer<3>();
      test_mat_inv_lehmer<4>();
      test_mat_inv_lehmer<5>();
      test_mat_inv_lehmer<6>();
    }

#undef TWO
#undef ONE
#undef ZERO
}; // class LinAlgTest

LinAlgTest<double> tagged_linalg_test_d("LinearAlgebra-test double");
LinAlgTest<float>  tagged_linalg_test_f("LinearAlgebra-test float");
