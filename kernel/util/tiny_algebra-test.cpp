#include <kernel/util/tiny_algebra.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::Tiny;
using namespace FEAT::TestSystem;

using Math::sqr;
using Math::cub;

/**
 * \brief Tagged-Test for tiny algebra vector and matrix.
 *
 * \author Peter Zajac
 */
template<typename DataType_>
class TinyAlgebraTest :
  public TaggedTest<Archs::None, DataType_>
{
private:
  // machine exactness
  const DataType_ _eps;

  // initialise Lehmer-Matrix; see http://en.wikipedia.org/wiki/Lehmer_matrix
  template<int n_, int sm_, int sn_>
  static void _init_lehmer_mat(Matrix<DataType_, n_, n_, sm_, sn_>& a)
  {
    for(int i(0) ; i < n_ ; ++i)
    {
      for(int j(0) ; j < n_ ; ++j)
      {
        a(i,j) = DataType_(Math::min(i, j) + 1) / DataType_(Math::max(i, j) + 1);
      }
    }
  }

  // initialise Lehmer-Matrix inverse
  template<int n_, int sm_, int sn_>
  static void _init_lehmer_inv(Matrix<DataType_, n_, n_, sm_, sn_>& a)
  {
    if(n_ == 1)
    {
      a(0,0) = DataType_(1);
      return;
    }

    DataType_ b;
    a = DataType_(0);
    a(0,0) = DataType_(4) / DataType_(3);
    a(0,1) = b = -DataType_(2) / DataType_(3);
    for(int i(1) ; i < n_ - 1 ; ++i)
    {
      a(i,i-1) = b;
      a(i,i  ) = DataType_(4*cub(i+1)) / DataType_(4*sqr(i+1) - 1);
      a(i,i+1) = b = -DataType_((i+2)*(i+1)) / DataType_(2*i + 3);
    }
    a(n_-1,n_-2) = b;
    a(n_-1,n_-1) = DataType_(sqr(n_)) / DataType_(2*n_ - 1);
  }

public:
  // constructor
  explicit TinyAlgebraTest(const std::string& test_name) :
    TaggedTest<Archs::None, DataType_>(test_name),
    _eps(Math::eps<DataType_>())
  {
  }

  virtual ~TinyAlgebraTest()
  {
  }

  template<int n_>
  void test_mat_inv_lehmer() const
  {
    // set tolerance
    const DataType_ tol(Math::pow(_eps, DataType_(0.75)));

    // initialise a Lehmer-Matrix
    Matrix<DataType_, n_, n_> a, b, c;
    _init_lehmer_mat(a);
    _init_lehmer_inv(b);

    // invert matrix and subtract analytic inverse
    c.set_inverse(a);

    // compute norm
    DataType_ def = DataType_(0);
    for(int i(0); i < n_; ++i)
    {
      for(int j(0); j < n_; ++j)
      {
        def = Math::max(def, Math::abs(c[i][j] - b[i][j]));
      }
    }

    // check norm of error
    TEST_CHECK_EQUAL_WITHIN_EPS(def, DataType_(0), tol);
  }

  template<int n_>
  void test_mat_det_lehmer() const
  {
    // set tolerance
    const DataType_ tol(Math::pow(_eps, DataType_(0.75)));

    // initialise a Lehmer-Matrix
    Matrix<DataType_, n_, n_> a;
    _init_lehmer_mat(a);

    // reference matrix determinants
    static const DataType_ dets[] =
    {
      DataType_(0),
      DataType_(1),
      DataType_(3) / DataType_(4),
      DataType_(5) / DataType_(12),
      DataType_(35) / DataType_(192),
      DataType_(21) / DataType_(320),
      DataType_(77) / DataType_(3840),
      DataType_(143) / DataType_(26880),
      DataType_(143) / DataType_(114688),
      DataType_(2431) / DataType_(9289728)
    };

    // check determinat
    TEST_CHECK_EQUAL_WITHIN_EPS(a.det(), dets[n_], tol);
  }


  template<int n_>
  void test_mat_cof_lehmer() const
  {
    // Set tolerance
    const DataType_ tol(Math::pow(_eps, DataType_(0.75)));

    // Initialise a Lehmer-Matrix
    Matrix<DataType_, n_, n_> a, b, c;
    _init_lehmer_mat(a);
    _init_lehmer_inv(b);

    // Set C to the cofactor matrix and subtract the analytic cofactor values. Since the Lehmer matrix is invertible,
    // Cof(A) = 1/det(A) A^(-1)
    c.set_cofactor(a);

    DataType_ my_det = a.det();

    // compute norm
    DataType_ def = DataType_(0);
    for(int i(0); i < n_; ++i)
    {
      for(int j(0); j < n_; ++j)
      {
        def = Math::max(def, Math::abs(c[i][j] - my_det*b[i][j]));
      }
    }

    // check norm of error
    TEST_CHECK_EQUAL_WITHIN_EPS(def, DataType_(0), tol);
  }

  virtual void run() const override
  {
    // test matrix inversion
    test_mat_inv_lehmer<1>(); // specialised
    test_mat_inv_lehmer<2>(); // specialised
    test_mat_inv_lehmer<3>(); // specialised
    test_mat_inv_lehmer<4>(); // specialised
    test_mat_inv_lehmer<5>(); // specialised
    test_mat_inv_lehmer<6>(); // specialised
    test_mat_inv_lehmer<7>(); // generic
    test_mat_inv_lehmer<8>(); // generic
    test_mat_inv_lehmer<9>(); // generic

    // test matrix determinant calculation
    test_mat_det_lehmer<1>(); // specialised
    test_mat_det_lehmer<2>(); // specialised
    test_mat_det_lehmer<3>(); // specialised
    test_mat_det_lehmer<4>(); // specialised
    test_mat_det_lehmer<5>(); // specialised
    test_mat_det_lehmer<6>(); // specialised
    test_mat_det_lehmer<7>(); // generic
    test_mat_det_lehmer<8>(); // generic
    test_mat_det_lehmer<9>(); // generic

    test_mat_cof_lehmer<2>(); // specialised
    test_mat_cof_lehmer<3>(); // specialised
    test_mat_cof_lehmer<4>(); // specialised
    test_mat_cof_lehmer<5>(); // specialised
    test_mat_cof_lehmer<6>(); // specialised
    test_mat_cof_lehmer<7>(); // generic
    test_mat_cof_lehmer<8>(); // generic
    test_mat_cof_lehmer<9>(); // generic
  }
};

TinyAlgebraTest<float> tagged_tiny_test_f("TinyAlgebraTest<float>");
TinyAlgebraTest<double> tagged_tiny_test_d("TinyAlgebraTest<double>");
#ifdef FEAT_HAVE_QUADMATH
TinyAlgebraTest<__float128> tagged_tiny_test_f128("TinyAlgebraTest<__float128>");
#endif // FEAT_HAVE_QUADMATH
