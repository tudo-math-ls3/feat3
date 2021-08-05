#include <area51/ccnd_fiber/tensor_operations.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::Tiny;
using namespace FEAT::TestSystem;
using namespace FEAT::CCND_FIBER;

using Math::sqr;
using Math::cub;

/**
 * \brief Tagged-Test for tensor operations for fourth moment tensors.
 *
 * \author Maximilian Esser
 */
template<typename DataType_>
class TensorOperatorTest :
  public TaggedTest<Archs::None, DataType_>
{
private:
  // machine exactness
  const DataType_ _eps;

public:
  // constructor
  explicit TensorOperatorTest(const std::string& test_name) :
  TaggedTest<Archs::None, DataType_>(test_name),
  _eps(Math::eps<DataType_>())
  {
  }

  virtual ~TensorOperatorTest()
  {
  }

  void test_symmetric_tensor_contraction_2D() const
  {
    // Set tolerance
    const DataType_ tol(Math::pow(_eps, DataType_(0.75))*DataType_(10));

    //initialize the symmetric tensor as vector
    Vector<DataType_, 16> t(DataType_(0));
    //write values in
    t(0) = DataType_(0.3); //A1111
    t(15) = DataType_(-0.334); //A2222

    t(0*8 + 0*4 + 0*2 + 1) = t(0*8 + 0*4 + 1*2 + 0) = t(0*8 + 1*4 + 0*2 + 0) = t(1*8 + 0*4 + 0*2 + 0) = DataType_(5.7); //A1112
    t(1*8 + 1*4 + 1*2 + 0) = t(1*8 + 1*4 + 0*2 + 1) = t(1*8 + 0*4 + 1*2 + 1) = t(0*8 + 1*4 + 1*2 + 1) = DataType_(2.082); //A1222

    t(0*8 + 0*4 + 1*2 + 1) = t(0*8 + 1*4 + 0*2 + 1) = t(1*8 + 0*4 + 0*2 + 1) = t(1*8 + 1*4 + 0*2 + 0) = t(1*8 + 0*4 + 1*2 + 0) = t(0*8 + 1*4 + 1*2 + 0) = DataType_(-1.77); //A1122, A1212 A2112 A2121 A2211 A1221

    //test if some value was not initialized
    for(int i(0); i < 16; ++i)
    {
      TEST_CHECK_MSG(Math::abs(t(i)) > tol, "Value " + stringify(i) + " was not intialized.");
    }

    //now initialize the symmetric reduced vector...
    Vector<DataType_, 5> t_sym(DataType_(0));
    t_sym(0) = t(0); //A1111
    t_sym(1) = t(15); //A2222
    t_sym(2) = t(0*8 + 0*4 + 0*2 + 1); //A1112
    t_sym(3) = t(1*8 + 1*4 + 1*2 + 0); //A2221
    t_sym(4) = t(0*8 + 0*4 + 1*2 + 1); //A1122


    //now we test our symmetric and non symmetric implementations against each other

    //define two vectors, who represent our matrix... should not be symmetric in any way
    Vector<DataType_, 2> a,b;
    a(0) = DataType_(1.2);
    a(1) = DataType_(-3.3);
    b(0) = DataType_(9.1679);
    b(1) = DataType_(2.003);


    //define our two matrices, in which we will write the results
    Matrix<DataType_, 2, 2> mat_A, mat_B;

    //since t is symmetric, contraction_ij should be the same for all i and j
    {
      mat_A.format();
      mat_B.format();

      set_tensor4_outer_product_contraction_12(mat_A, a, b, t);
      set_tensor4_outer_product_contraction_13(mat_B, a, b, t);
      for(int i(0); i < 2; ++i)
      {
        for(int j(0); j < 2; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      set_tensor4_outer_product_contraction_14(mat_B, a, b, t);
      for(int i(0); i < 2; ++i)
      {
        for(int j(0); j < 2; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      set_tensor4_outer_product_contraction_23(mat_B, a, b, t);
      for(int i(0); i < 2; ++i)
      {
        for(int j(0); j < 2; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      set_tensor4_outer_product_contraction_24(mat_B, a, b, t);
      for(int i(0); i < 2; ++i)
      {
        for(int j(0); j < 2; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      set_tensor4_outer_product_contraction_34(mat_B, a, b, t);
      for(int i(0); i < 2; ++i)
      {
        for(int j(0); j < 2; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }
    }

    //now test if the symmetric version is the same as our non symmetric version for a symmetric tensor
    {
      mat_A.format();
      mat_B.format();

      set_tensor4_outer_product_contraction_12(mat_A, a, b, t);
      mat_B = fourth_order_contraction(a, b, t_sym, DataType_(1));
      for(int i(0); i < 2; ++i)
      {
        for(int j(0); j < 2; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }
    }

    //test whether order for symmteric variant is irrelevant
    {
      mat_A.format();
      mat_B.format();

      mat_A = fourth_order_contraction(b, a, t_sym, DataType_(1.36));
      mat_B = fourth_order_contraction(a, b, t_sym, DataType_(1.36));
      for(int i(0); i < 2; ++i)
      {
        for(int j(0); j < 2; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }
    }

    //done

  }


  void test_symmetric_tensor_contraction_3D() const
  {
    // Set tolerance
    const DataType_ tol(Math::pow(_eps, DataType_(0.75))*DataType_(10));

    //initialize the symmetric tensor as vector
    Vector<DataType_, 81> t(DataType_(0));
    //write values in
    t(0) = DataType_(0.3); //A1111
    t(27 + 9 + 3 + 1) = DataType_(-0.334); //A2222
    t(80) = DataType_(0.243); //A3333

    t(0*27 + 0*9 + 0*3 + 1) = t(0*27 + 0*9 + 1*3 + 0) = t(0*27 + 1*9 + 0*3 + 0) = t(1*27 + 0*9 + 0*3 + 0) = DataType_(5.7); //A1112
    t(0*27 + 0*9 + 0*3 + 2) = t(0*27 + 0*9 + 2*3 + 0) = t(0*27 + 2*9 + 0*3 + 0) = t(2*27 + 0*9 + 0*3 + 0) = DataType_(1.3176); //A1113
    t(1*27 + 1*9 + 1*3 + 0) = t(1*27 + 1*9 + 0*3 + 1) = t(1*27 + 0*9 + 1*3 + 1) = t(0*27 + 1*9 + 1*3 + 1) = DataType_(2.082); //A1222
    t(1*27 + 1*9 + 1*3 + 2) = t(1*27 + 1*9 + 2*3 + 1) = t(1*27 + 2*9 + 1*3 + 1) = t(2*27 + 1*9 + 1*3 + 1) = DataType_(0.854); //A2223
    t(2*27 + 2*9 + 2*3 + 0) = t(2*27 + 2*9 + 0*3 + 2) = t(2*27 + 0*9 + 2*3 + 2) = t(0*27 + 2*9 + 2*3 + 2) = DataType_(-1.302); //A1333
    t(2*27 + 2*9 + 2*3 + 1) = t(2*27 + 2*9 + 1*3 + 2) = t(2*27 + 1*9 + 2*3 + 2) = t(1*27 + 2*9 + 2*3 + 2) = DataType_(-2.93); //A2333

    t(0*27 + 0*9 + 1*3 + 1) = t(0*27 + 1*9 + 0*3 + 1) = t(1*27 + 0*9 + 0*3 + 1) = t(1*27 + 1*9 + 0*3 + 0) = t(1*27 + 0*9 + 1*3 + 0) = t(0*27 + 1*9 + 1*3 + 0) = DataType_(0.49965); //A1122, A1212 A2112 A2121 A2211 A1221
    t(0*27 + 0*9 + 2*3 + 2) = t(0*27 + 2*9 + 0*3 + 2) = t(2*27 + 0*9 + 0*3 + 2) = t(2*27 + 2*9 + 0*3 + 0) = t(2*27 + 0*9 + 2*3 + 0) = t(0*27 + 2*9 + 2*3 + 0) = DataType_(-1.77); //A1133
    t(1*27 + 1*9 + 2*3 + 2) = t(1*27 + 2*9 + 1*3 + 2) = t(2*27 + 1*9 + 1*3 + 2) = t(2*27 + 2*9 + 1*3 + 1) = t(2*27 + 1*9 + 2*3 + 1) = t(1*27 + 2*9 + 2*3 + 1) = DataType_(0.9534); //A2233

    t(0*27 + 0*9 + 1*3 + 2) = t(0*27 + 0*9 + 2*3 + 1) = t(0*27 + 1*9 + 0*3 + 2) = t(0*27 + 1*9 + 2*3 + 0) = t(1*27 + 0*9 + 2*3 + 0) = t(1*27 + 0*9 + 0*3 + 2) //A1123, A1132, A1213, A1231, A2131, A2113
      = t(1*27 + 2*9 + 0*3 + 0) = t(2*27 + 0*9 + 0*3 + 1) = t(2*27 + 1*9 + 0*3 + 0) = t(2*27 + 0*9 + 1*3 + 0) = t(0*27 + 2*9 + 1*3 + 0) = t(0*27 + 2*9 + 0*3 + 1) //A2311, A3112, A3211, A3121, A1321, A1312
      = DataType_(1.635);

    t(1*27 + 1*9 + 0*3 + 2) = t(1*27 + 1*9 + 2*3 + 0) = t(1*27 + 0*9 + 1*3 + 2) = t(1*27 + 0*9 + 2*3 + 1) = t(0*27 + 1*9 + 2*3 + 1) = t(0*27 + 1*9 + 1*3 + 2) //A2213
      = t(0*27 + 2*9 + 1*3 + 1) = t(2*27 + 1*9 + 1*3 + 0) = t(2*27 + 0*9 + 1*3 + 1) = t(2*27 + 1*9 + 0*3 + 1) = t(1*27 + 2*9 + 0*3 + 1) = t(1*27 + 2*9 + 1*3 + 0)
      = DataType_(-0.178);

    t(2*27 + 2*9 + 0*3 + 1) = t(2*27 + 2*9 + 1*3 + 0) = t(2*27 + 0*9 + 2*3 + 1) = t(2*27 + 0*9 + 1*3 + 2) = t(0*27 + 2*9 + 1*3 + 2) = t(0*27 + 2*9 + 2*3 + 1) //A3312
      = t(0*27 + 1*9 + 2*3 + 2) = t(1*27 + 2*9 + 2*3 + 0) = t(1*27 + 0*9 + 2*3 + 2) = t(1*27 + 2*9 + 0*3 + 2) = t(2*27 + 1*9 + 0*3 + 2) = t(2*27 + 1*9 + 2*3 + 0)
      = DataType_(3.418);

    //test if some value was not initialized
    for(int i(0); i < 81; ++i)
    {
      TEST_CHECK_MSG(Math::abs(t(i)) > tol, "Value " + stringify(i) + " was not intialized.");
    }

    //now initialize the symmetric reduced vector...
    Vector<DataType_, 15> t_sym(DataType_(0));
    t_sym(0) = t(0); //A1111
    t_sym(1) = t(27+9+3+1); //A2222
    t_sym(2) = t(80); //A3333
    t_sym(3) = t(0*27 + 0*9 + 0*3 + 1); //A1112
    t_sym(4) = t(0*27 + 0*9 + 0*3 + 2); //A1113
    t_sym(5) = t(1*27 + 1*9 + 1*3 + 0); //A2221
    t_sym(6) = t(1*27 + 1*9 + 1*3 + 2); //A2223
    t_sym(7) = t(2*27 + 2*9 + 2*3 + 0); //A3331
    t_sym(8) = t(2*27 + 2*9 + 2*3 + 1); //A3332
    t_sym(9) = t(0*27 + 0*9 + 1*3 + 1); //A1122
    t_sym(10) = t(0*27 + 0*9 + 2*3 + 2); //A1133
    t_sym(11) = t(1*27 + 1*9 + 2*3 + 2); //A2233
    t_sym(12) = t(0*27 + 0*9 + 1*3 + 2); //A1123
    t_sym(13) = t(1*27 + 1*9 + 0*3 + 2); //A2213
    t_sym(14) = t(2*27 + 2*9 + 0*3 + 1); //A3312

    //now we test our symmetric and non symmetric implementations against each other

    //define two vectors, who represent our matrix... should not be symmetric in any way
    Vector<DataType_, 3> a,b;
    a(0) = DataType_(1.2);
    a(1) = DataType_(-3.3);
    a(2) = DataType_(-9.23);
    b(0) = DataType_(1.1679);
    b(1) = DataType_(2.003);
    b(2) = DataType_(-4.7);

    //define our two matrices, in which we will write the results
    Matrix<DataType_, 3, 3> mat_A, mat_B;

    //first of all, the add and set method should result in the same matrix, if formated beforehand
    {
      mat_A.format();
      mat_B.format();

      add_tensor4_outer_product_contraction_12(mat_A, a, b, t);
      set_tensor4_outer_product_contraction_12(mat_B, a, b, t);
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      mat_A.format();
      mat_B.format();

      add_tensor4_outer_product_contraction_symmetric(mat_A, a, b, t_sym);
      set_tensor4_outer_product_contraction_symmetric(mat_B, a, b, t_sym);
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }
    }

    //since t is symmetric, contraction_ij should be the same for alle i and j
    {
      mat_A.format();
      mat_B.format();

      set_tensor4_outer_product_contraction_12(mat_A, a, b, t);
      set_tensor4_outer_product_contraction_13(mat_B, a, b, t);
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      set_tensor4_outer_product_contraction_14(mat_B, a, b, t);
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      set_tensor4_outer_product_contraction_23(mat_B, a, b, t);
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      set_tensor4_outer_product_contraction_24(mat_B, a, b, t);
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }

      set_tensor4_outer_product_contraction_34(mat_B, a, b, t);
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }
    }

    //now test if the symmetric version is the same as our non symmetric version for a symmetric tensor
    {
      mat_A.format();
      mat_B.format();

      set_tensor4_outer_product_contraction_12(mat_A, a, b, t);
      set_tensor4_outer_product_contraction_symmetric(mat_B, a, b, t_sym);
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }
    }

    //also test the contraction function
    {
      set_tensor4_outer_product_contraction_symmetric(mat_A, a, b, t_sym, DataType_(1.34));
      mat_B = fourth_order_contraction(a, b, t_sym, DataType_(1.34));
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }
    }

    //test whether order for symmteric variant is irrelevant
    {
      mat_A.format();
      mat_B.format();

      mat_A = fourth_order_contraction(b, a, t_sym, DataType_(1.36));
      mat_B = fourth_order_contraction(a, b, t_sym, DataType_(1.36));
      for(int i(0); i < 3; ++i)
      {
        for(int j(0); j < 3; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(mat_A[i][j], mat_B[i][j], tol);
        }
      }
    }
    //done
  }

  virtual void run() const override
  {
    //test symmetric tensor4 contraction in 2D and 3D
    test_symmetric_tensor_contraction_2D();
    test_symmetric_tensor_contraction_3D();
  }
};

TensorOperatorTest<float> tagged_tensor_operator_test_f("TensorOperatorTest<float>");
TensorOperatorTest<double> tagged_tensor_operator_test_d("TensorOperatorTest<double>");
#ifdef FEAT_HAVE_QUADMATH
TensorOperatorTest<__float128> tagged_tensor_operator_test_f128("TensorOperatorTest<__float128>");
#endif // FEAT_HAVE_QUADMATH
