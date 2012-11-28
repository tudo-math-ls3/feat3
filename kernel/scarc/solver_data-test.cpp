#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/solver_data.hpp>
#include <kernel/util/cpp11_smart_pointer.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;


template<typename Tag_, typename Algo_, typename DataType_>
class SolverDataTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    SolverDataTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("SolverDataTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      DenseVector<Tag_, DataType_> x(1000, DataType_(1)), b(1000, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T(1000, 1000);
      for(Index i(0) ; i < 1000 ; ++i)
        T(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A(T);

      SolverData<> sd0(A, x, b);
      TEST_CHECK_EQUAL(sd0.stored_sys(), A);
      TEST_CHECK_EQUAL(sd0.stored_sol(), x);
      TEST_CHECK_EQUAL(sd0.stored_rhs(), b);
      TEST_CHECK_EQUAL(sd0.stored_temp().size(), 0);

      SolverData<> sd1(A, x, b, 2);
      TEST_CHECK_EQUAL(sd1.stored_sys(), A);
      TEST_CHECK_EQUAL(sd1.stored_sol(), x);
      TEST_CHECK_EQUAL(sd1.stored_rhs(), b);
      TEST_CHECK_EQUAL(sd1.stored_temp().size(), 2);
      TEST_CHECK_EQUAL(sd1.stored_temp().at(0).size(), x.size());
      TEST_CHECK_EQUAL(sd1.stored_temp().at(1).size(), x.size());

      PreconditionedSolverData<> psd(A, A, x, b, 3);
      TEST_CHECK_EQUAL(psd.stored_sys(), A);
      TEST_CHECK_EQUAL(psd.stored_prec, A);
      TEST_CHECK_EQUAL(psd.stored_sol(), x);
      TEST_CHECK_EQUAL(psd.stored_rhs(), b);
      TEST_CHECK_EQUAL(psd.stored_temp().size(), 3);
      TEST_CHECK_EQUAL(psd.stored_temp().at(0).size(), x.size());
      TEST_CHECK_EQUAL(psd.stored_temp().at(1).size(), x.size());
      TEST_CHECK_EQUAL(psd.stored_temp().at(2).size(), x.size());

      SolverData<> sd2(psd);
      TEST_CHECK_EQUAL(sd2.stored_sys(), A);
      TEST_CHECK_EQUAL(sd2.stored_sol(), x);
      TEST_CHECK_EQUAL(sd2.stored_rhs(), b);
      TEST_CHECK_EQUAL(sd2.stored_temp().size(), 3);
      TEST_CHECK_EQUAL(sd2.stored_temp().at(0).size(), x.size());
      TEST_CHECK_EQUAL(sd2.stored_temp().at(1).size(), x.size());
      TEST_CHECK_EQUAL(sd2.stored_temp().at(2).size(), x.size());

      //-------------------------------------------------------------------------
      MultiLevelSolverData<> mlsd(A, b, x);
    }
};
SolverDataTest<Mem::Main, Algo::Generic,  double> sf_cpu_double("StorageType: std::vector, DataType: double");
