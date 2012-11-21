#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/solver_functor.hpp>
#include <kernel/util/cpp11_smart_pointer.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;
using namespace FEAST::TestSystem;


template<typename Tag_, typename DataType_>
class SolverFunctorTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    SolverFunctorTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("SolverFunctorTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      SparseMatrixCOO<Tag_, DataType_> T(1000, 1000);

      for(Index i(0) ; i < 1000 ; ++i)
        T(i, i, DataType_(1));

      SparseMatrixCSR<Tag_, DataType_> A(T);
      DenseVector<Tag_, DataType_> b(1000, DataType_(2));
      DenseVector<Tag_, DataType_> x(1000, DataType_(1));
      DenseVector<Tag_, DataType_> c(1000);
      DenseVector<Tag_, DataType_> d(1000);

      DefectFunctor<Algo::Generic, DenseVector<Tag_, DataType_>, SparseMatrixCSR<Tag_, DataType_> > f(d, b, A, x);
      f.execute();

      Defect<Algo::Generic>::value(c, b, A, x);

      TEST_CHECK_EQUAL(d, c);
    }
};
SolverFunctorTest<Mem::Main, double> sf_cpu_double("StorageType: std::vector, DataType: double");
