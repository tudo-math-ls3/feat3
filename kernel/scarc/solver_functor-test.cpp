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


template<typename Tag_, typename Algo_, typename DataType_>
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
      SparseMatrixCOO<Mem::Main, DataType_> T(1000, 1000);

      for(Index i(0) ; i < 1000 ; ++i)
        T(i, i, DataType_(1));

      SparseMatrixCSR<Tag_, DataType_> A(T);
      DenseVector<Tag_, DataType_> b(1000, DataType_(2));
      DenseVector<Tag_, DataType_> x(1000, DataType_(1));
      DenseVector<Tag_, DataType_> d(1000);

      //reference solution
      DenseVector<Tag_, DataType_> c(1000);
      Defect<Algo_>::value(c, b, A, x);

      DefectFunctor<Algo_, DenseVector<Tag_, DataType_>, SparseMatrixCSR<Tag_, DataType_> > f(d, b, A, x);
      f.execute();

      TEST_CHECK_EQUAL(d, c);

      //----------------------------------------------------------------------------------------------------------

      DenseVector<Tag_, DataType_> dn(1000);
      CompoundFunctor<> cf;
      cf.add_functor(new DefectFunctor<Algo_, DenseVector<Tag_, DataType_>, SparseMatrixCSR<Tag_, DataType_> >(dn, b, A, x));
      cf.add_functor(new DefectFunctor<Algo_, DenseVector<Tag_, DataType_>, SparseMatrixCSR<Tag_, DataType_> >(dn, b, A, x));
      cf.add_functor(new DefectFunctor<Algo_, DenseVector<Tag_, DataType_>, SparseMatrixCSR<Tag_, DataType_> >(dn, b, A, x));
      cf.execute();
      TEST_CHECK_EQUAL(d, c);
    }
};
SolverFunctorTest<Mem::Main, Algo::Generic, double> sf_cpu_double("StorageType: std::vector, DataType: double");
