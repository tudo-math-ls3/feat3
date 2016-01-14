#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/solver_functor.hpp>
#include <kernel/scarc/solver_pattern.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;


template<typename Tag_, typename DataType_>
class SolverPatternTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    SolverPatternTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("SolverPatternTest<" + tag + ">")
    {
    }

    virtual void run() const override
    {
      DenseVector<Tag_, DataType_> x(100, DataType_(1)), b(100, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T(100, 100);
      for(Index i(0) ; i < 100 ; ++i)
        T(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A(T);
      SparseMatrixCSR<Tag_, DataType_> P;
      P.clone(A);

      PreconditionedSolverData<> data(std::move(A), std::move(P), std::move(x), std::move(b),
                        SolverPatternGeneration<Richardson>::min_num_temp_vectors(),
                        SolverPatternGeneration<Richardson>::min_num_temp_scalars());
      std::shared_ptr<SolverFunctorBase<DenseVector<Tag_, DataType_> > > solver(SolverPatternGeneration<Richardson>::execute(data, 20, 1e-8));
      TEST_CHECK_EQUAL(solver->type_name(), "CompoundSolverFunctor[DefectFunctor, NormFunctor2, IterateFunctor[CompoundSolverFunctor[ProductFunctor, SumFunctor, DefectFunctor, NormFunctor2, DivFunctor]]]");
      solver->execute();

      DenseVector<Tag_, DataType_> x1(100, DataType_(1)), b1(100, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T1(100, 100);
      for(Index i(0) ; i < 100 ; ++i)
        T1(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A1(T1);
      SparseMatrixCSR<Tag_, DataType_> P1;
      P1.clone(A1);

      //---------------------------------------------------------------------------------------------------------------------------------------------
      DenseVector<Tag_, DataType_> dummy;
      PreconditionedSolverData<> data3(std::move(A1), std::move(P1), std::move(x1), std::move(b1),
                         SolverPatternGeneration<RichardsonLayer>::min_num_temp_vectors(),
                         SolverPatternGeneration<RichardsonLayer>::min_num_temp_scalars());
      std::shared_ptr<SolverFunctorBase<DenseVector<Tag_, DataType_> > > solver3(SolverPatternGeneration<RichardsonLayer>::execute(data3, dummy, 20, 1e-8));
      TEST_CHECK_EQUAL(solver3->type_name(), "CompoundSolverFunctor[DefectFunctor, NormFunctor2, IterateFunctor[CompoundSolverFunctor[ProductFunctor, SumFunctor, DefectFunctor, NormFunctor2, DivFunctor]]]");
      TEST_CHECK_THROWS(solver3->execute(), ScaRCError);

    }
};
SolverPatternTest<Mem::Main, double> sf_cpu_double("StorageType: std::vector, DataType: double");
