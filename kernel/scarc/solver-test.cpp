#include <kernel/base_header.hpp>
#ifdef SERIAL
#include <test_system/test_system.hpp>

#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/scarc/solver_functor.hpp>
#include <kernel/scarc/solver_data.hpp>
#include <kernel/scarc/solver_pattern.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;

template<typename Tag_, typename DataType_>
class SolverTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    SolverTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("SolverTest<" + tag + ">")
    {
    }

    virtual void run() const
    {

      SparseMatrixCOO<Mem::Main, DataType_> fcoo(10, 10);
      for (unsigned long row(0) ; row < fcoo.rows() ; ++row)
      {
        for (unsigned long col(0) ; col < fcoo.columns() ; ++col)
        {
          if(row == col)
            fcoo(row, col, DataType_(2));
          else if((row == col+1) || (row+1 == col))
            fcoo(row, col, DataType_(-1));
        }
      }
      SparseMatrixELL<Tag_, DataType_> A(fcoo);

      DenseVector<Tag_, DataType_> x_ref(A.rows(), DataType_(2));
      DenseVector<Tag_, DataType_> x(A.rows(), DataType_(0));
      DenseVector<Tag_, DataType_> b(A.rows());
      A.apply(b, x_ref);

      SparseMatrixCOO<Tag_, DataType_> P_proxy(x.size(), x.size());
      for(Index i(0) ; i < x.size() ; ++i)
      {
        P_proxy(i, i, DataType_(0.75 * (1. / A(i, i))));
      }
      SparseMatrixELL<Tag_, DataType_> P(P_proxy);

      PreconditionedSolverData<DataType_, Tag_, DenseVector, SparseMatrixELL, SparseMatrixELL >data(std::move(A), std::move(P), std::move(x), std::move(b),
                        SolverPatternGeneration<Richardson>::min_num_temp_vectors(),
                        SolverPatternGeneration<Richardson>::min_num_temp_scalars());
      std::shared_ptr<SolverFunctorBase<DenseVector<Tag_, DataType_> > > solver(SolverPatternGeneration<Richardson>::execute(data, 2000, 1e-8));
      solver->execute();

      for(Index i(0) ; i < data.sol().size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(i), x_ref(i), 1e-6);

      //synchronised version must at least (uses other reduction mechanism) deliver the same
      SparseMatrixELL<Tag_, DataType_> A1(fcoo);
      DenseVector<Tag_, DataType_> b1(A1.rows());
      A1.apply(b1, x_ref);
      DenseVector<Tag_, DataType_> x1(A1.rows(), DataType_(0));
      SparseMatrixELL<Tag_, DataType_> P1(P_proxy);

      SynchronisedPreconditionedSolverData<DataType_,
                                           Tag_,
                                           DenseVector,
                                           VectorMirror,
                                           SparseMatrixELL,
                                           SparseMatrixELL >data1(std::move(A1), std::move(P1), std::move(x1), std::move(b1),
                        SolverPatternGeneration<Richardson>::min_num_temp_vectors(),
                        SolverPatternGeneration<Richardson>::min_num_temp_scalars());
      std::shared_ptr<SolverFunctorBase<DenseVector<Tag_, DataType_> > > solver1(SolverPatternGeneration<Richardson>::execute(data1, 2000, 1e-8));
      solver1->execute();

      for(Index i(0) ; i < data1.sol().size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(data1.sol()(i), x_ref(i), 1e-6);
      std::cout << "M7" << std::endl;
    }
};
SolverTest<Mem::Main, double> sf_cpu_double("ELL double");
#endif // SERIAL
