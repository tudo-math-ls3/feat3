#define SERIAL
#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/util/cpp11_smart_pointer.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/product_matvec.hpp>
#include <kernel/scarc/solver_functor.hpp>
#include <kernel/scarc/solver_data.hpp>
#include <kernel/scarc/solver_pattern.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;

template<typename Tag_, typename Algo_, typename DataType_>
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
      ProductMatVec<Algo_>::value(b, A, x_ref);

      SparseMatrixCOO<Tag_, DataType_> P_proxy(x.size(), x.size());
      for(Index i(0) ; i < x.size() ; ++i)
      {
        P_proxy(i, i, DataType_(0.75 * (1. / A(i, i))));
      }
      SparseMatrixELL<Tag_, DataType_> P(P_proxy);

      PreconditionedSolverData<DataType_, Tag_, DenseVector, SparseMatrixELL, SparseMatrixELL<Tag_, DataType_> >data(A, P, x, b,
                        SolverPatternGeneration<Richardson, Algo_>::min_num_temp_vectors(),
                        SolverPatternGeneration<Richardson, Algo_>::min_num_temp_scalars());
      std::shared_ptr<SolverFunctorBase<DenseVector<Tag_, DataType_> > > solver(SolverPatternGeneration<Richardson, Algo_>::execute(data, 2000, 1e-8));
      solver->execute();

      for(Index i(0) ; i < x.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(i), x_ref(i), 1e-6);


      //synchronised version must at least (uses other reduction mechanism) deliver the same
      DenseVector<Tag_, DataType_> x1(A.rows(), DataType_(0));
      SynchronisedPreconditionedSolverData<DataType_,
                                           Tag_,
                                           DenseVector,
                                           VectorMirror,
                                           SparseMatrixELL,
                                           SparseMatrixELL<Tag_, DataType_> >data1(A, P, x1, b,
                        SolverPatternGeneration<Richardson, Algo_>::min_num_temp_vectors(),
                        SolverPatternGeneration<Richardson, Algo_>::min_num_temp_scalars());
      std::shared_ptr<SolverFunctorBase<DenseVector<Tag_, DataType_> > > solver1(SolverPatternGeneration<Richardson, Algo_>::execute(data1, 2000, 1e-8));
      solver1->execute();

      for(Index i(0) ; i < x1.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(data1.sol()(i), x_ref(i), 1e-6);
    }
};
SolverTest<Mem::Main, Algo::Generic,  double> sf_cpu_double("ELL double");
