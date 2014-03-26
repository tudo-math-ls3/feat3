#include <kernel/base_header.hpp>
//#ifdef SERIAL
#include <test_system/test_system.hpp>

#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/scarc/scarc_functor.hpp>
#include <kernel/scarc/scarc_data.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;

template<typename Tag_, typename Algo_, typename DataType_>
class ScaRCFunctorTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    ScaRCFunctorTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("ScaRCFunctorTest<" + tag + ">")
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
      A.template apply<Algo_>(b, x_ref);

      SparseMatrixCOO<Tag_, DataType_> P_proxy(x.size(), x.size());
      for(Index i(0) ; i < x.size() ; ++i)
      {
        P_proxy(i, i, DataType_(0.75 * (1. / A(i, i))));
      }
      SparseMatrixELL<Tag_, DataType_> P(P_proxy);

      //synchronised version must at least (uses other reduction mechanism) deliver the same
      SparseMatrixELL<Tag_, DataType_> A1(fcoo);
      DenseVector<Tag_, DataType_> b1(A1.rows());
      A1.template apply<Algo_>(b1, x_ref);
      DenseVector<Tag_, DataType_> x1(A1.rows(), DataType_(0));
      SparseMatrixELL<Tag_, DataType_> P1(P_proxy);

      std::cout << A1 << std::endl;

      auto pdata1(std::shared_ptr<ScaRCDataBase<DataType_,
          Tag_,
          DenseVector<Tag_, DataType_>,
          SparseMatrixELL<Tag_, DataType_>,
          std::vector> >(new SynchronisedPreconditionedScaRCData<DataType_,
            Tag_,
            DenseVector<Tag_, DataType_>,
            VectorMirror<Tag_, DataType_>,
            SparseMatrixELL<Tag_, DataType_>,
            SparseMatrixELL<Tag_, DataType_> >(std::move(A1), std::move(P1), std::move(x1), std::move(b1))));

      std::cout <<  pdata1->sol() << std::endl;

      /*auto solver1(std::shared_ptr<ScaRCFunctorBase<DataType_, Tag_, DenseVector, SparseMatrixELL> >(new ScaRCFunctorRichardson<DataType_, Tag_, DenseVector, SparseMatrixELL>(pdata1)));
      solver1->execute();

      std::cout <<  pdata1->sol() << std::endl;

      for(Index i(0) ; i < pdata1->sol().size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(pdata1->sol()(i), x_ref(i), 1e-6);
      std::cout << "M7" << std::endl;*/
    }
};
ScaRCFunctorTest<Mem::Main, Algo::Generic,  double> sf_cpu_double("ELL double");
//#endif // SERIAL
