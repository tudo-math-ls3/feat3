#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/transposition.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Algo_,
  typename SM_>
class TranspositionTest
  : public TaggedTest<typename SM_::MemType, typename SM_::DataType, Algo_>
{

public:
  typedef typename SM_::MemType MemType;
  typedef typename SM_::DataType DataType;

  TranspositionTest()
    : TaggedTest<MemType, DataType, Algo_>("transposition_test" + SM_::name())
  {
  }

  virtual void run() const
  {
    for (Index size(2) ; size < 3e2 ; size*=3)
    {
      SparseMatrixCOO<Mem::Main, DataType> a_local(size, size + 2);

      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DataType(2));
          else if(row == col+1)
            a_local(row, col, DataType(-1));
          else if(row+1 == col)
            a_local(row, col, DataType(-3));
        }
      }
      SM_ a;
      a.convert(a_local);

      SM_ b;
      b = Transposition<Algo_>::value(a);

      for (Index i(0) ; i < a.rows() ; ++i)
      {
        for (Index j(0) ; j < a.columns() ; ++j)
        {
          TEST_CHECK_EQUAL(b(j, i), a(i, j));
        }
      }

      a.transpose(b);

      for (Index i(0) ; i < a.rows() ; ++i)
      {
        for (Index j(0) ; j < a.columns() ; ++j)
        {
          TEST_CHECK_EQUAL(b(j, i), a(i, j));
        }
      }

      a.transpose();

      TEST_CHECK_EQUAL(a, b);
    }
  }
};
TranspositionTest<Algo::Generic, SparseMatrixCOO<Mem::Main, float, Index> > coo_cpu_transposition_test_float;
TranspositionTest<Algo::Generic, SparseMatrixCOO<Mem::Main, double, Index> > coo_cpu_transposition_test_double;
TranspositionTest<Algo::Generic, SparseMatrixCSR<Mem::Main, float, Index> > csr_cpu_transposition_test_float;
TranspositionTest<Algo::Generic, SparseMatrixCSR<Mem::Main, double, Index> > csr_cpu_transposition_test_double;
TranspositionTest<Algo::Generic, SparseMatrixELL<Mem::Main, float, Index> > ell_cpu_transposition_test_float;
TranspositionTest<Algo::Generic, SparseMatrixELL<Mem::Main, double, Index> > ell_cpu_transposition_test_double;
