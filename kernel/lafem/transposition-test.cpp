#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/transposition.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename SM_>
class TranspositionTest
  : public TaggedTest<Mem_, DT_, Algo_>
{

public:

  TranspositionTest()
    : TaggedTest<Mem_, DT_, Algo_>("transposition_test" + SM_::name())
  {
  }

  virtual void run() const
  {
    for (Index size(2) ; size < 3e2 ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
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
    }
  }
};
TranspositionTest<Mem::Main, Algo::Generic, float, SparseMatrixCOO<Mem::Main, float> > coo_cpu_transposition_test_float;
TranspositionTest<Mem::Main, Algo::Generic, double, SparseMatrixCOO<Mem::Main, double> > coo_cpu_transposition_test_double;
TranspositionTest<Mem::Main, Algo::Generic, float, SparseMatrixCSR<Mem::Main, float> > csr_cpu_transposition_test_float;
TranspositionTest<Mem::Main, Algo::Generic, double, SparseMatrixCSR<Mem::Main, double> > csr_cpu_transposition_test_double;
TranspositionTest<Mem::Main, Algo::Generic, float, SparseMatrixELL<Mem::Main, float> > ell_cpu_transposition_test_float;
TranspositionTest<Mem::Main, Algo::Generic, double, SparseMatrixELL<Mem::Main, double> > ell_cpu_transposition_test_double;
