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
  typename Tag_,
  typename BType_,
  typename DT_>
class TranspositionTest
  : public TaggedTest<Tag_, DT_>
{

public:

  TranspositionTest()
    : TaggedTest<Tag_, DT_>("transposition_test")
  {
  }

  virtual void run() const
  {
    for (Index size(2) ; size < 3e2 ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
      for (unsigned long row(0) ; row < a_local.rows() ; ++row)
      {
        for (unsigned long col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
        }
      }
      SparseMatrixCSR<Tag_, DT_> a(a_local);

      SparseMatrixCSR<Tag_, DT_> b;
      b = Transposition<Tag_, BType_>::value(a);

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
TranspositionTest<Mem::Main, Algo::Generic, float> cpu_transposition_test_float;
TranspositionTest<Mem::Main, Algo::Generic, double> cpu_transposition_test_double;
