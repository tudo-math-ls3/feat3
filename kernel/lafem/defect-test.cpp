#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/defect.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class CSRDefectTest
  : public TaggedTest<Arch_, DT_>
{

public:

  CSRDefectTest()
    : TaggedTest<Arch_, DT_>("defect_test")
  {
  }

  virtual void run() const
  {
    const DT_ eps = std::pow(std::numeric_limits<DT_>::epsilon(), DT_(0.8));
    const DT_ pi = DT_(2) * std::acos(DT_(0));

    for (Index size(2) ; size < 3e2 ; size*=2)
    {
      SparseMatrixCOO<Arch_, DT_> ac(size, size + 2);
      DenseVector<Arch_, DT_> b(size + 2);
      DenseVector<Arch_, DT_> rhs(size);
      DenseVector<Arch_, DT_> c(size, 4711);
      DenseVector<Arch_, DT_> ref(size);
      for (unsigned long row(0) ; row < ac.rows() ; ++row)
      {
        for (unsigned long col(0) ; col < ac.columns() ; ++col)
        {
          if(row == col)
            ac(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            ac(row, col, DT_(-1));
        }
      }
      SparseMatrixCSR<Arch_, DT_> a(ac);
      for (Index i(0) ; i < size ; ++i)
      {
        b(i, std::sin(pi * DT_(i) / DT_(size-1)));
        rhs(i, b(i) - DT_(5));
      }

      Defect<Arch_, BType_>::value(c, rhs, a, b);

      DT_ dev(DT_(0));
      for (Index i(0) ; i < size ; ++i)
      {
        DT_ sum(c(i) - rhs(i));
        for (Index j(0) ; j < a.columns() ; ++j)
          sum += ac(i, j) * b(j);
        dev = std::max(dev, std::abs(sum));
      }

      TEST_CHECK(dev <= eps);
    }
  }
};
CSRDefectTest<Archs::CPU, Archs::Generic, float> csrv_defect_test_float;
CSRDefectTest<Archs::CPU, Archs::Generic, double> csr_defect_test_double;
