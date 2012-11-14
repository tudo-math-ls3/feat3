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
  typename Algo_,
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
      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
      DenseVector<Mem::Main, DT_> b_local(size + 2);
      DenseVector<Mem::Main, DT_> rhs_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
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
      for (Index i(0) ; i < size ; ++i)
      {
        b_local(i, std::sin(pi * DT_(i) / DT_(size-1)));
        rhs_local(i, b_local(i) - DT_(5));
      }

      SparseMatrixCSR<Arch_, DT_> a(a_local);
      DenseVector<Arch_, DT_> b(size + 2);
      copy(b, b_local);
      DenseVector<Arch_, DT_> rhs(size);
      copy(rhs, rhs_local);
      DenseVector<Arch_, DT_> c(size, 4711);

      Defect<Algo_>::value(c, rhs, a, b);
      copy(result_local, c);

      DT_ dev(DT_(0));
      for (Index i(0) ; i < size ; ++i)
      {
        DT_ sum(result_local(i) - rhs_local(i));
        for (Index j(0) ; j < a_local.columns() ; ++j)
          sum += a_local(i, j) * b_local(j);
        dev = std::max(dev, std::abs(sum));
      }

      TEST_CHECK(dev <= eps);
    }
  }
};
CSRDefectTest<Mem::Main, Algo::Generic, float> csrv_defect_test_float;
CSRDefectTest<Mem::Main, Algo::Generic, double> csrv_defect_test_double;
#ifdef FEAST_BACKENDS_CUDA
CSRDefectTest<Mem::CUDA, Algo::CUDA, float> cuda_csrv_defect_test_float;
CSRDefectTest<Mem::CUDA, Algo::CUDA, double> cuda_csrv_defect_test_double;
#endif
