#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/container_pool.hpp>


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the container pool class.
*
* \test test description missing
*
* \tparam Mem_
* description missing
*
* \tparam DT_
* description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Mem_,
  typename DT_,
  typename CT_>
class ContainerPoolTest
  : public TaggedTest<Mem_, DT_>
{

public:

  ContainerPoolTest()
    : TaggedTest<Mem_, DT_>("container_pool_test " + CT_::type_name())
  {
  }

  virtual void run() const
  {
    std::vector<CT_> list;
    CT_ z;
    {
      CT_ a(10);
      CT_ y(a);
      std::cout<<a.size()<<" "<<y.size()<<std::endl;
      //ContainerPool<CT_>::instance()->push_back(a);
      list.push_back(a);
      z = a.clone();
    }
    CT_ b(8);
    //ContainerPool<CT_>::instance()->push_back(b);
    list.push_back(b);

    //CT_ c = ContainerPool<CT_>::instance()->at(0);
    //CT_ d(ContainerPool<CT_>::instance()->at(1));
    CT_ c = list.at(0);
    CT_ d = list.at(1);

    TEST_CHECK_EQUAL(c.size(), z.size());
    TEST_CHECK_EQUAL(d.size(), b.size());
  }
};
ContainerPoolTest<Mem::Main, float, DenseVector<Mem::Main, float> > cpu_dv_container_pool_test_float;
ContainerPoolTest<Mem::Main, double, DenseVector<Mem::Main, double> > cpu_dv_container_pool_test_double;
#ifdef FEAST_BACKENDS_CUDA
ContainerPoolTest<Mem::CUDA, float, DenseVector<Mem::CUDA, float> > gpu_dv_container_pool_test_float;
ContainerPoolTest<Mem::CUDA, double, DenseVector<Mem::CUDA, double> > gpu_dv_container_pool_test_double;
#endif
ContainerPoolTest<Mem::Main, float, SparseMatrixCOO<Mem::Main, float> > cpu_coo_container_pool_test_float;
ContainerPoolTest<Mem::Main, double, SparseMatrixCOO<Mem::Main, double> > cpu_coo_container_pool_test_double;
#ifdef FEAST_BACKENDS_CUDA
ContainerPoolTest<Mem::CUDA, float, SparseMatrixCOO<Mem::CUDA, float> > gpu_coo_container_pool_test_float;
ContainerPoolTest<Mem::CUDA, double, SparseMatrixCOO<Mem::CUDA, double> > gpu_coo_container_pool_test_double;
#endif
