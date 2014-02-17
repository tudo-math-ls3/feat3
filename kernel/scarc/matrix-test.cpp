#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/matrix.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_, typename DataType_>
class MatrixTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    MatrixTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("MatrixTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      ScaRC::DynamicAOSMatrix<> m1(2000, 3000);
      ScaRC::DynamicAOSMatrix<> m2(m1);

      m2.add_block(m1, 100, 100);

      TEST_CHECK_EQUAL(m2.rows(), 2000ul);
      TEST_CHECK_EQUAL(m2.columns(), 3000ul);
      TEST_CHECK_EQUAL(m2.num_blocks(), 1ul);
      TEST_CHECK_EQUAL(m1.num_blocks(), 0ul);
    }
};
MatrixTest<Mem::Main, double> attribute_test_cpu_v_ulong_float("StorageType: std::vector, DataType: double");
