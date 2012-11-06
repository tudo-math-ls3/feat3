#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/vector.hpp>
#include <kernel/util/cpp11_smart_pointer.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_, typename DataType_>
class VectorTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    VectorTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("VectorTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      ScaRC::DynamicVector<> v1(2000);
      ScaRC::DynamicVector<> v2(v1);

      v2.add_block(&v1, 100);

      TEST_CHECK_EQUAL(v2.size(), 2000ul);
      TEST_CHECK_EQUAL(v2.num_blocks(), 1ul);
      TEST_CHECK_EQUAL(v1.num_blocks(), 0ul);
    }
};
VectorTest<Archs::None, unsigned long> attribute_test_cpu_v_ulong_float("StorageType: std::vector, DataType: ulong");
