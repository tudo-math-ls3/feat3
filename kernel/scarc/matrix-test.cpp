#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/matrix.hpp>
#include <kernel/util/cpp11_smart_pointer.hpp>
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
      ScaRC::Matrix<> m;
    }
};
MatrixTest<Archs::None, unsigned long> attribute_test_cpu_v_ulong_float("StorageType: std::vector, DataType: ulong");
