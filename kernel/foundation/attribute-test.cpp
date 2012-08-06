#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/attribute.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_, typename DataType_, template<typename, typename> class ST_>
class AttributeTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    AttributeTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("AttributeTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      Foundation::Attribute<DataType_, ST_> a(2000);
      Foundation::AttributeBase<ST_> * ab;
    }
};
AttributeTest<Archs::None, unsigned long, std::vector> attribute_test_cpu_v_v("std::vector");
