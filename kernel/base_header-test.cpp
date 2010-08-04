#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

using namespace TestSystem;
using namespace Feast;

template <typename Tag_, typename DT_>
class BaseHeaderTest :
  public TaggedTest<Tag_, DT_>
{
  public:
    BaseHeaderTest(const std::string & id) :
      TaggedTest<Tag_, DT_>(id)
    {
    }

    virtual void run() const
    {
      TEST_CHECK_EQUAL(nullptr, NULL);
    }
};
BaseHeaderTest<Nil, Nil> base_header_test ("BaseHeaderTest");
