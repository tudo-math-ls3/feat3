#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/dense_data_wrapper.hpp>
#include <kernel/hornet/dense_vector.hpp>
#include <kernel/archs.hpp>

#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_,
  typename DataType_,
  template<typename, typename> class ContType_>
class DenseDataWrapperTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    DenseDataWrapperTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("DenseDataWrapperTest<" + tag + ">")
    {
    }

    void run() const
    {
      Foundation::DenseDataWrapper<15u, Tag_, DataType_, ContType_> test;
    }
};
DenseDataWrapperTest<Archs::CPU, double, DenseVector> ddw_test_DV("Hornet::DenseVector<double>");
