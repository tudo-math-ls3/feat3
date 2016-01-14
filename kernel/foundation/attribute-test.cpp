#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/attribute.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_, typename DataType1_, typename DataType2_, template<typename, typename> class ST_>
class AttributeTest:
  public TaggedTest<Tag_, DataType1_>
{
  public:
    AttributeTest(const std::string & tag) :
      TaggedTest<Tag_, DataType1_>("AttributeTest<" + tag + ">")
    {
    }

    virtual void run() const override
    {
      std::vector<std::shared_ptr<Foundation::AttributeBase<ST_> > > attrs;

      attrs.push_back(std::shared_ptr<Foundation::AttributeBase<ST_> >(new Foundation::Attribute<DataType1_, ST_>));
      attrs.push_back(std::shared_ptr<Foundation::AttributeBase<ST_> >(new Foundation::Attribute<DataType2_, ST_>));
      for(Index j(0) ; j < 1000 ; ++j)
      {
        ((Foundation::Attribute<DataType1_, ST_>*)(attrs.at(0).get()))->get_data().push_back(DataType1_(DataType1_(5) + j));
        ((Foundation::Attribute<DataType2_, ST_>*)(attrs.at(1).get()))->get_data().push_back(DataType2_(DataType2_(5) * j));
      }

      for(Index j(0) ; j < 1000 ; ++j)
      {
        TEST_CHECK_EQUAL(((Foundation::Attribute<DataType1_, ST_>*)attrs.at(0).get())->at(j), DataType1_(5 + j));
        TEST_CHECK_EQUAL(((Foundation::Attribute<DataType2_, ST_>*)attrs.at(1).get())->at(j), DataType2_(5 * j));
        TEST_CHECK_EQUAL((*((Foundation::Attribute<DataType1_, ST_>*)attrs.at(0).get()))[j], DataType1_(5 + j));
        TEST_CHECK_EQUAL((*((Foundation::Attribute<DataType2_, ST_>*)attrs.at(1).get()))[j], DataType2_(5 * j));
      }
    }
};
AttributeTest<Mem::Main, Index, float, std::vector> attribute_test_cpu_v_ulong_float("StorageType: std::vector, DataTypes: ulong, float");
AttributeTest<Mem::Main, double, float, std::vector> attribute_test_cpu_v_double_float("StorageType: std::vector, DataTypes: double, float");
AttributeTest<Mem::Main, double, int, std::vector> attribute_test_cpu_v_double_int("StorageType: std::vector, DataTypes: double, int");
AttributeTest<Mem::Main, Index, float, std::deque> attribute_test_cpu_d_ulong_float("StorageType: std::deque, DataTypes: ulong, float");
AttributeTest<Mem::Main, double, float, std::deque> attribute_test_cpu_d_double_float("StorageType: std::deque, DataTypes: double, float");
AttributeTest<Mem::Main, double, int, std::deque> attribute_test_cpu_d_double_int("StorageType: std::deque, DataTypes: double, int");
