#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/dense_data_wrapper.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/archs.hpp>
#include <algorithm>

#include<deque>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;


template<typename Tag_,
  typename DataType_,
  typename Algo_,
  template<typename, typename, typename> class ContType_>
class DenseDataWrapperTest:
  public TaggedTest<Tag_, DataType_, Algo_>
{
  public:
    DenseDataWrapperTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_, Algo_>("DenseDataWrapperTest<" + tag + ">")
    {
    }

    void run() const
    {
      Foundation::DenseDataWrapper<15u, Tag_, DataType_, ContType_> test;
      test.push_back(3); test.push_back(1); test.push_back(0); test.push_back(2);
      std::sort(test.begin(), test.end());

      TEST_CHECK_EQUAL(test.at(0), DataType_(0));
      TEST_CHECK_EQUAL(test.at(1), DataType_(1));
      TEST_CHECK_EQUAL(test.at(2), DataType_(2));
      TEST_CHECK_EQUAL(test.at(3), DataType_(3));

      test.erase(test.begin() + 3);
      TEST_CHECK_EQUAL(test.size(), 3ul);
      TEST_CHECK_EQUAL(test.at(0), DataType_(0));
      TEST_CHECK_EQUAL(test.at(1), DataType_(1));
      TEST_CHECK_EQUAL(test.at(2), DataType_(2));

      test.insert(test.begin() + 1, DataType_(118));
      TEST_CHECK_EQUAL(test.size(), 4ul);
      TEST_CHECK_EQUAL(test.at(0), DataType_(0));
      TEST_CHECK_EQUAL(test.at(1), DataType_(118));
      TEST_CHECK_EQUAL(test.at(2), DataType_(1));
      TEST_CHECK_EQUAL(test.at(3), DataType_(2));


      std::vector<Foundation::DenseDataWrapper<15u, Tag_, DataType_, ContType_> > test1;
      test1.push_back(std::move(test));

      std::vector<Foundation::DenseDataWrapper<15u, Tag_, DataType_, ContType_> > test2(std::move(test1));
    }
};
DenseDataWrapperTest<Mem::Main, double, Algo::Generic, DenseVector> ddw_test_DV("DenseVector<double>");
