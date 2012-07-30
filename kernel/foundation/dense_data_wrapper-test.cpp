#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/dense_data_wrapper.hpp>
#include <kernel/archs.hpp>

#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;

//Test container
template<typename DT_>
class TestArrayClass
{
  public:
    TestArrayClass(Index size) :
      _size(size),
      _data(new DT_[size])
    {
    }

    ~TestArrayClass()
    {
      delete[] _data;
    }

    DT_ & operator[] (Index i)
    {
      return _data[i];
    }

    Index size()
    {
      return _size;
    }

  private:
    Index _size;
    DT_ * _data;

};

template<typename Tag_, typename DataType_, template<typename> class ContType_ = TestArrayClass>
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
      Foundation::DenseDataWrapper<15u, DataType_, ContType_> test;
    }
};
DenseDataWrapperTest<Archs::None, double> ddw_test("Array<double>");
