#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/buffer.hpp>
#include <kernel/archs.hpp>
#include <vector>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_, typename DataType_>
class BufferTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    BufferTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("BufferTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      std::vector<std::shared_ptr<Foundation::SharedArrayBase > > buffers;
      buffers.push_back(Foundation::BufferedSharedArray<DataType_>::create(20));

      for(Index i(0) ; i < 20 ; ++i)
        (*(Foundation::BufferedSharedArray<DataType_>*)((buffers.at(0).get())))[i] = i;

      for(Index i(0) ; i < 20 ; ++i)
        TEST_CHECK_EQUAL((*(Foundation::BufferedSharedArray<DataType_>*)((buffers.at(0).get())))[i], i);

      //------------------------------------------------------------------------------------------------
      Foundation::BufferedData<> b;
      b.get().push_back(Foundation::BufferedSharedArray<DataType_>::create(20));
      for(Index i(0) ; i < 20 ; ++i)
        (*(Foundation::BufferedSharedArray<DataType_>*)((b.get().at(0).get())))[i] = i;

      for(Index i(0) ; i < 20 ; ++i)
        TEST_CHECK_EQUAL((*(Foundation::BufferedSharedArray<DataType_>*)((b.get().at(0).get())))[i], i);
    }
};
BufferTest<Archs::None, Index> buffer_test_cpu_ulong("ulong");
