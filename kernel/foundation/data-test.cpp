#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/data.hpp>
#include<kernel/archs.hpp>

#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class PDataTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    PDataTest(const std::string & tag) :
      TaggedTest<Tag_, Index>("PDataTest<" + tag + ">")
    {
    }

    virtual void run() const override
    {
      PData<
        Dim1D,
        Topology<IndexType_, OT_, IT_>,
        OT_,
        Mesh,
        double> p;
    }
};
PDataTest<Mem::Main, Index, std::vector, std::vector<Index> > ref_test_cpu_v_v("std::vector, std::vector");
PDataTest<Mem::Main, Index, std::vector, std::deque<Index> > ref_test_cpu_v_d("std::vector, std::deque");
PDataTest<Mem::Main, Index, std::deque, std::vector<Index> > ref_test_cpu_d_v("std::deque, std::vector");
PDataTest<Mem::Main, Index, std::deque, std::deque<Index> > ref_test_cpu_d_d("std::deque, std::deque");
