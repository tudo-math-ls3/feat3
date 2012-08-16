#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/topology.hpp>
#include <kernel/foundation/mapped_topology.hpp>
#include <kernel/foundation/dense_data_wrapper.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;


template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MappedTopologyTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MappedTopologyTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MappedTopologyTest<" + tag + ">")
    {
    }

    void run() const
    {
      Foundation::MappedTopology<std::string, Index, OT_, IT_> mt;
      mt.insert(1, "test");

      std::string b(mt.find(1));
      std::string c(mt.at(0));
    }
};
MappedTopologyTest<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > mappedtopology_test_cpu_v_v("std::vector, std::vector");
MappedTopologyTest<Archs::None, unsigned long, std::deque, std::vector<unsigned long> > mappedtopology_test_cpu_d_v("std::deque, std::vector");
MappedTopologyTest<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > mappedtopology_test_cpu_v_d("std::vector, std::deque");
MappedTopologyTest<Archs::None, unsigned long, std::deque, std::deque<unsigned long> > mappedtopology_test_cpu_d_d("std::deque, std::deque");

MappedTopologyTest<Archs::CPU, unsigned long, std::vector, Foundation::DenseDataWrapper<15, Archs::CPU, unsigned long, DenseVector> > mappedtopology_test_cpu_v_ddwdv("std::vector, DV");
