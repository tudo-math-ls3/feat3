#ifndef SERIAL
#define SERIAL

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/topology.hpp>
#include <kernel/foundation/dense_data_wrapper.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::LAFEM;


template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class TopologyTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    TopologyTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("TopologyTest<" + tag + ">")
    {
    }

    void run() const
    {
      Foundation::Topology<> t;

      t.push_back();

      TEST_CHECK_EQUAL(t.size(), 1ul);

      t.at(0).push_back(1123);
      t[0].push_back(878);
      TEST_CHECK_EQUAL(t.at(0).at(0), 1123ul);
      TEST_CHECK_EQUAL(t[0][1], 878ul);

      TEST_CHECK_EQUAL(t.size(), 1ul);

      Foundation::Topology<IndexType_, OT_, IT_> t2;
      t2.push_back();

      TEST_CHECK_EQUAL(t2.size(), 1ul);

      t2.at(0).push_back(1123);
      t2[0].push_back(878);
      TEST_CHECK_EQUAL(t2.at(0).at(0), 1123ul);
      TEST_CHECK_EQUAL(t2[0][1], 878ul);

      TEST_CHECK_EQUAL(t2.size(), 1ul);

      TEST_CHECK_EQUAL(t2.get_history().size(), 1ul);

      t2.erase(0);
      TEST_CHECK_EQUAL(t2.get_history().size(), 2ul);
      TEST_CHECK_EQUAL(t2.size(), 0ul);
    }
};
TopologyTest<Archs::None, Index, std::vector, std::vector<Index> > topology_test_cpu_v_v("std::vector, std::vector");
TopologyTest<Archs::None, Index, std::deque, std::vector<Index> > topology_test_cpu_d_v("std::deque, std::vector");
TopologyTest<Archs::None, Index, std::vector, std::deque<Index> > topology_test_cpu_v_d("std::vector, std::deque");
TopologyTest<Archs::None, Index, std::deque, std::deque<Index> > topology_test_cpu_d_d("std::deque, std::deque");

TopologyTest<Mem::Main, Index, std::vector, Foundation::DenseDataWrapper<15, Mem::Main, Index, DenseVector> > topology_test_cpu_v_ddwdv("std::vector, DV");

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class TopologyCommTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    TopologyCommTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("TopologyCommTest<" + tag + ">")
    {
    }

    void run() const
    {
      Foundation::Topology<> t;
      t.push_back();
      t.at(0).push_back(42);
      t.at(0).push_back(47);
      t.push_back();
      t.at(1).push_back(52);
      t.at(1).push_back(57);

      Foundation::Topology<>::buffer_type_ sendbuf(t.buffer());
      Foundation::Topology<>::buffer_type_ recvbuf(t.buffer());

      t.to_buffer(sendbuf);

      t.send_recv(
          sendbuf,
          0,
          recvbuf,
          0);

      t.from_buffer(recvbuf);

      TEST_CHECK_EQUAL(t.at(0).at(0), 42);
      TEST_CHECK_EQUAL(t.at(0).at(1), 47);
      TEST_CHECK_EQUAL(t.at(1).at(0), 52);
      TEST_CHECK_EQUAL(t.at(1).at(1), 57);
    }
};
/*TopologyCommTest<Archs::None, Index, std::vector, std::vector<Index> > topology_commtest_cpu_v_v("std::vector, std::vector");
TopologyCommTest<Archs::None, Index, std::deque, std::vector<Index> > topology_commtest_cpu_d_v("std::deque, std::vector");
TopologyCommTest<Archs::None, Index, std::vector, std::deque<Index> > topology_commtest_cpu_v_d("std::vector, std::deque");
TopologyCommTest<Archs::None, Index, std::deque, std::deque<Index> > topology_commtest_cpu_d_d("std::deque, std::deque");

TopologyCommTest<Mem::Main, Index, std::vector, Foundation::DenseDataWrapper<15, Mem::Main, Index, DenseVector> > topology_commtest_cpu_v_ddwdv("std::vector, DV");*/
#endif // SERIAL
