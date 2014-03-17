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


template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class TopologyTest:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    TopologyTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_, Algo_>("TopologyTest<" + tag + ">")
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

      //TEST_CHECK_EQUAL(t2.get_history().size(), 1ul);

      t2.erase(0);
      //TEST_CHECK_EQUAL(t2.get_history().size(), 2ul);
      TEST_CHECK_EQUAL(t2.size(), 0ul);


      Foundation::Topology<IndexType_, OT_, IT_> t3;
      t3.push_back();
      t3.at(0).push_back(0);
      t3.at(0).push_back(3);
      t3.at(0).push_back(1);
      t3.at(0).push_back(4);

      t3.at(0).erase(t3.at(0).begin() + 2);
      TEST_CHECK_EQUAL(t3.at(0).size(), 3ul);
      TEST_CHECK_EQUAL(t3.at(0).at(0), 0ul);
      TEST_CHECK_EQUAL(t3.at(0).at(1), 3ul);
      TEST_CHECK_EQUAL(t3.at(0).at(2), 4ul);

      t3.at(0).insert(t3.at(0).begin() + 2, IndexType_(5));
      TEST_CHECK_EQUAL(t3.at(0).size(), 4ul);
      TEST_CHECK_EQUAL(t3.at(0).at(0), 0ul);
      TEST_CHECK_EQUAL(t3.at(0).at(1), 3ul);
      TEST_CHECK_EQUAL(t3.at(0).at(2), 5ul);
      TEST_CHECK_EQUAL(t3.at(0).at(3), 4ul);
    }
};
TopologyTest<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > topology_test_cpu_v_v("std::vector, std::vector");
TopologyTest<Mem::Main, Index, Algo::Generic, std::deque, std::vector<Index> > topology_test_cpu_d_v("std::deque, std::vector");
TopologyTest<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > topology_test_cpu_v_d("std::vector, std::deque");
TopologyTest<Mem::Main, Index, Algo::Generic, std::deque, std::deque<Index> > topology_test_cpu_d_d("std::deque, std::deque");
TopologyTest<Mem::Main, Index, Algo::Generic, std::vector, Foundation::DenseDataWrapper<15, Mem::Main, Index, DenseVector> > topology_test_cpu_v_ddwdv("std::vector, DV");
