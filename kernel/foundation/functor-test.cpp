#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/archs.hpp>
#include<deque>
#include <kernel/foundation/functor.hpp>
#include <kernel/foundation/topology.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_= Archs::None, typename IndexType_ = Index>
class FunctorTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    FunctorTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("FunctorTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      std::vector<IndexType_> vector;
      vector.push_back(IndexType_(0));
      vector.push_back(IndexType_(42));
      vector.push_back(IndexType_(3));

      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func(vector, IndexType_(2), IndexType_(3));
      TEST_CHECK_EQUAL(func.get_position(), IndexType_(2));
      TEST_CHECK_EQUAL(func.get_value(), IndexType_(3));

      TEST_CHECK_THROWS(func.execute(), FunctorError);

      TEST_CHECK_EQUAL(vector.at(IndexType_(2)), IndexType_(3));

      func.undo();

      TEST_CHECK_EQUAL(vector.size(), 2);

      TEST_CHECK_THROWS(func.undo(), FunctorError);

      //-------------------------------------------
      Foundation::Topology<> topology;
      topology.push_back();
      topology.push_back();
      topology.push_back();
      topology.erase();
      topology.erase();
      topology.erase();

      TEST_CHECK_EQUAL(topology.get_history().size(), 6ul);
      TEST_CHECK_EQUAL(topology.get_history().at(0)->name(), "push_back()");
      TEST_CHECK_EQUAL(topology.get_history().at(1)->name(), "push_back()");
      TEST_CHECK_EQUAL(topology.get_history().at(2)->name(), "push_back()");
      TEST_CHECK_EQUAL(topology.get_history().at(3)->name(), "erase()");
      TEST_CHECK_EQUAL(topology.get_history().at(4)->name(), "erase()");
      TEST_CHECK_EQUAL(topology.get_history().at(5)->name(), "erase()");
    }
};
FunctorTest<> simple_functor_test("None, Index");
