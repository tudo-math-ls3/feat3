#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/archs.hpp>
#include<deque>
#include <kernel/foundation/functor.hpp>
#include <kernel/foundation/topology.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_= Mem::Main, typename IndexType_ = Index>
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
      //Testing with simple stl container
      std::vector<IndexType_> vector_0;
      std::vector<std::vector<IndexType_> > vector_vector_0;

      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_0(vector_0, IndexType_(0), IndexType_(42));
      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_1(vector_0, IndexType_(1), IndexType_(43));
      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_2(vector_0, IndexType_(2), IndexType_(44));

      TEST_CHECK_THROWS(func_0.undo(), Foundation::FunctorError);
      TEST_CHECK_THROWS(func_1.undo(), Foundation::FunctorError);
      TEST_CHECK_THROWS(func_2.undo(), Foundation::FunctorError);
      func_0.execute();
      func_1.execute();
      func_2.execute();
      TEST_CHECK_THROWS(func_0.execute(), Foundation::FunctorError);
      TEST_CHECK_THROWS(func_1.execute(), Foundation::FunctorError);
      TEST_CHECK_THROWS(func_2.execute(), Foundation::FunctorError);
      TEST_CHECK_EQUAL(vector_0.at(0), IndexType_(42));
      TEST_CHECK_EQUAL(vector_0.at(1), IndexType_(43));
      TEST_CHECK_EQUAL(vector_0.at(2), IndexType_(44));

      Foundation::EraseFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_3(vector_0, IndexType_(0), IndexType_(vector_0.at(0)));
      Foundation::EraseFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_4(vector_0, IndexType_(1), IndexType_(vector_0.at(1)));
      Foundation::EraseFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_5(vector_0, IndexType_(2), IndexType_(vector_0.at(2)));
      TEST_CHECK_THROWS(func_5.undo(), Foundation::FunctorError);
      TEST_CHECK_THROWS(func_4.undo(), Foundation::FunctorError);
      TEST_CHECK_THROWS(func_3.undo(), Foundation::FunctorError);
      func_5.execute();
      func_4.execute();
      func_3.execute();
      TEST_CHECK_THROWS(func_3.execute(), Foundation::FunctorError);
      TEST_CHECK_THROWS(func_4.execute(), Foundation::FunctorError);
      TEST_CHECK_THROWS(func_5.execute(), Foundation::FunctorError);
      TEST_CHECK_EQUAL(vector_0.size(), IndexType_(0));

      //testing with complex STL container
      std::vector<IndexType_> vector_1;
      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_6(vector_1, IndexType_(0), IndexType_(42));
      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_7(vector_1, IndexType_(1), IndexType_(43));
      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_8(vector_1, IndexType_(2), IndexType_(44));
      func_6.execute();
      func_7.execute();
      func_8.execute();

      Foundation::PushBackFunctor<std::vector<std::vector<IndexType_> >, IndexType_, std::vector<IndexType_> > func_9(vector_vector_0, IndexType_(0), vector_1);
      func_9.execute();
      TEST_CHECK_EQUAL(vector_vector_0.at(0).at(0), IndexType_(42));
      TEST_CHECK_EQUAL(vector_vector_0.at(0).at(1), IndexType_(43));
      TEST_CHECK_EQUAL(vector_vector_0.at(0).at(2), IndexType_(44));
      Foundation::EraseFunctor<std::vector<std::vector<IndexType_> >, IndexType_, std::vector<IndexType_> > func_10(vector_vector_0, IndexType_(0), vector_vector_0.at(0));
      func_10.execute();
      TEST_CHECK_EQUAL(vector_vector_0.size(), IndexType_(0));

      //testing compound functor with simple STL container
      Foundation::CompoundFunctor<> cfunc_0;
      cfunc_0.add_functor(new Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_>(vector_0, IndexType_(0), IndexType_(42)));
      cfunc_0.add_functor(new Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_>(vector_0, IndexType_(1), IndexType_(43)));
      cfunc_0.add_functor(new Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_>(vector_0, IndexType_(2), IndexType_(44)));
      cfunc_0.execute();
      TEST_CHECK_EQUAL(vector_0.at(0), IndexType_(42));
      TEST_CHECK_EQUAL(vector_0.at(1), IndexType_(43));
      TEST_CHECK_EQUAL(vector_0.at(2), IndexType_(44));
      cfunc_0.undo();
      TEST_CHECK_EQUAL(vector_0.size(), IndexType_(0));

      //testing compound functor with complex STL container
      Foundation::CompoundFunctor<> cfunc_1;
      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_11(vector_0, IndexType_(0), IndexType_(42));
      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_12(vector_0, IndexType_(1), IndexType_(43));
      Foundation::PushBackFunctor<std::vector<IndexType_>, IndexType_, IndexType_> func_13(vector_0, IndexType_(2), IndexType_(44));
      func_11.execute();
      func_12.execute();
      func_13.execute();
      cfunc_1.add_functor(new Foundation::PushBackFunctor<std::vector<std::vector<IndexType_> >, IndexType_, std::vector<IndexType_> >(vector_vector_0, IndexType_(0), vector_0));
      cfunc_1.execute();
      TEST_CHECK_EQUAL(vector_vector_0.at(0).at(0), IndexType_(42));
      TEST_CHECK_EQUAL(vector_vector_0.at(0).at(1), IndexType_(43));
      TEST_CHECK_EQUAL(vector_vector_0.at(0).at(2), IndexType_(44));
      cfunc_1.undo();
      TEST_CHECK_EQUAL(vector_vector_0.size(), IndexType_(0));

    }
};
FunctorTest<> simple_functor_test("None, Index");
