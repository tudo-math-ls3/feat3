#include <test_system/test_system.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

int main(int argc, char** argv)
{
  int result(EXIT_SUCCESS);

  if(argc > 1)
  {
    std::list<String> labels;
    for(int i(1) ; i < argc ; ++i)
    {
      labels.push_back(argv[i]);
    }
    for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
        i != i_end ; )
    {
      if((find(labels.begin(), labels.end(), (*i)->get_memory_name()) == labels.end()) &&
          (find(labels.begin(), labels.end(), (*i)->get_prec_name()) == labels.end()))
      {
        i = TestList::instance()->erase(i);
        continue;
      }
      ++i;
    }
  }

  size_t list_size(TestList::instance()->size());
  size_t tests_passed(0u);
  size_t tests_failed(0u);
  unsigned long iterator_index(1);
  for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
      i != i_end ; )
  {
    CONTEXT("When running test case '" + (*i)->id() + ":");
    try
    {
      std::cout << "(" << iterator_index << "/" << list_size << ") " << (*i)->id()
        << " [Memory: " << (*i)->get_memory_name() << "]"
        << " [Algo: " << (*i)->get_algo_name() << "]"
        << " [Precision: "<< (*i)->get_prec_name() << "]"
        << std::endl;
      (*i)->run();
      std::cout << "PASSED" << std::endl;
      ++tests_passed;
    }
    catch (TestFailedException & e)
    {
      std::cout << "FAILED: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl;
      result = EXIT_FAILURE;
      ++tests_failed;
    }
    catch (InternalError & e)
    {
      std::cout << "FAILED with InternalError: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl
        << stringify(e.message()) << std::endl;
      result = EXIT_FAILURE;
      ++tests_failed;
    }
    catch (std::exception & e)
    {
      std::cout << "FAILED with unknown Exception: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl
        << std::endl;
      result = EXIT_FAILURE;
      ++tests_failed;
    }
    i = TestList::instance()->erase(i);
    iterator_index++;
  }

  for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
      i != i_end ; )
  {
    i = TestList::instance()->erase(i);
  }

  if(result == EXIT_SUCCESS)
  {
    std::cout << "All " << list_size << " tests PASSED!" << std::endl;
  }
  else
  {
    std::cout << tests_passed << " of " << list_size << " tests PASSED, "
      << tests_failed << " tests FAILED!" << std::endl;
  }

  return result;
}
