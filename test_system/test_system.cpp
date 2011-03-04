// includes, system
#include <cstdlib>
#include <iostream>
#include <algorithm>

// includes, FEAST
#include <test_system/test_system.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/exception.hpp>

using namespace FEAST;
using namespace TestSystem;

int main(int argc, char** argv)
{
  int result(EXIT_SUCCESS);

  std::string mpc("mpiproccount");
  if(argc == 2 && std::string(argv[1]) == mpc.c_str())
  {
    return (*TestList::instance()->begin_tests())->mpi_proc_count();
  }

  if(argc > 1)
  {
    for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
        i != i_end ; ++i)
    {
      if((*i)->mpi_proc_count() != (*TestList::instance()->begin_tests())->mpi_proc_count())
      {
        std::cout << "mpi_proc_count missmatch!"<<std::endl;
        result = EXIT_FAILURE;
        return (result);
      }
    }

    std::list<std::string> labels;
    for(int i(1) ; i < argc ; ++i)
    {
      labels.push_back(argv[i]);
    }
    for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
        i != i_end ; )
    {
      if((find(labels.begin(), labels.end(), (*i)->get_tag_name()) == labels.end()) &&
         (find(labels.begin(), labels.end(), (*i)->get_prec_name()) == labels.end()))
      {
        i = TestList::instance()->erase(i);
        continue;
      }
      ++i;
    }
  }

  size_t list_size(TestList::instance()->size());
  unsigned long iterator_index(1);
  for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
      i != i_end ; )
  {
    CONTEXT("When running test case '" + (*i)->id() + "':");
    try
    {
      std::cout << "(" << iterator_index << "/" << list_size << ") " << (*i)->id() + " [Backend: "
        << (*i)->get_tag_name() << "]" << " [Precision: "<< (*i)->get_prec_name() << "]" << std::endl;
      (*i)->run();
      std::cout << "PASSED" << std::endl;
    }
    catch (TestFailedException & e)
    {
      std::cout << "FAILED: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl;
      result = EXIT_FAILURE;
    }
    catch (InternalError & e)
    {
      std::cout << "FAILED: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl
                << stringify(e.message()) << std::endl;
      result = EXIT_FAILURE;
    }
    i = TestList::instance()->erase(i);
    iterator_index++;
  }

  return result;
}
