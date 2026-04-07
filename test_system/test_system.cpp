// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/runtime.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/time_stamp.hpp>

#include <cstring>

using namespace FEAT;
using namespace FEAT::TestSystem;

int main(int argc, char** argv)
{
  // Initialse FEAT runtime
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  std::cout << "CTEST_FULL_OUTPUT\n";

  int result(EXIT_SUCCESS);

  TimeStamp t_start;

  if(argc > 1)
  {
    bool all_filter(false);
    std::list<String> labels;
    for(int i(1) ; i < argc ; ++i)
    {
      if (0 == strcmp(argv[i], "and"))
      {
        all_filter = true;
        continue;
      }

      //discard any unused parameters, marked by "--"
      if (strlen(argv[i]) > 1 && argv[i][0] == '-' && argv[i][1] == '-')
        continue;

      labels.push_back(argv[i]);
    }

    for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
        i != i_end ; )
    {
      // any_filter: run test if any one label matches
      if (all_filter == false)
      {
        if((find(labels.begin(), labels.end(), (*i)->get_preferred_backend_name()) == labels.end()) &&
            (find(labels.begin(), labels.end(), (*i)->get_datatype_name()) == labels.end()) &&
            (find(labels.begin(), labels.end(), (*i)->get_index_name()) == labels.end()) )
        {
          i = TestList::instance()->erase(i);
          continue;
        }
        ++i;
      }
      // all_filter: run test if every single label matches
      else
      {
        bool deleted(false);
        for (auto li(labels.begin()) ; li != labels.end() ; ++li)
        {
          if(*li != (*i)->get_preferred_backend_name() &&
              *li != (*i)->get_datatype_name() &&
              *li != (*i)->get_index_name())
          {
            i = TestList::instance()->erase(i);
            deleted = true;
            break;
          }
        }
        if (!deleted)
          ++i;
      }
    }
  }

  size_t list_size(TestList::instance()->size());
  size_t tests_passed(0u);
  size_t tests_failed(0u);
  unsigned long iterator_index(1);
  for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
      i != i_end ; )
  {
    try
    {
      std::cout << "(" << iterator_index << "/" << list_size << ") " << (*i)->id()
        << " [Preferred Backend: " << (*i)->get_preferred_backend_name() << "]"
        << " [Precision: "<< (*i)->get_datatype_name() << "]"
        << " [Indexing: "<< (*i)->get_index_name() << "]"
        << "\n";
      Backend::set_preferred_backend((*i)->get_preferred_backend());
      TimeStamp t1;
      (*i)->run();
      TimeStamp t2;
      std::cout << "PASSED in " << t2.elapsed_string(t1) <<" seconds\n";
      ++tests_passed;
    }
    catch (TestFailedException & e)
    {
      std::cout << "FAILED: " << (*i)->id() << "\n" << stringify(e.what()) << "\n";
      result = EXIT_FAILURE;
      ++tests_failed;
    }
    /*
    catch(std::exception& e)
    {
      std::cout << "CRASHED: " << (*i)->id() << "\n" << stringify(e.what()) << "\n";
      result = EXIT_FAILURE;
      ++tests_failed;
    }
    catch(...)
    {
      std::cout << "CRASHED: " << (*i)->id() << "\ncaught unknown exception\n";
      result = EXIT_FAILURE;
      ++tests_failed;
    }
    */

    i = TestList::instance()->erase(i);
    iterator_index++;
  }

  for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
      i != i_end ; )
  {
    i = TestList::instance()->erase(i);
  }

  std::cout << "\n";
  if(result == EXIT_SUCCESS)
  {
    std::cout << "All " << list_size << " tests PASSED!\n";
  }
  else
  {
    std::cout << tests_passed << " of " << list_size << " tests PASSED, "
      << tests_failed << " tests FAILED!\n";
  }

  TimeStamp t_end;
  std::cout << "Total Runtime: " << t_end.elapsed_string(t_start) << " seconds\n";

  return result;
}
