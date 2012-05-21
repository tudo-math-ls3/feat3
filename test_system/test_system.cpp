#include <test_system/test_system.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

#ifdef PARALLEL
/* ********************************************************************************************* */
/*   -  P A R A L L E L   T E S T   D R I V E R  *  P A R A L L E L   T E S T   D R I V E R   -  */
/* ********************************************************************************************* */
int main(int argc, char** argv)
{
  int result(EXIT_SUCCESS);
  init_mpi(argc, argv);
  int rank(-1);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc == 2)
  {
    String mpc("mpiproccount");
    String argv1(argv[1]);
    if(0 == mpc.compare(argv1))
    {
      return (*TestList::instance()->begin_tests())->mpi_proc_count();
    }
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

    std::list<String> labels;
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
    CONTEXT("When running test case '" + (*i)->id() + "' on mpi process " + stringify(rank) + ":");
    try
    {
      if (rank == 0)
        std::cout << "(" << iterator_index << "/" << list_size << ") " << (*i)->id() + " [Backend: "
          << (*i)->get_tag_name() << "]" << " [Precision: "<< (*i)->get_prec_name() << "]" << std::endl;
      (*i)->run();
      if (rank == 0)
        std::cout << "PASSED" << std::endl;
    }
    catch (TestFailedException & e)
    {
      if (rank == 0)
        std::cout << "FAILED: " << (*i)->id() << std::endl << stringify(e.what()) << " on mpi process " << stringify(rank) << std::endl;
      result = EXIT_FAILURE;
    }
    catch (InternalError & e)
    {
        std::cout << "FAILED with InternalError: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl
          << stringify(e.message()) << " on mpi process " << stringify(rank) << std::endl;
      result = EXIT_FAILURE;
    }
    catch (std::exception & e)
    {
      std::cout << "FAILED with unknown Exception: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl
        << " on mpi process " << stringify(rank) << std::endl;
      result = EXIT_FAILURE;
    }
    i = TestList::instance()->erase(i);
    iterator_index++;
  }

  for(TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
      i != i_end ; )
  {
    i = TestList::instance()->erase(i);
  }

  if (result == EXIT_SUCCESS && rank == 0)
    std::cout<<"All " << list_size << " tests PASSED!" << std::endl;
  finalise_mpi();
  return result;
}
#else
/* ************************************************************************************* */
/*   -  S E R I A L   T E S T   D R I V E R  *  S E R I A L   T E S T   D R I V E R   -  */
/* ************************************************************************************* */
int main(int /*argc*/, char** /*argv*/)
{
  int result(EXIT_SUCCESS);

  size_t list_size = TestList::instance()->size();
  unsigned long iterator_index(1);

  TestList::Iterator it(TestList::instance()->begin_tests());
  TestList::Iterator jt(TestList::instance()->end_tests());

  // Check whether the test list is empty; this might happen if a test can only run in parallel mode.
  if(it == jt)
  {
    std::cout << "PASSED" << std::endl;
    return result;
  }

  while(it != jt)
  {
    try
    {
      std::cout << "(" << iterator_index << "/" << list_size << ") " << (*it)->id() + " [Backend: "
        << (*it)->get_tag_name() << "]" << " [Precision: "<< (*it)->get_prec_name() << "]" << std::endl;
      (*it)->run();
      std::cout << "PASSED" << std::endl;
    }
    catch (TestFailedException & e)
    {
      std::cout << "FAILED: " << (*it)->id() << std::endl << stringify(e.what()) << std::endl;
      result = EXIT_FAILURE;
    }
    it = TestList::instance()->erase(it);
    ++iterator_index;
  }

  if (result == EXIT_SUCCESS)
    std::cout<<"All " << list_size << " tests PASSED!" << std::endl;

  return result;
}
#endif // PARALLEL
