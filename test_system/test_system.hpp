#pragma once
#ifndef TEST_SYSTEM_TEST_SYSTEM_HPP
/// Header guard
#define TEST_SYSTEM_TEST_SYSTEM_HPP 1

// includes, system
#include <string>
#include <exception>
#include <list>
#include <typeinfo>
#include <cstdlib>
#include <iostream>
#include <algorithm>

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/mpi_utils.hpp>

/**
* \file
*
* Implementation of Test and related classes.
*/

/// TestSystem namespace
namespace FEAST
{
  namespace TestSystem
  {
    // Forwared declaration
    class BaseTest;

    /// exception thrown by the check method in BaseTest
    class TestFailedException
      : public std::exception
    {

      private:

        /// failure message
        const std::string _message;


      public:

        /**
         * \brief CTOR
         *
         * \param[in] function
         * The function that trows the exception.
         *
         * \param[in] file
         * The file the exception was thrown in.
         *
         * \param[in] line
         * The line the exception was thrown in.
         *
         * \param[in] message
         * A message describing the exception.
         */
        TestFailedException(
            const char* const function,
            const char* const file,
            const long line,
            const std::string & message) throw ()
          : _message(stringify(file) + ":" + stringify(line) + ": in " + stringify(function) + ": " + message )
        {
        }

        /// DTOR
        virtual ~TestFailedException() throw ()
        {
        }

        /// description
        const char* what() const throw ()
        {
          return _message.c_str();
        }
    };

    /// list of all instantiated tests
    class TestList
      : public InstantiationPolicy<TestList, Singleton>
    {

      private:

        /// internal STL list representation of TestList
        std::list<BaseTest*> _tests;

        /// default CTOR
        TestList()
        {
        }


      public:

        friend TestList* InstantiationPolicy<TestList, Singleton>::instance();

        /// TestList forward iterator.
        typedef std::list<BaseTest*>::iterator Iterator;

        /**
         * \brief adds a test to the TestList
         *
         * \param[in] test
         * the test that will be added
         */
        void register_test(BaseTest* const test)
        {
          _tests.push_back(test);
        }

        /// returns an iterator to the front of the TestList
        Iterator begin_tests()
        {
          return _tests.begin();
        }

        /// returns an iterator beyond the last element of the TestList
        Iterator end_tests()
        {
          return _tests.end();
        }

        /// returns the size of the TestList
        size_t size()
        {
          return _tests.size();
        }

        /// removes iterator target from the TestList
        Iterator erase (Iterator i)
        {
          return _tests.erase(i);
        }
    };


    /**
     * \brief base class for all Tests
     *
     * \author Dirk Ribbrock
     */
    class BaseTest
    {

      protected:

        /// test description String
        const std::string _id;
        /// architecture description String
        std::string _tag_name;
        /// precision description String
        std::string _prec_name;


      public:

        /**
         * \brief CTOR
         *
         * \param[in] id
         * the testcase's id string
         */
        BaseTest(const std::string& id)
          : _id(id),
          _tag_name(TypeTraits<Nil>::name()),
          _prec_name(TypeTraits<Nil>::name())
      {
        TestList::instance()->register_test(this);
      }

        /// DTOR
        virtual ~BaseTest() {}

        /// returns our id string
        virtual const std::string id() const
        {
          return _id;
        }

        /// returns the mpi proc count to use
        virtual unsigned long mpi_proc_count() const
        {
          return 1;
        }

        /// utility method used bei TEST_CHECK_*
        virtual void check(const char * const function, const char * const file,
            const long line, bool was_ok, const std::string & message) const
        {
          if(! was_ok)
            throw TestFailedException(function, file, line, message);
        }

        /**
         * \brief runs the test case
         *
         * Called by unittest framework only.
         */
        virtual void run() const = 0;

        /// returns our target platform
        virtual std::string get_tag_name()
        {
          return _tag_name;
        }

        /// returns our target platform
        virtual std::string get_prec_name()
        {
          return _prec_name;
        }

        /// utility class used by TEST_CHECK_EQUAL
        struct TwoVarHolder
        {
          /// result of comparison
          bool result;
          /// string representation of first parameter
          std::string s_a;
          /// string representation of second parameter
          std::string s_b;

          /**
           * \brief CTOR
           *
           * \param[in] a
           * First value to compare
           *
           * \param[in] b
           * Second value to compare
           */
          template<
            typename T1_,
                     typename T2_>
                       TwoVarHolder(
                           T1_ a,
                           T2_ b)
                       : result(a == b),
                       s_a(),
                       s_b()
          {
            if(!result)
            {
              s_a = stringify(a);
              s_b = stringify(b);
            }
          }
        };

        /// utility class used by TEST_CHECK_NOT_EQUAL
        struct TwoVarHolder2
        {
          /// result of comparison
          bool result;
          /// string representation of first parameter
          std::string s_a;
          /// string representation of second parameter
          std::string s_b;

          /**
           * \brief Constructor
           *
           * \param[in] a
           * First value to compare
           *
           * \param[in] b
           * Second value to compare
           */
          template<
            typename T1_,
                     typename T2_>
                       TwoVarHolder2(
                           T1_ a,
                           T2_ b)
                       : result(a == b),
                       s_a(),
                       s_b()
          {
            if(result)
            {
              s_a = stringify(a);
              s_b = stringify(b);
            }
          }
        };

        /// utility class used by TEST_CHECK_EQUAL_WITHIN_EPS.
        struct WithinEpsCalculator
        {
          /// result of comparison
          bool result;
          /// string representation of first parameter
          std::string s_a;
          /// string representation of second parameter
          std::string s_b;
          /// string representation of eps
          std::string s_diff;

          /**
           * \brief Constructor
           *
           * \param[in] a
           * first value to compare
           *
           * \param[in] b
           * second value to compare
           *
           * \param[in] c
           * maximum epsilon the values a and b are allowed to differ from each other
           */
          template<
            typename T1_,
                     typename T2_,
                     typename T3_>
                       WithinEpsCalculator(
                           T1_ a,
                           T2_ b,
                           T3_ c)
                       : s_a(),
                       s_b()
          {
            if(a >= b)
            {
              result = ((a - b) <= c);
              if(!result)
              {
                s_diff = stringify(a - b);
                s_a = stringify(a);
                s_b = stringify(b);
              }
            }
            else
            {
              result = ((b - a) <= c);
              if(!result)
              {
                s_diff = stringify(b - a);
                s_a = stringify(a);
                s_b = stringify(b);
              }
            }
          }
        };
    };

    /**
     * \brief abstract base class for all tagged test classes
     *
     * \author Dirk Ribbrock
     */
    template<
      typename Tag_,
               typename DataType_>
                 class TaggedTest
                 : public BaseTest
                 {
                   public:
                     /**
                      * \brief CTOR
                      *
                      * \param[in] id
                      * the testcase's id string
                      */
                     TaggedTest(const std::string & id)
                       : BaseTest(id)
                     {
                       _tag_name = TypeTraits<Tag_>::name();
                       _prec_name = TypeTraits<DataType_>::name();
                     };

                     /// DTOR
                     virtual ~TaggedTest() {}
                 };
  }
}

// let Doxygen ignore the following block
// \cond
// define __PRETTY_FUNCTION if not defined
#ifndef __PRETTY_FUNCTION__
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif
// \endcond

/// checks if a == b
#define TEST_CHECK_EQUAL(a, b) \
  do { \
    try { \
      BaseTest::TwoVarHolder test_h(a, b); \
      check(__PRETTY_FUNCTION__, __FILE__, __LINE__, test_h.result, \
          this->_id + "\n" +  "Expected '" #a "' to equal \n'" + test_h.s_b + \
          "'\nbut got\n'" + test_h.s_a + "'"); \
    } catch (const TestFailedException &) { \
      throw; \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception "+ FEAST::stringify(test_e.what()) + \
          " inside a TEST_CHECK_EQUAL block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL block"); \
    } \
  } while (false)

/// checks if a != b
#define TEST_CHECK_NOT_EQUAL(a, b) \
  do { \
    try { \
      BaseTest::TwoVarHolder2 test_h(a, b); \
      check(__PRETTY_FUNCTION__, __FILE__, __LINE__, !test_h.result, \
          this->_id + "\n" +  "Expected '" #a "' that is'" + test_h.s_a + \
          "' to equal not '" + test_h.s_b + "'"); \
    } catch (const TestFailedException &) { \
      throw; \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception "+ FEAST::stringify(test_e.what()) + \
          " inside a TEST_CHECK_NOT_EQUAL block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_NOT_EQUAL block"); \
    } \
  } while (false)

/// checks if stringify(a) == stringify(b)
#define TEST_CHECK_STRINGIFY_EQUAL(a, b) \
  do { \
    try { \
      std::string s_a(FEAST::stringify(a)); \
      std::string s_b(FEAST::stringify(b)); \
      check(__PRETTY_FUNCTION__, __FILE__, __LINE__, s_a == s_b, \
          this->_id + "\n" +  "Expected '" #a "' to equal '" + s_b + \
          "'\nbut got\n'" + s_a + "'"); \
    } catch (const TestFailedException &) { \
      throw; \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception  "+ FEAST::stringify(test_e.what()) + \
          " inside a TEST_CHECK_STRINGIFY_EQUAL block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_STRINGIFY_EQUAL block"); \
    } \
  } while (false)

/// checks if a is true
#define TEST_CHECK(a) \
  do { \
    try { \
      check(__PRETTY_FUNCTION__, __FILE__, __LINE__, a, \
          this->_id + "\n" +  "Check '" #a "' failed"); \
    } catch (const TestFailedException &) { \
      throw; \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception "+ FEAST::stringify(test_e.what()) + \
          " inside a TEST_CHECK block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK block"); \
    } \
  } while (false)

/// checks if a throws an exception of type b
#define TEST_CHECK_THROWS(a, b) \
  do { \
    try { \
      try { \
        a; \
        check(__PRETTY_FUNCTION__, __FILE__, __LINE__, false, \
            this->_id + "\n" +  "Expected exception of type '" #b "' not thrown"); \
      } catch (b &) { \
        TEST_CHECK(true); \
      } \
    } catch (const TestFailedException &) { \
      throw; \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception "+ FEAST::stringify(test_e.what()) + \
          " inside a TEST_CHECK_THROWS block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_THROWS block"); \
    } \
  } while (false)

/// checks if a - b < epsilon
#define TEST_CHECK_EQUAL_WITHIN_EPS(a, b, eps) \
  do { \
    try { \
      BaseTest::WithinEpsCalculator calc(a, b, eps); \
      check(__PRETTY_FUNCTION__, __FILE__, __LINE__, calc.result,  \
          this->_id + "\n" + "Expected '|" #a " - " #b \
          "|' < '" + FEAST::stringify(eps) + "' but was '" + calc.s_diff +"'"); \
    } catch (const TestFailedException & test_e) { \
      throw;  \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception  "+ FEAST::stringify(test_e.what()) + \
          " inside a TEST_CHECK_EQUAL_WITHIN_EPS block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL_WITHIN_EPS block"); \
    } \
  } while (false)

/// runs the given test with pre- and postprocessing
#define TEST(pre, test, post) \
  do { \
    pre; \
    try { \
      test; \
    } catch (const TestFailedException & test_e) { \
      post; \
      throw; } \
    post; \
  } while (false)




using namespace FEAST;
using namespace FEAST::TestSystem;

int main(int argc, char** argv)
{
  int result(EXIT_SUCCESS);
  init_mpi(argc, argv);
  int rank(-1);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc == 2)
  {
    std::string mpc("mpiproccount");
    std::string argv1(argv[1]);
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
        std::cout << "FAILED: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl;
      result = EXIT_FAILURE;
    }
    catch (InternalError & e)
    {
      if (rank == 0)
        std::cout << "FAILED: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl
          << stringify(e.message()) << std::endl;
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

  finalise_mpi();
  return result;
}

#endif // TEST_SYSTEM_TEST_SYSTEM_HPP
