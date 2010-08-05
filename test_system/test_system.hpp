#pragma once
#ifndef TEST_SYSTEM_TEST_SYSTEM_HPP
#define TEST_SYSTEM_TEST_SYSTEM_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/stringify.hpp>

#include <string>
#include <exception>
#include <list>
#include <typeinfo>

/**
 * \file
 *
 * Implementation of Test and related classes.
 *
 **/

using namespace Feast;

namespace TestSystem
{
  // Forwared declaration
  class BaseTest;

  /**
   * Exception thrown by the check method in BaseTest
   */
  class TestFailedException :
    public std::exception
  {
    private:
      const std::string _message;

    public:
      /**
       * Constructor.
       */
      TestFailedException(const char* const function, const char* const file,
          const long line, const std::string & message) throw ()
      {
      }

      /**
       * Destructor.
       */
      virtual ~TestFailedException() throw ()
      {
      }

      /**
       * Description.
       */
      const char* what() const throw ()
      {
        return _message.c_str();
      }
  };

  /**
   * List of all instantiated tests
   */
  class TestList
  {
    private:
      static std::list<BaseTest*> _tests;

      TestList()
      {
      }

    public:
      typedef std::list<BaseTest*>::iterator Iterator;

      static TestList* instance()
      {
        static TestList result;

        return &result;
      }

      void register_test(BaseTest* const test)
      {
        _tests.push_back(test);
      }

      Iterator begin_tests() const
      {
        return _tests.begin();
      }

      Iterator end_tests() const
      {
        return _tests.end();
      }

      unsigned long size()
      {
        return _tests.size();
      }

      Iterator erase (Iterator i)
      {
        return _tests.erase(i);
      }

  };

  /**
   * \brief Baseclass for all Tests
   *
   * \author Dirk Ribbrock
   */
  class BaseTest
  {
    protected:
      const std::string _id;
      std::string _tag_name;
      std::string _prec_name;

    public:
      /**
       * Constructor.
       *
       * \param id The testcase's id string.
       */
      BaseTest(const std::string& id) :
        _id(id),
        _tag_name("NONE"),
        _prec_name("NONE")
      {
        TestList::instance()->register_test(this);
      }

      /// Destructor.
      virtual ~BaseTest() {}

      /// Returns our id string.
      virtual const std::string id() const
      {
        return _id;
      }

      /// Utility method used bei TEST_CHECK_*
      virtual void check(const char * const function, const char * const file,
          const long line, bool was_ok, const std::string & message) const
      {
        if (! was_ok)
          throw TestFailedException(function, file, line, message);
      }

      /**
       * \brief Runs the test case.
       *
       * Called by unittest framework only.
       */
      virtual void run() const = 0;

      /// Returns our target platform.
      virtual std::string get_tag_name()
      {
        return _tag_name;
      }

      /// Returns our target platform.
      virtual std::string get_prec_name()
      {
        return _prec_name;
      }

      /**
       * Utility class used by TEST_CHECK_EQUAL.
       */
      struct TwoVarHolder
      {
        bool result;
        std::string s_a;
        std::string s_b;

        template <typename T1_, typename T2_>
          TwoVarHolder(T1_ a, T2_ b) :
            result(a == b),
            s_a(),
            s_b()
        {
          if (!result)
          {
            s_a = stringify(a);
            s_b = stringify(b);
          }
        }
      };

      /**
       * Utility class used by TEST_CHECK_NOT_EQUAL.
       */
      struct TwoVarHolder2
      {
        bool result;
        std::string s_a;
        std::string s_b;

        template <typename T1_, typename T2_>
          TwoVarHolder2(T1_ a, T2_ b) :
            result(a == b),
            s_a(),
            s_b()
        {
          if (result)
          {
            s_a = stringify(a);
            s_b = stringify(b);
          }
        }
      };

      /**
       * Utility class used by TEST_CHECK_EQUAL_WITHIN_EPS.
       */
      struct WithinEpsCalculator
      {
        bool result;
        std::string s_a;
        std::string s_b;
        std::string s_diff;

        template <typename T1_, typename T2_, typename T3_>
          WithinEpsCalculator(T1_ a, T2_ b, T3_ c) :
            s_a(),
            s_b()
        {
          if (a >= b)
          {
            result = ((a - b) <= c);
            if (!result)
            {
              s_diff = stringify(a - b);
              s_a = stringify(a);
              s_b = stringify(b);
            }
          }
          else
          {
            result = ((b - a) <= c);
            if (!result)
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
   * \brief Abstract Baseclass for all tagged test classes.
   *
   * \author Dirk Ribbrock
   */
  template <typename Tag_, typename DataType_>
    class TaggedTest : public BaseTest
  {
    public:
      /**
       * Constructor.
       *
       * \param id The testcase's id string.
       */
      TaggedTest(const std::string & id):
        BaseTest(id)
    {
      /// \todo Use tag provided name.
      //_tag_name = Tag_::name;
      _tag_name = stringify(typeid(Tag_).name());
      _prec_name = stringify(typeid(DataType_).name());
    };

      /// Destructor.
      virtual ~TaggedTest() {}
  };
}

/**
 * Check that a == b.
 */
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
          "Test threw unexpected exception "+ Feast::stringify(test_e.what()) + \
          " inside a TEST_CHECK_EQUAL block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL block"); \
    } \
  } while (false)

/**
 * Check that a != b.
 */
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
          "Test threw unexpected exception "+ Feast::stringify(test_e.what()) + \
          " inside a TEST_CHECK_NOT_EQUAL block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_NOT_EQUAL block"); \
    } \
  } while (false)

/**
 * Check that stringify(a) == stringify(b).
 */
#define TEST_CHECK_STRINGIFY_EQUAL(a, b) \
  do { \
    try { \
      std::string s_a(Feast::stringify(a)); \
      std::string s_b(Feast::stringify(b)); \
      check(__PRETTY_FUNCTION__, __FILE__, __LINE__, s_a == s_b, \
          this->_id + "\n" +  "Expected '" #a "' to equal '" + s_b + \
          "'\nbut got\n'" + s_a + "'"); \
    } catch (const TestFailedException &) { \
      throw; \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception  "+ Feast::stringify(test_e.what()) + \
          " inside a TEST_CHECK_STRINGIFY_EQUAL block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_STRINGIFY_EQUAL block"); \
    } \
  } while (false)

/**
 * Check that a is true.
 */
#define TEST_CHECK(a) \
  do { \
    try { \
      check(__PRETTY_FUNCTION__, __FILE__, __LINE__, a, \
          this->_id + "\n" +  "Check '" #a "' failed"); \
    } catch (const TestFailedException &) { \
      throw; \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception "+ Feast::stringify(test_e.what()) + \
          " inside a TEST_CHECK block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK block"); \
    } \
  } while (false)

/**
 * Check that a throws an exception of type b.
 */
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
          "Test threw unexpected exception "+ Feast::stringify(test_e.what()) + \
          " inside a TEST_CHECK_THROWS block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_THROWS block"); \
    } \
  } while (false)

/**
 * Check that a - b < epsilon.
 */
#define TEST_CHECK_EQUAL_WITHIN_EPS(a, b, eps) \
  do { \
    try { \
      BaseTest::WithinEpsCalculator calc(a, b, eps); \
      check(__PRETTY_FUNCTION__, __FILE__, __LINE__, calc.result,  \
          this->_id + "\n" + "Expected '|" #a " - " #b \
          "|' < '" + stringify(eps) + "' but was '" + calc.s_diff +"'"); \
    } catch (const TestFailedException & test_e) { \
      throw;  \
    } catch (const std::exception & test_e) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected exception  "+ Feast::stringify(test_e.what()) + \
          " inside a TEST_CHECK_EQUAL_WITHIN_EPS block"); \
    } catch (...) { \
      throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL_WITHIN_EPS block"); \
    } \
  } while (false)

/**
 * Run the given test with pre- and postprocessing
 */
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


#endif // TEST_SYSTEM_TEST_SYSTEM_HPP
