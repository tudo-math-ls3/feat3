// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef TEST_SYSTEM_TEST_SYSTEM_HPP
/// Header guard
#define TEST_SYSTEM_TEST_SYSTEM_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/util/exception.hpp>

// includes, system
#include <string>
#include <exception>
#include <list>
#include <typeinfo>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cmath>

#define CHECK_INTERNAL(was_ok, message)\
        if(! (was_ok))\
          throw FEAT::TestSystem::TestFailedException(__func__, __FILE__, __LINE__, message);

/**
* \file
*
* Implementation of Test and related classes.
*/

namespace FEAT
{
  /// TestSystem namespace
  namespace TestSystem
  {
    // Forwared declaration
    class UnitTest;

    /// exception thrown by the check method in UnitTest
    class TestFailedException
      : public std::exception
    {

    private:

      /// failure message
      const String _message;


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
          const String & message)
        : _message(stringify(file) + ":" + stringify(line) + ": in " + stringify(function) + ": " + message )
      {
      }

      /// DTOR
      virtual ~TestFailedException() throw ()
      {
      }

      /// description
      virtual const char* what() const throw() override
      {
        return _message.c_str();
      }
    }; // class TestFailedException

    /// list of all instantiated tests
    class TestList
    {
    private:
      class DeleteOnDestruction
      {
      private:
        /// Pointer to our managed object
        TestList * * const _ptr;

      public:
        /// Constructor
        explicit DeleteOnDestruction(TestList * * const ptr)
          : _ptr(ptr)
        {
        }

        /// Destructor
        ~DeleteOnDestruction()
        {
          TestList::_delete(*_ptr);
          *_ptr = nullptr;
        }
      };

      /// Returns a pointer to our instance pointer.
      static TestList * * _instance_ptr()
      {
        // The Microsoft compiler warns that creating a local static object is not thread-safe.
        // The corresponding warning is enabled by default, so we'll only disable it for the following
        // code and restore the warning state afterwards.
#ifdef FEAT_COMPILER_MICROSOFT
#  pragma warning(push)
#  pragma warning(disable: 4640)
#endif

        static TestList * instance(nullptr);
        static DeleteOnDestruction delete_instance(&instance);

#ifdef FEAT_COMPILER_MICROSOFT
#  pragma warning(pop)
#endif

        return &instance;
      }

      static void _delete(TestList * const ptr)
      {
        delete ptr;
      }

      /// internal STL list representation of TestList
      std::list<UnitTest*> _tests;

      /// default CTOR
      TestList()
      {
      }

      TestList(const TestList&) = delete;
      TestList& operator=(const TestList&) = delete;

    public:
      /**
       * \brief Returns the instance.
       *
       * \returns
       * A pointer to the singleton.
       */
      static TestList * instance()
      {
        TestList * * instance_ptr = _instance_ptr();

        if (nullptr == *instance_ptr)
        {
          /// \todo Make thread safe
          //static Mutex m;
          //Lock l(m);

          instance_ptr = _instance_ptr();

          if (nullptr == *instance_ptr)
          {
            *instance_ptr = new TestList;
          }
        }
        return *instance_ptr;
      }

      /// TestList forward iterator.
      typedef std::list<UnitTest*>::iterator Iterator;

      /**
        * \brief adds a test to the TestList
        *
        * \param[in] test
        * the test that will be added
        */
      void register_test(UnitTest* const test)
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
    }; // class TestList


    /**
     * \brief base class for all Tests
     *
     * \author Dirk Ribbrock
     */
    class UnitTest
    {
    protected:

      /// test description String
      const String _id;
      /// precision description String
      String _datatype_name;
      /// index type description String
      String _index_name;
      /// preferred compute intensive backend
      PreferredBackend _preferred_backend;


    public:

      /**
        * \brief CTOR
        *
        * \param[in] id_in
        * the testcase's id string
        */
      explicit UnitTest(const String& id_in, const String datatype_name = "none", const String index_name = "none", PreferredBackend preferred_backend = PreferredBackend::generic)
        : _id(id_in),
        _datatype_name(datatype_name),
        _index_name(index_name),
        _preferred_backend(preferred_backend)
      {
        TestList::instance()->register_test(this);
      }

      /// DTOR
      virtual ~UnitTest() {}

      /// returns our id string
      virtual const String id() const
      {
        return _id;
      }

      /// utility method used bei TEST_CHECK_THROWS
      virtual void check(const char * const function, const char * const file,
          const long line, bool was_ok, const String & message) const
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
      virtual String get_datatype_name()
      {
        return _datatype_name;
      }

      /// returns our target platform
      virtual String get_index_name()
      {
        return _index_name;
      }

      virtual PreferredBackend get_preferred_backend() const
      {
        return _preferred_backend;
      }

      virtual String get_preferred_backend_name() const
      {
        return stringify(_preferred_backend);
      }
    }; // class UnitTest

  } // namespace TestSystem
} // namespace FEAT
/// checks if a == b
#define TEST_CHECK_EQUAL(a, b) \
  do { \
    CHECK_INTERNAL((a)==(b), this->_id + "\n" +  "Expected '" #a "' to equal \n'" + FEAT::stringify(b) + "'\nbut got\n'" + FEAT::stringify(a) + "'")\
  } while (false)

/// checks if a != b
#define TEST_CHECK_NOT_EQUAL(a, b) \
  do { \
    CHECK_INTERNAL(!((a)==(b)), this->_id + "\n" +  "Expected '" #a "' that is'" + FEAT::stringify(a) + "' to equal not '" + FEAT::stringify(b) + "'")\
  } while (false)

/// checks if a <= x <= b
#define TEST_CHECK_IN_RANGE(x, a, b) \
  do { \
    CHECK_INTERNAL(((a) <= (x)) && ((x) <= (b)), this->_id + "\n" +  "Expected '" #x "' that is'" + FEAT::stringify(x) + "' to be in range [" + FEAT::stringify(a) + "," + FEAT::stringify(b) + "]")\
  } while (false)

/// checks if stringify(a) == stringify(b)
#define TEST_CHECK_STRINGIFY_EQUAL(a, b) \
  do { \
    String s_a(FEAT::stringify(a)); \
    String s_b(FEAT::stringify(b)); \
    CHECK_INTERNAL(s_a == s_b, this->_id + "\n" +  "Expected '" #a "' to equal '" + s_b + "'\nbut got\n'" + s_a + "'")\
  } while (false)

/// checks if a is true
#define TEST_CHECK(a) \
  do { \
    CHECK_INTERNAL(a, this->_id + "\n" +  "Check '" #a "' failed") \
  } while (false)

/// checks if \c a is true; prints \c msg if \c a is false
#define TEST_CHECK_MSG(a,msg) \
  do { \
    CHECK_INTERNAL(a, this->_id + "\n" + (msg))\
  } while (false)

/// checks if a throws an exception of type b
#define TEST_CHECK_THROWS(a, b) \
  do { \
    try { \
      try { \
        a; \
        this->check(__func__, __FILE__, __LINE__, false, \
            this->_id + "\n" +  "Expected exception of type '" #b "' not thrown"); \
      } catch (b &) { \
        TEST_CHECK(true); \
      } \
    } catch (const FEAT::TestSystem::TestFailedException &) { \
      throw; \
    } catch (const std::exception & test_e) { \
      throw FEAT::TestSystem::TestFailedException(__func__, __FILE__, __LINE__, \
          "Test threw unexpected exception "+ FEAT::stringify(test_e.what()) + \
          " inside a TEST_CHECK_THROWS block"); \
    } catch (...) { \
      throw FEAT::TestSystem::TestFailedException(__func__, __FILE__, __LINE__, \
          "Test threw unexpected unknown exception inside a TEST_CHECK_THROWS block"); \
    } \
  } while (false)

/// checks if |a - b| <= epsilon
#define TEST_CHECK_EQUAL_WITHIN_EPS(a, b, eps) \
  do { \
    CHECK_INTERNAL(((a) <= ((b) + (eps))) && ((b) <= ((a) + (eps))), \
        this->_id + "\n" + "Expected '|" #a " - " #b \
        "|' <= '" + FEAT::stringify(eps) + "' but was '" + \
        FEAT::stringify((a) < (b) ? (b) - (a) : (a) - (b)) + "'" \
        + ", with " #a "=" + FEAT::stringify(a) + " and " #b "=" + FEAT::stringify(b))\
  } while (false)

/// runs the given test with pre- and postprocessing
#define TEST(pre, test, post) \
  do { \
    pre; \
    try { \
      test; \
    } catch (const FEAT::TestSystem::TestFailedException & test_e) { \
      post; \
      throw; } \
    post; \
  } while (false)

#endif // TEST_SYSTEM_TEST_SYSTEM_HPP
