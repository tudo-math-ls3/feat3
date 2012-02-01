#pragma once

// includes, FEAST
#include <kernel/util/string.hpp>

// includes, system
#include <map>

namespace FEAST
{
  /**
   * \brief String map class
   *
   * This class implements a map of Strings with case insensitive comparison.
   * This class can be used to store strings for variable substitution.
   *
   * \author Peter Zajac
   */
  class StringMap
    : public std::map<String, String, String::NoCaseLess>
  {
  public:
    /// base class typedef
    typedef std::map<String, String, String::NoCaseLess> BaseClass;

  public:
    /// default constructor
    StringMap()
      : BaseClass()
    {
    }

    /// virtual destructor
    virtual ~StringMap()
    {
    }
  }; // class StringMap
} // namespace FEAST