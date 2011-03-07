#pragma once
#ifndef STRING_UTILS_HPP
#define STRING_UTILS_HPP

// includes, system
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <cassert>

/**
* \file collection of various string utilities
*
* \author Dirk Ribbrock
* \author Dominik Goeddeke
*/

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief converting an item to a string
  *
  * \tparam T_
  * type to be stringified
  *
  * \param[in] item
  * the item to stringify
  */
  template<typename T_>
  inline std::string stringify(const T_ & item)
  {
    std::ostringstream s;
    s << item;
    return s.str();
  }


  /**
  * \brief converting an item to a string (overload for std::string)
  *
  * \param[in] item
  * the item to stringify
  */
  inline std::string stringify(const std::string & item)
  {
    return item;
  }


  /**
  * \brief converting an item to a string (overload for char)
  *
  * \param[in] item
  * the item to stringify
  */
  inline std::string stringify(const char & item)
  {
    return std::string(1, item);
  }


  /**
  * \brief converting an item to a string (overload for unsigned char)
  *
  * \param[in] item
  * the item to stringify
  */
  inline std::string stringify(const unsigned char & item)
  {
    return std::string(1, item);
  }


  /**
  * \brief converting an item to a string (overload for bool)
  *
  * \param[in] item
  * the item to stringify
  */
  inline std::string stringify(const bool & item)
  {
    return item ? "true" : "false";
  }


  /**
  * \brief converts an item to a string (overload for char *, which isn't a screwup like other pointers)
  *
  * This function is selected when using the syntax \code stringify("some string") \endcode .
  *
  * \param[in] item
  * the item to stringify
  */
  inline std::string stringify(const char * const item)
  {
    return std::string(item);
  }


  /**
  * \brief concatenates strings between begin and end iterator
  *
  * \param[in] begin
  * first element
  *
  * \param[in] end
  * last element
  *
  * \param[in] delimiter
  * string to use as delimiter
  *
  * \todo also implement for std::vector iterator
//COMMENT_HILMAR: is there a clever way to provide this functionality for all STL containers at once?
  *
  * \author Dirk Ribbrock
  */
  std::string join_strings(
    std::list<std::string>::const_iterator begin,
    std::list<std::string>::const_iterator end,
    const std::string& delimiter)
  {
    std::string result;

    if (begin != end)
    {
      while (true)
      {
        result += *begin;
        if (++begin == end)
        {
          break;
        }
        result += delimiter;
      }
    }
// COMMENT_HILMAR: Do not add linebreak. The user can do this by himself if he wants.
//      result +="\n";
    return result;
  }
} // namespace FEAST

#endif //  #ifndef STRING_UTILS_HPP
