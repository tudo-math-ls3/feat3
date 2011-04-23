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


  /**
  * \brief trims leading and trailing white spaces from a string
  *
  * This function simply removes all sorts of leading and trailing whitespace from the given string and returns the
  * result in a new string object.
  *
  * \note basically taken from Bruce Eckel, Thinking in C++, Volume 2
  *
  * \param[in] str
  * string to be trimmed
  *
  * \return trimmed string
  *
  * \author Hilmar Wobker
  */
  std::string trim(const std::string& str)
  {
    if(str.length() == 0)
    {
      return str;
    }
    size_t start = str.find_first_not_of(" \a\b\f\n\r\t\v");
    size_t end = str.find_last_not_of(" \a\b\f\n\r\t\v");
    // test if there are any non-whitespace characters at all
    if(start == std::string::npos)
    {
      return "";
    }
    return std::string(str, start, end - start + 1);
  }


  /**
  * \brief turns a string into its upper case variant
  *
  * \note basically taken from Bruce Eckel, Thinking in C++, Volume 2
  *
  * \param[in] str
  * string to be turned to upper case
  *
  * \return upper case string
  *
  * \author Hilmar Wobker
  */
  inline std::string upper_case(const std::string& str)
  {
    std::string str_up(str);
    for(size_t i = 0; i < str.length(); ++i)
    {
      str_up[i] = std::toupper(str_up[i]);
    }
    return str_up;
  }


  /**
  * \brief turns a string into its lower case variant
  *
  * \note basically taken from Bruce Eckel, Thinking in C++, Volume 2
  *
  * \param[in] str
  * string to be turned to lower case
  *
  * \return lower case string
  *
  * \author Hilmar Wobker
  */
  inline std::string lower_case(const std::string& str)
  {
    std::string str_low(str);
    for(size_t i = 0; i < str.length(); ++i)
    {
      str_low[i] = std::tolower(str_low[i]);
    }
    return str_low;
  }

} // namespace FEAST

#endif //  #ifndef STRING_UTILS_HPP
