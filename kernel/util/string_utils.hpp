#pragma once
#ifndef STRING_UTILS_HPP
#define STRING_UTILS_HPP

#include <sstream>
#include <string>

/**
* \brief collection of various string utilities
*
* \author Dirk Ribbrock
* \author Dominik Goeddeke
*/
class StringUtils
{
  public:
  /* *****************
  * member functions *
  *******************/

  /**
   * Convert item to a string.
   * \param[in] item
   * The item to stringify
   */
  template <typename T_>
  static inline std::string stringify(const T_ & item)
  {
    std::ostringstream s;
    s << item;
    return s.str();
  }

  /**
   * Convert item to a string (overload for std::string).
   * \param[in] item
   * The item to stringify
   */
  static inline std::string stringify(const std::string & item)
  {
    return item;
  }

  /**
   * Convert item to a string (overload for char).
   * \param[in] item
   * The item to stringify
   */
  static inline std::string stringify(const char & item)
  {
    return std::string(1, item);
  }

  /**
   * Convert item to a string (overload for unsigned char).
   * \param[in] item
   * The item to stringify
   */
  static inline std::string stringify(const unsigned char & item)
  {
    return std::string(1, item);
  }

  /**
   * Convert item to a string (overload for bool).
   * \param[in] item
   * The item to stringify
   */
  static inline std::string stringify(const bool & item)
  {
    return item ? "true" : "false";
  }

  /**
   * Convert item to a string (overload for char *, which isn't a
   * screwup like other pointers).
   * \param[in] item
   * The item to stringify
   */
  static inline std::string stringify(const char * const item)
  {
    return std::string(item);
  }

};
#endif //  #ifndef STRING_UTILS_HPP
