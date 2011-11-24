#pragma once
#ifndef KERNEL_UTIL_STRING_HPP
#define KERNEL_UTIL_STRING_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#include <string>
#include <sstream>
#include <locale>

#ifdef FEAST_COMPILER_MICROSOFT
#  include <string.h> // for _stricmp
#endif

namespace FEAST
{
  /**
   * \brief String class implementation.
   *
   * This class inherits from the STL class \c std::string, so in consequence any algorithm working
   * with \c std::string will also work with this class.
   *
   * \author Peter Zajac
   */
  class String
    : virtual public std::string
  {
  public:
    /// default constructor
    inline String()
      : std::string()
    {
    }

    /// CTOR
    inline String(const char* str)
      : std::string(str)
    {
    }

    /// CTOR
    inline String(const char* str, size_type count)
      : std::string(str, count)
    {
    }

    /// CTOR
    inline String(const std::string& str)
      : std::string(str)
    {
    }

    /// CTOR
    inline String(const std::string& str, size_type offset, size_type count = npos)
      : std::string(str, offset, count)
    {
    }

    /// CTOR
    inline String(size_type count, char c)
      : std::string(count, c)
    {
    }

    /**
     * \brief Trims the string.
     *
     * This method removes any leading and trailing characters contained in the character set from the string.
     *
     * \param charset
     * The character set which is to be trimmed from the string.
     *
     * \returns \p *this
     */
    String& trim(const String& charset = " \a\b\f\n\r\t\v")
    {
      if(empty() || charset.empty())
        return *this;

      // remove trailing chars
      size_type pos = find_last_not_of(charset);
      if(pos == npos)
      {
        // the whole string consists only of characters to be trimmed, so clear the whole string
        clear();
      }
      else
      {
        // remove trailing chars
        erase(pos+1, npos);

        // remove leading chars
        pos = find_first_not_of(charset);
        if(pos > 0)
        {
          erase(0, pos);
        }
      }

      return *this;
    }

    /**
     * \brief Converts the string to upper case.
     *
     * \returns \p *this
     */
    String& make_upper_case()
    {
      for(size_type i = 0; i < size(); ++i)
      {
#ifdef FEAST_COMPILER_MICROSOFT
        (*this)[i] = std::toupper((*this)[i], std::locale::classic());
#else
        (*this)[i] = std::toupper((*this)[i]);
#endif
      }
      return *this;
    }

    /**
     * \brief Converts the string to lower case.
     *
     * \returns \p *this
     */
    String& make_lower_case()
    {
      for(size_type i = 0; i < size(); ++i)
      {
#ifdef FEAST_COMPILER_MICROSOFT
        (*this)[i] = std::tolower((*this)[i], std::locale::classic());
#else
        (*this)[i] = std::tolower((*this)[i]);
#endif
      }
      return *this;
    }

    /**
     * \brief Compares two strings without regard to case.
     *
     * \param[in] other
     * The string that is to be compared to \p this.
     *
     * \returns
     *  - 0 if both strings are equal without regard to case.
     *  - -1, if \p this is less than \p other
     *  - +1, if \p this is greater than \p other
     */
    int compare_no_case(const String& other) const
    {
#ifdef FEAST_COMPILER_MICROSOFT
      // The MS C library offers a function for this task.
      return _stricmp(c_str(), other.c_str());
#else
      // reference implementation
      int k;
      size_type n1 = size();
      size_type n2 = other.size();
      size_type n = std::min(n1, n2);

      // loop over all characters and compare them
      for(size_t i = 0; i < n; ++i)
      {
        k = int(std::tolower((*this)[i])) - int(std::tolower(other[i]));
        if(k < 0)
        {
          return -1;
        }
        else if(k > 0)
        {
          return 1;
        }
      }

      // If we come out here, then both strings are identical for the first n characters.
      // Now let's check whether their length is equal, too.
      if(n1 < n2)
      {
        return -1;
      }
      else if(n1 > n2)
      {
        return 1;
      }
      else
      {
        return 0;
      }
#endif
    }

    /**
     * \brief Appends a set of strings.
     *
     * This method appends a set of strings given by two input iterators, separating them by a delimiter.
     *
     * \param[in] first, last
     * Two forward iterators representing the string set.
     *
     * \param[in] delimiter
     * A string acting as a delimiter.
     *
     * \returns \p *this
     */
    template<typename Iterator_>
    String& join(
      Iterator_ first,
      Iterator_ last,
      const String& delimiter = "")
    {
      Iterator_ it(first);
      while(it != last)
      {
        append(*it);
        if(++it == last)
        {
          return *this;
        }
        append(delimiter);
      }
      return *this;
    }

    /**
     * \brief Appends a set of strings.
     *
     * This method appends a set of strings given by a container, e.g. <c>std::list<String></c>, separating each
     * entry in the container by a delimiter.
     *
     * \param[in] container
     * A string container.
     *
     * \param[in] delimiter
     * A string acting as a delimiter.
     *
     * \returns \p *this
     */
    template<typename Container_>
    String& join(
      const Container_& container,
      const String& delimiter = "")
    {
      return join(container.cbegin(), container.cend(), delimiter);
    }
  }; // class String

  /**
   * \brief Converts an item into a String.
   *
   * \tparam T_
   * The type of the item to be converted.
   *
   * \param[in] item
   * The item to be converted.
   *
   * \returns
   * A String representation of the item.
   */
  template<typename T_>
  inline String stringify(const T_& item)
  {
    std::ostringstream oss;
    oss << item;
    return oss.str();
  }

  /// \cond internal
  template<>
  inline String stringify<const char*>(const char* const& item)
  {
    return String(item);
  }

  template<>
  inline String stringify<std::string>(const std::string& item)
  {
    return item;
  }

  template<>
  inline String stringify<char>(const char& item)
  {
    return String(1, item);
  }

  template<>
  inline String stringify<bool>(const bool& item)
  {
    return String(item ? "true" : "false");
  }
  /// \endcond
} // namespace FEAST

#endif // KERNEL_UTIL_STRING_HPP
