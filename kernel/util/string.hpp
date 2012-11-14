#pragma once
#ifndef KERNEL_UTIL_STRING_HPP
#define KERNEL_UTIL_STRING_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#include <string>
#include <sstream>
#include <locale>
#include <vector>

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
    : public std::string
  {
  public:
    /**
     * \brief A class providing case-insensitive String comparison.
     *
     * \details
     * This class provides an STL-conforming function object for case-insensitive comparison of Strings
     * which can be used in associative containers as e.g. \c std::map or \c std::set.
     *
     * \note
     * This class does not implement the String comparison itself but calls the String::compare_no_case()
     * member function for the dirty work instead.
     *
     * \author Peter Zajac
     */
    class NoCaseLess
      : public std::binary_function<String, String, bool>
    {
    public:
      /**
       * \brief Compares two Strings without regard to case.
       *
       * \param[in] left, right
       * The two Strings which are to be compared.
       *
       * \returns
       * \c true, if \p left is less than \p right without regard to case, otherwise \c false.
       */
      bool operator()(const String& left, const String& right) const
      {
        return left.compare_no_case(right) < 0;
      }
    }; // class NoCaseLess

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
     * \brief Returns a reference to the first character in the string.
     * \returns
     * A reference to the first character in the string.
     */
    reference front()
    {
      return at(0);
    }

    /** \copydoc front() */
    const_reference front() const
    {
      return at(0);
    }

    /**
     * \brief Returns a reference to the last character in the string.
     * \returns
     * A reference to the last character in the string.
     */
    reference back()
    {
      return at(size() - 1);
    }

    /** \copydoc back() */
    const_reference back() const
    {
      return at(size() - 1);
    }

    /**
     * \brief Inserts a character at the front of the string.
     *
     * \param[in] value
     * The character to the pushed.
     */
    void push_front(char value)
    {
      insert(0, 1, value);
    }

    /// Removes the first character from the string.
    void pop_front()
    {
      erase(0, 1);
    }

    /// Removes the last character from the string.
    void pop_back()
    {
      erase(size() - 1, 1);
    }

    /**
     * \brief Trims the front of the string.
     *
     * This method removes any leading characters contained in the character set from the string.
     *
     * \param[in] charset
     * The character set which is to be trimmed from the front of the string.
     *
     * \returns
     * The front-trimmed string.
     */
    String trim_front(const String& charset = " \a\b\f\n\r\t\v") const
    {
      // find first character not to be trimmed
      size_type pos = find_first_not_of(charset);
      if(pos == npos)
        return String();
      else
        return substr(pos);
    }

    /**
     * \brief Trims the back of the string.
     *
     * This method removes any trailing characters contained in the character set from the string.
     *
     * \param[in] charset
     * The character set which is to be trimmed from the back of the string.
     *
     * \returns
     * The back-trimmed string.
     */
    String trim_back(const String& charset = " \a\b\f\n\r\t\v") const
    {
      // find last character not to be trimmed
      size_type pos = find_last_not_of(charset);
      if(pos == npos)
        return String();
      else
        return substr(0, pos + 1);
    }

    /**
     * \brief Trims the string.
     *
     * This method removes any leading and trailing characters contained in the character set from the string.
     *
     * \note If you want to trim \e this string, use trim_me() instead.
     *
     * \param[in] charset
     * The character set which is to be trimmed from the string.
     *
     * \returns
     * The trimmed string.
     */
    String trim(const String& charset = " \a\b\f\n\r\t\v") const
    {
      // trim front and back
      return trim_front(charset).trim_back(charset);
    }

    /**
     * \brief Trims this string.
     *
     * This method removes any leading and trailing characters contained in the character set from this string.
     *
     * \note If you want to have the trimmed string without modifying \e this string, use trim() instead.
     *
     * \param[in] charset
     * The character set which is to be trimmed from the string.
     *
     * \returns \p *this
     */
    String& trim_me(const String& charset = " \a\b\f\n\r\t\v")
    {
      return (*this = trim(charset));
    }

    /**
     * \brief Forks the string by a delimiter charset.
     *
     * This function separates the string into substrings, where two substrings are separated a delimiter
     * substring consisting only of delimiter charset characters.
     *
     * This function is frequently used by parsers, which fork strings by whitespace characters.
     *
     * <b>Example:</b>\n
     * When using the default whitespace delimiter character set, the string " 5  42\t7 " will be forked
     * into 3 substrings: "5", "42" and "7".
     *
     * \param[out] words
     * A vector of Strings which receives the substrings.
     *
     * \param[in] charset
     * The character set which is to be treated as the delimiter charset.
     *
     * \returns
     * <c>words.size()</c>
     */
    size_t fork_by_charset(
      std::vector<String>& words,
      const String& charset = " \a\b\f\n\r\t\v") const
    {
      words.clear();
      if(empty() || charset.empty())
      {
        return 0;
      }

      // find first occurance of fork substring
      size_t off1(find_first_not_of(charset));
      if(off1 == npos)
      {
        // only delimiter characters; nothing to be extracted
        return 0;
      }

      // go forking
      while(off1 != npos)
      {
        // find next occurance of delimiter string
        size_t off2(find_first_of(charset, off1));

        // add next fork substring to vector
        if(off2 == npos)
        {
          // extract last substring
          words.push_back(substr(off1));
          return words.size();
        }
        else
        {
          // extract next substring
          words.push_back(substr(off1, (off2 == npos ? npos : off2 - off1)));
        }

        // find next occurance of fork substring
        off1 = find_first_not_of(charset, off2);
      }

      // okay
      return words.size();
    }

    /**
     * \brief Forks the string by a delimiter substring.
     *
     * This function separates the string into substrings, where the substrings are separated by a delimiter
     * string.
     *
     * <b><Example:</b>\n
     * When using "," as a delimiter string, the input string " ,5,,3" will be forked into 4 substrings:
     * " ", "5", "" and "3".
     *
     * \param[out] words
     * A vector of Strings which receives the substrings.
     *
     * \param[in] delimiter
     * The string that is to be treated as a delimiter.
     *
     * \returns
     * <c>words.size()</c>
     */
    size_t fork_by_string(
      std::vector<String>& words,
      const String& delimiter = " ") const
    {
      words.clear();
      if(empty() || delimiter.empty())
        return 0u;

      // find first occurance of delimiter substring
      size_t off1(find(delimiter));
      words.push_back(substr(0, off1));
      if(off1 == npos)
        return words.size();

      // go forking
      const size_t dellen(delimiter.size());
      while(off1 != npos)
      {
        // increase leading offset by delimiter length
        off1 += dellen;

        // find next substring occurance
        size_t off2 = find(delimiter, off1);
        if(off2 == npos)
        {
          // extract last substring
          words.push_back(substr(off1));
          return words.size();
        }
        else
        {
          // extract next substring
          words.push_back(substr(off1, off2 - off1));
        }

        // update offset
        off1 = off2;
      }

      // okay
      return words.size();
    }

    /**
     * \brief Converts the string to upper case.
     *
     * \returns
     * The upper-case string.
     */
    String upper() const
    {
      String str;
      str.reserve(size());
      for(const_iterator it(begin()); it != end(); ++it)
      {
#ifdef FEAST_COMPILER_MICROSOFT
        str.push_back(std::toupper(*it, std::locale::classic()));
#else
        str.push_back(std::toupper(*it));
#endif
      }
      return str;
    }

    /**
     * \brief Converts the string to lower case.
     *
     * \returns
     * The lower-case string.
     */
    String lower() const
    {
      String str;
      str.reserve(size());
      for(const_iterator it(begin()); it != end(); ++it)
      {
#ifdef FEAST_COMPILER_MICROSOFT
        str.push_back(std::tolower(*it, std::locale::classic()));
#else
        str.push_back(std::tolower(*it));
#endif
      }
      return str;
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
     * \brief Parses the string.
     *
     * \param[out] t
     * The value that is to be parsed.
     *
     * \returns
     * \c true, if the parse succeeds, otherwise \c false.
     */
    template<typename T_>
    bool parse(T_& t) const
    {
      std::istringstream iss(trim());
      iss >> t;
      return !iss.fail();
    }

    /// \cond internal
    bool parse(bool& b) const
    {
      if(trim().compare_no_case("true") == 0)
      {
        b = true;
        return true;
      }
      if(trim().compare_no_case("false") == 0)
      {
        b = false;
        return true;
      }
      return false;
    }

    bool parse(std::string& s) const
    {
      s.assign(*this);
      return true;
    }
    /// \endcond

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
  inline String stringify(const char* item)
  {
    return String(item);
  }

  inline String stringify(const std::string& item)
  {
    return item;
  }

  inline String stringify(char item)
  {
    return String(1, item);
  }

  inline String stringify(bool item)
  {
    return String(item ? "true" : "false");
  }
  /// \endcond
} // namespace FEAST

#endif // KERNEL_UTIL_STRING_HPP
