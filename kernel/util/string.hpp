// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_STRING_HPP
#define KERNEL_UTIL_STRING_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>

// includes, system
#include <locale>
#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <iomanip>
#include <cstddef>

#ifdef FEAT_COMPILER_MICROSOFT
#  include <string.h> // for _stricmp
#endif

#ifndef __CUDACC__
#ifdef FEAT_HAVE_QUADMATH
extern "C"
{
#    include <quadmath.h>
}
#endif // FEAT_HAVE_QUADMATH
#endif // __CUDACC__

namespace FEAT
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

#ifndef __CUDACC__
    /// CTOR
    inline String(std::string&& str)
      : std::string(str)
    {
    }
#endif

    /// copy CTOR
    inline String(const String& str)
      : std::string(str)
    {
    }

#ifndef __CUDACC__
    /// move CTOR
    inline String(String&& str)
      : std::string(str)
    {
    }
#endif

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

    /// \cond internal
    String& operator=(const std::string& str)
    {
      std::string::operator=(str);
      return *this;
    }

#ifndef __CUDACC__
    String& operator=(std::string&& str)
    {
      std::string::operator=(str);
      return *this;
    }
#endif

    String& operator=(const String& str)
    {
      std::string::operator=(str);
      return *this;
    }

#ifndef __CUDACC__
    String& operator=(String&& str)
    {
      std::string::operator=(str);
      return *this;
    }
#endif

    String& operator=(const char* s)
    {
      std::string::operator=(s);
      return *this;
    }

    String& append(const std::string& str)
    {
      std::string::append(str);
      return *this;
    }

    String& append(const std::string& str, size_type pos, size_type n)
    {
      std::string::append(str, pos, n);
      return *this;
    }

    String& append(const char* s, size_type n)
    {
      std::string::append(s, n);
      return *this;
    }

    String& append(const char* s)
    {
      std::string::append(s);
      return *this;
    }

    String& append(size_type n, char c)
    {
      std::string::append(n, c);
      return *this;
    }

    String& operator+=(const std::string& str)
    {
      return append(str);
    }

    String& operator+=(const char* s)
    {
      return append(s);
    }

    String substr(size_type pos = size_type(0), size_type n = npos) const
    {
      return String(std::string::substr(pos, n));
    }
    /// \endcond

    /**
     * \brief Returns a null-terminated char string containing all white-space characters
     */
    static const char* whitespaces()
    {
      return " \a\b\f\n\r\t\v";
    }

    /**
     * \brief Inserts a character at the front of the string.
     *
     * \param[in] value
     * The character to the pushed.
     */
    void push_front(char value)
    {
      insert(size_type(0), size_type(1), value);
    }

    /// Removes the first character from the string.
    void pop_front()
    {
      erase(size_type(0), size_type(1));
    }

    /// Removes the last character from the string.
    void pop_back()
    {
      erase(size() - size_type(1), size_type(1));
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
    String trim_front(const String & charset) const
    {
      // find first character not to be trimmed
      size_type pos = find_first_not_of(charset);
      if(pos == npos)
        return String();
      else
        return substr(pos);
    }

    /**
     * \brief Trims the front of the string of all white-spaces.
     *
     * \returns
     * The front-trimmed string.
     */
    String trim_front() const
    {
      return trim_front(whitespaces());
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
    String trim_back(const String & charset) const
    {
      // find last character not to be trimmed
      size_type pos = find_last_not_of(charset);
      if(pos == npos)
        return String();
      else
        return substr(size_type(0), pos + size_type(1));
    }

    /**
     * \brief Trims the back of the string of all white-spaces.
     *
     * \returns
     * The back-trimmed string.
     */
    String trim_back() const
    {
      return trim_back(whitespaces());
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
    String trim(const String & charset) const
    {
      // trim front and back
      return trim_front(charset).trim_back(charset);
    }

    /**
     * \brief Trims the string of all white-spaces.
     *
     * \note If you want to trim \e this string, use trim_me() instead.
     *
     * \returns
     * The trimmed string.
     */
    String trim() const
    {
      return trim(whitespaces());
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
    String& trim_me(const String & charset)
    {
      return (*this = trim(charset));
    }

    /**
     * \brief Trims this string of all white-spaces.
     *
     * \note If you want to have the trimmed string without modifying \e this string, use trim() instead.
     *
     * \returns \p *this
     */
    String& trim_me()
    {
      return (*this = trim());
    }

    /**
     * \brief Pads the front of the string up to a desired length.
     *
     * This function returns a string that is front-padded with a specific character up to a desired minimum length.
     * If the length of \c this already has the desired minimum lengh, this function returns \c *this.
     *
     * \note This function is virtually the counter-part of #trim_front() and #trunc_front().
     *
     * \param[in] len
     * The desired (minimum) length of the string.
     *
     * \param[in] c
     * The character to be used for padding.
     *
     * \returns
     * The padded string.
     */
    String pad_front(size_type len, char c = ' ') const
    {
      size_type l(length());
      return (l < len) ? String(len - l, c).append(*this) : *this;
    }

    /**
     * \brief Pads the back of the string up to a desired length.
     *
     * This function returns a string that is back-padded with a specific character up to a desired minimum length.
     * If the length of \c this already has the desired minimum lengh, this function returns \c *this.
     *
     * \note This function is virtually the counter-part of #trim_back() and #trunc_back().
     *
     * \param[in] len
     * The desired (minimum) length of the string.
     *
     * \param[in] c
     * The character to be used for padding.
     *
     * \returns
     * The padded string.
     */
    String pad_back(size_type len, char c = ' ') const
    {
      size_type l(length());
      return (l < len) ? String(*this).append(len - l, c) : *this;
    }

    /**
     * \brief Truncates the front of the string to a given maximum length.
     *
     * \param[in] len
     * The desired maximum length of the string.
     *
     * \returns
     * The truncated string.
     */
    String trunc_front(size_type len) const
    {
      return length() <= len ? *this : substr(length() - len);
    }

    /**
     * \brief Truncates the back of the string to a given maximum length.
     *
     * \param[in] len
     * The desired maximum length of the string.
     *
     * \returns
     * The truncated string.
     */
    String trunc_back(size_type len) const
    {
      return length() <= len ? *this : substr(size_type(0), len);
    }

    /**
     * \brief Splits the string by a delimiter charset.
     *
     * This function separates the string into substrings, where two substrings are separated a delimiter
     * substring consisting only of delimiter charset characters.
     *
     * This function is frequently used by parsers, which split strings by whitespace characters.
     *
     * <b>Example:</b>\n
     * When using the default whitespace delimiter character set, the string " 5  42\t7 " will be split
     * into 3 substrings: "5", "42" and "7".
     *
     * \param[in] charset
     * The character set which is to be treated as the delimiter charset.
     *
     * \returns
     * A deque of the sub-strings that resulted from the splitting.
     */
    std::deque<String> split_by_charset(const String & charset) const
    {
      std::deque<String> words;
      if(empty() || charset.empty())
      {
        return words;
      }

      // find first occurrence of split substring
      size_type off1(find_first_not_of(charset));
      if(off1 == npos)
      {
        // only delimiter characters; nothing to be extracted
        return words;
      }

      // go splitting
      while(off1 != npos)
      {
        // find next occurrence of delimiter string
        size_type off2(find_first_of(charset, off1));

        // add next split substring to vector
        if(off2 == npos)
        {
          // extract last substring
          words.push_back(substr(off1));
          return words;
        }
        else
        {
          // extract next substring
          words.push_back(substr(off1, (off2 == npos ? npos : off2 - off1)));
        }

        // find next occurrence of split substring
        off1 = find_first_not_of(charset, off2);
      }

      // okay
      return words;
    }

    /**
     * \brief Splits the string by white-spaces.
     *
     * \see #split_by_charset()
     *
     * \returns
     * A deque of the sub-strings that resulted from the splitting.
     */
    std::deque<String> split_by_whitespaces() const
    {
      return split_by_charset(whitespaces());
    }

    /**
     * \brief Splits the string by a delimiter substring.
     *
     * This function separates the string into substrings, where the substrings are separated by a delimiter
     * string.
     *
     * <b>Example:</b>\n
     * When using "," as a delimiter string, the input string " ,5,,3" will be split into 4 substrings:
     * " ", "5", "" and "3".
     *
     * \param[in] delimiter
     * The string that is to be treated as a delimiter.
     *
     * \returns
     * A deque of the sub-strings that resulted from the splitting.
     */
    std::deque<String> split_by_string(const String & delimiter) const
    {
      std::deque<String> words;
      if(empty() || delimiter.empty())
        return words;

      // find first occurrence of delimiter substring
      size_type off1(find(delimiter));
      words.push_back(substr(0, off1));
      if(off1 == npos)
        return words;

      // go splitting
      const size_type del_len(delimiter.size());
      while(off1 != npos)
      {
        // increase leading offset by delimiter length
        off1 += del_len;

        // find next substring occurrence
        size_type off2 = find(delimiter, off1);
        if(off2 == npos)
        {
          // extract last substring
          words.push_back(substr(off1));
          return words;
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
      return words;
    }

    /**
     * \brief Replaces all occurrences of a substring by another substring.
     *
     * This function will replace all occurrences of a substring \c find_string by \c replace_string.
     *
     * \note This function is not successive, i.e. it will \b not replace an occurrence of the find-string
     * if this occurrence is a result of a previous replacement, e.g. for the combination \c *this = "aaa",
     * \c find_string = "aa" and \c replace_string = "pa" the resulting string will be "paa" and not "ppa".
     *
     * \param[in] find_string
     * The substring that is to be searched for. Must not be empty.
     *
     * \param[in] replace_string
     * The substring that is to be replaced by.
     *
     * \returns
     * The total number of replacements made.
     */
    size_type replace_all(const String & find_string, const String & replace_string)
    {
      size_type flen(find_string.size());
      size_type rlen(replace_string.size());
      if(flen <= 0)
        return size_type(0);

      // find first occurrence of find string
      size_type pos(find(find_string));
      size_type counter(size_type(0));
      while(pos != npos)
      {
        // replace substring
        replace(pos, flen, replace_string);

        // increment counter
        ++counter;

        // find next occurrence
        pos = find(find_string, pos + rlen);
      }

      // return replacement count
      return counter;
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
#ifdef FEAT_COMPILER_MICROSOFT
        str.push_back(std::toupper(*it, std::locale::classic()));
#else
        std::locale loc;
        str.push_back(std::toupper<char>(*it, loc));
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
#ifdef FEAT_COMPILER_MICROSOFT
        str.push_back(std::tolower(*it, std::locale::classic()));
#else
        std::locale loc;
        str.push_back(std::tolower<char>(*it, loc));
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
#ifdef FEAT_COMPILER_MICROSOFT
      // The MS C library offers a function for this task.
      return _stricmp(c_str(), other.c_str());
#else
      // reference implementation
      std::locale loc;
      size_type n1 = size();
      size_type n2 = other.size();
      size_type n = std::min(n1, n2);

      // loop over all characters and compare them
      for(size_type i = 0; i < n; ++i)
      {
        int k = int(std::tolower<char>((*this)[i], loc)) - int(std::tolower<char>(other[i], loc));
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
     * \brief Checks whether this string starts with another string.
     *
     * \param[in] head
     * The string that the front of this string is to be checked against.
     *
     * \returns
     * \c true, if \c this starts with \p head, otherwise \c false
     */
    bool starts_with(const String& head) const
    {
      // every string starts with an empty string
      if(head.empty())
        return true;

      // check size
      if(this->size() <  head.size())
        return false;

      // compare head
      return (this->compare(std::size_t(0), head.size(), head) == 0);
    }

    /**
     * \brief Checks whether this string ends with another string.
     *
     * \param[in] tail
     * The string that the back of this string is to be checked against.
     *
     * \returns
     * \c true, if \c this ends with \p tail, otherwise \c false
     */
    bool ends_with(const String& tail) const
    {
      // every string ends with an empty string
      if(tail.empty())
        return true;

      // check size
      if(this->size() <  tail.size())
        return false;

      // compare tail
      return (this->compare(this->size() - tail.size(), tail.size(), tail) == 0);
    }

    /**
     * \brief Checks whether this string starts with a specified character.
     *
     * \param[in] head
     * The character that the front of this string is to be checked against.
     *
     * \returns
     * \c true, if \c this starts with \p head, otherwise \c false
     */
    bool starts_with(const char head) const
    {
      return (this->empty() ? false : this->front() == head);
    }

    /**
     * \brief Checks whether this string ends with a specified character.
     *
     * \param[in] tail
     * The character that the back of this string is to be checked against.
     *
     * \returns
     * \c true, if \c this ends with \p tail, otherwise \c false
     */
    bool ends_with(const char tail) const
    {
      return (this->empty() ? false : this->back() == tail);
    }

    /**
     * \brief Parses the string and stores its value in the provided variable.
     *
     * \param[out] t
     * The parsed value, if the parse succeeded.
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

    // This one is really required, as otherwise some compilers choose the
    // generic template instead of the overload for std::string above...
    bool parse(String& s) const
    {
      s.assign(*this);
      return true;
    }

#ifndef __CUDACC__
#ifdef FEAT_HAVE_QUADMATH
    bool parse(__float128& x) const
    {
      if(this->empty())
        return false;

      const char* nptr(this->c_str());
      char* endptr(nullptr);
      x = strtoflt128(nptr, &endptr);
      // Note: According to the C-Standard (ISO/IEC 9899:1190 (E), 7.10.1.4 The strod function),
      // the 'strtod' function sets 'endptr' to 'nptr' if the conversion fails, so we simply
      // hope that the quadmath function for __float128 honors this convention...
      return (nptr != endptr);
    }
#endif // FEAT_HAVE_QUADMATH
#endif // __CUDACC__
    /// \endcond
  }; // class String

  /// \cond internal
  inline String operator+(const String& a, const String& b)
  {
    return String(a).append(b);
  }

  inline String operator+(const char* a, const String& b)
  {
    return String(a).append(b);
  }

  inline String operator+(const String& a, const char* b)
  {
    return String(a).append(b);
  }

  inline String operator+(const String& a, char c)
  {
    return String(a).append(String::size_type(1), c);
  }

  inline String operator+(char c, const String& b)
  {
    return String(String::size_type(1), c).append(b);
  }
  /// \endcond

  /**
   * \brief Joins a sequence of strings.
   *
   * This functions joins a sequence of items given by two input iterators,
   * separating them by a delimiter string.
   *
   * \param[in] first, last
   * Two forward iterators representing the item sequence.
   * Each item referenced by the iterators is converted via the stringify function.
   *
   * \param[in] delimiter
   * A string acting as a delimiter.
   *
   * \returns A String containing the stringified items.
   */
  template<typename Iterator_>
  String stringify_join(
    Iterator_ first,
    Iterator_ last,
    const String & delimiter = "")
  {
    if(first == last)
      return String();

    Iterator_ it(first);
    String str = stringify(*it);

    for(++it; it != last; ++it)
      str.append(delimiter).append(stringify(*it));

    return str;
  }

  /**
   * \brief Joins a sequence of strings.
   *
   * This method joins a sequence of items given by a container, e.g. <c>std::list<double></c>,
   * separating each item in the container by a delimiter string.
   *
   * \param[in] container
   * A container representing the item sequence.
   * Each item referenced of the container is converted via the stringify function.
   *
   * \param[in] delimiter
   * A string acting as a delimiter.
   *
   * \returns A String containing the stringified items.
   */
  template<typename Container_>
  String stringify_join(
    const Container_& container,
    const String & delimiter = "")
  {
    return stringify_join(container.cbegin(), container.cend(), delimiter);
  }

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

#ifndef __CUDACC__
  inline String stringify(std::nullptr_t)
  {
    return String("nullptr");
  }

#ifdef FEAT_HAVE_QUADMATH
  inline String stringify(__float128 value)
  {
    // get buffer length
    int len = ::quadmath_snprintf(nullptr, 0, "%Qg", value);
    // allocate buffer
    std::vector<char> buffer(len+16);
    // print to buffer
    quadmath_snprintf(buffer.data(), buffer.size(), "%Qg", value);
    // convert buffer to string
    return String(buffer.data());
  }
#endif // FEAT_HAVE_QUADMATH
#endif // __CUDACC__
  /// \endcond

  /**
   * \brief Prints a floating point value to a string in scientific notation.
   *
   * \tparam DataType_
   * The data type of the value to be printed. Is silently assumed to be a standard floating-point datatype,
   * i.e. either \c float, \c double or <c>long double</c>.
   *
   * \note
   * There exists an overload of this function for the \c __float128 data type if FEAT is
   * configured and build with support for the \c quadmath library.
   *
   * \param[in] value
   * The floating-point value that is to be printed.
   *
   * \param[in] precision
   * The precision that is to be used, i.e. the number of mantissa digits to be printed.\n
   * If set to 0, the default (compiler-dependent) precision will be used.
   *
   * \param[in] width
   * The width that is to be used, i.e. the total number of characters to be printed.\n
   * If set to 0, the default (compiler-dependent) width will be used.
   *
   * \param[in] sign
   * Specifies whether non-negative numbers are to be prefixed by '+'.
   *
   * \returns
   * A String containing the value of \p value in scientific notation.
   *
   * \author Peter Zajac
   */
  template<typename DataType_>
  inline String stringify_fp_sci(DataType_ value, int precision = 0, int width = 0, bool sign = false)
  {
    std::ostringstream oss;
    oss << std::scientific;
    if(precision > 0)
      oss << std::setprecision(precision);
    if(width > 0)
      oss << std::setw(width);
    if(sign)
      oss << std::showpos;
    oss << value;
    return oss.str();
  }

  /**
   * \brief Prints a floating point value to a string in fixed-point notation.
   *
   * \tparam DataType_
   * The data type of the value to be printed. Is silently assumed to be a standard floating-point datatype,
   * i.e. either \c float, \c double or <c>long double</c>.
   *
   * \note
   * There exists an overload of this function for the \c __float128 data type if FEAT is
   * configured and build with support for the \c quadmath library.
   *
   * \param[in] value
   * The floating-point value that is to be printed.
   *
   * \param[in] precision
   * The precision that is to be used, i.e. the number of mantissa digits to be printed.\n
   * If set to 0, the default (compiler-dependent) precision will be used.
   *
   * \param[in] width
   * The width that is to be used, i.e. the total number of characters to be printed.\n
   * If set to 0, the default (compiler-dependent) width will be used.
   *
   * \param[in] sign
   * Specifies whether non-negative numbers are to be prefixed by '+'.
   *
   * \returns
   * A String containing the value of \p value in fixed-point notation.
   *
   * \author Peter Zajac
   */
  template<typename DataType_>
  inline String stringify_fp_fix(DataType_ value, int precision = 0, int width = 0, bool sign = false)
  {
    std::ostringstream oss;
    oss << std::fixed;
    if(precision > 0)
      oss << std::setprecision(precision);
    if(width > 0)
      oss << std::setw(width);
    if(sign)
      oss << std::showpos;
    oss << value;
    return oss.str();
  }

#ifndef __CUDACC__
#ifdef FEAT_HAVE_QUADMATH
  inline String stringify_fp_sci(__float128 value, int precision = 0, int width = 0, bool sign = false)
  {
    String format("%");
    if(sign)
      format.append("+");
    if(width > 0)
      format.append(stringify(width));
    if(precision > 0)
    {
      format.append(".");
      format.append(stringify(precision));
    }
    format.append("Qe");
    // get buffer length
    int len = ::quadmath_snprintf(nullptr, 0, format.c_str(), value);
    // allocate buffer
    std::vector<char> buffer(len+16);
    // print to buffer
    quadmath_snprintf(buffer.data(), buffer.size(), format.c_str(), value);
    // convert buffer to string
    return String(buffer.data());
  }

  inline String stringify_fp_fix(__float128 value, int precision = 0, int width = 0, bool sign = false)
  {
    String format("%");
    if(sign)
      format.append("+");
    if(width > 0)
      format.append(stringify(width));
    if(precision > 0)
    {
      format.append(".");
      format.append(stringify(precision));
    }
    format.append("Qf");
    // get buffer length
    int len = ::quadmath_snprintf(nullptr, 0, format.c_str(), value);
    // allocate buffer
    std::vector<char> buffer(len+16);
    // print to buffer
    quadmath_snprintf(buffer.data(), buffer.size(), format.c_str(), value);
    // convert buffer to string
    return String(buffer.data());
  }

  inline std::istream& operator>>(std::istream& is, __float128& x)
  {
    String buffer;
    // try to parse a string
    if(!(is >> buffer).fail())
    {
      if(!buffer.parse(x))
      {
        // parse failed, so put back the string
        for(std::size_t i(0); i < buffer.size(); ++i)
          is.putback(buffer.at(i));

        // set failbit
        is.setstate(std::ios_base::failbit);
      }
    }
    return is;
  }
#endif // FEAT_HAVE_QUADMATH
#endif // __CUDACC__
} // namespace FEAT

// operator<<(__float128) must reside in the std namespace because gcc/icc search the namespace of ostream for a fitting op<< first
// and would otherwise only find op<< with conversions of __float128 to int/double/long etc, which would lead to ambigue overloads, too.
#ifndef __CUDACC__
#ifdef FEAT_HAVE_QUADMATH
namespace std
{
  inline std::ostream& operator<<(std::ostream& os, __float128 x)
  {
    // write to stream
    return (os << FEAT::stringify(x));
  }
}
#endif // FEAT_HAVE_QUADMATH
#endif // __CUDACC__

#endif // KERNEL_UTIL_STRING_HPP
