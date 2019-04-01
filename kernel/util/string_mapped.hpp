// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_STRING_MAPPED_HPP
#define KERNEL_UTIL_STRING_MAPPED_HPP 1

// includes, FEAT
#include <kernel/util/string.hpp>

// includes, system
#include <map>
#include <iostream>

namespace FEAT
{
  /**
   * \brief String-Mapped class template
   *
   * This class implements an auxiliary class that can be used to look up a value of
   * a specific type from a map of strings. The primary use of this class is to parse
   * enumeration values from the set of command line arguments by using the
   * SimpleArgParser class.
   *
   * \tparam ValueType_
   * The type of the value that is to be mapped. Must be copy-assignable.
   *
   * \author Peter Zajac
   */
  template<typename ValueType_>
  class StringMapped
  {
  private:
    /// the value
    ValueType_& _value;
    /// the string map
    const std::map<String, ValueType_>& _value_map;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] value
     * A reference to the value to be looked up.
     *
     * \param[in] value_map
     * A const reference to the value map to be used for lookup.
     */
    explicit StringMapped(ValueType_& value, const std::map<String, ValueType_>& value_map) :
      _value(value),
      _value_map(value_map)
    {
    }

    /**
     * \brief Looks up the value in the map from a given string
     *
     * This function parses the input string by looking up the corresponding value
     * in the underlying string-value map and saves the corresponding value in this
     * object's variable reference. If the string was not found in the underlying map,
     * the internal value of this object is left unmodified.
     *
     * \param[in] str
     * The string that is to be looked up.
     *
     * \returns \c true, if \p string represents a valid key in the underlying map,
     * or \c false, if the \p string was not found in the map.
     */
    bool lookup(const String& str)
    {
      // try to find
      auto it = _value_map.find(str);
      if(it == _value_map.end())
        return false;
      // okay, assign value
      _value = it->second;
      return true;
    }

    /**
     * \brief Input stream extraction operator
     *
     * \param[in,out] is
     * The input stream to extract from.
     *
     * \param[in] val
     * A rvalue-reference to the string mapped value to extract.
     *
     * \returns
     * \c is
     */
    friend std::istream& operator>>(std::istream& is, StringMapped& val)
    {
      String str;
      // Try to fetch a string from the stream
      if(!(is >> str).fail())
      {
        // Try to look up the string; if that fails, set the failbit
        if(!val.lookup(str))
        {
          is.setstate(std::ios::failbit);
        }
      }
      return is;
    }

    /// \cond internal
    friend std::istream& operator>>(std::istream& is, StringMapped&& val)
    {
      String str;
      // Try to fetch a string from the stream
      if(!(is >> str).fail())
      {
        // Try to look up the string; if that fails, set the failbit
        if(!val.lookup(str))
        {
          is.setstate(std::ios::failbit);
        }
      }
      return is;
    }
    /// \endcond

    /**
     * \brief Output stream insertion operator
     *
     * \param[in,out] os
     * The output stream to insert to.
     *
     * \param[in] val
     * A rvalue-reference to the string mapped value to insert.
     *
     * \returns
     * \c os
     */
    friend std::ostream& operator<<(std::ostream& os, const StringMapped& val)
    {
      return (os << (val._value));
    }

    /// \cond internal
    friend std::ostream& operator<<(std::ostream& os, StringMapped&& val)
    {
      return (os << (val._value));
    }
    /// \endcond
  }; // class StringMapped

  /**
   * \brief Creates a StringMapped object.
   *
   * \param[in] value
   * A reference to the value to be looked up.
   *
   * \param[in] value_map
   * A const reference to the value map to be used for lookup.
   *
   * \returns
   * A StringMapped<ValueType_> object.
   */
  template<typename ValueType_>
  inline StringMapped<ValueType_> string_mapped(ValueType_& value, const std::map<String, ValueType_>& value_map)
  {
    return StringMapped<ValueType_>(value, value_map);
  }

  /**
   * \brief Looks up a StringMapped value.
   *
   * \param[in] value
   * A reference to the value to be looked up.
   *
   * \param[in] value_map
   * A const reference to the value map to be used for lookup.
   *
   * \param[in] str
   * The string to be looked up.
   *
   * \returns
   * \c true, if \c str was found in the value map, otherwise \c false.
   */
  template<typename ValueType_>
  inline bool string_mapped_lookup(ValueType_& value, const std::map<String, ValueType_>& value_map, const String& str)
  {
    return string_mapped(value, value_map).lookup(str);
  }
} // namespace FEAT

#endif // KERNEL_UTIL_STRING_MAPPED_HPP
