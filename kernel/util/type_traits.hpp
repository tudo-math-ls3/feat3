#pragma once
#ifndef KERNEL_UTIL_TYPE_TRAITS_HPP
#define KERNEL_UTIL_TYPE_TRAITS_HPP 1

// includes, FEAST
#include <kernel/util/string.hpp>

// includes, system
#include <typeinfo>

namespace FEAST
{
  /**
  * \brief basic TypeTraits struct
  *
  * \author Dirk Ribbrock
  */
  template<typename DT_>
  struct TypeTraits
  {
    /// returns a string identifying the datatype
    static String name()
    {
      //return stringify(typeid(DT_).name());
      return DT_::name;
    }
  };

  /**
  * \brief TypeTraits specialisation for float
  *
  * \author Dirk Ribbrock
  */
  template<>
  struct TypeTraits<float>
  {
    /// returns a string identifying the datatype
    static String name()
    {
      return "float";
    }
  };

  /**
  * \brief TypeTraits specialisation for double
  *
  * \author Dirk Ribbrock
  */
  template<>
  struct TypeTraits<double>
  {
    /// returns a string identifying the datatype
    static String name()
    {
      return "double";
    }
  };

  /**
  * \brief TypeTraits specialisation for unsigned long
  *
  * \author Dirk Ribbrock
  */
  template<>
  struct TypeTraits<unsigned long>
  {
    /// returns a string identifying the datatype
    static String name()
    {
      return "unsigned long";
    }
  };

  /**
  * \brief TypeTraits specialisation for unsigned int
  *
  * \author Dirk Ribbrock
  */
  template<>
  struct TypeTraits<unsigned int>
  {
    /// returns a string identifying the datatype
    static String name()
    {
      return "unsigned int";
    }
  };

  /**
  * \brief TypeTraits specialisation for signed int
  *
  * \author Dirk Ribbrock
  */
  template<>
  struct TypeTraits<int>
  {
    /// returns a string identifying the datatype
    static String name()
    {
      return "int";
    }
  };
} // FEAST

#endif // KERNEL_UTIL_TYPE_TRAITS_HPP
