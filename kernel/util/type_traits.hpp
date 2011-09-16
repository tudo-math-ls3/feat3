#pragma once
#ifndef KERNEL_UTIL_TYPE_TRAITS_HPP
#define KERNEL_UTIL_TYPE_TRAITS_HPP 1

// includes, FEAST
#include <kernel/util/string_utils.hpp>

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
    static std::string name()
    {
      return stringify(typeid(DT_).name());
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
    static std::string name()
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
    static std::string name()
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
    static std::string name()
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
    static std::string name()
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
    static std::string name()
    {
      return "int";
    }
  };

  /**
  * \brief TypeTraits specialisation for Nil tag class
  *
  * \author Dirk Ribbrock
  */
  template<>
  struct TypeTraits<Nil>
  {
    /// returns a string identifying the datatype
    static std::string name()
    {
      return "Nil";
    }
  };
}

#endif // KERNEL_UTIL_TYPE_TRAITS_HPP
