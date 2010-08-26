#pragma once
#ifndef UTIL_TYPE_TRAITS_HHP
/// Header guard
#define UTIL_TYPE_TRAITS_HHP 1

#include <typeinfo>
#include <kernel/util/stringify.hpp>

namespace Feast
{
  /**
   * \brief Basic TypeTraits struct
   *
   * \author Dirk Ribbrock
   */
  template <typename DT_>
  struct TypeTraits
  {
    /// Return a string idetifying the Datatype
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
  template <>
  struct TypeTraits<float>
  {
    /// Return a string idetifying the Datatype
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
  template <>
  struct TypeTraits<double>
  {
    /// Return a string idetifying the Datatype
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
  template <>
  struct TypeTraits<unsigned long>
  {
    /// Return a string idetifying the Datatype
    static std::string name()
    {
      return "unsigned long";
    }
  };
}

#endif //UTIL_TYPE_TRAITS_HHP
