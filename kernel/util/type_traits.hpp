#pragma once
#ifndef KERNEL_UTIL_TYPE_TRAITS_HPP
#define KERNEL_UTIL_TYPE_TRAITS_HPP 1

// includes, FEAST
#include <kernel/util/string.hpp>

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
      return DT_::name();
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>float</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<float>
  {
    /// dummy enum
    enum
    {
      /// this type is not integral
      is_int = 0,
      /// this type is floating
      is_float = 1,
      /// this type is not boolean
      is_bool = 0,
      /// this type is signed
      is_signed = 1
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "float";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>double</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<double>
  {
    /// dummy enum
    enum
    {
      /// this type is not integral
      is_int = 0,
      /// this type is floating
      is_float = 1,
      /// this type is not boolean
      is_bool = 0,
      /// this type is signed
      is_signed = 1
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "double";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>long double</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<long double>
  {
    /// dummy enum
    enum
    {
      /// this type is not integral
      is_int = 0,
      /// this type is floating
      is_float = 1,
      /// this type is not boolean
      is_bool = 0,
      /// this type is signed
      is_signed = 1
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "long double";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>unsigned int</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<unsigned int>
  {
    /// dummy enum
    enum
    {
      /// this type is  integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is unsigned
      is_signed = 0
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "unsigned int";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>signed int</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<signed int>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is signed
      is_signed = 1
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "signed int";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>unsigned char</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<unsigned char>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is unsigned
      is_signed = 0
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "unsigned char";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>signed char</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<signed char>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is signed
      is_signed = 1
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "signed char";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>unsigned short</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<unsigned short>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is unsigned
      is_signed = 0
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "unsigned short";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>signed short</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<signed short>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is signed
      is_signed = 1
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "signed short";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>unsigned long</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<unsigned long>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is unsigned
      is_signed = 0
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "unsigned long";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>signed long</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<signed long>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is signed
      is_signed = 1
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "signed long";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>unsigned long long</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<unsigned long long>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is unsigned
      is_signed = 0
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "unsigned long long";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>signed long long</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<signed long long>
  {
    /// dummy enum
    enum
    {
      /// this type is integral
      is_int = 1,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 0,
      /// this type is signed
      is_signed = 1
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "signed long long";
    }
  };

  /**
   * \brief TypeTraits specialisation for <c>bool</c>
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<bool>
  {
    /// dummy enum
    enum
    {
      /// this type is not integral
      is_int = 0,
      /// this type is not floating
      is_float = 0,
      /// this type is not boolean
      is_bool = 1,
      /// this type is unsigned
      is_signed = 0
    };

    /// returns a string identifying the datatype
    static String name()
    {
      return "bool";
    }
  };

  /**
   * \brief TypeTraits specialisation for Nil tag class
   *
   * \author Dirk Ribbrock
   * \author Peter Zajac
   */
  template<>
  struct TypeTraits<Nil>
  {
    /// returns a string identifying the datatype
    static String name()
    {
      return "Nil";
    }
  };
} // FEAST

#endif // KERNEL_UTIL_TYPE_TRAITS_HPP
