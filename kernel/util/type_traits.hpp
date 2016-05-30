#pragma once
#ifndef KERNEL_UTIL_TYPE_TRAITS_HPP
#define KERNEL_UTIL_TYPE_TRAITS_HPP 1

// includes, FEAT
#include <kernel/util/string.hpp>

namespace FEAT
{
  /**
   * \brief Type namespace
   *
   * This namespace encapsulates the TypeTraits class template as well as tag classes used by it.
   */
  namespace Type
  {
    /// Tag class for any data type not matching any other type class
    class AnotherClass {};
    /// Tag class for integral data types
    class IntegralClass {};
    /// Tag class for floating data types
    class FloatingClass {};
    /// Tag class for the one and only boolean data type
    class BooleanClass {};

    /**
     * \brief basic Type Traits struct
     *
     * \author Dirk Ribbrock
     */
    template<typename DT_>
    struct Traits
    {
      /// this type is of another class
      typedef AnotherClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return DT_::name();
      }
    };

    /**
     * \brief Type Traits specialisation for <c>float</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<float>
    {
      /// this type is not integral
      static constexpr bool is_int = false;
      /// this type is floating
      static constexpr bool is_float = true;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of floating class
      typedef FloatingClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "float";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>double</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<double>
    {
      /// this type is not integral
      static constexpr bool is_int = false;
      /// this type is floating
      static constexpr bool is_float = true;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of floating class
      typedef FloatingClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "double";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>long double</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<long double>
    {
      /// this type is not integral
      static constexpr bool is_int = false;
      /// this type is floating
      static constexpr bool is_float = true;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of floating class
      typedef FloatingClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "long double";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>unsigned int</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<unsigned int>
    {
      /// this type is  integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is unsigned
      static constexpr bool is_signed = false;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "unsigned int";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>signed int</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<signed int>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "signed int";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>unsigned char</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<unsigned char>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is unsigned
      static constexpr bool is_signed = false;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "unsigned char";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>signed char</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<signed char>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "signed char";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>unsigned short</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<unsigned short>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is unsigned
      static constexpr bool is_signed = false;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "unsigned short";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>signed short</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<signed short>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "signed short";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>unsigned long</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<unsigned long>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is unsigned
      static constexpr bool is_signed = false;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "unsigned long";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>signed long</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<signed long>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "signed long";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>unsigned long long</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<unsigned long long>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is unsigned
      static constexpr bool is_signed = false;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "unsigned long long";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>signed long long</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<signed long long>
    {
      /// this type is integral
      static constexpr bool is_int = true;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "signed long long";
      }
    };

    /**
     * \brief Type Traits specialisation for <c>bool</c>
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<>
    struct Traits<bool>
    {
      /// this type is not integral
      static constexpr bool is_int = false;
      /// this type is not floating
      static constexpr bool is_float = false;
      /// this type is not boolean
      static constexpr bool is_bool = true;
      /// this type is unsigned
      static constexpr bool is_signed = false;

      /// this type is of boolean class
      typedef BooleanClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "bool";
      }
    };

#ifdef FEAT_HAVE_QUADMATH
    /**
     * \brief Type Traits specialisation for <c>__float128</c>
     *
     * \author Peter Zajac
     */
    template<>
    struct Traits<__float128>
    {
      /// this type is not integral
      static constexpr bool is_int = false;
      /// this type is floating
      static constexpr bool is_float = true;
      /// this type is not boolean
      static constexpr bool is_bool = false;
      /// this type is signed
      static constexpr bool is_signed = true;

      /// this type is of floating class
      typedef FloatingClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "__float128";
      }
    };
#endif // FEAT_HAVE_QUADMATH
  } // namespace Type
} // namespace FEAT

#endif // KERNEL_UTIL_TYPE_TRAITS_HPP
