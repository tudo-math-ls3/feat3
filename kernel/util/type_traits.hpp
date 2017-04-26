#pragma once
#ifndef KERNEL_UTIL_TYPE_TRAITS_HPP
#define KERNEL_UTIL_TYPE_TRAITS_HPP 1

// includes, FEAT
#include <kernel/util/string.hpp>

#include <typeindex>

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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(float)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(double)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(long double)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(unsigned int)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(signed int)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(unsigned char)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(signed char)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(unsigned short)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(signed short)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(unsigned long)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(signed long)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(unsigned long long)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(signed long long)).hash_code();
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

      static uint64_t hash_code()
      {
        return (uint64_t)std::type_index(typeid(bool)).hash_code();
      }
    };

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
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

      static uint64_t hash_code()
      {
        return uint64_t(12345);
      }
    };
#endif // FEAT_HAVE_QUADMATH && !__CUDA__CC
  } // namespace Type
} // namespace FEAT

#endif // KERNEL_UTIL_TYPE_TRAITS_HPP
