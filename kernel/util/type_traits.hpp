#pragma once
#ifndef KERNEL_UTIL_TYPE_TRAITS_HPP
#define KERNEL_UTIL_TYPE_TRAITS_HPP 1

// includes, FEAST
#include <kernel/util/string.hpp>

namespace FEAST
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

      /// this type is of floating class
      typedef FloatingClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "float";
      }

      /// returns the items value in double precision
      static double to_double(float val)
      {
        return (double)val;
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

      /// this type is of floating class
      typedef FloatingClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "double";
      }

      /// returns the items value in double precision
      static double to_double(double val)
      {
        return val;
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

      /// this type is of integral class
      typedef IntegralClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "unsigned long";
      }

      /// returns the items value in double precision
      static double to_double(unsigned long val)
      {
        return (double)val;
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

      /// this type is of boolean class
      typedef BooleanClass TypeClass;

      /// returns a string identifying the datatype
      static String name()
      {
        return "bool";
      }
    };
  } // namespace Type
} // namespace FEAST

#endif // KERNEL_UTIL_TYPE_TRAITS_HPP
