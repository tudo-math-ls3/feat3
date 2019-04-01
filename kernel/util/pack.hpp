// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_PACK_HPP
#define KERNEL_UTIL_PACK_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/type_traits.hpp>

// includes, system
#include <cstdint>
#include <vector>

// includes, thirdparty
#ifdef FEAT_HAVE_HALFMATH
#include <half.hpp>
#endif // FEAT_HAVE_HALFMATH
#ifdef FEAT_HAVE_ZLIB
#include <zlib.h>
#endif // FEAT_HAVE_ZLIB

namespace FEAT
{
  /**
   * \brief Data Array Pack namespace
   *
   * This namespace encapsulates various functions than provide means of packing and unpacking
   * data arrays as well as compression and decompression of buffers using the third-party
   * library zlib.
   */
  namespace Pack
  {
    /*inline bool little_endian()
    {
      union { u8 c[4]; u32 i; } data;
      data.i = 0x12345678;
      return (data.c[0] == 0x78);
    }*/

    /// 8-bit signed integer type
    typedef std::int8_t i8;
    /// 16-bit signed integer type
    typedef std::int16_t i16;
    /// 32-bit signed integer type
    typedef std::int32_t i32;
    /// 64-bit signed integer type
    typedef std::int64_t i64;

    /// 8-bit unsigned integer type
    typedef std::uint8_t u8;
    /// 16-bit unsigned integer type
    typedef std::uint16_t u16;
    /// 32-bit unsigned integer type
    typedef std::uint32_t u32;
    /// 64-bit unsigned integer type
    typedef std::uint64_t u64;

#if defined(FEAT_HAVE_HALFMATH) || defined(DOXYGEN)
#define FEAT_HAVE_PACK_TYPE_F16 1
    /// 16-bit floating point type
    typedef half_float::half f16;
    // make sure this really has only 16 bits; this may be violated by a stupid compiler
    static_assert(sizeof(f16) == 2, "16-bit float is assumed to have exactly 16 bits");
#endif

    /// 32-bit floating point type
    typedef float f32;
    /// 64-bit floating point type
    typedef double f64;

#if defined(FEAT_HAVE_QUADMATH) || defined(DOYGEN)
#define FEAT_HAVE_PACK_TYPE_F128 1
    /// 128-bit floating piont type
    typedef __float128 f128;
    // make sure this really has only 16 bytes
    static_assert(sizeof(f128) == 16, "128-bit float is assumed to have exactly 128 bits");
#endif

    /**
     * \brief Type enumeration
     */
    enum class Type : u16
    {
      None    = 0x0000, //< none (for comparison)
      Mask_T  = 0x0FFF, //< raw type mask
      Mask_Z  = 0x8000, //< zlib compression mask
    //Mask_E  = 0x4000, //< swap endianness mask

      // raw signed integer types
      I8      = 0x0011, //<  8-bit unsigned integer
      I16     = 0x0012, //< 16-bit unsigned integer
      I32     = 0x0013, //< 32-bit unsigned integer
      I64     = 0x0014, //< 64-bit unsigned integer

      // raw unsigned integer types
      U8      = 0x0021, //<  8-bit unsigned integer
      U16     = 0x0022, //< 16-bit unsigned integer
      U32     = 0x0023, //< 32-bit unsigned integer
      U64     = 0x0024, //< 64-bit unsigned integer

      // raw floating point types
    //F8      = 0x0031, //<   8-bit floating point (quarter precision)
      F16     = 0x0032, //<  16-bit floating point (half precision)
      F32     = 0x0033, //<  32-bit floating point (single precision)
      F64     = 0x0034, //<  64-bit floating point (double precision)
      F128    = 0x0035, //< 128-bit floating point (quadruple precision)

      // zlib compressed signed integer types
      ZI8     = 0x8011, //<  8-bit unsigned integer
      ZI16    = 0x8012, //< 16-bit unsigned integer
      ZI32    = 0x8013, //< 32-bit unsigned integer
      ZI64    = 0x8014, //< 64-bit unsigned integer

      // zlib compressed unsigned integer types
      ZU8     = 0x8021, //<  8-bit unsigned integer
      ZU16    = 0x8022, //< 16-bit unsigned integer
      ZU32    = 0x8023, //< 32-bit unsigned integer
      ZU64    = 0x8024, //< 64-bit unsigned integer

      // zlib compressed floating point types
    //ZF8     = 0x8031, //<   8-bit floating point (quarter precision)
      ZF16    = 0x8032, //<  16-bit floating point (half precision)
      ZF32    = 0x8033, //<  32-bit floating point (single precision)
      ZF64    = 0x8034, //<  64-bit floating point (double precision)
      ZF128   = 0x8035  //< 128-bit floating point (quadruple precision)
    }; // enum class Type

    /// bit-wise AND operator for Pack::Type
    inline Pack::Type operator&(Pack::Type a, Pack::Type b)
    {
      return (Pack::Type)(((Pack::u16)a) & ((Pack::u16)b));
    }

    /// bit-wise OR operator for Pack::Type
    inline Pack::Type operator|(Pack::Type a, Pack::Type b)
    {
      return (Pack::Type)(((Pack::u16)a) | ((Pack::u16)b));
    }

    /// stream output operator for Pack::Type
    inline std::ostream& operator<<(std::ostream& os, Pack::Type t)
    {
      switch(t)
      {
    //case Pack::Type::F8:    return os << "F8";
      case Pack::Type::F16:   return os << "F16";
      case Pack::Type::F32:   return os << "F32";
      case Pack::Type::F64:   return os << "F64";
      case Pack::Type::F128:  return os << "F128";
      case Pack::Type::I8:    return os << "I8";
      case Pack::Type::I16:   return os << "I16";
      case Pack::Type::I32:   return os << "I32";
      case Pack::Type::I64:   return os << "I64";
      case Pack::Type::U8:    return os << "U8";
      case Pack::Type::U16:   return os << "U16";
      case Pack::Type::U32:   return os << "U32";
      case Pack::Type::U64:   return os << "U64";
    //case Pack::Type::ZF8:   return os << "ZF8";
      case Pack::Type::ZF16:  return os << "ZF16";
      case Pack::Type::ZF32:  return os << "ZF32";
      case Pack::Type::ZF64:  return os << "ZF64";
      case Pack::Type::ZF128: return os << "ZF128";
      case Pack::Type::ZI8:   return os << "ZI8";
      case Pack::Type::ZI16:  return os << "ZI16";
      case Pack::Type::ZI32:  return os << "ZI32";
      case Pack::Type::ZI64:  return os << "ZI64";
      case Pack::Type::ZU8:   return os << "ZU8";
      case Pack::Type::ZU16:  return os << "ZU16";
      case Pack::Type::ZU32:  return os << "ZU32";
      case Pack::Type::ZU64:  return os << "ZU64";
      default:
        return os << "???";
      }
    }

    /// stream input operator for Pack::Type
    inline std::istream& operator>>(std::istream& is, Pack::Type& t)
    {
      String s;
      if((is >> s).fail())
        return is;

    //if(s.compare_no_case("F8") == 0)    t = Pack::Type::F8; else
      if(s.compare_no_case("F16") == 0)   t = Pack::Type::F16; else
      if(s.compare_no_case("F32") == 0)   t = Pack::Type::F32; else
      if(s.compare_no_case("F64") == 0)   t = Pack::Type::F64; else
      if(s.compare_no_case("F128") == 0)  t = Pack::Type::F128; else
      if(s.compare_no_case("I8") == 0)    t = Pack::Type::I8; else
      if(s.compare_no_case("I16") == 0)   t = Pack::Type::I16; else
      if(s.compare_no_case("I32") == 0)   t = Pack::Type::I32; else
      if(s.compare_no_case("I64") == 0)   t = Pack::Type::I64; else
      if(s.compare_no_case("U8") == 0)    t = Pack::Type::U8; else
      if(s.compare_no_case("U16") == 0)   t = Pack::Type::U16; else
      if(s.compare_no_case("U32") == 0)   t = Pack::Type::U32; else
      if(s.compare_no_case("U64") == 0)   t = Pack::Type::U64; else
    //if(s.compare_no_case("ZF8") == 0)    t = Pack::Type::ZF8; else
      if(s.compare_no_case("ZF16") == 0)   t = Pack::Type::ZF16; else
      if(s.compare_no_case("ZF32") == 0)   t = Pack::Type::ZF32; else
      if(s.compare_no_case("ZF64") == 0)   t = Pack::Type::ZF64; else
      if(s.compare_no_case("ZF128") == 0)  t = Pack::Type::ZF128; else
      if(s.compare_no_case("ZI8") == 0)    t = Pack::Type::ZI8; else
      if(s.compare_no_case("ZI16") == 0)   t = Pack::Type::ZI16; else
      if(s.compare_no_case("ZI32") == 0)   t = Pack::Type::ZI32; else
      if(s.compare_no_case("ZI64") == 0)   t = Pack::Type::ZI64; else
      if(s.compare_no_case("ZU8") == 0)    t = Pack::Type::ZU8; else
      if(s.compare_no_case("ZU16") == 0)   t = Pack::Type::ZU16; else
      if(s.compare_no_case("ZU32") == 0)   t = Pack::Type::ZU32; else
      if(s.compare_no_case("ZU64") == 0)   t = Pack::Type::ZU64; else
        is.setstate(std::ios_base::failbit);

      return is;
    }

    inline std::size_t element_size(const Pack::Type type)
    {
      return std::size_t(1) << (((int)type & 0xF) - 1);
    }

    /// \cond internal
    namespace Intern
    {
      /// auxiliary helper class: swaps the bytes of an 8/16/32/64/128 byte type
      template<int bytes_>
      struct SwapHelper;

      template<>
      struct SwapHelper<1>
      {
        static void swap(void*)
        {
          // nothing to do here
        }
      };

      template<>
      struct SwapHelper<2>
      {
        static void swap(void* x)
        {
          u16& v = *reinterpret_cast<u16*>(x);
          v = u16(((v & 0x00FF) << 8) | ((v & 0xFF00) >> 8));
        }
      };

      template<>
      struct SwapHelper<4>
      {
        static void swap(void* x)
        {
          u32& v = *reinterpret_cast<u32*>(x);
          v = ((v & 0x000000FF) << 24)
            | ((v & 0x0000FF00) <<  8)
            | ((v & 0x00FF0000) >>  8)
            | ((v & 0xFF000000) >> 24);
        }
      };

      template<>
      struct SwapHelper<8>
      {
        static void swap(void* x)
        {
          u64& v = *reinterpret_cast<u64*>(x);
          v = ((v & 0x00000000000000FFull) << 56)
            | ((v & 0x000000000000FF00ull) << 40)
            | ((v & 0x0000000000FF0000ull) << 24)
            | ((v & 0x00000000FF000000ull) <<  8)
            | ((v & 0x000000FF00000000ull) >>  8)
            | ((v & 0x0000FF0000000000ull) >> 24)
            | ((v & 0x00FF000000000000ull) >> 40)
            | ((v & 0xFF00000000000000ull) >> 56);
        }
      };

      template<>
      struct SwapHelper<16>
      {
        static void swap(void* x)
        {
          // swap two 64-bit ints
          u64* vv = reinterpret_cast<u64*>(x);
          u64 t = vv[0];
          vv[0] = vv[1];
          vv[1] = t;
          // the the bytes within
          SwapHelper<8>::swap(&vv[0]);
          SwapHelper<8>::swap(&vv[1]);
        }
      };

      /**
       * \brief Swaps the bytes of an object.
       */
      template<typename X_>
      inline X_ xswap(X_ x)
      {
        SwapHelper<int(sizeof(X_))>::swap(&x);
        return x;
      }

      template<typename X_, typename T_>
      static std::size_t xencode(void* buf, const T_* t, const std::size_t count, bool swap_bytes)
      {
        X_* x = reinterpret_cast<X_*>(buf);
        if(swap_bytes)
        {
          for(std::size_t i(0); i < count; ++i)
          {
            x[i] = xswap(static_cast<X_>(t[i]));
          }
        }
        else
        {
          for(std::size_t i(0); i < count; ++i)
          {
            x[i] = static_cast<X_>(t[i]);
          }
        }
        return count * sizeof(X_);
      }

      template<typename X_, typename T_>
      static std::size_t xdecode(T_* t, const void* buf, const std::size_t count, bool swap_bytes)
      {
        const X_* x = reinterpret_cast<const X_*>(buf);
        if(swap_bytes)
        {
          for(std::size_t i(0); i < count; ++i)
          {
            t[i] = static_cast<T_>(xswap(x[i]));
          }
        }
        else
        {
          for(std::size_t i(0); i < count; ++i)
          {
            t[i] = static_cast<T_>(x[i]);
          }
        }
        return count * sizeof(X_);
      }

      template<typename Tclass_, bool signed_>
      struct TypeHelper;

      // specialisation for floating point types
      template<>
      struct TypeHelper<FEAT::Type::FloatingClass, true>
      {
        template<typename T_>
        static std::size_t encode(void* buf, const T_* t, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
        {
          switch(pack_type)
          {
#ifdef FEAT_HAVE_PACK_TYPE_F16
          case Pack::Type::F16:  return xencode<Pack::f16>(buf, t, count, swap_bytes);
#endif
          case Pack::Type::F32:  return xencode<Pack::f32>(buf, t, count, swap_bytes);
          case Pack::Type::F64:  return xencode<Pack::f64>(buf, t, count, swap_bytes);
#ifdef FEAT_HAVE_PACK_TYPE_F128
          case Pack::Type::F128: return xencode<Pack::f128>(buf, t, count, swap_bytes);
#endif
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "invalid data type conversion");
          }
        }

        template<typename T_>
        static std::size_t decode(T_* t, const void* buf, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
        {
          switch(pack_type)
          {
#ifdef FEAT_HAVE_PACK_TYPE_F16
          case Pack::Type::F16:  return xdecode<Pack::f16>(t, buf, count, swap_bytes);
#endif
          case Pack::Type::F32:  return xdecode<Pack::f32>(t, buf, count, swap_bytes);
          case Pack::Type::F64:  return xdecode<Pack::f64>(t, buf, count, swap_bytes);
#ifdef FEAT_HAVE_PACK_TYPE_F128
          case Pack::Type::F128: return xdecode<Pack::f128>(t, buf, count, swap_bytes);
#endif
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "invalid data type conversion");
          }
        }
      };

      // specialisation for signed integer types
      template<>
      struct TypeHelper<FEAT::Type::IntegralClass, true>
      {
        template<typename T_>
        static std::size_t encode(void* buf, const T_* t, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
        {
          switch(pack_type)
          {
          case Pack::Type::I8:   return xencode<Pack::i8> (buf, t, count, swap_bytes);
          case Pack::Type::I16:  return xencode<Pack::i16>(buf, t, count, swap_bytes);
          case Pack::Type::I32:  return xencode<Pack::i32>(buf, t, count, swap_bytes);
          case Pack::Type::I64:  return xencode<Pack::i64>(buf, t, count, swap_bytes);
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "invalid data type conversion");
          }
        }

        template<typename T_>
        static std::size_t decode(T_* t, const void* buf, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
        {
          switch(pack_type)
          {
          case Pack::Type::I8:   return xdecode<Pack::i8> (t, buf, count, swap_bytes);
          case Pack::Type::I16:  return xdecode<Pack::i16>(t, buf, count, swap_bytes);
          case Pack::Type::I32:  return xdecode<Pack::i32>(t, buf, count, swap_bytes);
          case Pack::Type::I64:  return xdecode<Pack::i64>(t, buf, count, swap_bytes);
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "invalid data type conversion");
          }
        }
      };

      // specialisation for unsigned integer types
      template<>
      struct TypeHelper<FEAT::Type::IntegralClass, false>
      {
        template<typename T_>
        static std::size_t encode(void* buf, const T_* t, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
        {
          switch(pack_type)
          {
          case Pack::Type::U8:   return xencode<Pack::u8> (buf, t, count, swap_bytes);
          case Pack::Type::U16:  return xencode<Pack::u16>(buf, t, count, swap_bytes);
          case Pack::Type::U32:  return xencode<Pack::u32>(buf, t, count, swap_bytes);
          case Pack::Type::U64:  return xencode<Pack::u64>(buf, t, count, swap_bytes);
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "invalid data type conversion");
          }
        }

        template<typename T_>
        static std::size_t decode(T_* t, const void* buf, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
        {
          switch(pack_type)
          {
          case Pack::Type::U8:   return xdecode<Pack::u8> (t, buf, count, swap_bytes);
          case Pack::Type::U16:  return xdecode<Pack::u16>(t, buf, count, swap_bytes);
          case Pack::Type::U32:  return xdecode<Pack::u32>(t, buf, count, swap_bytes);
          case Pack::Type::U64:  return xdecode<Pack::u64>(t, buf, count, swap_bytes);
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "invalid data type conversion");
          }
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Encodes an array into a packed buffer without compression support.
     *
     * \attention
     * This function does not encode compressed pack types.
     *
     * \note
     * This function exists primarily for code outsourcing; it is recommended
     * that you use the encode() function instead.
     *
     * \param[out] buf
     * A pointer to the output buffer. Must not be \c nullptr.
     *
     * \param[in] src
     * A pointer to the input array that is to be packed. Must not be \c nullptr.
     *
     * \param[in] count
     * The number of input array elements.
     *
     * \param[in] pack_type
     * The desired output buffer pack type.
     *
     * \param[in] swap_bytes
     * Specifies whether to swap the pack type bytes.
     *
     * \returns
     * The total number of bytes written into the output buffer.
     */
    template<typename T_>
    static std::size_t encode_raw(void* buf, const T_* src, const std::size_t count,
      const Pack::Type pack_type, bool swap_bytes)
    {
      if(count <= std::size_t(0))
        return std::size_t(0);

      XASSERT(buf != nullptr);
      XASSERT(src != nullptr);

      XASSERTM((pack_type & Pack::Type::Mask_T) != Pack::Type::None, "invalid pack type");
      XASSERTM((pack_type & Pack::Type::Mask_Z) != Pack::Type::Mask_Z, "cannot encode compressed type");

      typedef typename FEAT::Type::Traits<T_> TypeTraits;

      return Intern::TypeHelper<typename TypeTraits::TypeClass, TypeTraits::is_signed>::
        encode(buf, src, count, pack_type, swap_bytes);
    }

    /**
     * \brief Decodes an array from a packed buffer without compression support.
     *
     * \attention
     * This function does not decode compressed pack types.
     *
     * \note
     * This function exists primarily for code outsourcing; it is recommended
     * that you use the encode() function instead.
     *
     * \param[out] dst
     * A pointer to the output array that is to be unpacked. Must not be \c nullptr.
     *
     * \param[in] buf
     * A pointer to the input buffer. Must not be \c nullptr.
     *
     * \param[in] count
     * The number of output array elements.
     *
     * \param[in] pack_type
     * The input buffer pack type.
     *
     * \param[in] swap_bytes
     * Specifies whether to swap the pack type bytes.
     *
     * \returns
     * The total number of bytes consumed from the input buffer.
     */
    template<typename T_>
    static std::size_t decode_raw(T_* dst, const void* buf, const std::size_t count,
      const Pack::Type pack_type, bool swap_bytes)
    {
      if(count <= std::size_t(0))
        return std::size_t(0);

      XASSERT(dst != nullptr);
      XASSERT(buf != nullptr);

      XASSERTM((pack_type & Pack::Type::Mask_T) != Pack::Type::None, "invalid pack type");
      XASSERTM((pack_type & Pack::Type::Mask_Z) != Pack::Type::Mask_Z, "cannot decode compressed type");

      typedef typename FEAT::Type::Traits<T_> TypeTraits;

      return Intern::TypeHelper<typename TypeTraits::TypeClass, TypeTraits::is_signed>::
        decode(dst, buf, count, pack_type, swap_bytes);
    }

    /**
     * \brief Computes the estimated (upper bound) pack buffer size for an array
     *
     * \param[in] count
     * The number of array elements to be packed.
     *
     * \param[in] type
     * The packed type of the array elements.
     *
     * \returns
     * An estimated (upper bound) array buffer size.
     */
    inline std::size_t estimate_size(const std::size_t count, const Pack::Type type)
    {
      // compute raw element size
      std::size_t raw_size = count * Pack::element_size(type & Pack::Type::Mask_T);

#ifdef FEAT_HAVE_ZLIB
      if((type & Pack::Type::Mask_Z) != Pack::Type::None)
      {
        return std::size_t(::compressBound(uLong(raw_size)));
      }
#else // no FEAT_HAVE_ZLIB
      XASSERTM((type & Pack::Type::Mask_Z) == Pack::Type::None, "cannot estimate compressed size; zlib not available");
#endif // FEAT_HAVE_ZLIB

      // return raw size
      return raw_size;
    }

    /**
     * \brief Encodes an array into a packed buffer
     *
     * \param[out] buf
     * A pointer to the output buffer. Must not be \c nullptr.
     *
     * \param[in] src
     * A pointer to the input array that is to be packed. Must not be \c nullptr.
     *
     * \param[in] buf_size
     * The size of the output buffer \p buf in bytes. Must be at least as big as the estimated buffer size.
     *
     * \param[in] count
     * The number of input array elements.
     *
     * \param[in] pack_type
     * The desired output buffer pack type.
     *
     * \param[in] swap_bytes
     * Specifies whether to swap the pack type bytes.
     *
     * \returns
     * The total number of bytes written into the output buffer.
     */
    template<typename T_>
    static std::size_t encode(void* buf, const T_* src, const std::size_t buf_size, const std::size_t count,
      const Pack::Type pack_type, bool swap_bytes)
    {
      if(count <= std::size_t(0))
        return std::size_t(0);

      XASSERT(buf != nullptr);
      XASSERT(src != nullptr);

      // get raw type
      const Pack::Type raw_type = pack_type & Pack::Type::Mask_T;

      // no compression required?
      if((pack_type & Pack::Type::Mask_Z) == Pack::Type::None)
        return encode_raw(buf, src, count, raw_type, swap_bytes);

#ifdef FEAT_HAVE_ZLIB
      // ensure that the temporary buffer is large enough
      std::size_t raw_size = element_size(raw_type) * count;
      std::vector<char> tmp(raw_size);

      // encode into temporary buffer
      encode_raw(tmp.data(), src, count, raw_type, swap_bytes);

      // cast buffer lengths to zlib types
      uLongf dl = static_cast<uLongf>(buf_size);
      uLong  sl = static_cast<uLong> (raw_size);

      // cast buffer pointers to zlib types
            Bytef* dbuf = reinterpret_cast<      Bytef*>(buf);
      const Bytef* sbuf = reinterpret_cast<const Bytef*>(tmp.data());

      // compress via zlib
      if(::compress(dbuf, &dl, sbuf, sl) != Z_OK)
        throw INTERNAL_ERROR("zlib compression error");

      // return final destination buffer usage as reported by zlib
      return static_cast<std::size_t>(dl);
#else // no FEAT_HAVE_ZLIB
      (void)buf_size;
      throw InternalError(__func__, __FILE__ , __LINE__, "cannot encode compressed type; zlib not available");
#endif // FEAT_HAVE_ZLIB
    }

    /**
     * \brief Decodes an array from a packed buffer
     *
     * \param[out] dst
     * A pointer to the output array that is to be unpacked. Must not be \c nullptr.
     *
     * \param[in] buf
     * A pointer to the input buffer. Must not be \c nullptr.
     *
     * \param[in] count
     * The number of output array elements.
     *
     * \param[in] buf_size
     * The size of the input buffer \p buf in bytes.
     *
     * \param[in] pack_type
     * The input buffer pack type.
     *
     * \param[in] swap_bytes
     * Specifies whether to swap the pack type bytes.
     *
     * \returns
     * The total number of bytes consumed from the input buffer.
     */
    template<typename T_>
    static std::size_t decode(T_* dst, const void* buf, const std::size_t count, const std::size_t buf_size,
      const Pack::Type pack_type, bool swap_bytes)
    {
      if(count <= std::size_t(0))
        return std::size_t(0);

      XASSERT(dst != nullptr);
      XASSERT(buf != nullptr);

      // get raw type
      const Pack::Type raw_type = pack_type & Pack::Type::Mask_T;

      // no decompression required?
      if((pack_type & Pack::Type::Mask_Z) == Pack::Type::None)
        return decode_raw(dst, buf, count, raw_type, swap_bytes);

#ifdef FEAT_HAVE_ZLIB
      // ensure that the temporary buffer is large enough
      std::size_t raw_size = element_size(raw_type) * count + 1024u; // + 1KB
      std::vector<char> tmp(raw_size);

      // cast buffer lengths to zlib types
      uLongf dl = static_cast<uLongf>(raw_size);
      uLong  sl = static_cast<uLong> (buf_size);

      // cast buffer pointers to zlib types
            Bytef* dbuf = reinterpret_cast<      Bytef*>(tmp.data());
      const Bytef* sbuf = reinterpret_cast<const Bytef*>(buf);

      // decompress via zlib
      if(::uncompress2(dbuf, &dl, sbuf, &sl) != Z_OK)
        throw INTERNAL_ERROR("zlib decompression error");

      // decode from temporary buffer
      decode_raw(dst, tmp.data(), count, raw_type, swap_bytes);

      // return final source buffer usage as reported by zlib
      return static_cast<std::size_t>(sl);
#else // no FEAT_HAVE_ZLIB
      (void)buf_size;
      throw InternalError(__func__, __FILE__ , __LINE__, "cannot decode compressed type; zlib not available");
#endif // FEAT_HAVE_ZLIB
    }
  } // namespace Pack
} // namespace FEAT

#endif // KERNEL_UTIL_PACK_HPP
