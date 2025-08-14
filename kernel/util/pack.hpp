// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/util/string.hpp>

// includes, system
#include <cstdint>
#include <ostream>
#include <istream>

/**
 * \brief Data Array Pack namespace
 *
 * This namespace encapsulates various functions than provide means of packing
 * and unpacking data arrays as well as compression and decompression of buffers
 * using the third-party library zlib or lossy compression and decompression of
 * buffers using the third-party library zfp.
 */
namespace FEAT::Pack
{
  /// 8-bit signed integer type
  using i8 = std::int8_t;
  /// 16-bit signed integer type
  using i16 = std::int16_t;
  /// 32-bit signed integer type
  using i32 = std::int32_t;
  /// 64-bit signed integer type
  using i64 = std::int64_t;

  /// 8-bit unsigned integer type
  using u8 = std::uint8_t;
  /// 16-bit unsigned integer type
  using u16 = std::uint16_t;
  /// 32-bit unsigned integer type
  using u32 = std::uint32_t;
  /// 64-bit unsigned integer type
  using u64 = std::uint64_t;

#if defined(FEAT_HAVE_HALFMATH) || defined(DOXYGEN)
#define FEAT_HAVE_PACK_TYPE_F16 1
  /// 16-bit floating point type
  using f16 = Half;
  // make sure this really has only 16 bits; this may be violated by a stupid
  // compiler
  static_assert(sizeof(f16) == 2, "16-bit float is assumed to have exactly 16 bits");
#endif

  /// 32-bit floating point type
  using f32 = float;
  /// 64-bit floating point type
  using f64 = double;

#if defined(FEAT_HAVE_QUADMATH) || defined(DOYGEN)
#define FEAT_HAVE_PACK_TYPE_F128 1
  /// 128-bit floating piont type
  using f128 = __float128;
  // make sure this really has only 16 bytes
  static_assert(sizeof(f128) == 16, "128-bit float is assumed to have exactly 128 bits");
#endif

/**
 * \brief bitmask for zfp header
 *
 * //0x1u HEADER_MAGIC -> for version control--> only checks codec... in last 8
 * bits of first 32 bits...
 * //0x2u HEADER_META
 * //0x4u HEADER_MODE
 */
#ifdef FEAT_HAVE_ZFP
  static constexpr uint zfp_header_mask = 7u;
#endif // FEAT_HAVE_ZFP
  /**
   * \brief Type enumeration
   */
  enum class Type : u16
  {
    None = 0x0000,   //< none (for comparison)
    Mask_T = 0x0FFF, //< raw type mask
    Mask_Z = 0x8000, //< zlib compression mask
    Mask_P = 0x2000, //< zfp compression mask
                     // Mask_E  = 0x4000, //< swap endianness mask

    Type_I = 0x0010, //< signed integer type
    Type_U = 0x0020, //< unsigend integer type
    Type_F = 0x0030, //< floating point type

    // raw signed integer types
    I8 = 0x0011,  //<  8-bit unsigned integer
    I16 = 0x0012, //< 16-bit unsigned integer
    I32 = 0x0013, //< 32-bit unsigned integer
    I64 = 0x0014, //< 64-bit unsigned integer

    // raw unsigned integer types
    U8 = 0x0021,  //<  8-bit unsigned integer
    U16 = 0x0022, //< 16-bit unsigned integer
    U32 = 0x0023, //< 32-bit unsigned integer
    U64 = 0x0024, //< 64-bit unsigned integer

    // raw floating point types
    // F8 = 0x0031, //<   8-bit floating point (quarter precision)
    F16 = 0x0032,  //<  16-bit floating point (half precision)
    F32 = 0x0033,  //<  32-bit floating point (single precision)
    F64 = 0x0034,  //<  64-bit floating point (double precision)
    F128 = 0x0035, //< 128-bit floating point (quadruple precision)

    // zlib compressed signed integer types
    ZI8 = 0x8011,  //<  8-bit unsigned integer
    ZI16 = 0x8012, //< 16-bit unsigned integer
    ZI32 = 0x8013, //< 32-bit unsigned integer
    ZI64 = 0x8014, //< 64-bit unsigned integer

    // zlib compressed unsigned integer types
    ZU8 = 0x8021,  //<  8-bit unsigned integer
    ZU16 = 0x8022, //< 16-bit unsigned integer
    ZU32 = 0x8023, //< 32-bit unsigned integer
    ZU64 = 0x8024, //< 64-bit unsigned integer

    // zlib compressed floating point types
    // ZF8 = 0x8031, //<   8-bit floating point (quarter precision)
    ZF16 = 0x8032,  //<  16-bit floating point (half precision)
    ZF32 = 0x8033,  //<  32-bit floating point (single precision)
    ZF64 = 0x8034,  //<  64-bit floating point (double precision)
    ZF128 = 0x8035, //< 128-bit floating point (quadruple precision)

    // zfp compressed floating point types
    // PF8 = 0x2031, //<   8-bit floating point (quarter precision)
    // PF16 = 0x2032, //<  16-bit floating point (half precision)
    PF32 = 0x2033, //<  32-bit floating point (single precision)
    PF64 = 0x2034  //<  64-bit floating point (double precision)
    // PF128 = 0x2035  //< 128-bit floating point (quadruple precision)

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
      // case Pack::Type::F8:    return os << "F8";
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
    // case Pack::Type::ZF8:   return os << "ZF8";
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
    case Pack::Type::PF32:  return os << "PF32";
    case Pack::Type::PF64:  return os << "PF64";
    default:                return os << "???";
    }
  }

  /// stream input operator for Pack::Type
  inline std::istream& operator>>(std::istream& is, Pack::Type& t)
  {
    String s;
    if((is >> s).fail())
    {
      return is;
    }

    // if(s.compare_no_case("F8") == 0)    t = Pack::Type::F8; else
    if(s.compare_no_case("F16") == 0)
      t = Pack::Type::F16;
    else if(s.compare_no_case("F32") == 0)
      t = Pack::Type::F32;
    else if(s.compare_no_case("F64") == 0)
      t = Pack::Type::F64;
    else if(s.compare_no_case("F128") == 0)
      t = Pack::Type::F128;
    else if(s.compare_no_case("I8") == 0)
      t = Pack::Type::I8;
    else if(s.compare_no_case("I16") == 0)
      t = Pack::Type::I16;
    else if(s.compare_no_case("I32") == 0)
      t = Pack::Type::I32;
    else if(s.compare_no_case("I64") == 0)
      t = Pack::Type::I64;
    else if(s.compare_no_case("U8") == 0)
      t = Pack::Type::U8;
    else if(s.compare_no_case("U16") == 0)
      t = Pack::Type::U16;
    else if(s.compare_no_case("U32") == 0)
      t = Pack::Type::U32;
    else if(s.compare_no_case("U64") == 0)
      t = Pack::Type::U64;
    else
      // if(s.compare_no_case("ZF8") == 0)    t = Pack::Type::ZF8; else
    if(s.compare_no_case("ZF16") == 0)
      t = Pack::Type::ZF16;
    else if(s.compare_no_case("ZF32") == 0)
      t = Pack::Type::ZF32;
    else if(s.compare_no_case("ZF64") == 0)
      t = Pack::Type::ZF64;
    else if(s.compare_no_case("ZF128") == 0)
      t = Pack::Type::ZF128;
    else if(s.compare_no_case("ZI8") == 0)
      t = Pack::Type::ZI8;
    else if(s.compare_no_case("ZI16") == 0)
      t = Pack::Type::ZI16;
    else if(s.compare_no_case("ZI32") == 0)
      t = Pack::Type::ZI32;
    else if(s.compare_no_case("ZI64") == 0)
      t = Pack::Type::ZI64;
    else if(s.compare_no_case("ZU8") == 0)
      t = Pack::Type::ZU8;
    else if(s.compare_no_case("ZU16") == 0)
      t = Pack::Type::ZU16;
    else if(s.compare_no_case("ZU32") == 0)
      t = Pack::Type::ZU32;
    else if(s.compare_no_case("ZU64") == 0)
      t = Pack::Type::ZU64;
    else if(s.compare_no_case("PF32") == 0)
      t = Pack::Type::PF32;
    else if(s.compare_no_case("PF64") == 0)
      t = Pack::Type::PF64;
    else
      is.setstate(std::ios_base::failbit);

    return is;
  }

  /**
   * \brief Returns the size of a Pack::Type element in bytes
   *
   * \param[in] type
   * The type whose element size is to be determined
   *
   * \returns
   * The size of the type in bytes or 0, if \p type is not a valid data type
   */
  inline std::size_t element_size(const Pack::Type type)
  {
    // Lower 4 bits contain the size of the data type.
    // Upper bits contain the type of compression.
    constexpr auto size_mask = 0xFULL;
    return std::size_t(1) << ((static_cast<u16>(type) & size_mask) - 1);
  }

  /**
   * \brief Deduct the (raw) Pack::Type from a given data type \p T_
   *
   * Instantiated for i{8,16,32,64}, u{8,16,32,64}, f{16,32,64,128}.
   *
   * \tparam T_
   * The type whose Pack::Type value is to be determined.
   *
   * \returns
   * The (raw) Pack::Type value for the given type \p T_ or Pack::Type::None,
   * if \p T_ does not represent a type from the Pack::Type enum.
   */
  template<typename T_>
  Pack::Type deduct_type();

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
  std::size_t estimate_size(std::size_t count, Pack::Type type, double tolerance = 1e-16);

  /**
   * \brief Encodes an array into a packed buffer
   *
   * Instantiated for i{8,16,32,64}, u{8,16,32,64}, f{16,32,64,128}.
   *
   * \param[out] buf
   * A pointer to the output buffer. Must not be \c nullptr.
   *
   * \param[in] src
   * A pointer to the input array that is to be packed. Must not be \c nullptr.
   *
   * \param[in] buf_size
   * The size of the output buffer \p buf in bytes. Must be at least as big as the
   * estimated buffer size.
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
   * \param[in] tolerance
   * Optional paramter: The desired maximum error for lossy compression.
   *
   * \returns
   * The total number of bytes written into the output buffer.
   */
  template<typename T_>
  std::size_t encode(
    void* buf,
    const T_* src,
    std::size_t buf_size,
    std::size_t count,
    Pack::Type pack_type,
    bool swap_bytes,
    double tolerance = 1e-16);

  /**
   * \brief Decodes an array from a packed buffer
   *
   * Instantiated for i{8,16,32,64}, u{8,16,32,64}, f{16,32,64,128}.
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
  std::size_t decode(
    T_* dst,
    void* buf,
    std::size_t count,
    std::size_t buf_size,
    Pack::Type pack_type,
    bool swap_bytes);
} // namespace FEAT::Pack
