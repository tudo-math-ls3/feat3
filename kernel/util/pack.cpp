// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/pack.hpp>

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/type_traits.hpp>

// includes, system
#include <cstring>

// includes, thirdparty
#ifdef FEAT_HAVE_ZLIB
#include <zlib.h>
#endif // FEAT_HAVE_ZLIB
#ifdef FEAT_HAVE_ZFP
FEAT_DISABLE_WARNINGS
#include "zfp.h"
FEAT_RESTORE_WARNINGS
#endif // FEAT_HAVE_ZFP

namespace FEAT::Pack
{
  /// auxiliary helper class: swaps the bytes of an 8/16/32/64/128 byte type
  template<int bytes_>
  struct SwapHelper;

  template<>
  struct SwapHelper<1>
  {
    static void swap(void* /*unused*/)
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
      v = u16(((v & 0x00FFULL) << 8) | ((v & 0xFF00ULL) >> 8));
    }
  };

  template<>
  struct SwapHelper<4>
  {
    static void swap(void* x)
    {
      u32& v = *reinterpret_cast<u32*>(x);
      v = u32(((v & 0x000000FFULL) << 24) | ((v & 0x0000FF00ULL) <<  8) |
              ((v & 0x00FF0000ULL) >>  8) | ((v & 0xFF000000ULL) >> 24));
    }
  };

  template<>
  struct SwapHelper<8>
  {
    static void swap(void* x)
    {
      u64& v = *reinterpret_cast<u64*>(x);
      v = ((v & 0x00000000000000FFULL) << 56) | ((v & 0x000000000000FF00ULL) << 40) |
          ((v & 0x0000000000FF0000ULL) << 24) | ((v & 0x00000000FF000000ULL) << 8) |
          ((v & 0x000000FF00000000ULL) >> 8) | ((v & 0x0000FF0000000000ULL) >> 24) |
          ((v & 0x00FF000000000000ULL) >> 40) | ((v & 0xFF00000000000000ULL) >> 56);
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

  /**
   * \brief Encodes an array by converting each of its elements to a desired type
   *
   * Effectively, this function performs <c>((X_*)buf)[i] = X_(src[i])</c> for all
   * i < count.
   *
   * \tparam X_
   * The desired type of the encoded destination array buffer \p buf
   *
   * \tparam T_
   * The type of the source array \p src whose elements are to be converted.
   *
   * \param[out] buf
   * The destination buffer that receives the array elements of type \p X_.
   * This buffer is assumed to be at least <c>sizeof(X_)*count</c> bytes in size.
   *
   * \param[in] src
   * The source array whose elements of type \p T_ are to be converted
   *
   * \param[in] count
   * The number of elements in the source array \p src.
   *
   * \param[in] swap_bytes
   * Specifies whether the bytes in the destination buffer \p buf are to be
   * swapped after conversion.
   *
   * \returns
   * The number bytes written to the destination buffer, i.e.
   * <c>sizeof(X_)*count</c>.
   */
  template<typename X_, typename T_>
  static std::size_t xencode(void* buf, const T_* src, const std::size_t count, bool swap_bytes)
  {
    X_* x = reinterpret_cast<X_*>(buf);
    if(swap_bytes)
    {
      for(std::size_t i(0); i < count; ++i)
      {
        x[i] = xswap(static_cast<X_>(src[i]));
      }
    }
    else
    {
      for(std::size_t i(0); i < count; ++i)
      {
        x[i] = static_cast<X_>(src[i]);
      }
    }
    return count * sizeof(X_);
  }

  /**
   * \brief Encodes an array by converting each of its elements to a desired type
   *
   * Effectively, this function performs <c>dest[i] = T_(((X_*)buf)[i])</c> for
   * all i < count.
   *
   * \tparam X_
   * The type of the encoded source array buffer \p buf whose elements are to be
   * converted
   *
   * \tparam T_
   * The type of the destination array \p dest
   *
   * \param[out] dest
   * The destination array that receives the converted array elements
   *
   * \param[in] buf
   * The source buffer whose elements of type \p X_ are to be converted.
   * This buffer is assumed to be at least <c>sizeof(X_)*count</c> bytes in size.
   *
   * \param[in] count
   * The number of elements in the destination array \p dest.
   *
   * \param[in] swap_bytes
   * Specifies whether the bytes in the source buffer \p buf are to be swapped
   * before conversion.
   *
   * \returns
   * The number bytes read from the source buffer, i.e. <c>sizeof(X_)*count</c>.
   */
  template<typename X_, typename T_>
  static std::size_t xdecode(T_* dest, const void* buf, const std::size_t count, bool swap_bytes)
  {
    const X_* x = reinterpret_cast<const X_*>(buf);
    if(swap_bytes)
    {
      for(std::size_t i(0); i < count; ++i)
      {
        dest[i] = static_cast<T_>(xswap(x[i]));
      }
    }
    else
    {
      for(std::size_t i(0); i < count; ++i)
      {
        dest[i] = static_cast<T_>(x[i]);
      }
    }
    return count * sizeof(X_);
  }

  template<typename Tclass_, bool signed_>
  struct TypeHelper
  {
    template<typename T_>
    static Pack::Type deduct()
    {
      return Pack::Type::None;
    }
  };

  // specialization for floating point types
  template<>
  struct TypeHelper<FEAT::Type::FloatingClass, true>
  {
    template<typename T_>
    static Pack::Type deduct()
    {
      switch(sizeof(T_))
      {
      // case 1: return Pack::Type::F8;
      case 2: return Pack::Type::F16;
      case 4: return Pack::Type::F32;
      case 8: return Pack::Type::F64;
      case 16: return Pack::Type::F128;
      default: return Pack::Type::None;
      }
    }

    template<typename T_>
    static std::size_t
    encode(void* buf, const T_* t, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
    {
      switch(pack_type)
      {
#ifdef FEAT_HAVE_PACK_TYPE_F16
      case Pack::Type::F16: return xencode<Pack::f16>(buf, t, count, swap_bytes);
#endif
      case Pack::Type::F32: return xencode<Pack::f32>(buf, t, count, swap_bytes);
      case Pack::Type::F64: return xencode<Pack::f64>(buf, t, count, swap_bytes);
#ifdef FEAT_HAVE_PACK_TYPE_F128
      case Pack::Type::F128: return xencode<Pack::f128>(buf, t, count, swap_bytes);
#endif
      default: XABORTM("invalid data type conversion");
      }
    }

    template<typename T_>
    static std::size_t
    decode(T_* t, const void* buf, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
    {
      switch(pack_type)
      {
#ifdef FEAT_HAVE_PACK_TYPE_F16
      case Pack::Type::F16: return xdecode<Pack::f16>(t, buf, count, swap_bytes);
#endif
      case Pack::Type::F32: return xdecode<Pack::f32>(t, buf, count, swap_bytes);
      case Pack::Type::F64: return xdecode<Pack::f64>(t, buf, count, swap_bytes);
#ifdef FEAT_HAVE_PACK_TYPE_F128
      case Pack::Type::F128: return xdecode<Pack::f128>(t, buf, count, swap_bytes);
#endif
      default: XABORTM("invalid data type conversion");
      }
    }
  };

  // specialization for signed integer types
  template<>
  struct TypeHelper<FEAT::Type::IntegralClass, true>
  {
    template<typename T_>
    static Pack::Type deduct()
    {
      switch(sizeof(T_))
      {
      case 1: return Pack::Type::I8;
      case 2: return Pack::Type::I16;
      case 4: return Pack::Type::I32;
      case 8: return Pack::Type::I64;
      default: return Pack::Type::None;
      }
    }

    template<typename T_>
    static std::size_t
    encode(void* buf, const T_* t, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
    {
      switch(pack_type)
      {
      case Pack::Type::I8: return xencode<Pack::i8>(buf, t, count, swap_bytes);
      case Pack::Type::I16: return xencode<Pack::i16>(buf, t, count, swap_bytes);
      case Pack::Type::I32: return xencode<Pack::i32>(buf, t, count, swap_bytes);
      case Pack::Type::I64: return xencode<Pack::i64>(buf, t, count, swap_bytes);
      default: XABORTM("invalid data type conversion");
      }
    }

    template<typename T_>
    static std::size_t
    decode(T_* t, const void* buf, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
    {
      switch(pack_type)
      {
      case Pack::Type::I8: return xdecode<Pack::i8>(t, buf, count, swap_bytes);
      case Pack::Type::I16: return xdecode<Pack::i16>(t, buf, count, swap_bytes);
      case Pack::Type::I32: return xdecode<Pack::i32>(t, buf, count, swap_bytes);
      case Pack::Type::I64: return xdecode<Pack::i64>(t, buf, count, swap_bytes);
      default: XABORTM("invalid data type conversion");
      }
    }
  };

  // specialization for unsigned integer types
  template<>
  struct TypeHelper<FEAT::Type::IntegralClass, false>
  {
    template<typename T_>
    static Pack::Type deduct()
    {
      switch(sizeof(T_))
      {
      case 1: return Pack::Type::U8;
      case 2: return Pack::Type::U16;
      case 4: return Pack::Type::U32;
      case 8: return Pack::Type::U64;
      default: return Pack::Type::None;
      }
    }

    template<typename T_>
    static std::size_t
    encode(void* buf, const T_* t, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
    {
      switch(pack_type)
      {
      case Pack::Type::U8: return xencode<Pack::u8>(buf, t, count, swap_bytes);
      case Pack::Type::U16: return xencode<Pack::u16>(buf, t, count, swap_bytes);
      case Pack::Type::U32: return xencode<Pack::u32>(buf, t, count, swap_bytes);
      case Pack::Type::U64: return xencode<Pack::u64>(buf, t, count, swap_bytes);
      default: XABORTM("invalid data type conversion");
      }
    }

    template<typename T_>
    static std::size_t
    decode(T_* t, const void* buf, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
    {
      switch(pack_type)
      {
      case Pack::Type::U8: return xdecode<Pack::u8>(t, buf, count, swap_bytes);
      case Pack::Type::U16: return xdecode<Pack::u16>(t, buf, count, swap_bytes);
      case Pack::Type::U32: return xdecode<Pack::u32>(t, buf, count, swap_bytes);
      case Pack::Type::U64: return xdecode<Pack::u64>(t, buf, count, swap_bytes);
      default: XABORTM("invalid data type conversion");
      }
    }
  };

  template<typename T_>
  Pack::Type deduct_type()
  {
    using TypeTraits = typename FEAT::Type::Traits<T_>;
    return TypeHelper<typename TypeTraits::TypeClass, TypeTraits::is_signed>::template deduct<T_>();
  }

  /**
   * \brief Encodes an array into a packed buffer without compression support.
   *
   * \attention
   * This function does not encode compressed pack types.
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
  static std::size_t
  encode_raw(void* buf, const T_* src, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
  {
    if(count <= std::size_t(0))
    {
      return std::size_t(0);
    }

    XASSERT(buf != nullptr);
    XASSERT(src != nullptr);

    XASSERTM((pack_type & Pack::Type::Mask_T) != Pack::Type::None, "invalid pack type");
    XASSERTM((pack_type & Pack::Type::Mask_Z) != Pack::Type::Mask_Z, "cannot encode compressed type");

    using TypeTraits = typename FEAT::Type::Traits<T_>;

    return TypeHelper<typename TypeTraits::TypeClass, TypeTraits::is_signed>::encode(
      buf,
      src,
      count,
      pack_type,
      swap_bytes);
  }

  /**
   * \brief Decodes an array from a packed buffer without compression support.
   *
   * \attention
   * This function does not decode compressed pack types.
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
  static std::size_t
  decode_raw(T_* dst, const void* buf, const std::size_t count, const Pack::Type pack_type, bool swap_bytes)
  {
    if(count <= std::size_t(0))
    {
      return std::size_t(0);
    }

    XASSERT(dst != nullptr);
    XASSERT(buf != nullptr);

    XASSERTM((pack_type & Pack::Type::Mask_T) != Pack::Type::None, "invalid pack type");
    XASSERTM((pack_type & Pack::Type::Mask_Z) != Pack::Type::Mask_Z, "cannot decode compressed type");

    using TypeTraits = typename FEAT::Type::Traits<T_>;

    return TypeHelper<typename TypeTraits::TypeClass, TypeTraits::is_signed>::decode(
      dst,
      buf,
      count,
      pack_type,
      swap_bytes);
  }

  /**
   * \brief Computes the estimated (upper bound) pack buffer size for an array for
   * lossy compression
   *
   * \warning
   * Lossless estimation and lossy estimation of buffer size are different!
   *
   * \param[in] src
   * A pointer to the input array that is to be packed. Must not be \c nullptr.
   *
   * \param[in] count
   * The number of input array elements.
   *
   * \param[in] tolerance
   * The desired maximum error of compressed data.
   *
   * \returns
   * An estimated (upper bound) array buffer size.
   *
   */
  static std::size_t lossy_estimate_size(const std::size_t count, const Pack::Type type, const double tolerance)
  {
    if(count <= std::size_t(0))
    {
      return std::size_t(0);
    }
#ifdef FEAT_HAVE_ZFP
    // get type T_, int would be handleable... but why lossy save an int array?:
    zfp_type t = zfp_type::zfp_type_none;
    switch(type)
    {
    case Pack::Type::PF32: t = zfp_type::zfp_type_float; break;
    case Pack::Type::PF64: t = zfp_type::zfp_type_double; break;
    default: XABORTM("cannot encode compressed type; data typ not available");
    }
    XASSERT(t != zfp_type::zfp_type_none);

    // declare needed zfp internal variables:
    zfp_field* field;
    zfp_stream* zfp;
    std::size_t bytes;
    std::size_t x_size, y_size;
    void* temp_arr = nullptr;

    // calculate needed representation size, e.g. 4 times (size/4 +1):
    x_size = 4u;
    y_size = count / x_size + 1;

    // check if we lose significant bits through type conversion...
    if((sizeof(std::size_t) > sizeof(uint)) && (count > std::numeric_limits<uint>::max()))
      XABORTM("cannot encode compressed type; array size to big for internal data");

    // initialize zfp field as 2D variant, that src is not right type does not
    // matter at this point...( memory has to be freed in the end):
    field = zfp_field_2d(temp_arr, t, (uint)x_size, (uint)y_size);
    // open stream(has to be freed)
    zfp = zfp_stream_open(NULL);
    // set error tolerance of stream:
    zfp_stream_set_accuracy(zfp, tolerance);
    // get a conservativ estimate for the buffer size(header size included)
    bytes = zfp_stream_maximum_size(zfp, field);

    // free allocated data:
    zfp_field_free(field);
    zfp_stream_close(zfp);

    // return compression size
    return bytes;

#else  // no FEAT_HAVE_ZFP
    (void)type;
    (void)tolerance;
    XABORTM("cannot encode compressed type; zfp not available");
    return std::size_t(0);
#endif // FEAT_HAVE_ZFP*/
  }

  std::size_t estimate_size(const std::size_t count, const Pack::Type type, const double tolerance)
  {
    if((type & Pack::Type::Mask_P) != Pack::Type::None)
    {
      return lossy_estimate_size(count, type, tolerance);
    }
    // compute raw element size
    std::size_t raw_size = count * Pack::element_size(type & Pack::Type::Mask_T);

#ifdef FEAT_HAVE_ZLIB
    if((type & Pack::Type::Mask_Z) != Pack::Type::None)
    {
      return std::size_t(::compressBound(uLong(raw_size)));
    }
#else  // no FEAT_HAVE_ZLIB
    XASSERTM((type & Pack::Type::Mask_Z) == Pack::Type::None, "cannot estimate compressed size; zlib not available");
#endif // FEAT_HAVE_ZLIB

    // return raw size
    return raw_size;
  }

  template<typename T_>
  std::size_t encode(
    void* buf,
    const T_* src,
    const std::size_t buf_size,
    const std::size_t count,
    const Pack::Type pack_type,
    bool swap_bytes,
    double tolerance)
  {
    if(count <= std::size_t(0))
    {
      return std::size_t(0);
    }

    XASSERT(buf != nullptr);
    XASSERT(src != nullptr);

    // get raw type
    const Pack::Type raw_type = pack_type & Pack::Type::Mask_T;

    // no compression required?
    if((pack_type & (Pack::Type::Mask_Z | Pack::Type::Mask_P)) == Pack::Type::None)
    {
      return encode_raw(buf, src, count, raw_type, swap_bytes);
    }

    // zlib compression?
    if((pack_type & Pack::Type::Mask_Z) == Pack::Type::Mask_Z)
    {
      return lossless_encode(buf, src, buf_size, count, pack_type, swap_bytes);
    }

    // zfp compression?
    if((pack_type & Pack::Type::Mask_P) == Pack::Type::Mask_P)
    {
      return lossy_encode(buf, src, buf_size, count, pack_type, swap_bytes, tolerance);
    }

    return std::size_t(0);
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
   * \returns
   * The total number of bytes written into the output buffer.
   */
  template<typename T_>
  static std::size_t lossless_encode(
    void* buf,
    const T_* src,
    const std::size_t buf_size,
    const std::size_t count,
    const Pack::Type pack_type,
    bool swap_bytes)
  {
    // get raw type
    const Pack::Type raw_type = pack_type & Pack::Type::Mask_T;

#ifdef FEAT_HAVE_ZLIB
    // ensure that the temporary buffer is large enough
    std::size_t raw_size = element_size(raw_type) * count;
    std::vector<char> tmp(raw_size);

    // encode into temporary buffer
    encode_raw(tmp.data(), src, count, raw_type, swap_bytes);

    // cast buffer lengths to zlib types
    uLongf dl = static_cast<uLongf>(buf_size);
    uLong sl = static_cast<uLong>(raw_size);

    // cast buffer pointers to zlib types
    Bytef* dbuf = reinterpret_cast<Bytef*>(buf);
    const Bytef* sbuf = reinterpret_cast<const Bytef*>(tmp.data());

    // compress via zlib
    if(::compress(dbuf, &dl, sbuf, sl) != Z_OK)
      XABORTM("zlib compression error");

    // return final destination buffer usage as reported by zlib
    return static_cast<std::size_t>(dl);
#else  // no FEAT_HAVE_ZLIB
    (void)buf_size;
    (void)buf;
    (void)src;
    (void)raw_type;
    (void)swap_bytes;
    (void)pack_type;
    (void)count;
    XABORTM("cannot encode compressed type; zlib not available");
    return std::size_t(0);
#endif // FEAT_HAVE_ZLIB
  }

  template<typename T_>
  std::size_t decode(
    T_* dst,
    void* buf,
    const std::size_t count,
    const std::size_t buf_size,
    const Pack::Type pack_type,
    bool swap_bytes)
  {
    if(count <= std::size_t(0))
    {
      return std::size_t(0);
    }

    XASSERT(dst != nullptr);
    XASSERT(buf != nullptr);

    // get raw type
    const Pack::Type raw_type = pack_type & Pack::Type::Mask_T;

    // no compression required?
    if((pack_type & (Pack::Type::Mask_Z | Pack::Type::Mask_P)) == Pack::Type::None)
    {
      return decode_raw(dst, buf, count, raw_type, swap_bytes);
    }

    // zlib compression?
    if((pack_type & Pack::Type::Mask_Z) == Pack::Type::Mask_Z)
    {
      return lossless_decode(dst, buf, count, buf_size, pack_type, swap_bytes);
    }

    // zfp compression?
    if((pack_type & Pack::Type::Mask_P) == Pack::Type::Mask_P)
    {
      return lossy_decode(dst, buf, count, buf_size, pack_type, swap_bytes);
    }

    return std::size_t(0);
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
  static std::size_t lossless_decode(
    T_* dst,
    const void* buf,
    const std::size_t count,
    const std::size_t buf_size,
    const Pack::Type pack_type,
    bool swap_bytes)
  {
    // get raw type
    const Pack::Type raw_type = pack_type & Pack::Type::Mask_T;

#ifdef FEAT_HAVE_ZLIB
    // ensure that the temporary buffer is large enough
    std::size_t raw_size = element_size(raw_type) * count + 1024u; // + 1KB
    std::vector<char> tmp(raw_size);

    // cast buffer lengths to zlib types
    uLongf dl = static_cast<uLongf>(raw_size);
    uLong sl = static_cast<uLong>(buf_size);

    // cast buffer pointers to zlib types
    Bytef* dbuf = reinterpret_cast<Bytef*>(tmp.data());
    const Bytef* sbuf = reinterpret_cast<const Bytef*>(buf);

    // decompress via zlib
    if(::uncompress2(dbuf, &dl, sbuf, &sl) != Z_OK)
      XABORTM("zlib decompression error");

    // decode from temporary buffer
    decode_raw(dst, tmp.data(), count, raw_type, swap_bytes);

    // return final source buffer usage as reported by zlib
    return static_cast<std::size_t>(sl);
#else  // no FEAT_HAVE_ZLIB
    (void)buf_size;
    (void)buf;
    (void)dst;
    (void)raw_type;
    (void)swap_bytes;
    (void)pack_type;
    (void)count;
    XABORTM("cannot decode compressed type; zlib not available");
    return std::size_t(0);
#endif // FEAT_HAVE_ZLIB
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
   * The desired maximum error of compressed data.
   *
   * \warning
   * Tolerance can be exceeded in fringe cases. If data set is highly
   * uncontinuous, error testing is recommended.
   *
   * \returns
   * The total number of bytes written into the output buffer.
   *
   */
  template<typename T_>
  static std::size_t lossy_encode(
    void* buf,
    T_* src,
    const std::size_t buf_size,
    const std::size_t count,
    const Pack::Type pack_type,
    bool swap_bytes,
    const double tolerance)
  {
    // get raw type
    const Pack::Type raw_type = pack_type & Pack::Type::Mask_T;
    // Test if zfp is loaded:
#ifdef FEAT_HAVE_ZFP
    // get type T_ :
    zfp_type t = zfp_type::zfp_type_none;
    switch(pack_type)
    {
    case Pack::Type::PF32: t = zfp_type::zfp_type_float; break;
    case Pack::Type::PF64: t = zfp_type::zfp_type_double; break;
    default: XABORTM("cannot encode compressed type; data typ not available");
    }
    XASSERT(t != zfp_type::zfp_type_none);
    // check if header can hold size of the array, in case it cant (more than 1
    // petabyte of double data) give out error for now limitation is only for the
    // header, not for the underlying zfp_array, so there could be a workaround if
    // needed...
    XASSERT(count <= u64(2e14));

    // cast src to the needed type:
    //  ensure that the temporary buffer is large enough
    std::size_t raw_size = element_size(raw_type) * count + 3 * element_size(raw_type);
    std::vector<char> tmp(raw_size);

    // encode into temporary buffer
    encode_raw(tmp.data(), src, count, raw_type, swap_bytes);

    // declare needed zfp internal variables:
    zfp_field* field;
    zfp_stream* zfp;
    std::size_t real_bytes;
    bitstream* stream;
    std::size_t x_size, y_size;

    // check if we lose significant bits through type conversion...
    if((sizeof(std::size_t) > sizeof(uint)) && (count > std::numeric_limits<uint>::max()))
      XABORTM("cannot encode compressed type; array size to big for internal "
              "data structure");
    // calculate needed representation size and rest:
    x_size = 4u;
    if(count % 4u == 0)
      y_size = count / x_size;
    else
      y_size = count / x_size + 1;

    // initialize zfp structures:
    field = zfp_field_2d(tmp.data(), t, (uint)x_size, (uint)y_size);
    zfp = zfp_stream_open(NULL);
    zfp_stream_set_accuracy(zfp, tolerance);
    stream = stream_open(buf, buf_size);
    zfp_stream_set_bit_stream(zfp, stream);
    stream_rewind(stream);
    // write header to stream
    zfp_write_header(zfp, field, zfp_header_mask);
    // compress
    real_bytes = zfp_compress(zfp, field);

    // free allocated data:
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);

    return real_bytes;
#else
    (void)buf_size;
    (void)buf;
    (void)src;
    (void)raw_type;
    (void)swap_bytes;
    (void)pack_type;
    (void)count;
    (void)tolerance;
    XABORTM("cannot encode compressed type; zfp not available");
    return std::size_t(0);
#endif // FEAT_HAVE_ZFP
  }

  /**
   * \brief Decodes an array from a packed buffer with lossy saved compression
   * data.
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
  static std::size_t lossy_decode(
    T_* dst,
    void* buf,
    const std::size_t count,
    const std::size_t buf_size,
    const Pack::Type pack_type,
    bool swap_bytes)
  {
    // get raw type
    const Pack::Type raw_type = pack_type & Pack::Type::Mask_T;
#ifdef FEAT_HAVE_ZFP
    // get type T_ :
    zfp_type t = zfp_type::zfp_type_none;
    switch(pack_type)
    {
    case Pack::Type::PF32: t = zfp_type::zfp_type_float; break;
    case Pack::Type::PF64: t = zfp_type::zfp_type_double; break;
    default: XABORTM("cannot encode compressed type; data typ not available");
    }
    XASSERT(t != zfp_type::zfp_type_none);

    // ensure that the temporary buffer is large enough
    std::size_t raw_size = element_size(raw_type) * count + 1024u; // + 1KB
    std::vector<char> tmp(raw_size);

    // declare intern zfp variables
    zfp_field* field;
    zfp_stream* zfp;
    bitstream* stream;
    std::size_t array_size;
    std::size_t real_bytes;
    // allocating memory
    field = zfp_field_alloc();
    zfp = zfp_stream_open(NULL);
    // aligning stream with buf:
    stream = stream_open(buf, buf_size);
    // read in Codec:
    stream_rseek(stream, 24u);
    uint codec = (uint)stream_read_bits(stream, 8u);
    XASSERTM(
      codec == zfp_codec_version,
      "cannot decode compressed type; zfp version not compatible \n "
      "zfp_Systemcodec: " +
        stringify(zfp_codec_version) + "\n zfp_Compressed codec: " + stringify(codec));
    stream_rewind(stream);
    // aligning zfp_stream with bitsream
    zfp_stream_set_bit_stream(zfp, stream);
    // read header
    zfp_read_header(zfp, field, zfp_header_mask);
    array_size = zfp_field_size(field, NULL);
    // check if T_ and type to be decompressed is same:
    if(t != field->type)
      XABORTM("cannot decode compressed type; given array and saved data do not "
              "have the same type");
    // checking if buffersize of given array is enough to take up saved array....
    if(array_size - 3 > count)
      XABORTM("cannot decode compressed data; given count is too small for "
              "decompressed data!");
    // alligning field and tmp
    field->data = tmp.data();
    // decompress data
    real_bytes = zfp_decompress(zfp, field);

    // decode from temporary buffer
    decode_raw(dst, tmp.data(), count, raw_type, swap_bytes);

    // free memory
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);

    return real_bytes;
#else
    (void)buf_size;
    (void)buf;
    (void)dst;
    (void)count;
    (void)raw_type;
    (void)swap_bytes;
    (void)pack_type;
    XABORTM("cannot decode compressed type; zfp not available");
    return std::size_t(0);
#endif // FEAT_HAVE_ZFP
  }

  /////////////////////////////////////////////////////////////
  // Explicit template instantiations for supported data types
  /////////////////////////////////////////////////////////////

  template Pack::Type deduct_type<i8>();
  template Pack::Type deduct_type<i16>();
  template Pack::Type deduct_type<i32>();
  template Pack::Type deduct_type<i64>();

  template Pack::Type deduct_type<u8>();
  template Pack::Type deduct_type<u16>();
  template Pack::Type deduct_type<u32>();
  template Pack::Type deduct_type<u64>();

#if defined(FEAT_HAVE_PACK_TYPE_F16)
  template Pack::Type deduct_type<f16>();
#endif
  template Pack::Type deduct_type<f32>();
  template Pack::Type deduct_type<f64>();

#if defined(FEAT_HAVE_PACK_TYPE_F128)
  template Pack::Type deduct_type<f128>();
#endif

  template std::size_t encode<i8>(void*, const i8*, std::size_t, std::size_t, Pack::Type, bool, double);
  template std::size_t encode<i16>(void*, const i16*, std::size_t, std::size_t, Pack::Type, bool, double);
  template std::size_t encode<i32>(void*, const i32*, std::size_t, std::size_t, Pack::Type, bool, double);
  template std::size_t encode<i64>(void*, const i64*, std::size_t, std::size_t, Pack::Type, bool, double);

  template std::size_t encode<u8>(void*, const u8*, std::size_t, std::size_t, Pack::Type, bool, double);
  template std::size_t encode<u16>(void*, const u16*, std::size_t, std::size_t, Pack::Type, bool, double);
  template std::size_t encode<u32>(void*, const u32*, std::size_t, std::size_t, Pack::Type, bool, double);
  template std::size_t encode<u64>(void*, const u64*, std::size_t, std::size_t, Pack::Type, bool, double);

  template std::size_t encode<long long>(void*, const long long*, std::size_t, std::size_t, Pack::Type, bool, double);

#if defined(FEAT_HAVE_PACK_TYPE_F16)
  template std::size_t encode<f16>(void*, const f16*, std::size_t, std::size_t, Pack::Type, bool, double);
#endif
  template std::size_t encode<f32>(void*, const f32*, std::size_t, std::size_t, Pack::Type, bool, double);
  template std::size_t encode<f64>(void*, const f64*, std::size_t, std::size_t, Pack::Type, bool, double);

#if defined(FEAT_HAVE_PACK_TYPE_F128)
  template std::size_t encode<f128>(void*, const f128*, std::size_t, std::size_t, Pack::Type, bool, double);
#endif

  template std::size_t decode<i8>(i8*, void*, std::size_t, std::size_t, Pack::Type, bool);
  template std::size_t decode<i16>(i16*, void*, std::size_t, std::size_t, Pack::Type, bool);
  template std::size_t decode<i32>(i32*, void*, std::size_t, std::size_t, Pack::Type, bool);
  template std::size_t decode<i64>(i64*, void*, std::size_t, std::size_t, Pack::Type, bool);

  template std::size_t decode<u8>(u8*, void*, std::size_t, std::size_t, Pack::Type, bool);
  template std::size_t decode<u16>(u16*, void*, std::size_t, std::size_t, Pack::Type, bool);
  template std::size_t decode<u32>(u32*, void*, std::size_t, std::size_t, Pack::Type, bool);
  template std::size_t decode<u64>(u64*, void*, std::size_t, std::size_t, Pack::Type, bool);

  template std::size_t decode<long long>(long long*, void*, std::size_t, std::size_t, Pack::Type, bool);

#if defined(FEAT_HAVE_PACK_TYPE_F16)
  template std::size_t decode<f16>(f16*, void*, std::size_t, std::size_t, Pack::Type, bool);
#endif
  template std::size_t decode<f32>(f32*, void*, std::size_t, std::size_t, Pack::Type, bool);
  template std::size_t decode<f64>(f64*, void*, std::size_t, std::size_t, Pack::Type, bool);

#if defined(FEAT_HAVE_PACK_TYPE_F128)
  template std::size_t decode<f128>(f128*, void*, std::size_t, std::size_t, Pack::Type, bool);
#endif
} // namespace FEAT::Pack
