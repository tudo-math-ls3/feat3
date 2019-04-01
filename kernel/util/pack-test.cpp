// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/pack.hpp>

#include <vector>

using namespace FEAT;
using namespace FEAT::TestSystem;

class PackTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  PackTest() :
    TaggedTest<Archs::None, Archs::None>("PackTest")
  {
  }

  template<typename DT_>
  void test_pack_float(const Pack::Type pack_type, const DT_ tol, bool swap_bytes = false) const
  {
    static constexpr std::size_t n = 37;
    const DT_ scal = Math::pi<DT_>() / DT_(2*n);
    DT_ idata[n], odata[n];

    // fill array with sine values
    for(std::size_t i(0); i < n; ++i)
      idata[i] = Math::sin(DT_(i+1) * scal);

    // estimate buffer size and create buffer vector
    const std::size_t est_size = Pack::estimate_size(n, pack_type);
    std::vector<char> pbuf(est_size + std::size_t(64u));

    // pack array into buffer and validate size
    const std::size_t dst_size = Pack::encode(pbuf.data(), idata, pbuf.size(), n, pack_type, swap_bytes);
    TEST_CHECK(dst_size <= est_size);

    // unpack array from buffer and validate size
    const std::size_t src_size = Pack::decode(odata, pbuf.data(), n, pbuf.size(), pack_type, swap_bytes);
    TEST_CHECK_EQUAL(src_size, dst_size);

    // compare with input array
    for(std::size_t i(0); i < n; ++i)
    {
      DT_ abs_err = Math::abs(odata[i] - idata[i]);
      DT_ rel_err = abs_err / Math::abs(idata[i]);
      TEST_CHECK_IN_RANGE(rel_err, DT_(0), tol);
    }
  }

  template<typename IT_>
  void test_pack_int(const Pack::Type pack_type, bool swap_bytes = false) const
  {
    // Note: do not increase, because (n/2)^2 has to fit into an
    // 8-bit signed integer due to the formula used below!
    static constexpr std::size_t n = 23;
    IT_ idata[n], odata[n];

    // fill array with square values
    for(std::size_t i(0); i < n; ++i)
    {
      IT_ j = IT_(i) - IT_(n) / 2;
      idata[i] = (j < IT_(0) ? -(j*j) : +(j*j));
    }

    // estimate buffer size and create buffer vector
    const std::size_t est_size = Pack::estimate_size(n, pack_type);
    std::vector<char> pbuf(est_size + std::size_t(64u));

    // pack array into buffer and validate size
    const std::size_t dst_size = Pack::encode(pbuf.data(), idata, pbuf.size(), n, pack_type, swap_bytes);
    TEST_CHECK(dst_size <= est_size);

    // unpack array from buffer and validate size
    const std::size_t src_size = Pack::decode(odata, pbuf.data(), n, pbuf.size(), pack_type, swap_bytes);
    TEST_CHECK_EQUAL(src_size, dst_size);

    // compare with input array
    for(std::size_t i(0); i < n; ++i)
    {
      TEST_CHECK_EQUAL(idata[i], odata[i]);
    }
  }

  virtual void run() const override
  {
    // test element sizes
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::F16), 2);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::F32), 4);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::F64), 8);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::F128), 16);

    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::I8), 1);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::I16), 2);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::I32), 4);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::I64), 8);

    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::U8), 1);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::U16), 2);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::U32), 4);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::U64), 8);

    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZF16), 2);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZF32), 4);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZF64), 8);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZF128), 16);

    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZI8), 1);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZI16), 2);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZI32), 4);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZI64), 8);

    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZU8), 1);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZU16), 2);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZU32), 4);
    TEST_CHECK_EQUAL(Pack::element_size(Pack::Type::ZU64), 8);

    // test array packing
    test_pack_float<float> (Pack::Type::F64, 1E-7F); // float -> double
    test_pack_float<double>(Pack::Type::F32, 1E-7);  // double -> float

    test_pack_float<float> (Pack::Type::F32, 1E-7F, true); // float -> swap(float)
    test_pack_float<double>(Pack::Type::F64, 1E-15, true); // double -> swap(double)

#ifdef FEAT_HAVE_PACK_TYPE_F16
    test_pack_float<float> (Pack::Type::F16, 1E-3F); // float -> half
    test_pack_float<double>(Pack::Type::F64, 1E-3);  // double -> half
#endif // FEAT_HAVE_PACK_TYPE_F16
#ifdef FEAT_HAVE_PACK_TYPE_F128
    test_pack_float<double>    (Pack::Type::F128, 1E-15);  // double -> quad
    test_pack_float<__float128>(Pack::Type::F64 , 1E-15Q); // quad -> double
    test_pack_float<__float128>(Pack::Type::F128, 1E-33Q, true); // quad -> swap(quad)
#endif // FEAT_HAVE_PACK_TYPE_F128

    test_pack_int<int>      (Pack::Type::I16); // int -> int16 (shrink)
    test_pack_int<int>      (Pack::Type::I32); // int -> int32
    test_pack_int<int>      (Pack::Type::I64); // int -> int64 (expand)

    test_pack_int<long long>(Pack::Type::I16); // long long -> int16
    test_pack_int<long long>(Pack::Type::I32); // long long -> int32
    test_pack_int<long long>(Pack::Type::I64); // long long -> int64

    test_pack_int<int>      (Pack::Type::I8 , true);  // int -> swap(int8)
    test_pack_int<int>      (Pack::Type::I16, true);  // int -> swap(int16)
    test_pack_int<int>      (Pack::Type::I32, true);  // int -> swap(int32)
    test_pack_int<int>      (Pack::Type::I64, true);  // int -> swap(int64)

#ifdef FEAT_HAVE_ZLIB
    // a subset of the above types should suffice here
    test_pack_float<float> (Pack::Type::ZF32, 1E-7F); // float -> float
    test_pack_float<float> (Pack::Type::ZF64, 1E-7F); // float -> double
    test_pack_float<double>(Pack::Type::ZF32, 1E-7);  // double -> float
    test_pack_float<double>(Pack::Type::ZF64, 1E-15); // double -> double
    test_pack_int<int>(Pack::Type::ZI16, false); // int -> int16 (shrink)
    test_pack_int<int>(Pack::Type::ZI32, true ); // int -> int32 (swap)
    test_pack_int<int>(Pack::Type::ZI64, false); // int -> int64 (expand)
#endif // FEAT_HAVE_ZLIB
  }
} pack_test;
