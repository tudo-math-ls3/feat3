#pragma once
#ifndef KERNEL_UTIL_RANDOM_HPP
#define KERNEL_UTIL_RANDOM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#include <stdint.h>

namespace FEAST
{
  /**
   * \brief Pseudo-Random Number Generator
   *
   * This class implements a 32-bit XOR-Shift pseudo-random number generator with a 128-bit cycle.
   *
   * The intention of this class is to provide a pseudo-random number generator which generates (warning:
   * paradoxon ahead) a reproducible random number stream independent of the platform and the compiler in use.
   *
   * \warning This class is (intentionally) not thread-safe -- in a multi-threaded application each thread should
   * use its own random number generator.
   *
   * \see <b>G. Marsaglia:</b> <em>Xorshift RNGs</em>; Journal of Statistical Software, Vol. 8 (Issue 14), 2003\n
   *      http://www.jstatsoft.org/v08/i14
   * \see http://en.wikipedia.org/wiki/Xorshift
   *
   * \note The implementation in this class corresponds (with the exception of the default seed value) to the
   * <c>xor128()</c> example in Section 4 of the <em>Xorshift RNGs</em> article mentioned above.
   *
   * \author Peter Zajac
   */
  class Random
  {
  private:
    /// \cond internal
    uint32_t _s, _x, _y, _z;

    // Returns the next number in the stream.
    inline uint32_t _get()
    {
      uint32_t t = _s ^ (_s << 11);
      _s = _x;
      _x = _y;
      _y = _z;
      return _z = _z ^ (_z >> 19) ^ t ^ (t >> 8);
    }
    /// \endcond

  public:
    /**
     * \brief CTOR
     *
     * \param[in] seed
     * The seed for the random number generator.
     */
    explicit Random(uint32_t seed = 428147976u) :
      _s(seed),
      _x(362436069u),
      _y(521288629u),
      _z(88675123u)
    {
    }

    /**
     * \brief Integer extraction operator
     *
     * This function extracts the next random integer in the stream.
     *
     * \param[out] x
     * The value to be extracted.
     *
     * \returns \p *this
     */
    Random& operator>>(uint8_t& x)
    {
      x = uint8_t(_get() & 0xFFu);
      return *this;
    }

    /** \copydoc operator>>(uint8_t&) */
    Random& operator>>(uint16_t& x)
    {
      x = uint16_t(_get() & 0xFFFFu);
      return *this;
    }

    /** \copydoc operator>>(uint8_t&) */
    Random& operator>>(uint32_t& x)
    {
      x = _get();
      return *this;
    }

    /** \copydoc operator>>(uint8_t&) */
    Random& operator>>(uint64_t& x)
    {
      uint32_t* t = reinterpret_cast<uint32_t*>(&x);
      t[0] = _get();
      t[1] = _get();
      return *this;
    }

    /** \copydoc operator>>(uint8_t&) */
    Random& operator>>(int8_t& x)
    {
      return this->operator>>(*reinterpret_cast<uint8_t*>(&x));
    }

    /** \copydoc operator>>(uint8_t&) */
    Random& operator>>(int16_t& x)
    {
      return this->operator>>(*reinterpret_cast<uint16_t*>(&x));
    }

    /** \copydoc operator>>(uint8_t&) */
    Random& operator>>(int32_t& x)
    {
      return this->operator>>(*reinterpret_cast<uint32_t*>(&x));
    }

    /** \copydoc operator>>(uint8_t&) */
    Random& operator>>(int64_t& x)
    {
      return this->operator>>(*reinterpret_cast<uint64_t*>(&x));
    }

    /**
     * \brief Floating-Point extraction operator
     *
     * This function extracts the next random floating point number in the stream.
     * The value extracted by this operator will be in range [0,1].
     *
     * \param[out] x
     * The value to be extracted.
     *
     * \returns \p *this
     */
    Random& operator>>(double& x)
    {
      static const double d = 1.0 / double(0xFFFFFFFFul);
      x = double(_get()) * d;
      return *this;
    }

    /** \copydoc operator>>(double&) */
    Random& operator>>(long double& x)
    {
      double d(0.0);
      this->operator>>(d);
      x = long double(d);
      return *this;
    }

    /** \copydoc operator>>(double&) */
    Random& operator>>(float& x)
    {
      double d(0.0);
      this->operator>>(d);
      x = float(d);
      return *this;
    }
  }; // class Random

} // namespace FEAST

#endif // KERNEL_UTIL_RANDOM_HPP
