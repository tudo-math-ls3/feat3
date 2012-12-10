#pragma once
#ifndef KERNEL_UTIL_RANDOM_HPP
#define KERNEL_UTIL_RANDOM_HPP 1

// includes, FEAST
#include <kernel/util/type_traits.hpp>

// includes, system
#include <algorithm>
#include <stdint.h>
#include <limits>

namespace FEAST
{
  /// \cond interal
  namespace Intern
  {
    // Random number class; this class's gen() function takes a Random object reference and
    // uses its next() function to generate a random number of the specified type T_.
    template<typename T_, typename TypeClass_ = typename Type::Traits<T_>::TypeClass>
    class RandomNumber;
  } // namespace Intern
  /// \endcond

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
  public:
    /// seed type
    typedef uint32_t SeedType;

    /// default seed enumeration
    enum
    {
      /// default s-seed value
      def_seed_s = 428147976,
      /// default x-seed value
      def_seed_x = 362436069,
      /// default y-seed value
      def_seed_y = 521288629,
      /// default z-seed value
      def_seed_z = 88675123
    };

  private:
    /// the rng's working values
    uint32_t _s, _x, _y, _z;

  public:
    /**
     * \brief CTOR
     *
     * \param[in] seed
     * The seed for the random number generator.
     */
    explicit Random(SeedType seed = SeedType(def_seed_s)) :
      _s(seed),
      _x(uint32_t(def_seed_x)),
      _y(uint32_t(def_seed_y)),
      _z(uint32_t(def_seed_z))
    {
    }

    /**
     * \brief Returns the next number in the stream.
     *
     * This function returns the next unsigned 32-bit integer in the random number stream and advances the
     * stream.
     */
    uint32_t next()
    {
      uint32_t t = _s ^ (_s << 11);
      _s = _x;
      _x = _y;
      _y = _z;
      return _z = _z ^ (_z >> 19) ^ t ^ (t >> 8);
    }

    /**
     * \brief Extraction operator
     *
     * This function extracts the next random number in the stream.
     * - If the type of the value to be extracted is integral, then the value will be within the full range, i.e.
     *   - for a signed <em>n</em>-bit integer, it will hold <em>-2^(n-1) <= x < 2^(n-1)</em>
     *   - for an unsigned <em>n</em>-bit integer, it will hol that <em>0 <= x < 2^n</em>
     * - If the type of the value to be extracted is a floating point type, then the value will be within the
     *   closed interval [0,1].
     *
     * \tparam T_
     * The type of the value to be extracted. Must be either of integral or floating point type.
     *
     * \param[out] x
     * The reference to a variable receiving the extracted value.
     *
     * \returns \p *this
     */
    template<typename T_>
    Random& operator>>(T_& x)
    {
      x = Intern::RandomNumber<T_>::gen(*this);
      return *this;
    }

    /**
     * \brief Ranged evaluation operator.
     *
     * This operator returns a random number in the range passed to this operator.
     *
     * \tparam T_
     * The type of the value to be extracted. Is determined automatically by \p a and \p b.
     *
     * \param[in] a, b
     * The minimum and maximum values for the range of the random number to be generated.
     * \note If \p b < \p a, then \p a and \p b are swapped.
     * \note If \p a = \p b, then the returned value will be equal to \p a, however, the
     *       generator will still advance its position within the random number stream.
     *
     * \returns
     * A random value in the range [a, b].
     */
    template<typename T_>
    T_ operator()(T_ a, T_ b)
    {
      return Intern::RandomNumber<T_>::gen_ranged(*this, std::min(a, b), std::max(a, b));
    }
  }; // class Random

  /// \cond internal
  namespace Intern
  {
    // Random integer class; this class is the base class for the integer specialisation
    // of the RandomNumber class defined below.
    template<typename T_, size_t num_bytes_ = sizeof(T_)>
    class RandomInteger;

    // specialisation for 8-bit integers
    template<typename T_>
    class RandomInteger<T_, 1>
    {
    public:
      static T_ gen(Random& rng)
      {
        return (T_)(rng.next() & 0xFFu);
      }
    };

    // specialisation for 16-bit integers
    template<typename T_>
    class RandomInteger<T_, 2>
    {
    public:
      static T_ gen(Random& rng)
      {
        return (T_)(rng.next() & 0xFFFFu);
      }
    };

    // specialisation for 32-bit integers
    template<typename T_>
    class RandomInteger<T_, 4>
    {
    public:
      static T_ gen(Random& rng)
      {
        return (T_)(rng.next());
      }
    };

    // specialisation for 64-bit integers
    template<typename T_>
    class RandomInteger<T_, 8>
    {
    public:
      static T_ gen(Random& rng)
      {
        // the rng works with 32-bit ints, so we need to generate a pair of those
        // and concatenate them to obtain a random 64-bit integer...
        union
        {
          uint32_t t[2];
          uint64_t tq;
        } x;
        x.t[0] = rng.next();
        x.t[1] = rng.next();
        return (T_)x.tq;;
      }
    };

    // specialisation for integers
    template<typename T_>
    class RandomNumber<T_, Type::IntegralClass> :
      public RandomInteger<T_>
    {
    public:
      static T_ gen_ranged(Random& rng, T_ a, T_ b)
      {
        // call gen() prior to the if-clause to advance the rng in any case
        T_ x(RandomInteger<T_>::gen(rng));
        if(a < b)
          return a + x % (b - a + 1);
        else
          return a;
      }
    };

    // specialisation for floating point numbers
    template<typename T_>
    class RandomNumber<T_, Type::FloatingClass>
    {
    public:
      static T_ gen(Random& rng)
      {
        // the rng generates numbers in range [0,2^32-1], so divide by that maximum range
        static const double d = 1.0 / double(0xFFFFFFFFu);
        return T_(double(rng.next()) * d);
      }

      static T_ gen_ranged(Random& rng, T_ a, T_ b)
      {
        // machine exactness; we do not want to divide by anything smaller than that...
        static const T_ eps = std::numeric_limits<T_>::epsilon();

        // call gen() prior to the if-clause to advance the rng in any case
        T_ x(gen(rng));
        if(a + eps < b)
          return a + x / (b - a);
        else
          return a;
      }
    };
  } // namespacce Intern
  /// \endcond
} // namespace FEAST

#endif // KERNEL_UTIL_RANDOM_HPP
