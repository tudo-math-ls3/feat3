#include <kernel/util/random.hpp>
#include <test_system/test_system.hpp>
#include <cstdint>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the Random class.
 *
 * \test Tests the Random class.
 *
 * \author Peter Zajac
 */
class RandomTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  RandomTest() :
    TaggedTest<Archs::None, Archs::None>("RandomTest")
  {
  }

  virtual ~RandomTest()
  {
  }

  template<typename T_>
  void testi(Random& rng, T_ a, T_ b) const
  {
    // get an arbitrary int
    T_ x;
    rng >> x;

    // get a ranged int
    x = rng(a, b);
    TEST_CHECK((a <= x) && (x <= b));
  }

  template<typename T_>
  void testf(Random& rng, T_ a, T_ b) const
  {
    // get an arbitrary float in range [0,1]
    T_ x;
    rng >> x;
    TEST_CHECK((T_(0) <= x) && (x <= T_(1)));

    // get a ranged float
    x = rng(a, b);
    TEST_CHECK((a <= x) && (x <= b));
  }

  void testb(Random& rng) const
  {
    // get a bool
    bool b;
    rng >> b;

    // get a 'true'
    b = rng(true, true);
    TEST_CHECK(b == true);

    // get a 'false'
    b = rng(false,false);
    TEST_CHECK(b == false);
  }

  virtual void run() const override
  {
    // create an RNG
    Random rng;

    // test 10 times
    for(int i(0); i < 10; ++i)
    {
      // test signed int generations
      testi(rng, std::int8_t(-71), std::int8_t(39));
      testi(rng, std::int16_t(-11549), std::int16_t(17537));
      testi(rng, std::int32_t(-715823882), std::int32_t(1431651765));
      testi(rng, std::int64_t(-715823882), std::int64_t(1431655765));

      // test unsigned int generations
      testi(rng, std::uint8_t(20u), std::uint8_t(190u));
      testi(rng, std::uint16_t(300u), std::uint16_t(57000u));
      testi(rng, std::uint32_t(5000u), std::uint32_t(3894167296u));
      testi(rng, std::uint64_t(70000u), std::uint64_t(13446744073709551616ull));

      // test float extractions
      testf(rng, 0.707f, 3.14f); // float
      testf(rng, 0.707107, 3.14159265); // double
      testf(rng, 0.707107l, 3.14159265l); // long double

      // test bool extraction
      testb(rng);
    }

    // create time seeded RNG
    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rngt(seed);

    // test 10 times
    for(int i(0); i < 10; ++i)
    {
      // test signed int generations
      testi(rngt, std::int8_t(-71), std::int8_t(39));
      testi(rngt, std::int16_t(-11549), std::int16_t(17537));
      testi(rngt, std::int32_t(-715823882), std::int32_t(1431651765));
      testi(rngt, std::int64_t(-715823882), std::int64_t(1431655765));

      // test unsigned int generations
      testi(rngt, std::uint8_t(20u), std::uint8_t(190u));
      testi(rngt, std::uint16_t(300u), std::uint16_t(57000u));
      testi(rngt, std::uint32_t(5000u), std::uint32_t(3894167296u));
      testi(rngt, std::uint64_t(70000u), std::uint64_t(13446744073709551616ull));

      // test float extractions
      testf(rngt, 0.707f, 3.14f); // float
      testf(rngt, 0.707107, 3.14159265); // double
      testf(rngt, 0.707107l, 3.14159265l); // long double

      // test bool extraction
      testb(rngt);
    }
  }
} random_test;
