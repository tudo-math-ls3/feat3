// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/slotmap.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Util;

struct DestructorCounter
{
  int* counter;

  DestructorCounter() : counter(nullptr)
  {
  }

  DestructorCounter(int* ref) : counter(ref)
  {
  }

  ~DestructorCounter()
  {
    if(counter != nullptr)
    {
      *counter += 1;
    }
  }

  DestructorCounter(const DestructorCounter& other) : counter(other.counter)
  {
  }

  DestructorCounter(DestructorCounter&& other) noexcept : counter(other.counter)
  {
    other.counter = nullptr;
  }

  DestructorCounter& operator=(const DestructorCounter& other)
  {
    if(this == &other)
    {
      return *this;
    }

    counter = other.counter;

    return *this;
  }

  DestructorCounter& operator=(DestructorCounter&& other)
  {
    this->counter = other.counter;
    other.counter = nullptr;
    return *this;
  }
};

struct CustomKey
{
  CustomKey() : index(0), generation(0)
  {
  }
  CustomKey(std::uint32_t i, std::uint32_t g) : index(i), generation(g)
  {
  }

  std::uint32_t index;
  std::uint32_t generation;
};

class SlotMapTest : public UnitTest
{
public:
  SlotMapTest() : UnitTest("SlotMapTest")
  {
  }

  virtual ~SlotMapTest()
  {
  }

  virtual void run() const override
  {
    test_insertion();
    test_deletion();
    test_move_semantics();
    test_capacity();
    test_custom_key();
    test_contains_key();
    test_replace();
    test_secondarymap();
    test_iterators();
  }

  void test_iterators() const
  {
    SlotMap<int> slotmap;
    std::vector<SlotMapKey> keys;

    keys.push_back(slotmap.insert(1));
    keys.push_back(slotmap.insert(3));
    keys.push_back(slotmap.insert(5));
    keys.push_back(slotmap.insert(7));
    keys.push_back(slotmap.insert(9));

    Index count = 0;
    for(auto&& entry : slotmap)
    {
      // Keys point towards correct values
      TEST_CHECK_EQUAL(entry.value, slotmap[entry.key]);

      // Mutate value
      entry.value *= 2;

      count++;
    }

    // Iterator covers correct count of elements
    TEST_CHECK_EQUAL(count, slotmap.size());

    for(auto entry : slotmap)
    {
      // Entries should all be even after mutation in previous loop
      TEST_CHECK(entry.value % 2 == 0);
    }

    // Create some erase-insert churn
    SlotMapKey key = keys[0];
    for(int i = 0; i < 20; i++)
    {
      slotmap.erase(key);
      key = slotmap.insert(i);
    }

    for(auto entry : slotmap)
    {
      // Entries still point towards correct values
      TEST_CHECK_EQUAL(entry.value, slotmap[entry.key]);
    }

    // Deletion during iteration should work just like for std::unordered_maps
    auto it = slotmap.begin();
    while(it != slotmap.end())
    {
      slotmap.erase((*it).key);
    }
    TEST_CHECK_EQUAL(slotmap.size(), 0);
  }

  void test_secondarymap() const
  {
    SlotMap<int> slotmap;
    SecondaryMap<char> secondary;

    SlotMapKey a = slotmap.insert(5);
    SlotMapKey b = slotmap.insert(10);
    secondary.insert(a, 'a');

    TEST_CHECK_EQUAL(secondary[a], 'a');
    TEST_CHECK_EQUAL(secondary.contains_key(b), false);

    a.generation++;
    TEST_CHECK_EQUAL(secondary.contains_key(a), false);
  }

  void test_replace() const
  {
    SlotMap<int> slotmap;
    SlotMapKey a = slotmap.insert(5);
    SlotMapKey b = slotmap.replace(a, 10);
    // New key should lead to new value
    TEST_CHECK_EQUAL(slotmap[b], 10);
    // Old key should be deleted
    TEST_CHECK(!slotmap.contains_key(a));
  }

  void test_contains_key() const
  {
    SlotMap<int> slotmap;
    SlotMapKey a = slotmap.insert(5);
    // SlotMap should contain just inserted key
    TEST_CHECK(slotmap.contains_key(a));
    slotmap.erase(a);
    // SlotMap should not contain erased key
    TEST_CHECK(!slotmap.contains_key(a));
    // SlotMap should not contain a made up key
    TEST_CHECK((!slotmap.contains_key(SlotMapKey(10, 10))));
  }

  void test_custom_key() const
  {
    SlotMap<int, CustomKey> slotmap;

    CustomKey a = slotmap.insert(5);
    TEST_CHECK_EQUAL(slotmap[a], 5);

    slotmap.erase(a);
    TEST_CHECK(!slotmap.contains_key(a));
  }

  void test_capacity() const
  {
    SlotMap<int> slotmap;
    slotmap.reserve(16);
    TEST_CHECK(slotmap.capacity() >= 16);
    slotmap.reserve(64);
    TEST_CHECK(slotmap.capacity() >= 64);
  }

  void test_move_semantics() const
  {
    int counter = 0;
    SlotMap<DestructorCounter> slotmap;
    SlotMapKey keys[10];

    for(int i = 0; i < 10; i++)
    {
      keys[i] = slotmap.insert(DestructorCounter(&counter));
    }

    for(int i = 0; i < 10; i++)
    {
      slotmap.erase(keys[i]);
    }

    TEST_CHECK_EQUAL(counter, 10);
  }

  void test_deletion() const
  {
    SlotMap<int> slotmap;

    for(int i = 0; i < 64; i++)
    {
      slotmap.erase(slotmap.insert(i));
    }
    // Memory should be reused
    TEST_CHECK(slotmap.capacity() == 1);

    SlotMapKey key = slotmap.insert(1);
    slotmap.erase(key);
    // Keys should be invalid after erasure
    TEST_CHECK(!slotmap.contains_key(key));
  }

  void test_insertion() const
  {
    SlotMap<int> slotmap;

    SlotMapKey a = slotmap.insert(5);
    TEST_CHECK_EQUAL(slotmap[a], 5);

    SlotMapKey b = slotmap.insert(10);
    TEST_CHECK_EQUAL(slotmap[b], 10);
  }
} slotmaptest;
