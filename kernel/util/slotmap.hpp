// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, system
#include <algorithm>
#include <optional>
#include <vector>

// includes, feat
#include <kernel/util/assertion.hpp>

namespace FEAT::Util
{
  /**
   * \brief Default key type for SlotMap
   *
   * This is the default key type handed out by the SlotMap implementation.
   * See the SlotMap documentation for more details.
     *
   * \author Markus Muegge
   */
  struct SlotMapKey
  {
    SlotMapKey() : index(0), generation(0)
    {
    }
    SlotMapKey(std::uint32_t i, std::uint32_t g) : index(i), generation(g)
    {
    }

    std::uint32_t index;
    std::uint32_t generation;
  };

  /**
   * \brief Return type for \c SlotMapIterator
   *
   * Contains a reference to data stored in the \c SlotMap and the corresponding key.
   *
   * \author Markus Muegge
   */
  template<typename V, typename Key = SlotMapKey>
  struct KeyValuePair
  {
    KeyValuePair(V& v, Key k) : value(v), key(k)
    {
    }

    V& value;
    Key key;
  };

  /**
   * \brief Return type for \c ConstSlotMapIterator
   *
   * Contains a const reference to data stored in the \c SlotMap and the corresponding key.
   *
   * \author Markus Muegge
   */
  template<typename V, typename Key = SlotMapKey>
  struct ConstKeyValuePair
  {
    ConstKeyValuePair(const V& v, Key k) : value(v), key(k)
    {
    }

    const V& value;
    Key key;
  };

  /**
   * \brief High-performance associative container
   *
   * \tparam V
   * Type of value to be stored.
   * \tparam Key
   * Type of keys returned by this SlotMap. SlotMapKey by default.
   *
   * A SlotMap is a high-performance associative container. It:
   * * allows for O(1) insertion, deletion, lookup
   * * stores its elements densely
   * * avoids the ABA problem, i.e. the user can notice that data stored at a key has changed
   *
   * To do so it hands out its own keys, instead of letting the user provide a key like
   * a traditional hashmap would do. Keys consist of an index into the SlotMap used to
   * access elements, and a generation used to detect if the slot has been reused since
   * the key has been issued.
   *
   * # Implementation
   * The implementation consists of three dynamic arrays (std::vectors) called
   * slots, data, erase, as well as head and tail indices for a list of free slots.
   *
   * The data array is where the inserted elements are stored.
   *
   * The slots array provides the indirection required to keep keys stable while guaranteeing
   * to store the elements densely. Each entry in the slots array consists of an index and a
   * generation. The generation indicates how often this slot has been reused. This is what
   * allows detecting if an old key is used to access new data. The index either points into
   * the data array (if the slot is currently occupied) or to the next free slot (if the slot
   * is currently not occupied). The keys given out by the SlotMap point towards slots. This way
   * elements in the data array can be shifted around and only indices in the slots array need
   * to be fixed, while the keys stay valid.
   *
   * The erase array contains the index of the slot pointing towards itself (and the associated
   * element in the data array). In short slots[erase[i]].index == i. It is used to quickly find
   * affected slots during deletion of a key, allowing O(1) deletions.
   *
   * The freelist tracks the start and end of the linked list of unoccupied slots established.
   * As an invariant the SlotMap always keeps at least one Slot unoccupied
   * to guarantee that freelist.head and freelist.tail point towards a valid slot.
   * Note that keeping track of the tail is not strictly required.
   * A stack of unoccupied slots could be used instead, but using a list
   * increments the slots generations more evenly.
   *
   * Inserting a new value works as follows
   * * Retrieve the next unoccupied slot from the freelist
   * * Advance the freelist head or push a new slot to maintain freelist invariant
   * * Point retrieved slot to new element
   * * Push value into data array
   * * Push slot index into erase array
   * * Return key pointing to retrieved slot to user. Keys generation matches that of the slot.
   *
   * Deleting a value works as follows
   * * Retrieve the slot for the given key, if the key is valid
   * * Increment the slots generation. This invalidates all prior keys.
   * * Delete the data and erase entries the slot is pointing to
   * * Move the back elements of data and erase into the now free positions
   * * Fix the slot pointing to the previously last element using the erase array
   * * Add the now unoccupied slot to the back of the freelist
   *
   * See https://www.youtube.com/watch?v=SHaAR7XPtNU for more details.
   *
   * # Custom Keys
   * The SlotMap returns keys of type SlotMapKey per default. For more type safety
   * when using multiple SlotMaps simultaneously a custom key type can be provided
   * as a second template argument.
   *
   * A valid key type requires a constructor accepting the index and the generation
   * in that order, a public member index and a public member generation.
   *
   * \author Markus Muegge
   */
  template<typename V, typename Key = SlotMapKey>
  class SlotMap
  {
  public:
    using ValueType = V;
    using KeyType = Key;

    /**
     * \brief Construct \c SlotMap with initial capacity zero.
     */
    SlotMap() : freelist{0, 0}
    {
      // Create initial slot, so that freelist is always valid
      _slots.emplace_back(Slot());
    }

    /**
     * \brief Insert \c value into the \c SlotMap
     *
     * \param[in] value
     * The value to be inserted.
     *
     * \returns A key pointing to the stored value.
     */
    Key insert(V&& value)
    {
      std::uint32_t idx = freelist.head;

      // Adjust freelist
      if(idx == freelist.tail)
      {
        // Used up the last Slot. Push a new one to keep freelist valid
        _slots.push_back(Slot());
        freelist.head = freelist.tail = (std::uint32_t)_slots.size() - 1;
      }
      else
      {
        // Move head to next element of freelist
        freelist.head = _slots[idx].index;
      }

      Slot& slot = _slots[idx];

      // Point slot to new data
      slot.index = (std::uint32_t)_data.size();

      // Push new data
      _data.push_back(std::move(value));

      // Update erase list
      // Slot for new data is _slots[idx]
      _erase.push_back(idx);

      return Key(idx, slot.generation);
    }

    /**
     * \brief Insert a value into the \c SlotMap
     *
     * \param[in] value
     * The value to be inserted.
     *
     * \returns A key pointing to the stored value.
     */
    Key insert(const V& value)
    {
      std::uint32_t idx = freelist.head;

      if(idx == freelist.tail)
      {
        // Used up the last Slot. Push a new one to keep freelist valid
        _slots.push_back(Slot());
        freelist.head = freelist.tail = (std::uint32_t)_slots.size() - 1;
      }
      else
      {
        // Move head to next element of freelist
        freelist.head = _slots[idx].index;
      }

      Slot& slot = _slots[idx];

      // Point slot to new data
      slot.index = static_cast<std::uint32_t>(_data.size());

      // Push new data
      _data.push_back(value);

      // Update erase list
      // Slot for new data is _slots[idx]
      _erase.push_back(idx);

      return Key(idx, slot.generation);
    }

    /**
     * \brief Remove a key from the \c SlotMap
     *
     * Will invalidate all further uses of the key.
     *
     * \param[in] key The key to remove
     */
    void erase(const Key& key)
    {
      if(!contains_key(key))
      {
        return;
      }
      // Retrieve Slot for key
      Slot& slot = try_get_slot(key);

      // Increment the slots generation
      // Previous keys are no longer valid
      slot.generation++;

      // Exchange positions of to-be-deleted element and last element
      // of _data and _erase for O(1) deletion later
      std::swap(_data[slot.index], _data[_data.size() - 1]);
      std::swap(_erase[slot.index], _erase[_erase.size() - 1]);

      // Fix Slot of moved element
      _slots[_erase[slot.index]].index = slot.index;

      // Actually delete elements
      _data.pop_back();
      _erase.pop_back();

      // Slot can now be reused. Add slot to freelist
      _slots[freelist.tail].index = key.index;
      freelist.tail = key.index;
      _slots[freelist.tail].index = freelist.tail;
    }

    /**
     * \brief Replace the value stored at \c key with \c value
     *
     * \param[in] key The key to replace
     * \param[in] value The new value to store
     * \returns The key of the new value
     *
     * Shortcut to erasing a key and then inserting a new value.
     * The new value is not guaranteed to be stored at the same memory location
     * as the previous one.
     */
    Key replace(Key& key, V&& value)
    {
      erase(key);
      return insert(value);
    }

    /**
     * \brief Replace the value stored at \c key with \c value
     *
     * \param[in] key The key to replace
     * \param[in] value The new value to store
     * \returns The key of the new value
     *
     * Shortcut to erasing a key and then inserting a new value.
     * The new value is not guaranteed to be stored at the same memory location
     * as the previous one.
     */
    Key replace(Key& key, V& value)
    {
      erase(key);
      return insert(value);
    }

    /**
     * \brief Check if given key is still valid.
     *
     * \param[in] key The key to check.
     * \returns \c true if the key is valid, \c false otherwise.

     * A key is valid if it points to an existing slot and matches that slots generation.
    */
    bool contains_key(const Key& key) const
    {
      return (key.index < _slots.size() && _slots[key.index].generation == key.generation);
    }

    /**
     * \brief Access Operator
     *
     * \param[in] key The key to retrieve
     * \returns A reference to the data pointed to by the key
     *
     * Aborts if the given key is not valid.
     */
    const V& operator[](const Key& key) const
    {
      const Slot& slot = try_get_slot(key);
      ASSERTM(slot.index < _data.size(), "Invalid slot index.");
      return _data[slot.index];
    }

    /**
     * \brief Access Operator
     *
     * \param[in] key The key to retrieve
     * \returns A reference to the data pointed to by the key
     *
     * Aborts if the given key is not valid.
     */
    V& operator[](const Key& key)
    {
      Slot& slot = try_get_slot(key);
      ASSERTM(slot.index < _data.size(), "Invalid slot index.");
      return _data[slot.index];
    }

    /**
     * \brief Reserve at least \c capacity space in the \c SlotMap
     *
     * \param[in] capacity Minimum capacity
     */
    void reserve(Index capacity)
    {
      _slots.reserve(capacity + 1);
      _data.reserve(capacity);
      _slots.reserve(capacity);
    }

    /**
     * \brief Retrieve current capacity of the \c SlotMap
     *
     * \returns Capacity of the \c SlotMap
     */
    Index capacity() const
    {
      return _data.capacity();
    }

    /**
     * \brief Retrieve current size of the \c SlotMap
     *
     * \returns Size of the \c SlotMap
     */
    Index size() const
    {
      return _data.size();
    }

    V* data()
    {
      return _data.data();
    }

    const V* data() const
    {
      return _data.data();
    }

    std::size_t bytes() const
    {
      return sizeof(freelist) + _slots.size() * sizeof(Slot) + _data.size() * sizeof(V) +
             _erase.size() * sizeof(Index);
    }

    /**
     * \brief Iterator type for the \c SlotMap
     *
     * Iterates over the data array in order.
     * Returns pairs consisting of references to the data and the corresponding keys.
     *
     * Note that this provides a more idiomatic api than forwarding the std::vector
     * iterator. Doing so would expose a pointer into the data array to the user
     * and allow access via indices, bypassing the slotmaps keys. The user
     * can still retrieve a pointer from the passed out references, but doing
     * so is a more obvious code smell.
     *
     * \author Markus Muegge
     */
    class SlotMapIterator
    {
    public:
      SlotMapIterator(SlotMap<V, Key>& map, Index i) : _slotmap(map), _index(i)
      {
      }

      bool operator==(const SlotMapIterator& rhs)
      {
        // Iterators are identical if they refer to the same slotmap and point to the same indices
        return &(this->_slotmap) == &(rhs._slotmap) && this->_index == rhs._index;
      }
      bool operator!=(const SlotMapIterator& rhs)
      {
        return !(*this == rhs);
      }

      SlotMapIterator& operator++()
      {
        _index++;
        return *this;
      }

      KeyValuePair<V, Key> operator*()
      {
        // Retrieve slot for current index via erase table
        std::uint32_t slot_index = _slotmap._erase[_index];
        Slot& slot = _slotmap._slots[slot_index];
        return KeyValuePair<V, Key>(_slotmap._data[slot.index], Key(slot_index, slot.generation));
      }

    private:
      SlotMap<V, Key>& _slotmap;
      Index _index;
    };

    /**
     * \brief Const iterator type for the \c SlotMap
     *
     * See \c SlotMapIterator for details
       *
     * \author Markus Muegge
     */
    class ConstSlotMapIterator
    {
    public:
      ConstSlotMapIterator(const SlotMap<V, Key>& map, Index i) : _slotmap(map), _index(i)
      {
      }

      bool operator==(const ConstSlotMapIterator& rhs)
      {
        // Iterators are identical if they refer to the same slotmap and point to the same indices
        return &(this->_slotmap) == &(rhs._slotmap) && this->_index == rhs._index;
      }
      bool operator!=(const ConstSlotMapIterator& rhs)
      {
        return !(*this == rhs);
      }

      ConstSlotMapIterator& operator++()
      {
        _index++;
        return *this;
      }

      ConstKeyValuePair<V, Key> operator*()
      {
        // Retrieve slot for current index via erase table
        std::uint32_t slot_index = _slotmap._erase[_index];
        const Slot& slot = _slotmap._slots[slot_index];
        return ConstKeyValuePair<V, Key>(_slotmap._data[slot.index], Key(slot_index, slot.generation));
      }

    private:
      const SlotMap<V, Key>& _slotmap;
      Index _index;
    };

    SlotMapIterator begin()
    {
      return SlotMapIterator(*this, 0);
    }

    SlotMapIterator end()
    {
      return SlotMapIterator(*this, _data.size());
    }

    ConstSlotMapIterator begin() const
    {
      return ConstSlotMapIterator(*this, 0);
    }

    ConstSlotMapIterator end() const
    {
      return ConstSlotMapIterator(*this, _data.size());
    }

  private:
    /**
     * \brief Internal data structure for slots
     */
    struct Slot
    {
      Slot() : index(0), generation(1)
      {
      }
      Slot(std::uint32_t i, std::uint32_t g) : index(i), generation(g)
      {
      }

      std::uint32_t index;
      std::uint32_t generation;
    };

    /**
     * \brief Helper method to retrieve a slot
     *
     * \param[in] key Key
     * \returns Slot pointed to by key
     *
     * Aborts if given key is not valid, returns reference to slot otherwise.
     */
    Slot& try_get_slot(const Key& key)
    {
      ASSERTM(key.index < _slots.size(), "Invalid SlotMap key. Index exceeds number of slots.");
      Slot& slot = _slots[key.index];
      ASSERTM(key.generation == slot.generation, "Invalid SlotMap key. Generation mismatch.");
      return slot;
    }

    const Slot& try_get_slot(const Key& key) const
    {
      ASSERTM(key.index < _slots.size(), "Invalid SlotMap key. Index exceeds number of slots.");
      const Slot& slot = _slots[key.index];
      ASSERTM(key.generation == slot.generation, "Invalid SlotMap key. Generation mismatch.");
      return slot;
    }

    struct
    {
      std::uint32_t head;
      std::uint32_t tail;
    } freelist;

    std::vector<Slot> _slots;
    std::vector<V> _data;
    std::vector<std::uint32_t> _erase;
  };

  template<typename V, typename Key = SlotMapKey>
  class SecondaryMap
  {
  public:
    void insert(const Key& key, const V& value)
    {
      if(key.index >= _slots.size())
      {
        Index extend_by = 1 + key.index - _slots.size();
        _slots.insert(_slots.end(), extend_by, Slot());
      }

      Slot& slot = _slots[key.index];
      slot.value = std::optional<V>(value);
      slot.generation = key.generation;
    }

    bool contains_key(const Key& key)
    {
      return key.index < _slots.size() && _slots[key.index].generation == key.generation &&
             _slots[key.index].value.has_value();
    }

    /**
     * \brief Access Operator
     *
     * \param[in] key The key to retrieve
     * \returns A reference to the data pointed to by the key
     *
     * Aborts if the given key is not valid.
     */
    const V& operator[](const Key& key) const
    {
      ASSERTM(key.index < _slots.size(), "Invalid SecondaryMap key. Index exceeds number of slots.");
      const Slot& slot = _slots[key.index];
      ASSERTM(slot.value.has_value(), "Invalid SecondaryMap key. Slot is unoccupied.");
      ASSERTM(key.generation == slot.generation, "Invalid SecondaryMap key. Generation mismatch.");
      return slot.value.value();
    }

    /**
     * \brief Access Operator
     *
     * \param[in] key The key to retrieve
     * \returns A reference to the data pointed to by the key
     *
     * Aborts if the given key is not valid.
     */
    V& operator[](const Key& key)
    {
      ASSERTM(key.index < _slots.size(), "Invalid SecondaryMap key. Index exceeds number of slots.");
      Slot& slot = _slots[key.index];
      ASSERTM(slot.value.has_value(), "Invalid SecondaryMap key. Slot is unoccupied.");
      ASSERTM(key.generation == slot.generation, "Invalid SecondaryMap key. Generation mismatch.");
      return slot.value.value();
    }

  private:
    struct Slot
    {
      Slot() : generation(1), value(std::nullopt)
      {
      }
      explicit Slot(V&& v) : generation(1), value(std::move(v))
      {
      }

      std::uint32_t generation;
      std::optional<V> value;
    };

    std::vector<Slot> _slots;
  };
} // namespace FEAT::Util
