#pragma once
#ifndef KERNEL_BASE_MESH_ITEM_HPP
#define KERNEL_BASE_MESH_ITEM_HPP

// includes, system
#include <sstream> // for std::ostream
#include <cassert>  // for assert()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/constants.hpp>

namespace FEAST
{
  namespace BaseMesh
  {
    /**
    * \brief Base class for every item that can be part of a mesh (vertex, edges, face, volume)
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    */
    class Item
    {

    private:

      /// index
      global_index_t _index;

      /// number
      global_index_t _number;

    public:

      /// CTOR 1
      Item()
      {
        _index = Constants::MAX_INDEX;
        _number = Constants::MAX_NUMBER;
      }

      /// CTOR 2
      Item(global_index_t index, global_index_t number)
        : _index(index),
          _number(number)
      {
      }

      /// DTOR (must be virtual)
      virtual ~Item()
      {
      }

      /// returns the index of this item
      inline global_index_t index() const
      {
        return _index;
      }

      /// sets the index of this item
      inline void set_index(const global_index_t index)
      {
        assert(index < Constants::MAX_INDEX);
        _index = index;
      }

      /// returns the number of this item
      inline global_index_t number() const
      {
        return _number;
      }

      /// sets the number of this item
      inline void set_number(const global_index_t number)
      {
        assert(number < Constants::MAX_NUMBER);
        _number = number;
      }

      /// prints the index of this item to the given stream
      inline void print_index(std::ostream& stream) const
      {
        if (index() != Constants::MAX_INDEX)
        {
          stream << index();
        }
        else
        {
          // for debugging purpose, print the address of the object if no index is set yet
          stream << (long)this;
          //stream << "-";
        }
        if (number() != Constants::MAX_NUMBER)
        {
          stream << "/" << number();
        }
      }

      /// returns the index of this item as string
      inline std::string print_index() const
      {
        std::ostringstream oss;
        print_index(oss);
        return oss.str();
      }
    };  // class Item
  } // namespace BaseMesh
} // namespace FEAST
#endif // KERNEL_BASE_MESH_ITEM_HPP
