/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR below.
 * Generally, Peter wanted to take a deeper look at the base mesh implementation.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_BASE_MESH_ITEM_HPP
#define KERNEL_BASE_MESH_ITEM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/constants.hpp>

// includes, system
#include <sstream> // for std::ostringstream

namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
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

      /**
      * \brief index of the item
      *
      * The index is set when the item is added to the base mesh. It will never be changed again until destruction
      * of the item. (Currently, the index is simply the position in the corresponding vector in the BaseMesh class
      * (_vertices, _edges, ...)
      */
      Index _index;

      /**
      * \brief number of the item
      *
      * While all items have an index, only *active* items have a number. When a mesh is created, indices and numbers
      * usually are equal, but as soon as one cell is subdivided, indices and numbers will differ.
      *
      * \sa BaseMesh2D::set_cell_numbers()
      */
// COMMENT_HILMAR: Brauchen wirklich alle items eine number? Oder nur die Zellen groesster Dimension? Wenn letzters,
// dann wird die Variable in die CellData<...> Klasse verschoben.
      Index _number;

    public:

      /// CTOR 1
      Item()
      {
        CONTEXT("BaseMesh::Item::Item()");
        _index = Constants::MAX_INDEX;
        _number = Constants::MAX_NUMBER;
      }

      /// CTOR 2
      Item(Index index, Index number)
        : _index(index),
          _number(number)
      {
        CONTEXT("BaseMesh::Item::Item()");
      }

      /// DTOR
      virtual ~Item()
      {
        CONTEXT("BaseMesh::Item::~Item()");
      }

      /// returns the index of this item
      inline Index index() const
      {
        CONTEXT("BaseMesh::Item::index()");
        return _index;
      }

      /// sets the index of this item
      inline void set_index(Index const index)
      {
        CONTEXT("BaseMesh::Item::set_index()");
        ASSERT(index < Constants::MAX_INDEX, "Index " + stringify(index) + " must not exceed Constants::MAX_INDEX.");
        _index = index;
      }

      /// returns the number of this item
      inline Index number() const
      {
        CONTEXT("BaseMesh::Item::number()");
        return _number;
      }

      /// sets the number of this item
      inline void set_number(Index const number)
      {
        CONTEXT("BaseMesh::Item::set_number()");
        ASSERT(number < Constants::MAX_INDEX, "Number " + stringify(number) + " must not exceed Constants::MAX_INDEX.");
        _number = number;
      }

      /// unsets the number of this item, i.e. sets it to Constants::MAX_NUMBER
      inline void unset_number()
      {
        CONTEXT("BaseMesh::Item::unset_number()");
        _number = Constants::MAX_NUMBER;
      }

      /// prints the index of this item to the given stream
      inline void print_index(std::ostream& stream) const
      {
        CONTEXT("BaseMesh::Item::print_index()");
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
      inline String print_index() const
      {
        CONTEXT("BaseMesh::Item::print_index()");
        std::ostringstream oss;
        print_index(oss);
        return oss.str();
      }
    };  // class Item
  } // namespace BaseMesh
} // namespace FEAST
#endif // KERNEL_BASE_MESH_ITEM_HPP
