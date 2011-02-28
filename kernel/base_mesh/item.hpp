#pragma once
#ifndef KERNEL_BASE_MESH_ITEM_HPP
#define KERNEL_BASE_MESH_ITEM_HPP

// includes, system
#include <sstream> // for std::ostringstream
#include <cassert>  // for assert()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/constants.hpp>

/// FEAST namespace
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
      index_t_glob _index;

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
      index_t_glob _number;

    public:

      /// CTOR 1
      Item()
      {
        CONTEXT("BaseMesh::Item::Item()");
        _index = Constants::MAX_INDEX;
        _number = Constants::MAX_NUMBER;
      }

      /// CTOR 2
      Item(index_t_glob index, index_t_glob number)
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
      inline index_t_glob index() const
      {
        CONTEXT("BaseMesh::Item::index()");
        return _index;
      }

      /// sets the index of this item
      inline void set_index(index_t_glob const index)
      {
        CONTEXT("BaseMesh::Item::set_index()");
        assert(index < Constants::MAX_INDEX);
        _index = index;
      }

      /// returns the number of this item
      inline index_t_glob number() const
      {
        CONTEXT("BaseMesh::Item::number()");
        return _number;
      }

      /// sets the number of this item
      inline void set_number(index_t_glob const number)
      {
        CONTEXT("BaseMesh::Item::set_number()");
        assert(number < Constants::MAX_NUMBER);
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
      inline std::string print_index() const
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
