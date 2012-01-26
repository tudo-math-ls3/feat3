#pragma once
#ifndef KERNEL_UTIL_INDEX_TABLE_HPP
#define KERNEL_UTIL_INDEX_TABLE_HPP 1

// includes, FEAST
#include <kernel/util/assertion.hpp>

namespace FEAST
{
  /**
   * \brief Index table implementation
   *
   * \todo detailed description
   *
   * \author Peter Zajac
   */
  class IndexTable
  {
  public:
    /// ImageIterator typedef for Adjactor interface implementation
    typedef const Index* ImageIterator;

  private:
    /// number of rows in the table
    Index _num_rows;
    /// number of columns in the table
    Index _num_cols;
    /// table index array
    Index* _table;

  public:
    /// Standard constructor
    IndexTable() :
      _num_rows(0),
      _num_cols(0),
      _table(nullptr)
    {
      CONTEXT("IndexTable::IndexTable()");
    }

    /**
     * \brief Constructor
     *
     * \param[in] num_rows, num_cols
     * The dimensions of the table.
     */
    IndexTable(
      Index num_rows,
      Index num_cols)
       :
      _num_rows(num_rows),
      _num_cols(num_cols),
      _table(nullptr)
    {
      CONTEXT("IndexTable::IndexTable(Index,Index)");
      ASSERT(num_rows > 0, "num_rows must be greater than 0");
      ASSERT(num_cols > 0, "num_cols must be greater than 0");
      _table = new Index[num_rows * num_cols];
    }

    /// virtual destructor
    virtual ~IndexTable()
    {
      CONTEXT("IndexTable::~IndexTable()");
      if(_table != nullptr)
      {
        delete [] _table;
      }
    }

    /**
     * \brief Creates the index table.
     *
     * \param[in] num_rows, num_cols
     * The dimensions of the index table.
     */
    void create(
      Index num_rows,
      Index num_cols)
    {
      CONTEXT("IndexTable::create()");
      ASSERT(num_rows > 0, "num_rows must be greater than 0");
      ASSERT(num_cols > 0, "num_cols must be greater than 0");
      if(_table != nullptr)
      {
        delete [] _table;
      }
      _num_rows = num_rows;
      _num_cols = num_cols;
      _table = new Index[num_rows * num_cols];
    }

    /**
     * \brief Returns a table entry.
     *
     * \param[in] row
     * The row index of the table entry.
     *
     * \param[in] col
     * The column index of the table entry.
     *
     * \returns
     * A reference to the table entry.
     */
    inline Index& at(
      Index row,
      Index col)
    {
      CONTEXT("IndexTable::at()");
      ASSERT(row < _num_rows, "row index out of range");
      ASSERT(col < _num_cols, "column index out of range");
      return _table[row*_num_cols + col];
    }

    /** \copydoc at() */
    inline const Index& at(
      Index row,
      Index col) const
    {
      CONTEXT("IndexTable::at() [const]");
      ASSERT(row < _num_rows, "row index out of range");
      ASSERT(col < _num_cols, "column index out of range");
      return _table[row * _num_cols + col];
    }

    /**
     * \brief Returns the number of rows.
     * \returns
     * The number of rows in the table.
     */
    inline Index get_num_rows() const
    {
      CONTEXT("IndexTable::get_num_rows()");
      return _num_rows;
    }

    /**
     * \brief Returns the number of columns.
     * \returns
     * The number of columns in the table.
     */
    inline Index get_num_cols() const
    {
      CONTEXT("IndexTable::get_num_cols()");
      return _num_cols;
    }

    /**
     * \brief Returns the table's index array.
     * \returns
     * The table's index array.
     */
    inline Index* get_table()
    {
      CONTEXT("IndexTable::get_table()");
      return _table;
    }

    /** \copydoc get_table() */
    inline const Index* get_table() const
    {
      CONTEXT("IndexTable::get_table()");
      return _table;
    }

    /**
     * \brief Direct access operator.
     *
     * \param[in] row
     * The row that is to be returned.
     *
     * \returns
     * A pointer to the row of the table.
     */
    inline Index* operator[](Index row)
    {
      CONTEXT("IndexTable::operator[]()");
      ASSERT(row < _num_rows, "row index out of range");
      return &_table[row * _num_cols];
    }

    /** \copydoc operator[]() */
    inline const Index* operator[](Index row) const
    {
      CONTEXT("IndexTable::operator[]()");
      ASSERT(row < _num_rows, "row index out of range");
      return &_table[row * _num_cols];
    }

    /* ******************************************************************* */
    /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
    /* ******************************************************************* */
    /** \copydoc Adjactor::num_nodes_domain() */
    inline Index num_nodes_domain() const
    {
      CONTEXT("IndexTable::num_nodes_domain()");
      return _num_rows;
    }

    /** \copydoc Adjactor::num_nodes_image() */
    inline Index num_nodes_image() const
    {
      CONTEXT("IndexTable::num_nodes_image()");
      return _num_cols;
    }

    /** \copydoc Adjactor::image_begin() */
    inline ImageIterator image_begin(Index domain_node) const
    {
      CONTEXT("IndexTable::image_begin()");
      ASSERT(domain_node < _num_rows, "Domain node index out of range");
      return &_table[domain_node * _num_cols];
    }

    /** \copydoc Adjactor::image_end() */
    inline ImageIterator image_end(Index domain_node) const
    {
      CONTEXT("IndexTable::image_end()");
      ASSERT(domain_node < _num_rows, "Domain node index out of range");
      return &_table[(domain_node+1) * _num_cols];
    }
  }; // class IndexTable
} // namespace FEAST

#endif // KERNEL_UTIL_INDEX_TABLE_HPP
