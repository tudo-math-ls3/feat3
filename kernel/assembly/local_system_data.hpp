#pragma once
#ifndef KERNEL_ASSEMBLY_LOCAL_SYSTEM_DATA_HPP
#define KERNEL_ASSEMBLY_LOCAL_SYSTEM_DATA_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Local Vector data class template
     *
     * This class is used by various assembly algorithms for the handling of local vectors.
     *
     * \tparam LocalVector_
     * The underlying local vector type. Shall be a Tiny::Vector object.
     *
     * \tparam Mapping_
     * The underlying dof-mapping type.
     *
     * \author Peter Zajac
     */
    template<
      typename LocalVector_,
      typename Mapping_>
    class LocalVectorData :
      public LocalVector_
    {
    public:
      /// data type
      typedef typename LocalVector_::DataType DataType;

    protected:
      /// dof-mapping reference
      const Mapping_& _mapping;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mapping
       * A reference to the dof-mapping object that is to be used.
       */
      explicit LocalVectorData(const Mapping_& mapping) :
        _mapping(mapping)
      {
      }

      /**
       * \brief Returns the number of entries in the local vector.
       */
      Index get_num_entries() const
      {
        return _mapping.get_num_local_dofs();
      }

      /**
       * \brief Returns the number of contributions for a local entry.
       *
       * \param[in] entry
       * The index of the entry whose contribution count is to be returned.
       *
       * \returns
       * The number of contributions for \p entry.
       */
      Index get_num_contribs(Index entry) const
      {
        return _mapping.get_num_contribs(entry);
      }

      /**
       * \brief Returns the mapped index for an entry contribution.
       *
       * \param[in] entry
       * The index of the local entry that is to be mapped.
       *
       * \param[in] contrib
       * The index of the contribution that is to be mapped.
       *
       * \returns
       * The mapped index of the entry contribution.
       */
      Index get_index(Index entry, Index contrib) const
      {
        return _mapping.get_index(entry, contrib);
      }

      /**
       * \brief Returns the mapped weight for an entry contribution.
       *
       * \param[in] entry
       * The index of the local entry that is to be mapped.
       *
       * \param[in] contrib
       * The index of the contribution that is to be mapped.
       *
       * \returns
       * The mapped weight of the entry contribution.
       */
      DataType get_weight(Index entry, Index contrib) const
      {
        return DataType(_mapping.get_weight(entry, contrib));
      }
    };

    /**
     * \brief Local Matrix data class template
     *
     * This class is used by various assembly algorithms for the handling of local matrices.
     *
     * \tparam LocalMatrix_
     * The underlying local matrix type. Shall be a Tiny::Matrix object.
     *
     * \tparam RowMapping_
     * The underlying dof-mapping type for the rows.
     *
     * \tparam ColMapping_
     * the underlying dof-mapping type for the columns.
     *
     * \author Peter Zajac
     */
    template<
      typename LocalMatrix_,
      typename RowMapping_,
      typename ColMapping_>
    class LocalMatrixData :
      public LocalMatrix_
    {
    public:
      /// data type
      typedef typename LocalMatrix_::DataType DataType;

    protected:
      /// row dof-mapping reference
      const RowMapping_& _row_map;
      /// column dof-mapping reference
      const ColMapping_& _col_map;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] row_map, col_map
       * References to the row and column dof-mapping objects that are to be used.
       */
      explicit LocalMatrixData(const RowMapping_& row_map, const ColMapping_& col_map) :
        _row_map(row_map),
        _col_map(col_map)
      {
      }

      /**
       * \brief Returns the number of rows in the local matrix.
       */
      Index get_num_rows() const
      {
        return _row_map.get_num_local_dofs();
      }

      /**
       * \brief Returns the number of columns in the local matrix.
       */
      Index get_num_cols() const
      {
        return _col_map.get_num_local_dofs();
      }

      /**
       * \brief Returns the number of contributions for a local row.
       *
       * \param[in] row
       * The index of the row whose contribution count is to be returned.
       *
       * \returns
       * The number of contributions for \p row.
       */
      Index get_num_row_contribs(Index row) const
      {
        return _row_map.get_num_contribs(row);
      }

      /**
       * \brief Returns the number of contributions for a local column.
       *
       * \param[in] col
       * The index of the column whose contribution count is to be returned.
       *
       * \returns
       * The number of contributions for \p col.
       */
      Index get_num_col_contribs(Index col) const
      {
        return _col_map.get_num_contribs(col);
      }

      /**
       * \brief Returns the mapped index for a row contribution.
       *
       * \param[in] row
       * The index of the local row that is to be mapped.
       *
       * \param[in] contrib
       * The index of the contribution that is to be mapped.
       *
       * \returns
       * The mapped index of the row contribution.
       */
      Index get_row_index(Index row, Index contrib) const
      {
        return _row_map.get_index(row, contrib);
      }

      /**
       * \brief Returns the mapped index for a column contribution.
       *
       * \param[in] col
       * The index of the local column that is to be mapped.
       *
       * \param[in] contrib
       * The index of the contribution that is to be mapped.
       *
       * \returns
       * The mapped index of the column contribution.
       */
      Index get_col_index(Index col, Index contrib) const
      {
        return _col_map.get_index(col, contrib);
      }

      /**
       * \brief Returns the mapped weight for a row contribution.
       *
       * \param[in] row
       * The index of the local row that is to be mapped.
       *
       * \param[in] contrib
       * The index of the contribution that is to be mapped.
       *
       * \returns
       * The mapped weight of the row contribution.
       */
      DataType get_row_weight(Index row, Index contrib) const
      {
        return DataType(_row_map.get_weight(row, contrib));
      }

      /**
       * \brief Returns the mapped weight for an column contribution.
       *
       * \param[in] col
       * The index of the local column that is to be mapped.
       *
       * \param[in] contrib
       * The index of the contribution that is to be mapped.
       *
       * \returns
       * The mapped weight of the column contribution.
       */
      DataType get_col_weight(Index col, Index contrib) const
      {
        return DataType(_col_map.get_weight(col, contrib));
      }
    }; // class LocalMatrixData<...>
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_LOCAL_SYSTEM_DATA_HPP
