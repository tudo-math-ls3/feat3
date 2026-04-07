// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Matrix Mirror Buffer class template.
     *
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     *
     * This class holds all matrix mirror informations needed on the corresponding other rank. \n \n
     * The data layout resembles the common csr layout, with the exception of multiple data entries per non zero element. \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _indices[0]: row start index (including matrix end index)\n
     * _indices[1]: column index per non zero element \n
     * _scalar_index[0]: row count \n
     * _scalar_index[1]: column count \n
     * _scalar_index[2]: non zero element count (used elements) \n
     * _scalar_index[3]: entries per nonzero count
     *
     * \author Dirk Ribbrock
     */
    template<typename DT_, typename IT_>
    class MatrixMirrorBuffer : public Container<DT_, IT_>
    {
    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our value type
      typedef DT_ ValueType;
      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = MatrixMirrorBuffer<DT2_, IT2_>;

      /// this typedef lets you create a buffer container with different Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional buffer.
       */
      MatrixMirrorBuffer() = default;

      /**
       * \brief Basic Constructor
       *
       * \param[in] rows_in The row count.
       * \param[in] columns_in The column count.
       * \param[in] used_elements_in The amount of non zero elements in the created matrix mirror buffer.
       * \param[in] entries_per_nonzero_in The amount of data entries per nonzero matrix entry.
       *
       * Creates a new MatrixMirrorBuffer with the given dimensions.
       */
      explicit MatrixMirrorBuffer(Index rows_in, Index columns_in, Index used_elements_in, Index entries_per_nonzero_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(used_elements_in);
        this->_scalar_index.push_back(entries_per_nonzero_in);

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * (num_rows() + 1)));
        this->_indices_size.push_back(num_rows() + 1);

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * num_nzes()));
        this->_indices_size.push_back(num_nzes());

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * (num_nzes() * entries_per_nonzero())));
        this->_elements_size.push_back(num_nzes() * entries_per_nonzero());
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix mirror buffer from
       *
       * Creates a matrix mirror buffer based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit MatrixMirrorBuffer(const Adjacency::Graph & graph, Index entries_per_nonzero_in) :
        MatrixMirrorBuffer(graph.get_num_nodes_domain(), graph.get_num_nodes_image(), graph.get_num_indices(), entries_per_nonzero_in)
      {
        const Index num_nnze = graph.get_num_indices();
        const Index num_rows = graph.get_num_nodes_domain();
        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());

        Memory::TypedView<IT_> prow_ptr = this->row_ptr_view_w();
        Memory::TypedView<IT_> pcol_idx = this->col_idx_view_w();

        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i <= num_rows; ++i)
          prow_ptr[i] = IT_(dom_ptr[i]);

        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < num_nnze; ++i)
          pcol_idx[i] = IT_(img_idx[i]);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count.
       * \param[in] columns_in The column count.
       * \param[in] entries_per_nonzero_in The amount of data entries per nonzero matrix entry.
       * \param[in] col_idx_in Vector with column indices.
       * \param[in] val_in Vector with non zero elements.
       * \param[in] row_ptr_in Vector with start indices of all num_rows into the val/col_idx arrays.
       * Note that this vector must also contain the end index of the last row and thus has a size of row_count + 1.
       *
       * Creates a new MatrixMirrorBuffer with the given dimensions and content.
       */
      explicit MatrixMirrorBuffer(const Index rows_in, const Index columns_in, Index entries_per_nonzero_in,
        DenseVector<IT_, IT_> & col_idx_in, DenseVector<DT_, IT_> & val_in, DenseVector<IT_, IT_> & row_ptr_in)
      {
        XASSERT(val_in.size() % entries_per_nonzero_in == 0);
        XASSERT(val_in.size() == col_idx_in.size() * entries_per_nonzero_in);

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(col_idx_in.size());
        this->_scalar_index.push_back(entries_per_nonzero_in);

        this->_indices.push_back(row_ptr_in.elements_arbiter().attach());
        this->_indices_size.push_back(row_ptr_in.size());
        this->_indices.push_back(col_idx_in.elements_arbiter().attach());
        this->_indices_size.push_back(col_idx_in.size());
        this->_elements.push_back(val_in.elements_arbiter().attach());
        this->_elements_size.push_back(val_in.size());
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source buffer.
       *
       * Moves another buffer to this buffer.
       */
      MatrixMirrorBuffer(MatrixMirrorBuffer && other) :
        Container<DT_, IT_>(std::forward<MatrixMirrorBuffer>(other))
      {
      }

      /**
       * \brief Assignment move operator
       *
       * \param[in] other The source buffer.
       *
       * Moves another buffer to the target buffer.
       */
      MatrixMirrorBuffer& operator= (MatrixMirrorBuffer && other)
      {
        this->move(std::forward<MatrixMirrorBuffer>(other));
        return *this;
      }

      /** \brief Clone operation
       *
       * Create a clone of this container.
       *
       * \param[in] clone_mode The actual cloning procedure.
       *
       */
      MatrixMirrorBuffer clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        MatrixMirrorBuffer t;
        t.clone(*this, clone_mode);
        return t;
      }

      /** \brief Clone operation
       *
       * Create a clone of another container.
       *
       * \param[in] other The source container to create the clone from.
       * \param[in] clone_mode The actual cloning procedure.
       *
       */
      template<typename DT2_, typename IT2_>
      void clone(const MatrixMirrorBuffer<DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Weak)
      {
        Container<DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source MatrixMirrorBuffer.
       *
       * Use source matrix mirror buffers content as content of current matrix mirror buffer
       */
      template <typename DT2_, typename IT2_>
      void convert(const MatrixMirrorBuffer<DT2_, IT2_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count.
       */
      Index num_rows() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(0);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      Index num_cols() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      Index num_nzes() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve entries per non zero element count.
       *
       * \returns Entries per non zero element count.
       */
      Index entries_per_nonzero() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve total length of value array.
       *
       * That's simply the product of non-zero count and number of entries per non-zero.
       *
       * \returns Total length of value array
       */
      Index num_nzes_raw() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(2) * this->_scalar_index.at(3);
      }


      /// Returns a reference to the row-pointer array arbiter
      Memory::Arbiter& row_ptr_arbiter()
      {
        return this->_indices.at(0);
      }

      /// Returns a reference to the row-pointer array arbiter
      const Memory::Arbiter& row_ptr_arbiter() const
      {
        return this->_indices.at(0);
      }

      /// Returns a reference to the column-index array arbiter
      Memory::Arbiter& col_idx_arbiter()
      {
        return this->_indices.at(1);
      }

      /// Returns a reference to the column-index array arbiter
      const Memory::Arbiter& col_idx_arbiter() const
      {
        return this->_indices.at(1);
      }

      /// Returns a reference to the values array arbiter
      Memory::Arbiter& val_arbiter()
      {
        return this->_elements.front();
      }

      /// Returns a reference to the values array arbiter
      const Memory::Arbiter& val_arbiter() const
      {
        return this->_elements.front();
      }

      Memory::TypedView<DT_> val_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<DT_> val_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<DT_> val_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<IT_>(this->_elements.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<DT_> val_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<IT_>(this->_elements.at(0).view(loc, acc));
      }

      Memory::TypedView<IT_> row_ptr_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<IT_> row_ptr_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<IT_> row_ptr_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<IT_> row_ptr_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, acc));
      }

      Memory::TypedView<IT_> col_idx_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::read));
      }

      Memory::TypedView<IT_> col_idx_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::write));
      }

      Memory::TypedView<IT_> col_idx_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<IT_> col_idx_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, acc));
      }
    }; // class MatrixMirrorBuffer<...>
  } // namespace LAFEM
} // namespace FEAT
