// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_MATRIX_MIRROR_BUFFER_HPP
#define KERNEL_LAFEM_MATRIX_MIRROR_BUFFER_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/adjacency/graph.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Matrix Mirror Buffer class template.
     *
     * \tparam Mem_ The \ref FEAT::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     *
     * This class holds all matrix mirror informations needed on the corresponding other rank. \n \n
     * The data layout resembles the common csr layout, with the exception of multiple data entries per non zero element. \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _indices[0]: column index per non zero element \n
     * _indices[1]: row start index (including matrix end index)\n
     * _scalar_index[0]: container size
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_index[4]: entries per nonzero count
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class MatrixMirrorBuffer : public Container<Mem_, DT_, IT_>
    {
    public:

      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Our value type
      typedef DT_ ValueType;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = MatrixMirrorBuffer<Mem2_, DT2_, IT2_>;

      /// this typedef lets you create a buffer container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional buffer.
       */
      explicit MatrixMirrorBuffer() :
        Container<Mem_, DT_, IT_> (0)
      {
      }

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
      MatrixMirrorBuffer(Index rows_in, Index columns_in, Index used_elements_in, Index entries_per_nonzero_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(used_elements_in);
        this->_scalar_index.push_back(entries_per_nonzero_in);

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(used_elements()));
        this->_indices_size.push_back(used_elements());

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(rows() + 1));
        this->_indices_size.push_back(rows() + 1);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(used_elements() * entries_per_nonzero()));
        this->_elements_size.push_back(used_elements() * entries_per_nonzero());
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix mirror buffer from
       *
       * Creates a matrix mirror buffer based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit MatrixMirrorBuffer(const Adjacency::Graph & graph, Index entries_per_nonzero_in) :
        Container<Mem_, DT_, IT_>(0)
      {
        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nnze = graph.get_num_indices();

        // create temporary vectors
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vrow_ptr(num_rows+1);
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vcol_idx(num_nnze);
        LAFEM::DenseVector<Mem::Main, DT_, IT_> vdata(num_nnze * entries_per_nonzero_in, DT_(0));

        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());
        IT_ * prow_ptr(vrow_ptr.elements());
        IT_ * pcol_idx(vcol_idx.elements());

        // build row-end
        prow_ptr[0] = IT_(dom_ptr[0]);
        for(Index i(0); i < num_rows; ++i)
          prow_ptr[i+1] = IT_(dom_ptr[i+1]);

        // build col-idx
        for(Index i(0); i < num_nnze; ++i)
          pcol_idx[i] = IT_(img_idx[i]);

        // build the matrix
        this->assign(MatrixMirrorBuffer<Mem::Main, DT_, IT_>(num_rows, num_cols, entries_per_nonzero_in, vcol_idx, vdata, vrow_ptr));
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count.
       * \param[in] columns_in The column count.
       * \param[in] entries_per_nonzero_in The amount of data entries per nonzero matrix entry.
       * \param[in] col_ind_in Vector with column indices.
       * \param[in] val_in Vector with non zero elements.
       * \param[in] row_ptr_in Vector with start indices of all rows into the val/col_ind arrays.
       * Note that this vector must also contain the end index of the last row and thus has a size of row_count + 1.
       *
       * Creates a new MatrixMirrorBuffer with the given dimensions and content.
       */
      explicit MatrixMirrorBuffer(const Index rows_in, const Index columns_in, Index entries_per_nonzero_in,
                               DenseVector<Mem_, IT_, IT_> & col_ind_in, DenseVector<Mem_, DT_, IT_> & val_in, DenseVector<Mem_, IT_, IT_> & row_ptr_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        XASSERT(val_in.size() % entries_per_nonzero_in == 0);
        XASSERT(val_in.size() == col_ind_in.size() * entries_per_nonzero_in);

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(col_ind_in.size());
        this->_scalar_index.push_back(entries_per_nonzero_in);

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(col_ind_in.elements());
        this->_indices_size.push_back(col_ind_in.size());
        this->_indices.push_back(row_ptr_in.elements());
        this->_indices_size.push_back(row_ptr_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix mirror buffer from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit MatrixMirrorBuffer(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        deserialize<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source buffer.
       *
       * Moves another buffer to this buffer.
       */
      MatrixMirrorBuffer(MatrixMirrorBuffer && other) :
        Container<Mem_, DT_, IT_>(std::forward<MatrixMirrorBuffer>(other))
      {
      }

      /**
       * \brief Assignment move operator
       *
       * \param[in] other The source buffer.
       *
       * Moves another buffer to the target buffer.
       */
      MatrixMirrorBuffer & operator= (MatrixMirrorBuffer && other)
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
      template<typename Mem2_, typename DT2_, typename IT2_>
      void clone(const MatrixMirrorBuffer<Mem2_, DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Weak)
      {
        Container<Mem_, DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source MatrixMirrorBuffer.
       *
       * Use source matrix mirror buffers content as content of current matrix mirror buffer
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const MatrixMirrorBuffer<Mem2_, DT2_, IT2_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Deserialization of complete container entity.
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialize(std::vector<char> input)
      {
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_csr, input);
      }

      /**
       * \brief Serialization of complete container entity.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialize for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialize()
      {
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_csr);
      }

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count.
       */
      template <Perspective = Perspective::native>
      Index rows() const
      {
        return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      template <Perspective = Perspective::native>
      Index columns() const
      {
        return this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      template <Perspective = Perspective::native>
      Index used_elements() const
      {
        return this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve entries per non zero element count.
       *
       * \returns Entries per non zero element count.
       */
      template <Perspective = Perspective::native>
      Index entries_per_nonzero() const
      {
        return this->_scalar_index.at(4);
      }

      /**
       * \brief Retrieve total length of value array.
       *
       * That's simply the product of non-zero count and number of entries per non-zero.
       *
       * \returns Total length of value array
       */
      template <Perspective = Perspective::native>
      Index val_size() const
      {
        return this->_scalar_index.at(3) * this->_scalar_index.at(4);
      }

      /**
       * \brief Retrieve column indices array.
       *
       * \returns Column indices array.
       */
      IT_ * col_ind()
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(0);
      }

      IT_ const * col_ind() const
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(0);
      }

      /**
       * \brief Retrieve non zero element array.
       *
       * \returns Non zero element array.
       */
      DT_ * val()
      {
        if (this->_elements.size() == 0)
          return nullptr;

        return this->_elements.at(0);
      }

      DT_ const * val() const
      {
        if (this->_elements.size() == 0)
          return nullptr;

        return this->_elements.at(0);
      }

      /**
       * \brief Retrieve row start index array.
       *
       * \returns Row start index array.
       */
      IT_ * row_ptr()
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(1);
      }

      IT_ const * row_ptr() const
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(1);
      }

    }; //MatrixMirrorBuffer

  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_MATRIX_MIRROR_BUFFER_HPP
