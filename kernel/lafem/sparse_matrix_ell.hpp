#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_ELL_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_ELL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/matrix_base.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/scale_row_col.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/adjacency/graph.hpp>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief ELL based sparse matrix.
     *
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the ELL-C format.\n\n
     * Data survey: \n
     * _elements[0]: val     - raw non zero number values, stored in ELL-C storage format \n
     * _indices[0]:  col_ind - column index per non zero element, stored in same format as val \n
     * _indices[1]:  cs      - starting offset of each chunk (including matrix end index)\n
     * _indices[2]:  cl      - length of the longset row in each chunk
     * _indices[3]:  rl      - length of each row
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: chunk size (C) \n
     * _scalar_index[4]: num_of_chunks rounded up division of row-counter by chunk size \n
     * _scalar_index[5]: size of val- and row-arrays \n
     * _scalar_index[6]: non zero element count (used elements) \n
     * _scalar_dt[0]: zero element
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseMatrixELL : public Container<Mem_, DT_, IT_>, public MatrixBase
    {
    private:
      void _read_from_mtx(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        _read_from_mtx(file);
        file.close();
      }

      void _read_from_mtx(std::istream& file)
      {
        std::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >

        Index trows, tcols, ue(0);
        String line;
        std::getline(file, line);
        const bool general((line.find("%%MatrixMarket matrix coordinate real general") != String::npos) ? true : false);
        const bool symmetric((line.find("%%MatrixMarket matrix coordinate real symmetric") != String::npos) ? true : false);

        if (symmetric == false && general == false)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a compatible mtx-file");
        }

        while(!file.eof())
        {
          std::getline(file,line);
          if (file.eof())
            throw InternalError(__func__, __FILE__, __LINE__, "Input-file is empty");

          String::size_type begin(line.find_first_not_of(" "));
          if (line.at(begin) != '%')
            break;
        }
        {
          String::size_type begin(line.find_first_not_of(" "));
          line.erase(0, begin);
          String::size_type end(line.find_first_of(" "));
          String srow(line, 0, end);
          trows = Index(atol(srow.c_str()));
          line.erase(0, end);

          begin = line.find_first_not_of(" ");
          line.erase(0, begin);
          end = line.find_first_of(" ");
          String scol(line, 0, end);
          tcols = Index(atol(scol.c_str()));
          line.erase(0, end);
        }

        while(!file.eof())
        {
          std::getline(file, line);
          if (file.eof())
            break;

          String::size_type begin(line.find_first_not_of(" "));
          line.erase(0, begin);
          String::size_type end(line.find_first_of(" "));
          String srow(line, 0, end);
          IT_ row((IT_)atol(srow.c_str()));
          --row;
          line.erase(0, end);

          begin = line.find_first_not_of(" ");
          line.erase(0, begin);
          end = line.find_first_of(" ");
          String scol(line, 0, end);
          IT_ col((IT_)atol(scol.c_str()));
          --col;
          line.erase(0, end);

          begin = line.find_first_not_of(" ");
          line.erase(0, begin);
          end = line.find_first_of(" ");
          String sval(line, 0, end);
          DT_ tval((DT_)atof(sval.c_str()));

          entries[IT_(row)].insert(std::pair<IT_, DT_>(col, tval));
          ++ue;
          if (symmetric == true && row != col)
          {
            entries[IT_(col)].insert(std::pair<IT_, DT_>(row, tval));
            ++ue;
          }
        }

        const IT_ tC(IT_(this->_C()));
        const Index tnum_of_chunks((trows + tC - Index(1)) / tC);

        LAFEM::DenseVector<Mem::Main, IT_, IT_> tcl(tnum_of_chunks, IT_(0));
        IT_ * ptcl(tcl.elements());
        LAFEM::DenseVector<Mem::Main, IT_, IT_> tcs(tnum_of_chunks + 1);
        IT_ * ptcs(tcs.elements());
        LAFEM::DenseVector<Mem::Main, IT_, IT_> trl(trows);
        IT_ * ptrl(trl.elements());

        for (IT_ i(0); i < trows; ++i)
        {
          ptrl[i] = IT_(entries[i].size());
          if (ptrl[i] > ptcl[i/tC])
            ptcl[i/tC] = IT_(ptrl[i]);
        }

        ptcs[0] = IT_(0);
        for (Index i(0); i < tnum_of_chunks; ++i)
        {
          ptcs[i+1] = ptcs[i] + IT_(tC) * ptcl[i];
        }

        Index tval_size = Index(ptcs[tnum_of_chunks]);

        LAFEM::DenseVector<Mem::Main, IT_, IT_> tcol_ind(tval_size);
        IT_ * ptcol_ind(tcol_ind.elements());
        LAFEM::DenseVector<Mem::Main, DT_, IT_> tval(tval_size);
        DT_ * ptval(tval.elements());

        IT_ row_idx(0);
        for (auto row : entries)
        {
          IT_ idx(ptcs[row_idx/tC] + row_idx%tC);

          for (auto col : row.second)
          {
            ptcol_ind[idx] = col.first;
            ptval    [idx] = col.second;
            idx += tC;
          }
          for (; idx < ptcs[row_idx/tC + 1]; idx += tC)
          {
            ptcol_ind[idx] = IT_(0);
            ptval    [idx] = DT_(0);
          }
          ++row_idx;
        }

        _init_padded_elements<Mem::Main, Mem::Main>(ptval, ptcol_ind, ptcs, trows, tnum_of_chunks, tC);

        this->assign(SparseMatrixELL<Mem::Main, DT_, IT_>(trows, tcols, ue, tval, tcol_ind, tcs, tcl, trl, tC));
      }

      void _read_from_ell(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        _read_from_ell(file);
        file.close();
      }

      void _read_from_ell(std::istream& file)
      {
        this->template _deserialise<double, uint64_t>(FileMode::fm_ell, file);
      }

      Index & _size()
      {
        return this->_scalar_index.at(0);
      }
      Index & _rows()
      {
        return this->_scalar_index.at(1);
      }
      Index & _columns()
      {
        return this->_scalar_index.at(2);
      }
      Index & _C()
      {
        return this->_scalar_index.at(3);
      }
      Index & _num_of_chunks()
      {
        return this->_scalar_index.at(4);
      }
      Index & _val_size()
      {
        return this->_scalar_index.at(5);
      }
      Index & _used_elements()
      {
        return this->_scalar_index.at(6);
      }

      /// initialise padded elements of val- and col_ind-array if rows modulo C != 0
      template <typename Mem2_, typename Mem3_>
      void _init_padded_elements(DT_ * const pval, IT_ * const pcol_ind, const IT_ * const pcs,
                                 const Index trows, const Index tnum_of_chunks, const Index tC)
      {
        if (trows % tC != 0)
        {
          for (Index i(pcs[tnum_of_chunks-1] + trows - (tnum_of_chunks-1) * tC); i < pcs[tnum_of_chunks]; i+=tC)
          {
            Util::MemoryPool<Mem2_>::instance()->set_memory(pval     + i, DT_(0), tnum_of_chunks * tC - trows);
            Util::MemoryPool<Mem3_>::instance()->set_memory(pcol_ind + i, IT_(0), tnum_of_chunks * tC - trows);
          }
        }
      }

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indexype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Compatible L-vector type
      typedef DenseVector<Mem_, DT_, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<Mem_, DT_, IT_> VectorTypeR;
      /// Our used layout type
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_ell;
      /// ImageIterator class for Adjactor interface implementation
      class ImageIterator
      {
      private:
        const IT_* _ai;
        IT_ _c;

      public:
        ImageIterator() : _ai(nullptr), _c(IT_(0)) {}

        ImageIterator(const IT_* ai, IT_ c) : _ai(ai), _c(c) {}

        ImageIterator& operator=(const ImageIterator& other)
        {
          _ai = other._ai;
          _c = other._c;
          return *this;
        }

        bool operator!=(const ImageIterator& other) const
        {
          return this->_ai != other._ai;
        }

        ImageIterator& operator++()
        {
          _ai += _c;
          return *this;
        }

        Index operator*() const
        {
          return Index(*_ai);
        }
      };
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = class SparseMatrixELL<Mem2_, DT2_, IT2_>;

      /**
       * \brief Constructor
       *
       * \param[in] C chunk size (default = 32).
       *
       * Creates an empty non dimensional matrix.
       */
      explicit SparseMatrixELL(const Index C_in = 32) :
        Container<Mem_, DT_, IT_> (0)
      {
        CONTEXT("When creating SparseMatrixELL");
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(C_in);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] layout The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixELL(const SparseLayout<Mem_, IT_, layout_id> & layout_in) :
        Container<Mem_, DT_, IT_> (layout_in._scalar_index.at(0))
      {
        CONTEXT("When creating SparseMatrixELL");
        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          Util::MemoryPool<Mem_>::instance()->increase_memory(i);

        this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_val_size()));
        this->_elements_size.push_back(_val_size());

        Util::MemoryPool<Mem_>::instance()->set_memory(val(), DT_(0), _val_size());
      }

      /**
       * \brief Constructor
       *
       * \param[in] other The source matrix.
       * \param[in] C chunk size (default = 32).
       *
       * Creates a ELL matrix based on the source matrix.
       */
      template <typename MT_>
      explicit SparseMatrixELL(const MT_ & other, const Index C_in = 32) :
        SparseMatrixELL(C_in)
      {
        CONTEXT("When creating SparseMatrixELL");

        convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows The row count of the created matrix.
       * \param[in] columns The column count of the created matrix.
       * \param[in] used_elements number of non zero elements.
       * \param[in] val Vector with data values.
       * \param[in] col_ind Vector with column indices.
       * \param[in] cs starting-offset of each chunk.
       * \param[in] cl length of the longest row in each chunk.
       * \param[in] rl length of each row.
       * \param[in] C chunk size (default = 32).
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixELL(const Index rows_in, const Index columns_in,
                               const Index used_elements_in,
                               DenseVector<Mem_, DT_, IT_> & val_in,
                               DenseVector<Mem_, IT_, IT_> & col_ind_in,
                               DenseVector<Mem_, IT_, IT_> & cs_in,
                               DenseVector<Mem_, IT_, IT_> & cl_in,
                               DenseVector<Mem_, IT_, IT_> & rl_in,
                               const Index C_in = 32) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        CONTEXT("When creating SparseMatrixELL");

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(C_in);
        this->_scalar_index.push_back((rows_in + C_in - Index(1)) / C_in);
        this->_scalar_index.push_back(val_in.size());
        this->_scalar_index.push_back(used_elements_in);
        this->_scalar_dt.push_back(DT_(0));

        ASSERT(val_in.size() == col_ind_in.size(), "Error: val- and col-arrays must have the same size!");
        ASSERT(cs_in.size() == num_of_chunks() + 1, "Error: cs-array-size must match to row-count and chunk size!");
        ASSERT(cl_in.size() == num_of_chunks(), "Error: cl-array-size must match to row-count and chunk size!");
        ASSERT(rl_in.size() == rows(), "Error: rl-array-size must match to row-count!");

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(col_ind_in.elements());
        this->_indices_size.push_back(col_ind_in.size());
        this->_indices.push_back(cs_in.elements());
        this->_indices_size.push_back(cs_in.size());
        this->_indices.push_back(cl_in.elements());
        this->_indices_size.push_back(cl_in.size());
        this->_indices.push_back(rl_in.elements());
        this->_indices_size.push_back(rl_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix from
       * \param[in] C chunk size (default = 32)
       *
       * Creates a ELL matrix based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit SparseMatrixELL(const Adjacency::Graph & graph, const Index C_in = 32) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixELL");

        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nnze = graph.get_num_indices();
        Index num_of_chunks_in = (num_rows + C_in - Index(1)) / C_in;

        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());

        // create temporary vector
        LAFEM::DenseVector<Mem::Main, IT_, IT_> tcl(num_of_chunks_in, IT_(0));
        IT_ * ptcl(tcl.elements());
        LAFEM::DenseVector<Mem::Main, IT_, IT_> tcs(num_of_chunks_in + 1);
        IT_ * ptcs(tcs.elements());
        LAFEM::DenseVector<Mem::Main, IT_, IT_> trl(num_rows);
        IT_ * ptrl(trl.elements());

        for (Index i(0); i < num_rows; ++i)
        {
          ptrl[i] = IT_(dom_ptr[i+1] - dom_ptr[i]);
          if (dom_ptr[i+1] - dom_ptr[i] > ptcl[i/C_in])
            ptcl[i/C_in] = IT_(dom_ptr[i+1] - dom_ptr[i]);
        }

        ptcs[0] = IT_(0);
        for (Index i(0); i < num_of_chunks_in; ++i)
        {
          ptcs[i+1] = ptcs[i] + IT_(C_in) * ptcl[i];
        }

        Index val_size_in = Index(ptcs[num_of_chunks_in]);

        LAFEM::DenseVector<Mem::Main, IT_, IT_> tcol_ind(val_size_in);
        IT_ * ptcol_ind(tcol_ind.elements());

        for (Index i(0); i < num_rows; ++i)
        {
          for (Index j(dom_ptr[i]), k(0); j < dom_ptr[i+1]; ++j, ++k)
          {
            ptcol_ind[ptcs[i/C_in] + i%C_in + k*C_in] = IT_(img_idx[j]);
          }
          for (Index k(ptrl[i]); k < ptcl[i/C_in]; ++k)
          {
            ptcol_ind[ptcs[i/C_in] + i%C_in + k*C_in] = IT_(0);
          }
        }

        LAFEM::DenseVector<Mem_, DT_, IT_> tval_mem(val_size_in);
        auto * ptval_mem(tval_mem.elements());

        _init_padded_elements<Mem_, Mem::Main>(ptval_mem, ptcol_ind, ptcs, num_rows, num_of_chunks_in, C_in);

        LAFEM::DenseVector<Mem_, IT_, IT_> tcol_ind_mem;
        tcol_ind_mem.convert(tcol_ind);
        LAFEM::DenseVector<Mem_, IT_, IT_> tcs_mem;
        tcs_mem.convert(tcs);
        LAFEM::DenseVector<Mem_, IT_, IT_> tcl_mem;
        tcl_mem.convert(tcl);
        LAFEM::DenseVector<Mem_, IT_, IT_> trl_mem;
        trl_mem.convert(trl);

        // build the matrix
        this->assign(SparseMatrixELL<Mem_, DT_, IT_>(num_rows, num_cols, num_nnze,
                                                     tval_mem, tcol_ind_mem, tcs_mem, tcl_mem, trl_mem, C_in));
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       * \param[in] C chunk size (default = 32).
       *
       * Creates a ELL matrix based on the source file.
       */
      explicit SparseMatrixELL(FileMode mode, String filename, const Index C_in = 32) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixELL");

        switch(mode)
        {
          case FileMode::fm_mtx:
            this->assign(SparseMatrixELL(C_in));
            _read_from_mtx(filename);
            break;
          case FileMode::fm_ell:
            _read_from_ell(filename);
            break;
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       * \param[in] C chunk size (default = 32).
       *
       * Creates a ELL matrix based on the source filestream.
       */
      explicit SparseMatrixELL(FileMode mode, std::istream& file, const Index C_in = 32) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixELL");

        switch(mode)
        {
          case FileMode::fm_mtx:
            this->assign(SparseMatrixELL(C_in));
            _read_from_mtx(file);
            break;
          case FileMode::fm_ell:
            _read_from_ell(file);
            break;
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "Filemode nggot supported!");
        }
      }

      /**
       * \brief Constructor
       *
       * \param[in] std::vector<char> A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixELL(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixELL");
        deserialise<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixELL(SparseMatrixELL && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixELL>(other))
      {
        CONTEXT("When moving SparseMatrixELL");
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixELL & operator= (SparseMatrixELL && other)
      {
        CONTEXT("When moving SparseMatrixELL");

        this->move(std::forward<SparseMatrixELL>(other));

        return *this;
      }

      InsertWeakClone( SparseMatrixELL );

      /** \brief Shallow copy operation
       *
       * Create a shallow copy of itself.
       *
       */
      SparseMatrixELL shared() const
      {
        CONTEXT("When sharing SparseMatrixELL");
        SparseMatrixELL r;
        r.assign(*this);
        return r;
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixELL<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting SparseMatrixELL");
        this->assign(other);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixBanded<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting SparseMatrixELL");

        const Index tC(this->_C());
        this->clear();

        this->_scalar_index.push_back(other.size());
        this->_scalar_index.push_back(other.rows());
        this->_scalar_index.push_back(other.columns());
        this->_scalar_index.push_back(tC);
        this->_scalar_index.push_back((_rows() + _C() - Index(1)) / _C());
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(other.used_elements());
        this->_scalar_dt.push_back(other.zero_element());

        SparseMatrixBanded<Mem::Main, DT_, IT_> cother;
        cother.convert(other);

        IT_ * tcl = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_num_of_chunks());
        Util::MemoryPool<Mem::Main>::instance()->set_memory(tcl, IT_(0), _num_of_chunks());
        IT_ * trl = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_rows());

        const DT_ * cval(cother.val());
        const IT_ * coffsets(cother.offsets());
        const Index cnum_of_offsets(cother.num_of_offsets());
        const Index crows(cother.rows());

        // Search first offset of the upper triangular matrix
        Index k(0);
        while (k < cnum_of_offsets && coffsets[k] + 1 < crows)
        {
          ++k;
        }

        for (Index i(k + 1); i > 0;)
        {
          --i;

          // iteration over all offsets of the upper triangular matrix
          for (Index j(cnum_of_offsets + 1); j > 0;)
          {
            --j;

            const Index start(Math::max(cother.start_offset(i),
                                        cother.end_offset(j) + 1));
            const Index end  (Math::min(cother.start_offset(i-1),
                                        cother.end_offset(j-1) + 1));
            if (start < end)
            {
              for (Index l(start/_C()); l <= (end-1)/_C(); ++l)
              {
                if (j-i > tcl[l])
                  tcl[l] = IT_(j - i);
              }
              for (Index l(start); l < end; ++l)
              {
                trl[l] = IT_(j - i);
              }
            }
          }
        }

        IT_ * tcs = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_num_of_chunks() + 1);

        tcs[0] = IT_(0);
        for (Index i(0); i < _num_of_chunks(); ++i)
        {
          tcs[i+1] = tcs[i] + IT_(_C()) * tcl[i];
        }

        _val_size() = tcs[_num_of_chunks()];

        DT_ * tval = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(_val_size());
        IT_ * tcol_ind = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_val_size());

        for (Index i(k + 1); i > 0;)
        {
          --i;

          // iteration over all offsets of the upper triangular matrix
          for (Index j(cnum_of_offsets + 1); j > 0;)
          {
            --j;

            // iteration over all rows which contain the offsets between offset i and offset j
            const Index start(Math::max(cother.start_offset(i),
                                        cother.end_offset(j) + 1));
            const Index end  (Math::min(cother.start_offset(i-1),
                                        cother.end_offset(j-1) + 1));
            for (Index l(start); l < end; ++l)
            {
              for (Index a(i); a < j; ++a)
              {
                tcol_ind[tcs[l/_C()] + l%_C() + (a - i) * _C()] = IT_(l + coffsets[a] + 1 - crows);
                tval[tcs[l/_C()] + l%_C() + (a - i) * _C()] = cval[a * crows + l];
              }
              for (Index a(j-i); a < tcl[l/_C()]; ++a)
              {
                tcol_ind[tcs[l/_C()] + l%_C() + a * _C()] = IT_(0);
                tval[tcs[l/_C()] + l%_C() + a * _C()] = DT_(0);
              }
            }
          }
        }

        _init_padded_elements<Mem::Main, Mem::Main>(tval, tcol_ind, tcs, _rows(), _num_of_chunks(), tC);

        this->_elements_size.push_back(_val_size());
        this->_indices_size.push_back(_val_size());
        this->_indices_size.push_back(_num_of_chunks() + 1);
        this->_indices_size.push_back(_num_of_chunks());
        this->_indices_size.push_back(_rows());
        if (std::is_same<Mem_, Mem::Main>::value)
        {
          this->_elements.push_back(tval);
          this->_indices.push_back(tcol_ind);
          this->_indices.push_back(tcs);
          this->_indices.push_back(tcl);
          this->_indices.push_back(trl);
        }
        else
        {
          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_val_size()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_val_size()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_num_of_chunks() + 1));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_num_of_chunks()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_rows()));

          Util::MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval, _val_size());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(0), tcol_ind, _val_size());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(1), tcs, _num_of_chunks() + 1);
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(2), tcl, _num_of_chunks());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(3), trl, _rows());
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tval);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcol_ind);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcs);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcl);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(trl);
        }
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixCSR<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting SparseMatrixELL");

        const Index tC(this->_C());
        this->clear();

        this->_scalar_index.push_back(other.size());
        this->_scalar_index.push_back(other.rows());
        this->_scalar_index.push_back(other.columns());
        this->_scalar_index.push_back(tC);
        this->_scalar_index.push_back((_rows() + _C() - Index(1)) / _C());
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(other.used_elements());
        this->_scalar_dt.push_back(other.zero_element());

        SparseMatrixCSR<Mem::Main, DT_, IT_> cother;
        cother.convert(other);

        IT_ * tcl = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_num_of_chunks());
        Util::MemoryPool<Mem::Main>::instance()->set_memory(tcl, IT_(0), _num_of_chunks());
        IT_ * trl = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_rows());

        const IT_ * const crow_ptr(cother.row_ptr());
        const IT_ * const ccol_ind(cother.col_ind());
        const DT_ * const cval(cother.val());

        for (Index i(0); i < _rows(); ++i)
        {
          trl[i] = crow_ptr[i+1] - crow_ptr[i];
          if (crow_ptr[i+1] - crow_ptr[i] > tcl[i/_C()])
            tcl[i/_C()] = crow_ptr[i+1] - crow_ptr[i];
        }

        IT_ * tcs = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_num_of_chunks() + 1);

        tcs[0] = IT_(0);
        for (Index i(0); i < _num_of_chunks(); ++i)
        {
          tcs[i+1] = tcs[i] + IT_(_C()) * tcl[i];
        }

        _val_size() = tcs[_num_of_chunks()];

        DT_ * tval = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(_val_size());
        IT_ * tcol_ind = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_val_size());

        for (Index i(0); i < _rows(); ++i)
        {
          for (Index j(crow_ptr[i]), k(0); j < crow_ptr[i+1]; ++j, ++k)
          {
            tcol_ind[tcs[i/_C()] + i%_C() + k*_C()] = ccol_ind[j];
            tval    [tcs[i/_C()] + i%_C() + k*_C()] = cval[j];
          }
          for (Index k(trl[i]); k < tcl[i/_C()]; ++k)
          {
            tcol_ind[tcs[i/_C()] + i%_C() + k*_C()] = IT_(0);
            tval    [tcs[i/_C()] + i%_C() + k*_C()] = DT_(0);
          }
        }

        _init_padded_elements<Mem::Main, Mem::Main>(tval, tcol_ind, tcs, _rows(), _num_of_chunks(), tC);

        this->_elements_size.push_back(_val_size());
        this->_indices_size.push_back(_val_size());
        this->_indices_size.push_back(_num_of_chunks() + 1);
        this->_indices_size.push_back(_num_of_chunks());
        this->_indices_size.push_back(_rows());
        if (std::is_same<Mem_, Mem::Main>::value)
        {
          this->_elements.push_back(tval);
          this->_indices.push_back(tcol_ind);
          this->_indices.push_back(tcs);
          this->_indices.push_back(tcl);
          this->_indices.push_back(trl);
        }
        else
        {
          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_val_size()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_val_size()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_num_of_chunks() + 1));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_num_of_chunks()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_rows()));

          Util::MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval, _val_size());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(0), tcol_ind, _val_size());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(1), tcs, _num_of_chunks() + 1);
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(2), tcl, _num_of_chunks());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(3), trl, _rows());
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tval);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcol_ind);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcs);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcl);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(trl);
        }
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixCOO<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting SparseMatrixELL");

        const Index tC(this->_C());
        this->clear();

        this->_scalar_index.push_back(other.size());
        this->_scalar_index.push_back(other.rows());
        this->_scalar_index.push_back(other.columns());
        this->_scalar_index.push_back(tC);
        this->_scalar_index.push_back((_rows() + _C() - Index(1)) / _C());
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(other.used_elements());
        this->_scalar_dt.push_back(other.zero_element());

        SparseMatrixCOO<Mem::Main, DT_, IT_> cother;
        cother.convert(other);

        IT_ * tcl = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_num_of_chunks());
        Util::MemoryPool<Mem::Main>::instance()->set_memory(tcl, IT_(0), _num_of_chunks());
        IT_ * trl = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_rows());
        Util::MemoryPool<Mem::Main>::instance()->set_memory(trl, IT_(0), _rows());

        const IT_ * const crow_ind(cother.row_indices());
        const IT_ * const ccol_ind(cother.column_indices());
        const DT_ * const cval(cother.val());

        {
          Index i(0);
          while (i < _used_elements())
          {
            const Index row(crow_ind[i]);
            IT_ counter(0);
            while (i < _used_elements() && crow_ind[i] == row)
            {
              ++counter;
              ++i;
            }
            trl[row] = counter;
            if (counter > tcl[row/_C()])
              tcl[row/_C()] = counter;
          }
        }

        IT_ * tcs = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_num_of_chunks() + 1);

        tcs[0] = IT_(0);
        for (Index i(0); i < _num_of_chunks(); ++i)
        {
          tcs[i+1] = tcs[i] + IT_(_C()) * tcl[i];
        }

        _val_size() = tcs[_num_of_chunks()];

        DT_ * tval = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(_val_size());
        IT_ * tcol_ind = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(_val_size());

        for (Index row(0), k(0); row < rows(); ++row)
        {
          for (Index i(0); i < trl[row]; ++i, ++k)
          {
            tcol_ind[tcs[row/_C()] + row%_C() + i*_C()] = ccol_ind[k];
            tval    [tcs[row/_C()] + row%_C() + i*_C()] = cval[k];
          }
          for (Index i(trl[row]); i < tcl[row/_C()]; ++i)
          {
            tcol_ind[tcs[row/_C()] + row%_C() + i*_C()] = IT_(0);
            tval    [tcs[row/_C()] + row%_C() + i*_C()] = DT_(0);
          }
        }

        _init_padded_elements<Mem::Main, Mem::Main>(tval, tcol_ind, tcs, _rows(), _num_of_chunks(), tC);

        this->_elements_size.push_back(_val_size());
        this->_indices_size.push_back(_val_size());
        this->_indices_size.push_back(_num_of_chunks() + 1);
        this->_indices_size.push_back(_num_of_chunks());
        this->_indices_size.push_back(_rows());
        if (std::is_same<Mem_, Mem::Main>::value)
        {
          this->_elements.push_back(tval);
          this->_indices.push_back(tcol_ind);
          this->_indices.push_back(tcs);
          this->_indices.push_back(tcl);
          this->_indices.push_back(trl);
        }
        else
        {
          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_val_size()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_val_size()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_num_of_chunks() + 1));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_num_of_chunks()));
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_rows()));

          Util::MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval, _val_size());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(0), tcol_ind, _val_size());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(1), tcs, _num_of_chunks() + 1);
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(2), tcl, _num_of_chunks());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(3), trl, _rows());
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tval);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcol_ind);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcs);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcl);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(trl);
        }
      }

      /**
       * \brief Conversion method
       *
       * \param[in] a The input matrix.
       *
       * Converts any matrix to SparseMatrixELL-format
       */
      template <typename MT_>
      void convert(const MT_ & a)
      {
        CONTEXT("When converting SparseMatrixELL");
        std::cerr<<"Warning: Generic matrix convert used!"<<std::endl;

        typename MT_::template ContainerType<Mem::Main, DT_, IT_> ta;
        ta.convert(a);

        const Index aC(this->C());
        const Index arows(ta.rows());
        const Index acolumns(ta.columns());
        const Index aused_elements(ta.used_elements());
        const Index anum_of_chunks((arows + aC - Index(1)) / aC);

        DenseVector<Mem::Main, IT_, IT_> arl(arows);
        IT_ * prl(arl.elements());
        DenseVector<Mem::Main, IT_, IT_> acl(anum_of_chunks, IT_(0));
        IT_ * pcl(acl.elements());
        DenseVector<Mem::Main, IT_, IT_> acs(anum_of_chunks + 1);
        IT_ * pcs(acs.elements());

        for (Index i(0); i < arows; ++i)
        {
          prl[i] = IT_(ta.get_length_of_line(i));
          if (prl[i] > pcl[i/_C()])
            pcl[i/_C()] = prl[i];
        }

        pcs[0] = IT_(0);
        for (Index i(0); i < anum_of_chunks; ++i)
        {
          pcs[i+1] = pcs[i] + IT_(aC) * pcl[i];
        }

        const Index aval_size(pcs[anum_of_chunks]);

        DenseVector<Mem::Main, DT_, IT_> aval(aval_size);
        DenseVector<Mem::Main, IT_, IT_> acol_ind(aval_size);

        DT_ * pval(aval.elements());
        IT_ * pcol_ind(acol_ind.elements());

        for (Index i(0); i < arows; ++i)
        {
          ta.set_line(i, pval + pcs[i/aC] + i%aC, pcol_ind + pcs[i/aC] + i%aC, 0, aC);

          for (Index k(prl[i]); k < pcl[i/_C()]; ++k)
          {
            pcol_ind[pcs[i/aC] + i%aC + k*aC] = IT_(0);
            pval    [pcs[i/aC] + i%aC + k*aC] = DT_(0);
          }
        }

        _init_padded_elements<Mem::Main, Mem::Main>(pval, pcol_ind, pcs, arows, anum_of_chunks, aC);

        SparseMatrixELL<Mem::Main, DT_, IT_> ta_ell(arows, acolumns, aused_elements, aval, acol_ind, acs, acl, arl, aC);
        this->assign(ta_ell);
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixELL & operator= (const SparseLayout<Mem_, IT_, layout_id> & layout_in)
      {
        CONTEXT("When assigning SparseMatrixELL");

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

        this->_elements.clear();
        this->_indices.clear();
        this->_elements_size.clear();
        this->_indices_size.clear();
        this->_scalar_index.clear();
        this->_scalar_dt.clear();

        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          Util::MemoryPool<Mem_>::instance()->increase_memory(i);

        this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_val_size()));
        this->_elements_size.push_back(_val_size());

        return *this;
      }

      /**
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] std::vector<char> A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialise(std::vector<char> input)
      {
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_ell, input);
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * \param[in] mode FileMode enum, describing the actual container specialisation.
       * \param[out] std::vector<char> A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAST::LAFEM::Container::_serialise for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialise()
      {
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_ell);
      }

      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        CONTEXT("When writing out SparseMatrixELL");

        switch(mode)
        {
          case FileMode::fm_ell:
            write_out_ell(filename);
            break;
          case FileMode::fm_mtx:
            write_out_mtx(filename);
            break;
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be written to.
       */
      void write_out(FileMode mode, std::ostream& file) const
      {
        CONTEXT("When writing out SparseMatrixELL");

        switch(mode)
        {
          case FileMode::fm_ell:
            write_out_ell(file);
            break;
          case FileMode::fm_mtx:
           write_out_mtx(file);
           break;
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Write out matrix to ell binary file.
       *
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out_ell(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        write_out_ell(file);
        file.close();
      }

      /**
       * \brief Write out matrix to ell binary file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_ell(std::ostream& file) const
      {
        if (! std::is_same<DT_, double>::value)
          std::cout<<"Warning: You are writing out an ell matrix that is not double precision!"<<std::endl;

        this->template _serialise<double, uint64_t>(FileMode::fm_ell, file);
      }

      /**
       * \brief Write out matrix to MatrixMarktet mtx file.
       *
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out_mtx(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        write_out_mtx(file);
        file.close();
      }

      /**
       * \brief Write out matrix to MatrixMarktet mtx file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_mtx(std::ostream& file) const
      {
        SparseMatrixELL<Mem::Main, DT_, IT_> temp;
        temp.convert(*this);

        file << "%%MatrixMarket matrix coordinate real general" << std::endl;
        file << temp.rows() << " " << temp.columns() << " " << temp.used_elements() << std::endl;

        const IT_ * const ptcs(temp.cs());
        const IT_ * const ptrl(temp.rl());
        const IT_ * const ptcol_ind(temp.col_ind());
        const DT_ * const ptval(temp.val());
        const IT_ tC(IT_(temp.C()));


        for (IT_ row(0) ; row < rows() ; ++row)
        {
          IT_ end(ptrl[row]);
          for (IT_ i(0) ; i < end ; ++i)
          {
            const IT_ idx(ptcs[row/tC] + row%tC + tC*i);
            file << row + 1 << " " << ptcol_ind[idx] + 1 << " " << std::scientific << ptval[idx] << std::endl;
          }
        }
      }

      /**
       * \brief Retrieve specific matrix element.
       *
       * \param[in] row The row of the matrix element.
       * \param[in] col The column of the matrix element.
       *
       * \returns Specific matrix element.
       */
      DT_ operator()(Index row, Index col) const
      {
        CONTEXT("When retrieving SparseMatrixELL element");

        ASSERT(row < rows(), "Error: " + stringify(row) + " exceeds sparse matrix ell row size " + stringify(rows()) + " !");
        ASSERT(col < columns(), "Error: " + stringify(col) + " exceeds sparse matrix ell column size " + stringify(columns()) + " !");

        const Index nchunk(Index(floor(row / float(C()))));
        Index start(Index(Util::MemoryPool<Mem_>::get_element(cs(), nchunk)));
        Index max(Index(Util::MemoryPool<Mem_>::get_element(rl(), row)));
        for (Index i(start + row - nchunk * C()), j(0) ; j < max && Index(Util::MemoryPool<Mem_>::get_element(col_ind(), i)) <= col ; i += C(), ++j)
        {
          if (Index(Util::MemoryPool<Mem_>::get_element(col_ind(), i)) == col)
            return Util::MemoryPool<Mem_>::get_element(val(), i);
        }
        return zero_element();
      }

      /**
       * \brief Retrieve convenient sparse matrix layout object.
       *
       * \return An object containing the sparse matrix layout.
       */
      SparseLayout<Mem_, IT_, layout_id> layout() const
      {
        return SparseLayout<Mem_, IT_, layout_id>(this->_indices, this->_indices_size, this->_scalar_index);
      }

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count.
       */
      const Index & rows() const
      {
        return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      const Index & columns() const
      {
        return this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve chunk size.
       *
       * \returns Chunk size.
       */
      const Index & C() const
      {
        return this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve number of chunks.
       *
       * \returns number of chunks.
       */
      const Index & num_of_chunks() const
      {
        return this->_scalar_index.at(4);
      }

      /**
       * \brief Retrieve size of val- and col-arrays
       *
       * \returns Size of val- and col-arrays
       */
      const Index & val_size() const
      {
        return this->_scalar_index.at(5);
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      const Index & used_elements() const override
      {
        return this->_scalar_index.at(6);
      }

      /**
       * \brief Retrieve column indices array.
       *
       * \returns Column indices array.
       */
      IT_ const * col_ind() const
      {
        return this->_indices.at(0);
      }

      /**
       * \brief Retrieve array with starting-offsets of each chunk.
       *
       * \returns array with starting-offsets of each chunk.
       */
      IT_ const * cs() const
      {
        return this->_indices.at(1);
      }

      /**
       * \brief Retrieve array with lengths of the longest row in each chunk.
       *
       * \returns array with lengths of the longest row in each chunk.
       */
      IT_ const * cl() const
      {
        return this->_indices.at(2);
      }

      /**
       * \brief Retrieve array with lengths of each row.
       *
       * \returns array with lengths of each row.
       */
      IT_ const * rl() const
      {
        return this->_indices.at(3);
      }

      /**
       * \brief Retrieve non zero element array.
       *
       * \returns Non zero element array.
       */
      DT_ * val()
      {
        return this->_elements.at(0);
      }

      DT_ const * val() const
      {
        return this->_elements.at(0);
      }

      /**
       * \brief Retrieve non zero element.
       *
       * \returns Non zero element.
       */
      DT_ zero_element() const
      {
        return this->_scalar_dt.at(0);
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SparseMatrixELL";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      void copy(const SparseMatrixELL & x)
      {
        this->_copy_content(x);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      template <typename Mem2_>
      void copy(const SparseMatrixELL<Mem2_, DT_, IT_> & x)
      {
        this->_copy_content(x);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow y + \alpha x\f$
       *
       * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
       *
       * \param[in] x The first summand matrix to be scaled.
       * \param[in] y The second summand matrix
       * \param[in] alpha A scalar to multiply x with.
       */
      template <typename Algo_>
      void axpy(
                const SparseMatrixELL & x,
                const SparseMatrixELL & y,
                DT_ alpha = DT_(1))
      {
        if (x.rows() != y.rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.used_elements() != y.used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");
        if (x.C() != y.C())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix chunk size do not match!");
        if (x.C() != this->C())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix chunk size do not match!");

        // check for special cases
        // r <- x + y
        if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
          Arch::Sum<Mem_, Algo_>::value(this->val(), x.val(), y.val(), this->val_size());
        // r <- y - x
        else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          Arch::Difference<Mem_, Algo_>::value(this->val(), y.val(), x.val(), this->val_size());
        // r <- y
        else if (Math::abs(alpha) < Math::eps<DT_>())
          this->copy(y);
        // r <- y + alpha*x
        else
          Arch::Axpy<Mem_, Algo_>::dv(this->val(), alpha, x.val(), y.val(), this->val_size());
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha x\f$
       *
       * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      template <typename Algo_>
      void scale(const SparseMatrixELL & x, const DT_ alpha)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
        if (x.C() != this->C())
          throw InternalError(__func__, __FILE__, __LINE__, "Chunk size does not match!");

        Arch::Scale<Mem_, Algo_>::value(this->val(), x.val(), alpha, this->val_size());
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      template <typename Algo_>
      DT_ norm_frobenius() const
      {
        return Arch::Norm2<Mem_, Algo_>::value(this->val(), this->val_size());
      }

      /**
       * \brief Calculate \f$this^\top \f$
       *
       * \return The transposed matrix
       */
      SparseMatrixELL transpose() const
      {
        SparseMatrixELL x_t;
        x_t.transpose(*this);
        return x_t;
      }

      /**
       * \brief Calculate \f$this \leftarrow x^\top \f$
       *
       * \param[in] x The matrix to be transposed.
       */
      void transpose(const SparseMatrixELL & x)
      {
        SparseMatrixELL<Mem::Main, DT_, IT_> tx;
        tx.convert(x);

        const Index txrows(tx.rows());
        const Index txcolumns(tx.columns());
        const Index txused_elements(tx.used_elements());
        const Index txC(tx.C());

        const DT_ * ptxval(tx.val());
        const IT_ * ptxcol_ind(tx.col_ind());
        const IT_ * ptxcs(tx.cs());
        const IT_ * ptxrl(tx.rl());

        const Index tC(this->C());
        const Index tnum_of_chunks((txcolumns + tC - Index(1)) / tC);

        DenseVector<Mem::Main, IT_, IT_> trl(txcolumns, IT_(0));
        IT_ * ptrl(trl.elements());

        for (Index i(0); i < txrows; ++i)
        {
          for (Index k(0), j(ptxcs[i/txC] + i%txC); k < ptxrl[i]; ++k, j+= txC)
          {
            ++ptrl[ptxcol_ind[j]];
          }
        }

        DenseVector<Mem::Main, IT_, IT_> tcl(tnum_of_chunks, IT_(0));
        IT_ * ptcl(tcl.elements());
        DenseVector<Mem::Main, IT_, IT_> tcs(tnum_of_chunks + 1);
        IT_ * ptcs(tcs.elements());

        for (Index i(0); i < txcolumns; ++i)
        {
          if (ptcl[i/tC] < ptrl[i])
          {
            ptcl[i/tC] = ptrl[i];
          }
          ptrl[i] = IT_(0);
        }

        ptcs[0] = IT_(0);
        for (Index i(0); i < tnum_of_chunks; ++i)
        {
          ptcs[i+1] = ptcs[i] + IT_(tC) * ptcl[i];
        }

        const Index tval_size(ptcs[tnum_of_chunks]);

        DenseVector<Mem::Main, IT_, IT_> tcol_ind(tval_size);
        DenseVector<Mem::Main, DT_, IT_> tval(tval_size);

        IT_ * ptcol_ind(tcol_ind.elements());
        DT_ * ptval(tval.elements());

        for (Index i(0); i < txrows; ++i)
        {
          for (Index k(0), j(ptxcs[i/txC] + i%txC); k < ptxrl[i]; ++k, j+= txC)
          {
            const IT_ row(ptxcol_ind[j]);
            ptcol_ind[ptcs[row/tC] + row%tC + ptrl[row] * tC] = IT_(i);
            ptval    [ptcs[row/tC] + row%tC + ptrl[row] * tC] = ptxval[j];
            ++ptrl[row];
          }
        }

        SparseMatrixELL<Mem::Main, DT_, IT_> t(txcolumns, txrows, txused_elements, tval, tcol_ind, tcs, tcl, trl, tC);
        this->assign(t);
      }

      /**
       * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_i\f$
       *
       * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
       *
       * \param[in] x The matrix whose rows are to be scaled.
       * \param[in] s The vector to the scale the rows by.
       */
      template<typename Algo_>
      void scale_rows(const SparseMatrixELL & x, const DenseVector<Mem_,DT_,IT_> & s)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
        if (s.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        Arch::ScaleRows<Mem_, Algo_>::ell(this->val(), x.val(), this->col_ind(), this->cs(),
                                          this->cl(), this->rl(), s.elements(), this->C(), rows());
      }

      /**
       * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_j\f$
       *
       * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
       *
       * \param[in] x The matrix whose columns are to be scaled.
       * \param[in] s The vector to the scale the columns by.
       */
      template<typename Algo_>
      void scale_cols(const SparseMatrixELL & x, const DenseVector<Mem_,DT_,IT_> & s)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
        if (s.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        Arch::ScaleCols<Mem_, Algo_>::ell(this->val(), x.val(), this->col_ind(), this->cs(),
                                          this->cl(), this->rl(), s.elements(), this->C(), rows());
      }


      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      template<typename Algo_>
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        Arch::ProductMatVec<Mem_, Algo_>::ell(r.elements(), this->val(), this->col_ind(), this->cs(), this->cl(),
                                              x.elements(), this->C(), this->rows());
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$, global version.
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] gate The gate base pointer
       */
      template<typename Algo_>
      void apply(DenseVector<Mem_,DT_, IT_>& r,
        const DenseVector<Mem_, DT_, IT_>& x,
        Arch::ProductMat0Vec1GatewayBase<Mem_, Algo_, DenseVector<Mem_, DT_, IT_>, SparseMatrixELL<Mem_, DT_, IT_> >* gate)
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        gate->value(r, *this, x);
      }

      /**
       * \brief Calculate \f$r \leftarrow y + \alpha this\cdot x \f$
       *
       * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      template<typename Algo_>
      void apply(
                 DenseVector<Mem_,DT_, IT_>& r,
                 const DenseVector<Mem_, DT_, IT_>& x,
                 const DenseVector<Mem_, DT_, IT_>& y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Arch::Defect<Mem_, Algo_>::ell(r.elements(), y.elements(), this->val(), this->col_ind(),
                                         this->cs(), this->cl(), x.elements(), this->C(), this->rows());
        }
        // r <- y
        else if (Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_, Algo_>::ell(r.elements(), alpha, x.elements(), y.elements(), this->val(),
                                       this->col_ind(), this->cs(), this->cl(), this->C(), this->rows());
        }
      }
      ///@}

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(this->rows());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->columns());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        return this->rl()[row];
      }

      /// \cond internal

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DT_ * const pval_set, IT_ * const pcol_set,
                    const Index col_start, const Index stride_in = 1) const
      {
        const Index tC(this->C());
        const Index length(this->rl()[row]);
        const IT_ * pcol_ind(this->col_ind() + this->cs()[row/tC] + row%tC);
        const DT_ * pval    (this->val()     + this->cs()[row/tC] + row%tC);

        for (Index i(0); i < length; ++i)
        {
          pval_set[i * stride_in] =     pval[i * tC];
          pcol_set[i * stride_in] = pcol_ind[i * tC] + IT_(col_start);
        }
      }
      /// \endcond

      /* ******************************************************************* */
      /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
      /* ******************************************************************* */
    public:
      /** \copydoc Adjactor::get_num_nodes_domain() */
      const inline Index & get_num_nodes_domain() const
      {
        return rows();
      }

      /** \copydoc Adjactor::get_num_nodes_image() */
      const inline Index & get_num_nodes_image() const
      {
        return columns();
      }

      /** \copydoc Adjactor::image_begin() */
      inline ImageIterator image_begin(Index domain_node) const
      {
        ASSERT(domain_node < rows(), "Domain node index out of range");
        return ImageIterator(&col_ind()[cs()[domain_node/C()] + domain_node%C()], C());
      }

      /** \copydoc Adjactor::image_end() */
      inline ImageIterator image_end(Index domain_node) const
      {
        CONTEXT("Graph::image_end()");
        ASSERT(domain_node < rows(), "Domain node index out of range");
        return ImageIterator(&col_ind()[cs()[domain_node/C()] + domain_node%C() + rl()[domain_node] * C()], C());
      }


      /**
       * \brief SparseMatrixELL comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_> friend bool operator== (const SparseMatrixELL & a, const SparseMatrixELL<Mem2_, DT_, IT_> & b)
      {
        CONTEXT("When comparing SparseMatrixELLs");

        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;
        if (a.used_elements() != b.used_elements())
          return false;
        if (a.zero_element() != b.zero_element())
          return false;
        if (a.C() != b.C())
          return false;
        if (a.val_size() != b.val_size())
          return false;

        if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && a.get_indices().size() == 0 && b.get_elements().size() == 0 && b.get_indices().size() == 0)
          return true;

        SparseMatrixELL<Mem::Main, DT_, IT_> ta;
        SparseMatrixELL<Mem::Main, DT_, IT_> tb;
        ta.convert(a);
        tb.convert(b);

        const IT_ * const tacs(ta.cs());
        const IT_ * const tacl(ta.cl());
        const IT_ * const tarl(ta.rl());
        const IT_ * const tacol(ta.col_ind());
        const DT_ * const taval(ta.val());
        const IT_ * const tbcs(tb.cs());
        const IT_ * const tbcl(tb.cl());
        const IT_ * const tbrl(tb.rl());
        const IT_ * const tbcol(tb.col_ind());
        const DT_ * const tbval(tb.val());

        const Index taC(ta.C());
        const Index tarows(ta.rows());

        for (Index i(0); i < ta.num_of_chunks(); ++i)
        {
          if (tacs[i] != tbcs[i] || tacl[i] != tbcl[i])
            return false;
        }

        for (Index i(0); i < ta.rows(); ++i)
        {
          if (tarl[i] != tbrl[i])
            return false;
        }

        for (Index i(0); i < tarows; ++i)
        {
          for (Index j(0); j < tarl[i]; ++j)
          {
            if (tacol[tacs[i/taC] + i%taC + j*taC] != tbcol[tacs[i/taC] + i%taC + j*taC] || taval[tacs[i/taC] + i%taC + j*taC] != tbval[tacs[i/taC] + i%taC + j*taC])
              return false;
          }
        }

        return true;
      }

      /**
       * \brief SparseMatrixELL streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const SparseMatrixELL & b)
      {
        lhs << "[" << std::endl;
        for (Index i(0) ; i < b.rows() ; ++i)
        {
          lhs << "[";
          for (Index j(0) ; j < b.columns() ; ++j)
          {
            lhs << "  " << b(i, j);
          }
          lhs << "]" << std::endl;
        }
        lhs << "]" << std::endl;

        return lhs;
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_MATRIX_ELL_HPP
