#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/matrix_base.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/adjacency/graph.hpp>

#include <fstream>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief sparse banded matrix
     *
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its diagonal entries
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: vector of offsets (bottom-left diagonal has offset 0,
     *                                 main diagonal has offset rows - 1 and
     *                                 top-right diagonal has offset row + columns - 1)\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_index[4]: number of offsets \n
     * _scalar_dt[0]: zero element
     *
     * This class saves a sparse-matrix with a banded structure. For each diagonal of
     * the matrix with non-zero elements there must be reserved memory for the whole
     * diagonal. For faster access on the matrix-elements each diagonal get the virtual
     * length of the row-count of the matrix. They are enlarged to the left and right
     * side of the matrix as shown in the following layout.
     *       +--                  --+
     *  \    | \           \      \ |
     *   \   |\ \           \      \|
     *    \  | \ \           \      |
     *     \ |  \ \           \     |\
     *      \|   \ \           \    | \
     *       |    \ \           \   |  \
     *       |\    \ \           \  |   \
     *       +--                  --+
     * To get the position of the diagonals in the matrix, the matching offsets are
     * saved from left to right in the offsets-array.
     * - The first diagonal is the one at the bottom-left and gets the offset = 1,
     * - the main diaognal has the offset = rows - 1
     * - and the last offset at the top-right has the offset = rows + columns - 1.
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseMatrixBanded : public Container<Mem_, DT_, IT_>, public MatrixBase
    {
    private:
      void _read_from_bm(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        _read_from_bm(file);
        file.close();
      }

      void _read_from_bm(std::istream& file)
      {
        this->template _deserialise<double, uint64_t>(FileMode::fm_bm, file);
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

      Index & _used_elements()
      {
        return this->_scalar_index.at(3);
      }

      Index & _num_of_offsets()
      {
        return this->_scalar_index.at(4);
      }

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Compatible L-vector type
      typedef DenseVector<MemType, DataType, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<MemType, DataType, IT_> VectorTypeR;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit SparseMatrixBanded() :
        Container<Mem_, DT_, IT_> (0)
      {
        CONTEXT("When creating SparseMatrixBanded");
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] val_in The vector with non zero elements.
       * \param[in] offsets_in The vector of offsets.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixBanded(const Index rows_in, const Index columns_in,
                                  DenseVector<Mem_, DT_, IT_> & val_in,
                                  DenseVector<Mem_, IT_, IT_> & offsets_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        if (val_in.size() != rows_in * offsets_in.size())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Size of values does not match to number of offsets and row count!");
        }

        CONTEXT("When creating SparseMatrixBanded");
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);

        Index tused_elements(0);

        for (Index i(0); i < offsets_in.size(); ++i)
        {
          const Index toffset(offsets_in(i));

          if (toffset + 2 > rows_in + columns_in)
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Offset out of matrix!");
          }

          tused_elements += columns_in + Math::min(rows_in, columns_in + rows_in - toffset - 1) - Math::max(columns_in + rows_in - toffset - 1, columns_in);
        }

        this->_scalar_index.push_back(tused_elements);
        this->_scalar_index.push_back(offsets_in.size());
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(offsets_in.elements());
        this->_indices_size.push_back(offsets_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a banded matrix based on the source filestream.
       */
      explicit SparseMatrixBanded(FileMode mode, std::istream& file) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixBanded");

        switch(mode)
        {
          case FileMode::fm_bm:
            _read_from_bm(file);
            break;
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
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
      explicit SparseMatrixBanded(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixBanded");
        deserialise<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixBanded(SparseMatrixBanded && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixBanded>(other))
      {
        CONTEXT("When moving SparseMatrixBanded");
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixBanded & operator= (SparseMatrixBanded && other)
      {
        CONTEXT("When moving SparseMatrixBanded");

        this->move(std::forward<SparseMatrixBanded>(other));

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a deep copy of itself.
       *
       * \param[in] clone_indices Should we create a deep copy of the index arrays, too ?
       *
       * \return A deep copy of itself.
       *
       */
      SparseMatrixBanded clone(bool clone_indices = false) const
      {
        SparseMatrixBanded t;
        t.clone(*this, clone_indices);
        return t;
      }

      using Container<Mem_, DT_, IT_>::clone;

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixBanded<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting SparseMatrixBanded");
        this->assign(other);
      }

      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        CONTEXT("When writing out SparseMatrixBanded");

        switch(mode)
        {
          case FileMode::fm_bm:
            write_out_bm(filename);
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
        CONTEXT("When writing out SparseMatrixBanded");

        switch(mode)
        {
          case FileMode::fm_bm:
            write_out_bm(file);
            break;
          default:
            throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Write out matrix to banded binary file.
       *
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out_bm(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        write_out_bm(file);
        file.close();
      }

      /**
       * \brief Write out matrix to banded binary file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_bm(std::ostream& file) const
      {
        if (! std::is_same<DT_, double>::value)
          std::cout<<"Warning: You are writing out a banded matrix that is not double precision!"<<std::endl;

        this->template _serialise<double, uint64_t>(FileMode::fm_bm, file);
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
        CONTEXT("When retrieving SparseMatrixBanded element");

        ASSERT(row < rows(), "Error: " + stringify(row) + " exceeds sparse matrix banded row size " + stringify(rows()) + " !");
        ASSERT(col < columns(), "Error: " + stringify(col) + " exceeds sparse matrix banded column size " + stringify(columns()) + " !");

        const Index trows(this->_scalar_index.at(1));

        for (Index i(0); i < this->_scalar_index.at(4); ++i)
        {
          const Index toffset(Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i));
          if (row + toffset + 1 == col + trows)
          {
            return Util::MemoryPool<Mem_>::get_element(this->_elements.at(0), i * trows + row);
          }
        }
        return zero_element();
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
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      const Index & used_elements() const override
      {
        return this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve number of offsets.
       *
       * \returns Number of offsets.
       */
      const Index & num_of_offsets() const
      {
        return this->_scalar_index.at(4);
      }

      /**
       * \brief Retrieve element array.
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
       * \brief Retrieve offsets array.
       *
       * \returns offsets array.
       */
      IT_ * offsets()
      {
        return this->_indices.at(0);
      }

      IT_ const * offsets() const
      {
        return this->_indices.at(0);
      }

      /**
       * \brief Retrieve non zero element.
       *
       * \returns Non zero element.
       */
      const DT_ zero_element() const
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
        return "SparseMatrixBanded";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      void copy(const SparseMatrixBanded & x)
      {
        this->_copy_content(x);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      template <typename Mem2_>
      void copy(const SparseMatrixBanded<Mem2_, DT_, IT_> & x)
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
                const SparseMatrixBanded & x,
                const SparseMatrixBanded & y,
                const DT_ alpha = DT_(1))
      {
        if (x.rows() != y.rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.num_of_offsets() != y.num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.num_of_offsets() != this->num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.used_elements() != y.used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

        // check for special cases
        // r <- x + y
        if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
          Arch::Sum<Mem_, Algo_>::value(this->val(), x.val(), y.val(), this->rows() * this->num_of_offsets());
        // r <- y - x
        else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          Arch::Difference<Mem_, Algo_>::value(this->val(), y.val(), x.val(), this->rows() * this->num_of_offsets());
        // r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          this->copy(y);
        // r <- y + alpha*x
        else
          Arch::Axpy<Mem_, Algo_>::dv(this->val(), alpha, x.val(), y.val(), this->rows() * this->num_of_offsets());
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha x \f$
       *
       * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      template <typename Algo_>
      void scale(const SparseMatrixBanded & x, const DT_ alpha)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.num_of_offsets() != this->num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

        Arch::Scale<Mem_, Algo_>::value(this->val(), x.val(), alpha, this->rows() * this->num_of_offsets());
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      template <typename Algo_>
      DT_ norm_frobenius() const
      {
        return Arch::Norm2<Mem_, Algo_>::value(this->val(), this->rows() * this->num_of_offsets());
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
      void apply(DenseVector<Mem_,DT_, IT_>& r, const DenseVector<Mem_, DT_, IT_>& x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        Arch::ProductMatVec<Mem_, Algo_>::banded(r.elements(),
                                                 this->val(),
                                                 this->offsets(),
                                                 x.elements(),
                                                 this->num_of_offsets(),
                                                 this->rows(),
                                                 this->columns());
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      template<typename Algo_>
      void apply(DenseVector<Mem_,DT_, IT_>& r,
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
          Arch::Defect<Mem_, Algo_>::banded(r.elements(),
                                            y.elements(),
                                            this->val(),
                                            this->offsets(),
                                            x.elements(),
                                            this->num_of_offsets(),
                                            this->rows(),
                                            this->columns());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_, Algo_>::banded(r.elements(),
                                          y.elements(),
                                          alpha,
                                          this->val(),
                                          this->offsets(),
                                          x.elements(),
                                          this->num_of_offsets(),
                                          this->rows(),
                                          this->columns());
        }
      }
      ///@}
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
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_bm, input);
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
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_bm);
      }

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

      /// Returns first row-index of the diagonal matching to the offset i
      Index start_offset(const Index i) const
      {
        return Arch::Intern::ProductMatVecBanded::start_offset(i, offsets(), rows(), columns(), num_of_offsets());
      }

      /// Returns last row-index of the diagonal matching to the offset i
      Index end_offset(const Index i) const
      {
        return Arch::Intern::ProductMatVecBanded::end_offset(i, offsets(), rows(), columns(), num_of_offsets());
      }

      /**
       * \brief SparseMatrixBanded comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_> friend bool operator== (const SparseMatrixBanded & a, const SparseMatrixBanded<Mem2_, DT_, IT_> & b)
      {
        CONTEXT("When comparing SparseMatrixBandeds");

        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;
        if (a.num_of_offsets() != b.num_of_offsets())
          return false;
        if (a.used_elements() != b.used_elements())
          return false;
        if (a.zero_element() != b.zero_element())
          return false;

        if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && a.get_indices().size() == 0 && b.get_elements().size() == 0 && b.get_indices().size() == 0)
          return true;

        IT_ * offsets_a;
        IT_ * offsets_b;
        DT_ * val_a;
        DT_ * val_b;

        if(std::is_same<Mem::Main, Mem_>::value)
        {
          offsets_a = (IT_*)a.offsets();
          val_a = (DT_*)a.val();
        }
        else
        {
          offsets_a = new IT_[a.num_of_offsets()];
          Util::MemoryPool<Mem_>::instance()->template download<IT_>(offsets_a, a.offsets(), a.num_of_offsets());
          val_a = new DT_[a.num_of_offsets() * a.rows()];
          Util::MemoryPool<Mem_>::instance()->template download<DT_>(val_a, a.val(), a.num_of_offsets() * a.rows());
        }
        if(std::is_same<Mem::Main, Mem2_>::value)
        {
          offsets_b = (IT_*)b.offsets();
          val_b = (DT_*)b.val();
        }
        else
        {
          offsets_b = new IT_[b.num_of_offsets()];
          Util::MemoryPool<Mem2_>::instance()->template download<IT_>(offsets_b, b.offsets(), b.num_of_offsets());
          val_b = new DT_[b.num_of_offsets() * b.rows()];
          Util::MemoryPool<Mem2_>::instance()->template download<DT_>(val_b, b.val(), b.num_of_offsets() * b.rows());
        }

        bool ret(true);

        for (Index i(0); i < a.num_of_offsets(); ++i)
        {
          if (offsets_a[i] != offsets_b[i])
          {
            ret = false;
            break;
          }
        }

        for (Index i(0) ; i < a.num_of_offsets() * a.rows() ; ++i)
        {
          if (val_a[i] != val_b[i])
          {
            ret = false;
            break;
          }
        }

        if(! std::is_same<Mem::Main, Mem_>::value)
        {
          delete[] offsets_a;
          delete[] val_a;
        }
        if(! std::is_same<Mem::Main, Mem2_>::value)
        {
          delete[] offsets_b;
          delete[] val_b;
        }

        return ret;
      }

      /**
       * \brief SparseMatrixBanded streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const SparseMatrixBanded & b)
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

#endif // KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP
