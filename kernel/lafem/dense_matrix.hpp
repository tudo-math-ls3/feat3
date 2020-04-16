// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_DENSE_MATRIX_HPP
#define KERNEL_LAFEM_DENSE_MATRIX_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/apply.hpp>
#include <kernel/lafem/arch/product_matmat.hpp>
#include <kernel/lafem/dense_vector.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Dense data matrix class template.
     *
     * \tparam Mem_ The \ref FEAT::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a matrix of continuous data in memory. \n\n
     * Data survey: \n
     * _elements[0]: raw number values \n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class DenseMatrix : public Container<Mem_, DT_, IT_>
    {
    private:
      Index & _rows()
      {
        return _rows();
      }

      Index & _columns()
      {
        return _columns();
      }

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Compatible L-vector type
      typedef DenseVector<Mem_, DT_, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<Mem_, DT_, IT_> VectorTypeR;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = DenseMatrix<Mem2_, DT2_, IT2_>;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit DenseMatrix() :
        Container<Mem_, DT_, IT_> (0)
      {
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       *
       * Creates a matrix with given dimensions.
       */
      explicit DenseMatrix(Index rows_in, Index columns_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        XASSERT(rows_in != Index(0) && columns_in != Index(0));
        this->_scalar_index.at(0) = rows_in * columns_in;
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(this->_scalar_index.at(0)));
        this->_elements_size.push_back(this->_scalar_index.at(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] value The value, each element will be set to.
       *
       * Creates a matrix with given dimensions and value.
       */
      explicit DenseMatrix(Index rows_in, Index columns_in, DT_ value) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        XASSERT(rows_in != Index(0) && columns_in != Index(0));
        this->_scalar_index.at(0) = rows_in * columns_in;
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(this->_scalar_index.at(0)));
        this->_elements_size.push_back(this->_scalar_index.at(0));
        MemoryPool<Mem_>::set_memory(this->_elements.at(0), value, this->_scalar_index.at(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit DenseMatrix(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        deserialise<DT2_, IT2_>(input);
      }

      //Just a test:
      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       */
      explicit DenseMatrix(FileMode mode, String filename) :
      Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source filestream.
       */
      explicit DenseMatrix(FileMode mode, std::istream& file) :
      Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, file);
      }
      //end of test
      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      DenseMatrix(DenseMatrix && other) :
        Container<Mem_, DT_, IT_>(std::forward<DenseMatrix>(other))
      {
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      DenseMatrix & operator= (DenseMatrix && other)
      {
        this->move(std::forward<DenseMatrix>(other));

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a clone of this container.
       *
       * \param[in] clone_mode The actual cloning procedure.
       * \returns The created clone.
       *
       */
      DenseMatrix clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        DenseMatrix t;
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
      void clone(const DenseMatrix<Mem2_, DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Deep)
      {
        Container<Mem_, DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const DenseMatrix<Mem2_, DT2_, IT2_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Get a pointer to the data array.
       *
       * \returns Pointer to the data array.
       */
      DT_ * elements()
      {
        return this->_elements.at(0);
      }

      DT_ const * elements() const
      {
        return this->_elements.at(0);
      }

      /**
       * \brief Retrieve specific matrix element.
       *
       * \param[in] row The row of the matrix element.
       * \param[in] col The column of the matrix element.
       *
       * \returns Specific matrix element.
       */
      const DT_ operator()(Index row, Index col) const
      {
        ASSERT(row < this->rows());
        ASSERT(col < this->columns());
        return MemoryPool<Mem_>::get_element(this->_elements.at(0), row * this->columns() + col);
      }

      /**
       * \brief Set specific matrix element.
       *
       * \param[in] row The row of the matrix element.
       * \param[in] col The column of the matrix element.
       * \param[in] value The value to be set.
       */
      void operator()(Index row, Index col, DT_ value)
      {
        ASSERT(row < this->rows());
        ASSERT(col < this->columns());
        MemoryPool<Mem_>::set_memory(this->_elements.at(0) + row * this->columns() + col, value);
      }

      /**
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialise(std::vector<char> input)
      {
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_dm, input);
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * \returns A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialise for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialise()
      {
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_dm);
      }

      /**
       * \brief Read in matrix from file
       *
       * \param[in] mode The used file format
       * \param[in] filename The file that shall be read in
       */
      //begin of test
      void read_from(FileMode mode, String filename)
      {
        std::ios_base::openmode bin = std::ifstream::in | std::ifstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ifstream::in;
        std::ifstream file(filename.c_str(), bin);
        if(! file.is_open())
          XABORTM("Unable to open Vector file " + filename);
        read_from(mode, file);
        file.close();
      }

      /**
      *\brief Read in matrix from file
      *
      * \param[in] mode The used file format
      * \param[in] file The stream that shall be writen to.
      */
      void read_from(FileMode mode, std::istream& file)
      {
        switch(mode)
        {
          case FileMode::fm_mtx:
          {
            std::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >

            IT_ trows;
            Index tcols, ue(0); //ue is not needed...
            String line;
            std::getline(file, line); // !!? Test on overflow error... could be an enormous matrix... !??
            //for now, just symmetric and general matrices...
            const bool general((line.find("%%MatrixMarket matrix coordinate real general") != String::npos) ? true : false);
            const bool symmetric((line.find("%%MatrixMarket matrix coordinate real symmetric") != String::npos) ? true : false);
            if (symmetric == false && general == false)
            {
              XABORTM("Input-file is not a compatible mtx-file");
            }

            while(!file.eof())
            {
              std::getline(file,line);
              if (file.eof())
                XABORTM("Input-file is empty");

              String::size_type begin(line.find_first_not_of(" "));
              if (line.at(begin) != '%')
                break;
            }
            //Read in number of rows and columns
            {
              String::size_type begin(line.find_first_not_of(" "));
              line.erase(0, begin);
              String::size_type end(line.find_first_of(" "));
              String srow(line, 0, end);
              trows = IT_(atol(srow.c_str()));
              line.erase(0, end);

              begin = line.find_first_not_of(" ");
              line.erase(0, begin);
              end = line.find_first_of(" ");
              String scol(line, 0, end);
              tcols = Index(atol(scol.c_str()));
              line.erase(0, end);
            }
            //Read in row, column and value of lines:
            while(!file.eof())
            {
              std::getline(file, line);
              if(file.eof())
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
            //create temp Matrix, write in values, assign temp to this
            {
              DenseMatrix<Mem::Main, DT_, IT_> swapMat(Index(trows), tcols);
              for(auto row : entries)
              {
                for(auto col : row.second)
                {
                  swapMat(Index(row.first),Index(col.first), col.second);
                }
              }
              this->assign(swapMat);
            }
            break;

          }
          case FileMode::fm_dm:
          case FileMode::fm_binary:
            this->template _deserialise<double, uint64_t>(FileMode::fm_dm, file);
            break;
          default:
            XABORTM("Filemode not supported!");
        }
      }

      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        std::ios_base::openmode bin = std::ofstream::out | std::ofstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ofstream::out;
        std::ofstream file(filename.c_str(), bin);
        if(! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);
        write_out(mode, file);
        file.close();
      }

      /**
       * \brief Write outmatrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be written to.
       */
      void write_out(FileMode mode, std::ostream& file) const
      {
        switch(mode)
        {
          case FileMode::fm_mtx:
          {
            DenseMatrix<Mem::Main, DT_, IT_> temp;
            temp.convert(*this);
            file << "%%MatrixMarket matrix coordinate real general" << std::endl;
            file << temp.rows() << " " << temp.columns() << " " << temp.used_elements() << std::endl;

            for(IT_ row(0) ; row < rows(); ++row)
            {
              for(IT_ col(0) ; col < columns() ; ++col)
              {
                if(temp(row, col) != DT_(0))
                  file << row + 1 << " " << col + 1 << " " << std::scientific << temp(row, col) << std::endl;
              }
            }
            break;
          }
          case FileMode::fm_dm:
          case FileMode::fm_binary:
            this->template _serialise<double, uint64_t>(FileMode::fm_dm, file);
            break;
          default:
            XABORTM("Filemode not supported!");
        }
      }
      //end of test

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count.
       */
      Index rows() const
      {
        return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      Index columns() const
      {
        return this->_scalar_index.at(2);
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "DenseMatrix";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const DenseMatrix & x, bool full = false)
      {
        this->_copy_content(x, full);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      template <typename Mem2_>
      void copy(const DenseMatrix<Mem2_, DT_, IT_> & x, bool full = false)
      {
        this->_copy_content(x, full);
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

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index /*row*/) const
      {
        return this->columns();
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow \alpha x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const DenseMatrix & x, const DT_ alpha)
      {
        XASSERTM(x.rows() == this->rows(), "Row count does not match!");
        XASSERTM(x.columns() == this->columns(), "Column count does not match!");
        XASSERTM(x.used_elements() == this->used_elements(), "Nonzero count does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->size());
        Arch::Scale<Mem_>::value(this->elements(), x.elements(), alpha, this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements() * 2);
        DT_ result = Arch::Norm2<Mem_>::value(this->elements(), this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements() * 2);
        Arch::Apply<Mem_>::dense(r.elements(), DT_(1), DT_(0), r.elements(), this->elements(),
                                 x.elements(), this->rows(), this->columns());

        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        if(Math::abs(alpha) < Math::eps<DT_>())
        {
          r.copy(y);
          return;
        }

        Statistics::add_flops( (this->used_elements() + this->rows()) * 2 );
        Arch::Apply<Mem_>::dense(r.elements(), alpha, DT_(1), y.elements(), this->elements(),
                                 x.elements(), this->rows(), this->columns());

        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ this \leftarrow x \cdot y \f$
       */
      void multiply(DenseMatrix & x, DenseMatrix & y)
      {
        XASSERTM(x.columns() == y.rows(), "dimension mismatch!");
        XASSERTM(this->rows() == x.rows(), "dimension mismatch!");
        XASSERTM(this->columns() == y.columns(), "dimension mismatch!");

        TimeStamp ts_start;
        Statistics::add_flops(x.used_elements() * y.columns()*2);

        Arch::ProductMatMat<Mem_>::dense(this->elements(), x.elements(),
                                         y.elements(), this->rows(), this->columns(), x.columns());

        TimeStamp ts_stop;
        Statistics::add_time_blas3(ts_stop.elapsed(ts_start));

      }

      /// Invert the matrix insitu
      void invert()
      {
        XASSERTM(this->rows() == this->columns(), "matrix must be square!");

        TimeStamp ts_start;
        Statistics::add_flops(this->used_elements() * this->columns()*2);

        DenseMatrix<Mem::Main, DT_, IT_> m_main;
        m_main.convert(*this);
        IT_ * temp = new IT_[this->rows()];
        Math::invert_matrix((IT_)this->rows(), (IT_)this->rows(), m_main.elements(), temp);
        delete[] temp;
        this->convert(m_main);

        TimeStamp ts_stop;
        Statistics::add_time_blas3(ts_stop.elapsed(ts_start));
      }

      /// Create an inverse of the current matrix
      DenseMatrix inverse() const
      {
        DenseMatrix result;
        result.clone(*this);
        result.invert();
        return result;
      }
      ///@}

      /// \cond internal

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DT_ * const pval_set, IT_ * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        const DT_ * pval(this->elements());

        for (Index i(0); i < columns(); ++i)
        {
          pval_set[i * stride] = pval[columns() * row + i];
          pcol_set[i * stride] = IT_(i) + IT_(col_start);
        }
      }

      void set_line_reverse(const Index row, DT_ * const pval_set, const Index stride = 1)
      {
        const DT_ * pval(this->elements());

        for (Index i(0); i < columns(); ++i)
        {
          pval_set[i * stride] = pval[columns() * row + i];
        }
      }
      /// \endcond

      /**
       * \brief DenseMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_> friend bool operator== (const DenseMatrix & a, const DenseMatrix<Mem2_, DT_, IT_> & b)
      {
        if (a.size() != b.size())
          return false;
        if (a.get_elements().size() != b.get_elements().size())
          return false;
        if (a.get_indices().size() != b.get_indices().size())
          return false;
        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;

        if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && b.get_elements().size() == 0)
          return true;

        bool ret(true);

        DT_ * ta;
        DT_ * tb;

        if(std::is_same<Mem::Main, Mem_>::value)
        {
          ta = const_cast<DT_*>(a.elements());
        }
        else
        {
          ta = new DT_[a.size()];
          MemoryPool<Mem_>::template download<DT_>(ta, a.elements(), a.size());
        }

        if(std::is_same<Mem::Main, Mem2_>::value)
        {
          tb = const_cast<DT_*>(b.elements());
        }
        else
        {
          tb = new DT_[b.size()];
          MemoryPool<Mem2_>::template download<DT_>(tb, b.elements(), b.size());
        }

        for (Index i(0) ; i < a.size() ; ++i)
        {
          if (ta[i] != tb[i])
          {
            ret = false;
            break;
          }
        }

        if(! std::is_same<Mem::Main, Mem_>::value)
          delete[] ta;
        if(! std::is_same<Mem::Main, Mem2_>::value)
          delete[] tb;

        return ret;
      }

      /**
       * \brief DenseMatrix streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const DenseMatrix & b)
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
} // namespace FEAT

#endif // KERNEL_LAFEM_DENSE_MATRIX_HPP
