// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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
#include <kernel/lafem/arch/transpose.hpp>
#include <kernel/lafem/dense_vector.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Dense data matrix class template.
     *
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
    template <typename DT_, typename IT_ = Index>
    class DenseMatrix : public Container<DT_, IT_>
    {
    private:
      Index & _rows()
      {
        return this->_scalar_index.at(1);
      }

      Index & _columns()
      {
        return this->_scalar_index.at(2);
      }

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Compatible L-vector type
      typedef DenseVector<DT_, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<DT_, IT_> VectorTypeR;
      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = DenseMatrix<DT2_, IT2_>;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit DenseMatrix() :
        Container< DT_, IT_> (0)
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
        Container<DT_, IT_>(rows_in * columns_in)
      {
        XASSERT(rows_in != Index(0) && columns_in != Index(0));
        this->_scalar_index.at(0) = rows_in * columns_in;
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);

        this->_elements.push_back(MemoryPool::template allocate_memory<DT_>(this->_scalar_index.at(0)));
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
        Container<DT_, IT_>(rows_in * columns_in)
      {
        XASSERT(rows_in != Index(0) && columns_in != Index(0));
        this->_scalar_index.at(0) = rows_in * columns_in;
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_elements.push_back(MemoryPool::template allocate_memory<DT_>(this->_scalar_index.at(0)));
        this->_elements_size.push_back(this->_scalar_index.at(0));
        MemoryPool::set_memory(this->_elements.at(0), value, this->_scalar_index.at(0));
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
        Container<DT_, IT_>(0)
      {
        deserialize<DT2_, IT2_>(input);
      }

      //Just a test:
      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       */
      explicit DenseMatrix(FileMode mode, String filename) :
      Container<DT_, IT_>(0)
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
      Container<DT_, IT_>(0)
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
        Container<DT_, IT_>(std::forward<DenseMatrix>(other))
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
      template<typename DT2_, typename IT2_>
      void clone(const DenseMatrix<DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Deep)
      {
        Container<DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename DT2_, typename IT2_>
      void convert(const DenseMatrix<DT2_, IT2_> & other)
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
        MemoryPool::synchronize();
        return this->elements()[row * this->columns() + col];
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
        MemoryPool::set_memory(this->_elements.at(0) + row * this->columns() + col, value);
        MemoryPool::synchronize();
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
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_dm, input);
      }

      /**
       * \brief Serialization of complete container entity.
       *
       * \param[in] config LAFEM::SerialConfig, a struct describing the serialize configuration.
       * \note the corresponding configure flags 'zlib' and/or 'zfp' need to be added in the build-id at the configure call.
       *
       * \returns A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialize for details.
       */

      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialize(const LAFEM::SerialConfig& config = LAFEM::SerialConfig()) const
      {
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_dm, config);
      }

      /**
       * \brief Read in matrix from file
       *
       * \param[in] mode The used file format
       * \param[in] filename The file that shall be read in
       */
      //begin of test
      void read_from(FileMode mode, const String& filename)
      {
        std::ios_base::openmode bin = std::ifstream::in | std::ifstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ifstream::in;
        std::ifstream file(filename.c_str(), bin);
        if(! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);
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
        this->clear();

        switch(mode)
        {
          case FileMode::fm_mtx:
          {

            Index trows, tcols;
            String line;
            std::getline(file, line); // !!? Test on overflow error... could be an enormous matrix... !??
            //for now, just array real general (aka dense) matrices
            const bool array_format((line.find("%%MatrixMarket matrix array real general") != String::npos) ? true : false);
            if (array_format == false)
            {
              XABORTM("Input-file is not a compatible array real mtx-file");
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
              trows = Index(atol(srow.c_str()));
              line.erase(0, end);

              begin = line.find_first_not_of(" ");
              line.erase(0, begin);
              end = line.find_first_of(" ");
              String scol(line, 0, end);
              tcols = Index(atol(scol.c_str()));
              line.erase(0, end);
            }

            DenseMatrix<DT_, IT_> result(Index(trows), tcols);
            Index i(0);

            //Read in value of lines:
            while(!file.eof())
            {
              std::getline(file, line);
              if(file.eof())
                break;

              String::size_type begin = line.find_first_not_of(" ");
              line.erase(0, begin);
              String::size_type end = line.find_first_of(" ");
              String sval(line, 0, end);
              DT_ tval((DT_)atof(sval.c_str()));

              Index row(i / tcols);
              Index col(i % tcols);
              result(row, col, tval);

              ++i;
            }
            XASSERTM(i == trows * tcols, "Dense MTX file did not contain enough entries!");

            this->move(std::move(result));
            break;

          }
          case FileMode::fm_dm:
          case FileMode::fm_binary:
            this->template _deserialize<double, std::uint64_t>(FileMode::fm_dm, file);
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
      void write_out(FileMode mode, const String& filename) const
      {
        std::ios_base::openmode bin = std::ofstream::out | std::ofstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ofstream::out;
        std::ofstream file;
        char* buff = nullptr;
        if(mode == FileMode::fm_mtx)
        {
          buff = new char[LAFEM::FileOutStreamBufferSize];
          file.rdbuf()->pubsetbuf(buff, LAFEM::FileOutStreamBufferSize);
        }
        file.open(filename.c_str(), bin);
        if(! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);
        write_out(mode, file);
        file.close();
        delete[] buff;
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
            file << "%%MatrixMarket matrix array real general\n";
            file << this->rows() << " " << this->columns() << " " << this->used_elements() << "\n";

            for(IT_ row(0) ; row < rows(); ++row)
            {
              for(IT_ col(0) ; col < columns() ; ++col)
              {
                file << stringify_fp_sci((*this)(row, col)) << "\n";
              }
            }
            break;
          }
          case FileMode::fm_dm:
          case FileMode::fm_binary:
            this->template _serialize<double, std::uint64_t>(FileMode::fm_dm, file);
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
        Arch::Scale::value(this->elements(), x.elements(), alpha, this->used_elements());

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
        DT_ result = Arch::Norm2::value(this->elements(), this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x + this\f$
       *
       * \param[in] x The first summand vector to be scaled.
       * \param[in] y The second summand vector.
       * \param[in] alpha A scalar to multiply x with.
       */
      void axpy(
        const DenseMatrix & x,
        const DT_ alpha = DT_(1))
      {
        XASSERTM(x.size() == this->size(), "Vector size does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->size() * 2);
        Arch::Axpy::value(this->elements(), alpha, x.elements(), this->size());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<DT_, IT_> & r, const DenseVector<DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements() * 2);
        Arch::Apply::dense(r.elements(), DT_(1), DT_(0), r.elements(), this->elements(),
                                 x.elements(), this->rows(), this->columns());

        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
      * \brief Calculate \f$ r \leftarrow this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      */
      void apply_transposed(DenseVector<DT_, IT_> & r, const DenseVector<DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements() * 2);
        Arch::Apply::dense_transposed(r.elements(), DT_(1), DT_(0), r.elements(), this->elements(),
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
                 DenseVector<DT_, IT_> & r,
                 const DenseVector<DT_, IT_> & x,
                 const DenseVector<DT_, IT_> & y,
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
        Arch::Apply::dense(r.elements(), alpha, DT_(1), y.elements(), this->elements(),
                                 x.elements(), this->rows(), this->columns());

        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      * \param[in] y The summand vector.
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(
        DenseVector<DT_, IT_> & r,
        const DenseVector<DT_, IT_> & x,
        const DenseVector<DT_, IT_> & y,
        const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->columns(), "Vector size of y does not match!");

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        if(Math::abs(alpha) < Math::eps<DT_>())
        {
          r.copy(y);
          return;
        }

        Statistics::add_flops( (this->used_elements() + this->rows()) * 2 );
        Arch::Apply::dense_transposed(r.elements(), alpha, DT_(1), y.elements(), this->elements(),
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

        Arch::ProductMatMat::dense(this->elements(), DT_(1.0), DT_(0.0), x.elements(),
                                         y.elements(), this->elements(), this->rows(), this->columns(), x.columns());

        TimeStamp ts_stop;
        Statistics::add_time_blas3(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ this \leftarrow x \cdot y \f$
       */
      void multiply(SparseMatrixCSR<DT_, IT_> & x, DenseMatrix & y)
      {
        XASSERTM(x.columns() == y.rows(), "dimension mismatch!");
        XASSERTM(this->rows() == x.rows(), "dimension mismatch!");
        XASSERTM(this->columns() == y.columns(), "dimension mismatch!");

        TimeStamp ts_start;
        //Statistics::add_flops(x.used_elements() * y.columns()*2);

        Arch::ProductMatMat::dsd(this->elements(), DT_(1.0), DT_(0.0), x.val(), x.col_ind(), x.row_ptr(), x.used_elements(),
                                         y.elements(), this->rows(), this->columns(), x.columns());

        TimeStamp ts_stop;
        Statistics::add_time_blas3(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x y + \beta~ this\f$
       *
       * \param[in] x The first matrix to be scaled with alpha.
       * \param[in] y The second matrix to be multiplied with x.
       * \param[in] alpha A scalar to multiply x with.
       * \param[in] beta A scalar to multiply z with.
       */
      void multiply(
        const DenseMatrix & x,
        const DenseMatrix & y,
        const DT_ alpha = DT_(1),
        const DT_ beta = DT_(1))
      {
        XASSERTM(x.columns() == y.rows(), "dimension mismatch!");
        XASSERTM(this->rows() == x.rows(), "dimension mismatch!");
        XASSERTM(this->columns() == y.columns(), "dimension mismatch!");

        TimeStamp ts_start;

        Statistics::add_flops(this->size() * 2);
        Arch::ProductMatMat::dense(this->elements(), alpha, beta, x.elements(),
                                         y.elements(), this->elements(), this->rows(), this->columns(), x.columns());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x y + \beta~ this\f$
       *
       * \param[in] x The first matrix to be scaled with alpha.
       * \param[in] y The second matrix to be multiplied with x.
       * \param[in] alpha A scalar to multiply x with.
       * \param[in] beta A scalar to multiply this with.
       */
      void multiply(
        const SparseMatrixCSR<DT_, IT_> & x,
        const DenseMatrix & y,
        const DT_ alpha = DT_(1),
        const DT_ beta = DT_(1))
      {
        XASSERTM(x.columns() == y.rows(), "dimension mismatch!");
        XASSERTM(this->rows() == x.rows(), "dimension mismatch!");
        XASSERTM(this->columns() == y.columns(), "dimension mismatch!");

        TimeStamp ts_start;

        Statistics::add_flops(x.used_elements() * y.columns()*2);
        Arch::ProductMatMat::dsd(this->elements(), alpha, beta, x.val(), x.col_ind(), x.row_ptr(), x.used_elements(),
                                         y.elements(), this->rows(), this->columns(), x.columns());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /// Invert the matrix insitu
      void invert()
      {
        XASSERTM(this->rows() == this->columns(), "matrix must be square!");

        TimeStamp ts_start;
        Statistics::add_flops(this->used_elements() * this->columns()*2);

        IT_ * temp = new IT_[this->rows()];
        Math::invert_matrix((IT_)this->rows(), (IT_)this->rows(), this->elements(), temp);
        delete[] temp;

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

      /**
       * \brief Calculate \f$this^\top \f$
       *
       * \return The transposed matrix
       */
      DenseMatrix transpose() const
      {
        DenseMatrix x_t;
        x_t.transpose(*this);
        return x_t;
      }

      /**
       * \brief Calculate \f$this \leftarrow x^\top \f$
       *
       * \param[in] x The matrix to be transposed.
       *
       * \warning This obviously flips the row- and column count of the matrix
       */
      void transpose(const DenseMatrix & x)
      {
        if (rows() == x.columns() && columns() == x.rows())
        {
          Arch::Transpose::value(this->elements(), x.elements(), x.rows(), x.columns());
        }
        else
        {
          DenseMatrix r(x.columns(), x.rows());
          Arch::Transpose::value(r.elements(), x.elements(), x.rows(), x.columns());
          this->move(std::move(r));
        }
      }

      /**
       * \brief Calculate \f$this^\top \f$ inplace
       *
       * \warning This obviously flips the row- and column count of the matrix
       */
      void transpose_inplace()
      {
        Arch::Transpose::value(this->elements(), this->elements(), this->rows(), this->columns());

        Index t(this->rows());
        this->_rows() = this->columns();
        this->_columns() = t;
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
      friend bool operator== (const DenseMatrix & a, const DenseMatrix<DT_, IT_> & b)
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

        ta = const_cast<DT_*>(a.elements());
        tb = const_cast<DT_*>(b.elements());

        for (Index i(0) ; i < a.size() ; ++i)
        {
          if (ta[i] != tb[i])
          {
            ret = false;
            break;
          }
        }

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
        lhs << "[\n";
        for (Index i(0) ; i < b.rows() ; ++i)
        {
          lhs << "[";
          for (Index j(0) ; j < b.columns() ; ++j)
          {
            lhs << "  " << stringify(b(i, j));
          }
          lhs << "]\n";
        }
        lhs << "]\n";

        return lhs;
      }
    };
  } // namespace LAFEM
} // namespace FEAT
