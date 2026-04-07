// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/arch/axpy_dense.hpp>
#include <kernel/lafem/arch/matmatmult_dense_dense.hpp>
#include <kernel/lafem/arch/matmatmult_sparse_dense.hpp>
#include <kernel/lafem/arch/matvecmult_dense_dense.hpp>
#include <kernel/lafem/arch/norm2_dense.hpp>
#include <kernel/lafem/arch/scale_dense.hpp>
#include <kernel/lafem/arch/transpose_dense.hpp>

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
     * _scalar_index[0]: row count \n
     * _scalar_index[1]: column count
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename DT_, typename IT_ = Index>
    class DenseMatrix : public Container<DT_, IT_>
    {
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

    protected:
      Index & _rows()
      {
        return this->_scalar_index.at(0);
      }

      Index & _cols()
      {
        return this->_scalar_index.at(1);
      }

    public:
      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      DenseMatrix() = default;

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       *
       * Creates a matrix with given dimensions.
       */
      explicit DenseMatrix(Index rows_in, Index columns_in)
      {
        XASSERT(rows_in != Index(0) && columns_in != Index(0));

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * rows_in * columns_in));
        this->_elements_size.push_back(rows_in * columns_in);
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
        DenseMatrix(rows_in, columns_in)
      {
        this->format(value);
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit DenseMatrix(std::vector<char> input)
      {
        this->deserialize<DT2_, IT2_>(input);
      }

      //Just a test:
      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       */
      explicit DenseMatrix(FileMode mode, String filename)
      {
        this->read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source filestream.
       */
      explicit DenseMatrix(FileMode mode, std::istream& file)
      {
        this->read_from(mode, file);
      }

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
       * \param[in] file The stream that shall be written to.
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
            //Read in number of num_rows and num_cols
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

            DenseMatrix<DT_, IT_> tmp(Index(trows), tcols);
            Index i(0);
            Memory::TypedView<DT_> tmp_view(tmp.elements_view_w());

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

              tmp_view[i] = tval;

              ++i;
            }
            XASSERTM(i == trows * tcols, "Dense MTX file did not contain enough entries!");

            tmp_view.release();
            this->move(std::move(tmp));

          }
          break;

        case FileMode::fm_dm:
          [[fallthrough]];

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
          buff = new char[LAFEM::file_out_stream_buffer_size];
          file.rdbuf()->pubsetbuf(buff, LAFEM::file_out_stream_buffer_size);
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
            file << this->num_rows() << " " << this->num_cols() << " " << this->num_nzes() << "\n";
            file << std::scientific << std::setprecision(Type::Traits<DT_>::format_precision);

            Memory::TypedView<DT_> elem_view(this->elements_view_r());
            Index tsize = this->num_rows() * this->num_cols();
            for (Index i(0); i < tsize; ++i)
            {
              file << elem_view(i) << "\n";
            }
            break;
          }

        case FileMode::fm_dm:
          [[fallthrough]];

        case FileMode::fm_binary:
          this->template _serialize<double, std::uint64_t>(FileMode::fm_dm, file);
          break;

        default:
          XABORTM("Filemode not supported!");
        }
      }
      //end of test

      /// Retrieve matrix row count
      Index num_rows() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(0);
      }

      /// Retrieve matrix column count
      Index num_cols() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(1);
      }

      /// Retrieve matrix non-zero element count
      Index num_nzes() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(0) * this->_scalar_index.at(1);
      }

      /// Retrieve matrix element count
      Index size() const
      {
        return this->num_nzes();
      }

      /// Checks whether the vector is empty, i.e. if it has size 0
      bool empty() const
      {
        return this->_elements_size.empty() || (this->_elements_size.at(0) <= Index(0));
      }

      /// Returns a reference to the element array arbiter
      Memory::Arbiter& elements_arbiter()
      {
        return this->_elements.front();
      }

      /// Returns a reference to the element array arbiter
      const Memory::Arbiter& elements_arbiter() const
      {
        return this->_elements.front();
      }

      Memory::TypedView<DT_> elements_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<DT_> elements_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<DT_> elements_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<DT_> elements_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, acc));
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
        return VectorTypeL(this->num_rows());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->num_cols());
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
        XASSERTM(x.num_rows() == this->num_rows(), "Row count does not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Column count does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::AxpyDense::template exec<DT_>(this->elements_arbiter(), alpha, x.elements_arbiter(), this->num_nzes());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size());
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DataType result = Arch::Norm2Dense::template exec<DT_>(this->elements_arbiter(), this->num_nzes());

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
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
      void axpy(const DenseMatrix & x, const DT_ alpha = DT_(1))
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Matrix row count does not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Matrix column count does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::AxpyDense::template exec<DT_>(this->elements_arbiter(), alpha, x.elements_arbiter(), this->num_nzes());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size() * 2);
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
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");

        if(this->empty())
          return;

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultDenseDense::template exec<DT_>(r.elements_arbiter(), DT_(1), Memory::Arbiter(),
          this->elements_arbiter(), x.elements_arbiter(), this->num_rows(), this->num_cols(), false);

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
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
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");

        if(this->empty())
          return;

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultDenseDense::template exec<DT_>(r.elements_arbiter(), DT_(1), Memory::Arbiter(),
          this->elements_arbiter(), x.elements_arbiter(), this->num_rows(), this->num_cols(), true);

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
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
      void apply(DenseVector<DT_, IT_> & r, const DenseVector<DT_, IT_> & x, const DenseVector<DT_, IT_> & y, const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows(), "Vector size of y does not match!");

        if(this->empty())
          return;

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultDenseDense::template exec<DT_>(r.elements_arbiter(), alpha, y.elements_arbiter(),
          this->elements_arbiter(), x.elements_arbiter(), this->num_rows(), this->num_cols(), false);

        TimeStamp ts_stop;

        Statistics::add_flops( (this->num_nzes() + this->num_rows()) * 2 );
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
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_cols(), "Vector size of y does not match!");

        if(this->empty())
          return;

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultDenseDense::template exec<DT_>(r.elements_arbiter(), alpha, y.elements_arbiter(),
          this->elements_arbiter(), x.elements_arbiter(), this->num_rows(), this->num_cols(), true);

        TimeStamp ts_stop;

        Statistics::add_flops( (this->num_nzes() + this->num_rows()) * 2 );
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Retrieve the maximum relative difference of this matrix and another one
       * y.max_rel_diff(x) returns  \f$ \max_{0\leq i < n}\frac{|x_i-y_i|}{\max{|x_i|+|y_i|, eps}} \f$
       *
       * \return The largest relative difference.
       */
      DT_ max_rel_diff(const DenseMatrix& x) const
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Matrix row count does not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Matrix column count does not match!");
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DataType max_rel_diff = Arch::MaxRelDiffDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), this->num_nzes());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return max_rel_diff;
      }

      /**
       * \brief Checks if the structural layout of this matrix matches that of another matrix.
       * This excludes comparison of the actual data values.
       *
       * \param[in] x The matrix to compare this matrix to
       *
       * \returns true if the layouts match, false otherwise.
       */
      bool same_layout(const DenseMatrix& x) const
      {
        if (this->num_rows() != x.num_rows())
          return false;
        if (this->num_cols() != x.num_cols())
          return false;

        return true;
      }

      /**
       * \brief Calculate \f$ this \leftarrow x \cdot y \f$
       */
      void set_product(DenseMatrix & x, DenseMatrix & y)
      {
        XASSERTM(x.num_cols() == y.num_rows(), "dimension mismatch!");
        XASSERTM(this->num_rows() == x.num_rows(), "dimension mismatch!");
        XASSERTM(this->num_cols() == y.num_cols(), "dimension mismatch!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::MatMatMultDenseDense::template exec<DT_>(this->elements_arbiter(), DT_(1), x.elements_arbiter(),
          y.elements_arbiter(), Memory::Arbiter(), this->num_rows(), this->num_cols(), x.num_cols());

        TimeStamp ts_stop;
        Statistics::add_flops(x.num_nzes() * y.num_cols()*2);
        Statistics::add_time_blas3(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ this \leftarrow x \cdot y \f$
       *
       * \todo count flops properly
       */
      void set_product(SparseMatrixCSR<DT_, IT_> & x, DenseMatrix & y)
      {
        XASSERTM(x.num_cols() == y.num_rows(), "dimension mismatch!");
        XASSERTM(this->num_rows() == x.num_rows(), "dimension mismatch!");
        XASSERTM(this->num_cols() == y.num_cols(), "dimension mismatch!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::MatMatMultSparseDense::template exec<DT_, IT_>(this->elements_arbiter(), DT_(1), x.val_arbiter(),
          x.col_idx_arbiter(), x.row_ptr_arbiter(), x.num_nzes(), y.elements_arbiter(), Memory::Arbiter(),
          this->num_rows(), this->num_cols(), x.num_cols());

        TimeStamp ts_stop;
        //Statistics::add_flops(x.num_nzes() * y.num_cols()*2);
        Statistics::add_time_blas3(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x y + \beta~ this\f$
       *
       * \param[in] x The first matrix to be scaled with alpha.
       * \param[in] y The second matrix to be multiplied with x.
       * \param[in] alpha A scalar to multiply x with.
       */
      void add_product(
        const DenseMatrix & x,
        const DenseMatrix & y,
        const DT_ alpha = DT_(1))
      {
        XASSERTM(x.num_cols() == y.num_rows(), "dimension mismatch!");
        XASSERTM(this->num_rows() == x.num_rows(), "dimension mismatch!");
        XASSERTM(this->num_cols() == y.num_cols(), "dimension mismatch!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::MatMatMultDenseDense::template exec<DT_>(this->elements_arbiter(), alpha, x.elements_arbiter(),
          y.elements_arbiter(), this->elements_arbiter(), this->num_rows(), this->num_cols(), x.num_cols());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size() * 2);
        Statistics::add_time_blas3(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x y + \beta~ this\f$
       *
       * \param[in] x The first matrix to be scaled with alpha.
       * \param[in] y The second matrix to be multiplied with x.
       * \param[in] alpha A scalar to multiply x with.
       * \param[in] beta A scalar to multiply this with.
       */
      void add_product(
        const SparseMatrixCSR<DT_, IT_> & x,
        const DenseMatrix & y,
        const DT_ alpha = DT_(1))
      {
        XASSERTM(x.num_cols() == y.num_rows(), "dimension mismatch!");
        XASSERTM(this->num_rows() == x.num_rows(), "dimension mismatch!");
        XASSERTM(this->num_cols() == y.num_cols(), "dimension mismatch!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::MatMatMultSparseDense::template exec<DT_, IT_>(this->elements_arbiter(), alpha, x.val_arbiter(),
          x.col_idx_arbiter(), x.row_ptr_arbiter(), x.num_nzes(), y.elements_arbiter(), this->elements_arbiter(),
          this->num_rows(), this->num_cols(), x.num_cols());

        TimeStamp ts_stop;

        Statistics::add_flops(x.num_nzes() * y.num_cols()*2);
        Statistics::add_time_blas3(ts_stop.elapsed(ts_start));
      }

      /// Invert the matrix insitu
      void invert()
      {
        XASSERTM(this->num_rows() == this->num_cols(), "matrix must be square!");
        XASSERT(!this->empty());

        Memory::TypedView<DT_> val = this->elements_view_rw();
        std::vector<IT_> temp(this->num_rows());

        TimeStamp ts_start;

        Math::invert_matrix((IT_)this->num_rows(), (IT_)this->num_rows(), val.get_w(), temp.data());
        /// \todo deal with matrix inversion error

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * this->num_cols()*2);
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
        DenseMatrix r(this->num_cols(), this->num_rows());
        Arch::TransposeDense::template exec<DT_>(r.elements_arbiter(), this->elements_arbiter(), this->num_rows(), this->num_cols());
        return r;
      }

      ///@}

      /**
       * \brief DenseMatrix streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const DenseMatrix & b)
      {
        const Memory::TypedView<DT_> val = b.elements_view_r();

        lhs << "[\n";
        const Index m = b.num_rows();
        const Index n = b.num_cols();
        for (Index i(0) ; i < m ; ++i)
        {
          lhs << "[";
          for (Index j(0) ; j < n ; ++j)
          {
            lhs << "  " << val(i*n + j);
          }
          lhs << "]\n";
        }
        lhs << "]\n";

        return lhs;
      }
    };
  } // namespace LAFEM
} // namespace FEAT
