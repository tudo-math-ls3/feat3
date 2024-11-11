// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/lafem/base.hpp>
#include <kernel/util/exception.hpp>

#include <utility>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Grid-Transfer operator class template
     *
     * \tparam Matrix_
     * The matrix type that is to be used for the prolongation and restriction matrices.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_>
    class Transfer
    {
    public:
      /// the data type
      typedef typename Matrix_::DataType DataType;
      /// the index type
      typedef typename Matrix_::IndexType IndexType;

      /// the internal matrix type
      typedef Matrix_ MatrixType;
      /// the compatible vector type
      typedef typename Matrix_::VectorTypeL VectorType;

      /// Our 'base' class type
      template <typename DT2_, typename IT2_>
      using TransferTypeByDI = Transfer<typename Matrix_::template ContainerTypeByDI<DT2_, IT2_>>;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

    protected:
      /// the internal prolongation matrix
      Matrix_ _mat_prol;
      /// the internal restriction matrix
      Matrix_ _mat_rest;
      /// the internal truncation matrix
      Matrix_ _mat_trunc;

    public:
      /// standard constructor
      Transfer()
      {
      }

      /**
       * \brief Creates the transfer from given prolongation and restriction matrices
       *
       * \param[in] mat_prol
       * The prolongation matrix.
       *
       * \param[in] mat_rest
       * The restriction matrix.
       */
      explicit Transfer(Matrix_&& mat_prol, Matrix_&& mat_rest) :
        _mat_prol(std::forward<Matrix_>(mat_prol)),
        _mat_rest(std::forward<Matrix_>(mat_rest)),
        _mat_trunc()
      {
      }


      /**
       * \brief Creates the transfer from given prolongation, restriction and reduction matrices
       *
       * \param[in] mat_prol
       * The prolongation matrix.
       *
       * \param[in] mat_rest
       * The restriction matrix.
       *
       * \param[in] mat_trunc
       * The truncation matrix.
       */
      explicit Transfer(Matrix_&& mat_prol, Matrix_&& mat_rest, Matrix_&& mat_trunc) :
        _mat_prol(std::forward<Matrix_>(mat_prol)),
        _mat_rest(std::forward<Matrix_>(mat_rest)),
        _mat_trunc(std::forward<Matrix_>(mat_trunc))
      {
      }

      /// move-constructor
      Transfer(Transfer&& other) :
        _mat_prol(std::forward<Matrix_>(other._mat_prol)),
        _mat_rest(std::forward<Matrix_>(other._mat_rest)),
        _mat_trunc(std::forward<Matrix_>(other._mat_trunc))
      {
      }

      /// virtual destructor
      virtual ~Transfer()
      {
      }

      /// move-assignment operator
      Transfer& operator=(Transfer&& other)
      {
        if(this == &other)
          return *this;

        _mat_prol = std::forward<Matrix_>(other._mat_prol);
        _mat_rest = std::forward<Matrix_>(other._mat_rest);
        _mat_trunc = std::forward<Matrix_>(other._mat_trunc);

        return *this;
      }

      /// container conversion function
      template<typename Matrix2_>
      void convert(const Transfer<Matrix2_>& other)
      {
        if((void*)this == (void*)&other)
          return;

        _mat_prol.convert(other.get_mat_prol());
        _mat_rest.convert(other.get_mat_rest());
        _mat_trunc.convert(other.get_mat_trunc());
      }

      /**
       * \brief Creates a clone of this object.
       *
       * \param[in] clone_mode
       * The desired clone mode
       *
       * \returns A clone of this object
       */
      Transfer clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        return Transfer(_mat_prol.clone(clone_mode), _mat_rest.clone(clone_mode), _mat_trunc.clone(clone_mode));
      }

      /**
       * \brief Returns the internal data size in bytes.
       *
       * \returns The internal data size in bytes.
       */
      std::size_t bytes() const
      {
        return _mat_prol.bytes() + _mat_rest.bytes() + _mat_trunc.bytes();
      }

      void compile()
      {
        // nothing to do here...
      }

      /// \cond internal
      Matrix_& get_mat_prol()
      {
        return _mat_prol;
      }

      const Matrix_& get_mat_prol() const
      {
        return _mat_prol;
      }

      Matrix_& get_mat_rest()
      {
        return _mat_rest;
      }

      const Matrix_& get_mat_rest() const
      {
        return _mat_rest;
      }

      Matrix_& get_mat_trunc()
      {
        return _mat_trunc;
      }

      const Matrix_& get_mat_trunc() const
      {
        return _mat_trunc;
      }
      /// \endcond

      /**
       * \brief Checks whether this transfer is a ghost-operator.
       */
      bool is_ghost() const
      {
        return false;
      }

      /**
       * \brief Applies the truncation operator
       *
       * \param[in] vec_fine
       * The fine-mesh primal vector to be truncated.
       *
       * \param[out] vec_coarse
       * The truncated coarse-mesh primal vector.
       *
       * \returns \c true
       */
      bool trunc(const VectorType& vec_fine, VectorType& vec_coarse) const
      {
        _mat_trunc.apply(vec_coarse, vec_fine);
        return true;
      }

      /**
       * \brief Sends the truncation for a ghost operator
       *
       * \attention
       * This function is defined only for compatibility with the Global::Transfer class
       * and must not be called, as this function will throw an exception otherwise.
       *
       * \param[in] vec_fine
       * The fine-mesh primal vector to be truncated.
       */
      bool trunc_send(const VectorType& DOXY(vec_fine)) const
      {
        XABORTM("This function must not be called");
      }

      /**
       * \brief Applies the restriction operator
       *
       * \param[in] vec_fine
       * The fine-mesh dual vector to be restricted.
       *
       * \param[out] vec_coarse
       * The restricted coarse-mesh dual vector.
       *
       * \returns \c true
       */
      bool rest(const VectorType& vec_fine, VectorType& vec_coarse) const
      {
        _mat_rest.apply(vec_coarse, vec_fine);
        return true;
      }

      /**
       * \brief Sends the restriction for a ghost operator
       *
       * \attention
       * This function is defined only for compatibility with the Global::Transfer class
       * and must not be called, as this function will throw an exception otherwise.
       *
       * \param[in] vec_fine
       * The fine-mesh dual vector to be restricted.
       */
      bool rest_send(const VectorType& DOXY(vec_fine)) const
      {
        XABORTM("This function must not be called");
      }

      /**
       * \brief Applies the prolongation operator
       *
       * \param[out] vec_fine
       * The prolonged fine-mesh primal vector.
       *
       * \param[in] vec_coarse
       * The coarse-mesh primal vector to be prolonged.
       *
       * \returns \c true
       */
      bool prol(VectorType& vec_fine, const VectorType& vec_coarse) const
      {
        _mat_prol.apply(vec_fine, vec_coarse);
        return true;
      }

      /**
       * \brief Receives the prolongation for a ghost operator
       *
       * \attention
       * This function is defined only for compatibility with the Global::Transfer class
       * and must not be called, as this function will throw an exception otherwise.
       *
       * \param[out] vec_fine
       * The prolonged fine-mesh primal vector.
       */
      bool prol_recv(VectorType& DOXY(vec_fine)) const
      {
        XABORTM("This function must not be called");
      }

      /**
       * \brief Cancels the prolongation.
       *
       * \attention
       * This function is defined only for compatibility with the Global::Transfer class
       * and must not be called, as this function will throw an exception otherwise.
       */
      void prol_cancel() const
      {
        XABORTM("This function must not be called");
      }
    }; // class Transfer<...>
  } // namespace LAFEM
} // namespace FEAT
