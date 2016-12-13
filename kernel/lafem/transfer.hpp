#pragma once
#ifndef KERNEL_LAFEM_TRANSFER_HPP
#define KERNEL_LAFEM_TRANSFER_HPP 1

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
      /// the memory type
      typedef typename Matrix_::MemType MemType;
      /// the data type
      typedef typename Matrix_::DataType DataType;
      /// the index type
      typedef typename Matrix_::IndexType IndexType;

      /// the internal matrix type
      typedef Matrix_ MatrixType;
      /// the compatible vector type
      typedef typename Matrix_::VectorTypeL VectorType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_, typename IT2_>
      using TransferTypeByMDI = class Transfer<typename Matrix_::template ContainerTypeByMDI<Mem2_, DT2_, IT2_>>;

    protected:
      /// the internal prolongation matrix
      Matrix_ _mat_prol;
      /// the internal restriction matrix
      Matrix_ _mat_rest;

    public:
      /// standard constrctor
      Transfer()
      {
      }

      /**
       * \brief Creates the transfer from given prolongation and restiction matrices
       *
       * \param[in] mat_prol
       * The prolongation matrix.
       *
       * \param[in] mat_rest
       * The restriction matrix.
       */
      explicit Transfer(Matrix_&& mat_prol, Matrix_&& mat_rest) :
        _mat_prol(std::forward<Matrix_>(mat_prol)),
        _mat_rest(std::forward<Matrix_>(mat_rest))
      {
      }

      /// move-constructor
      Transfer(Transfer&& other) :
        _mat_prol(std::forward<Matrix_>(other._mat_prol)),
        _mat_rest(std::forward<Matrix_>(other._mat_rest))
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
        return Transfer(_mat_prol.clone(clone_mode), _mat_rest.clone(clone_mode));
      }

      /**
       * \brief Returns the internal data size in bytes.
       *
       * \returns The internal data size in bytes.
       */
      std::size_t bytes() const
      {
        return _mat_prol.bytes() + _mat_rest.bytes();
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
      /// \endcond

      /**
       * \brief Checks whether this transfer is a ghost-operator.
       */
      bool is_ghost() const
      {
        return false;
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
       * \param[in] vec_fine
       * The fine-mesh dual vector to be restricted.
       */
      bool rest_send(const VectorType&) const
      {
        throw InternalError(__func__, __FILE__, __LINE__, "This function must not be called");
      }

      /**
       * \brief Applies the prolongation operator
       *
       * \param[out] vec_fine
       * The prolongated fine-mesh primal vector.
       *
       * \param[in] vec_coarse
       * The coarse-mesh primal vector to be prolongated.
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
       * \param[out] vec_fine
       * The prolongated fine-mesh primal vector.
       */
      bool prol_recv(VectorType& DOXY(vec_fine)) const
      {
        throw InternalError(__func__, __FILE__, __LINE__, "This function must not be called");
      }

      /**
       * \brief Cancels the prolongation.
       */
      void prol_cancel() const
      {
        throw InternalError(__func__, __FILE__, __LINE__, "This function must not be called");
      }
    }; // class Transfer<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_TRANSFER_HPP
