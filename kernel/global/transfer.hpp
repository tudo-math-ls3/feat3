// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GLOBAL_TRANSFER_HPP
#define KERNEL_GLOBAL_TRANSFER_HPP 1

#include <kernel/global/muxer.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/lafem/transfer.hpp>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Global grid-transfer operator class template
     *
     * \author Peter Zajac
     */
    template<typename LocalTransfer_, typename Mirror_>
    class Transfer
    {
    public:
      /// our local transfer
      typedef LocalTransfer_ LocalTransferType;
      /// our internal local matrix type
      typedef typename LocalTransfer_::MatrixType LocalMatrixType;
      /// our internal local vector type
      typedef typename LocalTransfer_::VectorType LocalVectorType;

      /// our global vector type
      typedef Global::Vector<LocalVectorType, Mirror_> VectorType;

      /// our coarse grid multiplexer type
      typedef Global::Muxer<LocalVectorType, Mirror_> MuxerType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_, typename IT2_>
      using TransferTypeByMDI = Transfer<
        typename LocalTransfer_::template TransferTypeByMDI<Mem2_, DT2_, IT2_>,
        typename Mirror_::template MirrorType<Mem2_, DT2_, IT2_> >;

      static constexpr bool is_global = true;
      static constexpr bool is_local = false;

    public:
      /// the coarse-level multiplexer
      MuxerType* _coarse_muxer;
      /// the local transfer operator
      LocalTransfer_ _transfer;
      /// a temporary local vector
      mutable LocalVectorType _vec_tmp;

    public:
      /// standard constructor
      Transfer() :
        _coarse_muxer(nullptr)
      {
      }

      /// move-constructor
      Transfer(Transfer&& other) :
        _coarse_muxer(other._coarse_muxer),
        _transfer(std::forward<LocalTransfer_>(other._transfer)),
        _vec_tmp(std::forward<LocalVectorType>(other._vec_tmp))
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] coarse_muxer
       * The coarse-level muxer.
       *
       * \param[in] args
       * The arguments that are passed to the local transfer operator constructor.
       */
      template<typename... Args_>
      explicit Transfer(MuxerType* coarse_muxer, Args_&&... args) :
        _coarse_muxer(coarse_muxer),
        _transfer(std::forward<Args_>(args)...),
        _vec_tmp(_transfer.get_mat_rest().create_vector_l())
      {
      }

      /// virtual destructor
      virtual ~Transfer()
      {
      }

      /// move-assign operator
      Transfer& operator=(Transfer&& other)
      {
        if(this == &other)
          return *this;

        _coarse_muxer = other._coarse_muxer;
        _transfer = std::forward<LocalTransfer_>(other._transfer);
        _vec_tmp = std::forward<LocalVectorType>(other._vec_tmp);

        return *this;
      }

      /// container conversion function
      template<typename LocalTransfer2_, typename Mirror2_>
      void convert(MuxerType* coarse_muxer, const Transfer<LocalTransfer2_, Mirror2_>& other)
      {
        if((void*)this == (void*)&other)
          return;

        _coarse_muxer = coarse_muxer;
        _transfer.convert(other._transfer);
        _vec_tmp.convert(other._vec_tmp);
      }

      /**
       * \brief Creates a clone of this object.
       *
       * \param[in] clone_mode
       * The desired clone mode
       *
       * \returns A clone of this object
       */
      Transfer clone(LAFEM::CloneMode clone_mode = LAFEM::CloneMode::Weak) const
      {
        return Transfer(_coarse_muxer, _transfer.clone(clone_mode));
      }

      /// \returns the local transfer operator
      LocalTransfer_& local()
      {
        return _transfer;
      }

      /// \returns the local transfer operator
      const LocalTransfer_& local() const
      {
        return _transfer;
      }

      /// \returns the internal data size in bytes
      std::size_t bytes() const
      {
        return _transfer.bytes() + _vec_tmp.bytes();
      }

      void compile()
      {
        _transfer.compile();
        _vec_tmp = _transfer.get_mat_rest().create_vector_l();
      }

      /// \cond internal
      LocalVectorType& get_vec_temp()
      {
        return _vec_tmp;
      }

      LocalMatrixType& get_mat_prol()
      {
        return _transfer.get_mat_prol();
      }

      const LocalMatrixType& get_mat_prol() const
      {
        return _transfer.get_mat_prol();
      }

      LocalMatrixType& get_mat_rest()
      {
        return _transfer.get_mat_rest();
      }

      const LocalMatrixType& get_mat_rest() const
      {
        return _transfer.get_mat_rest();
      }

      LocalMatrixType& get_mat_trunc()
      {
        return _transfer.get_mat_trunc();
      }

      const LocalMatrixType& get_mat_trunc() const
      {
        return _transfer.get_mat_trunc();
      }
      /// \endcond

      /**
       * \brief Checks whether this transfer is a ghost-operator.
       */
      bool is_ghost() const
      {
        return (_coarse_muxer != nullptr) && _coarse_muxer->is_ghost();
      }

      /**
       * \brief Applies the truncation operator
       *
       * \param[in] vec_fine
       * The fine-mesh dual vector to be truncated.
       *
       * \param[out] vec_coarse
       * The truncated coarse-mesh dual vector.
       *
       * \returns \c true
       */
      bool trunc(const VectorType& vec_fine, VectorType& vec_coarse) const
      {
        if((_coarse_muxer == nullptr) || (!_coarse_muxer->is_child()))
        {
          _transfer.trunc(vec_fine.local(), vec_coarse.local());
        }
        else
        {
          XASSERT(_coarse_muxer->is_parent());
          _transfer.trunc(vec_fine.local(), _vec_tmp);
          _coarse_muxer->join(_vec_tmp, vec_coarse.local());
        }
        vec_coarse.sync_0();
        return true;
      }

      /**
       * \brief Sends the truncation for a ghost operator
       *
       * \param[in] vec_fine
       * The fine-mesh dual vector to be truncated.
       */
      bool trunc_send(const VectorType& vec_fine) const
      {
        XASSERT(_coarse_muxer != nullptr);
        XASSERT(_coarse_muxer->is_ghost());
        _transfer.trunc(vec_fine.local(), _vec_tmp);
        _coarse_muxer->join_send(_vec_tmp);
        return true;
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
        if((_coarse_muxer == nullptr) || (!_coarse_muxer->is_child()))
        {
          _transfer.rest(vec_fine.local(), vec_coarse.local());
        }
        else
        {
          XASSERT(_coarse_muxer->is_parent());
          _transfer.rest(vec_fine.local(), _vec_tmp);
          _coarse_muxer->join(_vec_tmp, vec_coarse.local());
        }
        vec_coarse.sync_0();
        return true;
      }

      /**
       * \brief Sends the restriction for a ghost operator
       *
       * \param[in] vec_fine
       * The fine-mesh dual vector to be restricted.
       */
      bool rest_send(const VectorType& vec_fine) const
      {
        XASSERT(_coarse_muxer != nullptr);
        XASSERT(_coarse_muxer->is_ghost());
        _transfer.rest(vec_fine.local(), _vec_tmp);
        _coarse_muxer->join_send(_vec_tmp);
        return true;
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
        if((_coarse_muxer == nullptr) || (!_coarse_muxer->is_child()))
        {
          _transfer.prol(vec_fine.local(), vec_coarse.local());
        }
        else
        {
          XASSERT(_coarse_muxer->is_child());
          XASSERT(_coarse_muxer->is_parent());
          _coarse_muxer->split(_vec_tmp, vec_coarse.local());
          _transfer.prol(vec_fine.local(), _vec_tmp);
        }
        vec_fine.sync_0();
        return true;
      }

      /**
       * \brief Cancels the prolongation.
       */
      void prol_cancel() const
      {
        XASSERTM(false, "this function must not be called");
      }

      /**
       * \brief Receives the prolongation for a ghost operator
       *
       * \param[out] vec_fine
       * The prolonged fine-mesh primal vector.
       */
      bool prol_recv(VectorType& vec_fine) const
      {
        XASSERT(_coarse_muxer != nullptr);
        XASSERT(_coarse_muxer->is_ghost());
        _coarse_muxer->split_recv(_vec_tmp);
        _transfer.prol(vec_fine.local(), _vec_tmp);
        vec_fine.sync_0();
        return true;
      }
    }; // class Transfer<..>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_TRANSFER_HPP
