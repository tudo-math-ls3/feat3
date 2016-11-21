#pragma once
#ifndef CONTROL_BLOCKED_BASIC_HPP
#define CONTROL_BLOCKED_BASIC_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/comm_base.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>

namespace FEAT
{
  namespace Control
  {
    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename BlockedMatrix_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_> >
    class BlockedBasicSystemLevel
    {
    public:
      static_assert(std::is_same<MemType_, typename BlockedMatrix_::MemType>::value, "MemType mismatch!");
      static_assert(std::is_same<DataType_, typename BlockedMatrix_::DataType>::value, "DataType mismatch!");
      static_assert(std::is_same<IndexType_, typename BlockedMatrix_::IndexType>::value, "IndexType mismatch!");

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;
      static constexpr int dim = dim_;

      /// define local blocked matrix type
      typedef BlockedMatrix_ LocalBlockedMatrix;

      /// define local system matrix type
      typedef BlockedMatrix_ LocalSystemMatrix;

      /// define local system vector type
      typedef typename LocalSystemMatrix::VectorTypeR LocalSystemVector;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define system gate
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;

      /// define global system vector type
      typedef Global::Vector<LocalSystemVector, SystemMirror> GlobalSystemVector;

      /// define global system matrix type
      typedef Global::Matrix<LocalSystemMatrix, SystemMirror, SystemMirror> GlobalSystemMatrix;

      /* ***************************************************************************************** */

      /// our system gate
      SystemGate gate_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;

      /// CTOR
      BlockedBasicSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys)
      {
      }

      virtual ~BlockedBasicSystemLevel()
      {
      }

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const BlockedBasicSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        gate_sys.convert(other.gate_sys);
        matrix_sys.convert(&gate_sys, &gate_sys, other.matrix_sys);
      }
    }; // class BlockedBasicSystemLevel<...>


    template<
      int dim_,
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename BlockedMatrix_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_> >
    class BlockedUnitFilterSystemLevel :
      public BlockedBasicSystemLevel<dim_, MemType_, DataType_, IndexType_, BlockedMatrix_>
    {
    public:
      typedef BlockedBasicSystemLevel<dim_, MemType_, DataType_, IndexType_, BlockedMatrix_> BaseClass;

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define local filter type
      typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename Mem2_, typename DT2_, typename IT2_, typename BlockedMatrix2_>
      using BaseType = class BlockedUnitFilterSystemLevel<dim_, Mem2_, DT2_, IT2_, BlockedMatrix2_>;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      BlockedUnitFilterSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->matrix_sys.bytes () + filter_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source BlockedUnitFilterSystemLevel content as content of current BlockedUnitFilterSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const BlockedUnitFilterSystemLevel<dim_, M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }
    }; // class BlockedUnitFilterSystemLevel<...>

    template<typename SystemLevel_, //typename BlockedMatrix_ = typename SystemLevel_::LocalBlockedMatrix>
      typename BlockedMatrix_ = LAFEM::SparseMatrixBWrappedCSR<
        typename SystemLevel_::MemType,
        typename SystemLevel_::DataType,
        typename SystemLevel_::IndexType,
        SystemLevel_::dim>
      >
    class BlockedBasicTransferLevel
    {
    public:
      /// define system mirror type
      typedef typename SystemLevel_::SystemMirror SystemMirror;

      /// our local transfer matrix type
      typedef BlockedMatrix_ LocalSystemTransferMatrix;

      /// our global transfer matrix type
      typedef Global::Matrix<LocalSystemTransferMatrix, SystemMirror, SystemMirror> GlobalSystemTransferMatrix;

      /// Our class base type
      template <typename SystemLevel2_>
      using BaseType = class BlockedBasicTransferLevel<SystemLevel2_>;

      /// our global transfer matrices
      GlobalSystemTransferMatrix prol_sys;

      /// \copydoc ScalarBasicTransferLevel::prol_sys
      GlobalSystemTransferMatrix rest_sys;

      BlockedBasicTransferLevel()
      {
      }

      /// CTOR
      explicit BlockedBasicTransferLevel(SystemLevel_& lvl_coarse, SystemLevel_& lvl_fine) :
        prol_sys(&lvl_fine.gate_sys, &lvl_coarse.gate_sys),
        rest_sys(&lvl_coarse.gate_sys, &lvl_fine.gate_sys)
      {
      }

      virtual ~BlockedBasicTransferLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return prol_sys.bytes() + rest_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source BlockedBasicTransferLevel content as content of current BlockedBasicTransferLevel.
       *
       * \warning The provided SystemLevels must already be converted to the matching
       * configuration, as they contain the used gateways.
       *
       */
      template <typename SL_, typename SM_>
      void convert(SystemLevel_ & lvl_coarse , SystemLevel_ & lvl_fine, const BlockedBasicTransferLevel<SL_, SM_> & other)
      {
        prol_sys.convert(&lvl_fine.gate_sys, &lvl_coarse.gate_sys, other.prol_sys);
        rest_sys.convert(&lvl_coarse.gate_sys, &lvl_fine.gate_sys, other.rest_sys);
      }
    }; // class BlockedBasicTransferLevel<...>

  } // namespace Control
} // namespace FEAT

#endif // CONTROL_BLOCKED_BASIC_HPP
