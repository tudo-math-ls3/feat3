#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_BWRAPPEDCSR_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_BWRAPPEDCSR_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Wraps a SparseMatrixCSR to SparseMatrixBCSR
     *
     * \tparam Mem_
     * The memory architecture
     *
     * \tparam DT_
     * Floating point datatype
     *
     * \tparam IT_
     * Index type
     *
     * \tparam BlockSize_
     * Size of the matrice's blocks
     *
     * This class acts as a wrapper around a standard SparseMatrixBCSR  pretending to be a SparseMatrixBCSR.
     * This is required i.e. for the inter-level transfer of vectors in multigrid. This operation is just a
     * multiplication with the transfer matrix, and SparseMatrixCSR can be multiplied with a DenseVectorBlocked
     * in the obvious and desired fashion. The only pitfall is that the corresponding routine only accepts inputs
     * of type TransferMatrix::VectorTypeR.
     *
     * \author Jordi Paul
     *
     */
    template<typename Mem_, typename DT_, typename IT_, int BlockSize_>
    class SparseMatrixBWrappedCSR : public LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>
    {
    public:
      /// The memory architecture
      typedef Mem_ MemType;
      /// The floating point datatype
      typedef DT_ DataType;
      /// The index type
      typedef IT_ IndexType;

      /// The block height the SparseMatrixBWrappedCSR pretends to have
      static constexpr int BlockHeight = BlockSize_;
      /// The block width the SparseMatrixBWrappedCSR pretends to have
      static constexpr int BlockWidth = BlockSize_;

      /// The real type of the underlying matrix
      typedef LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> BaseClass;
      /// What this matrix pretends to be
      typedef LAFEM::SparseMatrixBCSR<Mem_, DT_, IT_, BlockSize_, BlockSize_> PretendType;

      /// Vector type accepted for multiplication form the left
      typedef typename PretendType::VectorTypeL VectorTypeL;
      /// Vector type accepted for multiplication form the right
      typedef typename PretendType::VectorTypeR VectorTypeR;

      /// Inherit base class constructors
      using BaseClass::BaseClass;

      /**
       * \brief Empty standard constructor
       */
      SparseMatrixBWrappedCSR()
      {
      }

      /**
       * \brief From-baseclass move constructor
       *
       * \param[in] other
       * The object to construct from
       *
       */
      explicit SparseMatrixBWrappedCSR(BaseClass&& other) :
        BaseClass(std::forward<BaseClass>(other))
      {
      }

      /**
       * \brief Returns a new compatible vector for left-multiplication
       *
       * \returns An empty left-vector
       *
       */
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(this->rows());
      }

      /**
       * \brief Returns a new compatible vector for right-multiplication
       *
       * \returns An empty right-vector
       *
       */
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->columns());
      }

      /**
       * \brief Clone operation
       *
       * \param[in] mode
       * The LAFEM::CloneMode to use
       *
       */
      SparseMatrixBWrappedCSR clone(CloneMode mode = CloneMode::Weak) const
      {
        return SparseMatrixBWrappedCSR(BaseClass::clone(mode));
      }
    };
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_SPARSE_MATRIX_BWRAPPEDCSR_HPP
