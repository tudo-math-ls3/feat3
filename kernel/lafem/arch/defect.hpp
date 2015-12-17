#pragma once
#ifndef KERNEL_LAFEM_ARCH_DEFECT_HPP
#define KERNEL_LAFEM_ARCH_DEFECT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <typeinfo>


namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct Defect;

      template <>
      struct Defect<Mem::Main>
      {
        template <typename DT_, typename IT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
          csr_generic(r, rhs, val, col_ind, row_ptr, x, rows, columns, used_elements);
        }

        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
#ifdef FEAST_BACKENDS_MKL
          csr_mkl(r, rhs, val, col_ind, row_ptr, x, rows, columns, used_elements);
#else
          csr_generic(r, rhs, val, col_ind, row_ptr, x, rows, columns, used_elements);
#endif
        }

#if defined(FEAST_HAVE_QUADMATH) && !defined(__CUDACC__)
        static void csr(__float128 * r, const __float128 * const rhs, const __float128 * const val, const Index * const col_ind, const Index * const row_ptr, const __float128 * const x, const Index rows, const Index columns, const Index used_elements)
        {
          csr_generic(r, rhs, val, col_ind, row_ptr, x, rows, columns, used_elements);
        }
#endif

        template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
        static void csrb(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
#ifdef FEAST_BACKENDS_MKL
          if (typeid(IT_) == typeid(unsigned long) && BlockHeight_ == BlockWidth_)
            csrb_mkl(r, rhs, val, (unsigned long*)col_ind, (unsigned long*)row_ptr, x, rows, columns, used_elements, BlockHeight_);
          else
            csrb_generic<DT_, IT_, BlockHeight_, BlockWidth_>(r, rhs, val, col_ind, row_ptr, x, rows, columns, used_elements);
#else
          csrb_generic<DT_, IT_, BlockHeight_, BlockWidth_>(r, rhs, val, col_ind, row_ptr, x, rows, columns, used_elements);
#endif
        }

#if defined(FEAST_HAVE_QUADMATH) && !defined(__CUDACC__)
        template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
        static void csrb(__float128 * r, const __float128 * const rhs, const __float128 * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const __float128 * const x, const Index rows, const Index columns, const Index used_elements)
        {
          csrb_generic<__float128, IT_, BlockHeight_, BlockWidth_>(r, rhs, val, col_ind, row_ptr, x, rows, columns, used_elements);
        }
#endif

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs, const IT_ * const cl, const DT_ * const x, const Index C, const Index rows)
        {
          ell_generic(r, rhs, val, col_ind, cs, cl, x, C, rows);
        }

        template <typename DT_, typename IT_>
        static void coo(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const row_ptr, const IT_ * const col_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
          coo_generic(r, rhs, val, row_ptr, col_ptr, x, rows, columns, used_elements);
        }

        template <typename DT_>
        static void coo(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
#ifdef FEAST_BACKENDS_MKL
          coo_mkl(r, rhs, val, row_ptr, col_ptr, x, rows, columns, used_elements);
#else
          coo_generic(r, rhs, val, row_ptr, col_ptr, x, rows, columns, used_elements);
#endif
        }

#if defined(FEAST_HAVE_QUADMATH) && !defined(__CUDACC__)
        static void coo(__float128 * r, const __float128 * const rhs, const __float128 * const val, const Index * const row_ptr, const Index * const col_ptr, const __float128 * const x, const Index rows, const Index columns, const Index used_elements)
        {
          coo_generic(r, rhs, val, row_ptr, col_ptr, x, rows, columns, used_elements);
        }
#endif

        template <typename DT_, typename IT_>
        static void banded(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns)
        {
          banded_generic(r, rhs, val, offsets, x, num_of_offsets, rows, columns);
        }

        template <typename DT_>
        static void dense(DT_ * r, const DT_ * const rhs, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
        {
          dense_generic(r, rhs, val, x, rows, columns);
        }

        template <typename DT_, typename IT_>
        static void csr_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
        static void csrb_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_>
        static void ell_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs, const IT_ * const cl, const DT_ * const x, const Index C, const Index rows);

        template <typename DT_, typename IT_>
        static void coo_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const row_ptr, const IT_ * const col_ptr, const DT_ * const x, const Index rows, const Index, const Index used_elements);

        template <typename DT_, typename IT_>
        static void banded_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns);

        template <typename DT_>
        static void dense_generic(DT_ * r, const DT_ * const rhs, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns);

        static void csr_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const float * const x, const Index rows, const Index columns, const Index);
        static void csr_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const double * const x, const Index rows, const Index columns, const Index);

        static void csrb_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const float * const x, const Index rows, const Index columns, const Index, const int blocksize);
        static void csrb_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const double * const x, const Index rows, const Index columns, const Index, const int blocksize);

        static void coo_mkl(float * r, const float * const rhs, const float * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const float * const x, const Index rows, const Index columns, const Index used_elements);
        static void coo_mkl(double * r, const double * const rhs, const double * const val, const unsigned long * const row_ptr, const unsigned long * const col_ptr, const double * const x, const Index rows, const Index columns, const Index used_elements);
      };

      extern template void Defect<Mem::Main>::csr_generic(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main>::csr_generic(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main>::csr_generic(double *, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main>::csr_generic(double *, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void Defect<Mem::Main>::ell_generic(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      extern template void Defect<Mem::Main>::ell_generic(double *, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      extern template void Defect<Mem::Main>::ell_generic(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      extern template void Defect<Mem::Main>::ell_generic(double *, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      extern template void Defect<Mem::Main>::coo_generic(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main>::coo_generic(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main>::coo_generic(double *, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main>::coo_generic(double *, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      extern template void Defect<Mem::Main>::banded_generic(float *, const float * const, const float * const, const Index * const, const float * const, const Index, const Index, const Index);
      extern template void Defect<Mem::Main>::banded_generic(double *, const double * const, const double * const, const Index * const, const double * const, const Index, const Index, const Index);

      extern template void Defect<Mem::Main>::dense_generic(float *, const float * const, const float * const, const float * const, const Index, const Index);
      extern template void Defect<Mem::Main>::dense_generic(double *, const double * const, const double * const, const double * const, const Index, const Index);

      template <>
      struct Defect<Mem::CUDA>
      {
        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_>
        static void csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index);

        template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
        static void csrb(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
        {
          static_assert(std::is_same<IT_, unsigned int>::value, "cuda bcsr only supports unsigned int!");
          static_assert(BlockHeight_ == BlockWidth_, "cuda bcsr only supports squared matrix blocks!");
          csrb_intern(r, rhs, val, col_ind, row_ptr, x, rows, columns, used_elements, BlockHeight_);
        }

        template <typename DT_>
        static void csrb_intern(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index, const int blocksize);

        template <typename DT_, typename IT_>
        static void ell(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs, const IT_ * const cl, const DT_ * const x, const Index C, const Index rows);

        template <typename DT_, typename IT_>
        static void banded(DT_ * r, const DT_ * const rhs, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/defect_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_DEFECT_HPP
