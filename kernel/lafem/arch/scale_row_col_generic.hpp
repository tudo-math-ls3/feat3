#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {

      template <typename DT_, typename IT_>
      void ScaleRows<Mem::Main>::csr_generic(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/,
                                                    const IT_ * const row_ptr, const DT_ * const x,
                                                    const Index rows, const Index, const Index)
      {
        for (Index row(0) ; row < rows ; ++row)
        {
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i(row_ptr[row]) ; i < end ; ++i)
          {
            r[i] = a[i] * x[row];
          }
        }
      }

      template <typename DT_, typename IT_>
      void ScaleRows<Mem::Main>::coo_generic(DT_ * r, const DT_ * const a, const IT_ * const row_idx,
                                                    const IT_ * const /*col_idx*/, const DT_ * const x,
                                                    const Index, const Index, const Index used_elements)
      {
        for (Index i(0) ; i < used_elements ; ++i)
        {
          r[i] = a[i] * x[row_idx[i]];
        }
      }

      template <typename DT_, typename IT_>
      void ScaleRows<Mem::Main>::ell_generic(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/,
                                                    const IT_ * const cs, const IT_ * const cl, const IT_ * const /*rl*/,
                                                    const DT_ * const x, const Index C, const Index rows)
      {
        for (Index i(0) ; i < rows/C ; ++i)
        {
          for (Index j(0) ; j < cl[i] ; ++j)
          {
            for (Index k(0); k < C; ++k)
            {
              r[cs[i]+j*C+k] = a[cs[i]+j*C+k] * x[i*C+k];
            }
          }
        }

        Index i(rows/C);
        {
          for (Index k(0) ; k < rows%C ; ++k)
          {
            for (Index j(0) ; j < cl[i] ; ++j)
            {
              r[cs[i]+j*C+k] = a[cs[i]+j*C+k] * x[i*C+k];
            }
          }
        }
      }


      // ***********************************************

      template <typename DT_, typename IT_>
      void ScaleCols<Mem::Main>::csr_generic(DT_ * r, const DT_ * const a, const IT_ * const col_ind,
                                                    const IT_ * const row_ptr, const DT_ * const x,
                                                    const Index rows, const Index, const Index)
      {
        for (Index row(0) ; row < rows ; ++row)
        {
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i(row_ptr[row]) ; i < end ; ++i)
          {
            r[i] = a[i] * x[col_ind[i]];
          }
        }
      }

      template <typename DT_, typename IT_>
      void ScaleCols<Mem::Main>::coo_generic(DT_ * r, const DT_ * const a, const IT_ * const /*row_idx*/,
                                                    const IT_ * const col_idx, const DT_ * const x, const Index,
                                                    const Index, const Index used_elements)
      {
        for (Index i(0) ; i < used_elements ; ++i)
        {
          r[i] = a[i] * x[col_idx[i]];
        }
      }

      template <typename DT_, typename IT_>
      void ScaleCols<Mem::Main>::ell_generic(DT_ * r, const DT_ * const a, const IT_ * const col_ind,
                                                    const IT_ * const cs, const IT_ * const cl, const IT_ * const /*rl*/,
                                                    const DT_ * const x, const Index C, const Index rows)
      {
        for (Index i(0) ; i < rows/C ; ++i)
        {
          for (Index j(0) ; j < cl[i] ; ++j)
          {
            for (Index k(0); k < C; ++k)
            {
              r[cs[i]+j*C+k] = a[cs[i]+j*C+k] * x[col_ind[cs[i]+j*C+k]];
            }
          }
        }

        Index i(rows/C);
        {
          for (Index k(0) ; k < rows%C ; ++k)
          {
            for (Index j(0) ; j < cl[i] ; ++j)
            {
              r[cs[i]+j*C+k] = a[cs[i]+j*C+k] * x[col_ind[cs[i]+j*C+k]];
            }
          }
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP
