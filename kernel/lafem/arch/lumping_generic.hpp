#pragma once
#ifndef KERNEL_LAFEM_ARCH_LUMPING_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_LUMPING_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_LUMPING_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_>
      void Lumping<Mem::Main>::csr_generic(DT_ * lump, const DT_ * const val, const IT_ * const /*col_ind*/,
        const IT_ * const row_ptr, const Index rows)
      {
        for (Index row(0); row < rows; row++)
        {
          Index end = row_ptr[row + 1];
          lump[row] = DT_(0);
          for (Index col = row_ptr[row]; col < end; col++)
          {
            lump[row] += val[col];
          }
        }
      }

      template <typename DT_, typename IT_>
      void Lumping<Mem::Main>::ell_generic(DT_ * lump, const DT_ * const val, const IT_ * const /*col_ind*/,
        const IT_ * const cs, const IT_ * const /*cl*/, const Index C, const Index rows)
      {
        for (Index row(0) ; row < rows ; ++row)
        {
          const Index chunk(row / C);
          const Index local_row(row % C);
          const Index chunk_end(cs[chunk+1]);

          lump[row] = DT_(0);
          for (Index pcol(cs[chunk] + local_row) ; pcol < chunk_end ; pcol+=C)
          {
            lump[row] += val[pcol];
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_LUMPING_GENERIC_HPP
