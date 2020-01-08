// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/lafem/arch/scale_row_col.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template void ScaleRows<Mem::Main>::csr_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      template void ScaleRows<Mem::Main>::csr_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      template void ScaleRows<Mem::Main>::csr_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      template void ScaleRows<Mem::Main>::csr_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      template void ScaleRows<Mem::Main>::coo_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      template void ScaleRows<Mem::Main>::coo_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      template void ScaleRows<Mem::Main>::coo_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      template void ScaleRows<Mem::Main>::coo_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      template void ScaleRows<Mem::Main>::ell_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      template void ScaleRows<Mem::Main>::ell_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      template void ScaleRows<Mem::Main>::ell_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      template void ScaleRows<Mem::Main>::ell_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

      template void ScaleCols<Mem::Main>::csr_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      template void ScaleCols<Mem::Main>::csr_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      template void ScaleCols<Mem::Main>::csr_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      template void ScaleCols<Mem::Main>::csr_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      template void ScaleCols<Mem::Main>::coo_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
      template void ScaleCols<Mem::Main>::coo_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
      template void ScaleCols<Mem::Main>::coo_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
      template void ScaleCols<Mem::Main>::coo_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

      template void ScaleCols<Mem::Main>::ell_generic(float *, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
      template void ScaleCols<Mem::Main>::ell_generic(double *, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
      template void ScaleCols<Mem::Main>::ell_generic(float *, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
      template void ScaleCols<Mem::Main>::ell_generic(double *, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
