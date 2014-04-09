// includes, FEAST
#include <kernel/lafem/arch/scale_row_col.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template void ScaleRows<Mem::Main, Algo::Generic>::csr(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index, const Index);
      template void ScaleRows<Mem::Main, Algo::Generic>::csr(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index, const Index);

      template void ScaleCols<Mem::Main, Algo::Generic>::csr(float *, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index, const Index);
      template void ScaleCols<Mem::Main, Algo::Generic>::csr(double *, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index, const Index);
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST
