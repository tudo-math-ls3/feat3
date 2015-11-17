#include <kernel/util/kahan_summation.hpp>

using namespace FEAST;

KahanAccumulation FEAST::KahanSum(KahanAccumulation accumulation, double value)
{
  KahanAccumulation result;
  double y = value - accumulation.correction;
  double t = accumulation.sum + y;
  result.correction = (t - accumulation.sum) - y;
  result.sum = t;
  return result;
}
