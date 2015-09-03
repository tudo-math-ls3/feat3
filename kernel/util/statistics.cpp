#include <kernel/util/statistics.hpp>

using namespace FEAST;

// static member initialisation
Index Statistics::_flops = Index(0);

void Statistics::add_flops(Index flops)
{
  _flops += flops;
}

Index Statistics::get_flops()
{
  return _flops;
}

String Statistics::get_formated_flops(double seconds)
{
  double flops((double)_flops);
  flops /= seconds;
  flops /= 1000.; // kilo
  flops /= 1000.; // mega
  flops /= 1000.; // giga
  return stringify(flops) + " GFlop/s";
}

void Statistics::reset_flops()
{
  _flops = Index(0);
}
