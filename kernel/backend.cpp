// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/backend.hpp>
#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

namespace FEAT
{
  /// the static variable for the preferred backend
  PreferredBackend Backend::_preferred_backend = PreferredBackend::generic;

  void Backend::set_preferred_backend(PreferredBackend preferred_backend)
  {
    _preferred_backend = preferred_backend;
    if (preferred_backend != PreferredBackend::cuda)
    {
#ifdef FEAT_HAVE_CUDA
      FEAT::Util::cuda_synchronize();
#endif
    }
  }

  PreferredBackend Backend::get_preferred_backend()
  {
    return _preferred_backend;
  }

  std::ostream & operator<< (std::ostream & left, PreferredBackend value)
  {
    switch (value)
    {
    case PreferredBackend::generic:
      left << "generic";
      break;

    case PreferredBackend::mkl:
      left << "mkl";
      break;

    case PreferredBackend::cuda:
      left << "cuda";
      break;

    default:
      left << "unknown preferred backend";
      break;
    }

    return left;
  }

} // namespace FEAT
