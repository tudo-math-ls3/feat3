// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
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

  /// do not warn for missing backends by default
  bool Backend::_enable_warn_missing = false;

  void Backend::set_preferred_backend(PreferredBackend preferred_backend)
  {
    // do we switch from CUDA to something else?
#ifdef FEAT_HAVE_CUDA
    if((PreferredBackend::cuda == _preferred_backend) && (preferred_backend != _preferred_backend))
    {
      FEAT::Util::cuda_synchronize();
    }
#endif

    _preferred_backend = preferred_backend;
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

  void Backend::set_warn_missing(bool enabled)
  {
    Backend::_enable_warn_missing = enabled;
  }

  void Backend::warn_missing(const char* operation, PreferredBackend preferred_backend, PreferredBackend fallback_backend)
  {
    if(!Backend::_enable_warn_missing)
      return;

    const char* pbs = "???";
    const char* fbs = "???";

    switch(preferred_backend)
    {
    case PreferredBackend::generic:
      pbs = "generic";
      break;
    case PreferredBackend::mkl:
      pbs = "mkl";
      break;
    case PreferredBackend::cuda:
      pbs = "cuda";
      break;
    }

    switch(fallback_backend)
    {
    case PreferredBackend::generic:
      fbs = "generic";
      break;
    case PreferredBackend::mkl:
      fbs = "mkl";
      break;
    case PreferredBackend::cuda:
      fbs = "cuda";
      break;
    }

    fprintf(stderr, "[WARNING] Operation '%s' misses preferred Backend '%s', using fallback backend '%s' instead\n",
      operation, pbs, fbs);
  }

} // namespace FEAT
