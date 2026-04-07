// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>

// includes, system
#include <iostream>

namespace FEAT
{
  /// The backend that shall be used in all compute heavy calculations
  enum class PreferredBackend
  {
    generic = 0, /**< generic c++ code */
    mkl, /**< intel mkl blas library */
    cuda /**< nvidia cuda gpgpu support */
  };

  /**
   * \brief Backend support class
   *
   * \author Dirk Ribbrock, Peter Zajac
   */
  class Backend
  {
  private:
    /// the currently preferred backend
    static PreferredBackend _preferred_backend;

    /// missing backend warnings enabled/disabled
    static bool _enable_warn_missing;

  public:
    /// set new preferred backend
    static void set_preferred_backend(PreferredBackend preferred_backend);

    /// get current preferred backend
    static PreferredBackend get_preferred_backend();

    /// enables or disables missing backend warning
    static void set_warn_missing(bool enabled);

    /// issues a warning about unsupported backend usage
    static void warn_missing(const char* operation, PreferredBackend preferred_backend, PreferredBackend fallback_backend = PreferredBackend::generic);

    /// issues a warning about unsupported backend usage
    static void warn_missing(const char* operation)
    {
      warn_missing(operation, get_preferred_backend());
    }
  }; // class Backend

  std::ostream & operator<< (std::ostream & left, FEAT::PreferredBackend backend);
} // namespace FEAT
