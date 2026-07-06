// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/omp_util.hpp>
#include <kernel/backend.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/solver/iterative.hpp>



namespace FEAT::Solver::CSRHelper
{
  // The matrix/filter/vector data is passed as Memory::Arbiter objects; the kernels acquire the
  // device (cuda) views they need themselves, instead of receiving pre-resolved raw device pointers.
  template<typename DT_, typename IT_>
  struct CRSMat
  {
    const Memory::Arbiter* val;
    const Memory::Arbiter* row;
    const Memory::Arbiter* col;
    Index num_rows;
    Index num_cols;
    Index num_nnz;
  };

  template<typename DT_, typename IT_>
  struct UnitFilterInfo
  {
    const Memory::Arbiter* val;
    const Memory::Arbiter* ind;
    Index used_elements;
  };

  template<typename DT_, typename IT_>
  struct GMRESApplyInfo
  {
    CRSMat<DT_, IT_> mat_sys;
    CRSMat<DT_, IT_> mat_vanka;
    UnitFilterInfo<DT_, IT_> def_filter;
    // krylov basis vectors: one arbiter per vector (read-write inside the kernel)
    const std::vector<Memory::Arbiter*>* vec_v;
    Memory::Arbiter* vec_w;
    Global::Gate<LAFEM::DenseVector<DT_, IT_>, LAFEM::VectorMirror<DT_, IT_>>* gate;
    Index size;

    DT_ res_scale;
    DT_ abs_tol;
    DT_ rel_tol;
    Index krylov_dim;
    Index max_iter;
    Index min_iter;
    Index plot_interval;
    FEAT::Solver::PlotMode plot_mode;
    bool prefetch_data;
  };

  template<typename DT_>
  struct GMRESStatus
  {
    DT_ def_norm;
    DT_ def_init;
    Index num_iter;
  };
}
