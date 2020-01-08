// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/solver/solver_factory.hpp>

/// \compilerhack ICC < 18 fail to link this with undefined reference to `__must_be_linked_with_icc_or_xild' error
#if not (defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1800))
namespace FEAT
{
  namespace Solver
  {
    using MST1_ = Solver::MatrixStock<
      Global::Matrix<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Filter<LAFEM::UnitFilter<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Transfer<LAFEM::Transfer<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>>, LAFEM::VectorMirror<Mem::Main, double, Index>>>;

    template std::shared_ptr<SolverBase<MST1_::VectorType>> SolverFactory::create_scalar_solver_by_section<
      MST1_,
      MST1_::VectorType>
      (MST1_&, PropertyMap*, const String&, size_t);

    template std::shared_ptr<SolverBase<MST1_::VectorType::LocalVectorType>> SolverFactory::create_scalar_solver_by_section<
      MST1_,
      MST1_::VectorType::LocalVectorType>
      (MST1_&, PropertyMap*, const String&, size_t);
  }
}
#endif //FEAT_COMPILER_INTEL
