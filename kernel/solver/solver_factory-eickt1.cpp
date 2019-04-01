// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/solver/solver_factory.hpp>

namespace FEAT
{
  namespace Solver
  {
    using MST1_ = Solver::MatrixStock<
      Global::Matrix<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Filter<LAFEM::UnitFilter<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Transfer<LAFEM::Transfer<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>>, LAFEM::VectorMirror<Mem::Main, double, Index>>>;

    template std::shared_ptr<SolverBase<MST1_::VectorType::LocalVectorType>>  SolverFactory::create_scalar_solver<
      MST1_,
      MST1_::VectorType::LocalVectorType>
      (MST1_&, PropertyMap*, const String&, std::size_t);

    template std::shared_ptr<SolverBase<MST1_::VectorType>>  SolverFactory::create_scalar_solver<
      MST1_,
      MST1_::VectorType>
      (MST1_&, PropertyMap*, const String&, std::size_t);
  }
}
