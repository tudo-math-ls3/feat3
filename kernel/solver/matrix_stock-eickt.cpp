// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/solver/matrix_stock.hpp>

namespace FEAT
{
  namespace Solver
  {
    template class MatrixStock<
      Global::Matrix<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Filter<LAFEM::UnitFilter<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      //Global::Matrix<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>
      Global::Transfer<LAFEM::Transfer<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>>, LAFEM::VectorMirror<Mem::Main, double, Index>>
        >;

    template class MatrixStock<
      Global::Matrix<LAFEM::SparseMatrixBCSR<Mem::Main, double, Index, 2, 2>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Filter<LAFEM::UnitFilterBlocked<Mem::Main, double, Index, 2>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      //Global::Matrix<LAFEM::SparseMatrixBCSR<Mem::Main, double, Index, 2, 2>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>
      Global::Transfer<LAFEM::Transfer<LAFEM::SparseMatrixBCSR<Mem::Main, double, Index, 2, 2>>, LAFEM::VectorMirror<Mem::Main, double, Index>>
        >;
  }
}
