// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    //forward declarations
    template <typename DT_, typename IT_>
    class DenseVector;

    template <typename DT_, typename IT_>
    class SparseVector;

    template <typename DT_, typename IT_, int BlockSize_>
    class DenseVectorBlocked;

    template <typename DT_, typename IT_, int BlockSize_>
    class SparseVectorBlocked;

    template <typename DT_, typename IT_>
    class DenseMatrix;

    template <typename DT_, typename IT_>
    class SparseMatrixCSR;

    template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
    class SparseMatrixBCSR;

    template <typename DT_, typename IT_>
    class SparseMatrixBanded;

    template <typename DT_, typename IT_>
    class SparseMatrixCSCR;

    template<typename DT_, typename IT_>
    class VectorMirror;

    template<typename DT_, typename IT_>
    class MatrixMirror;

  } // namespace LAFEM
} // namespace FEAT
