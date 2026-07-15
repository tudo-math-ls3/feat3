// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/solver/base.hpp>

#include <array>
#include <vector>

namespace FEAT::Solver::RichardsonHelper
{
  // All matrix/vector/mirror data is passed as Memory::Arbiter objects; the device driver acquires
  // the cuda views it needs itself, instead of receiving pre-resolved raw device pointers. This
  // mirrors the design of FEAT::Solver::CSRHelper used by the fused AmaVankaGMRES driver, but keeps
  // the meta-block (k x l) BCSR structure intact instead of flattening it into a single scalar CSR.

  /**
   * \brief A single (k x l) BCSR meta-block, referenced by its data arbiters.
   *
   * A grid cell may be "empty" (\c val == nullptr or \c num_nnz == 0), in which case the host driver
   * skips it entirely (no kernel launch). This is how the structurally-zero blocks of a saddle-point
   * system (and the numerically-empty pressure-pressure block of the AmaVanka inverse) cost nothing,
   * without the kernel ever special-casing the matrix structure.
   *
   * The size fields (\c num_rows, \c num_cols, \c num_nnz) are *block* counts, matching the
   * LAFEM::SparseMatrixBCSR num_rows()/num_cols()/num_nzes() semantics expected by
   * LAFEM::Arch::MatVecMultBCSRDense::exec_cuda_impl.
   */
  template<typename DT_, typename IT_>
  struct MetaBlock
  {
    const Memory::Arbiter* val = nullptr;
    const Memory::Arbiter* row = nullptr;
    const Memory::Arbiter* col = nullptr;
    /// persistent cuSPARSE SpMV/BSRMV workspace of this block (owned by the LAFEM matrix), or nullptr
    /// if the chosen SpMV path needs none; passed through to MatVecMultBCSRDense::exec_cuda_impl
    void* workspace = nullptr;
    int block_height = 0;
    int block_width = 0;
    Index num_rows = 0;
    Index num_cols = 0;
    Index num_nnz = 0;

    bool empty() const { return (val == nullptr) || (num_nnz == Index(0)); }
  };

  /**
   * \brief A whole meta matrix as a constexpr NR_ x NC_ grid of (possibly empty) BCSR blocks.
   *
   * The grid dimensions are compile-time known (the meta layout is fixed), so the container uses a
   * std::array. Blocks are stored row-major: \c blocks[r*NC_ + c] is the (r,c) meta-block, which maps
   * source vector component \c c to target vector component \c r.
   */
  template<typename DT_, typename IT_, int NR_, int NC_>
  struct MetaMatrix
  {
    static constexpr int num_row_blocks = NR_;
    static constexpr int num_col_blocks = NC_;
    std::array<MetaBlock<DT_, IT_>, std::size_t(NR_) * std::size_t(NC_)> blocks;

    MetaBlock<DT_, IT_>& at(int r, int c) { return blocks[std::size_t(r) * std::size_t(NC_) + std::size_t(c)]; }
    const MetaBlock<DT_, IT_>& at(int r, int c) const { return blocks[std::size_t(r) * std::size_t(NC_) + std::size_t(c)]; }
  };

  /// \brief Layout of one meta-vector component (velocity => block_size == dim, pressure => 1).
  struct MetaVectorInfo
  {
    int block_size = 1;
    /// number of *blocks* (nodes); the raw scalar length is num_entries * block_size
    Index num_entries = 0;

    Index size_raw() const { return num_entries * Index(block_size); }
  };

  /// \brief Flat unit-filter (defect/correction) info for a single meta component (phase 1: only
  /// UnitFilter / UnitFilterBlocked DOFs are collected; Slip/Mean filters are dropped).
  template<typename DT_, typename IT_>
  struct UnitFilterInfo
  {
    const Memory::Arbiter* val = nullptr;
    const Memory::Arbiter* ind = nullptr;
    Index used_elements = 0;
  };

  /**
   * \brief Per-meta-row (per-component) gate for the hand-rolled synchronization.
   *
   * One GateWrapper per meta component/row. A component whose \c ranks is empty needs no
   * synchronization (e.g. a discontinuous pressure space) and its sync is a no-op — which lets the
   * driver overlap the halo exchange of the components that do need it with the compute of those that
   * do not. All components of a saddle-point system share the same scalar shared-node mirror indices,
   * but each carries its own here to keep the interface fully general.
   */
  template<typename DT_, typename IT_>
  struct GateWrapper
  {
    /// per-neighbor shared-node (scalar) mirror index arrays
    std::vector<const Memory::Arbiter*> mirror_idx;
    /// per-neighbor number of shared nodes (mirror indices)
    std::vector<Index> mirror_size;
    /// neighbor ranks (empty => empty gate)
    std::vector<int> ranks;
    /// this component's frequency array (type-1 -> type-0 weighting for sync_1)
    const Memory::Arbiter* freqs = nullptr;
    Dist::Comm* comm = nullptr;

    bool empty() const { return ranks.empty(); }
  };

  /**
   * \brief Everything the fused Richardson driver needs for one apply.
   *
   * \tparam NC_ number of meta components (== number of meta rows/cols of the square system).
   */
  template<typename DT_, typename IT_, int NC_>
  struct RichardsonApplyInfo
  {
    /// system matrix A (structure-agnostic; empty cells skipped)
    MetaMatrix<DT_, IT_, NC_, NC_> mat_sys;
    /// AmaVanka approximate-inverse matrix V
    MetaMatrix<DT_, IT_, NC_, NC_> mat_vanka;
    /// one defect/correction unit filter per component
    std::array<UnitFilterInfo<DT_, IT_>, std::size_t(NC_)> def_filters;
    /// per-component layout (block size + node count)
    std::array<MetaVectorInfo, std::size_t(NC_)> components;
    /// one gate per meta row (some may be empty)
    std::array<GateWrapper<DT_, IT_>, std::size_t(NC_)> gates;

    /// per-component element arbiters of the solution vector (read-write, iterated in place)
    std::array<Memory::Arbiter*, std::size_t(NC_)> vec_sol;
    /// per-component element arbiters of the right-hand-side vector (read-only)
    std::array<const Memory::Arbiter*, std::size_t(NC_)> vec_rhs;
    /// per-component element arbiters of the defect scratch vector
    std::array<Memory::Arbiter*, std::size_t(NC_)> vec_def;
    /// per-component element arbiters of the correction scratch vector
    std::array<Memory::Arbiter*, std::size_t(NC_)> vec_cor;

    DT_ omega = DT_(1);
    Index min_iter = 0;
    Index max_iter = 0;
    bool always_use_parallel_kernel = false;
  };

  template<typename DT_>
  struct RichardsonStatus
  {
    DT_ def_norm = DT_(0);
    DT_ def_init = DT_(0);
    Index num_iter = 0;
  };
} // namespace FEAT::Solver::RichardsonHelper
