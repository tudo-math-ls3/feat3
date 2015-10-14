#pragma once
#ifndef KERNEL_LAFEM_BASE_HPP
#define KERNEL_LAFEM_BASE_HPP 1

/**
 * \file
 * \brief LAFEM common type definitions.
 *
 * This file is the supplementary header for the LAFEM library kernel, which is included by all other FEAST header and source files.
 * It defines macros and data types which are frequently used in other files.
 */

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  /**
   * \brief LAFEM namespace
   */
  namespace LAFEM
  {
    /**
     * Supported File modes.
     */
    enum class FileMode
    {
      fm_exp = 0, /**< Exponential ascii */
      fm_dv, /**< Binary data */
      fm_mtx, /**< Matrix market ascii */
      fm_ell, /**< Binary ell data */
      fm_csr, /**< Binary csr data */
      fm_coo, /**< Binary coo data */
      fm_bm, /**< Binary banded data */
      fm_dm,  /**< Binary dense matrix data */
      fm_sv  /**< Binary sparse vector data */
    };

    /**
     * Supported clone modes.
     */
    enum class CloneMode
    {
      Shallow = 0, /**< Share index and data arrays */
      Layout, /**< Share index arrays, allocate new data array */
      Weak, /**< Share index arrays, allocate new data array and copy content */
      Deep /**< Allocate new index and data arrays and copy content */
    };

    /**
     * Supported layout ids.
     */
    enum class SparseLayoutId
    {
      lt_csr = 0,
      lt_coo,
      lt_ell,
      lt_banded
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_BASE_HPP
