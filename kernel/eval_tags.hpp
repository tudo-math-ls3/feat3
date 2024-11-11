// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>

namespace FEAT
{
  /**
   * \brief Trafo configuration tags enum
   *
   * This enumeration is used to specify the evaluation requirements for
   * transformation evaluators.
   * For this enumeration, the "binary" AND and OR operators are overloaded,
   * so that the single enumeration values can be combined and masked
   * using these operators.
   */
  enum class TrafoTags
  {
    none      = 0,
    /// specifies whether the trafo should supply domain point coordinates
    dom_point = 0x0001,
    /// specifies whether the trafo should supply image point coordinates
    img_point = 0x0002,
    /// specifies whether the trafo should supply jacobian matrices
    jac_mat   = 0x0004,
    /// specifies whether the trafo should supply inverse jacobian matrices
    jac_inv   = 0x0008,
    /// specifies whether the trafo should supply jacobian determinants
    jac_det   = 0x0010,
    /// specifies whether the trafo should supply hessian tensors
    hess_ten  = 0x0020,
    /// specifies whether the trafo should supply inverse hessian tensors
    hess_inv  = 0x0040,
    /// specifies whether the trafo should supply outer normal vectors
    /// note: The normal vector is always automatically defined by the Assembly::TraceAssembler
    /// class and does not have to be specifies explicitly
    normal    = 0x0080
  };

  /**
   * \brief Binary OR operator for TrafoTags values.
   *
   * \param[in] a, b
   * The two TrafoTags values that are to be OR'ed.
   *
   * \returns
   * The combined tags of \p a and \p b.
   */
  static constexpr inline TrafoTags operator|(TrafoTags a, TrafoTags b)
  {
    return static_cast<TrafoTags>(static_cast<int>(a) | static_cast<int>(b));
  }

  /**
   * \brief Binary AND operator for TrafoTags values.
   *
   * \param[in] a, b
   * The two TrafoTags values that are to be AND'ed.
   *
   * \returns
   * The masked tags of \p a and \p b.
   */
  static constexpr inline TrafoTags operator&(TrafoTags a, TrafoTags b)
  {
    return static_cast<TrafoTags>(static_cast<int>(a) & static_cast<int>(b));
  }

  /**
   * \brief bool conversion operator
   *
   * \param[in] a
   * The TrafoTags value that is to be checked.
   *
   * \returns
   * \c false, if \p a is TrafoTags::none, otherwise \c true.
   */
  static constexpr inline bool operator*(TrafoTags a)
  {
    // just check the interesting bits
    return (static_cast<int>(a) & 0xFF) != 0;
  }

  /**
   * \brief Space configuration tags enum
   *
   * This enumeration is used to specify the evaluation requirements for
   * space evaluators.
   * For this enumeration, the "binary" AND and OR operators are overloaded,
   * so that the single enumeration values can be combined and masked
   * using these operators.
   */
  enum class SpaceTags
  {
    none      = 0,
    /// specifies whether the space should supply basis function values
    value     = 0x0001,
    /// specifies whether the space should supply basis function gradients
    grad      = 0x0002,
    /// specifies whether the space should supply basis function hessians
    hess      = 0x0004,
    /// specifies whether the space should supply reference basis function values
    ref_value = 0x0008,
    /// specifies whether the space should supply reference basis function gradients
    ref_grad  = 0x0010,
    /// specifies whether the space should supply reference basis function hessians
    ref_hess  = 0x0020
  };

  /**
   * \brief Binary OR operator for SpaceTags values.
   *
   * \param[in] a, b
   * The two SpaceTags values that are to be OR'ed.
   *
   * \returns
   * The combined tags of \p a and \p b.
   */
  static constexpr inline SpaceTags operator|(SpaceTags a, SpaceTags b)
  {
    return (SpaceTags)(static_cast<int>(a) | static_cast<int>(b));
  }

  /**
   * \brief Binary AND operator for SpaceTags values.
   *
   * \param[in] a, b
   * The two SpaceTags values that are to be AND'ed.
   *
   * \returns
   * The masked tags of \p a and \p b.
   */
  static constexpr inline SpaceTags operator&(SpaceTags a, SpaceTags b)
  {
    return (SpaceTags)(static_cast<int>(a) & static_cast<int>(b));
  }

  /**
   * \brief bool conversion operator
   *
   * \param[in] a
   * The SpaceTags value that is to be checked.
   *
   * \returns
   * \c false, if \p a is SpaceTags::none, otherwise \c true.
   */
  static constexpr inline bool operator*(SpaceTags a)
  {
    // just check the interesting bits
    return (static_cast<int>(a) & 0x3F) != 0;
  }
} // namespace FEAT
