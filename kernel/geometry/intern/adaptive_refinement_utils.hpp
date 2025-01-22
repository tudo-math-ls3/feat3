// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <array>

namespace FEAT::Geometry::Intern
{
  /**
   * \brief Returns coefficients for vertex interpolation
   *
   * \tparam Shape_ Shape to interpolate in
   * \param[in] vertex_idx Vertex for which interpolation coefficients should be returned
   */
  template<typename Shape_>
  static std::array<Real, Shape::FaceTraits<Shape_, 0>::count>
  vertex_coefficients(const Tiny::Vector<Real, Shape_::dimension>& vertex)
  {
    static_assert(Shape_::dimension >= 1, "Can't determine vertex coefficients for shape of this dimension!");
    static_assert(Shape_::dimension <= 3, "Can't determine vertex coefficients for shape of this dimension!");

    static constexpr int size = Shape::FaceTraits<Shape_, 0>::count;

    if constexpr(Shape_::dimension == 1)
    {
      const Real x = vertex[0];

      return std::array<Real, size>{
        1.0 - x,
        x,
      };
    }

    if constexpr(Shape_::dimension == 2)
    {
      const Real x = vertex[0];
      const Real y = vertex[1];
      const Real xy = x * y;

      return std::array<Real, size>{1 - x - y + xy, x - xy, y - xy, xy};
    }

    if constexpr(Shape_::dimension == 3)
    {
      const Real x = vertex[0];
      const Real y = vertex[1];
      const Real z = vertex[2];
      const Real xy = x * y;
      const Real xz = x * z;
      const Real yz = y * z;
      const Real xyz = x * y * z;

      return std::array<Real, size>{
        1 - x - y - z + xy + xz + yz - xyz,
        x - xy - xz + xyz,
        y - xy - yz + xyz,
        xy - xyz,
        z - xz - yz + xyz,
        xz - xyz,
        yz - xyz,
        xyz};
    }
  }

  /**
   * \brief Linear interpolation for vertices
   *
   * \param[in] vertices Vertex coordinates
   * \param[in] coeffs Interpolation coefficients
   */
  template<typename VertexType_, std::size_t num_verts_>
  VertexType_ static interpolate(
    const std::array<VertexType_, num_verts_>& vertices,
    const std::array<Real, num_verts_>& coeffs)
  {
    using CoordType = typename VertexType_::ValueType;
    VertexType_ result(0);
    for(Index i = 0; i < num_verts_; i++)
    {
      result += CoordType(coeffs[i]) * vertices[i];
    }
    return result;
  }

  template<typename ArrayLike_>
  static ArrayLike_ rotate_arraylike_2d(const ArrayLike_& input)
  {
    ArrayLike_ result;

    result[0] = input[2];
    result[1] = input[0];
    result[2] = input[3];
    result[3] = input[1];

    return result;
  }

  template<typename ArrayLike_>
  static ArrayLike_ rotate_arraylike_xaxis(const ArrayLike_& input)
  {
    ArrayLike_ result;

    result[0] = input[4];
    result[1] = input[5];
    result[2] = input[0];
    result[3] = input[1];
    result[4] = input[6];
    result[5] = input[7];
    result[6] = input[2];
    result[7] = input[3];

    return result;
  }

  template<typename ArrayLike_>
  static ArrayLike_ rotate_arraylike_yaxis(const ArrayLike_& input)
  {
    ArrayLike_ result;

    result[0] = input[1];
    result[1] = input[5];
    result[2] = input[3];
    result[3] = input[7];
    result[4] = input[0];
    result[5] = input[4];
    result[6] = input[2];
    result[7] = input[6];

    return result;
  }

  template<typename ArrayLike_>
  static ArrayLike_ rotate_arraylike_zaxis(const ArrayLike_& input)
  {
    ArrayLike_ result;

    result[0] = input[2];
    result[1] = input[0];
    result[2] = input[3];
    result[3] = input[1];
    result[4] = input[6];
    result[5] = input[4];
    result[6] = input[7];
    result[7] = input[5];

    return result;
  }
} // namespace FEAT::Geometry::Intern
