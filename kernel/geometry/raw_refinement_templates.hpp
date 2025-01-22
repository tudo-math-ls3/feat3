// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/geometry/intern/adaptive_refinement_utils.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/util/assertion.hpp>

#include <array>
#include <vector>

namespace FEAT::Geometry
{
  template<typename Shape_, int num_coords_ = Shape_::dimension>
  struct RawEntity
  {
    static const constexpr int num_coords = num_coords_;
    static const constexpr int num_vertices = Shape::FaceTraits<Shape_, 0>::count;

    using VertexType = Tiny::Vector<Real, num_coords>;

    std::array<Tiny::Vector<Real, num_coords>, num_vertices> coords;

    template<typename... Vector>
    explicit RawEntity(Vector... vectors) : coords{vectors...}
    {
    }

    template<int dim_>
    RawEntity<typename Shape::FaceTraits<Shape_, dim_>::ShapeType, num_coords_> face(int idx) const
    {
      using Mapping = Intern::FaceIndexMapping<Shape_, dim_, 0>;
      using SubShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      static constexpr const int num_verts = Shape::FaceTraits<SubShape, 0>::count;

      RawEntity<SubShape, num_coords_> result;

      for(int i(0); i < num_verts; ++i)
      {
        result.coords.at(i) = coords[Mapping::map(idx, i)];
      }

      return result;
    }
  };

  template<typename Shape_, int num_coords_>
  std::ostream& operator<<(std::ostream& stream, const RawEntity<Shape_, num_coords_>& entity)
  {
    stream << "RawEntity<" << Shape_::name() << ", " << stringify(num_coords_) << "> { coords: [ ";
    for(const auto& vertex : entity.coords)
    {
      stream << stringify(vertex) << ", ";
    }
    stream << "]}";
    return stream;
  }

  template<typename Shape_>
  struct RawTemplate
  {
    using EntityType = RawEntity<Shape_>;

    using VertexType = typename EntityType::VertexType;

    std::vector<RawEntity<Shape_>> entities;

    template<typename... Vector>
    RawTemplate& add_entity(Vector... vectors)
    {
      entities.emplace_back(vectors...);
      return *this;
    }

    RawTemplate& add_entity(const std::array<VertexType, Shape::FaceTraits<Shape_, 0>::count>& vertices)
    {
      RawEntity<Shape_> entity;
      entity.coords = vertices;
      entities.push_back(entity);

      return *this;
    }

    RawTemplate& axis_aligned(const VertexType& a, const VertexType& b)
    {
      static_assert(Shape_::dimension <= 3);

      static constexpr int num_verts = Shape::FaceTraits<Shape_, 0>::count;

      static std::array<std::array<Real, 3>, 8> coeffs {{
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1},
      }};

      VertexType diff = b - a;

      std::array<VertexType, num_verts> vertices;
      for(int i(0); i < num_verts; ++i)
      {
        auto& c = coeffs[i];
        VertexType v = a;
        for(int j(0); j < Shape_::dimension; ++j)
        {
          v[j] += c[j] * diff[j];
        }
        for (int j(0); j < v.n; ++j)
        {
          ASSERTM(v[j] >= 0.0 && v[j] <= 1.0, "Invalid vertex in axis_aligned");
        }
        vertices[i] = v;
      }
      add_entity(vertices);

      return *this;
    }

    RawTemplate& grid(std::array<Index, Shape_::dimension> size, VertexType stepsize)
    {
      std::array<Index, Shape_::dimension> coords = {};
      _grid(size, stepsize, 0, coords);
      return *this;
    }

    /**
     * \brief Recursively applies a template to a face of this template
     *
     * Uses linear interpolation to map all entities of the given template into
     * the chosen face. The chosen face is removed from this template.
     */
    RawTemplate& recurse(Index face, const RawTemplate& tmplt)
    {
      const EntityType parent = entities[face];
      entities.erase(entities.begin() + face);

      for(const EntityType& entity : tmplt.entities)
      {
        EntityType mapped{};

        for(Index i(0); i < entity.coords.size(); ++i)
        {
          const auto coeffs = Intern::vertex_coefficients<Shape_>(entity.coords[i]);
          mapped.coords[i] = Intern::interpolate(parent.coords, coeffs);
        }
        entities.push_back(mapped);
      }

      return *this;
    }

  private:
      void _grid(
          const std::array<Index, Shape_::dimension>& size,
          VertexType stepsize,
          int dim,
          std::array<Index, Shape_::dimension>& coords)
      {
        if (dim == Shape_::dimension)
        {
          VertexType v;
          for (int i = 0; i < Shape_::dimension; ++i)
          {
            v[i] = coords[i] * stepsize[i];
          }
          axis_aligned(v, v + stepsize);
        }
        else
        {
          for(Index i(0); i < size[dim]; i++)
          {
            coords[dim] = i;
            _grid(size, stepsize, dim + 1, coords);
          }
        }
      }
  };


  template<typename Shape_, typename FnVert, typename FnEntity>
  static RawTemplate<Shape_> transform_template_vertices(FnVert vertex_transform, FnEntity entity_transform, RawTemplate<Shape_>& tmplt)
  {

    RawTemplate<Shape_> result;
    for(auto& entity : tmplt.entities)
    {
      RawEntity<Shape_> new_entity;

      for(int i(0); i < entity.num_vertices; ++i)
      {
        new_entity.coords.at(i) = vertex_transform(entity.coords.at(i));
      }

      entity_transform(new_entity);
      result.entities.push_back(new_entity);
    }
    return result;
  }

  /**
   * \brief Rotates a raw template 90 degrees counterclockwise around (0.5, 0.5)
   */
  template<typename Shape_>
  static RawTemplate<Shape_> rotate_template_2d(RawTemplate<Shape_>& tmplt)
  {
    static_assert(Shape_::dimension == 2, "rotate_template_2d called for non 2D template!");

    // Rotate vertex 90 degrees counterclockwise around (0.5, 0.5)
    auto rotate = [](Tiny::Vector<Real, 2> vertex)
    {
      return Tiny::Vector<Real, 2>{-vertex[1] + 1, vertex[0]};
    };

    auto rotate_indices = [](RawEntity<Shape_>& entity)
    {
      entity.coords = Intern::rotate_arraylike_2d(entity.coords);
    };

    return transform_template_vertices(rotate, rotate_indices, tmplt);
  }

  /**
   * \brief Rotates a raw template 90 degrees counterclockwise around the x-axis.
   *
   * Vertices are transformed such that the template is centered around the
   * origin, then the rotation is applied, then the vertices are transformed
   * back.
   */
  template<typename Shape_>
  static RawTemplate<Shape_> rotate_template_xaxis(RawTemplate<Shape_>& tmplt)
  {
    static_assert(Shape_::dimension == 3, "rotate_template_2d called for non 2D template!");

    auto rotate = [](Tiny::Vector<Real, 3> vertex)
    {
      return Tiny::Vector<Real, 3>{vertex[0], -vertex[2] + 1, vertex[1]};
    };

    auto rotate_indices = [](RawEntity<Shape_>& entity)
    {
      entity.coords = Intern::rotate_arraylike_xaxis(entity.coords);
    };

    return transform_template_vertices(rotate, rotate_indices, tmplt);
  }

  /**
   * \brief Rotates a raw template 90 degrees counterclockwise around the y-axis.
   *
   * Vertices are transformed such that the template is centered around the
   * origin, then the rotation is applied, then the vertices are transformed
   * back.
   */
  template<typename Shape_>
  static RawTemplate<Shape_> rotate_template_yaxis(RawTemplate<Shape_>& tmplt)
  {
    static_assert(Shape_::dimension == 3, "rotate_template_2d called for non 2D template!");

    auto rotate = [](Tiny::Vector<Real, 3> vertex)
    {
      return Tiny::Vector<Real, 3>{vertex[2], vertex[1], -vertex[0] + 1};
    };

    auto rotate_indices = [](RawEntity<Shape_>& entity)
    {
      entity.coords = Intern::rotate_arraylike_yaxis(entity.coords);
    };

    return transform_template_vertices(rotate, rotate_indices, tmplt);
  }

  /**
   * \brief Rotates a raw template 90 degrees counterclockwise around the z-axis.
   *
   * Vertices are transformed such that the template is centered around the
   * origin, then the rotation is applied, then the vertices are transformed
   * back.
   */
  template<typename Shape_>
  static RawTemplate<Shape_> rotate_template_zaxis(RawTemplate<Shape_>& tmplt)
  {
    static_assert(Shape_::dimension == 3, "rotate_template_2d called for non 2D template!");

    auto rotate = [](Tiny::Vector<Real, 3> vertex)
    {
      return Tiny::Vector<Real, 3>{-vertex[1] + 1, vertex[0], vertex[2]};
    };

    auto rotate_indices = [](RawEntity<Shape_>& entity)
    {
      entity.coords = Intern::rotate_arraylike_zaxis(entity.coords);
    };

    return transform_template_vertices(rotate, rotate_indices, tmplt);
  }


  template<typename TemplateMap, typename RefinementType_>
  int create_2dtemplate_rotations(TemplateMap& map, const RefinementType_& base_type)
  {
    // Count base template towards created templates
    int templates_created = 1;

    RefinementType_ type = base_type;
    auto tmplt = map[type];

    for(int i = 0; i < 4; i++)
    {
      if(map.find(type) == map.end())
      {
        map[type] = tmplt;
        templates_created += 1;
      }

      type = type.rotate_2d();
      tmplt = rotate_template_2d(tmplt);
    }

    return templates_created;
  }

  template<typename TemplateMap, typename RefinementType_>
  int create_3dtemplate_rotations(TemplateMap& map, const RefinementType_& base_type)
  {
    enum class Rotation
    {
      None,
      X,
      Y,
    };

    static constexpr std::array<Rotation, 6> rotations = {{
      Rotation::None,
      Rotation::X,
      Rotation::Y,
      Rotation::X,
      Rotation::Y,
      Rotation::X,
    }};

    // Count base template towards created templates
    int templates_created = 1;

    RefinementType_ outer_type = base_type;
    auto outer_template = map[outer_type];

    for(int outer_rotation(0); outer_rotation < 6; outer_rotation++)
    {
      const Rotation next_outer_rotation = rotations[outer_rotation];
      if (next_outer_rotation == Rotation::X)
      {
        outer_type = outer_type.rotate_xaxis();
        outer_template = rotate_template_xaxis(outer_template);
      }
      else if (next_outer_rotation == Rotation::Y)
      {
        outer_type = outer_type.rotate_yaxis();
        outer_template = rotate_template_yaxis(outer_template);
      }

      RefinementType_ inner_type = outer_type;
      auto inner_template = outer_template;
      for(int inner_rotation(0); inner_rotation < 4; inner_rotation++)
      {
        if(map.find(inner_type) == map.end())
        {
          map[inner_type] = inner_template;
          templates_created += 1;
        }

        inner_type = inner_type.rotate_zaxis();
        inner_template = rotate_template_zaxis(inner_template);
      }
    }
    return templates_created;
  }
} // namespace FEAT::Geometry
