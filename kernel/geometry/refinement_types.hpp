// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/geometry/intern/adaptive_refinement_utils.hpp>
#include <kernel/geometry/intern/congruency_mapping.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/geometry/intern/refinement_field.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/shape.hpp>

#include <bitset>

namespace FEAT::Geometry
{
  namespace Intern
  {
    template<typename Shape_>
    static std::bitset<Shape::FaceTraits<Shape_, 0>::count>
    orient_bitset(const std::bitset<Shape::FaceTraits<Shape_, 0>::count>& bits, int orientation)
    {
      using CongruencyMapping = CongruencyMapping<Shape_, 0>;
      static constexpr const int num_bits = Shape::FaceTraits<Shape_, 0>::count;

      std::bitset<num_bits> result = 0;

      for(int current_bit = 0; current_bit < num_bits; current_bit++)
      {
        if(bits[current_bit])
        {
          int mapped_bit = CongruencyMapping::map(orientation, current_bit);
          result[mapped_bit] = true;
        }
      }
      return result;
    }
  } // namespace Intern

  template<typename Shape_>
  struct VertexMarking
  {
    static constexpr int num_vertices = Shape::FaceTraits<Shape_, 0>::count;

    std::bitset<num_vertices> bits;

    template<typename S_>
    friend bool operator==(const VertexMarking<S_>&, const VertexMarking<S_>&);

    template<typename S_>
    friend std::ostream& operator<<(std::ostream&, const VertexMarking<S_>&);

    explicit VertexMarking(std::uint64_t markings) : bits(markings)
    {
    }

    explicit VertexMarking(std::bitset<num_vertices> markings) : bits(markings)
    {
    }

    explicit VertexMarking(const SubdivisionLevelTuple<num_vertices>& tuple) : bits(0)
    {
      for(int i = 0; i < num_vertices; i++)
      {
        bits[i] = (tuple[i] > 0);
      }
    }

    static VertexMarking all_marked()
    {
      return VertexMarking((1ULL << (std::size_t)num_vertices) - 1ULL);
    }

    bool covers(const VertexMarking& other) const
    {
      // We cover other, if no bits in other are set, that are not already set for us
      return (bits | other.bits) == bits;
    }

    bool is_vertex_marked(int vertex) const
    {
      ASSERT(vertex < num_vertices);
      return bits[vertex];
    }

    std::size_t distance(const VertexMarking& other) const
    {
      return (bits ^ other.bits).count();
    }
  };

  template<typename Shape_>
  bool operator==(const VertexMarking<Shape_>& lhs, const VertexMarking<Shape_>& rhs)
  {
    return lhs._markings == rhs._markings;
  }

  template<typename Shape_>
  bool operator!=(const VertexMarking<Shape_>& lhs, const VertexMarking<Shape_>& rhs)
  {
    return !(lhs == rhs);
  }

  template<typename Shape_>
  std::ostream& operator<<(std::ostream& stream, const VertexMarking<Shape_>& type)
  {
    stream << "VertexMarkings<" << Shape_::name() << "> { markings: " << type.bits << " }";
    return stream;
  }

  template<typename Shape_>
  class StandardRefinementType
  {
    static constexpr int num_vertices = Shape::FaceTraits<Shape_, 0>::count;

    std::bitset<num_vertices> _type;

    template<typename S_>
    friend bool operator==(const StandardRefinementType<S_>&, const StandardRefinementType<S_>&);

    template<typename S_>
    friend std::ostream& operator<<(std::ostream&, const StandardRefinementType<S_>&);

    template<typename T_>
    friend struct std::hash;

  public:
    explicit StandardRefinementType(std::uint64_t type) : _type(type)
    {
    }

    explicit StandardRefinementType(std::bitset<num_vertices> bits) : _type(bits)
    {
    }

    explicit StandardRefinementType(VertexMarking<Shape_> markings) : _type(markings.bits)
    {
    }

    explicit StandardRefinementType(
      const Intern::template RefinementFieldTuple<std::uint64_t, num_vertices>& markings) :
      _type(0)
    {
      for(int i = 0; i < num_vertices; i++)
      {
        _type[i] = (markings[i] > 0);
      }
    }

    bool is_full_refinement() const
    {
      return _type.all();
    }

    bool is_zero_refinement() const
    {
      return _type.none();
    }

    template<int dim_>
    StandardRefinementType<typename Shape::FaceTraits<Shape_, dim_>::ShapeType> face_type(int face) const
    {
      using FaceShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      using Mapping = Intern::FaceIndexMapping<Shape_, dim_, 0>;

      static constexpr const int num_verts = Shape::FaceTraits<FaceShape, 0>::count;

      std::bitset<num_verts> result;
      for(int i(0); i < num_verts; i++)
      {
        auto vertex = static_cast<Index>(Mapping::map(face, i));
        result[i] = _type[vertex];
      }

      return StandardRefinementType<FaceShape>(result);
    }

    StandardRefinementType rotate_2d() const
    {
      return StandardRefinementType(Intern::rotate_arraylike_2d(_type));
    }

    StandardRefinementType rotate_xaxis() const
    {
      return StandardRefinementType(Intern::rotate_arraylike_xaxis(_type));
    }

    StandardRefinementType rotate_yaxis() const
    {
      return StandardRefinementType(Intern::rotate_arraylike_yaxis(_type));
    }

    StandardRefinementType rotate_zaxis() const
    {
      return StandardRefinementType(Intern::rotate_arraylike_zaxis(_type));
    }

    StandardRefinementType orient(int orientation) const
    {
      return StandardRefinementType(Intern::orient_bitset<Shape_>(_type, orientation));
    }

    VertexMarking<Shape_> to_vertex_marking() const
    {
      return VertexMarking<Shape_>(_type);
    }

    Index to_number() const
    {
      return _type.to_ulong();
    }
  };

  template<typename Shape_>
  bool operator==(const StandardRefinementType<Shape_>& lhs, const StandardRefinementType<Shape_>& rhs)
  {
    return lhs._type == rhs._type;
  }

  template<typename Shape_>
  bool operator!=(const StandardRefinementType<Shape_>& lhs, const StandardRefinementType<Shape_>& rhs)
  {
    return !(lhs == rhs);
  }

  template<typename Shape_>
  std::ostream& operator<<(std::ostream& stream, const StandardRefinementType<Shape_>& type)
  {
    stream << "StandardRefinementType<" << Shape_::name() << "> { type: " << type._type << " }";
    return stream;
  }

  struct IsolatedPointVertexMarking
  {
    std::uint64_t level = 0;
    bool is_isolated = false;

    IsolatedPointVertexMarking() = default;

    IsolatedPointVertexMarking(std::uint64_t l, bool i) : level(l), is_isolated(i)
    {
    }
  };

  inline bool operator==(const IsolatedPointVertexMarking& a, const IsolatedPointVertexMarking& b)
  {
    return a.level == b.level && a.is_isolated == b.is_isolated;
  }

  inline bool operator!=(const IsolatedPointVertexMarking& a, const IsolatedPointVertexMarking& b)
  {
    return !(a == b);
  }

  inline std::ostream& operator<<(std::ostream& stream, const IsolatedPointVertexMarking& marking)
  {
    stream << "IsolatedPointVertexMarking{ level: " << unsigned(marking.level) << ", is_isolated: " << marking.is_isolated << " }";
    return stream;
  }

  template<typename Shape_>
  class IsolatedPointRefinementType
  {
    static constexpr int num_vertices = Shape::FaceTraits<Shape_, 0>::count;

    std::bitset<num_vertices> _type;
    bool _isolated;

    template<typename S_>
    friend bool operator==(const IsolatedPointRefinementType<S_>&, const IsolatedPointRefinementType<S_>&);

    template<typename S_>
    friend std::ostream& operator<<(std::ostream&, const IsolatedPointRefinementType<S_>&);

    template<typename T_>
    friend struct std::hash;

  public:
    explicit IsolatedPointRefinementType(std::uint64_t type, bool isolated) : _type(type), _isolated(isolated)
    {
    }

    explicit IsolatedPointRefinementType(std::bitset<num_vertices> bits, bool isolated) :
      _type(bits),
      _isolated(isolated)
    {
    }

    explicit IsolatedPointRefinementType(VertexMarking<Shape_> markings, bool isolated) :
      _type(markings.bits),
      _isolated(isolated)
    {
    }

    explicit IsolatedPointRefinementType(
      const Intern::template RefinementFieldTuple<IsolatedPointVertexMarking, num_vertices>& markings) :
      _type(0),
      _isolated(false)
    {
      for(Index i = 0; i < num_vertices; i++)
      {
        _type[i] = (markings[i].level > 0);
        _isolated = _isolated || markings[i].is_isolated;
      }
    }

    bool is_full_refinement() const
    {
      return _type.all();
    }

    bool is_zero_refinement() const
    {
      return _type.none();
    }

    template<int dim_>
    IsolatedPointRefinementType<typename Shape::FaceTraits<Shape_, dim_>::ShapeType> face_type(int face) const
    {
      using FaceShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      using Mapping = Intern::FaceIndexMapping<Shape_, dim_, 0>;

      static constexpr const int num_verts = Shape::FaceTraits<FaceShape, 0>::count;

      std::bitset<num_verts> result_bits;

      for(int i(0); i < num_verts; i++)
      {
        auto vertex = static_cast<Index>(Mapping::map(face, i));
        result_bits[i] = _type[vertex];
      }

      // Result type has same isolation as us, as long as we picked any marked vertices
      return IsolatedPointRefinementType<FaceShape>(result_bits, _isolated && result_bits.any());
    }

    IsolatedPointRefinementType rotate_2d() const
    {
      return IsolatedPointRefinementType(Intern::rotate_arraylike_2d(_type), _isolated);
    }

    IsolatedPointRefinementType rotate_xaxis() const
    {
      return IsolatedPointRefinementType(Intern::rotate_arraylike_xaxis(_type), _isolated);
    }

    IsolatedPointRefinementType rotate_yaxis() const
    {
      return IsolatedPointRefinementType(Intern::rotate_arraylike_yaxis(_type), _isolated);
    }

    IsolatedPointRefinementType rotate_zaxis() const
    {
      return IsolatedPointRefinementType(Intern::rotate_arraylike_zaxis(_type), _isolated);
    }

    IsolatedPointRefinementType orient(int orientation) const
    {
      return IsolatedPointRefinementType(Intern::orient_bitset<Shape_>(_type, orientation), _isolated);
    }

    VertexMarking<Shape_> to_vertex_marking() const
    {
      return VertexMarking<Shape_>(_type);
    }

    Index to_number() const
    {
      return _type.to_ulong();
    }
  };

  template<typename Shape_>
  bool operator==(const IsolatedPointRefinementType<Shape_>& lhs, const IsolatedPointRefinementType<Shape_>& rhs)
  {
    return lhs._type == rhs._type && lhs._isolated == rhs._isolated;
  }

  template<typename Shape_>
  bool operator!=(const IsolatedPointRefinementType<Shape_>& lhs, const IsolatedPointRefinementType<Shape_>& rhs)
  {
    return !(lhs == rhs);
  }

  template<typename Shape_>
  std::ostream& operator<<(std::ostream& stream, const IsolatedPointRefinementType<Shape_>& type)
  {
    stream << "IsolatedPointRefinementType<" << Shape_::name() << "> { type: " << type._type
           << ", isolated: " << stringify(type._isolated) << " }";
    return stream;
  }

} // namespace FEAT::Geometry

namespace std
{
  template<typename Shape_>
  struct hash<FEAT::Geometry::StandardRefinementType<Shape_>>
  {
    static constexpr int num_vertices = FEAT::Shape::template FaceTraits<Shape_, 0>::count;

    std::size_t operator()(const FEAT::Geometry::StandardRefinementType<Shape_>& type) const
    {
      return std::hash<std::bitset<num_vertices>>{}(type._type);
    }
  };

  template<typename Shape_>
  struct hash<FEAT::Geometry::IsolatedPointRefinementType<Shape_>>
  {
    static constexpr int num_vertices = FEAT::Shape::template FaceTraits<Shape_, 0>::count;

    std::size_t operator()(const FEAT::Geometry::IsolatedPointRefinementType<Shape_>& type) const
    {
      auto h1 = std::hash<std::bitset<num_vertices>>{}(type._type);
      auto h2 = std::hash<bool>{}(type._isolated);
      return h1 ^ (h2 << 1UL);
    }
  };
} // namespace std
