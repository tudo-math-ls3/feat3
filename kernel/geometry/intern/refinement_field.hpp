// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/geometry/index_set.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace FEAT::Geometry::Intern
{
  template<typename VertexMarking_, int size_>
  class RefinementFieldTuple
  {
    static_assert(size_ > 0, "Invalid RefinementFieldTuple size. Must be greater than 0.");
  public:
    /// Type of vertex marking
    using VertexMarkingType = VertexMarking_;

    /// Tuple size
    static constexpr std::size_t size = size_;

    /// Tuple values
    std::array<VertexMarkingType, size_> levels;

    /// access operator
    VertexMarkingType& operator[](Index idx)
    {
      ASSERT(idx < size);
      return levels[idx];
    }

    /// access operator
    const VertexMarkingType& operator[](Index idx) const
    {
      ASSERT(idx < size);
      return levels[idx];
    }
  };

  template<typename VertexMarking_>
  class RefinementField
  {
  public:
    using VertexMarkingType = VertexMarking_;

  private:
    std::vector<VertexMarkingType> _markings;

  public:
    explicit RefinementField(Index num_vertices) : _markings(num_vertices)
    {
    }

    RefinementField(Index num_vertices, VertexMarkingType default_marking) : _markings(num_vertices, default_marking)
    {
    }

    /// move-constructor
    RefinementField(RefinementField&& other) noexcept :
      _markings(std::forward<std::vector<VertexMarkingType>>(other._markings))
    {
    }

    /// move-assign operator
    RefinementField& operator=(RefinementField&& other) noexcept
    {
      if(this == &other)
      {
        return *this;
      }

      _markings = std::forward<std::vector<VertexMarkingType>>(other._markings);

      return *this;
    }

    /// deleted copy-constructor
    RefinementField(RefinementField& other) = delete;

    /// deleted copy-assign operator
    RefinementField& operator=(RefinementField& other) = delete;

    /// destructor
    ~RefinementField() = default;

    /// access operator
    VertexMarkingType& operator[](Index idx)
    {
      ASSERT(idx < _markings.size());
      return _markings[idx];
    }

    /// const access operator
    const VertexMarkingType& operator[](Index idx) const
    {
      ASSERT(idx < _markings.size());
      return _markings[idx];
    }

    std::size_t size() const
    {
      return _markings.size();
    }

    template<int n_>
    RefinementFieldTuple<VertexMarkingType, n_> get_tuple(const Geometry::IndexTuple<n_>& indices) const
    {
      RefinementFieldTuple<VertexMarkingType, n_> result = {};

      for(int i(0); i < n_; ++i)
      {
        result[i] = _markings[indices[i]];
      }
      return result;
    }
  };
} // namespace FEAT::Geometry::Intern
