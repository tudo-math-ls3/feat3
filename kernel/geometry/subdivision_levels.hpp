// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <functional>

namespace FEAT::Geometry
{
  /**
   * \brief Tuple-type for SubdivisionLevels
   *
   * Represents subdivision levels for a single mesh entity.
   *
   * \author Markus Muegge
   */
  template<int size_>
  struct SubdivisionLevelTuple
  {
    static_assert(size_ > 0, "invalid number of indices");

    /// Tuple size
    static constexpr int size = size_;

    /// Tuple values
    std::array<std::uint64_t, size_> levels;

    /// access operator
    Index& operator[](Index idx)
    {
      ASSERT(idx < size);
      return levels[idx];
    }

    /// access operator
    const Index& operator[](Index idx) const
    {
      ASSERT(idx < size);
      return levels[idx];
    }

    /**
     * \brief Decrement levels
     *
     * \returns A new tuple where each level has been decremented by one, to a minimum of zero
     */
    SubdivisionLevelTuple<size> decrement() const
    {
      SubdivisionLevelTuple<size> result;
      for(int i = 0; i < size; i++)
      {
        result[i] = levels[i] > 0 ? levels[i] - 1 : 0;
      }
      return result;
    }
  };

  /**
   * \brief Subdivision level markings for meshes
   *
   * Used to assign a subdivision level to each vertex of a mesh.
   *
   * \author Markus Muegge
   */
  class SubdivisionLevels
  {
    /// Debug output operator
    friend std::ostream& operator<<(std::ostream& /*unused*/, const SubdivisionLevels& /*unused*/);

  public:
    /**
     * \brief Constructor
     *
     * \param[in] num_vertices Number of vertices to mark
     */
    explicit SubdivisionLevels(Index num_vertices) : _levels(num_vertices)
    {
    }

    /**
     * \brief Constructor
     *
     * \param[in] num_vertices Number of vertices to mark
     * \param[in] default_level Default subdivision level for each vertex
     */
    SubdivisionLevels(Index num_vertices, std::uint64_t default_level) : _levels(num_vertices, default_level)
    {
    }

    /// move-constructor
    SubdivisionLevels(SubdivisionLevels&& other) noexcept :
      _levels(std::forward<std::vector<std::uint64_t>>(other._levels))
    {
    }

    /// move-assign operator
    SubdivisionLevels& operator=(SubdivisionLevels&& other) noexcept
    {
      // avoid self-move
      if(this == &other)
      {
        return *this;
      }

      _levels = std::forward<std::vector<std::uint64_t>>(other._levels);

      return *this;
    }

    /// deleted copy-constructor
    SubdivisionLevels(SubdivisionLevels& other) = delete;
    /// deleted copy-assign operator
    SubdivisionLevels& operator=(SubdivisionLevels& other) = delete;

    /// destructor
    ~SubdivisionLevels() = default;

    /// access operator
    std::uint64_t& operator[](Index idx)
    {
      ASSERT(idx < _levels.size());
      return _levels[idx];
    }

    /// access operator
    const std::uint64_t& operator[](Index idx) const
    {
      ASSERT(idx < _levels.size());
      return _levels[idx];
    }

    /**
     * \brief Returns the maximum refinement level assigned to any vertex
     */
    std::uint64_t maximum_level()
    {
      std::uint64_t result = 0;
      for(auto level : _levels)
      {
        result = std::max(result, level);
      }
      return result;
    }

    std::size_t size() const
    {
      return _levels.size();
    }

    /**
     * \brief Sets subdivision levels for each cell of the given mesh
     *
     * Vertices will be set to the maximum level of any cell they are contained in.
     *
     * \tparam MeshType
     * The type of the given mesh
     *
     * \param[in] mesh
     * The mesh for whose cells levels should be set
     * \param[in] f
     * A function assigning a level to each cell of \c mesh
     */
    template<typename MeshType>
    void set_per_element(const MeshType& mesh, std::function<std::uint64_t(Index)> f)
    {
      constexpr int cell_dim = MeshType::shape_dim;

      Index num_verts = mesh.get_num_entities(0);
      Index num_cells = mesh.get_num_entities(cell_dim);

      // Initialize levels to correct size
      _levels = std::vector<std::uint64_t>(num_verts, 0);

      for(Index cell(0); cell < num_cells; ++cell)
      {
        std::uint64_t cell_level = f(cell);
        set_level_for_element(mesh, cell, cell_level);
      }
    }

    /**
     * \brief Set subdivision level for a mesh element
     *
     * Assigns maximum of current level and element level to all vertices making up the chosen element.
     *
     * \param[in] mesh Mesh to take elements from
     * \param[in] element Element to set level for
     * \param[in] element_level Level to set
     */
    template<typename MeshType>
    void set_level_for_element(const MeshType& mesh, Index element, std::uint64_t element_level)
    {
      constexpr int cell_dim = MeshType::shape_dim;
      constexpr int verts_per_cell = Shape::FaceTraits<typename MeshType::ShapeType, 0>::count;

      auto& v_at_c = mesh.template get_index_set<cell_dim, 0>();
      for(Index vert(0); vert < verts_per_cell; ++vert)
      {
        Index vert_idx = v_at_c(element, vert);
        _levels[vert_idx] = std::max(_levels[vert_idx], element_level);
      }
    }

    /**
     * \brief Adds up two SubdivisionLevels
     */
    void merge(const SubdivisionLevels& other)
    {
      Index len = std::max(_levels.size(), other._levels.size());
      for(Index i(0); i < len; ++i)
      {
        _levels[i] += other._levels[i];
      }
    }

    /**
     * \brief Retrieve SubdivisionLevelTuple
     *
     * \param[in] indices Vertex-Indices to retrieve levels for
     */
    template<int n>
    SubdivisionLevelTuple<n> get_levels(const IndexTuple<n>& indices) const
    {
      SubdivisionLevelTuple<n> result = {};
      for(int i = 0; i < indices.num_indices; i++)
      {
        result[i] = _levels[indices[i]];
      }
      return result;
    }

  private:
    /// Level markings
    std::vector<std::uint64_t> _levels;
  };

  /// Debug output operator for SubdivisionLevels
  inline std::ostream& operator<<(std::ostream& os, const SubdivisionLevels& sdls)
  {
    os << "SDLs {" << "\n";
    for(auto sdl : sdls._levels)
    {
      os << "\t" << sdl << ",\n";
    }
    os << "}\n";
    return os;
  }
} // namespace FEAT::Geometry
