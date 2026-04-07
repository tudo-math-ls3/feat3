// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/util/memory_arbiter.hpp>

// includes, system
#include <vector>
#include <algorithm>
#include <numeric>

namespace FEAT
{
  namespace Adjacency
  {
    /**
     * \brief Datahandler for inverse coloring data
     *
     * This class is initialized by coloring data and also holds (if supported) coloring data on the gpu.
     *
     * \warning This class does not handle allocating and freeing device data by itself, i.e. you have to
     *          call the init_device and free_device function yourself.
     *
     */
    class ColoringDataHandler
    {
    public:
      /// vector of unified memory pointer
      Memory::Arbiter _coloring_map;
      /// vector of coloring sizes
      std::vector<Index> _coloring_map_offsets;
      /// vector of coloring sizes
      std::vector<Index> _coloring_map_sizes;

      explicit ColoringDataHandler() = default;

      /**
       * \brief Constructor receiving coloring data and optional hint
       *
       * \tparam ColoringType_ The type of the input coloring data.
       *
       * Constructs the coloring data handle from a coloring data array,
       * i.e. a vector of ints or a Coloring graph.
       *
       * \param[in] coloring Some array of integers, representing colors.
       * \param[in] hint Optionally, provide the numbers of colors.
       *                 If negative, this is deduced automatically.
       * \warning If hint is smaller then the maximum color number, this will
       *          lead to an allocation error.
       */
      template<typename ColoringType_>
      ColoringDataHandler(const ColoringType_& coloring, int hint = -1)
      {
        fill_color(coloring, hint);
      }

      ~ColoringDataHandler()
      {
        release_color();
      }

      ColoringDataHandler(const ColoringDataHandler&) = delete;

      ColoringDataHandler& operator=(const ColoringDataHandler&) = delete;


      ColoringDataHandler(ColoringDataHandler&& other) noexcept
      {
        _coloring_map = std::move(other._coloring_map);
        _coloring_map_offsets = std::move(other._coloring_map_offsets);
        _coloring_map_sizes = std::move(other._coloring_map_sizes);
        other._coloring_map.release();
        other._coloring_map_offsets.clear();
        other._coloring_map_sizes.clear();
      }

      ColoringDataHandler& operator=(ColoringDataHandler&& other) noexcept
      {
        if(this == &other)
          return *this;
        release_color();
        _coloring_map = std::move(other._coloring_map);
        _coloring_map_offsets = std::move(other._coloring_map_offsets);
        _coloring_map_sizes = std::move(other._coloring_map_sizes);
        other._coloring_map.release();
        other._coloring_map_offsets.clear();
        other._coloring_map_sizes.clear();
        return *this;
      }

      /// Returns the number of colors
      Index get_num_colors() const
      {
        return _coloring_map_sizes.size();
      }

      Index get_color_offset(Index k) const
      {
        return _coloring_map_offsets.at(k);
      }

      std::vector<Index>& get_color_offsets()
      {
        return _coloring_map_offsets;
      }

      const std::vector<Index>& get_color_offsets() const
      {
        return _coloring_map_offsets;
      }

      Index get_color_size(Index k) const
      {
        return _coloring_map_sizes.at(k);
      }

      std::vector<Index>& get_color_sizes()
      {
        return _coloring_map_sizes;
      }

      const std::vector<Index>& get_color_sizes() const
      {
        return _coloring_map_sizes;
      }

      Memory::TypedView<int> color_map_view(Memory::Location mem_loc, Memory::Access mem_acc)
      {
        return _coloring_map.typed_view<int>(mem_loc, mem_acc);
      }

      Memory::TypedView<int> color_map_view(Memory::Location mem_loc, Memory::Access mem_acc = Memory::Access::read) const
      {
        return _coloring_map.typed_view<int>(mem_loc, mem_acc);
      }

      /// Get max size of all colors
      Index get_max_color_size() const
      {
        return Index(std::accumulate(_coloring_map_sizes.begin(), _coloring_map_sizes.end(), Index(0), [](Index a, Index b){return std::max(a,b);}));
      }

      /**
       * \brief Fill in the coloring array
       *
       * \param[in] coloring Array mapping a dof to a color.
       * \param[in] hint Hint on the number of colors. If negative, this is deduced automatically.
       */
      void fill_color(const std::vector<int>& coloring, int hint = -1)
      {
        int num_colors = hint;
        if(hint < 0)
        {
          num_colors = *std::max_element(coloring.begin(), coloring.end()) + 1;
        }
        // fill tmp vector with coloring
        std::vector<std::vector<int>> tmp_vector;
        tmp_vector.resize(Index(num_colors));
        for(std::size_t i = 0; i < coloring.size(); ++i)
        {
          tmp_vector.at(std::size_t(coloring.at(i))).push_back(int(i));
        }

        _coloring_map = Memory::Arbiter(coloring.size() * sizeof(int));
        _coloring_map_offsets.resize(std::size_t(num_colors));
        _coloring_map_sizes.resize(std::size_t(num_colors));
        Memory::TypedView<int> coloring_map_view(_coloring_map.view(Memory::Location::main, Memory::Access::write));
        Index offset(0);
        for(Index i = 0; i < Index(num_colors); ++i)
        {
          Memory::memcopy_main(&coloring_map_view.get_w()[offset], tmp_vector.at(i).data(), tmp_vector.at(i).size() * sizeof(int));
          _coloring_map_offsets.at(i) = offset;
          _coloring_map_sizes.at(i) = tmp_vector.at(i).size();
          offset += _coloring_map_sizes.at(i);
        }
      }

      /**
       * \brief Fill in the coloring array
       *
       * \param[in] coloring The coloring graph.
       * \param[in] hint Hint on the number of colors. If negative, this is deduced automatically.
       *                 For this version, the value is not used.
       */
      void fill_color(const Coloring& coloring, int hint = -1)
      {
        int num_colors = int(coloring.get_num_colors());
        if(hint >= 0)
        {
          ASSERTM(num_colors == hint, "Hint and number of colors do not fit!");
        }

        // fill tmp vector with coloring
        std::vector<std::vector<int>> tmp_vector;
        tmp_vector.resize(Index(num_colors));
        for(std::size_t i = 0; i < coloring.size(); ++i)
        {
          tmp_vector.at(std::size_t(coloring[i])).push_back(int(i));
        }
        _coloring_map = Memory::Arbiter(coloring.size() * sizeof(int));
        _coloring_map_offsets.resize(std::size_t(num_colors));
        _coloring_map_sizes.resize(std::size_t(num_colors));
        Memory::TypedView<int> coloring_map_view(_coloring_map.view(Memory::Location::main, Memory::Access::write));
        Index offset(0);
        for(Index i = 0; i < Index(num_colors); ++i)
        {
          Memory::memcopy_main(&coloring_map_view.get_w()[offset], tmp_vector.at(i).data(), tmp_vector.at(i).size() * sizeof(int));
          _coloring_map_offsets.at(i) = offset;
          _coloring_map_sizes.at(i) = tmp_vector.at(i).size();
          offset += _coloring_map_sizes.at(i);
        }
      }

      void release_color()
      {
        _coloring_map.release();
        _coloring_map_sizes.clear();
        _coloring_map_offsets.clear();
      }

      bool initialized() const
      {
        return !_coloring_map_sizes.empty();
      }

      Index get_max_size() const
      {
        return std::accumulate(_coloring_map_sizes.begin(), _coloring_map_sizes.end(), Index(0), [](const Index& a, const Index& b){return std::max(a, b);});
      }
    }; // class ColoringDataHandler
  } // namespace Adjacency
} // namespace FEAT
