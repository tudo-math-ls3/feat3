// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/memory_pool.hpp>

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

// includes, system
#include <vector>
#include <algorithm>
#include <numeric>

namespace FEAT
{
  namespace Adjacency
  {
    // forward declaration
    class Graph;

    /**
     * \brief Coloring object implementation
     *
     * \todo detailed description
     *
     * \author Constantin Christof
     */
    class Coloring
    {
    protected:
      using IndexVector = std::vector < Index>;

      /// total number of colors used
      Index _num_colors;

      /**
       * \brief coloring vector
       *
       * Dimension: #_num_nodes
       * The coloring vector is defined as the second row of the following chart:
       *    node number :         0                   1           ...         num_nodes-1
       *    color number: color of node 0    color of node 1    ...  color of node num_nodes-1
       *
       */
      IndexVector _coloring;

    public:

      /**
       * \brief Default constructor.
       *
       * This constructor creates a new empty coloring object, but does not allocate any vectors.
       */
      Coloring();

      /**
       * \brief Allocation Constructor.
       *
       * This constructor creates a new coloring object and allocates
       * a coloring vector of length num_nodes.
       *
       * \note This constructor does not initialize the allocated vector -- it has to be initialized by the user
       * after construction.
       *
       * \param[in] num_nodes
       * The total number of nodes.
       *
       * \param[in] num_colors
       * The number of colors.
       *
       */
      Coloring(
        Index num_nodes,
        Index num_colors);

      /**
       * \brief Array Constructor.
       *
       * This constructor creates a new coloring object from a given array.
       *
       * \param[in] num_nodes
       * The total number of nodes.
       *
       * \param[in] coloring
       * The \transient coloring array
       */
      Coloring(
        Index num_nodes,
        Index* coloring);

      /**
      * \brief Vector Constructor.
      *
      * This constructor creates a new coloring object from a given vector.
      *
      * \param[in] num_colors
      * The total number of colors.
      *
      * \param[in] coloring
      * The \transient coloring vector
      */
      Coloring(
        Index num_colors,
        const IndexVector& coloring);

      /**
       * \brief Creation out of a given Graph
       *
       * This constructor creates a coloring object out of a graph. The returned coloring
       * satisfies the condition that adjacent nodes do not have the same color.
       *
       * \param[in] graph
       * The \transient graph to create the coloring from.
       */
      explicit Coloring(const Graph& graph);

      /**
       * \brief Creation out of a given Graph with a prescribed order
       *
       * This constructor creates a coloring object out of a graph. The returned coloring
       * satisfies the condition that adjacent nodes do not have the same color. The algorithm
       * proceeds through the nodes as given by the parameter order.
       *
       * \param[in] graph
       * The \transient graph to create the coloring from.
       *
       * \param[in] order
       * Permutation array that describes how the algorithm is supposed to proceed through the nodes.
       */
      explicit Coloring(const Graph& graph, const Index* order);

      /// move ctor
      Coloring(Coloring&& other) noexcept;

      /// move-assign operator
      Coloring& operator=(Coloring&& other) noexcept;

      /// virtual destructor
      virtual ~Coloring();

      /**
       * \brief Clones this coloring.
       *
       * \returns A deep-copy of this coloring.
       */
      Coloring clone() const
      {
        if (!_coloring.empty())
          //return Coloring(Index(_coloring.size()), _coloring.data());
          return Coloring(_num_colors, _coloring);
        else
          return Coloring();
      }

      /**
       * \brief Creates a color partition graph.
       *
       * This function creates a graph out of this coloring object.
       * The graph's image is the set of nodes, the domain is the set of colors.
       */
      Graph create_partition_graph() const;

      /**
       * \brief Returns the coloring array.
       * \returns The coloring array.
       */
      Index* get_coloring()
      {
        return _coloring.data();
      }

      /** \copydoc get_coloring() */
      const Index* get_coloring() const
      {
        return _coloring.data();
      }

      /**
       * \brief Returns the total number of nodes.
       *
       * \returns The total number of nodes.
       */
      Index get_num_nodes() const
      {
        return Index(_coloring.size());
      }

      /// \copydoc get_num_nodes()
      auto size() const -> decltype(_coloring.size())
      {
        return _coloring.size();
      }

      /**
       * \brief Returns the maximum color index.
       *
       * \returns The maximum color index.
       */
      Index get_max_color() const
      {
        return _num_colors - 1;
      }

      /**
       * \brief Returns the number of colors.
       *
       * \returns The number of colors.
       */
      Index get_num_colors() const
      {
        return _num_colors;
      }

      /// Checks whether the coloring is empty
      bool empty() const
      {
        return _coloring.empty();
      }

      /// Returns the color for a node i.
      Index& operator[](Index i)
      {
        ASSERT(i < _coloring.size());
        return _coloring[i];
      }

      /// Returns the color for a node i.
      const Index& operator[](Index i) const
      {
        ASSERT(i < _coloring.size());
        return _coloring[i];
      }
    }; // class Coloring

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
        std::vector<int*> _coloring_maps;
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
        _coloring_maps = std::move(other._coloring_maps);
        _coloring_map_sizes = std::move(other._coloring_map_sizes);
        other._coloring_maps.clear();
        other._coloring_map_sizes.clear();
      }

      ColoringDataHandler& operator=(ColoringDataHandler&& other) noexcept
      {
        if(this == &other)
          return *this;
        release_color();
        _coloring_maps = std::move(other._coloring_maps);
        _coloring_map_sizes = std::move(other._coloring_map_sizes);
        other._coloring_maps.clear();
        other._coloring_map_sizes.clear();
        return *this;
      }

      /// Returns the number of colors
      Index get_num_colors() const
      {
        return _coloring_map_sizes.size();
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

      /// Get the k-th color map
      int* get_color_map(Index k)
      {
        return _coloring_maps.at(k);
      }

      /// \copydoc get_color_map
      const int* get_color_map(Index k) const
      {
        return _coloring_maps.at(k);
      }

      /// Retrieve the color maps
      std::vector<int*>& get_coloring_maps()
      {
        return _coloring_maps;
      }

      /// \copydoc get_coloring_maps
      const std::vector<int*>& get_coloring_maps() const
      {
        return _coloring_maps;
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
        _coloring_maps.resize(std::size_t(num_colors));
        _coloring_map_sizes.resize(std::size_t(num_colors));
        for(Index i = 0; i < Index(num_colors); ++i)
        {
          _coloring_maps.at(i) = MemoryPool::allocate_memory<int>(tmp_vector.at(i).size());
          MemoryPool::copy(_coloring_maps.at(i), tmp_vector.at(i).data(), tmp_vector.at(i).size());
          _coloring_map_sizes.at(i) = tmp_vector.at(i).size();
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
        _coloring_maps.resize(std::size_t(num_colors));
        _coloring_map_sizes.resize(std::size_t(num_colors));
        for(Index i = 0; i < Index(num_colors); ++i)
        {
          _coloring_maps.at(i) = MemoryPool::allocate_memory<int>(tmp_vector.at(i).size());
          MemoryPool::copy(_coloring_maps.at(i), tmp_vector.at(i).data(), tmp_vector.at(i).size());
          _coloring_map_sizes.at(i) = Index(tmp_vector.at(i).size());
        }
      }

      void release_color()
      {
        for(Index i = 0; i < _coloring_maps.size(); ++i)
        {
          MemoryPool::release_memory(_coloring_maps.at(i));
        }
        _coloring_maps.clear();
        _coloring_map_sizes.clear();
      }

      bool initialized() const
      {
        return _coloring_maps.size() > 0u;
      }

      Index get_max_size() const
      {
        return std::accumulate(_coloring_map_sizes.begin(), _coloring_map_sizes.end(), Index(0), [](const Index& a, const Index& b){return std::max(a, b);});
      }
    }; // class ColoringDataHandler
  } // namespace Adjacency
} // namespace FEAT
