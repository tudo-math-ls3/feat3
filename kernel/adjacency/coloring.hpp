// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ADJACENCY_COLOURING_HPP
#define KERNEL_ADJACENCY_COLOURING_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

// includes, system
#include <vector>
#include <algorithm>

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
      Coloring(Coloring&& other);

      /// move-assign operator
      Coloring& operator=(Coloring&& other);

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
        std::vector<std::vector<int>> _coloring_maps;
        std::vector<void*> _d_coloring_maps;

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

      /// Returns the number of colors
      Index get_num_colors() const
      {
        return _coloring_maps.size();
      }

      /// Get the k-th color map
      std::vector<int>& get_color_map(Index k)
      {
        return _coloring_maps.at(k);
      }

      /// \copydoc get_color_map
      const std::vector<int>& get_color_map(Index k) const
      {
        return _coloring_maps.at(k);
      }

      /// Get the k-th device color map
      void* get_color_map_device(Index k)
      {
        return _d_coloring_maps.at(k);
      }

      /// \copydoc get_color_map_device
      const void* get_color_map_device(Index k) const
      {
        return _d_coloring_maps.at(k);
      }

      /// Retrieve the color maps
      std::vector<std::vector<int>>& get_coloring_maps()
      {
        return _coloring_maps;
      }

      /// \copydoc get_coloring_maps
      const std::vector<std::vector<int>>& get_coloring_maps() const
      {
        return _coloring_maps;
      }

      /// Retrieve the device  color maps
      std::vector<void*>& get_coloring_maps_device()
      {
        return _d_coloring_maps;
      }

      /// \copydoc get_coloring_maps_device
      const std::vector<void*>& get_coloring_maps_device() const
      {
        return _d_coloring_maps;
      }

      /// Initialize the device coloring data. Beforehand, fill_color should have been called.
      void init_device()
      {
        #ifdef FEAT_HAVE_CUDA
          _d_coloring_maps.resize(_coloring_maps.size());
          for(Index i = 0; i < _coloring_maps.size(); ++i)
          {
            _d_coloring_maps[i] = Util::cuda_malloc(_coloring_maps[i].size() * sizeof(int));
            Util::cuda_copy_host_to_device(_d_coloring_maps[i], _coloring_maps[i].data(), _coloring_maps[i].size() * sizeof(int));
          }
        #endif
      }

      /// Free the device data.
      void free_device()
      {
        #ifdef FEAT_HAVE_CUDA
          for(Index i = 0; i < _d_coloring_maps.size(); ++i)
          {
            Util::cuda_free(_d_coloring_maps[i]);
          }
          _d_coloring_maps.clear();
        #endif
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
          num_colors = *std::max_element(coloring.begin(), coloring.end());
        }
        _coloring_maps.resize(std::size_t(num_colors));
        for(std::size_t i = 0; i < coloring.size(); ++i)
        {
          _coloring_maps.at(std::size_t(coloring.at(i))).push_back(int(i));
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
        _coloring_maps.resize(std::size_t(num_colors));
        for(int i = 0; i < int(coloring.get_num_nodes()); ++i)
        {
          _coloring_maps.at(coloring[Index(i)]).push_back(i);
        }
      }

      bool initialized() const
      {
        return _d_coloring_maps.size() > 0u;
      }
    }; // class ColoringDataHandler
  } // namespace Adjacency
} // namespace FEAT

#endif // KERNEL_ADJACENCY_COLOURING_HPP
