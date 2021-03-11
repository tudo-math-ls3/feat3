// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/adjacency/graph.hpp>

// includes, system
#include <vector>
#include <utility>
#include <unordered_set>

namespace FEAT
{
  namespace Adjacency
  {
    // Default constructor.
    Coloring::Coloring() :
      _num_colors(0),
      _coloring()
    {
    }

    // Allocation Constructor.
    Coloring::Coloring(
      Index num_nodes,
      Index num_colors)
       :
      _num_colors(num_colors)
    {
      _coloring = IndexVector(num_nodes);
    }

    // Array Constructor
    Coloring::Coloring(
      Index num_nodes,
      Index* coloring)
    {

      std::unordered_set<Index> colors;
      for (Index i(0) ; i < num_nodes ; ++i)
      {
        colors.insert(coloring[i]);
      }
      _num_colors = Index(colors.size());
      colors.clear();

      _coloring =IndexVector(num_nodes);
      for(Index i(0); i< num_nodes; i++)
      {
        _coloring[i] = coloring[i];
      }

    }
    //Vector Constructor

    Coloring::Coloring(
      Index num_colors,
      const IndexVector& coloring) :
      _num_colors(num_colors),
      _coloring(coloring)
    {
    }
    // Creation out of a given Graph
    Coloring::Coloring(const Graph& graph) :
      Coloring()
    {
      // get the graph's data
      Index num_nodes_domain = graph.get_num_nodes_domain();
      const Index* domain_ptr = graph.get_domain_ptr();
      const Index* image_idx = graph.get_image_idx();

      if(num_nodes_domain == Index(0))
        return;

      // allocate the coloring vector
      _coloring = IndexVector(num_nodes_domain);

      // highest possible number of colors
      Index mnc = graph.degree() + 1;

      // number of used colors
      _num_colors = 0;

      // auxiliary vector storing temporary information about used colors
      std::vector<Index> col_aux(mnc);

      // auxiliary vector storing the number of color uses
      std::vector<Index> col_num(mnc, Index(0));

      // auxiliary variables
      Index lower_bound;
      Index upper_bound;
      Index min_color_uses;
      Index min_color_index;

      // loop over all nodes
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        // reset the auxiliary vector
        for(Index j(0); j < mnc; ++j)
        {
          col_aux[j] = 0;
        }

        // loop over all adjancies of the i-th node
        lower_bound = domain_ptr[i];
        upper_bound = domain_ptr[i+1];
        for(Index j(lower_bound); j < upper_bound; ++j)
        {
          if(image_idx[j] < i)
          {
            // mark the used color
            col_aux[_coloring[image_idx[j]]] = 1;
          }
        }

        // calculate new color:

        // minimum number of color uses so far
        min_color_uses = num_nodes_domain + 1;

        // color with min_color_uses
        min_color_index = num_nodes_domain + 1;

        // loop over all colors
        for(Index j(0); j < _num_colors; ++j)
        {
          if(col_aux[j] !=1)
          {
            if(col_num[j] < min_color_uses)
            {
              min_color_uses = col_num[j];
              min_color_index = j;
            }
          }
        }

        // check if an old color was found or if a new color is needed
        if(min_color_index != num_nodes_domain +1)
        {
          // if an old color is acceptable
          _coloring[i] = min_color_index;
          ++col_num[min_color_index];
        }
        else
        {
          // if a new color is needed
          _coloring[i] = _num_colors;
          ++col_num[_num_colors];
          ++_num_colors;
        }
      }
    }

    // Creation out of a given Graph with a prescribed order
    Coloring::Coloring(const Graph& graph, const Index* order) :
      Coloring()
    {
      // get the graph's data
      Index num_nodes_domain = graph.get_num_nodes_domain();
      const Index* domain_ptr = graph.get_domain_ptr();
      const Index* image_idx = graph.get_image_idx();

      if(num_nodes_domain == Index(0))
        return;

      // allocate the coloring vector
      _coloring = IndexVector(num_nodes_domain);

      // highest possible number of colors
      Index mnc = graph.degree() + 1;

      // number of used colors
      _num_colors = 0;

      // auxiliary vector storing temporary information about used colors
      std::vector<Index> col_aux(mnc);

      // auxiliary vector storing the number of color uses
      std::vector<Index> col_num(mnc, Index(0));

      // auxiliary variables
      Index lower_bound;
      Index upper_bound;
      Index min_color_uses;
      Index min_color_index;

      // auxiliary index
      Index i;

      // initialize coloring data
      for(Index j(0); j < num_nodes_domain; ++j)
      {
        _coloring[j] = num_nodes_domain + 1;
      }

      // loop over all nodes
      for(Index k(0); k < num_nodes_domain; ++k)
      {
        // get index of the current node
        i = order[k];

        // reset the auxiliary array
        for(Index j(0); j < mnc; ++j)
        {
          col_aux[j] = 0;
        }

        // loop over all adjancies of the i-th node
        lower_bound = domain_ptr[i];
        upper_bound = domain_ptr[i+1];
        for(Index j(lower_bound); j < upper_bound; ++j)
        {
          if(_coloring[image_idx[j]] != num_nodes_domain + 1)
          {
            // mark the used color
            col_aux[_coloring[image_idx[j]]] = 1;
          }
        }

        // calculate new color:

        // minimum number of color uses so far
        min_color_uses = num_nodes_domain + 1;

        // color with min_color_uses
        min_color_index = num_nodes_domain + 1;

        // loop over all colors
        for(Index j(0); j < _num_colors; ++j)
        {
          if(col_aux[j] !=1)
          {
            if(col_num[j] < min_color_uses)
            {
              min_color_uses = col_num[j];
              min_color_index = j;
            }
          }
        }

        // check if an old color was found or if a new color is needed
        if(min_color_index != num_nodes_domain +1)
        {
          // if an old color is acceptable
          _coloring[i] = min_color_index;
          ++col_num[min_color_index];
        }
        else
        {
          // if a new color is needed
          _coloring[i] = _num_colors;
          ++col_num[_num_colors];
          ++_num_colors;
        }
      }
    }

    // move ctor
    Coloring::Coloring(Coloring&& other) :
      _num_colors(other._num_colors),
      _coloring(std::forward<IndexVector>(other._coloring))
    {
      other._num_colors = Index(0);
      other._coloring.clear();
    }

    // move-assign operator
    Coloring& Coloring::operator=(Coloring&& other)
    {
      // avoid self-move
      if(this == &other)
        return *this;

      _num_colors = other._num_colors;
      _coloring = std::forward<IndexVector>(other._coloring);

      other._num_colors = Index(0);
      other._coloring.clear();

      return *this;
    }

    // virtual destructor
    Coloring::~Coloring()
    {}

    Graph Coloring::create_partition_graph() const
    {
      // allocate a new graph
      Graph graph(get_max_color() + 1, get_num_nodes(), get_num_nodes());

      // create domain vector
      Index* domain_ptr = graph.get_domain_ptr();
      domain_ptr[0] = Index(0);

      // create image vector
      Index* image_idx = graph.get_image_idx();

      // index counter
      Index idx_counter = Index(0);

      // loop over all colors
      for(Index i(0); i < _num_colors; ++i)
      {
        // loop over all nodes
        for(Index j(0); j < _coloring.size(); ++j)
        {
          // if node j has the color i
          if(i == _coloring[j])
          {
            image_idx[idx_counter] = j;
            ++idx_counter;
          }
        }
        domain_ptr[i+1] = idx_counter;
      }

      // return graph
      return graph;
    }

  } // namespace Adjacency
} // namespace FEAT
