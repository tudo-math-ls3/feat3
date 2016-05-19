// includes, FEAST
#include <kernel/adjacency/colouring.hpp>
#include <kernel/adjacency/graph.hpp>

// includes, system
#include <vector>
#include <utility>

namespace FEAST
{
  namespace Adjacency
  {
    // Default constructor.
    Colouring::Colouring() :
      _num_nodes(0),
      _num_colours(0),
      _colouring(nullptr)
    {
    }

    // Allocation Constructor.
    Colouring::Colouring(
      Index num_nodes,
      Index num_colours)
       :
      _num_nodes(num_nodes),
      _num_colours(num_colours),
      _colouring(nullptr)
    {
      _colouring = new Index[_num_nodes];
    }

    // Array Constructor
    Colouring::Colouring(
      Index num_nodes,
      Index* colouring)
    {
      _num_nodes = num_nodes;
      _colouring = new Index[_num_nodes];
      _num_colours = 0;
      for(Index i(0); i< _num_nodes; i++)
      {
        _colouring[i] = colouring[i];
        if(colouring[i] >= _num_colours)
        {
          ++_num_colours;
        }
      }
    }

    // Creation out of a given Graph
    Colouring::Colouring(const Graph& graph)
    {
      // get the graph's data
      Index num_nodes_domain = graph.get_num_nodes_domain();
      const Index* domain_ptr = graph.get_domain_ptr();
      const Index* image_idx = graph.get_image_idx();

      _num_nodes = num_nodes_domain;

      // allocate the colouring array
      _colouring = new Index[_num_nodes];

      // highest possible number of colours
      Index mnc = graph.degree() + 1;

      // number of used colours
      _num_colours = 0;

      // auxiliary vector storing temporary information about used colours
      std::vector<Index> col_aux(mnc);

      // auxiliary vector storing the number of colour uses
      std::vector<Index> col_num(mnc, Index(0));

      // auxiliary variables
      Index lower_bound;
      Index upper_bound;
      Index min_colour_uses;
      Index min_colour_index;

      // loop over all nodes
      for(Index i(0); i < num_nodes_domain; ++i)
      {
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
          if(image_idx[j] < i)
          {
            // mark the used colour
            col_aux[_colouring[image_idx[j]]] = 1;
          }
        }

        // calculate new colour:

        // minimum number of colour uses so far
        min_colour_uses = num_nodes_domain + 1;

        // colour with min_colour_uses
        min_colour_index = num_nodes_domain + 1;

        // loop over all colours
        for(Index j(0); j < _num_colours; ++j)
        {
          if(col_aux[j] !=1)
          {
            if(col_num[j] < min_colour_uses)
            {
              min_colour_uses = col_num[j];
              min_colour_index = j;
            }
          }
        }

        // check if an old colour was found or if a new colour is needed
        if(min_colour_index != num_nodes_domain +1)
        {
          // if an old colour is acceptable
          _colouring[i] = min_colour_index;
          ++col_num[min_colour_index];
        }
        else
        {
          // if a new colour is needed
          _colouring[i] = _num_colours;
          ++col_num[_num_colours];
          ++_num_colours;
        }
      }
    }

    // Creation out of a given Graph with a prescribed order
    Colouring::Colouring(const Graph& graph, const Index* order)
    {
      // get the graph's data
      Index num_nodes_domain = graph.get_num_nodes_domain();
      const Index* domain_ptr = graph.get_domain_ptr();
      const Index* image_idx = graph.get_image_idx();

      _num_nodes = num_nodes_domain;

      // allocate the colouring array
      _colouring = new Index[_num_nodes];

      // highest possible number of colours
      Index mnc = graph.degree() + 1;

      // number of used colours
      _num_colours = 0;

      // auxiliary vector storing temporary information about used colours
      std::vector<Index> col_aux(mnc);

      // auxiliary vector storing the number of colour uses
      std::vector<Index> col_num(mnc, Index(0));

      // auxiliary variables
      Index lower_bound;
      Index upper_bound;
      Index min_colour_uses;
      Index min_colour_index;

      // auxiliary index
      Index i;

      // initialise colouring data
      for(Index j(0); j < num_nodes_domain; ++j)
      {
        _colouring[j] = num_nodes_domain + 1;
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
          if(_colouring[image_idx[j]] != num_nodes_domain + 1)
          {
            // mark the used colour
            col_aux[_colouring[image_idx[j]]] = 1;
          }
        }

        // calculate new colour:

        // minimum number of colour uses so far
        min_colour_uses = num_nodes_domain + 1;

        // colour with min_colour_uses
        min_colour_index = num_nodes_domain + 1;

        // loop over all colours
        for(Index j(0); j < _num_colours; ++j)
        {
          if(col_aux[j] !=1)
          {
            if(col_num[j] < min_colour_uses)
            {
              min_colour_uses = col_num[j];
              min_colour_index = j;
            }
          }
        }

        // check if an old colour was found or if a new colour is needed
        if(min_colour_index != num_nodes_domain +1)
        {
          // if an old colour is acceptable
          _colouring[i] = min_colour_index;
          ++col_num[min_colour_index];
        }
        else
        {
          // if a new colour is needed
          _colouring[i] = _num_colours;
          ++col_num[_num_colours];
          ++_num_colours;
        }
      }
    }

    // move ctor
    Colouring::Colouring(Colouring&& other) :
      _num_nodes(other._num_nodes),
      _num_colours(other._num_colours),
      _colouring(other._colouring)
    {
      other._num_nodes = other._num_colours = Index(0);
      other._colouring = nullptr;
    }

    // move-assign operator
    Colouring& Colouring::operator=(Colouring&& other)
    {
      // avoid self-move
      if(this == &other)
        return *this;

      if(_colouring != nullptr)
        delete [] _colouring;

      _num_nodes = other._num_nodes;
      _num_colours = other._num_colours;
      _colouring = other._colouring;

      other._num_nodes = other._num_colours = Index(0);
      other._colouring = nullptr;

      return *this;
    }

    // virtual destructor
    Colouring::~Colouring()
    {
      if(_colouring != nullptr)
        delete [] _colouring;
    }

    Graph Colouring::create_partition_graph() const
    {
      // allocate a new graph
      Graph graph(get_max_colour() + 1, get_num_nodes(), get_num_nodes());

      // create domain array
      Index* domain_ptr = graph.get_domain_ptr();
      domain_ptr[0] = Index(0);

      // create image array
      Index* image_idx = graph.get_image_idx();

      // index counter
      Index idx_counter = Index(0);

      // loop over all colours
      for(Index i(0); i < _num_colours; ++i)
      {
        // loop over all nodes
        for(Index j(0); j < _num_nodes; ++j)
        {
          // if node j has the colour i
          if(i == _colouring[j])
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
} // namespace FEAST
