#pragma once
#ifndef KERNEL_GEOMETRY_PARTI_ITERATIVE_HPP
#define KERNEL_GEOMETRY_PARTI_ITERATIVE_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/mutable_priority_queue.hpp>
#include <kernel/util/time_stamp.hpp>

#include <map>
#include <vector>
#include <list>

namespace FEAT
{
  namespace Geometry
  {
    namespace Intern
    {
      struct PartiIterativeItem
      {
        Index rank; //to which center does this cell belong
        Index distance; //distance to the nearest center

        PartiIterativeItem() :
          distance(std::numeric_limits<Index>::max())
        {
        }
      };

      template<typename Shape_, int num_coords_, int stride_, typename Coord_>
        std::vector<Index> parti_iterative_distance(Index start, ConformalMesh<Shape_, num_coords_, stride_, Coord_> & mesh, Index num_ranks)
        {
        typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
        static constexpr int facet_dim = MeshType::shape_dim-1;
        const auto& facet_idx = mesh.template get_index_set<MeshType::shape_dim, facet_dim>();
        // nachbar zu zelle i sind neighbours[i][j] mit j = 0 bis faced_idx.get_num_indices() und neighbours != ~Index(0)
        auto& neighbours = mesh.get_neighbours();

        Index num_elems(mesh.get_num_elements());

        // the list that shall contain our computed distance results
        std::vector<Index> distances(num_elems, std::numeric_limits<Index>::max());
        distances.at(start) = Index(0);
        // the 'queue' of nodes that have not yet been visited
        // the priorities are sorted from max to min, thus we need to use the inverse distance as the key, resulting in the lowest distance being at the front of the queue
        FEAT::mutable_priority_queue<Index, Index> pending_nodes;
        for (Index i(0) ; i < num_elems ; ++i)
        {
          pending_nodes.insert(i, 0);
        }
        pending_nodes.update(start, std::numeric_limits<Index>::max());

        Index exploration_threshold = Math::max(
            Index(Math::pow(double(num_elems), double(1) / double(MeshType::shape_dim)) + 1),
            num_elems / num_ranks);
        exploration_threshold = Math::max(exploration_threshold, Index(2));


        while (pending_nodes.size() > 0)
        {
          Index next_node = pending_nodes.front_value();
          Index next_distance = std::numeric_limits<Index>::max() - pending_nodes.front_key();
          // remove current node from queue of unvisited nodes
          pending_nodes.pop();
          //abort exploration if we have already reached a hopefully large enough portion of the mesh
          if (next_distance != std::numeric_limits<Index>::max() && next_distance > exploration_threshold)
          {
            return distances;
          }

          //update all neighbours
          for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
          {
            Index other_cell(neighbours[next_node][Index(j)]);
            //update neighbour cell distance, if neighbour exists and neighbour is in the queue of pending nodes
            if (other_cell != ~Index(0) && pending_nodes.count(other_cell) > 0)
            {
              Index new_distance = distances.at(next_node) + 1;
              if (distances.at(other_cell) > new_distance)
              {
                distances.at(other_cell) = new_distance;
                pending_nodes.update(other_cell, std::numeric_limits<Index>::max() - new_distance);
              }
            }
          }
        }
        return distances;
      }

      template<typename Shape_, int num_coords_, int stride_, typename Coord_>
      class PartiIterativeIndividual
      {
        public: //TODO write getter for cells per rank and make all member variables private
        /// our mesh type
        typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
        /// our mesh
        MeshType& _mesh;
        /// number of elements in input mesh
        const Index _num_elems;
        /// number of desired ranks/patches
        const Index _num_ranks;
        /// Mapping of ranks to cells assigned to each rank
        std::vector<std::set<Index>> _cells_per_rank;
        /// List of cluster centers: cell id for each rank
        std::set<Index> _centers;
        /// Our communicator
        const Dist::Comm & _comm;
        /// List of boundary cells per rank
        std::vector<std::set<Index>> _boundary_cells;
        /// mapping of cell to rank number
        std::vector<Index> _rank_per_cell;
        /// boundary size per rank, must be updated via update_boundary_size method
        std::vector<Index> _boundary_size_per_rank;

        public:

        /**
         * \brief Constructor
         *
         * \param[in] mesh
         * The mesh that is to be partitioned (on some refined level)
         *
         * \param[in] num_ranks
         * The desired number of ranks.
         */
        explicit PartiIterativeIndividual(MeshType& mesh, const Dist::Comm & comm, Random & rng) :
          _mesh(mesh),
          _num_elems(mesh.get_num_elements()),
          _num_ranks(Index(comm.size())),
          _cells_per_rank(_num_ranks),
          _comm(comm),
          _boundary_cells(_num_ranks),
          _rank_per_cell(_num_elems),
          _boundary_size_per_rank(_num_ranks)
        {
          XASSERT(_num_ranks > Index(0));
          XASSERT(_num_elems >= _num_ranks);

          // List of cells with informations regarding the partitioning layout
          std::vector<Intern::PartiIterativeItem> items(_num_elems);

          bool bad_centers(true);
          while (bad_centers)
          {
            bad_centers = false;

            // find randomly num_ranks different cells as cluster centers
            while(_centers.size() < _num_ranks)
            {
              Index cell = rng(Index(0), _num_elems - 1);
              if (_centers.count(cell) == 0)
                _centers.insert(cell);
            }

            // create distance list from each center and assign items the distance and rank of the nearest center
            auto center_iterator = _centers.begin();
            for (Index rank(0) ; rank < _num_ranks ; ++rank, ++center_iterator)
            {
              auto distance_list = Intern::parti_iterative_distance(*center_iterator, mesh, _num_ranks);
              for (Index node(0) ; node < _num_elems ; ++node)
              {
                if(distance_list.at(node) < items.at(node).distance)
                {
                  items.at(node).distance = distance_list.at(node);
                  items.at(node).rank = rank;
                }
              }
            }

            //check for max uint values in items list, i.e. we did not reach every cell
            for (Index node(0) ; node < _num_elems ; ++node)
            {
              bad_centers &= (items.at(node).distance == std::numeric_limits<Index>::max());
            }
          }

          for (Index i(0) ; i < _num_elems ; ++i)
          {
            _cells_per_rank.at(items.at(i).rank).insert(i);
          }

          //setup cell->rank mapping
          for (Index rank(0) ; rank < _num_ranks ; ++rank)
          {
            for (auto cell : _cells_per_rank.at(rank))
            {
              _rank_per_cell.at(cell) = rank;
            }
          }


          static constexpr int facet_dim = MeshType::shape_dim-1;
          const auto& facet_idx = mesh.template get_index_set<MeshType::shape_dim, facet_dim>();
          // nachbar zu zelle i sind neighbours[i][j] mit j = 0 bis faced_idx.get_num_indices() und neighbours != ~Index(0)
          auto& neighbours = mesh.get_neighbours();

          //setup boundary cells
          for (Index rank(0) ; rank < _num_ranks ; ++rank)
          {
            for (auto cell : _cells_per_rank.at(rank))
            {
              for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
              {
                Index other_cell(neighbours[cell][Index(j)]);
                // if neighbour exists and is not in our own patch
                if (other_cell != ~Index(0) && _rank_per_cell.at(other_cell) != rank)
                {
                  _boundary_cells.at(rank).insert(cell);
                }
              }
            }
          }

          update_boundary_size();

          //TODO anzahl nachbar ranks fuer fitness noetig ???
        }

        /// mutate individuum - let cells switch to smaller patches
        void mutate(MeshType& mesh, Random & rng, Index mutation_max)
        {
          static constexpr int facet_dim = MeshType::shape_dim-1;
          const auto& facet_idx = mesh.template get_index_set<MeshType::shape_dim, facet_dim>();
          // nachbar zu zelle i sind neighbours[i][j] mit j = 0 bis faced_idx.get_num_indices() und neighbours != ~Index(0)
          auto& neighbours = mesh.get_neighbours();

          Index mutation_count(0);
          Index tries(0);
          //iterate over every cell at the boundaries
          while (tries < 10*mutation_max && mutation_count < mutation_max)
          {
            ++tries;

            Index rank(rng(Index(0), _num_ranks - 1));
            if (_boundary_cells.at(rank).size() == 0)
              continue;
            if (_cells_per_rank.at(rank).size() == 1)
              continue;

            Index rng_cell_idx = rng(Index(0), Index(_boundary_cells.at(rank).size() - 1));
            Index counter(0);
            Index cell(*(_boundary_cells.at(rank).begin()));
            //choose cell randomly from current boundary set
            /// \todo make it more efficient
            for (auto it(_boundary_cells.at(rank).begin()) ; counter != rng_cell_idx ; ++counter, ++it)
            {
              cell = *it;
            }

            // list of neighbour cells in other patches
            std::list<Index> trans_rank_neighbours;
            for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
            {
              Index other_cell(neighbours[cell][Index(j)]);
              // if neighbour exists and is not in our own patch
              if (other_cell != ~Index(0) && _rank_per_cell.at(other_cell) != rank)
              {
                trans_rank_neighbours.push_back(other_cell);
              }
            }

            //check if any neighbour has a smaller patch
            Index smallest_rank(rank);
            Index smallest_size = Index(_cells_per_rank.at(rank).size());
            for (auto neighbour : trans_rank_neighbours)
            {
              if (_cells_per_rank.at(_rank_per_cell.at(neighbour)).size() < smallest_size)
              {
                smallest_size = Index(_cells_per_rank.at(_rank_per_cell.at(neighbour)).size());
                smallest_rank = _rank_per_cell.at(neighbour);
              }
            }

            if (smallest_rank == rank)
              continue;

            //build up list of neighbours in the same (old) rank
            std::list<Index> in_rank_neighbours;
            for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
            {
              Index other_cell(neighbours[cell][Index(j)]);
              // if neighbour exists and is in our own patch
              if (other_cell != ~Index(0) && _rank_per_cell.at(other_cell) == rank)
              {
                in_rank_neighbours.push_back(other_cell);
              }
            }
            //check if we would isolate one of our neighbours completly and abort
            //i.e. our neighbour from the same rank has only us as a link to its patch
            bool neighbour_missing = false;
            for (auto neighbour : in_rank_neighbours)
            {
              bool other_neighbour = false;
              for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
              {
                Index other_cell(neighbours[neighbour][Index(j)]);
                if (other_cell != ~Index(0) && other_cell != cell && _rank_per_cell.at(other_cell) == rank)
                {
                  other_neighbour = true;
                  break;
                }
              }
              if (other_neighbour == false)
                neighbour_missing = true;
            }
            if (neighbour_missing)
              continue;

            //found new rank for our cell, update structures
            _rank_per_cell.at(cell) = smallest_rank;
            XASSERT(_cells_per_rank.at(rank).count(cell) > 0);
            _cells_per_rank.at(rank).erase(cell);
            _cells_per_rank.at(smallest_rank).insert(cell);
            XASSERT(_boundary_cells.at(rank).count(cell) > 0);
            _boundary_cells.at(rank).erase(cell);
            _boundary_cells.at(smallest_rank).insert(cell);

            //update neighbour cells in other rank if they loose their boundary property due to our patch switch
            for (auto neighbour : trans_rank_neighbours)
            {
              //update only cells in common new rank
              if (_rank_per_cell.at(neighbour) != smallest_rank)
                continue;

              bool boundary(false);
              for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
              {
                Index other_cell(neighbours[neighbour][Index(j)]);
                if (other_cell != ~Index(0) && _rank_per_cell.at(other_cell) != smallest_rank)
                {
                  boundary = true;
                  break;
                }
              }

              if (! boundary)
                _boundary_cells.at(smallest_rank).erase(neighbour);
            }

            //update our own ranks neighbours to become a boundary cell
            for (auto neighbour : in_rank_neighbours)
            {
              _boundary_cells.at(rank).insert(neighbour);
            }

            //we mutate successfully
            ++mutation_count;
          }
          update_boundary_size();
        }

        /// update boundary size per rank from _boundary_cells structure
        void update_boundary_size()
        {
          static constexpr int facet_dim = MeshType::shape_dim-1;
          const auto& facet_idx = _mesh.template get_index_set<MeshType::shape_dim, facet_dim>();
          auto& neighbours = _mesh.get_neighbours();

          for (Index rank(0) ; rank < _num_ranks ; ++rank)
          {
            Index sum(0);
            for (auto cell : _boundary_cells.at(rank))
            {
              for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
              {
                Index other_cell(neighbours[cell][Index(j)]);
                if (other_cell != ~Index(0) && _rank_per_cell.at(other_cell) != rank)
                {
                  ++sum;
                }
              }
            }
            _boundary_size_per_rank.at(rank) = sum;
          }
        }

        Index get_boundary_size() const
        {
          Index sum(0);
          for (Index rank(0) ; rank < _num_ranks ; ++rank)
          {
            sum += get_boundary_size(rank);
          }
          return sum;
        }

        Index get_boundary_size(Index rank) const
        {
          return _boundary_size_per_rank.at(rank);
        }

        Index get_boundary_deviation() const
        {
          Index min(std::numeric_limits<Index>::max());
          Index max(0);
          for (Index rank(0) ; rank < _num_ranks ; ++rank)
          {
            Index ranks_boundary_size(get_boundary_size(rank));
            min = Math::min(ranks_boundary_size, min);
            max = Math::max(ranks_boundary_size, max);
          }
          return max - min;
        }

        Index get_cell_deviation() const
        {
          Index min(std::numeric_limits<Index>::max());
          Index max(0);
          for (Index rank(0) ; rank < _num_ranks ; ++rank)
          {
            min = Math::min(Index(_cells_per_rank.at(rank).size()), min);
            max = Math::max(Index(_cells_per_rank.at(rank).size()), max);
          }
          return max - min;
        }
      };

      template<typename T_>
      bool parti_iterative_compare_cell_deviation_simple(const T_ & first, const T_ & second)
      {
        if (first.get_cell_deviation() != second.get_cell_deviation())
          return first.get_cell_deviation() < second.get_cell_deviation();
        else if (first.get_boundary_deviation() != second.get_boundary_deviation() )
          return first.get_boundary_deviation() < second.get_boundary_deviation();
        else
          return first.get_boundary_size() < second.get_boundary_size();
      }

      template<typename T_>
      bool parti_iterative_compare_cell_deviation(const T_ & first, const T_ & second)
      {
        if (first.get_cell_deviation() != second.get_cell_deviation())
          return first.get_cell_deviation() < second.get_cell_deviation() && first.get_boundary_size() <= second.get_boundary_size() + 1;
        else if (first.get_boundary_deviation() != second.get_boundary_deviation() )
          return first.get_boundary_deviation() < second.get_boundary_deviation() && first.get_boundary_size() <= second.get_boundary_size() + 1;
        else
          return first.get_boundary_size() < second.get_boundary_size();
      }
    }
    /** \brief Iterative-Partitioner class template declaration */
    template<typename Mesh_>
    class PartiIterative;

    /**
     * \brief Iterative-Partitioner class template specialisation for ConformalMesh
     *
     * The basic usage of this class is as follows:
     * -# Refine the mesh until it contains at least num_rank cells.
     * -# Create an object of this class and pass the to-be-partitioned mesh as well
     *    as the desired number of ranks/patches to the constructor.
     * -# Create the Elements-At-Rank graph using the #build_elems_at_rank() function.
     *
     * \author Dirk Ribbrock
     */
    template<typename Shape_, int num_coords_, int stride_, typename Coord_>
    class PartiIterative<ConformalMesh<Shape_, num_coords_, stride_, Coord_>>
    {
      private:
        /// number of elements in input mesh
        const Index _num_elems;
        /// number of desired ranks/patches
        const Index _num_ranks;
        /// Our communicator
        const Dist::Comm & _comm;
        /// Our population
        std::list<Geometry::Intern::PartiIterativeIndividual<Shape_, num_coords_, stride_, Coord_>> _population;

      public:
      /// our mesh type
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;

      /**
       * \brief Constructor
       *
       * \param[in] mesh
       * The mesh that is to be partitioned (on some refined level).
       *
       * \param[in] comm
       * The global communicator to be used.
       *
       * \param[in] time_init
       * The amount of seconds to be used for initial patch center search.
       *
       * \param[in] time_mutate
       * The amount of seconds to be used for partitioning optimisation via patch mutation.
       */
      explicit PartiIterative(MeshType& mesh, const Dist::Comm & comm, double time_init, double time_mutate) :
        _num_elems(mesh.get_num_elements()),
        _num_ranks(Index(comm.size())),
        _comm(comm)
      {
        Random::SeedType seed(Random::SeedType(time(nullptr)));
        Random rng(seed + Random::SeedType(comm.rank()));

        mesh.fill_neighbours();

        {
          Geometry::Intern::PartiIterativeIndividual<Shape_, num_coords_, stride_, Coord_> indi(mesh, _comm, rng);
          _population.push_back(indi);
        }

        TimeStamp at;
        while(at.elapsed_now() < time_init)
        {
          for (Index i(0) ; i < 10 ; ++i)
          {
            Geometry::Intern::PartiIterativeIndividual<Shape_, num_coords_, stride_, Coord_> indi(mesh, _comm, rng);
            _population.push_back(indi);
          }
          _population.sort(Geometry::Intern::parti_iterative_compare_cell_deviation_simple<typename decltype(_population)::value_type>);
          while (_population.size() > 10)
            _population.pop_back();
        }

        //std::cout<<"pre: "<<_population.back().get_cell_deviation()<<" " <<_population.back().get_boundary_deviation()<<" " << _population.back().get_boundary_size()<<std::endl;

        //mutate each of the 10 individuals and store the top 10 fittest of all in _population
        at.stamp();
        while(at.elapsed_now() < time_mutate)
        {
          std::list<typename decltype(_population)::value_type> mutations;
          for (auto it(_population.begin()) ; it != _population.end() ; ++it)
          {
            auto indi(*it);
            indi.mutate(mesh, rng, 5);
            mutations.push_back(indi);
          }
          mutations.sort(Geometry::Intern::parti_iterative_compare_cell_deviation<typename decltype(_population)::value_type>);
          _population.merge(mutations, Geometry::Intern::parti_iterative_compare_cell_deviation<typename decltype(_population)::value_type>);
          while (_population.size() > 20)
            _population.pop_back();
        }

        _population.sort(Geometry::Intern::parti_iterative_compare_cell_deviation<typename decltype(_population)::value_type>);
        while (_population.size() > 1)
          _population.pop_back();

        //std::cout<<_population.back().get_cell_deviation()<<" " <<_population.back().get_boundary_deviation()<<" " << _population.back().get_boundary_size()<<std::endl;
      }

      /**
       * \brief Returns the Elements-at-Rank graph of the partitioning.
       *
       * \returns
       * The Elements-at-Rank graph of the partitioning.
       */
      Adjacency::Graph build_elems_at_rank() const
      {
        auto& indi = _population.front();

        Index own_cell_deviation = indi.get_cell_deviation();
        Index lowest_cell_deviation(std::numeric_limits<Index>::max());
        Index own_boundary_deviation = indi.get_boundary_deviation();
        Index lowest_boundary_deviation(std::numeric_limits<Index>::max());
        Index own_boundary_size = indi.get_boundary_size();
        Index lowest_boundary_size(std::numeric_limits<Index>::max());
        Index lowest_rank(0);
        ///\todo use mpi_reduce_all with op_minloc
        Index * all_cell_deviation = new Index[_num_ranks];
        Index * all_boundary_deviation = new Index[_num_ranks];
        Index * all_boundary_size = new Index[_num_ranks];
        _comm.allgather(&own_cell_deviation, 1, all_cell_deviation, 1);
        _comm.allgather(&own_boundary_deviation, 1, all_boundary_deviation, 1);
        _comm.allgather(&own_boundary_size, 1, all_boundary_size, 1);
        for (Index i(0) ; i < _num_ranks ; ++i)
        {
          if (all_cell_deviation[i] < lowest_cell_deviation)
          {
            lowest_cell_deviation = all_cell_deviation[i];
            lowest_boundary_deviation = all_boundary_deviation[i];
            lowest_boundary_size = all_boundary_size[i];
            lowest_rank = i;
          }
          else if (all_cell_deviation[i] == lowest_cell_deviation && all_boundary_deviation[i] < lowest_boundary_deviation)
          {
            lowest_cell_deviation = all_cell_deviation[i];
            lowest_boundary_deviation = all_boundary_deviation[i];
            lowest_boundary_size = all_boundary_size[i];
            lowest_rank = i;
          }
          else if (all_cell_deviation[i] == lowest_cell_deviation && all_boundary_deviation[i] == lowest_boundary_deviation && all_boundary_size[i] < lowest_boundary_size)
          {
            lowest_cell_deviation = all_cell_deviation[i];
            lowest_boundary_deviation = all_boundary_deviation[i];
            lowest_boundary_size = all_boundary_size[i];
            lowest_rank = i;
          }
        }
        delete[] all_cell_deviation;
        delete[] all_boundary_deviation;
        delete[] all_boundary_size;

        int _rank = _comm.rank();
        if (_rank == (int)lowest_rank)
        {
          //std::cout<<"selected: " << _population.front().get_cell_deviation()<<" " <<_population.front().get_boundary_deviation()<<" " << _population.front().get_boundary_size()<<std::endl;
          Adjacency::Graph graph(_num_ranks, _num_elems, _num_elems);
          Index* ptr = graph.get_domain_ptr();
          Index* idx = graph.get_image_idx();

          // build pointer array
          ptr[0] = 0;
          for(Index i(0); i < _num_ranks; ++i)
          {
            ptr[i+1] = Index(indi._cells_per_rank.at(i).size()) + ptr[i];
          }

          // build index array
          Index counter(0);
          for (Index rank(0) ; rank < _num_ranks ; ++rank)
          {
            for (auto cell : indi._cells_per_rank.at(rank))
            {
              idx[counter] = cell;
              ++counter;
            }
          }
          _comm.bcast(graph.get_domain_ptr(), _num_ranks + 1, _rank);
          _comm.bcast(graph.get_image_idx(), _num_elems, _rank);
          return graph;
        }
        else
        {
          Index * domain_ptr = new Index[_num_ranks + 1];
          _comm.bcast(domain_ptr, _num_ranks + 1, (int)lowest_rank);
          Index * image_idx = new Index[_num_elems];
          _comm.bcast(image_idx, _num_elems, (int)lowest_rank);
          Adjacency::Graph graph(_num_ranks, _num_elems, _num_elems, domain_ptr, image_idx);
          delete[] domain_ptr;
          delete[] image_idx;
          return graph;
        }
      }
    }; // class PartiIterative
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PARTI_ITERATIVE_HPP
