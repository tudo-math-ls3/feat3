// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
        Index patch; //to which center does this cell belong
        Index distance; //distance to the nearest center

        PartiIterativeItem() :
          distance(std::numeric_limits<Index>::max())
        {
        }
      };

      template<typename Shape_, int num_coords_, typename Coord_>
        std::vector<Index> parti_iterative_distance(Index start, ConformalMesh<Shape_, num_coords_, Coord_> & mesh, Index num_patches)
        {
        typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
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
            num_elems / num_patches);
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
            Index other_cell(neighbours[next_node][j]);
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

      template<typename Shape_, int num_coords_, typename Coord_>
      class PartiIterativeIndividual
      {
        public: //TODO write getter for cells per rank and make all member variables private
        /// our mesh type
        typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
        /// our mesh
        MeshType& _mesh;
        /// number of desired patches
        const Index _num_patches;
        /// number of elements in input mesh
        const Index _num_elems;
        /// Mapping of patches to cells assigned to each patch
        std::vector<std::set<Index>> _cells_per_patch;
        /// List of cluster centers: cell id for each patch
        std::set<Index> _centers;
        /// List of boundary cells per patch
        std::vector<std::set<Index>> _boundary_cells;
        /// mapping of cell to patch number
        std::vector<Index> _patch_per_cell;
        /// boundary size per patch, must be updated via update_boundary_size method
        std::vector<Index> _boundary_size_per_patch;

        public:

        /**
         * \brief Constructor
         *
         * \param[in] mesh
         * The mesh that is to be partitioned (on some refined level)
         *
         * \param[in] num_patches
         * The number of desired patches
         */
        explicit PartiIterativeIndividual(MeshType& mesh, Random & rng, Index num_patches) :
          _mesh(mesh),
          _num_patches(num_patches),
          _num_elems(mesh.get_num_elements()),
          _cells_per_patch(_num_patches),
          _boundary_cells(_num_patches),
          _patch_per_cell(_num_elems),
          _boundary_size_per_patch(_num_patches)
        {
          XASSERT(_num_patches > Index(0));
          XASSERT(_num_elems >= _num_patches);

          // List of cells with informations regarding the partitioning layout
          std::vector<Intern::PartiIterativeItem> items;

          bool bad_centers(true);
          while (bad_centers)
          {
            bad_centers = false;
            _centers.clear();
            items.clear();
            items.resize(_num_elems);

            // find randomly num_patches different cells as cluster centers
            while(_centers.size() < _num_patches)
            {
              Index cell = rng(Index(0), _num_elems - 1);
              if (_centers.count(cell) == 0)
                _centers.insert(cell);
            }

            // create distance list from each center and assign items the distance and patch of the nearest center
            auto center_iterator = _centers.begin();
            for (Index patch(0) ; patch < _num_patches ; ++patch, ++center_iterator)
            {
              auto distance_list = Intern::parti_iterative_distance(*center_iterator, mesh, _num_patches);
              for (Index node(0) ; node < _num_elems ; ++node)
              {
                if(distance_list.at(node) < items.at(node).distance)
                {
                  items.at(node).distance = distance_list.at(node);
                  items.at(node).patch = patch;
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
            _cells_per_patch.at(items.at(i).patch).insert(i);
          }

          //setup cell->patch mapping
          for (Index patch(0) ; patch < _num_patches ; ++patch)
          {
            for (auto cell : _cells_per_patch.at(patch))
            {
              _patch_per_cell.at(cell) = patch;
            }
          }


          static constexpr int facet_dim = MeshType::shape_dim-1;
          const auto& facet_idx = mesh.template get_index_set<MeshType::shape_dim, facet_dim>();
          // nachbar zu zelle i sind neighbours[i][j] mit j = 0 bis faced_idx.get_num_indices() und neighbours != ~Index(0)
          auto& neighbours = mesh.get_neighbours();

          //setup boundary cells
          for (Index patch(0) ; patch < _num_patches ; ++patch)
          {
            for (auto cell : _cells_per_patch.at(patch))
            {
              for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
              {
                Index other_cell(neighbours[cell][j]);
                // if neighbour exists and is not in our own patch
                if (other_cell != ~Index(0) && _patch_per_cell.at(other_cell) != patch)
                {
                  _boundary_cells.at(patch).insert(cell);
                }
              }
            }
          }

          update_boundary_size();

          //TODO anzahl nachbar patches fuer fitness noetig ???
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

            Index patch(rng(Index(0), _num_patches - 1));
            if (_boundary_cells.at(patch).size() == 0)
              continue;
            if (_cells_per_patch.at(patch).size() == 1)
              continue;

            Index rng_cell_idx = rng(Index(0), Index(_boundary_cells.at(patch).size() - 1));
            Index counter(0);
            Index cell(*(_boundary_cells.at(patch).begin()));
            //choose cell randomly from current boundary set
            /// \todo make it more efficient
            for (auto it(_boundary_cells.at(patch).begin()) ; counter != rng_cell_idx ; ++counter, ++it)
            {
              cell = *it;
            }

            // list of neighbour cells in other patches
            std::list<Index> trans_patch_neighbours;
            for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
            {
              Index other_cell(neighbours[cell][j]);
              // if neighbour exists and is not in our own patch
              if (other_cell != ~Index(0) && _patch_per_cell.at(other_cell) != patch)
              {
                trans_patch_neighbours.push_back(other_cell);
              }
            }

            //check if any neighbour has a smaller patch
            Index smallest_patch(patch);
            Index smallest_size = Index(_cells_per_patch.at(patch).size());
            for (auto neighbour : trans_patch_neighbours)
            {
              if (_cells_per_patch.at(_patch_per_cell.at(neighbour)).size() < smallest_size)
              {
                smallest_size = Index(_cells_per_patch.at(_patch_per_cell.at(neighbour)).size());
                smallest_patch = _patch_per_cell.at(neighbour);
              }
            }

            if (smallest_patch == patch)
              continue;

            //build up list of neighbours in the same (old) patch
            std::list<Index> in_patch_neighbours;
            for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
            {
              Index other_cell(neighbours[cell][j]);
              // if neighbour exists and is in our own patch
              if (other_cell != ~Index(0) && _patch_per_cell.at(other_cell) == patch)
              {
                in_patch_neighbours.push_back(other_cell);
              }
            }
            //check if we would isolate one of our neighbours completly and abort
            //i.e. our neighbour from the same patch has only us as a link to its patch
            bool neighbour_missing = false;
            for (auto neighbour : in_patch_neighbours)
            {
              bool other_neighbour = false;
              for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
              {
                Index other_cell(neighbours[neighbour][j]);
                if (other_cell != ~Index(0) && other_cell != cell && _patch_per_cell.at(other_cell) == patch)
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

            //found new patch for our cell, update structures
            _patch_per_cell.at(cell) = smallest_patch;
            XASSERT(_cells_per_patch.at(patch).count(cell) > 0);
            _cells_per_patch.at(patch).erase(cell);
            _cells_per_patch.at(smallest_patch).insert(cell);
            XASSERT(_boundary_cells.at(patch).count(cell) > 0);
            _boundary_cells.at(patch).erase(cell);
            _boundary_cells.at(smallest_patch).insert(cell);

            //update neighbour cells in other patch if they loose their boundary property due to our patch switch
            for (auto neighbour : trans_patch_neighbours)
            {
              //update only cells in common new patch
              if (_patch_per_cell.at(neighbour) != smallest_patch)
                continue;

              bool boundary(false);
              for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
              {
                Index other_cell(neighbours[neighbour][j]);
                if (other_cell != ~Index(0) && _patch_per_cell.at(other_cell) != smallest_patch)
                {
                  boundary = true;
                  break;
                }
              }

              if (! boundary)
                _boundary_cells.at(smallest_patch).erase(neighbour);
            }

            //update our own patchs neighbours to become a boundary cell
            for (auto neighbour : in_patch_neighbours)
            {
              _boundary_cells.at(patch).insert(neighbour);
            }

            //we mutate successfully
            ++mutation_count;
          }
          update_boundary_size();
        }

        /// update boundary size per patch from _boundary_cells structure
        void update_boundary_size()
        {
          static constexpr int facet_dim = MeshType::shape_dim-1;
          const auto& facet_idx = _mesh.template get_index_set<MeshType::shape_dim, facet_dim>();
          auto& neighbours = _mesh.get_neighbours();

          for (Index patch(0) ; patch < _num_patches ; ++patch)
          {
            Index sum(0);
            for (auto cell : _boundary_cells.at(patch))
            {
              for (int j(0) ; j < facet_idx.get_num_indices() ; ++j)
              {
                Index other_cell(neighbours[cell][j]);
                if (other_cell != ~Index(0) && _patch_per_cell.at(other_cell) != patch)
                {
                  ++sum;
                }
              }
            }
            _boundary_size_per_patch.at(patch) = sum;
          }
        }

        Index get_boundary_size() const
        {
          Index sum(0);
          for (Index patch(0) ; patch < _num_patches ; ++patch)
          {
            sum += get_boundary_size(patch);
          }
          return sum;
        }

        Index get_boundary_size(Index patch) const
        {
          return _boundary_size_per_patch.at(patch);
        }

        Index get_boundary_deviation() const
        {
          Index min(std::numeric_limits<Index>::max());
          Index max(0);
          for (Index patch(0) ; patch < _num_patches ; ++patch)
          {
            Index patches_boundary_size(get_boundary_size(patch));
            min = Math::min(patches_boundary_size, min);
            max = Math::max(patches_boundary_size, max);
          }
          return max - min;
        }

        Index get_cell_deviation() const
        {
          Index min(std::numeric_limits<Index>::max());
          Index max(0);
          for (Index patch(0) ; patch < _num_patches ; ++patch)
          {
            min = Math::min(Index(_cells_per_patch.at(patch).size()), min);
            max = Math::max(Index(_cells_per_patch.at(patch).size()), max);
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
     * -# Refine the mesh until it contains at least num_patches cells.
     * -# Create an object of this class and pass the to-be-partitioned mesh as well
     *    as the desired number of patches and the available number of ranks in a communicator to the constructor.
     * -# Create the Elements-At-Rank graph using the #build_elems_at_rank() function.
     *
     * \author Dirk Ribbrock
     */
    template<typename Shape_, int num_coords_, typename Coord_>
    class PartiIterative<ConformalMesh<Shape_, num_coords_, Coord_>>
    {
      private:
        /// number of elements in input mesh
        const Index _num_elems;
        /// number of usable ranks
        const Index _num_ranks;
        /// Our communicator
        const Dist::Comm & _comm;
        /// Our population
        std::list<Geometry::Intern::PartiIterativeIndividual<Shape_, num_coords_, Coord_>> _population;
        /// number of desired patches
        const Index _num_patches;

      public:
      /// our mesh type
      typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;

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
      explicit PartiIterative(MeshType& mesh, const Dist::Comm & comm, Index num_patches, double time_init, double time_mutate) :
        _num_elems(mesh.get_num_elements()),
        _num_ranks(Index(comm.size())),
        _comm(comm),
        _num_patches(num_patches)
      {
        Random::SeedType seed(Random::SeedType(time(nullptr)));
        Random rng(seed + Random::SeedType(comm.rank()));

        mesh.fill_neighbours();

        {
          Geometry::Intern::PartiIterativeIndividual<Shape_, num_coords_, Coord_> indi(mesh, rng, _num_patches);
          _population.push_back(indi);
        }

        TimeStamp at;
        while(at.elapsed_now() < time_init)
        {
          for (Index i(0) ; i < 10 ; ++i)
          {
            Geometry::Intern::PartiIterativeIndividual<Shape_, num_coords_, Coord_> indi(mesh, rng, _num_patches);
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
          Adjacency::Graph graph(_num_patches, _num_elems, _num_elems);
          Index* ptr = graph.get_domain_ptr();
          Index* idx = graph.get_image_idx();

          // build pointer array
          ptr[0] = 0;
          for(Index i(0); i < _num_patches; ++i)
          {
            ptr[i+1] = Index(indi._cells_per_patch.at(i).size()) + ptr[i];
          }

          // build index array
          Index counter(0);
          for (Index rank(0) ; rank < _num_patches ; ++rank)
          {
            for (auto cell : indi._cells_per_patch.at(rank))
            {
              idx[counter] = cell;
              ++counter;
            }
          }
          _comm.bcast(graph.get_domain_ptr(), _num_patches + 1, _rank);
          _comm.bcast(graph.get_image_idx(), _num_elems, _rank);
          return graph;
        }
        else
        {
          Index * domain_ptr = new Index[_num_patches + 1];
          _comm.bcast(domain_ptr, _num_patches + 1, (int)lowest_rank);
          Index * image_idx = new Index[_num_elems];
          _comm.bcast(image_idx, _num_elems, (int)lowest_rank);
          Adjacency::Graph graph(_num_patches, _num_elems, _num_elems, domain_ptr, image_idx);
          delete[] domain_ptr;
          delete[] image_idx;
          return graph;
        }
      }
    }; // class PartiIterative
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PARTI_ITERATIVE_HPP
