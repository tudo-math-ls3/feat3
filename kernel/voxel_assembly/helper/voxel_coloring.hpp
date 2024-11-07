#pragma once

#include <kernel/base_header.hpp>
#include <kernel/adjacency/graph.hpp>

#include <vector>
namespace FEAT
{
  namespace VoxelAssembly
  {
    template<typename Shape_>
    struct UnitCubeColoring
    {
      static constexpr int dim = Shape_::dimension;

      static std::vector<int> create_coloring(int lvl)
      {
        int num_per_row = 1 << lvl;
        // std::cout << num_per_row << "\n";
        std::vector<int> vec;
        if constexpr(dim == 3)
        {
          vec.resize(std::size_t(num_per_row * num_per_row * num_per_row));
          for(int j = 0; j < num_per_row * num_per_row * num_per_row; ++j)
          {
            vec[std::size_t(j)] = j % 8;

            // std::cout << "\n";
          }
        }
        else
        {
          vec.resize(std::size_t(num_per_row * num_per_row));
          for(int j = 0; j < num_per_row * num_per_row; ++j)
          {
            vec[std::size_t(j)] = j % 4;

            // std::cout << "\n";
          }
        }

        return vec;
      }

      static std::vector<int> naive_coloring(int lvl)
      {
        int num_per_row = 1 << lvl;
        std::vector<int> vec(std::size_t(num_per_row * num_per_row));
        for(int j = 0; j < num_per_row * num_per_row; ++j)
        {
          vec[std::size_t(j)] = j;
          // std::cout << "\n";
        }
        return vec;
      }
    };


    template<typename Mesh_>
    bool test_coloring(const Mesh_& mesh, const std::vector<int>& coloring)
    {
      auto& cell_to_vert = mesh.template get_index_set<Mesh_::world_dim, 0>();
      Adjacency::Graph vert_to_cell(Adjacency::RenderType::transpose, cell_to_vert);
      Adjacency::Graph cell_neighbours(Adjacency::RenderType::injectify, cell_to_vert, vert_to_cell);

      for(Index i = 0; i < cell_neighbours.get_num_nodes_domain(); ++i)
      {
        for(auto ptr = cell_neighbours.image_begin(i); ptr != cell_neighbours.image_end(i); ++ptr)
        {
          if(*ptr != i)
          {
            if(coloring.at(*ptr) == coloring.at(i))
            {
              std::cout << "Same color for neighbour!";
              return false;
            }
          }
        }
      }
      return true;

    }
  }
}
