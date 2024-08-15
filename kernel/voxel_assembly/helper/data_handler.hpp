#pragma once
#ifndef FEAT_KERNEL_VOXEL_ASSEMBLY_HELPER_DATA_HANDLER_HPP
#define FEAT_KERNEL_VOXEL_ASSEMBLY_HELPER_DATA_HANDLER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/voxel_assembly/voxel_assembly_common.hpp>
#include <kernel/voxel_assembly/helper/cell_to_dof_helper.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/adjacency/coloring.hpp>

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

#include <vector>

namespace FEAT
{
  namespace VoxelAssembly
  {
    //Warning: for now only works with Q2standard space... we should probably use this class as specialzation for lagrange based spaces...
    /**
     * \brief Data handler for Lagrange based FE spaces
     *
     * This class handles the data extraction of the underlying mesh, space and color data
     * for voxel assembly.
     * \note Only works for Lagrange space because it is implicitly required that the mapping of the
     * first k dofs are matched with the actual node numeration.
     *
     * \tparam Space_ The underlying space type.
     * \tparam DT_ The datatype.
     * \tparam IT_ The indextype.
     */
    template<typename Space_, typename DT_, typename IT_>
    class LagrangeDataHandler
    {
    public:
      typedef Space_ SpaceType;
      typedef DT_ DataType;
      typedef IT_ IndexType;

      static constexpr int dim = SpaceType::world_dim;
    protected:
      /// Array mapping each cell index to ALL its local dofs in correct numbering.
      IndexType* _cell_to_dof;
      /// Array mapping each _cell_to_dof entry to it locally sorted position.
      IndexType* _cell_to_dof_sorter;
      /// Size of cell to dof array
      Index _cell_to_dof_size;

      /// Array of node coordinates.
      Tiny::Vector<DataType, dim>* _nodes;
      /// Size of node array
      Index _nodes_size;

      /// Datahanlder for the coloring data.
      Adjacency::ColoringDataHandler coloring_data;

    public:
      inline AssemblyMappingData<DataType, IndexType> get_assembly_field() const
      {
        return AssemblyMappingData<DataType, IndexType>{_cell_to_dof, _cell_to_dof_sorter, _cell_to_dof_size,
                                    (void*)_nodes, _nodes_size};
      }

      Index get_num_colors() const
      {
        return coloring_data.get_num_colors();
      }

      int* get_color_map(Index k)
      {
        return coloring_data.get_color_map(k);
      }

      const int* get_color_map(Index k) const
      {
        return coloring_data.get_color_map(k);
      }

      std::vector<int*>& get_coloring_maps()
      {
        return coloring_data.get_coloring_maps();
      }

      const std::vector<int*>& get_coloring_maps() const
      {
        return coloring_data.get_coloring_maps();
      }

      Index get_color_map_size(Index k) const
      {
        return coloring_data.get_color_size(k);
      }

      std::vector<Index>& get_color_map_sizes()
      {
        return coloring_data.get_color_sizes();
      }

      const std::vector<Index>& get_color_map_sizes() const
      {
        return coloring_data.get_color_sizes();
      }

    protected:
      //Container has to be sortable and a and b have to be sorted beforhand....
      template<typename ITX_>
      bool _contains_common_element(const ITX_* a, const ITX_* ae, const ITX_* b, const ITX_* be) const
      {
          std::vector<ITX_> intersection;

          std::set<ITX_> aa { a, ae };
          std::set<ITX_> bb { b, be };

          std::set_intersection(aa.begin(), aa.end(),
                                bb.begin(), bb.end(),
                                std::inserter(intersection, intersection.end()));

          return !intersection.empty();
      }

      bool _test_coloring() const
      {
        const auto& _coloring_maps = coloring_data.get_coloring_maps();
        const auto& _coloring_map_sizes = coloring_data.get_color_sizes();
        for(int i = 0; i < int(_coloring_maps.size()); ++i)
        {
          for(int l = 0; l < int(_coloring_map_sizes.at(std::size_t(i))); ++l)
          {
            int cell_a = _coloring_maps[std::size_t(i)][std::size_t(l)];
            const IndexType* a_b = &_cell_to_dof[std::size_t(cell_a*SpaceType::DofMappingType::dof_count)];
            const IndexType* a_e = &_cell_to_dof[std::size_t(cell_a*SpaceType::DofMappingType::dof_count)] + SpaceType::DofMappingType::dof_count;
            for(int j = l+1; j < int(_coloring_map_sizes[std::size_t(i)]); ++j)
            {
              int cell_b = _coloring_maps[std::size_t(i)][std::size_t(j)];
              const IndexType* b_b = &_cell_to_dof[std::size_t(cell_b*SpaceType::DofMappingType::dof_count)];
              const IndexType* b_e = &_cell_to_dof[std::size_t(cell_b*SpaceType::DofMappingType::dof_count)] + SpaceType::DofMappingType::dof_count;

              if(_contains_common_element(a_b, a_e, b_b, b_e))
              {
                std::cout << "Cell 1: ";
                for(int r = 0; r < SpaceType::DofMappingType::dof_count; ++r)
                {
                  std::cout << *(a_b + r) << " ";
                }
                std::cout << "\nCell 2: ";
                for(int r = 0; r < SpaceType::DofMappingType::dof_count; ++r)
                {
                  std::cout << *(b_b + r) << " ";
                }
                std::cout << "\n";
                XABORTM("Intersection in color " + stringify(i) + " between cells " + stringify(cell_a) + " " + stringify(cell_b));
              }
            }
          }

          for(int j = i+1; j < int(_coloring_maps.size()); ++j)
          {
            if(_contains_common_element(_coloring_maps.at(std::size_t(i)), _coloring_maps.at(std::size_t(i))+_coloring_map_sizes.at(std::size_t(i)), _coloring_maps.at(std::size_t(j)), _coloring_maps.at(std::size_t(j))+_coloring_map_sizes.at(std::size_t(j))))
            {
              std::cout << "Colors contain common element!" << "\n";
              XABORTM("I think you have misunderstood colors.");
            }
          }
        }

        return true;

      }

      void _fill_color(const std::vector<int>& coloring, int hint)
      {
        coloring_data.fill_color(coloring, hint);
        #ifdef DEBUG
        _test_coloring();
        #endif
      }

      void _fill_color(const Adjacency::Coloring& coloring, int hint)
      {
        coloring_data.fill_color(coloring, hint);
        #ifdef DEBUG
        _test_coloring();
        #endif
      }

    public:
      explicit LagrangeDataHandler() = default;

      template<typename ColoringType_>
      explicit LagrangeDataHandler(const SpaceType& space, const ColoringType_& coloring, int hint = -1) :
      _cell_to_dof(nullptr),
      _cell_to_dof_sorter(nullptr),
      _cell_to_dof_size(Index(0)),
      _nodes(nullptr),
      _nodes_size(Index(0))
      {
        if constexpr(std::is_same<ColoringType_, Adjacency::Coloring>::value)
        {
          ASSERTM(space.get_mesh().get_num_entities(dim) == coloring.get_num_nodes(), "Coloring and space do not fit!");
        }
        else
        {
          ASSERTM(space.get_mesh().get_num_entities(dim) == coloring.size(), "Coloring and space do not fit!");
        }
        _nodes_size = space.get_mesh().get_vertex_set().get_num_vertices();
        _nodes = MemoryPool::allocate_memory<Tiny::Vector<DataType, dim>>(_nodes_size);
        // copy the internal nodes array
        const auto* vertex_begin = (const typename SpaceType::MeshType::VertexSetType::CoordType*)space.get_mesh().get_vertex_set().begin();
        std::transform(vertex_begin, vertex_begin + _nodes_size*Index(dim), (DataType*)_nodes, [](const auto& a) ->DataType {return DataType(a);});

        // define our _cell_to_dof
        _cell_to_dof_size = space.get_mesh().get_num_entities(dim) * SpaceType::DofMappingType::dof_count;
        _cell_to_dof = MemoryPool::allocate_memory<IndexType>(_cell_to_dof_size);
        _cell_to_dof_sorter = MemoryPool::allocate_memory<IndexType>(_cell_to_dof_size);
        // for this iterate through our target_sets and parse them in
        VoxelAssembly::fill_cell_to_dof(_cell_to_dof, space);
        VoxelAssembly::fill_sorter(_cell_to_dof_sorter, _cell_to_dof, space);

        // for(int cell = 0; cell < space.get_mesh().get_num_elements(); ++cell)
        // {
        //   std::cout << "For cell " << cell << ": \n";
        //   IndexType* ctd = _cell_to_dof.data() + cell * SpaceType::DofMappingType::dof_count;
        //   for(int i = 0; i < SpaceType::DofMappingType::dof_count; ++i)
        //   {
        //     std::cout << ctd[i] << "  ";
        //   }
        //   std::cout << std::endl;
        // }

        _fill_color(coloring, hint);

      }

      LagrangeDataHandler(const LagrangeDataHandler&) = delete;

      LagrangeDataHandler& operator=(const LagrangeDataHandler&) = delete;

      LagrangeDataHandler(LagrangeDataHandler&& other) noexcept :
      _cell_to_dof(other._cell_to_dof),
      _cell_to_dof_sorter(other._cell_to_dof_sorter),
      _cell_to_dof_size(other._cell_to_dof_size),
      _nodes(other._nodes),
      _nodes_size(other._nodes_size),
      coloring_data(std::move(other.coloring_data))
      {
        MemoryPool::increase_memory(_cell_to_dof);
        MemoryPool::increase_memory(_cell_to_dof_sorter);
        MemoryPool::increase_memory(_nodes);
      }

      LagrangeDataHandler& operator=(LagrangeDataHandler&& other) noexcept
      {
        if(this == &other)
          return *this;
        coloring_data = std::move(other.coloring_data);
        MemoryPool::release_memory(_nodes);
        MemoryPool::release_memory(_cell_to_dof_sorter);
        MemoryPool::release_memory(_cell_to_dof);
        _cell_to_dof = other._cell_to_dof;
        _cell_to_dof_sorter = other._cell_to_dof_sorter;
        _cell_to_dof_size = other._cell_to_dof_size;
        _nodes = other._nodes;
        _nodes_size = other._nodes_size;
        MemoryPool::increase_memory(_cell_to_dof);
        MemoryPool::increase_memory(_cell_to_dof_sorter);
        MemoryPool::increase_memory(_nodes);
        return *this;
      }

      ~LagrangeDataHandler()
      {
        MemoryPool::release_memory(_nodes);
        MemoryPool::release_memory(_cell_to_dof_sorter);
        MemoryPool::release_memory(_cell_to_dof);
      }
    }; //class GPUDataHandler
  }
}

#endif // FEAT_KERNEL_VOXEL_ASSEMBLY_HELPER_DATA_HANDLER_HPP