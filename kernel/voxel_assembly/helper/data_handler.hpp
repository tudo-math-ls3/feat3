#pragma once
#ifndef FEAT_KERNEL_VOXEL_ASSEMBLY_HELPER_DATA_HANDLER_HPP
#define FEAT_KERNEL_VOXEL_ASSEMBLY_HELPER_DATA_HANDLER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
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
      /// Vector mapping each cell index to ALL its local dofs in correct numbering.
      std::vector<IndexType> _cell_to_dof;
      /// Vector mapping each _cell_to_dof entry to it locally sorted position.
      std::vector<IndexType> _cell_to_dof_sorter;

      /// Device pointer pointing to device data holding the same information as _cell_to_dof.
      void* _d_cell_to_dof = nullptr;
      /// Device pointer to _cell_to_dof_sorter representation.
      void* _d_cell_to_dof_sorter = nullptr;

      /// Vector of node coordinates.
      std::vector<Tiny::Vector<DataType, dim>> _nodes;
      /// device pointer to array of nodes. Node data is saved sequentially
      void* _d_nodes = nullptr;

      /// The resulting coloring mapping, i.e. each std::vector contains the cell indices of one color
      std::vector<std::vector<int>> _coloring_maps;
      /// Device pointer to the single coloring maps
      std::vector<void*> _d_coloring_maps;

    public:
      inline AssemblyMappingData<DataType, IndexType> get_host_assembly_field() const
      {
        return AssemblyMappingData<DataType, IndexType>{_cell_to_dof.data(), _cell_to_dof_sorter.data(), _cell_to_dof.size(),
                                    (void*)_nodes.data(), _nodes.size()};
      }

      inline AssemblyMappingData<DataType, IndexType> get_device_assembly_field() const
      {
        return AssemblyMappingData<DataType, IndexType>{(IT_*)_d_cell_to_dof, (IT_*)_d_cell_to_dof_sorter, _cell_to_dof.size(),
                                    _d_nodes, _nodes.size()};
      }

      Index get_num_colors() const
      {
        return _coloring_maps.size();
      }

      std::vector<int>& get_color_map(Index k)
      {
        return _coloring_maps.at(k);
      }

      const std::vector<int>& get_color_map(Index k) const
      {
        return _coloring_maps.at(k);
      }

      void* get_color_map_device(Index k)
      {
        return _d_coloring_maps.at(k);
      }

      const void* get_color_map_device(Index k) const
      {
        return _d_coloring_maps.at(k);
      }

      std::vector<std::vector<int>>& get_coloring_maps()
      {
        return _coloring_maps;
      }

      const std::vector<std::vector<int>>& get_coloring_maps() const
      {
        return _coloring_maps;
      }

      std::vector<void*>& get_coloring_maps_device()
      {
        return _d_coloring_maps;
      }

      const std::vector<void*>& get_coloring_maps_device() const
      {
        return _d_coloring_maps;
      }

    protected:
      inline void _init_device()
      {
        #ifdef FEAT_HAVE_CUDA
          if(_d_cell_to_dof != nullptr || _d_nodes != nullptr || _d_coloring_maps.size() > 0)
            XABORTM("Device Memory already initialized!");
          _d_cell_to_dof = Util::cuda_malloc(_cell_to_dof.size()*sizeof(IndexType));
          Util::cuda_copy_host_to_device(_d_cell_to_dof, (void*)_cell_to_dof.data(), _cell_to_dof.size()*sizeof(IndexType));

          _d_cell_to_dof_sorter = Util::cuda_malloc(_cell_to_dof_sorter.size()*sizeof(IndexType));
          Util::cuda_copy_host_to_device(_d_cell_to_dof_sorter, (void*)_cell_to_dof_sorter.data(), _cell_to_dof_sorter.size()*sizeof(IndexType));

          _d_nodes = Util::cuda_malloc(_nodes.size() * dim * sizeof(DataType));
          Util::cuda_copy_host_to_device(_d_nodes, (void*)_nodes.data(), _nodes.size() * dim * sizeof(DataType));

          _d_coloring_maps.resize(_coloring_maps.size());
          for(Index i = 0; i < _coloring_maps.size(); ++i)
          {
            _d_coloring_maps[i] = Util::cuda_malloc(_coloring_maps[i].size() * sizeof(int));
            Util::cuda_copy_host_to_device(_d_coloring_maps[i], _coloring_maps[i].data(), _coloring_maps[i].size() * sizeof(int));
          }
        #endif
      }

      inline void _free_device()
      {
        #ifdef FEAT_HAVE_CUDA
          for(Index i = 0; i < _d_coloring_maps.size(); ++i)
          {
            Util::cuda_free(_d_coloring_maps[i]);
          }
          _d_coloring_maps.clear();
          Util::cuda_free(_d_nodes);
          _d_nodes = nullptr;
          Util::cuda_free(_d_cell_to_dof_sorter);
          _d_cell_to_dof_sorter = nullptr;
          Util::cuda_free(_d_cell_to_dof);
          _d_cell_to_dof = nullptr;
        #endif
      }

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
        for(int i = 0; i < int(_coloring_maps.size()); ++i)
        {
          for(int l = 0; l < int(_coloring_maps[std::size_t(i)].size()); ++l)
          {
            int cell_a = _coloring_maps[std::size_t(i)][std::size_t(l)];
            const IndexType* a_b = &_cell_to_dof[std::size_t(cell_a*SpaceType::DofMappingType::dof_count)];
            const IndexType* a_e = &_cell_to_dof[std::size_t(cell_a*SpaceType::DofMappingType::dof_count)] + SpaceType::DofMappingType::dof_count;
            for(int j = l+1; j < int(_coloring_maps[std::size_t(i)].size()); ++j)
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
                std::cout << std::endl << "Cell 2: ";
                for(int r = 0; r < SpaceType::DofMappingType::dof_count; ++r)
                {
                  std::cout << *(b_b + r) << " ";
                }
                std::cout << std::endl;
                XABORTM("Intersection in color " + stringify(i) + " between cells " + stringify(cell_a) + " " + stringify(cell_b));
              }
            }
          }

          for(int j = i+1; j < int(_coloring_maps.size()); ++j)
          {
            if(_contains_common_element(&_coloring_maps.at(std::size_t(i)).front(), &_coloring_maps.at(std::size_t(i)).back()+1, &_coloring_maps.at(std::size_t(j)).front(), &_coloring_maps.at(std::size_t(j)).back()+1))
            {
              std::cout << "Colors contain common element!" << std::endl;
              XABORTM("I think you have misunderstood colors.");
            }
          }
        }

        return true;

      }

      void _fill_color(const std::vector<int>& coloring, int hint)
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
        #ifdef DEBUG
        _test_coloring();
        #endif

      }

      void _fill_color(const Adjacency::Coloring& coloring, int hint)
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
        #ifdef DEBUG
        _test_coloring();
        #endif

      }

    public:
      explicit LagrangeDataHandler() = default;

      template<typename ColoringType_>
      explicit LagrangeDataHandler(const SpaceType& space, const ColoringType_& coloring, int hint = -1) :
      _d_cell_to_dof{nullptr},
      _d_cell_to_dof_sorter{nullptr},
      _d_nodes{nullptr}
      {
        if constexpr(std::is_same<ColoringType_, Adjacency::Coloring>::value)
        {
          ASSERTM(space.get_mesh().get_num_entities(dim) == coloring.get_num_nodes(), "Coloring and space do not fit!");
        }
        else
        {
          ASSERTM(space.get_mesh().get_num_entities(dim) == coloring.size(), "Coloring and space do not fit!");
        }
        _cell_to_dof.resize(space.get_mesh().get_num_entities(dim) * SpaceType::DofMappingType::dof_count);
        _cell_to_dof_sorter.resize(space.get_mesh().get_num_entities(dim) * SpaceType::DofMappingType::dof_count);
        _nodes = space.get_mesh().get_vertex_set().template clone_internal_vector<DataType>();
        // define our _cell_to_dof
        //for this iterate through our target_sets and parse them in
        VoxelAssembly::fill_cell_to_dof(_cell_to_dof.data(), space);
        VoxelAssembly::fill_sorter(_cell_to_dof_sorter.data(), _cell_to_dof.data(), space);

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



#ifdef FEAT_HAVE_CUDA
        if(Backend::get_preferred_backend() == PreferredBackend::cuda)
          _init_device();
#endif
      }

      ~LagrangeDataHandler()
      {
        _free_device();
      }
    }; //class GPUDataHandler
  }
}

#endif // FEAT_KERNEL_VOXEL_ASSEMBLY_HELPER_DATA_HANDLER_HPP