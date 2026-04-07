// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/voxel_assembly/voxel_assembly_common.hpp>
#include <kernel/voxel_assembly/helper/cell_to_dof_helper.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/adjacency/coloring_data_handler.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/stop_watch.hpp>

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
      Memory::Arbiter _cell_to_dof;
      /// Array mapping each _cell_to_dof entry to it locally sorted position.
      Memory::Arbiter _cell_to_dof_sorter;
      /// Size of cell to dof array
      Index _cell_to_dof_size;

      /// Array of node coordinates.
      Memory::Arbiter _nodes;
      /// Size of node array
      Index _nodes_size;

    public:
      /// Datahandler for the coloring data.
      Adjacency::ColoringDataHandler coloring_data;

    public:
      double time_off_load_mesh, time_off_load_color;
      double time_init_mesh, time_init_color;

    public:
      inline AssemblyMappingData<DataType, IndexType> get_assembly_field() const
      {
        return AssemblyMappingData<DataType, IndexType>{_cell_to_dof, _cell_to_dof_sorter, _cell_to_dof_size, _nodes, _nodes_size};
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
        //const auto& _coloring_maps = coloring_data.get_coloring_maps();
        const Memory::TypedView<int> coloring_view = coloring_data.color_map_view(Memory::Location::main, Memory::Access::read);
        const Memory::TypedView<IndexType> cell_to_dof_view(_cell_to_dof.view(Memory::Location::main, Memory::Access::read));
        const auto& _coloring_map_offsets = coloring_data.get_color_offsets();
        const auto& _coloring_map_sizes = coloring_data.get_color_sizes();
        for(std::size_t i = 0; i < _coloring_map_sizes.size(); ++i)
        {
          for(std::size_t l = 0; l < _coloring_map_sizes.at(i); ++l)
          {
            //int cell_a = _coloring_maps[std::size_t(i)][std::size_t(l)];
            int cell_a = coloring_view[_coloring_map_offsets[i] + l];
            const IndexType* a_b = &cell_to_dof_view[std::size_t(cell_a*SpaceType::DofMappingType::dof_count)];
            const IndexType* a_e = a_b + SpaceType::DofMappingType::dof_count;
            for(std::size_t j = l+1; j < _coloring_map_sizes[i]; ++j)
            {
              //int cell_b = _coloring_maps[std::size_t(i)][std::size_t(j)];
              int cell_b = coloring_view[_coloring_map_offsets[i] + j];
              const IndexType* b_b = &cell_to_dof_view[std::size_t(cell_b*SpaceType::DofMappingType::dof_count)];
              const IndexType* b_e = b_b + SpaceType::DofMappingType::dof_count;

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

          for(std::size_t j = i+1; j < _coloring_map_sizes.size(); ++j)
          {
            //if(_contains_common_element(
            //  _coloring_maps.at(std::size_t(i)),
            //  _coloring_maps.at(std::size_t(i))+_coloring_map_sizes.at(std::size_t(i)),
            //  _coloring_maps.at(std::size_t(j)),
            //  _coloring_maps.at(std::size_t(j))+_coloring_map_sizes.at(std::size_t(j))))
            const int* cvi = &coloring_view.get_r()[_coloring_map_offsets[i]];
            const int* cvj = &coloring_view.get_r()[_coloring_map_offsets[j]];
            if(_contains_common_element(
              cvi, cvi + _coloring_map_sizes.at(i),
              cvj, cvj + _coloring_map_sizes.at(j)))
            {
              std::cout << "Colors contain common element!\n";
              XABORTM("I think you have misunderstood colors.");
            }
          }
        }

        return true;

      }

      void _fill_color(const std::vector<int>& coloring, int hint)
      {
        TimeStamp stamp1;
        coloring_data.fill_color(coloring, hint);
        #ifdef DEBUG
        _test_coloring();
        #endif
        TimeStamp stamp2;
        time_init_color = stamp2.elapsed(stamp1);
      }

      void _fill_color(const Adjacency::Coloring& coloring, int hint)
      {
        TimeStamp stamp1;
        coloring_data.fill_color(coloring, hint);
        #ifdef DEBUG
        _test_coloring();
        #endif
        TimeStamp stamp2;
        time_init_color = stamp2.elapsed(stamp1);
      }

    public:
      explicit LagrangeDataHandler() = default;

      template<typename ColoringType_>
      explicit LagrangeDataHandler(const SpaceType& space, const ColoringType_& coloring, int hint = -1) :
      _cell_to_dof(),
      _cell_to_dof_sorter(),
      _cell_to_dof_size(Index(0)),
      _nodes(),
      _nodes_size(Index(0)),
      time_off_load_mesh(0.),
      time_off_load_color(0.),
      time_init_mesh(0.),
      time_init_color(0.)
      {
        TimeStamp stamp1;
        if constexpr(std::is_same<ColoringType_, Adjacency::Coloring>::value)
        {
          ASSERTM(space.get_mesh().get_num_entities(dim) == coloring.get_num_nodes(), "Coloring and space do not fit!");
        }
        else
        {
          ASSERTM(space.get_mesh().get_num_entities(dim) == coloring.size(), "Coloring and space do not fit!");
        }
        _nodes_size = space.get_mesh().get_vertex_set().get_num_vertices();
        _nodes = Memory::Arbiter(sizeof(Tiny::Vector<DataType, dim>) * _nodes_size);
        // copy the internal nodes array
        const auto* vertex_begin = (const typename SpaceType::MeshType::VertexSetType::CoordType*)space.get_mesh().get_vertex_set().begin();
        std::transform(vertex_begin, vertex_begin + _nodes_size*Index(dim),
          _nodes.template typed_view<DataType>(Memory::Location::main, Memory::Access::write).get_w(),
          [](const auto& a) ->DataType {return DataType(a);});

        // define our _cell_to_dof
        _cell_to_dof_size = space.get_mesh().get_num_entities(dim) * SpaceType::DofMappingType::dof_count;
        _cell_to_dof = Memory::Arbiter(sizeof(IndexType) * std::size_t(_cell_to_dof_size));
        _cell_to_dof_sorter = Memory::Arbiter(sizeof(IndexType) * std::size_t(_cell_to_dof_size));
        // for this iterate through our target_sets and parse them in
        {
          Memory::TypedView<IndexType> ctdv(_cell_to_dof.view(Memory::Location::main, Memory::Access::write));
          Memory::TypedView<IndexType> ctdsv(_cell_to_dof_sorter.view(Memory::Location::main, Memory::Access::write));
          VoxelAssembly::fill_cell_to_dof(ctdv.get_w(), space);
          VoxelAssembly::fill_sorter(ctdsv.get_w(), ctdv.get_r(), space);
        }
        TimeStamp stamp2;
        time_init_mesh = stamp2.elapsed(stamp1);

        // for(int cell = 0; cell < space.get_mesh().get_num_elements(); ++cell)
        // {
        //   std::cout << "For cell " << cell << ": \n";
        //   IndexType* ctd = _cell_to_dof.data() + cell * SpaceType::DofMappingType::dof_count;
        //   for(int i = 0; i < SpaceType::DofMappingType::dof_count; ++i)
        //   {
        //     std::cout << ctd[i] << "  ";
        //   }
        //   std::cout << "\n";
        // }

        _fill_color(coloring, hint);

      }

      LagrangeDataHandler(const LagrangeDataHandler&) = delete;

      LagrangeDataHandler& operator=(const LagrangeDataHandler&) = delete;

      LagrangeDataHandler(LagrangeDataHandler&& other) noexcept :
      _cell_to_dof(std::move(other._cell_to_dof)),
      _cell_to_dof_sorter(std::move(other._cell_to_dof_sorter)),
      _cell_to_dof_size(other._cell_to_dof_size),
      _nodes(std::move(other._nodes)),
      _nodes_size(other._nodes_size),
      coloring_data(std::move(other.coloring_data)),
      time_off_load_mesh(other.time_off_load_mesh),
      time_off_load_color(other.time_off_load_color),
      time_init_mesh(other.time_init_mesh),
      time_init_color(other.time_init_color)
      {
      }

      LagrangeDataHandler& operator=(LagrangeDataHandler&& other) noexcept
      {
        if(this == &other)
          return *this;
        time_off_load_mesh = other.time_off_load_mesh;
        time_off_load_color = other.time_off_load_color;
        time_init_mesh = other.time_init_mesh;
        time_init_color = other.time_init_color;
        coloring_data = std::move(other.coloring_data);
        _cell_to_dof = std::move(other._cell_to_dof);
        _cell_to_dof_sorter = std::move(other._cell_to_dof_sorter);
        _cell_to_dof_size = other._cell_to_dof_size;
        _nodes = std::move(other._nodes);
        _nodes_size = other._nodes_size;
        return *this;
      }

      ~LagrangeDataHandler() = default;
    }; //class GPUDataHandler
  }
}
