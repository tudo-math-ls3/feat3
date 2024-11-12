// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/shape.hpp>

#include <algorithm>

namespace FEAT
{
  namespace VoxelAssembly
  {
    namespace Intern
    {
      /// Prescribes the number of dofs per dimension for a given space
      template<typename SpaceType_, int target_dim_>
      struct DofHelper;

      template<int target_dim>
      struct DofHelper<VoxelAssembly::Q2StandardQuad, target_dim>
      {
        static constexpr int dof_per_entity = 1;
      };

      template<int target_dim>
      struct DofHelper<VoxelAssembly::Q2StandardHexa, target_dim>
      {
        static constexpr int dof_per_entity = 1;
      };

      template<int num_indices>
      Geometry::IndexTuple<num_indices> sort_index_tuple(const Geometry::IndexTuple<num_indices>& tuple)
      {
        Geometry::IndexTuple<num_indices> tmp = tuple;
        Index* begin = tmp.indices;
        Index* end = tmp.indices + num_indices;
        std::sort(begin, end);
        return tmp;
      }
    }

    /**
     * \brief Computes a sequentialized dof mapping for the given space
     *
     * \warning Should only be called with tar_dim = prior_ents = offset = 0
     */
    template<typename Space_, typename IT_, int tar_dim = 0, int prior_ents = 0>
    void fill_cell_to_dof(IT_* ctd, const Space_& space, IT_ offset = 0)
    {
      constexpr int dim = Space_::shape_dim;
      constexpr int dof_per_cell = Space_::DofMappingType::dof_count;
      //get index_set dim -> tar_dim
      const auto& mesh = space.get_mesh();
      const auto& index_set = mesh.template get_index_set<dim, tar_dim>();
      constexpr int num_indices = Shape::FaceTraits<typename Space_::MeshType::ShapeType, tar_dim>::count;
      /// How many dofs do we have per target dim?
      constexpr int dof_per_entity = Intern::DofHelper<Space_, tar_dim>::dof_per_entity;
      for(Index ent = 0; ent < index_set.get_num_entities(); ++ent)
      {
        const auto ind_arr = index_set[ent];
        IT_* tmp_ctd = ctd + ent * dof_per_cell;
        for(Index i = 0; i < ind_arr.num_indices; ++i)
        {
          for(IT_ j = 0; j < IT_(dof_per_entity); ++j)
          {
            tmp_ctd[i*IT_(dof_per_entity) + j + IT_(prior_ents)] = offset + IT_(dof_per_entity) * IT_(ind_arr[int(i)]) + j;
          }
        }
      }
      offset += IT_(mesh.get_num_entities(tar_dim)*dof_per_entity);
      if constexpr(tar_dim < (dim-1))
      {
        constexpr int new_prior_ents = prior_ents + num_indices*dof_per_entity;
        fill_cell_to_dof<Space_, IT_, tar_dim+1, new_prior_ents>(ctd, space, offset);
      }
      else if constexpr(tar_dim == (dim-1))
      {
        constexpr int shape_dim_dofs_per = Intern::DofHelper<Space_, dim>::dof_per_entity;
        for(Index ent = 0; ent < index_set.get_num_entities(); ++ent)
        {
          IT_* tmp_ctd = ctd + (ent+1)*dof_per_cell - shape_dim_dofs_per;
          for(int j = 0; j < shape_dim_dofs_per; ++j)
            tmp_ctd[j] = IT_(ent)*IT_(shape_dim_dofs_per) + IT_(j) + offset;
        }
      }
      else
      {
        XABORTM("Thou shallt not arrive here!");
      }
    }

    template<typename Space_, typename IT_>
    void fill_sorter(IT_* sorter, const IT_* mapping, const Space_& space)
    {
      constexpr int num_loc_dofs = Space_::DofMappingType::dof_count;
      for(Index cell = 0; cell < space.get_mesh().get_num_entities(Space_::world_dim); ++cell)
      {
        //fill up local range
        const Index begin = cell*num_loc_dofs;
        const Index end = (cell+1)*num_loc_dofs;
        const IT_* loc_map = mapping + begin;
        std::generate(sorter+begin, sorter+end, [n=0]() mutable{return n++;});
        std::sort(sorter+begin, sorter+end, [loc_map](IT_ va, IT_ vb){return loc_map[va] < loc_map[vb];});
      }

      /* For testing...
      for(Index cell = 0; cell < space.get_mesh().get_num_entities(Space_::world_dim); ++cell)
      {
        //fill up local range
        const Index begin = cell*num_loc_dofs;
        const Index end = (cell+1)*num_loc_dofs;
        const Index* loc_map = mapping + begin;
        std::cout << "Mapping at cell " << cell << "\n";
        for(Index i = begin; i < end; ++i)
        {
          std::cout << loc_map[i] << " ";
        }
        std::cout << "\n";
        for(Index i = begin; i < end; ++i)
        {
          std::cout << sorter[i] << " ";
        }
        std::cout << "\n";
      } //*/
    }
  }
}
