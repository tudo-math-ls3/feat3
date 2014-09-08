#pragma once
#ifndef KERNEL_FOUNDATION_AURA_HPP
#define KERNEL_FOUNDATION_AURA_HPP 1

#include<vector>
#include<kernel/archs.hpp>
#include<kernel/foundation/halo.hpp>

namespace FEAST
{
  namespace Foundation
  {
    template<typename M_, typename A_, typename HT_>
    struct Aura
    {
    };

    template<typename HT_>
    struct Aura<Mem::Main, Algo::Generic, HT_>
    {
      template<typename MeshType_,
               typename WT_,
               template<typename, typename> class StorageT_>
      static HT_ value(const StorageT_<std::shared_ptr<HaloBase<MeshType_, WT_, StorageT_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_, WT_, StorageT_> > > >& halos)
      {
        //retrieve deltas of all halos
        StorageT_<unsigned, std::allocator<unsigned> > delta;
        for(auto& h_i : halos)
          delta.push_back(h_i->get_overlap());
        for(auto d_i : delta)
          ASSERT(d_i == delta.at(0), "Error: Divergent delta values of halos in parameter pack!");

        //retrieve mesh pointers of all halos
        StorageT_<typename HT_::mesh_type_*, std::allocator<typename HT_::mesh_type_*> > meshes;
        for(auto& h_i : halos)
          meshes.push_back(h_i->get_mesh());
        for(auto m_i : meshes)
          ASSERT(m_i == meshes.at(0), "Error: Divergent mesh pointer values of halos in parameter pack!");
        const typename HT_::mesh_type_* mesh(meshes.at(0));

        //retrieve other ranks of all halos
        StorageT_<typename HT_::mesh_type_::index_type_, std::allocator<typename HT_::mesh_type_::index_type_> > others;
        for(auto& h_i : halos)
          others.push_back(h_i->get_other());
        for(auto o_i : others)
          ASSERT(o_i == others.at(0), "Error: Divergent comm links of halos in parameter pack!");

        //retrieve levels of all halos
        StorageT_<PolytopeLevels, std::allocator<PolytopeLevels> > pl;
        for(auto& h_i : halos)
          pl.push_back(h_i->get_level());

        //retrieve single pointers to data of all halos
        StorageT_<typename HT_::compound_storage_type_, std::allocator<typename HT_::compound_storage_type_> > datas;
        for(auto& h_i : halos)
          datas.push_back(h_i->get_elements());

        HT_ result;
        result.reset_mesh(meshes.at(0));
        result.reset_other(others.at(0));

        PolytopeLevels target_level(result.get_level());
        typename HT_::compound_storage_type_ target_data;

        for(typename HT_::mesh_type_::index_type_ i(0) ; i < halos.size() ; ++i)
        {
          typename HT_::compound_storage_type_& halo_data(datas.at(i));
          if(target_level == pl.at(i))
          {
            std::sort(halo_data.begin(), halo_data.end());
            typename HT_::compound_storage_type_ tmp(halo_data.size() + target_data.size());
            auto it(std::set_union(target_data.begin(), target_data.end(), halo_data.begin(), halo_data.end(), tmp.begin()));
            tmp.resize(Index(it - tmp.begin()));
            target_data = std::move(tmp);
          }
          else if(target_level < pl.at(i))
          {
            for(auto j : halo_data)
            {
              auto adj(mesh->get_adjacent_polytopes(pl.at(i), target_level, j));

              std::sort(adj.begin(), adj.end());
              typename HT_::compound_storage_type_ tmp(adj.size() + target_data.size());
              auto it(std::set_union(target_data.begin(), target_data.end(), adj.begin(), adj.end(), tmp.begin()));
              tmp.resize(Index(it - tmp.begin()));
              target_data = std::move(tmp);
            }
          }
        }

        result.set_elements(std::move(target_data));
        return result;
      }

      /*template<typename... T_>
      static HT_ value(T_&... halos)
      {
        //retrieve deltas of all halos
        std::vector<unsigned, std::allocator<unsigned> > delta;
        int r0[sizeof...(halos)] = {(delta.push_back(halos.get_overlap()), 0)...}; r0;
        for(auto d_i : delta)
          ASSERT(d_i == delta.at(0), "Error: Divergent delta values of halos in parameter pack!");

        //retrieve mesh pointers of all halos
        std::vector<typename HT_::mesh_type_*> meshes;
        int r3[sizeof...(halos)] = {(meshes.push_back(halos.get_mesh()), 0)...}; r3;
        for(auto m_i : meshes)
          ASSERT(m_i == meshes.at(0), "Error: Divergent mesh pointer values of halos in parameter pack!");
        const typename HT_::mesh_type_* mesh(meshes.at(0));

        //retrieve other ranks of all halos
        std::vector<typename HT_::mesh_type_::index_type_, std::allocator<typename HT_::mesh_type_::index_type_> > others;
        int r2[sizeof...(halos)] = {(others.push_back(halos.get_other()), 0)...}; r2;
        for(auto o_i : others)
          ASSERT(o_i == others.at(0), "Error: Divergent comm links of halos in parameter pack!");

        //retrieve levels of all halos
        std::vector<PolytopeLevels, std::allocator<PolytopeLevels> > pl;
        int r1[sizeof...(halos)] = {(pl.push_back(halos.get_level()), 0)...}; r1;

        //retrieve single pointers to data of all halos
        std::vector<typename HT_::compound_storage_type_> datas;
        int r4[sizeof...(halos)] = {(datas.push_back(halos.get_elements()), 0)...}; r4;

        HT_ result;
        result.reset_mesh(meshes.at(0));
        result.reset_other(others.at(0));

        PolytopeLevels target_level(result.get_level());
        typename HT_::compound_storage_type_ target_data;

        for(typename HT_::mesh_type_::index_type_ i(0) ; i < sizeof...(halos) ; ++i)
        {
          typename HT_::compound_storage_type_& halo_data(datas.at(i));
          if(target_level == pl.at(i))
          {
            std::sort(halo_data.begin(), halo_data.end());
            typename HT_::compound_storage_type_ tmp(halo_data.size() + target_data.size());
            auto it(std::set_union(target_data.begin(), target_data.end(), halo_data.begin(), halo_data.end(), tmp.begin()));
            tmp.resize(Index(it - tmp.begin()));
            target_data = std::move(tmp);
          }
          else if(target_level < pl.at(i))
          {
            for(auto j : halo_data)
            {
              auto adj(mesh->get_adjacent_polytopes(pl.at(i), target_level, j));

              std::sort(adj.begin(), adj.end());
              typename HT_::compound_storage_type_ tmp(adj.size() + target_data.size());
              auto it(std::set_union(target_data.begin(), target_data.end(), adj.begin(), adj.end(), tmp.begin()));
              tmp.resize(Index(it - tmp.begin()));
              target_data = std::move(tmp);
            }
          }
        }

        result.set_elements(std::move(target_data));
        return result;
      }*/

    };

  }
}
#endif
