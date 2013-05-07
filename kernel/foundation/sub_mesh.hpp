#pragma once
#ifndef KERNEL_FOUNDATION_SUBMESH_HPP
#define KERNEL_FOUNDATION_SUBMESH_HPP 1


#include "kernel/foundation/mesh.hpp"
#include "kernel/foundation/halo.hpp"

namespace FEAST
{
namespace Foundation
{

  template<
    typename Dim_,
    typename TopologyType_,
    template <typename, typename> class OuterStorageType_
    >
  class SubMesh : public Mesh<Dim_, TopologyType_, OuterStorageType_>
  {
    public:
      typedef typename Mesh<Dim_, TopologyType_, OuterStorageType_>::topology_type_::storage_type_ map_storage_type_;

      SubMesh(const HaloBase<Mesh<Dim_, TopologyType_, OuterStorageType_>, OuterStorageType_>* proxy) :
        Mesh<Dim_, TopologyType_, OuterStorageType_>(),
        _proxy(proxy),
        _map(),
        _level(proxy->get_level())
      {
        TopologyType_ subsets(Index(_level) + 1);

        for(Index i(0) ; i < _proxy->size() ; ++i)
        {
          this->add_polytope(_level);
          subsets.at(Index(_level)).push_back(_proxy->get_element(i));
          _map.push_back(_proxy->get_element(i));

          Index l(0);
          while(l < Index(_level))
          {
            typename TopologyType_::storage_type_ adjacent_polytopes(_proxy->get_mesh()->get_adjacent_polytopes(PolytopeLevels(_level),
                                                                                                                PolytopeLevels(l),
                                                                                                                _proxy->get_element(i)
                                                                                                                ));

            Index rsize(subsets.at(l).size() + adjacent_polytopes.size());
            typename TopologyType_::storage_type_ result_union(rsize);

            std::sort(adjacent_polytopes.begin(), adjacent_polytopes.end());

            typename TopologyType_::storage_type_::iterator it;

            it = std::set_union(subsets.at(l).begin(), subsets.at(l).end(), adjacent_polytopes.begin(), adjacent_polytopes.end(), result_union.begin());

            result_union.resize(it - result_union.begin());

            subsets.at(l) = result_union;

            ++l;
          }
        }

        for(Index i(0) ; i < Index(_level) ; ++i)
        {
          for(Index j(0) ; j < subsets.at(i).size() ; ++j)
          {
            this->add_polytope(PolytopeLevels(i));
          }
        }

        typename TopologyType_::storage_type_ vertex_map;
        Index l(1);
        while (l <= Index(_level))
        {
          for(Index i(0) ; i < subsets.at(l).size() ; ++i)
          {
            typename TopologyType_::storage_type_ adjacent_polytopes(_proxy->get_mesh()->get_adjacent_polytopes(PolytopeLevels(l),
                                                                                                                pl_vertex,
                                                                                                                subsets.at(l).at(i) ));

            for(auto v : adjacent_polytopes)
            {
              if(std::find(vertex_map.begin(), vertex_map.end(), v) != vertex_map.end())
              {
                ///v has already been mapped
                typename TopologyType_::storage_type_::iterator pos(std::find(vertex_map.begin(), vertex_map.end(), v));
                this->add_adjacency(PolytopeLevels(l), pl_vertex, i, pos - vertex_map.begin());
              }
              else
              {
                ///v has not yet been mapped
                this->add_adjacency(PolytopeLevels(l), pl_vertex, i, vertex_map.size());
                vertex_map.push_back(v);
              }
            }
          }

          ++l;
        }
      }

      const PolytopeLevels get_level() const
      {
        return _level;
      }

      typename TopologyType_::storage_type_ get_map()
      {
        return _map;
      }

      const typename TopologyType_::storage_type_ get_map() const
      {
        return _map;
        }

        const HaloBase<Mesh<Dim_, TopologyType_, OuterStorageType_>, OuterStorageType_>* get_proxy() const
        {
          return _proxy;
        }

      private:
        const HaloBase<Mesh<Dim_, TopologyType_, OuterStorageType_>, OuterStorageType_>* _proxy;
        typename Mesh<Dim_, TopologyType_, OuterStorageType_>::topology_type_::storage_type_ _map;
        PolytopeLevels _level;
    };
  }
}

#endif
