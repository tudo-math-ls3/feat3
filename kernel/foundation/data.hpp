#pragma once
#ifndef KERNEL_FOUNDATION_DATA_HPP
#define KERNEL_FOUNDATION_DATA_HPP 1

#include<kernel/base_header.hpp>
#include <kernel/foundation/attribute.hpp>
#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/sub_mesh.hpp>
#include <kernel/foundation/halo.hpp>

namespace FEAST
{
  namespace Foundation
  {

    template<
      typename Dim_,
      typename t_,
      template <typename, typename> class os_,
      template <typename, typename, template<typename, typename> class > class MeshType_,
      typename DT_>
    struct PData
    {
      PData() :
        basemesh(),
        submesh(),
        comm_halos(),
        basemesh_boundaries(),
        boundaries(),
        attrs(),
        comm_cost(0),
        comp_cost(0)
      {
      }

      PData(const PData& other) :
        basemesh(other.basemesh),
        submesh(other.submesh),
        comm_halos(other.comm_halos),
        basemesh_boundaries(other.basemesh_boundaries),
        boundaries(other.boundaries),
        attrs(other.attrs),
        comm_cost(other.comm_cost),
        comp_cost(other.comp_cost)
      {
      }

      PData& operator=(const PData& other)
      {
        this->basemesh = other.basemesh;
        this->submesh = other.submesh;
        this->comm_halos = other.comm_halos;
        this->basemesh_boundaries = other.basemesh_boundaries;
        this->boundaries = other.boundaries;
        this->attrs = other.attrs;
        this->comm_cost = other.comm_cost;
        this->comp_cost = other.comp_cost;
      }

      Mesh<Dim_, t_, os_> basemesh;
      std::shared_ptr<SubMesh<Dim_, t_, os_> > submesh;
      os_<std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> > > > comm_halos;
      os_<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_>, std::allocator<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_> > > basemesh_boundaries;
      os_<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_>, std::allocator<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_> > > boundaries;
      os_<std::shared_ptr<AttributeBase<os_> >, std::allocator<std::shared_ptr<AttributeBase<os_> > > > attrs;
      DT_ comm_cost;
      DT_ comp_cost;
    };
  }
}
#endif
