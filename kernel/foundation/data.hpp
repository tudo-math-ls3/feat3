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
        submesh(),
        comm_halos(),
        boundaries(),
        attrs()
      {
      }

      PData(const PData& other) :
        submesh(other.submesh),
        comm_halos(other.comm_halos),
        boundaries(other.boundaries),
        attrs(other.attrs)
      {
      }

      PData& operator=(const PData& other)
      {
        this->submesh = other.submesh;
        this->comm_halos = other.comm_halos;
        this->boundaries = other.boundaries;
        this->attrs = other.attrs;
      }

      std::shared_ptr<SubMesh<Dim_, t_, os_> > submesh;
      os_<std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> > > > comm_halos;
      os_<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_>, std::allocator<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_> > > boundaries;
      os_<std::shared_ptr<AttributeBase<os_> >, std::allocator<std::shared_ptr<AttributeBase<os_> > > > attrs;
    };
  }
}
#endif
