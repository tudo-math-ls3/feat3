#pragma once
#ifndef KERNEL_FOUNDATION_DATA_HH
#define KERNEL_FOUNDATION_DATA_HH 1

#include<kernel/base_header.hpp>
#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/halo.hpp>

namespace FEAST
{
  namespace Foundation
  {
    template<typename MeshType_, typename HaloType_, typename TopologyType_, template<typename, typename> class StorageType_ = std::vector>
    struct PatchData
    {
      PatchData() :
        meshes_on_process_patch(StorageType_<MeshType_, std::allocator<MeshType_> >()),
        halos_on_process_patch(StorageType_<HaloType_, std::allocator<HaloType_> >()),
        mesh_halo_map(StorageType_<TopologyType_, std::allocator<TopologyType_> >())
      {
      }

      StorageType_<MeshType_, std::allocator<MeshType_> > meshes_on_process_patch;
      StorageType_<HaloType_, std::allocator<HaloType_> > halos_on_process_patch;
      StorageType_<TopologyType_, std::allocator<TopologyType_> > mesh_halo_map;
    };

    struct SimpleDataFillPolicy
    {
      template<typename MeshType_, typename HaloType_, typename TopologyType_>
      static void execute(PatchData<MeshType_, HaloType_, TopologyType_> & target, int rank)
      {
      }
    };

    struct FillDataFromFilesPolicy
    {
      ///TODO
      template<typename MeshType_, typename HaloType_, typename TopologyType_>
      static void execute(PatchData<MeshType_, HaloType_, TopologyType_> & target, int rank, const std::string filename)
      {
      }
    };
  }
}
#endif
