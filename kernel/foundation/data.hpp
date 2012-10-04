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
        mesh_halo_map(StorageType_<TopologyType_, std::allocator<TopologyType_> >()),
        functions_on_process_patch(typename MeshType_::attr_base_type_())
      {
      }

      StorageType_<MeshType_, std::allocator<MeshType_> > meshes_on_process_patch;
      StorageType_<HaloType_, std::allocator<HaloType_> > halos_on_process_patch;
      StorageType_<TopologyType_, std::allocator<TopologyType_> > mesh_halo_map;
      typename MeshType_::attr_base_type_ functions_on_process_patch;
    };

    struct SimpleDataFillPolicy
    {
      template<typename MeshType_, typename HaloType_, typename TopologyType_>
      static void execute(PatchData<MeshType_, HaloType_, TopologyType_> & target, int rank)
      {
        //take exactly 1 XOR 4 processes and initialise the meshes and halos
        if(rank == 0)
        {
          MeshType_ local_mesh(0, &target.functions_on_process_patch); //hide attributes here for now
          target.meshes_on_process_patch.push_back(local_mesh);

          target.functions_on_process_patch.push_back(std::shared_ptr<Foundation::AttributeBase<> >(new Foundation::Attribute<double>()));
          Foundation::MeshAttributeRegistration::execute(target.meshes_on_process_patch.at(0), pl_vertex);
          //add vertices
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_vertex);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_vertex);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_vertex);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_vertex);
          ((Foundation::Attribute<double>*)(target.functions_on_process_patch.at(0).get()))->get_data().push_back(double(0.));
        }
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
