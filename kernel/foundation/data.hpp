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
        mesh_halo_map(TopologyType_()),
        functions_on_process_patch(typename MeshType_::attr_base_type_())
      {
      }

      StorageType_<MeshType_, std::allocator<MeshType_> > meshes_on_process_patch;
      StorageType_<HaloType_, std::allocator<HaloType_> > halos_on_process_patch;
      TopologyType_ mesh_halo_map;
      typename MeshType_::attr_base_type_ functions_on_process_patch;
    };

    struct SimpleDataFillPolicy
    {
      template<typename MeshType_, typename HaloType_, typename TopologyType_>
      static void execute(PatchData<MeshType_, HaloType_, TopologyType_> & target, int rank)
      {
        //take exactly 1 XOR 4 processes and initialise the meshes and halos
        if(rank < 4 && rank >= 0)
        {
          MeshType_ local_mesh(rank, &target.functions_on_process_patch); //hide attributes here for now
          target.meshes_on_process_patch.push_back(local_mesh);

          target.functions_on_process_patch.push_back(std::shared_ptr<Foundation::AttributeBase<> >(new Foundation::Attribute<double>()));
          Foundation::MeshAttributeRegistration::execute(target.meshes_on_process_patch.at(0), pl_vertex);

          //add vertices
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_vertex);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_vertex);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_vertex);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_vertex);
          ((Foundation::Attribute<double>*)(target.functions_on_process_patch.at(0).get()))->get_data().push_back(double(rank));
          ((Foundation::Attribute<double>*)(target.functions_on_process_patch.at(0).get()))->get_data().push_back(double(rank));
          ((Foundation::Attribute<double>*)(target.functions_on_process_patch.at(0).get()))->get_data().push_back(double(rank));
          ((Foundation::Attribute<double>*)(target.functions_on_process_patch.at(0).get()))->get_data().push_back(double(rank));

          //add edges
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_edge);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_edge);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_edge);
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_edge);

          //add face
          target.meshes_on_process_patch.at(0).add_polytope(Foundation::pl_face);

          //create adjacencies
          /*  2    3
           *  *-1--*
           *  2    |
           *  |    3
           *  *--0-*
           *  0    1
           */
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_edge, 0, 0);
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_edge, 0, 2);
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_face, 0, 0);

          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_edge, 1, 0);
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_edge, 1, 3);
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_face, 1, 0);

          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_edge, 2, 1);
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_edge, 2, 2);
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_face, 2, 0);

          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_edge, 3, 1);
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_edge, 3, 3);
          target.meshes_on_process_patch.at(0).add_adjacency(pl_vertex, pl_face, 3, 0);

          //create halos, every patch has three logical halos (we use vertex-halos)
          if(rank == 0)
          {
            HaloType_ h1(target.meshes_on_process_patch.at(0), 1); //right
            h1.add_element_pair(3, 2);
            h1.add_element_pair(1, 0);

            HaloType_ h2(target.meshes_on_process_patch.at(0), 2); //down
            h2.add_element_pair(0, 2);
            h2.add_element_pair(1, 3);

            HaloType_ h3(target.meshes_on_process_patch.at(0), 3); //diag
            h3.add_element_pair(1, 2);

            target.halos_on_process_patch.push_back(h1);
            target.halos_on_process_patch.push_back(h2);
            target.halos_on_process_patch.push_back(h3);

            target.mesh_halo_map.push_back();
            target.mesh_halo_map.get_topology().at(0).push_back(0);
            target.mesh_halo_map.get_topology().at(0).push_back(1);
            target.mesh_halo_map.get_topology().at(0).push_back(2);
          }
          else if(rank == 1)
          {
            HaloType_ h1(target.meshes_on_process_patch.at(0), 0); //left
            h1.add_element_pair(0, 1);
            h1.add_element_pair(2, 3);

            HaloType_ h2(target.meshes_on_process_patch.at(0), 3); //down
            h2.add_element_pair(0, 2);
            h2.add_element_pair(1, 3);

            HaloType_ h3(target.meshes_on_process_patch.at(0), 2); //diag
            h3.add_element_pair(0, 3);

            target.halos_on_process_patch.push_back(h1);
            target.halos_on_process_patch.push_back(h2);
            target.halos_on_process_patch.push_back(h3);

            target.mesh_halo_map.push_back();
            target.mesh_halo_map.get_topology().at(0).push_back(0);
            target.mesh_halo_map.get_topology().at(0).push_back(1);
            target.mesh_halo_map.get_topology().at(0).push_back(2);
          }
          else if(rank == 2)
          {
            HaloType_ h1(target.meshes_on_process_patch.at(0), 0); //up
            h1.add_element_pair(2, 0);
            h1.add_element_pair(3, 1);

            HaloType_ h2(target.meshes_on_process_patch.at(0), 3); //right
            h2.add_element_pair(1, 0);
            h2.add_element_pair(3, 2);

            HaloType_ h3(target.meshes_on_process_patch.at(0), 1); //diag
            h3.add_element_pair(3, 0);

            target.halos_on_process_patch.push_back(h1);
            target.halos_on_process_patch.push_back(h2);
            target.halos_on_process_patch.push_back(h3);

            target.mesh_halo_map.push_back();
            target.mesh_halo_map.get_topology().at(0).push_back(0);
            target.mesh_halo_map.get_topology().at(0).push_back(1);
            target.mesh_halo_map.get_topology().at(0).push_back(2);
          }
          else if(rank == 3)
          {
            HaloType_ h1(target.meshes_on_process_patch.at(0), 1); //up
            h1.add_element_pair(2, 0);
            h1.add_element_pair(3, 1);

            HaloType_ h2(target.meshes_on_process_patch.at(0), 2); //left
            h2.add_element_pair(0, 1);
            h2.add_element_pair(2, 3);

            HaloType_ h3(target.meshes_on_process_patch.at(0), 0); //diag
            h3.add_element_pair(2, 1);

            target.halos_on_process_patch.push_back(h1);
            target.halos_on_process_patch.push_back(h2);
            target.halos_on_process_patch.push_back(h3);

            target.mesh_halo_map.push_back();
            target.mesh_halo_map.get_topology().at(0).push_back(0);
            target.mesh_halo_map.get_topology().at(0).push_back(1);
            target.mesh_halo_map.get_topology().at(0).push_back(2);
          }
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
