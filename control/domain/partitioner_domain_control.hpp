#pragma once
#ifndef CONTROL_DOMAIN_PARTITIONER_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_PARTITIONER_DOMAIN_CONTROL_HPP 1

#include <kernel/foundation/comm_base.hpp>
#include <kernel/foundation/pexecutor.hpp>
#include <kernel/foundation/pgraph.hpp>
#include <kernel/foundation/psynch.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_node.hpp>

#include <control/domain/domain_control.hpp>

namespace FEAST
{
  namespace Control
  {
    namespace Domain
    {
      /**
       * \brief Domain control for a distributed mesh
       *
       * Distributed means that each process has its own patch root mesh node.
       *
       * \author Markus Geveler
       */
      template<typename PartitionerExecutorT_, typename MeshType_>
      class PartitionerDomainControl : public Control::Domain::DomainControl<MeshType_>
      {
        public:
          typedef PartitionerExecutorT_ PartT;
          /// Our base class
          typedef Control::Domain::DomainControl<MeshType_> BaseClass;
          /// The the root meshnode will have
          typedef Geometry::RootMeshNode<MeshType_> MeshNodeType;

          /**
           * \brief Constructor that reads the base mesh from a file
           */
          explicit PartitionerDomainControl(int lvl_max, int lvl_min, Index min_elems_partitioner,
          const String& meshfile)
          {
            std::vector<Index> ranks, ctags;
            std::stringstream synchstream;

            //MASTER
            if(Foundation::Comm::rank() == 0)
            {
              std::ifstream ifs(meshfile.c_str(), std::ios::binary);
              if(!ifs.is_open())
              {
                throw FileNotFound(meshfile);
              }
              synchstream << ifs.rdbuf();
            }

            // Synchronise the mesh stream in parallel mode
#ifndef SERIAL
            //MASTER to ALL
            Foundation::PSynch<PartT>::exec(synchstream);
#endif

            // Create MeshFileReader using the synched stream
            Geometry::MeshFileReader mesh_file_reader(synchstream);

            _create(lvl_max, lvl_min, min_elems_partitioner, mesh_file_reader);

          }

          /**
           * \brief Constructor that parses the base mesh from a MeshFileReader
           *
           * \note Since we are talking parallel here, the std::istream the MeshFileReader uses has to have been
           * synchronised before it was passes to it.
           */
          explicit PartitionerDomainControl(int lvl_max, int lvl_min, Index min_elems_partitioner,
          Geometry::MeshFileReader& mesh_file_reader)
          {
            _create(lvl_max, lvl_min, min_elems_partitioner, mesh_file_reader);
          }

        private:
          /**
           * \brief Routine handling the actual creation of the object
           *
           * This is outsourced from the constructor because the mesh_file_reader might be created and synchronised
           * there.
           */
          void _create(int lvl_max, int lvl_min, const Index min_elems_partitioner,
          Geometry::MeshFileReader& mesh_file_reader)
          {
            // create root mesh node
            MeshNodeType* mesh_node = nullptr;
            std::vector<Index> ranks, ctags;
            {
              //ALL
              Geometry::MeshAtlas<MeshType_>* atlas(new Geometry::MeshAtlas<MeshType_>());
              Geometry::RootMeshNode<MeshType_>* base_mesh_node(new Geometry::RootMeshNode<MeshType_>(nullptr, atlas));

              // Try to parse the mesh file
#ifndef DEBUG
              try
#endif
              {
                mesh_file_reader.parse(*base_mesh_node, *atlas);
              }
#ifndef DEBUG
              catch(std::exception& exc)
              {
                std::cerr << "ERROR: " << exc.what() << std::endl;
                if(Foundation::Comm::rank() == 0)
                  std::cout << "FAILED: " << exc.what() << " when parsing." << std::endl;
              }
              catch(...)
              {
                std::cerr << "ERROR: unknown exception" << std::endl;
                if(Foundation::Comm::rank() == 0)
                  std::cout << "FAILED: unknown exception" << std::endl;
              }
#endif
              //base-mesh node creation

              //partitioner code
              MeshType_* root_mesh(nullptr);
              root_mesh = base_mesh_node->get_mesh();

              Index num_elements(root_mesh->get_num_entities(MeshType_::shape_dim));
              while(num_elements < min_elems_partitioner)
              {
                MeshNodeType* coarse_node = base_mesh_node;
                base_mesh_node = coarse_node->refine();
                delete coarse_node;

                num_elements = base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim);
              }
              while(num_elements < Foundation::Comm::size())
              {
                MeshNodeType* coarse_node = base_mesh_node;
                base_mesh_node = coarse_node->refine();
                delete coarse_node;

                num_elements = base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim);
              }

              root_mesh = base_mesh_node->get_mesh();

              //ALL
#ifndef SERIAL
              Index num_global_elements(root_mesh->get_num_entities(MeshType_::shape_dim));
              typename PartT::PGraphT global_dual(
                *root_mesh, num_global_elements, Foundation::Communicator(MPI_COMM_WORLD));

              /*for(Index i(0) ; i < Index(global_dual.get_num_vtx()) ; ++i)
                {
                std::cout << Comm::rank() << " global dual " << i << "->";
                for(Index j(Index(global_dual.get_xadj()[i])) ; j < Index(global_dual.get_xadj()[i+1]) ; ++j)
                std::cout << " " << global_dual.get_adjncy()[j];

                std::cout << std::endl;
                }*/

              /*std::cout << Comm::rank() << ": APP: vtxdist global " << std::endl;
                for(Index i(0) ; i < Comm::size() + 1 ; ++i)
                std::cout << Comm::rank() << ": APP: vtxdist global[" << i << "] " << global_dual.get_vtxdist()[i] << std::endl;*/

              //local input for k-way partitioning
              auto local_dual(global_dual.create_local());
              /*std::cout << Comm::rank() << " local dual xadj: " << std::endl;
                for(Index i(0) ; i < Index(local_dual->get_num_vtx() + 1); ++i)
                {
                std::cout << " " << local_dual->get_xadj()[i];
                if(i == Index(local_dual->get_num_vtx()))
                std::cout << std::endl;
                }*/

              /*for(Index i(0) ; i < Comm::size() + 1 ; ++i)
                std::cout << Comm::rank() << ": APP: vtxdist local[" << i << "] " << local_dual->get_vtxdist()[i] << std::endl;*/

              auto part(PartT::part(*((typename PartT::PGraphT*)local_dual.get())));

              auto synched_part(Foundation::PSynch<PartT>::exec(
                part, typename PartT::IndexType(base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim))));

              PartT::fill_comm_structs_global(synched_part, global_dual);

              /*for(Index i(0) ; i <  base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim); ++i)
                std::cout << Comm::rank() << ": APP: PART[" << i << "] " << synched_part.get()[i] << std::endl;

                for(auto cr_i : synched_part.get_comm_ranks())
                std::cout << Comm::rank() << ": APP: COMM_RANKS " << cr_i << std::endl;

                for(auto ct_i : synched_part.get_comm_tags())
                std::cout << Comm::rank() << ": APP: COMM_TAGS " << ct_i << std::endl;*/

              Adjacency::Graph ranks_at_elem(synched_part.rank_at_element());
#else
              Index* ptr_serial(new Index[base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim) + 1]);
              Index* part_serial(new Index[base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim)]);
              for(Index i(0) ; i < Index(base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim) + 1) ; ++i)
                ptr_serial[i] = i;
              for(Index i(0) ; i < Index(base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim)) ; ++i)
                part_serial[i] = 0;

              Adjacency::Graph ranks_at_elem(
                base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim),
                1,
                base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim),
                ptr_serial, nullptr, part_serial);
#endif
              /*for(Index i(0) ; i <  base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim) + 1; ++i)
                std::cout << Comm::rank() << ": APP: ranks_at_elem ptr[" << i << "] " << (ranks_at_elem.get_domain_ptr())[i] << std::endl;
                for(Index i(0) ; i <  base_mesh_node->get_mesh()->get_num_entities(MeshType_::shape_dim) ; ++i)
                std::cout << Comm::rank() << ": APP: ranks_at_elem idx[" << i << "] " << (ranks_at_elem.get_image_idx())[i] << std::endl;*/

#ifndef SERIAL
              ranks = synched_part.get_comm_ranks();
              ctags = synched_part.get_comm_tags();
#endif
              /*for(auto& ranks_i : ranks)
                std::cout << Comm::rank() << ": APP: COMM_ranks " << ranks_i << std::endl;

                for(auto& tags_i : ctags)
                std::cout << Comm::rank() << ": APP: COMM_tags " << tags_i << std::endl;*/

              // <<<<< partitioner code

              // create patch mesh node
#ifndef SERIAL
              mesh_node = base_mesh_node->extract_patch(
                Index(Foundation::Comm::rank(synched_part.get_comm())), ranks_at_elem, ranks);
#else
              mesh_node = base_mesh_node->extract_patch(Index(0), ranks_at_elem, ranks);
#endif

              // delete base-mesh node
              delete base_mesh_node;
            }

            // push layer
            this->_layers.push_back(new typename BaseClass::LayerType(std::move(ranks), std::move(ctags)));

            // adjust lvl_max and lvl_min
            lvl_max = Math::max(lvl_max, 0);
            if(lvl_min < 0)
              lvl_min = Math::max(lvl_max + lvl_min + 1, 0);
            else
              lvl_min = Math::min(lvl_min, lvl_max);

            int lvl = 2;
            // refine up to desired minimum level
            for(; lvl < lvl_min; ++lvl)
            {
              MeshNodeType* coarse_node = mesh_node;
              mesh_node = coarse_node->refine();
              delete coarse_node;
            }

            // add coarse mesh node
            this->_levels.push_back(new typename BaseClass::LevelType(lvl, mesh_node));

            // refine up to desired maximum level
            for(; lvl < lvl_max;)
            {
              MeshNodeType* coarse_node = mesh_node;
              mesh_node = coarse_node->refine();
              this->_levels.push_back(new typename BaseClass::LevelType(++lvl, mesh_node));
            }
          }
      };

    } // namespace Domain
  } // namespace Control
} // namespace FEAST

#endif // CONTROL_DOMAIN_PARTITIONER_DOMAIN_CONTROL_HPP
