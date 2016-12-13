#pragma once
#ifndef CONTROL_DOMAIN_DOMAIN_LEVEL_HPP
#define CONTROL_DOMAIN_DOMAIN_LEVEL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/geometry/mesh_node.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Domain
    {
      template<typename Mesh_>
      class DomainLevel
      {
      public:
        typedef Mesh_ MeshType;
        typedef Geometry::MeshPart<Mesh_> PartType;
        typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

      protected:
        int _level_index;
        std::shared_ptr<MeshNodeType> _mesh_node;

      public:
        explicit DomainLevel(int _lvl_idx, std::shared_ptr<MeshNodeType> node) :
          _level_index(_lvl_idx),
          _mesh_node(node)
        {
          XASSERT(node != nullptr);
        }

        DomainLevel(DomainLevel&& other) :
          _level_index(_level_index),
          _mesh_node(other._mesh_node)
        {
        }

        virtual ~DomainLevel()
        {
        }

        std::size_t bytes() const
        {
          return _mesh_node->bytes();
        }

        int get_level_index() const
        {
          return _level_index;
        }

        MeshNodeType* get_mesh_node()
        {
          return _mesh_node.get();
        }

        const MeshNodeType* get_mesh_node() const
        {
          return _mesh_node.get();
        }

        MeshType& get_mesh()
        {
          return *_mesh_node->get_mesh();
        }

        const MeshType& get_mesh() const
        {
          return *_mesh_node->get_mesh();
        }

        const PartType* find_halo_part(int rank) const
        {
          // try to fetch the corresponding mesh part node
          return this->_mesh_node->find_mesh_part(String("_halo:") + stringify(rank));
        }

        const PartType* find_patch_part(int rank) const
        {
          return this->_mesh_node->find_mesh_part(String("_patch:") + stringify(rank));
        }
      }; // class DomainLevel<...>

      template<typename Mesh_, typename Trafo_, typename Space_>
      class SimpleDomainLevel :
        public Domain::DomainLevel<Mesh_>
      {
      public:
        typedef Domain::DomainLevel<Mesh_> BaseClass;

        typedef Mesh_ MeshType;
        typedef Trafo_ TrafoType;
        typedef Space_ SpaceType;

        TrafoType trafo;
        SpaceType space;

      public:
        explicit SimpleDomainLevel(int lvl_idx, std::shared_ptr<Geometry::RootMeshNode<MeshType>> node) :
          BaseClass(lvl_idx, node),
          trafo(BaseClass::get_mesh()),
          space(trafo)
        {
        }
      }; // class SimpleDomainLevel<...>

      template<typename Mesh_, typename Trafo_, typename SpaceVelo_, typename SpacePres_>
      class StokesDomainLevel :
        public Domain::DomainLevel<Mesh_>
      {
      public:
        typedef Domain::DomainLevel<Mesh_> BaseClass;

        typedef Mesh_ MeshType;
        typedef Trafo_ TrafoType;
        typedef SpaceVelo_ SpaceVeloType;
        typedef SpacePres_ SpacePresType;

        TrafoType trafo;
        SpaceVeloType space_velo;
        SpacePresType space_pres;

      public:
        explicit StokesDomainLevel(int lvl_idx, std::shared_ptr<Geometry::RootMeshNode<MeshType>> node) :
          BaseClass(lvl_idx, node),
          trafo(BaseClass::get_mesh()),
          space_velo(trafo),
          space_pres(trafo)
        {
        }
      }; // class StokesDomainLevel<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_DOMAIN_DOMAIN_LEVEL_HPP
