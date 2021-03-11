// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_UNIT_CUBE_PATCH_GENERATOR_HPP
#define KERNEL_GEOMETRY_UNIT_CUBE_PATCH_GENERATOR_HPP 1

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>

#include <vector>
#include <stdexcept>

namespace FEAT
{
  namespace Geometry
  {
    /// \todo Documentation
    template<typename Mesh_>
    class UnitCubePatchGenerator;

    /// \cond internal
    template<typename Coord_>
    class UnitCubePatchGenerator<ConformalMesh<Shape::Hypercube<1>, 1, Coord_>>
    {
    public:
      typedef Shape::Hypercube<1> ShapeType;
      typedef Coord_ CoordType;
      typedef ConformalMesh<ShapeType, 1, Coord_> MeshType;
      typedef MeshPart<MeshType> PartType;
      typedef RootMeshNode<MeshType> MeshNodeType;

    private:
      class RootMeshFactory :
        public Geometry::UnitCubeFactory<MeshType>
      {
      public:
        typedef typename MeshType::VertexSetType VertexSetType;

      private:
        const int _i, _n;

      public:
        RootMeshFactory(int i, int n) :
          _i(i), _n(n)
        {
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          vertex_set[0][0] = CoordType(_i  ) / CoordType(_n);
          vertex_set[1][0] = CoordType(_i+1) / CoordType(_n);
        }
      };

      class Halo0Factory :
        public Geometry::Factory<PartType>
      {
      public:
        typedef Geometry::Factory<PartType> BaseClass;
        typedef typename BaseClass::AttributeSetContainer AttributeSetContainer;
        typedef typename BaseClass::IndexSetHolderType IndexSetHolderType;
        typedef typename BaseClass::TargetSetHolderType TargetSetHolderType;

      private:
        const int _k;

      public:
        explicit Halo0Factory(int k) :
          _k(k)
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          return Index(dim == 0 ? 1 : 0);
        }

        virtual void fill_attribute_sets(AttributeSetContainer&) override
        {
        }

        virtual void fill_index_sets(IndexSetHolderType*&) override
        {
        }

        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
        {
          target_set_holder.template get_target_set<0>()[0] = Index(_k);
        }

      };

      static PartType* create_halo0(int k)
      {
        Halo0Factory factory(k);
        return new PartType(factory);
      }

    public:
      static Index create(int rank, int nprocs, std::shared_ptr<MeshNodeType>& node, std::vector<int>& ranks)
      {
        MeshNodeType* nnode = nullptr;
        std::vector<Index> cranks, ctags;
        int lvl = create(rank, nprocs, nnode, cranks, ctags);
        node = std::shared_ptr<MeshNodeType>(nnode);
        for(auto it = cranks.begin(); it != cranks.end(); ++it)
          ranks.push_back(int(*it));
        return Index(lvl);
      }

      static int create(int rank, int nprocs, MeshNodeType*& node, std::vector<Index>& ranks, std::vector<Index>& ctags)
      {
        XASSERT(nprocs > 0);
        XASSERT((rank >= 0) && (rank < nprocs));

        // determine slice count and refinement level
        int level(0);
        int n(1);
        for(; n < nprocs; n *= 2, ++level) {}
        XASSERTM(n == nprocs, "number of processes must be a power of 2");

        // decompose rank to (i)
        const int ii(rank);

        // create root mesh node
        {
          RootMeshFactory factory(ii, n);
          node = new MeshNodeType(new MeshType(factory));
        }

        // left neighbor
        if(ii > 0)
        {
          ranks.push_back(Index(ii-1));
          ctags.push_back(Index(ii-1));
          node->add_halo(int(ranks.back()), create_halo0(0));
          node->add_mesh_part("bnd:0", nullptr);
        }
        else
        {
          node->add_mesh_part("bnd:0", create_halo0(0));
        }
        // right neighbor
        if(ii+1 < n)
        {
          ranks.push_back(Index(ii+1));
          ctags.push_back(Index(ii));
          node->add_halo(int(ranks.back()), create_halo0(1));
          node->add_mesh_part("bnd:1", nullptr);
        }
        else
        {
          node->add_mesh_part("bnd:1", create_halo0(1));
        }

        // return refinement level
        return level;
      }
    }; // UnitCubePatchGenerator<Hypercube<1>>

    template<typename Coord_>
    class UnitCubePatchGenerator<ConformalMesh<Shape::Hypercube<2>, 2, Coord_>>
    {
    public:
      typedef Shape::Hypercube<2> ShapeType;
      typedef Coord_ CoordType;
      typedef ConformalMesh<ShapeType, 2, Coord_> MeshType;
      typedef MeshPart<MeshType> PartType;
      typedef RootMeshNode<MeshType> MeshNodeType;

    private:
      class RootMeshFactory :
        public Geometry::UnitCubeFactory<MeshType>
      {
      public:
        typedef typename MeshType::VertexSetType VertexSetType;

      private:
        const int _i, _j, _n;

      public:
        RootMeshFactory(int i, int j, int n) :
          _i(i), _j(j), _n(n)
        {
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          vertex_set[0][0] = CoordType(_i  ) / CoordType(_n);
          vertex_set[0][1] = CoordType(_j  ) / CoordType(_n);
          vertex_set[1][0] = CoordType(_i+1) / CoordType(_n);
          vertex_set[1][1] = CoordType(_j  ) / CoordType(_n);
          vertex_set[2][0] = CoordType(_i  ) / CoordType(_n);
          vertex_set[2][1] = CoordType(_j+1) / CoordType(_n);
          vertex_set[3][0] = CoordType(_i+1) / CoordType(_n);
          vertex_set[3][1] = CoordType(_j+1) / CoordType(_n);
        }
      };

      class Halo0Factory :
        public Geometry::Factory<PartType>
      {
      public:
        typedef Geometry::Factory<PartType> BaseClass;
        typedef typename BaseClass::AttributeSetContainer AttributeSetContainer;
        typedef typename BaseClass::IndexSetHolderType IndexSetHolderType;
        typedef typename BaseClass::TargetSetHolderType TargetSetHolderType;

      private:
        const int _k;

      public:
        explicit Halo0Factory(int k) :
          _k(k)
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          return Index(dim == 0 ? 1 : 0);
        }

        virtual void fill_attribute_sets(AttributeSetContainer&) override
        {
        }

        virtual void fill_index_sets(IndexSetHolderType*&) override
        {
        }

        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
        {
          target_set_holder.template get_target_set<0>()[0] = Index(_k);
        }
      };

      class Halo1Factory :
        public Geometry::Factory<PartType>
      {
      public:
        typedef Geometry::Factory<PartType> BaseClass;
        typedef typename BaseClass::AttributeSetContainer AttributeSetContainer;
        typedef typename BaseClass::IndexSetHolderType IndexSetHolderType;
        typedef typename BaseClass::TargetSetHolderType TargetSetHolderType;

      private:
        const int _k;

      public:
        explicit Halo1Factory(int k) :
          _k(k)
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          return Index(dim == 0 ? 2 : (dim == 1 ? 1 : 0));
        }

        virtual void fill_attribute_sets(AttributeSetContainer&) override
        {
        }

        virtual void fill_index_sets(IndexSetHolderType*&) override
        {
        }

        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
        {
          typedef Geometry::Intern::FaceIndexMapping<Shape::Quadrilateral, 1, 0> Fim;
          target_set_holder.template get_target_set<0>()[0] = Index(Fim::map(_k, 0));
          target_set_holder.template get_target_set<0>()[1] = Index(Fim::map(_k, 1));
          target_set_holder.template get_target_set<1>()[0] = Index(_k);
        }
      };

      static PartType* create_halo0(int k)
      {
        Halo0Factory factory(k);
        return new PartType(factory);
      }

      static PartType* create_halo1(int k)
      {
        Halo1Factory factory(k);
        return new PartType(factory);
      }

    public:
      static Index create(int rank, int nprocs, std::shared_ptr<MeshNodeType>& node, std::vector<int>& ranks)
      {
        MeshNodeType* nnode = nullptr;
        std::vector<Index> cranks, ctags;
        int lvl = create(rank, nprocs, nnode, cranks, ctags);
        node = std::shared_ptr<MeshNodeType>(nnode);
        for(auto it = cranks.begin(); it != cranks.end(); ++it)
          ranks.push_back(int(*it));
        return Index(lvl);
      }

      static int create(int rank, int nprocs, MeshNodeType*& node, std::vector<Index>& ranks, std::vector<Index>& ctags)
      {
        XASSERT(nprocs > 0);
        XASSERT((rank >= 0) && (rank < nprocs));

        // determine slice count and refinement level
        int level(0);
        int n(1);
        for(; n*n < nprocs; n *= 2, ++level) {}
        XASSERTM(n*n == nprocs, "number of processes must be a power of 4");

        // decompose rank to (i,j)
        const int ii(rank % n);
        const int jj(rank / n);

        // comm tag offsets
        const Index cto_x0 = Index(0);
        const Index cto_x1 = cto_x0 + Index((n-1)*(n-1));
        const Index cto_h  = cto_x1 + Index((n-1)*(n-1));
        const Index cto_v  = cto_h  + Index(n*(n-1));

        // create root mesh node
        {
          RootMeshFactory factory(ii, jj, n);
          node = new MeshNodeType(new MeshType(factory));
        }

        // lower left neighbor
        if((ii > 0) && (jj > 0))
        {
          ranks.push_back(Index(n*(jj-1) + ii-1));
          ctags.push_back(cto_x0 + Index((n-1)*(jj-1) + ii-1));
          node->add_halo(int(ranks.back()), create_halo0(0));
        }
        // lower right neighbor
        if((ii+1 < n) && (jj > 0))
        {
          ranks.push_back(Index(n*(jj-1) + ii+1));
          ctags.push_back(cto_x1 + Index((n-1)*(jj-1) + ii));
          node->add_halo(int(ranks.back()), create_halo0(1));
        }
        // upper left neighbor
        if((ii > 0) && (jj+1 < n))
        {
          ranks.push_back(Index(n*(jj+1) + ii-1));
          ctags.push_back(cto_x1 + Index((n-1)*jj + ii-1));
          node->add_halo(int(ranks.back()), create_halo0(2));
        }
        // upper right neighbor
        if((ii+1 < n) && (jj+1 < n))
        {
          ranks.push_back(Index(n*(jj+1) + ii+1));
          ctags.push_back(cto_x0 + Index((n-1)*jj + ii));
          node->add_halo(int(ranks.back()), create_halo0(3));
        }

        // bottom neighbor
        if(jj > 0)
        {
          ranks.push_back(Index(n*(jj-1) + ii));
          ctags.push_back(cto_v + Index(n*(jj-1) + ii));
          node->add_halo(int(ranks.back()), create_halo1(0));
          node->add_mesh_part("bnd:0", nullptr);
        }
        else
        {
          node->add_mesh_part("bnd:0", create_halo1(0));
        }
        // top neighbor
        if(jj+1 < n)
        {
          ranks.push_back(Index(n*(jj+1) + ii));
          ctags.push_back(cto_v + Index(n*jj + ii));
          node->add_halo(int(ranks.back()), create_halo1(1));
          node->add_mesh_part("bnd:1", nullptr);
        }
        else
        {
          node->add_mesh_part("bnd:1", create_halo1(1));
        }
        // left neighbor
        if(ii > 0)
        {
          ranks.push_back(Index(n*(jj) + ii-1));
          ctags.push_back(cto_h + Index(n*jj + ii - 1));
          node->add_halo(int(ranks.back()), create_halo1(2));
          node->add_mesh_part("bnd:2", nullptr);
        }
        else
        {
          node->add_mesh_part("bnd:2", create_halo1(2));
        }
        // right neighbor
        if(ii+1 < n)
        {
          ranks.push_back(Index(n*(jj) + ii+1));
          ctags.push_back(cto_h + Index(n*jj + ii));
          node->add_halo(int(ranks.back()), create_halo1(3));
          node->add_mesh_part("bnd:3", nullptr);
        }
        else
        {
          node->add_mesh_part("bnd:3", create_halo1(3));
        }

        // return refinement level
        return level;
      }
    }; // UnitCubePatchGenerator<Hypercube<2>>


    template<typename Coord_>
    class UnitCubePatchGenerator<ConformalMesh<Shape::Hypercube<3>, 3, Coord_>>
    {
    public:
      typedef Shape::Hypercube<3> ShapeType;
      typedef Coord_ CoordType;
      typedef ConformalMesh<ShapeType, 3, Coord_> MeshType;
      typedef MeshPart<MeshType> PartType;
      typedef RootMeshNode<MeshType> MeshNodeType;

    private:
      class RootMeshFactory :
        public Geometry::UnitCubeFactory<MeshType>
      {
      public:
        typedef typename MeshType::VertexSetType VertexSetType;

      private:
        const int _i, _j, _k, _n;

      public:
        RootMeshFactory(int i, int j, int k, int n) :
          _i(i), _j(j), _k(k), _n(n)
        {
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          vertex_set[0][0] = CoordType(_i  ) / CoordType(_n);
          vertex_set[0][1] = CoordType(_j  ) / CoordType(_n);
          vertex_set[0][2] = CoordType(_k  ) / CoordType(_n);

          vertex_set[1][0] = CoordType(_i+1) / CoordType(_n);
          vertex_set[1][1] = CoordType(_j  ) / CoordType(_n);
          vertex_set[1][2] = CoordType(_k  ) / CoordType(_n);

          vertex_set[2][0] = CoordType(_i  ) / CoordType(_n);
          vertex_set[2][1] = CoordType(_j+1) / CoordType(_n);
          vertex_set[2][2] = CoordType(_k  ) / CoordType(_n);

          vertex_set[3][0] = CoordType(_i+1) / CoordType(_n);
          vertex_set[3][1] = CoordType(_j+1) / CoordType(_n);
          vertex_set[3][2] = CoordType(_k  ) / CoordType(_n);

          vertex_set[4][0] = CoordType(_i  ) / CoordType(_n);
          vertex_set[4][1] = CoordType(_j  ) / CoordType(_n);
          vertex_set[4][2] = CoordType(_k+1) / CoordType(_n);

          vertex_set[5][0] = CoordType(_i+1) / CoordType(_n);
          vertex_set[5][1] = CoordType(_j  ) / CoordType(_n);
          vertex_set[5][2] = CoordType(_k+1) / CoordType(_n);

          vertex_set[6][0] = CoordType(_i  ) / CoordType(_n);
          vertex_set[6][1] = CoordType(_j+1) / CoordType(_n);
          vertex_set[6][2] = CoordType(_k+1) / CoordType(_n);

          vertex_set[7][0] = CoordType(_i+1) / CoordType(_n);
          vertex_set[7][1] = CoordType(_j+1) / CoordType(_n);
          vertex_set[7][2] = CoordType(_k+1) / CoordType(_n);
        }
      };

      class Halo0Factory :
        public Geometry::Factory<PartType>
      {
      public:
        typedef Geometry::Factory<PartType> BaseClass;
        typedef typename BaseClass::AttributeSetContainer AttributeSetContainer;
        typedef typename BaseClass::IndexSetHolderType IndexSetHolderType;
        typedef typename BaseClass::TargetSetHolderType TargetSetHolderType;

      private:
        const int _p;

      public:
        explicit Halo0Factory(int p) :
          _p(p)
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          return Index(dim == 0 ? 1 : 0);
        }

        virtual void fill_attribute_sets(AttributeSetContainer&) override
        {
        }

        virtual void fill_index_sets(IndexSetHolderType*&) override
        {
        }

        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
        {
          target_set_holder.template get_target_set<0>()[0] = Index(_p);
        }
      };

      class Halo1Factory :
        public Geometry::Factory<PartType>
      {
      public:
        typedef Geometry::Factory<PartType> BaseClass;
        typedef typename BaseClass::AttributeSetContainer AttributeSetContainer;
        typedef typename BaseClass::IndexSetHolderType IndexSetHolderType;
        typedef typename BaseClass::TargetSetHolderType TargetSetHolderType;

      private:
        const int _p;

      public:
        explicit Halo1Factory(int p) :
          _p(p)
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          return Index(dim == 0 ? 2 : (dim == 1 ? 1 : 0));
        }

        virtual void fill_attribute_sets(AttributeSetContainer&) override
        {
        }

        virtual void fill_index_sets(IndexSetHolderType*&) override
        {
        }

        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
        {
          typedef Geometry::Intern::FaceIndexMapping<ShapeType, 1, 0> Fim;
          target_set_holder.template get_target_set<0>()[0] = Index(Fim::map(_p, 0));
          target_set_holder.template get_target_set<0>()[1] = Index(Fim::map(_p, 1));
          target_set_holder.template get_target_set<1>()[0] = Index(_p);
        }
      };

      class Halo2Factory :
        public Geometry::Factory<PartType>
      {
      public:
        typedef Geometry::Factory<PartType> BaseClass;
        typedef typename BaseClass::AttributeSetContainer AttributeSetContainer;
        typedef typename BaseClass::IndexSetHolderType IndexSetHolderType;
        typedef typename BaseClass::TargetSetHolderType TargetSetHolderType;

      private:
        const int _p;

      public:
        explicit Halo2Factory(int p) :
          _p(p)
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          return Index(dim == 0 ? 4 : (dim == 1 ? 4 : (dim == 2 ? 1 : 0)));
        }

        virtual void fill_attribute_sets(AttributeSetContainer&) override
        {
        }

        virtual void fill_index_sets(IndexSetHolderType*&) override
        {
        }

        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
        {
          typedef Geometry::Intern::FaceIndexMapping<ShapeType, 2, 0> Fim0;
          target_set_holder.template get_target_set<0>()[0] = Index(Fim0::map(_p, 0));
          target_set_holder.template get_target_set<0>()[1] = Index(Fim0::map(_p, 1));
          target_set_holder.template get_target_set<0>()[2] = Index(Fim0::map(_p, 2));
          target_set_holder.template get_target_set<0>()[3] = Index(Fim0::map(_p, 3));
          typedef Geometry::Intern::FaceIndexMapping<ShapeType, 2, 1> Fim1;
          target_set_holder.template get_target_set<1>()[0] = Index(Fim1::map(_p, 0));
          target_set_holder.template get_target_set<1>()[1] = Index(Fim1::map(_p, 1));
          target_set_holder.template get_target_set<1>()[2] = Index(Fim1::map(_p, 2));
          target_set_holder.template get_target_set<1>()[3] = Index(Fim1::map(_p, 3));
          target_set_holder.template get_target_set<2>()[0] = Index(_p);
        }
      };

      static PartType* create_halo0(int k)
      {
        Halo0Factory factory(k);
        return new PartType(factory);
      }

      static PartType* create_halo1(int k)
      {
        Halo1Factory factory(k);
        return new PartType(factory);
      }

      static PartType* create_halo2(int k)
      {
        Halo2Factory factory(k);
        return new PartType(factory);
      }

    public:
      static Index create(int rank, int nprocs, std::shared_ptr<MeshNodeType>& node, std::vector<int>& ranks)
      {
        MeshNodeType* nnode = nullptr;
        std::vector<Index> cranks, ctags;
        int lvl = create(rank, nprocs, nnode, cranks, ctags);
        node = std::shared_ptr<MeshNodeType>(nnode);
        for(auto it = cranks.begin(); it != cranks.end(); ++it)
          ranks.push_back(int(*it));
        return Index(lvl);
      }

      static int create(int rank, int nprocs, MeshNodeType*& node, std::vector<Index>& ranks, std::vector<Index>& ctags)
      {
        XASSERT(nprocs > 0);
        XASSERT((rank >= 0) && (rank < nprocs));

        // determine slice count and refinement level
        int level(0);
        int n(1);
        for(; n*n*n < nprocs; n *= 2, ++level) {}
        XASSERTM(n*n*n == nprocs, "number of processes must be a power of 8");

        // decompose rank to (i,j,k)
        const int ii(rank % n);
        const int jj((rank / n) % n);
        const int kk(rank / (n*n));

        // comm tag offsets
        const Index cto_x0  = Index(0);
        const Index cto_x1  = cto_x0 + Index((n-1)*(n-1)*(n-1));
        const Index cto_x2  = cto_x1 + Index((n-1)*(n-1)*(n-1));
        const Index cto_x3  = cto_x2 + Index((n-1)*(n-1)*(n-1));

        const Index cto_h0  = cto_x3 + Index((n-1)*(n-1)*(n-1));
        const Index cto_h1  = cto_h0 + Index(n*(n-1)*(n-1));

        const Index cto_v0  = cto_h1 + Index(n*(n-1)*(n-1));
        const Index cto_v1  = cto_v0 + Index(n*(n-1)*(n-1));

        const Index cto_z0  = cto_v1 + Index(n*(n-1)*(n-1));
        const Index cto_z1  = cto_z0 + Index(n*(n-1)*(n-1));

        const Index cto_f0  = cto_z1 + Index(n*(n-1)*(n-1));
        const Index cto_f1  = cto_f0 + Index((n-1)*n*n);
        const Index cto_f2  = cto_f1 + Index((n-1)*n*n);

        // create root mesh
        {
          RootMeshFactory factory(ii, jj, kk, n);
          node = new MeshNodeType(new MeshType(factory));
        }

        // lower/upper -> kk
        // south/north -> jj
        // west /east  -> ii
        // lower south west neighbor
        if((ii > 0) && (jj > 0) && (kk > 0))
        {
          ranks.push_back(Index(n*n*(kk-1) + n*(jj-1) + ii-1));
          ctags.push_back(cto_x0 + Index((n-1)*(n-1)*(kk-1) + (n-1)*(jj-1) + ii-1));
          node->add_halo(int(ranks.back()), create_halo0(0));
        }
        // lower south east neighbor
        if((ii+1 < n) && (jj > 0) && (kk > 0))
        {
          ranks.push_back(Index(n*n*(kk-1) + n*(jj-1) + ii+1));
          ctags.push_back(cto_x1 + Index((n-1)*(n-1)*(kk-1) + (n-1)*(jj-1) + ii));
          node->add_halo(int(ranks.back()), create_halo0(1));
        }
        // lower north west neighbor
        if((ii > 0) && (jj+1 < n) && (kk > 0))
        {
          ranks.push_back(Index(n*n*(kk-1) + n*(jj+1) + ii-1));
          ctags.push_back(cto_x2 + Index((n-1)*(n-1)*(kk-1) + (n-1)*jj + ii-1));
          node->add_halo(int(ranks.back()), create_halo0(2));
        }
        // lower north east neighbor
        if((ii+1 < n) && (jj+1 < n) && (kk > 0))
        {
          ranks.push_back(Index(n*n*(kk-1) + n*(jj+1) + ii+1));
          ctags.push_back(cto_x3 + Index((n-1)*(n-1)*(kk-1) + (n-1)*jj + ii));
          node->add_halo(int(ranks.back()), create_halo0(3));
        }

        // upper south west
        if((ii > 0) && (jj > 0) && (kk+1 < n))
        {
          ranks.push_back(Index(n*n*(kk+1) + n*(jj-1) + ii-1));
          ctags.push_back(cto_x3 + Index((n-1)*(n-1)*(kk) + (n-1)*(jj-1) + ii-1));
          node->add_halo(int(ranks.back()), create_halo0(4));
        }
        // upper south east neighbor
        if((ii+1 < n) && (jj > 0) && (kk+1 < n))
        {
          ranks.push_back(Index(n*n*(kk+1) + n*(jj-1) + ii+1));
          ctags.push_back(cto_x2 + Index((n-1)*(n-1)*(kk) + (n-1)*(jj-1) + ii));
          node->add_halo(int(ranks.back()), create_halo0(5));
        }
        // upper north west neighbor
        if((ii > 0) && (jj+1 < n) && (kk+1 < n))
        {
          ranks.push_back(Index(n*n*(kk+1) + n*(jj+1) + ii-1));
          ctags.push_back(cto_x1 + Index((n-1)*(n-1)*(kk) + (n-1)*jj + ii-1));
          node->add_halo(int(ranks.back()), create_halo0(6));
        }
        // upper north east neighbor
        if((ii+1 < n) && (jj+1 < n) && (kk+1 < n))
        {
          ranks.push_back(Index(n*n*(kk+1) + n*(jj+1) + ii+1));
          ctags.push_back(cto_x0 + Index((n-1)*(n-1)*(kk) + (n-1)*jj + ii));
          node->add_halo(int(ranks.back()), create_halo0(7));
        }

        // lower south neighbor
        if((jj > 0) && (kk > 0))
        {
          ranks.push_back(Index(n*n*(kk-1) + n*(jj-1) + ii));
          ctags.push_back(cto_h0 + Index((n-1)*n*(kk-1) + n*(jj-1) + ii));
          node->add_halo(int(ranks.back()), create_halo1(0));
        }
        // lower north neighbor
        if((jj + 1 < n) && (kk > 0))
        {
          ranks.push_back(Index(n*n*(kk-1) + n*(jj+1) + ii));
          ctags.push_back(cto_h1 + Index((n-1)*n*(kk-1) + n*jj + ii));
          node->add_halo(int(ranks.back()), create_halo1(1));
        }
        // upper south neighbor
        if((jj > 0) && (kk + 1 < n))
        {
          ranks.push_back(Index(n*n*(kk+1) + n*(jj-1) + ii));
          ctags.push_back(cto_h1 + Index((n-1)*n*(kk) + n*(jj-1) + ii));
          node->add_halo(int(ranks.back()), create_halo1(2));
        }
        // upper north neighbor
        if((jj + 1 < n) && (kk + 1 < n))
        {
          ranks.push_back(Index(n*n*(kk+1) + n*(jj+1) + ii));
          ctags.push_back(cto_h0 + Index((n-1)*n*(kk) + n*jj + ii));
          node->add_halo(int(ranks.back()), create_halo1(3));
        }

        // lower west neighbor
        if((ii > 0) && (kk > 0))
        {
          ranks.push_back(Index(n*n*(kk-1) + n*jj + ii-1));
          ctags.push_back(cto_v0 + Index((n-1)*n*(kk-1) + n*jj + ii-1));
          node->add_halo(int(ranks.back()), create_halo1(4));
        }
        // lower east neighbor
        if((ii + 1 < n) && (kk > 0))
        {
          ranks.push_back(Index(n*n*(kk-1) + n*jj + ii+1));
          ctags.push_back(cto_v1 + Index((n-1)*n*(kk-1) + n*jj + ii));
          node->add_halo(int(ranks.back()), create_halo1(5));
        }
        // upper west neighbor
        if((ii > 0) && (kk + 1 < n))
        {
          ranks.push_back(Index(n*n*(kk+1) + n*jj + ii-1));
          ctags.push_back(cto_v1 + Index((n-1)*n*(kk) + n*jj + ii-1));
          node->add_halo(int(ranks.back()), create_halo1(6));
        }
        // upper east neighbor
        if((ii + 1 < n) && (kk + 1 < n))
        {
          ranks.push_back(Index(n*n*(kk+1) + n*jj + ii+1));
          ctags.push_back(cto_v0 + Index((n-1)*n*(kk) + n*jj + ii));
          node->add_halo(int(ranks.back()), create_halo1(7));
        }

        // south west neighbor
        if((ii > 0) && (jj > 0))
        {
          ranks.push_back(Index(n*n*kk + n*(jj-1) + ii-1));
          ctags.push_back(cto_z0 + Index((n-1)*(n-1)*(kk) + (n-1)*(jj-1) + ii-1));
          node->add_halo(int(ranks.back()), create_halo1(8));
        }
        // south east neighbor
        if((ii + 1 < n) && (jj > 0))
        {
          ranks.push_back(Index(n*n*kk + n*(jj-1) + ii+1));
          ctags.push_back(cto_z1 + Index((n-1)*(n-1)*(kk) + (n-1)*(jj-1) + ii));
          node->add_halo(int(ranks.back()), create_halo1(9));
        }
        // north west neighbor
        if((ii > 0) && (jj + 1 < n))
        {
          ranks.push_back(Index(n*n*kk + n*(jj+1) + ii-1));
          ctags.push_back(cto_z1 + Index((n-1)*(n-1)*(kk) + (n-1)*jj + ii-1));
          node->add_halo(int(ranks.back()), create_halo1(10));
        }
        // north east neighbor
        if((ii + 1 < n) && (jj + 1 < n))
        {
          ranks.push_back(Index(n*n*kk + n*(jj+1) + ii+1));
          ctags.push_back(cto_z0 + Index((n-1)*(n-1)*(kk) + (n-1)*jj + ii));
          node->add_halo(int(ranks.back()), create_halo1(11));
        }

        // lower neighbor
        if(kk > 0)
        {
          ranks.push_back(Index(n*n*(kk-1) + n*jj + ii));
          ctags.push_back(cto_f0 + Index(n*n*(kk-1) + jj*n + ii));
          node->add_halo(int(ranks.back()), create_halo2(0));
          node->add_mesh_part("bnd:0", nullptr);
        }
        else
          node->add_mesh_part("bnd:0", create_halo2(0));
        // upper neighbor
        if(kk + 1 < n)
        {
          ranks.push_back(Index(n*n*(kk+1) + n*jj + ii));
          ctags.push_back(cto_f0 + Index(n*n*kk + jj*n + ii));
          node->add_halo(int(ranks.back()), create_halo2(1));
          node->add_mesh_part("bnd:1", nullptr);
        }
        else
          node->add_mesh_part("bnd:1", create_halo2(1));
        // south neighbor
        if(jj > 0)
        {
          ranks.push_back(Index(n*n*kk + n*(jj-1) + ii));
          ctags.push_back(cto_f1 + Index(n*n*kk + (jj-1)*n + ii));
          node->add_halo(int(ranks.back()), create_halo2(2));
          node->add_mesh_part("bnd:2", nullptr);
        }
        else
          node->add_mesh_part("bnd:2", create_halo2(2));
        // north neighbor
        if(jj + 1 < n)
        {
          ranks.push_back(Index(n*n*kk + n*(jj+1) + ii));
          ctags.push_back(cto_f1 + Index(n*n*kk + jj*n + ii));
          node->add_halo(int(ranks.back()), create_halo2(3));
          node->add_mesh_part("bnd:3", nullptr);
        }
        else
          node->add_mesh_part("bnd:3", create_halo2(3));
        // west neighbor
        if(ii > 0)
        {
          ranks.push_back(Index(n*n*kk + n*jj + ii-1));
          ctags.push_back(cto_f2 + Index(n*n*kk + jj*n + ii-1));
          node->add_halo(int(ranks.back()), create_halo2(4));
          node->add_mesh_part("bnd:4", nullptr);
        }
        else
          node->add_mesh_part("bnd:4", create_halo2(4));
        // east neighbor
        if(ii + 1 < n)
        {
          ranks.push_back(Index(n*n*kk + n*jj + ii+1));
          ctags.push_back(cto_f2 + Index(n*n*kk + jj*n + ii));
          node->add_halo(int(ranks.back()), create_halo2(5));
          node->add_mesh_part("bnd:5", nullptr);
        }
        else
          node->add_mesh_part("bnd:5", create_halo2(5));

        return level;
      }
    }; // UnitCubePatchGenerator<Hypercube<3>>
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_UNIT_CUBE_PATCH_GENERATOR_HPP
