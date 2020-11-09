// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_COMMON_FACTORIES_HPP
#define KERNEL_GEOMETRY_COMMON_FACTORIES_HPP 1

// includes, FEAT
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/geometry/shape_convert_factory.hpp>

#include <deque>
namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Unit-Cube mesh factory
     *
     * This class template implements the mesh factory interface which generates a simple 1-cell
     * 1D/2D/3D unit-cube mesh.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
    class UnitCubeFactory DOXY({});

    /// \cond internal
    template<typename Coord_>
    class UnitCubeFactory< ConformalMesh<Shape::Hypercube<1>, 1, Coord_> > :
      public Factory< ConformalMesh<Shape::Hypercube<1>, 1, Coord_> >
    {
      public:
        /// mesh type
        typedef ConformalMesh<Shape::Hypercube<1>, 1, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      public:
        UnitCubeFactory()
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          return Index(2 - dim);
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          vertex_set[0][0] = Coord_(0);
          vertex_set[1][0] = Coord_(1);
        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
        {
          IndexSet<2>& v_e(index_set_holder.template get_index_set<1,0>());
          v_e[0][0] = 0;
          v_e[0][1] = 1;
        }

    };

    template<typename Coord_>
    class UnitCubeFactory< ConformalMesh<Shape::Hypercube<2>, 2, Coord_> > :
      public Factory< ConformalMesh<Shape::Hypercube<2>, 2, Coord_> >
    {
      public:
        /// shape type
        typedef Shape::Hypercube<2> ShapeType;
        /// mesh type
        typedef ConformalMesh<Shape::Hypercube<2>, 2, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      public:
        UnitCubeFactory()
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          switch(dim)
          {
            case 0:
              return 4u;
            case 1:
              return 4u;
            case 2:
              return 1u;
            default:
              return 0u;
          }
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          for(Index i(0); i < 4u; ++i)
          {
            for(int j(0); j < 2; ++j)
            {
              vertex_set[i][j] = Coord_((i >> j) & 0x1);
            }
          }
        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
        {
          _fill_sub_index_set<1,0>(index_set_holder);
          _fill_cell_index_set<0>(index_set_holder);
          _fill_cell_index_set<1>(index_set_holder);
        }

      private:
        template<int cell_dim_, int face_dim_>
          static void _fill_sub_index_set(IndexSetHolderType& index_set_holder)
          {
            typedef typename Intern::FaceIndexMapping<ShapeType, cell_dim_, face_dim_> FimType;
            typename IndexSetHolderType::template IndexSet<cell_dim_, face_dim_>::Type&
              idx(index_set_holder.template get_index_set<cell_dim_, face_dim_>());

            for(int i(0); i < Shape::FaceTraits<ShapeType, cell_dim_>::count; ++i)
            {
              for(int j(0); j < idx.num_indices; ++j)
              {
                idx[Index(i)][j] = Index(FimType::map(i, j));
              }
            }
          }

        template<int face_dim_>
          static void _fill_cell_index_set(IndexSetHolderType& index_set_holder)
          {
            typename IndexSetHolderType::template IndexSet<2, face_dim_>::Type&
              idx(index_set_holder.template get_index_set<2, face_dim_>());
            for(int j(0); j < idx.num_indices; ++j)
            {
              idx[0][j] = Index(j);
            }
          }
    };

    template<typename Coord_>
    class UnitCubeFactory< ConformalMesh<Shape::Hypercube<3>, 3, Coord_> > :
      public Factory< ConformalMesh<Shape::Hypercube<3>, 3, Coord_> >
    {
      public:
        /// shape type
        typedef Shape::Hypercube<3> ShapeType;
        /// mesh type
        typedef ConformalMesh<Shape::Hypercube<3>, 3, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      public:
        UnitCubeFactory()
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          switch(dim)
          {
            case 0:
              return 8u;
            case 1:
              return 12u;
            case 2:
              return 6u;
            case 3:
              return 1u;
            default:
              return 0u;
          }
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          for(Index i(0); i < 8u; ++i)
          {
            for(int j(0); j < 3; ++j)
            {
              vertex_set[i][j] = Coord_((i >> j) & 0x1);
            }
          }
        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
        {
          _fill_sub_index_set<1,0>(index_set_holder);
          _fill_sub_index_set<2,0>(index_set_holder);
          _fill_sub_index_set<2,1>(index_set_holder);
          _fill_cell_index_set<0>(index_set_holder);
          _fill_cell_index_set<1>(index_set_holder);
          _fill_cell_index_set<2>(index_set_holder);
        }

      private:
        template<int cell_dim_, int face_dim_>
          static void _fill_sub_index_set(IndexSetHolderType& index_set_holder)
          {
            typedef typename Intern::FaceIndexMapping<ShapeType, cell_dim_, face_dim_> FimType;
            typename IndexSetHolderType::template IndexSet<cell_dim_, face_dim_>::Type&
              idx(index_set_holder.template get_index_set<cell_dim_, face_dim_>());

            for(int i(0); i < Shape::FaceTraits<ShapeType, cell_dim_>::count; ++i)
            {
              for(int j(0); j < idx.num_indices; ++j)
              {
                idx[Index(i)][j] = Index(FimType::map(i, j));
              }
            }
          }

        template<int face_dim_>
          static void _fill_cell_index_set(IndexSetHolderType& index_set_holder)
          {
            typename IndexSetHolderType::template IndexSet<3, face_dim_>::Type&
              idx(index_set_holder.template get_index_set<3, face_dim_>());
            for(int j(0); j < idx.num_indices; ++j)
            {
              idx[Index(0)][j] = Index(j);
            }
          }
    };
    /// \endcond

    /// \cond internal
    template<>
      class UnitCubeFactory< MeshPart<ConformalMesh<Shape::Hypercube<2> > > >:
      public Factory< MeshPart<ConformalMesh<Shape::Hypercube<2> > > >
      {
        public:
          typedef Shape::Hypercube<2> ShapeType;
          /// mesh part type
          typedef MeshPart<ConformalMesh<ShapeType>> MeshType;
          /// target set holder type
          typedef MeshType::TargetSetHolderType TargetSetHolderType;

        public:
          virtual Index get_num_entities(int dim) override
          {
            switch(dim)
            {
              case 0:
                return 4;
              case 1:
                return 4;
              default:
                return 0;
            }
          }

          virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
          {
            // set vertex indices
            TargetSet& vi(target_set_holder.get_target_set<0>());
            vi[0] = 0;
            vi[1] = 1;
            vi[2] = 2;
            vi[3] = 3;

            // set edge indices
            TargetSet& ei(target_set_holder.get_target_set<1>());
            ei[0] = 0;
            ei[1] = 1;
            ei[2] = 2;
            ei[3] = 3;
          }

      }; //UnitCubeFactory< MeshPart<ConformalMesh<Shape::Hypercube<2> > > >

    /// \endcond

    /**
     * \brief Unit-Cube mesh factory specialization for simplical meshes
     *
     * This class template implements the mesh factory interface which generates a simple 1-cell
     * 1D/2D/3D unit-cube mesh by using the UnitCubeFactories for hypercube meshes and the
     * ShapeConvertFactory.
     *
     * \author Jordi Paul
     */
    template<typename Coord_, int dim_>
    class UnitCubeFactory< ConformalMesh<Shape::Simplex<dim_>, dim_, Coord_> > :
      public Factory< ConformalMesh<Shape::Simplex<dim_>, dim_, Coord_> >
    {
      public:
        /// mesh type
        typedef ConformalMesh<Shape::Simplex<dim_>, dim_, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
      private:
        /// shape convert factory type
        typedef ShapeConvertFactory<MeshType> FactoryType;
        /// mesh type to convert from
        typedef ConformalMesh<Shape::Hypercube<dim_>, dim_, Coord_> GeneratorMeshType;

      public:
        UnitCubeFactory() :
          _generator_mesh(nullptr),
          _factory(nullptr)
        {
          UnitCubeFactory<GeneratorMeshType> cube_factory;
          _generator_mesh = new GeneratorMeshType(cube_factory);
          _factory = new FactoryType(*_generator_mesh);
        }

        virtual ~UnitCubeFactory()
        {
          delete _generator_mesh;
          delete _factory;
        }

        virtual Index get_num_entities(int dim) override
        {
          return _factory->get_num_entities(dim);
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          _factory->fill_vertex_set(vertex_set);
        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
        {
          _factory->fill_index_sets(index_set_holder);
        }

      private:
        GeneratorMeshType* _generator_mesh;
        FactoryType* _factory;

    };
    /// \endcond


    template<int dim_, typename Coord_>
    class UnitCubeFactory< StructuredMesh<dim_, dim_, Coord_> > :
      public Factory< StructuredMesh<dim_, dim_, Coord_> >
    {
    public:
      /// shape type
      typedef Shape::Hypercube<dim_> ShapeType;
      /// mesh type
      typedef StructuredMesh<dim_, dim_, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;

    protected:
      const Index _level;

    public:
      explicit UnitCubeFactory(Index level = Index(0)) :
        _level(level)
      {
      }

      virtual Index get_num_slices(int DOXY(dir)) override
      {
        return (Index(1) << _level); // = 2^level in each direction
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        // coordinate scaling factor
        const Coord_ sc = Coord_(1) / Coord_(Index(1) << _level); // = 2^{-level}

        // number of vertices in each direction
        const Index nx = (Index(1) << _level) + Index(1); // = 2^level + 1

        // total number of vertices = nx ^ dim
        const Index nv = vertex_set.get_num_vertices();

        for(Index i(0); i < nv; ++i)
        {
          Index k = i;
          auto& v = vertex_set[i];
          for(int j(0); j < dim_; ++j, k /= nx)
          {
            v[j] = Coord_(k % nx) * sc;
          }
        }
      }
    };

    /**
     * \brief Refine mesh factory
     *
     * \tparam Mesh_
     * Mesh type
     *
     * \tparam Factory_
     * Factory type to create the initial mesh
     *
     * This uses a factory to create a mesh and then refines that mesh.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_, template<typename> class Factory_>
    class RefineFactory DOXY({});

    /// \cond internal
    template<typename Shape_, int num_coords_, typename Coord_, template<typename> class Factory_>
    class RefineFactory< ConformalMesh<Shape_, num_coords_, Coord_>, Factory_ > :
      public Factory< ConformalMesh<Shape_, num_coords_, Coord_> >
    {
      public:
        typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
        typedef typename MeshType::VertexSetType VertexSetType;
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      private:
        typedef Factory<MeshType> MeshFactory;
        typedef Factory_<MeshType> MyFactoryType;
        typedef StandardRefinery<MeshType> Refinery;

        MeshType* _coarse_mesh;
        MeshFactory* _factory;

      public:
        template<typename... Arguments>
        explicit RefineFactory(Index num_refines, Arguments&&... args) :
          _coarse_mesh(nullptr),
          _factory(nullptr)
        {
          if(num_refines <= 0)
          {
            _factory = new MyFactoryType(std::forward<Arguments>(args)...);
            return;
          }

          // create coarse mesh
          MyFactoryType my_factory(std::forward<Arguments>(args)...);
          _coarse_mesh = new MeshType(my_factory);

          // create refinery
          _factory = new Refinery(*_coarse_mesh);

          // refine n-1 times;
          for(Index i(1); i < num_refines; ++i)
          {
            // backup old mesh
            MeshType* mesh_old = _coarse_mesh;
            // refine mesh
            _coarse_mesh = new MeshType(*_factory);
            // delete old factory
            delete _factory;
            // delete old coarse mesh
            delete mesh_old;
            // create new factory
            _factory = new Refinery(*_coarse_mesh);
          }
        }

        virtual ~RefineFactory()
        {
          if(_factory != nullptr)
          {
            delete _factory;
          }
          if(_coarse_mesh != nullptr)
          {
            delete _coarse_mesh;
          }
        }

        virtual Index get_num_entities(int dim) override
        {
          return _factory->get_num_entities(dim);
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          _factory->fill_vertex_set(vertex_set);
        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
        {
          _factory->fill_index_sets(index_set_holder);
        }
    };

    template<int shape_dim_, int num_coords_, typename Coord_, template<typename> class Factory_>
    class RefineFactory< StructuredMesh<shape_dim_, num_coords_, Coord_>, Factory_ > :
      public Factory< StructuredMesh<shape_dim_, num_coords_, Coord_> >
    {
      public:
        typedef StructuredMesh<shape_dim_, num_coords_, Coord_> MeshType;
        typedef typename MeshType::VertexSetType VertexSetType;

      private:
        typedef Factory<MeshType> MeshFactory;
        typedef Factory_<MeshType> MyFactoryType;
        typedef StandardRefinery<MeshType> Refinery;

        MeshType* _coarse_mesh;
        MeshFactory* _factory;

      public:
        template<typename... Arguments>
        explicit RefineFactory(Index num_refines, Arguments&&... args) :
          _coarse_mesh(nullptr),
          _factory(nullptr)
        {
          if(num_refines <= 0)
          {
            _factory = new MyFactoryType(std::forward<Arguments>(args)...);
            return;
          }

          // create coarse mesh
          MyFactoryType my_factory(std::forward<Arguments>(args)...);
          _coarse_mesh = new MeshType(my_factory);

          // create refinery
          _factory = new Refinery(*_coarse_mesh);

          // refine n-1 times;
          for(Index i(1); i < num_refines; ++i)
          {
            // backup old mesh
            MeshType* mesh_old = _coarse_mesh;
            // refine mesh
            _coarse_mesh = new MeshType(*_factory);
            // delete old factory
            delete _factory;
            // delete old coarse mesh
            delete mesh_old;
            // create new factory
            _factory = new Refinery(*_coarse_mesh);
          }
        }

        virtual ~RefineFactory()
        {
          if(_factory != nullptr)
          {
            delete _factory;
          }
          if(_coarse_mesh != nullptr)
          {
            delete _coarse_mesh;
          }
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          _factory->fill_vertex_set(vertex_set);
        }

        virtual Index get_num_slices(int dir) override
        {
          return _factory->get_num_slices(dir);
        }
    };
    /// \endcond

    template<typename MeshType_>
    using RefinedUnitCubeFactory = RefineFactory<MeshType_, Geometry::UnitCubeFactory >;

    /**
     * \brief Unit cube factory with star shaped mesh topology
     *
     * 2-------------------------------1------------------------------3
     * |.\                                                          /.|
     * |  .\                                                      /.  |
     * |    .\                                                  /.    |
     * |      .\                                              /.      |
     * |        11                     2                    10        |
     * |          .\                                      /.          |
     * |            .\                                  /.            |
     * |              .\                              /.              |
     * |               6---------------5--------------7               |
     * |               |                              |               |
     * |               |                              |               |
     * |               |                              |               |
     * |               |                              |               |
     * |               |                              |               |
     * |               |                              |               |
     * 2       3       6               4              7       1       3
     * |               |                              |               |
     * |               |                              |               |
     * |               |                              |               |
     * |               |                              |               |
     * |               |                              |               |
     * |               4---------------4--------------5               |
     * |               /.                             .\              |
     * |             /.                                 .\            |
     * |           /.                                     .\          |
     * |         8                                           9        |
     * |       /.                      0                      .\      |
     * |     /.                                                 .\    |
     * |   /.                                                     .\  |
     * | /.                                                         .\|
     * 0-------------------------------0------------------------------1
     *
     * \author Jordi Paul
     *
     **/
    template<typename Mesh_>
    class UnitStarCubeFactory DOXY({});

    /// \cond internal
    /*
     * \brief Specialization for Hypercube<2> meshes
     **/
    template<typename Coord_>
    class UnitStarCubeFactory< ConformalMesh<Shape::Hypercube<2>, 2, Coord_> > :
    public Factory< ConformalMesh<Shape::Hypercube<2>, 2, Coord_> >
    {
      public:
        /// shape type
        typedef Shape::Hypercube<2> ShapeType;
        /// mesh type
        typedef ConformalMesh<Shape::Hypercube<2>, 2, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      public:
        UnitStarCubeFactory()
        {
        }

        virtual Index get_num_entities(int dim) override
        {
          switch(dim)
          {
            case 0:
              return 8u;
            case 1:
              return 12u;
            case 2:
              return 5u;
            default:
              return 0u;
          }
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          vertex_set[0][0] = Coord_(0);
          vertex_set[0][1] = Coord_(0);

          vertex_set[1][0] = Coord_(1);
          vertex_set[1][1] = Coord_(0);

          vertex_set[2][0] = Coord_(0);
          vertex_set[2][1] = Coord_(1);

          vertex_set[3][0] = Coord_(1);
          vertex_set[3][1] = Coord_(1);

          vertex_set[4][0] = Coord_(0.25);
          vertex_set[4][1] = Coord_(0.25);

          vertex_set[5][0] = Coord_(0.75);
          vertex_set[5][1] = Coord_(0.25);

          vertex_set[6][0] = Coord_(0.25);
          vertex_set[6][1] = Coord_(0.75);

          vertex_set[7][0] = Coord_(0.75);
          vertex_set[7][1] = Coord_(0.75);

        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
        {
          IndexSet<4>& v_q = index_set_holder.template get_index_set<2,0>();

          v_q(0,0) = Index(0);
          v_q(0,1) = Index(1);
          v_q(0,2) = Index(4);
          v_q(0,3) = Index(5);

          v_q(1,0) = Index(5);
          v_q(1,1) = Index(1);
          v_q(1,2) = Index(7);
          v_q(1,3) = Index(3);

          v_q(2,0) = Index(6);
          v_q(2,1) = Index(7);
          v_q(2,2) = Index(2);
          v_q(2,3) = Index(3);

          v_q(3,0) = Index(0);
          v_q(3,1) = Index(4);
          v_q(3,2) = Index(2);
          v_q(3,3) = Index(6);

          v_q(4,0) = Index(4);
          v_q(4,1) = Index(5);
          v_q(4,2) = Index(6);
          v_q(4,3) = Index(7);

          IndexSet<4>& e_q = index_set_holder.template get_index_set<2,1>();

          e_q(0,0) = Index(0);
          e_q(0,1) = Index(4);
          e_q(0,2) = Index(8);
          e_q(0,3) = Index(9);

          e_q(1,0) = Index(9);
          e_q(1,1) = Index(10);
          e_q(1,2) = Index(7);
          e_q(1,3) = Index(3);

          e_q(2,0) = Index(5);
          e_q(2,1) = Index(1);
          e_q(2,2) = Index(11);
          e_q(2,3) = Index(10);

          e_q(3,0) = Index(8);
          e_q(3,1) = Index(11);
          e_q(3,2) = Index(2);
          e_q(3,3) = Index(6);

          e_q(4,0) = Index(4);
          e_q(4,1) = Index(5);
          e_q(4,2) = Index(6);
          e_q(4,3) = Index(7);

          IndexSet<2>& v_e = index_set_holder.template get_index_set<1,0>();

          v_e(0,0) = Index(0);
          v_e(0,1) = Index(1);

          v_e(1,0) = Index(2);
          v_e(1,1) = Index(3);

          v_e(2,0) = Index(0);
          v_e(2,1) = Index(2);

          v_e(3,0) = Index(1);
          v_e(3,1) = Index(3);

          v_e(4,0) = Index(4);
          v_e(4,1) = Index(5);

          v_e(5,0) = Index(6);
          v_e(5,1) = Index(7);

          v_e(6,0) = Index(4);
          v_e(6,1) = Index(6);

          v_e(7,0) = Index(5);
          v_e(7,1) = Index(7);

          v_e(8,0) = Index(0);
          v_e(8,1) = Index(4);

          v_e(9,0) = Index(1);
          v_e(9,1) = Index(5);

          v_e(10,0) = Index(7);
          v_e(10,1) = Index(3);

          v_e(11,0) = Index(6);
          v_e(11,1) = Index(2);
        }

    };
    /// \endcond

    /// \cond internal
    /*
     * \brief Specialization for simplical meshes
     *
     * This uses the UnitStarCubeFactory for Hypercubes and then the ShapeConvertFactory.
     *
     **/
    template<typename Coord_, int dim_>
    class UnitStarCubeFactory< ConformalMesh<Shape::Simplex<dim_>, dim_, Coord_> > :
    public Factory< ConformalMesh<Shape::Simplex<dim_>, dim_, Coord_> >
    {
      public:
        /// shape type
        typedef Shape::Simplex<dim_> ShapeType;
        /// mesh type
        typedef ConformalMesh<Shape::Simplex<dim_>, dim_, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
        /// shape convert factory type
        typedef ShapeConvertFactory<MeshType> FactoryType;
        /// mesh type to convert from
        typedef ConformalMesh<Shape::Hypercube<dim_>, dim_, Coord_> GeneratorMeshType;

      private:
        GeneratorMeshType* _generator_mesh;
        FactoryType* _factory;

      public:
        UnitStarCubeFactory() :
          _generator_mesh(nullptr),
          _factory(nullptr)
        {
          UnitStarCubeFactory<GeneratorMeshType> cube_factory;
          _generator_mesh = new GeneratorMeshType(cube_factory);
          _factory = new FactoryType(*_generator_mesh);
        }

        virtual ~UnitStarCubeFactory()
        {
          delete _generator_mesh;
          delete _factory;
        }

        virtual Index get_num_entities(int dim) override
        {
          return _factory->get_num_entities(dim);
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          _factory->fill_vertex_set(vertex_set);
        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
        {
          _factory->fill_index_sets(index_set_holder);
        }

    };
    /// \endcond

    /**
     * \brief Constructs a polyline mesh
     *
     * \tparam dim_
     * Dimension of the input points
     *
     * \tparam Coord_
     * Floating point type for mesh coordinates
     *
     * \author Jordi Paul
     *
     * This simply joins all points of a given std::deque together to a polygonal straight line graph.
     *
     */
    template<int dim_, typename Coord_>
    class PolylineFactory :
    public Factory< ConformalMesh<Shape::Hypercube<1>, dim_, Coord_> >
    {
      public:
        /// mesh type
        typedef ConformalMesh<Shape::Hypercube<1>, dim_, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      private:
        /// Reference to the set of points in the polyline
        std::deque<Tiny::Vector<Coord_, dim_>>& _points;

      public:

        /**
         * \brief From deque of Tiny::Vectors constructor
         */
        explicit PolylineFactory(std::deque<typename VertexSetType::VertexType>& points_) :
          _points(points_)
        {
          XASSERTM(!points_.empty(), "PolylineFactory constructor called on empty point set!");

        }

        /**
         * \copydoc Factory::get_num_entities()
         */
        virtual Index get_num_entities(int dimension) override
        {
          switch(dimension)
          {
            case 0:
              return Index(_points.size());
            case 1:
              return Index(_points.size())-Index(1);
            default:
              return 0;
          }
        }

        /**
         * \copydoc Factory::fill_vertex_set()
         */
        virtual void fill_vertex_set(VertexSetType& vertex_set) override
        {
          Index i(0);
          const auto& jt(_points.end());
          for(auto it(_points.begin()); it != jt; ++i)
          {
            vertex_set[i] = *it;
            it++;
          }
        }

        /**
         * \copydoc Factory::fill_index_sets()
         */
        virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
        {
          IndexSet<2>& v_e(index_set_holder.template get_index_set<1,0>());
          for(Index i(0); i < get_num_entities(1); ++i)
          {
            v_e[i][0] = i;
            v_e[i][1] = i + Index(1);
          }
        }
    };

    template<typename Mesh_>
    class UnitSphereFactory DOXY({});

    /// \cond internal

    template<typename Coord_>
    class UnitSphereFactory< ConformalMesh<Shape::Simplex<2>, 3, Coord_> > :
      public Factory< ConformalMesh<Shape::Simplex<2>, 3, Coord_> >
    {
    public:
      /// shape type
      typedef Shape::Simplex<2> ShapeType;
      /// mesh type
      typedef ConformalMesh<Shape::Simplex<2>, 3, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    public:
      UnitSphereFactory()
      {
      }

      virtual Index get_num_entities(int dim) override
      {
        switch(dim)
        {
          case 0:
            return 12u;
          case 1:
            return 30u;
          case 2:
            return 20u;
          default:
            return 0u;
        }
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        const Coord_ s5 = Math::sqrt(Coord_(5));
        const Coord_ va = Coord_(0);
        const Coord_ vb = Math::sqrt(Coord_(2) / (Coord_(5) + s5));
        const Coord_ vc = (s5 + Coord_(1)) / Math::sqrt(Coord_(10) + Coord_(2)*s5);
        vertex_set[ 0][0] =  va;
        vertex_set[ 0][1] = -vb;
        vertex_set[ 0][2] =  vc;
        vertex_set[ 1][0] =  vc;
        vertex_set[ 1][1] =  va;
        vertex_set[ 1][2] =  vb;
        vertex_set[ 2][0] =  vc;
        vertex_set[ 2][1] =  va;
        vertex_set[ 2][2] = -vb;
        vertex_set[ 3][0] = -vc;
        vertex_set[ 3][1] =  va;
        vertex_set[ 3][2] = -vb;
        vertex_set[ 4][0] = -vc;
        vertex_set[ 4][1] =  va;
        vertex_set[ 4][2] =  vb;
        vertex_set[ 5][0] = -vb;
        vertex_set[ 5][1] =  vc;
        vertex_set[ 5][2] =  va;
        vertex_set[ 6][0] =  vb;
        vertex_set[ 6][1] =  vc;
        vertex_set[ 6][2] =  va;
        vertex_set[ 7][0] =  vb;
        vertex_set[ 7][1] = -vc;
        vertex_set[ 7][2] =  va;
        vertex_set[ 8][0] = -vb;
        vertex_set[ 8][1] = -vc;
        vertex_set[ 8][2] =  va;
        vertex_set[ 9][0] =  va;
        vertex_set[ 9][1] = -vb;
        vertex_set[ 9][2] = -vc;
        vertex_set[10][0] =  va;
        vertex_set[10][1] =  vb;
        vertex_set[10][2] = -vc;
        vertex_set[11][0] =  va;
        vertex_set[11][1] =  vb;
        vertex_set[11][2] =  vc;
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
      {
        // vertices-at-edge
        auto& ve = index_set_holder.template get_index_set<1,0>();
        ve[ 0][0] =  1;
        ve[ 0][1] =  2;
        ve[ 1][0] =  2;
        ve[ 1][1] =  6;
        ve[ 2][0] =  1;
        ve[ 2][1] =  6;
        ve[ 3][0] =  1;
        ve[ 3][1] =  7;
        ve[ 4][0] =  2;
        ve[ 4][1] =  7;
        ve[ 5][0] =  3;
        ve[ 5][1] =  4;
        ve[ 6][0] =  4;
        ve[ 6][1] =  5;
        ve[ 7][0] =  3;
        ve[ 7][1] =  5;
        ve[ 8][0] =  3;
        ve[ 8][1] =  8;
        ve[ 9][0] =  4;
        ve[ 9][1] =  8;
        ve[10][0] =  5;
        ve[10][1] =  6;
        ve[11][0] =  5;
        ve[11][1] = 11;
        ve[12][0] =  6;
        ve[12][1] = 11;
        ve[13][0] =  6;
        ve[13][1] = 10;
        ve[14][0] =  5;
        ve[14][1] = 10;
        ve[15][0] =  9;
        ve[15][1] = 10;
        ve[16][0] =  2;
        ve[16][1] = 10;
        ve[17][0] =  2;
        ve[17][1] =  9;
        ve[18][0] =  3;
        ve[18][1] =  9;
        ve[19][0] =  3;
        ve[19][1] = 10;
        ve[20][0] =  7;
        ve[20][1] =  8;
        ve[21][0] =  8;
        ve[21][1] =  9;
        ve[22][0] =  7;
        ve[22][1] =  9;
        ve[23][0] =  0;
        ve[23][1] =  7;
        ve[24][0] =  0;
        ve[24][1] =  8;
        ve[25][0] =  0;
        ve[25][1] = 11;
        ve[26][0] =  0;
        ve[26][1] =  1;
        ve[27][0] =  1;
        ve[27][1] = 11;
        ve[28][0] =  4;
        ve[28][1] = 11;
        ve[29][0] =  0;
        ve[29][1] =  4;

        // vertices-at-face
        auto& vf = index_set_holder.template get_index_set<2,0>();
        vf[ 0][0] =  1;
        vf[ 0][1] =  2;
        vf[ 0][2] =  6;
        vf[ 1][0] =  1;
        vf[ 1][1] =  7;
        vf[ 1][2] =  2;
        vf[ 2][0] =  3;
        vf[ 2][1] =  4;
        vf[ 2][2] =  5;
        vf[ 3][0] =  4;
        vf[ 3][1] =  3;
        vf[ 3][2] =  8;
        vf[ 4][0] =  6;
        vf[ 4][1] =  5;
        vf[ 4][2] = 11;
        vf[ 5][0] =  5;
        vf[ 5][1] =  6;
        vf[ 5][2] = 10;
        vf[ 6][0] =  9;
        vf[ 6][1] = 10;
        vf[ 6][2] =  2;
        vf[ 7][0] = 10;
        vf[ 7][1] =  9;
        vf[ 7][2] =  3;
        vf[ 8][0] =  7;
        vf[ 8][1] =  8;
        vf[ 8][2] =  9;
        vf[ 9][0] =  8;
        vf[ 9][1] =  7;
        vf[ 9][2] =  0;
        vf[10][0] = 11;
        vf[10][1] =  0;
        vf[10][2] =  1;
        vf[11][0] =  0;
        vf[11][1] = 11;
        vf[11][2] =  4;
        vf[12][0] =  6;
        vf[12][1] =  2;
        vf[12][2] = 10;
        vf[13][0] =  1;
        vf[13][1] =  6;
        vf[13][2] = 11;
        vf[14][0] =  3;
        vf[14][1] =  5;
        vf[14][2] = 10;
        vf[15][0] =  5;
        vf[15][1] =  4;
        vf[15][2] = 11;
        vf[16][0] =  2;
        vf[16][1] =  7;
        vf[16][2] =  9;
        vf[17][0] =  7;
        vf[17][1] =  1;
        vf[17][2] =  0;
        vf[18][0] =  3;
        vf[18][1] =  9;
        vf[18][2] =  8;
        vf[19][0] =  4;
        vf[19][1] =  8;
        vf[19][2] =  0;

        Geometry::RedundantIndexSetBuilder<Shape::Triangle>::compute(index_set_holder);
      }
    }; // class UnitSphereFactory<ConformalMesh<Shape::Simplex<2>, ...>>

    template<typename Coord_>
    class UnitSphereFactory< ConformalMesh<Shape::Hypercube<2>, 3, Coord_> > :
      public Factory< ConformalMesh<Shape::Hypercube<2>, 3, Coord_> >
    {
    public:
      /// shape type
      typedef Shape::Hypercube<2> ShapeType;
      /// mesh type
      typedef ConformalMesh<Shape::Hypercube<2>, 3, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    public:
      UnitSphereFactory()
      {
      }

      virtual Index get_num_entities(int dim) override
      {
        switch(dim)
        {
          case 0:
            return 8u;
          case 1:
            return 12u;
          case 2:
            return 6u;
          default:
            return 0u;
        }
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        const Coord_ scale = Coord_(2) / Math::sqrt(Coord_(3));
        for(Index i(0); i < 8u; ++i)
        {
          for(int j(0); j < 3; ++j)
          {
            vertex_set[i][j] = scale * (Coord_((i >> j) & 0x1) - Coord_(0.5));
          }
        }
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
      {
        _fill_index_set<1,0>(index_set_holder);
        _fill_index_set<2,0>(index_set_holder);
        _fill_index_set<2,1>(index_set_holder);
      }

    private:
      template<int cell_dim_, int face_dim_>
      static void _fill_index_set(IndexSetHolderType& index_set_holder)
      {
        typedef typename Intern::FaceIndexMapping<Shape::Hypercube<3>, cell_dim_, face_dim_> FimType;
        auto& idx = index_set_holder.template get_index_set<cell_dim_, face_dim_>();

        for(int i(0); i < Shape::FaceTraits<Shape::Hypercube<3>, cell_dim_>::count; ++i)
        {
          for(int j(0); j < idx.num_indices; ++j)
          {
            idx[Index(i)][j] = Index(FimType::map(i, j));
          }
        }
      }
    }; // class UnitSphereFactory<ConformalMesh<Shape::Hypercube<2>, ...>>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_COMMON_FACTORIES_HPP
