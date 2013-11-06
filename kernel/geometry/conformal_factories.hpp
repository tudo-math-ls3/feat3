#pragma once
#ifndef KERNEL_GEOMETRY_CONFORMAL_FACTORIES_HPP
#define KERNEL_GEOMETRY_CONFORMAL_FACTORIES_HPP 1

// includes, FEAST
#include <kernel/geometry/cell_sub_set.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/structured_mesh.hpp>

namespace FEAST
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
    template<int stride_, typename Coord_>
    class UnitCubeFactory< ConformalMesh<Shape::Hypercube<1>, 1, stride_, Coord_> > :
      public Factory< ConformalMesh<Shape::Hypercube<1>, 1, stride_, Coord_> >
    {
    public:
      /// mesh type
      typedef ConformalMesh<Shape::Hypercube<1>, 1, stride_, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    public:
      UnitCubeFactory()
      {
      }

      virtual Index get_num_entities(int dim)
      {
        return Index(2 - dim);
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set)
      {
        vertex_set[0][0] = Coord_(0);
        vertex_set[1][0] = Coord_(1);
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
      {
        IndexSet<2>& v_e(index_set_holder.template get_index_set<1,0>());
        v_e[0][0] = 0;
        v_e[0][1] = 1;
      }
    };

    template<int stride_, typename Coord_>
    class UnitCubeFactory< ConformalMesh<Shape::Hypercube<2>, 2, stride_, Coord_> > :
      public Factory< ConformalMesh<Shape::Hypercube<2>, 2, stride_, Coord_> >
    {
    public:
      /// shape type
      typedef Shape::Hypercube<2> ShapeType;
      /// mesh type
      typedef ConformalMesh<Shape::Hypercube<2>, 2, stride_, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    public:
      UnitCubeFactory()
      {
      }

      virtual Index get_num_entities(int dim)
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

      virtual void fill_vertex_set(VertexSetType& vertex_set)
      {
        for(Index i(0); i < 4u; ++i)
        {
          for(int j(0); j < 2; ++j)
          {
            vertex_set[i][j] = Coord_((i >> j) & 0x1);
          }
        }
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
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

        for(Index i(0); i < Index(Shape::FaceTraits<ShapeType, cell_dim_>::count); ++i)
        {
          for(int j(0); j < idx.num_indices; ++j)
          {
            idx[i][j] = Index(FimType::map(int(i), j));
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

    template<int stride_, typename Coord_>
    class UnitCubeFactory< ConformalMesh<Shape::Hypercube<3>, 3, stride_, Coord_> > :
      public Factory< ConformalMesh<Shape::Hypercube<3>, 3, stride_, Coord_> >
    {
    public:
      /// shape type
      typedef Shape::Hypercube<3> ShapeType;
      /// mesh type
      typedef ConformalMesh<Shape::Hypercube<3>, 3, stride_, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    public:
      UnitCubeFactory()
      {
      }

      virtual Index get_num_entities(int dim)
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

      virtual void fill_vertex_set(VertexSetType& vertex_set)
      {
        for(Index i(0); i < 8u; ++i)
        {
          for(int j(0); j < 3; ++j)
          {
            vertex_set[i][j] = Coord_((i >> j) & 0x1);
          }
        }
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
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

        for(Index i(0); i < Index(Shape::FaceTraits<ShapeType, cell_dim_>::count); ++i)
        {
          for(int j(0); j < idx.num_indices; ++j)
          {
            idx[i][j] = Index(FimType::map(int(i), j));
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
          idx[0][j] = Index(j);
        }
      }
    };
    /// \endcond

    template<>
    class UnitCubeFactory< CellSubSet<Shape::Hypercube<2> > > :
      public Factory< CellSubSet<Shape::Hypercube<2> > >
    {
    public:
      typedef Shape::Hypercube<2> ShapeType;
      /// mesh type
      typedef CellSubSet<ShapeType> MeshType;
      /// target set holder type
      typedef MeshType::TargetSetHolderType TargetSetHolderType;

    public:
      virtual Index get_num_entities(int dim)
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

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder)
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
    };


    /**
     * \brief Refined Unit-Cube mesh factory
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
    class RefinedUnitCubeFactory DOXY({});

    /// \cond internal
    template<typename Shape_, int num_coords_, int stride_, typename Coord_>
    class RefinedUnitCubeFactory< ConformalMesh<Shape_, num_coords_, stride_, Coord_> > :
      public Factory< ConformalMesh<Shape_, num_coords_, stride_, Coord_> >
    {
    public:
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
      typedef typename MeshType::VertexSetType VertexSetType;
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    private:
      typedef Factory<MeshType> MeshFactory;
      typedef UnitCubeFactory<MeshType> CubeFactory;
      typedef StandardRefinery<MeshType> RefineFactory;

      MeshType* _coarse_mesh;
      MeshFactory* _factory;

    public:
      explicit RefinedUnitCubeFactory(Index num_refines) :
        _coarse_mesh(nullptr),
        _factory(nullptr)
      {
        if(num_refines <= 0)
        {
          _factory = new CubeFactory();
          return;
        }

        // create coarse mesh
        CubeFactory cube_factory;
        _coarse_mesh = new MeshType(cube_factory);

        // create refinery
        _factory = new RefineFactory(*_coarse_mesh);

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
          _factory = new RefineFactory(*_coarse_mesh);
        }
      }

      virtual ~RefinedUnitCubeFactory()
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

      virtual Index get_num_entities(int dim)
      {
        return _factory->get_num_entities(dim);
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set)
      {
        _factory->fill_vertex_set(vertex_set);
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
      {
        _factory->fill_index_sets(index_set_holder);
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONFORMAL_FACTORIES_HPP
