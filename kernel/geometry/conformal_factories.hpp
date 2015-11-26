#pragma once
#ifndef KERNEL_GEOMETRY_CONFORMAL_FACTORIES_HPP
#define KERNEL_GEOMETRY_CONFORMAL_FACTORIES_HPP 1

// includes, FEAST
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/geometry/shape_convert_factory.hpp>

#include <deque>
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

            for(int i(0); i < Shape::FaceTraits<ShapeType, cell_dim_>::count; ++i)
            {
              for(int j(0); j < idx.num_indices; ++j)
              {
                idx[Index(i)][Index(j)] = Index(FimType::map(i, j));
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
              idx[0][Index(j)] = Index(j);
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

            for(int i(0); i < Shape::FaceTraits<ShapeType, cell_dim_>::count; ++i)
            {
              for(int j(0); j < idx.num_indices; ++j)
              {
                idx[Index(i)][Index(j)] = Index(FimType::map(i, j));
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
              idx[Index(0)][Index(j)] = Index(j);
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

      }; //UnitCubeFactory< MeshPart<ConformalMesh<Shape::Hypercube<2> > > >

    /// \endcond

    /**
     * \brief Unit-Cube mesh factory specialisation for simplical meshes
     *
     * This class template implements the mesh factory interface which generates a simple 1-cell
     * 1D/2D/3D unit-cube mesh by using the UnitCubeFactories for hypercube meshes and the
     * ShapeConvertFactory.
     *
     * \author Jordi Paul
     */
    template<int stride_, typename Coord_, int dim_>
    class UnitCubeFactory< ConformalMesh<Shape::Simplex<dim_>, dim_, stride_, Coord_> > :
    public Factory< ConformalMesh<Shape::Simplex<dim_>, dim_, stride_, Coord_> >
    {
      public:
        /// mesh type
        typedef ConformalMesh<Shape::Simplex<dim_>, dim_, stride_, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
      private:
        /// shape convert factory type
        typedef ShapeConvertFactory<MeshType> FactoryType;
        /// mesh type to convert from
        typedef ConformalMesh<Shape::Hypercube<dim_>, dim_, stride_, Coord_> GeneratorMeshType;

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

      private:
        GeneratorMeshType* _generator_mesh;
        FactoryType* _factory;

    };
    /// \endcond

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
    template<typename Shape_, int num_coords_, int stride_, typename Coord_, template<typename> class Factory_>
    class RefineFactory< ConformalMesh<Shape_, num_coords_, stride_, Coord_>, Factory_ > :
    public Factory< ConformalMesh<Shape_, num_coords_, stride_, Coord_> >
    {
      public:
        typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
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
     * \brief Specialisation for Hypercube<2> meshes
     **/
    template<int stride_, typename Coord_>
    class UnitStarCubeFactory< ConformalMesh<Shape::Hypercube<2>, 2, stride_, Coord_> > :
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
        UnitStarCubeFactory()
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
              return 5u;
            default:
              return 0u;
          }
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set)
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

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
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
     * \brief Specialisation for simplical meshes
     *
     * This uses the UnitStarCubeFactory for Hypercubes and then the ShapeConvertFactory.
     *
     **/
    template<int stride_, typename Coord_, int dim_>
    class UnitStarCubeFactory< ConformalMesh<Shape::Simplex<dim_>, dim_, stride_, Coord_> > :
    public Factory< ConformalMesh<Shape::Simplex<dim_>, dim_, stride_, Coord_> >
    {
      public:
        /// shape type
        typedef Shape::Simplex<dim_> ShapeType;
        /// mesh type
        typedef ConformalMesh<Shape::Simplex<dim_>, dim_, stride_, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
        /// shape convert factory type
        typedef ShapeConvertFactory<MeshType> FactoryType;
        /// mesh type to convert from
        typedef ConformalMesh<Shape::Hypercube<dim_>, dim_, stride_, Coord_> GeneratorMeshType;

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

    /**
     * \brief Constructs a polyline mesh
     *
     * \tparam[dim_]
     * Dimension of the input points
     *
     * \tparam[stride_]
     * Padded dimension of the input points
     *
     * \tparam[Coord_]
     * Floating point type for mesh coordinates
     *
     * \author Jordi Paul
     *
     * This simply joins all points of a given std::deque together to a polygonal straight line graph.
     *
     */
    template<int dim_, int stride_, typename Coord_>
    class PolylineFactory :
    public Factory< ConformalMesh<Shape::Hypercube<1>, dim_, stride_, Coord_> >
    {
      public:
        /// mesh type
        typedef ConformalMesh<Shape::Hypercube<1>, dim_, stride_, Coord_> MeshType;
        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      private:
        /// Reference to the set of points in the polyline
        std::deque<Tiny::Vector<Coord_, dim_, stride_>>& _points;

      public:

        /**
         * \brief From deque-of-Tiny::Vectors constructor
         */
        PolylineFactory(std::deque<typename VertexSetType::VertexType>& points_) :
          _points(points_)
        {
          if(points_.empty())
            throw InternalError("PolylineFactory constructor called on empty point set!");

        }

        /**
         * \copydoc Factory::get_num_entities(int)
         */
        virtual Index get_num_entities(int dimension)
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

        virtual void fill_vertex_set(VertexSetType& vertex_set)
        {
          Index i(0);
          const auto& jt(_points.end());
          for(auto it(_points.begin()); it != jt; ++i)
          {
            vertex_set[i] = *it;
            it++;
          }
        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
        {
          IndexSet<2>& v_e(index_set_holder.template get_index_set<1,0>());
          for(Index i(0); i < get_num_entities(1); ++i)
          {
            v_e[i][0] = i;
            v_e[i][1] = i + Index(1);
          }
        }

    };

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONFORMAL_FACTORIES_HPP
