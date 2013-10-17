#pragma once
#ifndef KERNEL_GEOMETRY_MESH_READER_FACTORY_HPP
#define KERNEL_GEOMETRY_MESH_READER_FACTORY_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/conformal_sub_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>
#include <kernel/geometry/index_calculator.hpp>
#include <kernel/util/mesh_streamer.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct MeshStreamerIndexer;

      template<typename Shape_>
      struct MeshStreamerTargeter;
    }
    /// \endcond

    template<typename Mesh_>
    class MeshStreamerFactory DOXY({});

    /**
     * \brief MeshStreamerFactory implementation for ConformalMesh
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      int num_coords_,
      int stride_,
      typename Coord_>
    class MeshStreamerFactory< ConformalMesh<Shape_, num_coords_, stride_, Coord_> > :
      public Factory< ConformalMesh<Shape_, num_coords_, stride_, Coord_> >
    {
    public:
      /// mesh type
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    private:
      MeshStreamer& _mesh_reader;
      MeshStreamer::MeshDataContainer* _mesh_data;

    public:
      explicit MeshStreamerFactory(MeshStreamer& mesh_reader) :
        _mesh_reader(mesh_reader),
        _mesh_data(mesh_reader.get_mesh())
      {
      }

      virtual Index get_num_entities(int dim)
      {
        switch(dim)
        {
        case 0:
          return _mesh_data->vertex_count;
        case 1:
          return _mesh_data->edge_count;
        case 2:
          return _mesh_data->tria_count + _mesh_data->quad_count;
        case 3:
          return _mesh_data->tetra_count + _mesh_data->hexa_count;
        default:
          return 0;
        }
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set)
      {
        const Index num_vertices(Index(_mesh_data->coords.size()));
        for(Index i(0); i < num_vertices; ++i)
        {
          // get a reference to the corresponding vertex
          MeshStreamer::MeshDataContainer::CoordVec& vtx(_mesh_data->coords[i]);

          ASSERT(Index(vtx.size()) == Index(num_coords_), "Vertex coordinate count mismatch!");

          // copy vertex coordinates
          for(Index j(0); j < Omdex(num_coords_); ++j)
          {
            vertex_set[i][j] = Coord_(vtx[j]);
          }
        }
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
      {
        // call wrapper
        Intern::MeshStreamerIndexer<Shape_>::wrap(index_set_holder, _mesh_data->adjacencies);

        // build redundant index sets
        RedundantIndexSetBuilder<Shape_>::compute(index_set_holder);
      }
    }; // class MeshStreamerFactory<ConformalMesh<...>>


    /**
     * \brief MeshStreamerFactory implementation for CellSubSet
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class MeshStreamerFactory< CellSubSet<Shape_> > :
      public Factory< CellSubSet<Shape_> >
    {
    public:
      /// mesh type
      typedef CellSubSet<Shape_> MeshType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    private:
      MeshStreamer& _mesh_reader;
      String _name;
      MeshStreamer::BaseContainer* _target_data;

    public:
      explicit MeshStreamerFactory(MeshStreamer& mesh_reader, String name) :
        _mesh_reader(mesh_reader),
        _name(name),
        _target_data(nullptr)
      {
        MeshStreamer::MeshNode* root(_mesh_reader.get_root_mesh_node());
        ASSERT_(root != nullptr);

        // try to find a sub-mesh node
        MeshStreamer::MeshNode* sub_mesh_node(root->find_sub_mesh(name));
        if(sub_mesh_node != nullptr)
        {
          _target_data = &sub_mesh_node->mesh_data;
          return;
        }

        // try to find a cell-set node
        MeshStreamer::CellSetNode* cell_set_node(root->find_cell_set(name));
        if(cell_set_node != nullptr)
        {
          _target_data = &cell_set_node->cell_set;
          return;
        }

        // no child with 'name' found
        throw InternalError("No sub-mesh or cell-set found with name '" + name + "'");
      }

      virtual Index get_num_entities(int dim)
      {
        switch(dim)
        {
        case 0:
          return _target_data->vertex_count;
        case 1:
          return _target_data->edge_count;
        case 2:
          return _target_data->tria_count + _target_data->quad_count;
        case 3:
          return _target_data->tetra_count + _target_data->hexa_count;
        default:
          return 0;
        }
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder)
      {
        Intern::MeshStreamerTargeter<Shape_>::wrap(target_set_holder, _target_data->parent_indices);
      }
    }; // class MeshStreamerFactory<CellSubSet<...>>

    /**
     * \brief MeshStreamerFactory implementation for ConformalSubMesh
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename Coord_>
    class MeshStreamerFactory< ConformalSubMesh<Shape_, Coord_> > :
      public Factory< ConformalSubMesh<Shape_, Coord_> >
    {
    public:
      /// mesh typedef
      typedef ConformalSubMesh<Shape_, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index set holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    private:
      MeshStreamer& _mesh_reader;
      String _name;
      MeshStreamer::MeshDataContainer* _mesh_data;

    public:
      explicit MeshStreamerFactory(MeshStreamer& mesh_reader, String name) :
        _mesh_reader(mesh_reader),
        _name(name),
        _mesh_data(nullptr)
      {
        MeshStreamer::MeshNode* root(_mesh_reader.get_root_mesh_node());
        ASSERT_(root != nullptr);

        // try to find a sub-mesh node
        MeshStreamer::MeshNode* sub_mesh_node(root->find_sub_mesh(name));
        if(sub_mesh_node != nullptr)
        {
          _mesh_data = &sub_mesh_node->mesh_data;
          return;
        }

        // no child with 'name' found
        throw InternalError("No sub-mesh found with name '" + name + "'");
      }

      virtual Index get_num_entities(int dim)
      {
        switch(dim)
        {
        case 0:
          return _mesh_data->vertex_count;
        case 1:
          return _mesh_data->edge_count;
        case 2:
          return _mesh_data->tria_count + _mesh_data->quad_count;
        case 3:
          return _mesh_data->tetra_count + _mesh_data->hexa_count;
        default:
          return 0;
        }
      }

      virtual int get_num_coords()
      {
        return int(_mesh_data->coord_per_vertex);
      }

      virtual int get_vertex_stride()
      {
        return get_num_coords();
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set)
      {
        const int num_coords = get_num_coords();
        const Index num_vertices(Index(_mesh_data->coords.size()));

        for(Index i(0); i < num_vertices; ++i)
        {
          // get a reference to the corresponding vertex
          MeshStreamer::MeshDataContainer::CoordVec& vtx(_mesh_data->coords[i]);

          ASSERT(int(vtx.size()) == num_coords, "Vertex coordinate count mismatch!");

          // copy vertex coordinates
          for(int j(0); j < num_coords; ++j)
          {
            vertex_set[i][j] = Coord_(vtx[j]);
          }
        }
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
      {
        // call wrapper
        Intern::MeshStreamerIndexer<Shape_>::wrap(index_set_holder, _mesh_data->adjacencies);

        // build redundant index sets
        RedundantIndexSetBuilder<Shape_>::compute(index_set_holder);
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder)
      {
        Intern::MeshStreamerTargeter<Shape_>::wrap(target_set_holder, _mesh_data->parent_indices);
      }
    }; // class MeshStreamerFactory<ConformalSubMesh<...>>

    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct MeshStreamerIndexer
      {
        typedef std::vector< std::vector<Index> > AdjStack;
        typedef AdjStack AdjStackMatrix[4][4];
        enum
        {
          num_indices = Shape::FaceTraits<Shape_, 0>::count
        };
        typedef IndexSet<num_indices> IdxSet;

        static void wrap(IndexSetHolder<Shape_>& idx, AdjStackMatrix& adj)
        {
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType FacetType;
          MeshStreamerIndexer<FacetType>::wrap(idx, adj);
          apply(idx.template get_index_set<Shape_::dimension, 0>(), adj[0][Shape_::dimension]);
        }

        static void apply(IdxSet& idx, AdjStack& adj)
        {
          const Index num_cells(Index(adj.size()));
          for(Index i(0); i < num_cells; ++i)
          {
            std::vector<Index>& ic(adj[i]);
            ASSERT(Index(ic.size()) == Index(num_indices), "Index tuple size mismatch");

            // copy indices
            for(Index j(0); j < Index(num_indices); ++j)
            {
              idx(i,j) = ic[j];
            }
          }
        }
      };

      template<>
      struct MeshStreamerIndexer<Shape::Vertex>
      {
        template<typename T1_, typename T2_>
        static void wrap(T1_&, T2_&)
        {
          // dummy
        }
      };

      template<typename Shape_>
      struct MeshStreamerTargeter
      {
        typedef std::vector<Index> ParentIndices[4];
        static void wrap(TargetSetHolder<Shape_>& tsh, ParentIndices& pi)
        {
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension-1>::ShapeType FacetType;
          MeshStreamerTargeter<FacetType>::wrap(tsh, pi);
          apply(tsh.template get_target_set<Shape_::dimension>(), pi[Shape_::dimension]);
        }

        static void apply(TargetSet& trg, std::vector<Index>& pix)
        {
          const Index num_cells(Index(pix.size()));
          for(Index i(0); i < num_cells; ++i)
          {
            trg[i] = pix[i];
          }
        }
      };

      template<>
      struct MeshStreamerTargeter<Shape::Vertex>
      {
        typedef std::vector<Index> ParentIndices[4];
        static void wrap(TargetSetHolder<Shape::Vertex>& tsh, ParentIndices& pi)
        {
          apply(tsh.get_target_set<0>(), pi[0]);
        }

        static void apply(TargetSet& trg, std::vector<Index>& pix)
        {
          const Index num_cells(Index(pix.size()));
          for(Index i(0); i < num_cells; ++i)
          {
            trg[i] = pix[i];
          }
        }
      };
    }
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_READER_FACTORY_HPP
