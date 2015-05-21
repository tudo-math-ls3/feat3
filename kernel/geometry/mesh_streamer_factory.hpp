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

      template<Index dim_>
      struct NumEntitiesExtractor;
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
      /// num_entities information, might differ from the information in _mesh_data
      Index _num_entities[Shape_::dimension+1];

    public:
      explicit MeshStreamerFactory(MeshStreamer& mesh_reader) :
        _mesh_reader(mesh_reader),
        _mesh_data(mesh_reader.get_mesh())
      {
        // Parse preliminary num_entities from _mesh_data
        for(Index d(0); d <= Index(Shape_::dimension); ++d)
          _num_entities[d] = parse_num_entities(d);
      }

      /// \brief Parses num_entities information from _mesh_data
      virtual Index parse_num_entities(Index dim)
      {
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
      }

      virtual Index get_num_entities(int dim)
      {
        ASSERT( (dim >= 0) && (dim <= Shape_::dimension), "No num_entities for dimension!");
        return _num_entities[dim];
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
          for(Index j(0); j < Index(num_coords_); ++j)
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

        // Set entries in num_entities. This has to be called after using the RedundantIndexSetBuilder because
        // the information from _mesh_data might not have been complete: If the mesh file did not contain edge/face
        // data, the number of entities of the corresponding dimension was zero.
        Intern::NumEntitiesExtractor<Shape_::dimension>::set_num_entities(index_set_holder, _num_entities);
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
      /// num_entities information, might differ from the information in _mesh_data
      Index _num_entities[Shape_::dimension+1];

    public:
      explicit MeshStreamerFactory(MeshStreamer& mesh_reader, String name) :
        _mesh_reader(mesh_reader),
        _name(name),
        _mesh_data(nullptr)
      {
        // Parse preliminary num_entities from _mesh_data
        for(Index d(0); d <= Shape_::dimension; ++d)
          _num_entities[d] = parse_num_entities(d);

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

      /// \brief Parses num_entities information from _mesh_data
      virtual Index parse_num_entities(Index dim)
      {
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
      }

      virtual Index get_num_entities(int dim)
      {
        ASSERT( (dim >= 0) && (dim <= Shape_::dimension), "No num_entities for dimension!");
        return _num_entities[dim];
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

        // Set entries in num_entities. This has to be called after using the RedundantIndexSetBuilder because
        // the information from _mesh_data might not have been complete: If the mesh file did not contain edge/face
        // data, the number of entities of the corresponding dimension was zero.
        Intern::NumEntitiesExtractor<Shape_::dimension>::set_num_entities(index_set_holder, _num_entities);
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
        static constexpr int num_indices = Shape::FaceTraits<Shape_, 0>::count;
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

      /**
       * \brief Helper class for extracting num_entities information from an IndexSetHolder
       *
       * \tparam dim_
       * Maximum dimension of the shape to extract said information for.
       *
       * If the mesh file read by the MeshStreamer used by the MeshStreamerFactory does not provide complete
       * information about edge@cell, edge@face and/or face@cell, the number of edges/faces is not correct in
       * the _mesh_data object. After fill_index_sets() is called, the size of these index set in the corresponding
       * IndexSetHolder gives the right number of edges/faces.
       *
       * Because the get_index_set routine uses the dimensions as template parameters, the loop has to be realised
       * by template recursion and as the end of that recursion needs partial specialisation, a helper class.
       *
       * \author Jordi Paul
       *
       */
      template<Index dim_>
      struct NumEntitiesExtractor
      {
        /**
         * \brief Computes certain entries of num_entities from an IndexSetHolder
         *
         * \tparam IndexSetHolderType_
         * Type of the IndexSetHolder
         *
         * \param[in] ish
         * IndexSetHolder to set num_entities from
         *
         * \param[out] num_entities
         * num_entities[d] = Number of objects of dimension d
         *
         */
        template<typename IndexSetHolderType_>
        static void set_num_entities(IndexSetHolderType_& ish, Index* num_entities)
        {
          // The number of entities of dimension dim_ is the length of the vertex@shape[dim_] IndexSet
          num_entities[dim_] = ish.template get_index_set<dim_,0>().get_num_entities();
          // Recurse down
          NumEntitiesExtractor<dim_-1>::set_num_entities(ish, num_entities);
        }
      };

      /**
       * \brief Full specialisation of NumEntitiesExtractor as end of the template recursion
       *
       * As num_entities[0] = number of vertices and this information must be present and correct, this stops at 1.
       *
       */
      template<>
      struct NumEntitiesExtractor<1>
      {
        template<typename IndexSetHolderType_>
        static void set_num_entities(IndexSetHolderType_& ish, Index* num_entities)
        {
          num_entities[1] = ish.template get_index_set<1,0>().get_num_entities();
        }
      };
    }
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_READER_FACTORY_HPP
