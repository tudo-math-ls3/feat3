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
      template<Index dim_>
      struct MeshDataContainerUpdater;

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

        // Update _mesh_data with the new information
        Intern::MeshDataContainerUpdater<Shape_::dimension>::update_from_ish(_mesh_data, index_set_holder);
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
      /// num_entities information, might differ from the information in _mesh_data
      Index _num_entities[Shape_::dimension+1];

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
            // Parse preliminary num_entities from _mesh_data
            for(Index d(0); d <= Index(Shape_::dimension); ++d)
              _num_entities[d] = parse_num_entities(d);
            return;
          }

          // try to find a cell-set node
          MeshStreamer::CellSetNode* cell_set_node(root->find_cell_set(name));
          if(cell_set_node != nullptr)
          {
            _target_data = &cell_set_node->cell_set;
            // Parse preliminary num_entities from _mesh_data
            for(Index d(0); d <= Index(Shape_::dimension); ++d)
              _num_entities[d] = parse_num_entities(d);
            return;
          }

          // no child with 'name' found
          throw InternalError("No sub-mesh or cell-set found with name '" + name + "'");
        }

      virtual Index get_num_entities(int dim)
      {
        ASSERT( (dim >= 0) && (dim <= Shape_::dimension), "No num_entities for dimension!");
        return _num_entities[dim];
      }

        /// \brief Parses num_entities information from _target_data
      virtual Index parse_num_entities(Index dim)
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
       * \brief Wrapper class for writing back information to a MeshDataContainer
       *
       * \tparam dim_
       * Dimension for the highest dimensional shape, i.e. 3 for a Hypercube<3> mesh
       *
       * \author Jordi Paul
       *
       */
      template<Index dim_>
      struct MeshDataContainerUpdater
      {
        /**
         * \brief Updates data in a MeshDataContainer from an IndexSetHolder
         *
         * If a streamed mesh only has vertex@cell information, the other information like edge@face etc. is
         * generated by the RedundantIndexSetBuilder and can be written back to the MeshDataContainer so it can be
         * used by other methods.
         *
         * \tparam IndexSetHolderType_
         * Type for the IndexSetHolder to copy the new information from
         *
         * \param[in,out] mesh_data
         * MeshDataContainer to be updated
         *
         * \param[in] ish
         * IndexSetHolder containing the new information
         *
         *
         */
        template<typename IndexSetHolderType_>
        static void update_from_ish(MeshStreamer::MeshDataContainer* DOXY(mesh_data),
        const IndexSetHolderType_& DOXY(ish))
        {
        }
      };

      template<>
      struct MeshDataContainerUpdater<2>
      {
        template<typename IndexSetHolderType_>
        static void update_from_ish(
          MeshStreamer::MeshDataContainer* mesh_data, const IndexSetHolderType_& ish)
          {

            if(mesh_data->shape_type == mesh_data->st_tria_quad || mesh_data->shape_type == mesh_data->st_tetra_hexa)
              throw InternalError("No mixed meshes, please.");

            // The first relevant data is vertex@edge, meaning dimension 1
            if(mesh_data->edge_count == 0)
            {
              ASSERT(mesh_data->adjacencies[0][1].size() == 0, "_mesh_data contains no edges, but adjacency vector is not empty!");
              // Update vertex@edge information
              auto& vert_at_edge = ish.template get_index_set<1,0>();
              // Set number of edges
              mesh_data->edge_count = vert_at_edge.get_num_entities();

              for(Index edge(0); edge < mesh_data->edge_count; ++edge)
              {
                // idx will contain all indices for edge
                std::vector<Index> idx;
                for(Index i(0); i < vert_at_edge.num_indices; ++i)
                  idx.push_back(vert_at_edge(edge,i));

                // Add this edge to the global adjacency vector
                mesh_data->adjacencies[0][1].push_back(idx);
              }

              // Next: edge@cell
              ASSERT(mesh_data->adjacencies[1][2].size() == 0, "_mesh_data has no edges, but edge@cell adjacency information!");
              // Update edge@cell information
              auto& edge_at_cell = ish.template get_index_set<1,0>();

              // Do not set the number of shapes of dimension 2 here. Either this method is called in a 2d context
              // (then that number is correct anyway) or in a 3d context (then that number is set in the higher
              // dimensional) version

              for(Index cell(0); cell < mesh_data->tria_count + mesh_data->quad_count; ++cell)
              {
                // idx will contain all indices for edge
                std::vector<Index> idx;
                // The _index_bound for the cell to edges mapping is the number of edges
                for(Index edge(0); edge < edge_at_cell.num_indices; ++edge)
                  idx.push_back(edge_at_cell(cell,edge));
                // Add this edge to the global adjacency vector
                mesh_data->adjacencies[1][2].push_back(idx);
              }
            }
          } // MeshDataContainerUpdater<2>::update_from_ish()
      }; //MeshDataContainerUpdater<2>

      template<>
      struct MeshDataContainerUpdater<3>
      {
        template<typename IndexSetHolderType_>
        static void update_from_ish(MeshStreamer::MeshDataContainer* mesh_data, const IndexSetHolderType_& ish)
        {
          if(mesh_data->shape_type == mesh_data->st_tria_quad || mesh_data->shape_type == mesh_data->st_tetra_hexa)
            throw InternalError("No mixed meshes, please.");

          // Update 2d data first
          MeshDataContainerUpdater<2>::update_from_ish(mesh_data, ish);

          // Next: Everything face related
          // Check if the data is missing
          if(mesh_data->tria_count + mesh_data->quad_count == 0)
          {
            // vertex@face
            ASSERT(mesh_data->adjacencies[0][2].size() == 0, "_mesh_data has no faces, but vertex adjacency vector is not empty!");
            auto& vert_at_face = ish.template get_index_set<2,0>();

            // Set number of faces
            if( mesh_data->shape_type == mesh_data->st_tetra)
              mesh_data->tria_count = vert_at_face.get_num_entities();

            if( mesh_data->shape_type == mesh_data->st_hexa)
              mesh_data->quad_count = vert_at_face.get_num_entities();

            for(Index face(0); face < mesh_data->tria_count + mesh_data->quad_count; ++face)
            {
              // idx will contain all indices for face
              std::vector<Index> idx;
              for(Index i(0); i < vert_at_face.num_indices; ++i)
                idx.push_back(vert_at_face(face,i));

              // Add this edge to the global adjacency vector
              mesh_data->adjacencies[0][2].push_back(idx);
            }

            // edge@face was set by the 2d routine called earlier

            // Next: edge@cell
            ASSERT(mesh_data->adjacencies[1][3].size() == 0, "_mesh_data has no edges, but adjacency vector for faces is not empty!");

            auto& edge_at_cell= ish.template get_index_set<2,1>();

            // Update edge@cell information
            for(Index cell(0); cell < mesh_data->tetra_count + mesh_data->hexa_count; ++cell)
            {
              // idx will contain all indices for edge
              std::vector<Index> idx;
              // num_indices is the (local!) number of edges per cell
              for(Index edge(0); edge < edge_at_cell.num_indices; ++edge)
                idx.push_back(edge_at_cell(cell,edge));
              // Add this edge to the global adjacency vector
              mesh_data->adjacencies[1][3].push_back(idx);
            }

            // Last: face@cell
            ASSERT(mesh_data->adjacencies[2][3].size() == 0, "_mesh_data has no faces, but adjacency vector for faces is not empty!");
            auto& face_at_cell= ish.template get_index_set<3,2>();

            // Update face@cell information
            for(Index cell(0); cell < mesh_data->tetra_count + mesh_data->hexa_count; ++cell)
            {
              // idx will contain all indices for face
              std::vector<Index> idx;
              // num_indices is the (local!) number of faces per cell
              for(Index face(0); face < face_at_cell.num_indices; ++face)
                idx.push_back(face_at_cell(cell,face));
              // Add this face to the global adjacency vector
              mesh_data->adjacencies[2][3].push_back(idx);
            }

          } // handling faces

        } // MeshDataContainerUpdater<3>::update_from_ish()
      }; // MeshDataContainerUpdater<3>

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
