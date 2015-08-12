#pragma once
#ifndef KERNEL_GEOMETRY_MESH_STREAMER_FACTORY_HPP
#define KERNEL_GEOMETRY_MESH_STREAMER_FACTORY_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/index_calculator.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/util/mesh_streamer.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<int dim_>
      struct AttributeSetBuilder;

      template<int dim>
      struct ContainerFiller;

      template<typename Shape_>
      struct MeshStreamerIndexer;

      template<typename Shape_>
      struct MeshStreamerTargeter;

      template<int dim_>
      struct NumEntitiesExtractor;

      template<typename MeshType_>
      struct TypeConverter;
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
      /// a pointer to our mesh data container
      MeshStreamer::MeshDataContainer* _mesh_data;
      /// num_entities information, might differ from the information in _mesh_data
      Index _num_entities[Shape_::dimension+1];

    public:
      explicit MeshStreamerFactory(MeshStreamer& mesh_reader) :
        _mesh_data(mesh_reader.get_mesh())
      {
        // Parse preliminary num_entities from _mesh_data
        for(int d(0); d <= Shape_::dimension; ++d)
          _num_entities[d] = _mesh_data->num_entities[d];
      }

      explicit MeshStreamerFactory(MeshStreamer::MeshDataContainer* mesh_data) :
        _mesh_data(mesh_data)
      {
        ASSERT_(mesh_data != nullptr);
        for(int d(0); d <= Shape_::dimension; ++d)
          _num_entities[d] = _mesh_data->num_entities[d];
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

          ASSERT(vtx.size() == Index(num_coords_), "Vertex coordinate count mismatch: "+stringify(vtx.size())+ " to "+stringify(num_coords_));

          // copy vertex coordinates
          for(int j(0); j < num_coords_; ++j)
          {
            vertex_set[i][j] = Coord_(vtx[Index(j)]);
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
        FEAST::Intern::MeshDataContainerUpdater<Shape_::dimension>::update_from_ish(_mesh_data, index_set_holder);
      }

    }; // class MeshStreamerFactory<ConformalMesh<...>>

    /**
     * \brief MeshStreamerFactory implementation for MeshPart
     *
     * \author Jordi Paul
     */
    template<
      typename Shape_,
      int num_coords_,
      int stride_,
      typename Coord_>
    class MeshStreamerFactory< MeshPart<ConformalMesh<Shape_, num_coords_, stride_, Coord_> > > :
      public Factory< MeshPart<ConformalMesh<Shape_, num_coords_, stride_, Coord_> > >
    {
    public:
      /// Type of mesh the MeshPart refers to
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
      /// Type of the MeshPart
      typedef MeshPart<MeshType> MeshPartType;
      /// Vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// Attribute type
      typedef typename MeshPartType::AttributeType AttributeType;
      /// Attribute holder type
      typedef typename MeshPartType::AttributeHolderType AttributeHolderType;
      /// Index set holder type
      typedef typename MeshPartType::IndexSetHolderType IndexSetHolderType;
      /// Target set holder type
      typedef typename MeshPartType::TargetSetHolderType TargetSetHolderType;

    private:
      /// our mesh data container
      MeshStreamer::MeshDataContainer* _mesh_data;
      /// num_entities information, might differ from the information in _mesh_data
      Index _num_entities[Shape_::dimension+1];

    public:
      explicit MeshStreamerFactory(MeshStreamer& mesh_reader, const String& name) :
        _mesh_data(nullptr)
      {
        MeshStreamer::MeshNode* root(mesh_reader.get_root_mesh_node());
        ASSERT_(root != nullptr);

        // try to find a sub-mesh node
        MeshStreamer::MeshNode* sub_mesh_node(root->find_meshpart(name));
        if(sub_mesh_node != nullptr)
        {
          _mesh_data = &sub_mesh_node->mesh_data;
          // Parse preliminary num_entities from _mesh_data
          for(int d(0); d <= Shape_::dimension; ++d)
            _num_entities[d] = _mesh_data->num_entities[d];
          return;
        }

        // no child with 'name' found
        throw InternalError("No sub-mesh found with name '" + name + "'");
      }

      explicit MeshStreamerFactory(MeshStreamer::MeshDataContainer* mesh_data) :
        _mesh_data(mesh_data)
      {
        ASSERT_(mesh_data != nullptr);

        // Parse preliminary num_entities from _mesh_data
        for(int d(0); d <= Shape_::dimension; ++d)
          _num_entities[d] = _mesh_data->num_entities[d];
      }

      virtual Index get_num_entities(int dim) override
      {
        ASSERT( (dim >= 0) && (dim <= Shape_::dimension), "No num_entities for dimension!");
        return _num_entities[dim];
      }

      virtual void fill_attribute_sets(AttributeHolderType& attributes) override
      {
        Intern::AttributeSetBuilder<Shape_::dimension>::build(attributes, _mesh_data);
      }

      virtual void fill_index_sets(IndexSetHolderType*& index_set_holder) override
      {
        ASSERT(index_set_holder == nullptr, "fill_index_sets: index_set_holder != nullptr!");
        Intern::MeshStreamerIndexer<Shape_>::build(index_set_holder, _mesh_data, _num_entities);
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        Intern::MeshStreamerTargeter<Shape_>::wrap(target_set_holder, _mesh_data->parent_indices);
      }

      virtual String get_identifier() const override
      {
        return _mesh_data->name;
      }

      virtual String get_parent_identifier() const override
      {
        return _mesh_data->parent;
      }
    }; // class MeshStreamerFactory<MeshPart<...>>

    /**
     *
     */
    template<typename RootMeshType_>
    struct MeshWriter
    {
      static void write_mesh_part(MeshStreamer::MeshDataContainer& mesh_data, const MeshPart<RootMeshType_>& mesh_part)
      {
        // Fill BaseContainer stuff
        mesh_data.name = mesh_part.get_identifier();
        mesh_data.parent = mesh_part.get_parent_identifier();
        mesh_data.info = "Written by write_mesh_part";

        for(int i(0); i < 4; ++i)
          mesh_data.num_entities[i] = Index(0);

        for(int i(0); i < RootMeshType_::shape_dim+1; ++i)
          mesh_data.num_entities[i] = mesh_part.get_num_entities(i);

        Intern::TypeConverter<RootMeshType_>::fill_slices(mesh_data.slices, mesh_part);

        mesh_data.attribute_count = Index(mesh_part.get_num_attributes());
        // MeshPart does not have vertices
        mesh_data.coord_per_vertex = 0;

        // Write target sets to parent_indices
        Intern::ContainerFiller<RootMeshType_::shape_dim>::fill_parents(mesh_data.parent_indices, mesh_part.get_target_set_holder());

        // Fill MeshDataContainer stuff
        // This has to be empty and gets filled as needed by the calling function
        mesh_data.chart = "";
        mesh_data.mesh_type = mesh_data.convert_mesh_type(Intern::TypeConverter<RootMeshType_>::mesh_type());
        mesh_data.shape_type = mesh_data.convert_shape_type(Intern::TypeConverter<RootMeshType_>::shape_type());
        // Coords do not exist

        // If we have a topology, fill it
        if(mesh_part.has_topology())
          Intern::ContainerFiller<RootMeshType_::shape_dim-1>::fill_adjacencies(mesh_data.adjacencies, *(mesh_part.get_topology()));

        // If the MeshPart has any attributes, fill the containers
        if(mesh_part.get_num_attributes() > 0)
          write_attributes(mesh_data.attributes, mesh_part.get_attribute_holder());
      }

      static void write_mesh(MeshStreamer::MeshDataContainer& mesh_data, const RootMeshType_& mesh)
      {
        // Fill MeshStreamer::BaseContainer stuff
        mesh_data.name = mesh.get_identifier();
        mesh_data.parent = "";
        mesh_data.info = "Written by write_mesh";

        for(int i(0); i < 4; ++i)
          mesh_data.num_entities[i] = Index(0);

        for(int i(0); i < RootMeshType_::shape_dim+1; ++i)
          mesh_data.num_entities[i] = mesh.get_num_entities(i);

        Intern::TypeConverter<RootMeshType_>::fill_slices(mesh_data.slices, mesh);

        // A real mesh has no attributes
        mesh_data.attribute_count = Index(0);
        mesh_data.coord_per_vertex = Index(RootMeshType_::world_dim);

        // parent_indices is empty

        // Fill MeshStreamer::MeshDataContainer stuff
        mesh_data.mesh_type = mesh_data.convert_mesh_type(Intern::TypeConverter<RootMeshType_>::mesh_type());
        mesh_data.shape_type = mesh_data.convert_shape_type(Intern::TypeConverter<RootMeshType_>::shape_type());

        // Fill coords
        auto& vtx = mesh.get_vertex_set();
        for(Index i(0); i < mesh.get_num_entities(0); ++i)
        {
          std::vector<double> tmp;
          for(int d(0); d < int(mesh_data.coord_per_vertex); ++d)
            tmp.push_back(vtx[i][d]);

          mesh_data.coords.push_back(tmp);
        }

        Intern::ContainerFiller<RootMeshType_::shape_dim-1>::fill_adjacencies(mesh_data.adjacencies, mesh.get_index_set_holder());
      }

      template<typename ChartType_>
      static void write_chart(MeshStreamer::ChartContainer& chart_data, String chart_name, const ChartType_& chart)
      {
        // Fill generic MeshStreamer::ChartContainer stuff
        chart_data.info = "Written by write_chart";
        chart_data.name = chart_name;
        chart_data.type = chart.get_type();

        chart.write_data_container(chart_data);
      }

      template<typename MeshAttributeHolderType_>
      static void write_attributes(std::vector<MeshStreamer::AttributeContainer> attribute_data[4], const MeshAttributeHolderType_& mah)
      {
        Intern::ContainerFiller<MeshAttributeHolderType_::shape_dim>::fill_attributes(attribute_data, mah);
      }

    }; // class MeshWriter<RootMeshType_>

    /// \cond internal
    namespace Intern
    {
      template<int dim_>
      struct AttributeSetBuilder
      {
        template<typename AttributeHolderType_, typename MeshDataContainerType_>
        static void build(AttributeHolderType_& attribute_set_holder, const MeshDataContainerType_& mesh_data)
        {
          AttributeSetBuilder<dim_-1>::build(attribute_set_holder, mesh_data);

          auto& attributes = attribute_set_holder.template get_mesh_attributes<dim_>();
          AttributeSetBuilder<dim_>::wrap(attributes, mesh_data->attributes[dim_]);
        }

        template<typename AttributeSetType_, typename AttributeContainerType_>
        static void wrap(AttributeSetType_& attributes, const AttributeContainerType_& attribute_container)
        {
          typedef typename AttributeSetType_::value_type::CoordType DataType;

          for(auto& it:attribute_container)
          {
            const Index value_count(it.value_count);
            const int value_dim(it.value_dim);
            const String my_name(it.identifier);

            typename AttributeSetType_::value_type new_attribute(value_count, value_dim, 0, my_name);

            for(Index i(0); i < value_count; ++i)
            {
              const std::vector<double>& vals(it.values[i]);

              ASSERT(int(vals.size()) == value_dim, "Attribute value count mismatch!");

              // copy vertex coordinates
              for(int j(0); j < value_dim; ++j)
                new_attribute[i][j] = DataType(vals[size_t(j)]);
            }

            for(auto& jt:attributes)
            {
              if(jt.get_identifier() == my_name)
                throw InternalError("Attribute set already contains an attribute set with name "+my_name);
            }

            attributes.push_back(new_attribute);
          }

        } // static void
      }; // struct AttributeSetBuilder<dim_>

      template<>
      struct AttributeSetBuilder<0>
      {
        template<typename AttributeHolderType_, typename MeshDataContainerType_>
        static void build(AttributeHolderType_& attribute_set_holder, const MeshDataContainerType_& mesh_data)
        {
          auto& attributes = attribute_set_holder.template get_mesh_attributes<0>();
          AttributeSetBuilder<1>::wrap(attributes, mesh_data->attributes[0]);
        }
      };

      template<typename Shape_>
      struct MeshStreamerIndexer
      {
        typedef std::vector< std::vector<Index> > AdjStack;
        typedef AdjStack AdjStackMatrix[4][4];
        static constexpr int num_indices = Shape::FaceTraits<Shape_, 0>::count;
        typedef typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType FacetType;
        typedef IndexSet<num_indices> IdxSet;

        static void wrap(IndexSetHolder<Shape_>& idx, AdjStackMatrix& adj)
        {
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
            for(int j(0); j < num_indices; ++j)
            {
              idx(i,Index(j)) = ic[Index(j)];
            }
          }
        }

        template<typename ISH_>
        static void build(ISH_*& index_set_holder, MeshStreamer::MeshDataContainer* mesh_data, Index* num_entities)
        {
          if(index_set_holder != nullptr)
            return;

          // If _mesh_data contains adjacency information, fill index_set_holder
          if(mesh_data->adjacencies[0][Shape_::dimension].size() > 0)
          {
            // allocate a new index set holder
            index_set_holder = new ISH_(num_entities);

            // call wrapper
            wrap(*index_set_holder, mesh_data->adjacencies);

            // build redundant index sets
            RedundantIndexSetBuilder<Shape_>::compute(*index_set_holder);

            // Set entries in num_entities. This has to be called after using the RedundantIndexSetBuilder because
            // the information from _mesh_data might not have been complete: If the mesh file did not contain edge/face
            // data, the number of entities of the corresponding dimension was zero.
            Intern::NumEntitiesExtractor<Shape_::dimension>::set_num_entities(*index_set_holder, num_entities);

            // Update _mesh_data with the new information
            FEAST::Intern::MeshDataContainerUpdater<Shape_::dimension>::update_from_ish(mesh_data, *index_set_holder);
          }
          else
            MeshStreamerIndexer<FacetType>::build(index_set_holder, mesh_data, num_entities);
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
        template<typename ISH_>
        static void build(ISH_*&, MeshStreamer::MeshDataContainer*, const Index*)
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
      template<int dim_>
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

      /**
       * \brief Helper class that converts mesh types to strings
       *
       * \tparam RootMeshType_
       * Type of the root mesh.
       *
       */
      template<typename RootMeshType_>
      struct TypeConverter
      {
        /**
         * \brief Returns the ShapeType as String for the MeshStreamer
         */
        static String shape_type();

        /**
         * \brief Returns the MeshType (conformal or structured) as String for the MeshStreamer
         */
        static String mesh_type();

        /**
         * \brief Fills the slices std::vector the RootMeshType_ has slices
         *
         * \param[in,out] slices
         * MeshStreamer::BaseContainer slices.
         *
         * \param[in] mesh
         * The Geometry root mesh.
         *
         */
        static void fill_slices(std::vector<Index>& DOXY(slices), const RootMeshType_& DOXY(mesh));

        /**
         * \brief Fills the slices std::vector the MeshPart<RootMeshType_> has slices
         *
         * \param[in,out] slices
         * MeshStreamer::BaseContainer slices.
         *
         * \param[in] mesh
         * The Geometry::MeshPart<RootMeshType_>.
         *
         */
        static void fill_slices(std::vector<Index>& DOXY(slices), const MeshPart<RootMeshType_>& DOXY(slices));
      };

      template< int dim_, int num_coords_, int stride_, typename Coord_>
      struct TypeConverter<StructuredMesh<dim_, num_coords_, stride_, Coord_>>
      {
        typedef StructuredMesh<dim_, num_coords_, stride_, Coord_> MeshType;

        static String shape_type()
        {
          if(dim_ == 0)
            return "vertex";
          if(dim_ == 1)
            return "edge";
          if(dim_ == 2)
            return "quad";
          if(dim_ == 3)
            return "hexa";
          throw InternalError("Invalid dim_ parameter: "+stringify(dim_));
        }

        static String mesh_type()
        {
          return "structured";
        }

        static void fill_slices(std::vector<Index>& slices, const MeshType& mesh)
        {
          for(int d(0); d <= dim_; ++d)
            slices.push_back(mesh.get_num_slices(d));
        }

        static void fill_slices(std::vector<Index>& slices, const MeshPart<MeshType>& mesh)
        {
          for(int d(0); d <= dim_; ++d)
            slices.push_back(mesh.get_num_slices(d));
        }
      };

      template< int dim_, int num_coords_, int stride_, typename Coord_>
      struct TypeConverter<ConformalMesh<Shape::Simplex<dim_>, num_coords_, stride_, Coord_>>
      {
        typedef ConformalMesh<Shape::Simplex<dim_>, num_coords_, stride_, Coord_> MeshType;
        static String shape_type()
        {
          if(dim_ == 0)
            return "vertex";
          if(dim_ == 1)
            return "edge";
          if(dim_ == 2)
            return "tria";
          if(dim_ == 3)
            return "tetra";
          throw InternalError("Invalid dim_ parameter: "+stringify(dim_));
        }

        static String mesh_type()
        {
          return "conformal";
        }

        static void fill_slices(std::vector<Index>& DOXY(slices), const MeshType& DOXY(mesh))
        {
          // dummy as this mesh type does not have slices
        }

        static void fill_slices(std::vector<Index>& DOXY(slices), const MeshPart<MeshType>& DOXY(mesh))
        {
          // dummy as this mesh type does not have slices
        }
      };

      template< int dim_, int num_coords_, int stride_, typename Coord_>
      struct TypeConverter<ConformalMesh<Shape::Hypercube<dim_>, num_coords_, stride_, Coord_>>
      {
        typedef ConformalMesh<Shape::Hypercube<dim_>, num_coords_, stride_, Coord_> MeshType;
        static String shape_type()
        {
          if(dim_ == 0)
            return "vertex";
          if(dim_ == 1)
            return "edge";
          if(dim_ == 2)
            return "quad";
          if(dim_ == 3)
            return "hexa";
          throw InternalError("Invalid dim_ parameter: "+stringify(dim_));
        }

        static String mesh_type()
        {
          return "conformal";
        }

        static void fill_slices(std::vector<Index>& DOXY(slices), const MeshType& DOXY(mesh))
        {
          // dummy as this mesh type does not have slices
        }

        static void fill_slices(std::vector<Index>& DOXY(slices), const MeshPart<MeshType>& DOXY(mesh))
        {
          // dummy as this mesh type does not have slices
        }
      };

      /**
       * \brief Helper class that fills various containers from FEAST::MeshStreamer
       *
       * \tparam dim_
       * (Shape) Dimension the operations are carried out for
       *
       * \author Jordi Paul
       *
       */
      template<int dim_>
      struct ContainerFiller
      {
        /**
         * \brief Fills the parent_indices of a MeshStreamer::BaseContainer
         *
         * \tparam ParentSetType
         * std::vector<Index> parent_indices[4]
         *
         * \tparam TargetSetHolderType
         * Geometry::TargetSetHolder<ShapeType>, but ShapeType is unknown here.
         *
         * \param[in,out] ps
         * The BaseContainer's parent_indices
         *
         * \param[in] tsh
         * The MeshPart's TargetSetHolder
         */
        template<typename ParentSetType, typename TargetSetHolderType>
        static void fill_parents(ParentSetType& ps, const TargetSetHolderType& tsh)
        {
          // Recurse down
          ContainerFiller<dim_-1>::fill_parents(ps, tsh);

          // Add everything from our TargetSet
          auto& ts = tsh.template get_target_set<dim_>();
          for(Index shape(0); shape < ts.get_num_entities(); ++shape)
            ps[dim_].push_back(ts[shape]);

        }

        /**
         * \brief Fills the adjacencies of a MeshStreamer::MeshDataContainer
         *
         * \tparam AdjacenciesType
         * std::vector<std::vector<Index>> adjacencies
         *
         * \tparam IndexSetHolderType
         * Geometry::IndexSetHolder<ShapeType>, but ShapeType is unknown here.
         *
         * \param[in,out] adj
         * The MeshDataContainer's adjacencies
         *
         * \param[in] tsh
         * The Mesh's or MeshPart's IndexSetHolder
         */
        template<typename AdjacenciesType, typename IndexSetHolderType>
        static void fill_adjacencies(AdjacenciesType& adj, const IndexSetHolderType& ish)
        {
          ContainerFiller<dim_-1>::fill_adjacencies(adj, ish);

          // As we have dim_+1 adjacencies
          const int my_dim(dim_+1);

          // Recurse down
          auto& idx = ish.template get_index_set<my_dim,0>();

          // Add all vertex@shape (of dimension dim_+1) adjacencies
          for(Index shape(0); shape < idx.get_num_entities(); ++shape)
          {
            std::vector<Index> tmp;
            for(Index i(0); i < idx.num_indices; ++i)
              tmp.push_back(idx(shape, i));

            adj[0][my_dim].push_back(tmp);
          }
        }

        /**
         * \brief Fills the attributes of a MeshStreamer::BaseContainer
         *
         * \tparam AttributeType
         * std::vector<MeshStreamer::AttributeContainer> attributes
         *
         * \tparam MeshAttributeHolderType
         * Geometry::MeshAttributeHolder<ShapeType, DataType>, but ShapeType is unknown here.
         *
         * \param[in,out] ad
         * The BaseContainer's attributes
         *
         * \param[in] ash
         * The MeshPart's AttributeSetHolder
         */
        template<typename AttributeType, typename MeshAttributeHolderType>
        static void fill_attributes(AttributeType& ad, const MeshAttributeHolderType& ash)
        {
          ContainerFiller<dim_-1>::fill_attributes(ad, ash);

          // These are the real types, but they can be inferred, too
          //typedef Real ValueType;
          //typedef std::vector<ValueType> ValueVec;

          typedef typename std::remove_pointer<AttributeType>::type::value_type::ValueVec ValueVec;
          typedef typename ValueVec::value_type ValueType;

          auto& attribute_set = ash.template get_mesh_attributes<dim_>();

          // Add a new AttributeContainer for every Attribute in the AttributeSet to dimension dim_
          for(auto& it : attribute_set)
          {
            MeshStreamer::AttributeContainer att;
            att.identifier = it.get_identifier();
            att.value_dim = it.get_num_coords();
            att.value_count = it.get_num_vertices();

            for(Index i(0); i < att.value_count; ++i)
            {
              ValueVec tmp;

              for(int d(0); d < att.value_dim; ++d)
                tmp.push_back(ValueType(it[i][d]));

              att.values.push_back(tmp);
            }

            ad[dim_].push_back(att);
          }
        }
      };

      /**
       * \brief Explicit specialisation as end of the template recursion
       *
       */
      template<>
      struct ContainerFiller<-1>
      {
        template<typename ParentSetType, typename TargetSetHolderType>
        static void fill_parents(ParentSetType& DOXY(ps), const TargetSetHolderType& DOXY(tsh))
        {
        }

        template<typename AdjacenciesType, typename IndexSetHolderType>
        static void fill_adjacencies(AdjacenciesType& DOXY(adj), const IndexSetHolderType& DOXY(ish))
        {
        }

        template<typename AttributeType, typename MeshAttributeHolderType>
        static void fill_attributes(AttributeType& DOXY(ad), const MeshAttributeHolderType& DOXY(ash))
        {
        }
      };

    } // namespace intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_STREAMER_FACTORY_HPP
