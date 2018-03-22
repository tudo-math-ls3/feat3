#pragma once
#ifndef KERNEL_GEOMETRY_MESH_FILE_WRITER_HPP
#define KERNEL_GEOMETRY_MESH_FILE_WRITER_HPP 1

#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/partition_set.hpp>
#include <kernel/util/xml_scanner.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_, int dim_ = Shape_::dimension>
      class TopoWriteHelper
      {
      public:
        static void write_topology(std::ostream& os, const IndexSetHolder<Shape_>& ish, const String& sindent, bool skip_empty)
        {
          // recurse down
          TopoWriteHelper<Shape_, dim_-1>::write_topology(os, ish, sindent, skip_empty);

          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          const auto& index_set = ish.template get_index_set<dim_, 0>();

          // empty topology?
          // Note: this case is valid for MeshPart topologies
          if(skip_empty && (index_set.get_num_entities() == Index(0)))
            return;

          os << sindent << "<Topology dim=\"" << dim_ << "\">" << std::endl;
          for(Index i(0); i < index_set.get_num_entities(); ++i)
          {
            const auto& idx = index_set[i];
            os << sind << idx[0];
            for(int j(1); j < index_set.num_indices; ++j)
              os << ' ' << idx[j];
            os << std::endl;
          }
          os << sindent << "</Topology>" << std::endl;
        }
      };

      template<typename Shape_>
      class TopoWriteHelper<Shape_, 0>
      {
      public:
        static void write_topology(std::ostream&, const IndexSetHolder<Shape_>&, const String&, bool)
        {
          // nothing to do here
        }
      };

      template<typename Shape_, int dim_ = Shape_::dimension>
      class MappWriteHelper
      {
      public:
        static void write_mapping(std::ostream& os, const TargetSetHolder<Shape_>& tsh, const String& sindent, bool skip_empty)
        {
          // recurse down
          MappWriteHelper<Shape_, dim_-1>::write_mapping(os, tsh, sindent, skip_empty);

          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          const auto& target_set = tsh.template get_target_set<dim_>();

          // empty mapping?
          if(skip_empty && (target_set.get_num_entities() == Index(0)))
            return;

          os << sindent << "<Mapping dim=\"" << dim_ << "\">" << std::endl;
          for(Index i(0); i < target_set.get_num_entities(); ++i)
            os << sind << target_set[i] << std::endl;
          os << sindent << "</Mapping>" << std::endl;
        }
      };

      template<typename Shape_>
      class MappWriteHelper<Shape_, 0>
      {
      public:
        static void write_mapping(std::ostream& os, const TargetSetHolder<Shape_>& tsh, const String& sindent, bool skip_empty)
        {
          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          const auto& target_set = tsh.template get_target_set<0>();

          // empty mapping?
          if(skip_empty && (target_set.get_num_entities() == Index(0)))
            return;

          os << sindent << "<Mapping dim=\"0\">" << std::endl;
          for(Index i(0); i < target_set.get_num_entities(); ++i)
            os << sind << target_set[i] << std::endl;
          os << sindent << "</Mapping>" << std::endl;
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Mesh file writer class
     *
     * This class implements a writer which can export objects of type
     * MeshAtlas, RootMeshNode and PartitionSet (and all of their sub-objects)
     * into the XML-based FEAT mesh file format.
     *
     * \author Peter Zajac
     */
    class MeshFileWriter
    {
    protected:
      /// the output stream to write to
      std::ostream& _os;
      /// specifies whether indentation is to be used
      bool _indent;
      /// the current indentation string
      String _sindent;

      /// \cond internal
      template<int shape_dim_>
      static String aux_shape_string(const Shape::Hypercube<shape_dim_>&)
      {
        return "hypercube";
      }

      template<int shape_dim_>
      static String aux_shape_string(const Shape::Simplex<shape_dim_>&)
      {
        return "simplex";
      }

      template<typename Shape_, int num_coords_, typename Coord_>
      static String aux_meshtype_string(const ConformalMesh<Shape_, num_coords_, Coord_>&)
      {
        return String("conformal:") + aux_shape_string(Shape_()) + ":" + stringify(int(Shape_::dimension)) + ":" + stringify(num_coords_);
      }
      /// \endcond

      /// increases the current indent by two spaces
      void _push_indent()
      {
        if(_indent)
          _sindent.resize(_sindent.size()+2, ' ');
      }

      /// decreases the current indent by two spaces
      void _pop_indent()
      {
        if(_indent)
          _sindent.resize(_sindent.size()-2);
      }

    public:
      /**
       * \brief Creates a writer for a given output stream
       *
       * \param[in] os
       * The output stream to write to.
       *
       * \param[in] indent
       * Specifies whether indentation is to be used.
       * Can be set to \c false to decrease the final file size.
       */
      explicit MeshFileWriter(std::ostream& os, bool indent = true) :
        _os(os),
        _indent(indent)
      {
        if(_indent)
          _sindent.reserve(32);
      }

      /// virtual destructor
      virtual ~MeshFileWriter()
      {
      }

      /**
       * \brief Writes a chart into the file.
       *
       * \param[in] chart
       * The chart to be exported.
       *
       * \param[in] name
       * The name of the chart.
       */
      template<typename RootMesh_>
      void write_chart(const Atlas::ChartBase<RootMesh_>& chart, const String& name)
      {
        String sind(_sindent);
        if(_indent)
          sind.append("  ");

        _os << _sindent << "<Chart name=\"" << name << "\">" << std::endl;
        chart.write(_os, sind);
        _os << _sindent << "</Chart>" << std::endl;
      }

      /**
       * \brief Writes all charts of an atlas into the file
       *
       * \param[in] mesh_atlas
       * The atlas whose charts are to be exported.
       */
      template<typename RootMesh_>
      void write_atlas(const MeshAtlas<RootMesh_>& mesh_atlas)
      {
        // get all charts
        std::deque<String> names = mesh_atlas.get_chart_names();

        // loop over all charts
        for(auto it = names.begin(); it != names.end(); ++it)
        {
          const Atlas::ChartBase<RootMesh_>* chart = mesh_atlas.find_mesh_chart(*it);
          if(chart == nullptr)
            throw InternalError(String("Chart '") + (*it) + "' not found in atlas");

          write_chart(*chart, *it);
        }
      }

      /**
       * \brief Writes a (root) mesh into the file
       *
       * \param[in] mesh
       * The (root) mesh to be exported.
       */
      template<typename Shape_, int num_coords_, typename Coord_>
      void write_mesh(const ConformalMesh<Shape_, num_coords_, Coord_>& mesh)
      {
        typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
        _os << _sindent << "<Mesh type=\"" << aux_meshtype_string(mesh) << "\"";
        _os << " size=\"" << mesh.get_num_entities(0);
        for(int i(1); i <= MeshType::shape_dim; ++i)
           _os << " " << mesh.get_num_entities(i);
        _os << "\">" << std::endl;
        _push_indent();

        _write_vertex_set(mesh.get_vertex_set());
        Intern::TopoWriteHelper<Shape_>::write_topology(_os, mesh.get_index_set_holder(), _sindent, false);

        _pop_indent();
        _os << _sindent << "</Mesh>" << std::endl;
      }

      /**
       * \brief Writes a mesh-part into the file
       *
       * \param[in] meshpart
       * The mesh-part to be exported.
       *
       * \param[in] parent_name
       * The name of the parent of the mesh-part.
       *
       * \param[in] part_name
       * The name of the mesh-part
       *
       * \param[in] chart_name
       * The name of the chart of the mesh-part.
       */
      template<typename Mesh_>
      void write_meshpart(const MeshPart<Mesh_>& meshpart, const String& parent_name, const String& part_name, const String& chart_name)
      {
        typedef typename Mesh_::ShapeType ShapeType;

        _os << _sindent << "<MeshPart";
        _os << " name=\"" << part_name << "\"";
        _os << " parent=\"" << parent_name << "\"";
        if(!chart_name.empty())
          _os << " chart=\"" << chart_name << "\"";
        _os << " topology=\"" << (meshpart.has_topology() ? "full" : "none") << "\"";
        _os << " size=\"" << meshpart.get_num_entities(0);
        for(int i(1); i <= ShapeType::dimension; ++i)
           _os << " " << meshpart.get_num_entities(i);
        _os << "\">" << std::endl;
        _push_indent();

        // write mapping
        Intern::MappWriteHelper<ShapeType>::write_mapping(_os, meshpart.get_target_set_holder(), _sindent, true);

        // write topology
        if(meshpart.has_topology())
          Intern::TopoWriteHelper<ShapeType>::write_topology(_os, *meshpart.get_topology(), _sindent, true);

        // write attributes
        const auto& attrs = meshpart.get_mesh_attributes();
        for(auto it = attrs.begin(); it != attrs.end(); ++it)
          _write_attribute(*(it->second), it->first);

        _pop_indent();
        _os << _sindent << "</MeshPart>" << std::endl;
      }

      /**
       * \brief Writes a partition to the file
       *
       * \param[in] partition
       * The partition to be exported-
       */
      void write_partition(const Partition& partition)
      {
        _os << _sindent << "<Partition";
        if(!partition.get_name().empty())
          _os << " name=\"" << partition.get_name() << "\"";
        _os << " priority=\"" << partition.get_priority() << "\"";
        _os << " level=\"" << partition.get_level() << "\"";
        _os << " size=\"" << partition.get_num_patches() << ' ' << partition.get_num_elements() << "\"";
        _os << ">" << std::endl;
        _push_indent();

        // loop over all patches
        const Adjacency::Graph& graph = partition.get_patches();
        for(Index i(0); i < graph.get_num_nodes_domain(); ++i)
        {
          _os << _sindent << "<Patch rank=\"" << i << "\" size=\"" << graph.degree(i) << "\">" << std::endl;
          _push_indent();
          for(auto it = graph.image_begin(i); it != graph.image_end(i); ++it)
            _os << _sindent << (*it) << std::endl;
          _pop_indent();
          _os << _sindent << "</Patch>" << std::endl;
        }

        _pop_indent();
        _os << _sindent << "</Partition>" << std::endl;
      }

      /**
       * \brief Writes all partitions of a partition set to the file
       *
       * \param[in] part_set
       * The partition set whose partitions are to be exported.
       */
      void write_partition_set(const PartitionSet& part_set)
      {
        for(const auto& p : part_set.get_partitions())
          write_partition(p);
      }

      /**
       * \brief Writes a full domain to the file
       *
       * \param[in] mesh_node
       * A root mesh node whose mesh and child mesh-parts are to be exported.
       * May be \c nullptr if no mesh or mesh-parts are to be exported.
       *
       * \param[in] mesh_atlas
       * An mesh atlas whose charts are to be exported.
       * May be \c nullptr if no charts are to be exported.
       *
       * \param[in] part_set
       * A partition set whose partitions are to be exported.
       * May be \c nullptr if no partitions are to be exported.
       *
       * \param[in] skip_internal_meshparts
       * Specifies whether internal mesh-parts (e.g. comm halos or partition patches) are
       * to be exported or not. Defaults to \c true, i.e. internal mesh parts are not exported by default.
       */
      template<typename RootMesh_>
      void write(
        const RootMeshNode<RootMesh_>* mesh_node,
        const MeshAtlas<RootMesh_>* mesh_atlas,
        const PartitionSet* part_set,
        bool skip_internal_meshparts = true)
      {
        const RootMesh_* root_mesh(nullptr);
        if(mesh_node != nullptr)
          root_mesh = mesh_node->get_mesh();

        // print root markup
        _os << "<FeatMeshFile version=\"1\"";
        if(root_mesh != nullptr)
          _os << " mesh=\"" << aux_meshtype_string(*root_mesh) << "\"";
        _os << ">" << std::endl;

        // increase indent
        _push_indent();

        // write mesh atlas
        if(mesh_atlas != nullptr)
        {
          // write atlas
          write_atlas(*mesh_atlas);
        }

        // write mesh node
        if(mesh_node != nullptr)
        {
          // write mesh
          if(root_mesh != nullptr)
            write_mesh(*root_mesh);

          // write meshparts
          std::deque<String> part_names = mesh_node->get_mesh_part_names();
          for(auto it = part_names.begin(); it != part_names.end(); ++it)
          {
            // skip mesh parts starting with an underscore;
            // these are reserved for internal use and should not be exported
            // unless explicitly specified by the caller
            if(skip_internal_meshparts && it->starts_with('_'))
              continue;

            const MeshPart<RootMesh_>* meshpart = mesh_node->find_mesh_part(*it);
            String chart_name = mesh_node->find_mesh_part_chart_name(*it);
            write_meshpart(*meshpart, "root", *it, chart_name);
          }
        }

        // write partitions
        if(part_set != nullptr)
        {
          write_partition_set(*part_set);
        }

        _pop_indent();

        // write terminator
        _os << "</FeatMeshFile>" << std::endl;
      }

    protected:
      /**
       * \brief Writes a vertex-set of a mesh into the file
       *
       * \param[in] vertex_set
       */
      template<int num_coords_, typename Coord_>
      void _write_vertex_set(const VertexSet<num_coords_, Coord_>& vertex_set)
      {
        _os << _sindent << "<Vertices>" << std::endl;
        _push_indent();
        for(Index i(0); i < vertex_set.get_num_vertices(); ++i)
        {
          const auto& v = vertex_set[i];
          _os << _sindent << v[0];
          for(int j(1); j < num_coords_; ++j)
            _os << ' ' << v[j];
          _os << std::endl;
        }
        _pop_indent();
        _os << _sindent << "</Vertices>" << std::endl;
      }

      /**
       * \brief Writes an attribte of a mesh part into the file
       *
       * \param[in] attr
       * The attribute that is to be exported.
       *
       * \param[in] name
       * The name of the attribute.
       */
      template<typename Data_>
      void _write_attribute(const AttributeSet<Data_>& attr, const String& name)
      {
        _os << _sindent << "<Attribute";
        _os << " name=\"" << name << "\"";
        _os << " dim=\"" << attr.get_dimension() << "\"";
        _os << ">" << std::endl;
        _push_indent();
        for(Index i(0); i < attr.get_num_values(); ++i)
        {
          _os << _sindent << attr(i, 0);
          for(int j(1); j < attr.get_dimension(); ++j)
            _os << ' ' << attr(i, j);
          _os << std::endl;
        }
        _pop_indent();
        _os << _sindent << "</Attribute>" << std::endl;
      }
    }; // class MeshFileWriter
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_MESH_FILE_WRITER_HPP
