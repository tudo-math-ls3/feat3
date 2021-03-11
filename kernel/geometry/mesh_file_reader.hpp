// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_MESH_FILE_READER_HPP
#define KERNEL_GEOMETRY_MESH_FILE_READER_HPP 1

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/atlas/bezier.hpp>
#include <kernel/geometry/atlas/circle.hpp>
#include <kernel/geometry/atlas/extrude.hpp>
#include <kernel/geometry/atlas/surface_mesh.hpp>
#include <kernel/geometry/atlas/sphere.hpp>
#include <kernel/adjacency/dynamic_graph.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/xml_scanner.hpp>
#include <kernel/util/dist_file_io.hpp>

#include <deque>
#include <vector>
#include <memory>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief MeshNodeLinker Error class
     *
     * This exception is thrown by the MeshNodeLinker class if anything goes wrong.
     *
     * \author Peter Zajac
     */
    class MeshNodeLinkerError :
      public Exception
    {
    public:
      explicit MeshNodeLinkerError(const String& msg) :
        Exception(msg)
      {
      }
    }; // class MeshNodeLinkerError

    /**
     * \brief MeshNode liner class template
     *
     * This linker class is used by the MeshFileReader for various post-parse operations,
     * such as deducting meshpart topologies or linking meshparts to their corresponding charts.
     *
     * \author Peter Zajac
     */
    template<typename RootMesh_>
    class MeshNodeLinker
    {
    protected:
      /// our mesh node
      RootMeshNode<RootMesh_>& _mesh_node;
      /// our atlas
      MeshAtlas<RootMesh_>& _atlas;

      // link mesh-parts to charts
      std::deque<std::pair<String,String>> _meshpart_to_chart;
      // deduct mesh-part topologies
      std::deque<String> _meshpart_deduct_topo;

    public:
      /**
       * \brief Constructor
       *
       * \param[in,out] mesh_node
       * The mesh node that the linker should operate on.
       *
       * \param[in,out] atlas
       * The atlas that contains the charts.
       */
      explicit MeshNodeLinker(RootMeshNode<RootMesh_>& mesh_node, MeshAtlas<RootMesh_>& atlas) :
        _mesh_node(mesh_node),
        _atlas(atlas)
      {
      }

      /**
       * \brief Adds a task to link a meshpart to a chart.
       *
       * \param[in] meshpart
       * The name of the meshpart that is to be linked to a chart.
       *
       * \param[in] chart
       * The name of the chart that the meshpart is to be linked to.
       */
      void meshpart_link_to_chart(const String& meshpart, const String& chart)
      {
        _meshpart_to_chart.emplace_back(std::make_pair(meshpart, chart));
      }

      /**
       * \brief Adds a task to deduct a meshpart topology from the root mesh.
       *
       * \param[in] meshpart
       * The name of the meshpart whose topology is to be deducted from the root mesh.
       */
      void meshpart_deduct_topology(const String& meshpart)
      {
        _meshpart_deduct_topo.emplace_back(meshpart);
      }

      /**
       * \brief Executes the linker.
       *
       * This function executes all the tasks that have been added to the linker.
       */
      void execute()
      {
        // link meshparts to charts
        while(!_meshpart_to_chart.empty())
        {
          String mpart_name = _meshpart_to_chart.front().first;
          String chart_name = _meshpart_to_chart.front().second;

          // try to find the chart
          const Atlas::ChartBase<RootMesh_>* chart = _atlas.find_mesh_chart(chart_name);
          if(chart == nullptr)
          {
            String msg = String("Chart '") + chart_name + "' not found for meshpart '" + mpart_name + "'";
            throw MeshNodeLinkerError(msg);
          }

          // set the chart
          if(!_mesh_node.set_mesh_part_chart(mpart_name, chart_name, chart))
          {
            String msg = String("meshpart '" + mpart_name + "' not found");
            throw MeshNodeLinkerError(msg);
          }

          // done
          _meshpart_to_chart.pop_front();
        }

        // ger the root mesh (if it exists)
        const RootMesh_* root_mesh = _mesh_node.get_mesh();

        // deduct meshpart topologies
        while(!_meshpart_deduct_topo.empty())
        {
          String mpart_name = _meshpart_deduct_topo.front();
          MeshPart<RootMesh_>* mesh_part = _mesh_node.find_mesh_part(mpart_name);

          if(mesh_part == nullptr)
          {
            String msg = String("meshpart '" + mpart_name + "' not found");
            throw MeshNodeLinkerError(msg);
          }

          // make sure we have the root mesh
          if(root_mesh == nullptr)
          {
            String msg = String("Cannot deduct topology for meshpart '") + mpart_name + "'; no root mesh found";
            throw MeshNodeLinkerError(msg);
          }

          // Create the MeshPart's topology from the parent mesh's topology if told so
          mesh_part->deduct_topology(*(root_mesh->get_topology()));

          // done
          _meshpart_deduct_topo.pop_front();
        }
      }
    }; // class MeshNodeLinker

    // Note: All the basic Parser classes are declared as internal.
    /// \cond internal

    // Forward declaration
    template<typename RootMesh_, int world_dim>
    struct DimensionalChartHelper;

    template<typename RootMesh_>
    class ChartParser :
      public Xml::MarkupParser
    {
    protected:
      MeshAtlas<RootMesh_>& _atlas;
      String _chart_name;
      Atlas::ChartBase<RootMesh_>* _chart;

    public:
      explicit ChartParser(MeshAtlas<RootMesh_>& atlas) :
        _atlas(atlas), _chart_name(), _chart(nullptr)
      {
      }

      virtual ~ChartParser()
      {
        if(_chart != nullptr)
          delete _chart;
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("name", true); // mandatory name
        return true;
      }

      virtual void create(
        int iline, const String& sline, const String&,
        const std::map<String, String>& attrs, bool closed) override
      {
        if(closed)
          throw Xml::GrammarError(iline, sline, "Invalid closed markup");

        // get the chart name
        _chart_name = attrs.find("name")->second.trim();
        if(_chart_name.empty())
          throw Xml::GrammarError(iline, sline, "Empty chart name");

        // make sure that we don't have that chart yet
        XASSERTM(_atlas.find_mesh_chart(_chart_name) == nullptr,
          String("Chart '") + _chart_name + "' already exists in atlas");

        // okay, that's all there is to check
      }

      virtual void close(int iline, const String& sline) override
      {
        // make sure that we have a chart here
        if(_chart == nullptr)
          throw Xml::GrammarError(iline, sline, "Invalid empty chart");

        // insert chart into atlas
        _atlas.add_mesh_chart(_chart_name, _chart);
        _chart = nullptr;
      }

      virtual bool content(int, const String&) override
      {
        return false; // invalid content
      }

      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String& name) override
      {
        return DimensionalChartHelper<RootMesh_, RootMesh_::world_dim>::markup(name, _chart);
      }

    }; // class ChartParser<...>

    template<int num_coords_, typename Coord_>
    class VerticesParser :
      public Xml::MarkupParser
    {
    public:
      typedef VertexSet<num_coords_, Coord_> VertexSetType;

    protected:
      VertexSetType& _vertex_set;
      Index _read;

    public:
      explicit VerticesParser(VertexSetType& vertex_set) :
        _vertex_set(vertex_set),
        _read(0)
      {
      }

      virtual bool attribs(std::map<String,bool>&) const override
      {
        return true;
      }

      virtual void create(int iline, const String& sline, const String&, const std::map<String, String>&, bool closed) override
      {
        if(closed)
          throw Xml::GrammarError(iline, sline, "Invalid closed markup");
      }

      virtual void close(int iline, const String& sline) override
      {
        // ensure that we have read all vertices
        if(_read < _vertex_set.get_num_vertices())
          throw Xml::GrammarError(iline, sline, "Invalid terminator; expected point");
      }

      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
      {
        // no children allowed
        return nullptr;
      }

      virtual bool content(int iline, const String& sline) override
      {
        // make sure that we do not read more points than expected
        if(_read >= _vertex_set.get_num_vertices())
          throw Xml::ContentError(iline, sline, "Invalid content; expected terminator");

        // split line by whitespaces
        std::deque<String> scoords = sline.split_by_whitespaces();

        // check size
        if(scoords.size() != std::size_t(num_coords_))
          throw Xml::ContentError(iline, sline, "Invalid number of coordinates");

        // get our current vertex
        auto& vtx = _vertex_set[_read];

        // try to parse all coords
        for(int i(0); i < num_coords_; ++i)
        {
          if(!scoords.at(std::size_t(i)).parse(vtx[i]))
            throw Xml::ContentError(iline, sline, "Failed to parse vertex coordinate");
        }

        // okay, another point done
        ++_read;

        return true;
      }
    };

    namespace Intern
    {
      template<typename Shape_, int dim_ = Shape_::dimension>
      struct TopoParseHelper
      {
        template<typename TopoParser_>
        static void init_topo(TopoParser_& tp, IndexSetHolder<Shape_>& ish, int dim)
        {
          if(dim_ == dim)
            tp.init_idx_set(ish.template get_index_set<dim_, 0>(), dim_);
          else
            TopoParseHelper<Shape_, dim_-1>::init_topo(tp, ish, dim);
        }
      };

      template<typename Shape_>
      struct TopoParseHelper<Shape_, 0>
      {
        template<typename TopoParser_>
        static void init_topo(TopoParser_&, IndexSetHolder<Shape_>&, int)
        {
          XABORTM("Thou shall not arrive here");
        }
      };

      template<typename Shape_, int dim_ = Shape_::dimension>
      struct MappParseHelper
      {
        template<typename MappParser_>
        static void init_mapp(MappParser_& mp, TargetSetHolder<Shape_>& tsh, int dim)
        {
          if(dim_ == dim)
            mp.init_trg_set(tsh.template get_target_set<dim_>(), dim_);
          else
            MappParseHelper<Shape_, dim_-1>::init_mapp(mp, tsh, dim);
        }
      };

      template<typename Shape_>
      struct MappParseHelper<Shape_, 0>
      {
        template<typename MappParser_>
        static void init_mapp(MappParser_& mp, TargetSetHolder<Shape_>& tsh, int)
        {
          mp.init_trg_set(tsh.template get_target_set<0>(), 0);
        }
      };
    } // namespace Intern

    template<typename Shape_>
    class TopologyParser :
      public Xml::MarkupParser
    {
    public:
      typedef IndexSetHolder<Shape_> TopologyType;

    protected:
      TopologyType& _topology;
      std::deque<int>& _have_topology;
      Index* _indices;
      Index _my_dim;
      Index _num_idx;
      Index _bound;
      Index _count;
      Index _read;

    public:
      template<int num_idx_>
      void init_idx_set(IndexSet<num_idx_>& idx_set, int dim)
      {
        _indices = reinterpret_cast<Index*>(idx_set.get_indices());
        _my_dim = Index(dim);
        _num_idx = Index(num_idx_);
        _bound = idx_set.get_index_bound();
        _count = idx_set.get_num_entities();
        _read = Index(0);
      }

    public:
      explicit TopologyParser(TopologyType& topology, std::deque<int>& have_topology) :
        _topology(topology),
        _have_topology(have_topology),
        _indices(nullptr),
        _my_dim(0),
        _num_idx(0),
        _bound(0),
        _count(0),
        _read(0)
      {
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("dim", true);
        return true;
      }

      virtual void create(int iline, const String& sline, const String&, const std::map<String, String>& attrs, bool closed) override
      {
        if(closed)
          throw Xml::GrammarError(iline, sline, "Invalid closed markup");

        // parse dimension
        if(!attrs.find("dim")->second.parse(_my_dim))
          throw Xml::ContentError(iline, sline, "Failed to parse 'dim' attribute");

        // make sure that the dimension is not out-of-bounds
        if((_my_dim > Index(_have_topology.size())) || (_my_dim == Index(0)))
          throw Xml::ContentError(iline, sline, "Invalid topology dimension");

        // make sure that we did not have that topology yet
        if(_have_topology.at(_my_dim-1) != 0)
          throw Xml::ContentError(iline, sline, "Multiple topology for dimension " + stringify(_my_dim));
        _have_topology.at(_my_dim-1) = 1;

        // initialize the corresponding topology:
        Intern::TopoParseHelper<Shape_>::init_topo(*this, _topology, int(_my_dim));
      }

      virtual void close(int iline, const String& sline) override
      {
        // ensure that we have read all index tuples
        if(_read < _count)
          throw Xml::GrammarError(iline, sline, "Invalid terminator; expected index tuple");
      }

      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
      {
        // no children allowed
        return nullptr;
      }

      virtual bool content(int iline, const String& sline) override
      {
        // make sure that we do not read more points than expected
        if(_read >= _count)
          throw Xml::ContentError(iline, sline, "Invalid content; expected terminator");

        // split line by whitespaces
        std::deque<String> sidx = sline.split_by_whitespaces();

        // check size
        if(sidx.size() != std::size_t(_num_idx))
          throw Xml::ContentError(iline, sline, "Invalid number of indices");

        // parse
        for(Index i(0); i < _num_idx; ++i)
        {
          Index& idx = _indices[_read*_num_idx+i];

          // try to parse
          if(!sidx.at(i).parse(idx))
            throw Xml::ContentError(iline, sline, "Failed to parse index");

          // check for out-of-bounds
          if(idx >= _bound)
            throw Xml::ContentError(iline, sline, "Index out of bounds");
        }

        // okay, another one processed
        ++_read;

        return true;
      }
    };

    template<typename Shape_>
    class MappingParser :
      public Xml::MarkupParser
    {
    public:
      typedef TargetSetHolder<Shape_> MappingType;

    protected:
      MappingType& _mapping;
      std::deque<int>& _have_mapping;
      Index* _indices;
      Index _my_dim;
      Index _count;
      Index _read;

    public:
      void init_trg_set(TargetSet& trg_set, int dim)
      {
        _indices = trg_set.get_indices();
        _my_dim = Index(dim);
        _count = trg_set.get_num_entities();
      }

    public:
      explicit MappingParser(MappingType& mapping, std::deque<int>& have_mapping) :
        _mapping(mapping),
        _have_mapping(have_mapping),
        _indices(nullptr),
        _my_dim(0),
        _count(0),
        _read(0)
      {
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("dim", true);
        return true;
      }

      virtual void create(int iline, const String& sline, const String&, const std::map<String, String>& attrs, bool closed) override
      {
        if(closed)
          throw Xml::GrammarError(iline, sline, "Invalid closed markup");

        // parse dimension
        if(!attrs.find("dim")->second.parse(_my_dim))
          throw Xml::ContentError(iline, sline, "Failed to parse 'dim' attribute");

        // make sure that the dimension is not out-of-bounds
        if((_my_dim > Index(_have_mapping.size())))
          throw Xml::ContentError(iline, sline, "Invalid mapping dimension");

        // make sure that we did not have that topology yet
        if(_have_mapping.at(_my_dim) != 0)
          throw Xml::ContentError(iline, sline, "Multiple mapping for dimension " + stringify(_my_dim));
        _have_mapping.at(_my_dim) = 1;

        // initialize the corresponding mapping:
        Intern::MappParseHelper<Shape_>::init_mapp(*this, _mapping, int(_my_dim));
      }

      virtual void close(int iline, const String& sline) override
      {
        // ensure that we have read all index tuples
        if(_read < _count)
          throw Xml::GrammarError(iline, sline, "Invalid terminator; expected index");
      }

      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
      {
        // no children allowed
        return nullptr;
      }

      virtual bool content(int iline, const String& sline) override
      {
        // make sure that we do not read more points than expected
        if(_read >= _count)
          throw Xml::ContentError(iline, sline, "Invalid content; expected terminator");

        // try to parse index
        if(!sline.parse(_indices[_read]))
          throw Xml::ContentError(iline, sline, "Failed to parse index");

        // okay, another one processed
        ++_read;

        return true;
      }
    };

    template<typename MeshPart_, typename DataType_ = typename MeshPart_::AttributeDataType>
    class AttributeParser :
      public Xml::MarkupParser
    {
    public:
      typedef AttributeSet<DataType_> AttribType;

    protected:
      MeshPart_& _mesh_part;
      AttribType* _attrib;
      Index _my_dim;
      Index _count;
      Index _read;
      String _my_name;

    public:
      explicit AttributeParser(MeshPart_& mesh_part) :
        _mesh_part(mesh_part),
        _attrib(nullptr),
        _my_dim(0),
        _count(0),
        _read(0)
      {
      }

      virtual ~AttributeParser()
      {
        // Because adding the AttributeSet to a MeshPart passes ownership, we just need to nullify the pointer here
        if(_attrib != nullptr)
          _attrib = nullptr;
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("dim", true);
        attrs.emplace("name", true);
        return true;
      }

      virtual void create(int iline, const String& sline, const String&, const std::map<String, String>& attrs, bool closed) override
      {
        if(closed)
          throw Xml::GrammarError(iline, sline, "Invalid closed markup");

        // parse dimension
        if(!attrs.find("dim")->second.parse(_my_dim))
          throw Xml::ContentError(iline, sline, "Failed to parse 'dim' attribute");

        // make sure that the dimension is valid
        if(_my_dim == Index(0))
          throw Xml::ContentError(iline, sline, "Invalid attribute dimension");

        // fetch name
        _my_name = attrs.find("name")->second;

        // create mesh attribute
        _count = _mesh_part.get_num_entities(0);
        _attrib = new AttribType(_count, int(_my_dim));
      }

      virtual void close(int iline, const String& sline) override
      {
        // ensure that we have read all index tuples
        if(_read < _count)
          throw Xml::GrammarError(iline, sline, "Invalid terminator; expected index");

        // okay, add attribute to mesh part
        _mesh_part.add_attribute(_attrib, _my_name);
      }

      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
      {
        // no children allowed
        return nullptr;
      }

      virtual bool content(int iline, const String& sline) override
      {
        // make sure that we do not read more points than expected
        if(_read >= _count)
          throw Xml::ContentError(iline, sline, "Invalid content; expected terminator");

        // split line by whitespaces
        std::deque<String> svals = sline.split_by_whitespaces();

        // check size
        if(svals.size() != std::size_t(_my_dim))
          throw Xml::ContentError(iline, sline, "Invalid number of values");

        // try to parse all coords
        for(int i(0); i < int(_my_dim); ++i)
        {
          if(!svals.at(std::size_t(i)).parse(_attrib->operator()(_read,i)))
            throw Xml::ContentError(iline, sline, "Failed to parse value");
        }

        // okay, another one processed
        ++_read;

        return true;
      }
    };

    template<typename Mesh_>
    class MeshParser;

    template<typename Shape_, int num_coords_, typename Coord_>
    class MeshParser<ConformalMesh<Shape_, num_coords_, Coord_>> :
      public Xml::MarkupParser
    {
    public:
      typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;

    protected:
      RootMeshNode<MeshType>& _root_node;
      MeshType* _mesh;
      bool _have_verts;
      std::deque<int> _have_topology;
      std::vector<Index> _sizes;

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

    public:
      explicit MeshParser(RootMeshNode<MeshType>& root_node) :
        _root_node(root_node),
        _mesh(nullptr),
        _have_verts(false)
      {
      }

      virtual ~MeshParser()
      {
        // Note: if the mesh was parsed successfully, then the
        // pointer to the mesh object is passed to the root mesh node
        // and the "_mesh" pointer is set back to nullptr.
        if(_mesh)
          delete _mesh;
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("type", true);
        attrs.emplace("size", true);
        return true;
      }

      virtual void create(int iline, const String& sline, const String&, const std::map<String, String>& attrs, bool closed) override
      {
        if(closed)
          throw Xml::GrammarError(iline, sline, "Invalid closed markup");

        // fetch mesh type substrings
        std::deque<String> stype = attrs.find("type")->second.split_by_string(":");

        // validate mesh type
        if(stype.size() != std::size_t(4))
          throw Xml::ContentError(iline, sline, "Invalid mesh type attribue");
        if(stype.at(0) != "conformal")
          throw Xml::ContentError(iline, sline, "Invalid mesh type; expected 'conformal'");
        String sshape = aux_shape_string(Shape_());
        if(stype.at(1) != sshape)
          throw Xml::ContentError(iline, sline, "Invalid mesh shape type; expected '" + sshape + "'");
        int ishape_dim(0), iworld_dim(0);
        if(!stype.at(2).parse(ishape_dim))
          throw Xml::ContentError(iline, sline, "Failed to parse mesh shape dimension");
        if(ishape_dim != Shape_::dimension)
          throw Xml::ContentError(iline, sline, "Invalid mesh shape dimension");
        if(!stype.at(3).parse(iworld_dim))
          throw Xml::ContentError(iline, sline, "Failed to parse mesh world dimension");
        if(iworld_dim != num_coords_)
          throw Xml::ContentError(iline, sline, "Invalid mesh world dimension");

        // okay, the mesh type seems fine

        // split sizes
        std::deque<String> ssize = attrs.find("size")->second.split_by_whitespaces();
        if(ssize.size() != std::size_t(Shape_::dimension+1))
          throw Xml::ContentError(iline, sline, "Invalid mesh size count");

        // parse sizes
        _sizes.resize(ssize.size());
        for(std::size_t i(0); i <= std::size_t(Shape_::dimension); ++i)
        {
          if(!ssize.at(i).parse(_sizes.at(i)))
            throw Xml::ContentError(iline, sline, "Failed to parse mesh size");
        }

        // try to create the mesh
        _mesh = new MeshType(_sizes.data());

        // initialize
        _have_verts = false;
        _have_topology.resize(std::size_t(Shape_::dimension), 0);
      }

      virtual void close(int iline, const String& sline) override
      {
        // do we have coords?
        if(!_have_verts)
          throw Xml::GrammarError(iline, sline, "Mesh has no Vertices");

        // make sure we have the whole topology
        for(std::size_t i(0); i < _have_topology.size(); ++i)
        {
          if(_have_topology.at(i) == 0)
            throw Xml::GrammarError(iline, sline, "Missing topology for dimension " + stringify(i+1));
        }

        // okay, the mesh is complete

        // compute remaining topology now
        RedundantIndexSetBuilder<Shape_>::compute(*_mesh->get_topology());

        // add to mesh node
        _root_node.set_mesh(_mesh);

        // clean up our pointer, so that the object is not deleted
        // upon destruction of this parser object
        _mesh = nullptr;
      }

      virtual std::shared_ptr<Xml::MarkupParser> markup(int iline, const String& sline, const String& name) override
      {
        // vertices?
        if(name == "Vertices")
        {
          if(_have_verts)
            throw Xml::GrammarError(iline, sline, "More than one Vertices block in mesh");

          // okay, create a vertex parser
          _have_verts = true;
          return std::make_shared<VerticesParser<num_coords_, Coord_>>(_mesh->get_vertex_set());
        }

        // topology?
        if(name == "Topology")
        {
          // create a topology parser
          return std::make_shared<TopologyParser<Shape_>>(*_mesh->get_topology(), _have_topology);
        }

        // invalid
        return nullptr;
      }

      virtual bool content(int, const String&) override
      {
        return false; // no content allowed
      }
    };

    template<typename Mesh_>
    class MeshPartParser;

    template<typename Shape_, int num_coords_, typename Coord_>
    class MeshPartParser<ConformalMesh<Shape_, num_coords_, Coord_>> :
      public Xml::MarkupParser
    {
    public:
      typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
      typedef MeshPart<MeshType> MeshPartType;
      typedef Atlas::ChartBase<MeshType> ChartType;
      enum class TopoType
      {
        none,
        full,
        parent
      };

    protected:
      RootMeshNode<MeshType>& _root_node;
      MeshNodeLinker<MeshType>& _linker;
      MeshPartType* _mesh_part;
      ChartType* _chart;
      TopoType _topo_type;
      String _name, _parent, _chart_name;
      std::deque<int> _have_topology;
      std::deque<int> _have_mapping;
      std::vector<Index> _sizes;

    public:
      explicit MeshPartParser(RootMeshNode<MeshType>& root_node, MeshNodeLinker<MeshType>& linker) :
        _root_node(root_node),
        _linker(linker),
        _mesh_part(nullptr),
        _chart(nullptr),
        _topo_type(TopoType::none)
      {
      }

      virtual ~MeshPartParser()
      {
        if(_mesh_part != nullptr)
          delete _mesh_part;
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("name", true);
        attrs.emplace("parent", true);
        attrs.emplace("size", true);
        attrs.emplace("topology", true);
        attrs.emplace("chart", false);
        return true;
      }

      virtual void create(int iline, const String& sline, const String&, const std::map<String, String>& attrs, bool closed) override
      {
        if(closed)
          throw Xml::GrammarError(iline, sline, "Invalid closed markup");

        // fetch name and parent
        _name = attrs.find("name")->second;
        _parent = attrs.find("parent")->second;

        // make sure that the mesh part does not already exist in the mesh node
        if(_root_node.find_mesh_part(_name) != nullptr)
          throw Xml::ContentError(iline, sline, String("Mesh Part '") + _name + "' already exists in mesh node");

        // fetch chart
        {
          auto it = attrs.find("chart");
          if(it != attrs.end())
          {
            // try to find chart
            _chart_name = it->second;
            _linker.meshpart_link_to_chart(_name, _chart_name);
          }
        }

        // make sure the parent is root
        if(_parent != "root")
          throw Xml::ContentError(iline, sline, "Invalid mesh part parent; only 'root' is accepted");

        // get topology type
        String stopo_type = attrs.find("topology")->second;
        if(stopo_type == "none")
          _topo_type = TopoType::none;
        else if(stopo_type == "full")
          _topo_type = TopoType::full;
        else if(stopo_type == "parent")
        {
          _topo_type = TopoType::parent;
          _linker.meshpart_deduct_topology(_name);
        }
        else
          throw Xml::ContentError(iline, sline, "Invalid topology attribute '" + stopo_type + "'");

        // split sizes
        std::deque<String> ssize = attrs.find("size")->second.split_by_whitespaces();
        if(ssize.size() > std::size_t(Shape_::dimension+1))
          throw Xml::ContentError(iline, sline, "Invalid mesh part size count");

        // parse sizes
        _sizes.resize(std::size_t(Shape_::dimension+1), Index(0));
        for(std::size_t i(0); i < ssize.size(); ++i)
        {
          if(!ssize.at(i).parse(_sizes.at(i)))
            throw Xml::ContentError(iline, sline, "Failed to parse mesh part size");
        }

        // try to create the mesh part
        _mesh_part = new MeshPartType(_sizes.data(), _topo_type != TopoType::none);

        // initialize
        _have_topology.resize(std::size_t(Shape_::dimension), 0);
        _have_mapping.resize(std::size_t(Shape_::dimension+1), 0);
      }

      virtual void close(int iline, const String& sline) override
      {
        // make sure we have all mappings
        for(std::size_t i(0); i < _have_mapping.size(); ++i)
        {
          if((_have_mapping.at(i) == 0) && (_sizes.at(i) > Index(0)))
            throw Xml::GrammarError(iline, sline, "Missing mapping for dimension " + stringify(i));
        }

        if(_topo_type == TopoType::full)
        {
          // make sure we have the whole topology
          bool topo_all(true);
          for(std::size_t i(0); i < _have_topology.size(); ++i)
          {
            // Note: there may be no entities of this dimension here
            if(_sizes.at(i+1) > 0)
              topo_all = topo_all && (_have_topology.at(i) != 0);
          }

          // Make sure we have either all or none topology entires
          if(!topo_all)
            throw Xml::GrammarError(iline, sline, "Incomplete topology");

          // compute remaining topology now
          if(topo_all)
            RedundantIndexSetBuilder<Shape_>::compute(*_mesh_part->get_topology());
        }
        else if(_topo_type == TopoType::parent)
        {
          // Note: this has been outsourced to the MeshNodeLinker class!
          // Create the MeshPart's topology from the parent mesh's topology if told so
          //_mesh_part->deduct_topology(*(_root_node.get_mesh()->get_topology()));
          //_linker.meshpart_deduct_topology(_name);
        }

        // Finally, insert mesh part into root node
        _root_node.add_mesh_part(_name, _mesh_part/*, _chart_name, _chart*/);

        // Reset pointer
        _mesh_part = nullptr;
      }

      virtual std::shared_ptr<Xml::MarkupParser> markup(int iline, const String& sline, const String& name) override
      {
        // mapping?
        if(name == "Mapping")
        {
          // create mapping parser
          return std::make_shared<MappingParser<Shape_>>(_mesh_part->get_target_set_holder(), _have_mapping);
        }

        // topology?
        if(name == "Topology")
        {
          if(_topo_type == TopoType::none)
            throw Xml::ContentError(iline, sline, "Unexpected Topology found");

          // create a topology parser
          return std::make_shared<TopologyParser<Shape_>>(*_mesh_part->get_topology(), _have_topology);
        }

        // attribute?
        if(name == "Attribute")
        {
          return std::make_shared<AttributeParser<MeshPartType>>(*_mesh_part);
        }

        // invalid
        return nullptr;
      }

      virtual bool content(int, const String&) override
      {
        return false; // no content allowed
      }
    };

    class PatchParser :
      public Xml::MarkupParser
    {
    protected:
      Adjacency::DynamicGraph& _patches;
      Index _rank;
      Index _size;
      Index _read;

    public:
      explicit PatchParser(Adjacency::DynamicGraph& patches) :
        _patches(patches),
        _rank(~Index(0)),
        _size(~Index(0)),
        _read(Index(0))
      {
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("rank", true);
        attrs.emplace("size", true);
        return true;
      }

      virtual void create(int iline, const String& sline, const String&, const std::map<String, String>& attrs, bool) override
      {
        // try to parse rank and size
        if(!attrs.find("rank")->second.parse(_rank))
          throw Xml::ContentError(iline, sline, "Cannot parse patch rank index");
        if(!attrs.find("size")->second.parse(_size))
          throw Xml::ContentError(iline, sline, "Cannot parse patch size");

        // ensure that the rank is valid
        if(_rank >= _patches.get_num_nodes_domain())
          throw Xml::ContentError(iline, sline, "Patch rank index is out of bounds");
      }

      virtual void close(int iline, const String& sline) override
      {
        if(_read < _size)
          throw Xml::GrammarError(iline, sline, "Invalid terminator; expected index");
      }

      virtual bool content(int iline, const String& sline) override
      {
        // make sure that we do not read more points than expected
        if(_read >= _size)
          throw Xml::ContentError(iline, sline, "Invalid content; expected terminator");

        // try to parse element index
        Index elem(0);
        if(!sline.parse(elem))
          throw Xml::ContentError(iline, sline, "Failed to parse element index");

        if(elem >= _patches.get_num_nodes_image())
          throw Xml::ContentError(iline, sline, "Patch element index is out of bounds");

        // insert adjacency
        _patches.insert(_rank, elem);

        // okay, another one processed
        ++_read;

        return true;
      }

      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
      {
        return nullptr;
      }
    }; // class PatchParser

    class PartitionParser :
      public Xml::MarkupParser
    {
    protected:
      PartitionSet& _part_set;
      String _name;
      int _prio;
      int _level;
      int _num_ranks, _num_elems;
      Adjacency::DynamicGraph _patches;

    public:
      explicit PartitionParser(PartitionSet& part_set) :
        _part_set(part_set),
        _name(),
        _prio(0),
        _level(0),
        _num_ranks(0),
        _num_elems(0),
        _patches()
      {
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("size", true);
        attrs.emplace("name", false);
        attrs.emplace("priority", false);
        attrs.emplace("level", false);
        return true;
      }

      virtual void create(int iline, const String& sline, const String&, const std::map<String, String>& attrs, bool) override
      {
        // split sizes
        std::deque<String> ssize = attrs.find("size")->second.split_by_whitespaces();
        if(ssize.size() != std::size_t(2))
          throw Xml::ContentError(iline, sline, "Invalid partition size");

        // parse sizes
        if(!ssize.front().parse(_num_ranks))
          throw Xml::ContentError(iline, sline, "Failed to parse partition rank count");
        if(!ssize.back().parse(_num_elems))
          throw Xml::ContentError(iline, sline, "Failed to parse partition element count");

        // fetch name
        {
          auto it = attrs.find("name");
          if(it != attrs.end())
            _name = it->second;
        }

        // fetch priority
        {
          auto it = attrs.find("priority");
          if(it != attrs.end())
          {
            if(!it->second.parse(_prio))
              throw Xml::ContentError(iline, sline, "Cannot parse partition priority");
          }
        }

        // fetch level
        {
          auto it = attrs.find("level");
          if(it != attrs.end())
          {
            if(!it->second.parse(_level))
              throw Xml::ContentError(iline, sline, "Cannot parse partition level");
            if(_level < 0)
              throw Xml::ContentError(iline, sline, "Invalid negative partition level");
          }
        }

        // okay, let's create a new partition
        _patches = Adjacency::DynamicGraph(Index(_num_ranks), Index(_num_elems));
      }

      virtual void close(int, const String&) override
      {
        // add our partition to the set
        _part_set.add_partition(Partition(_patches, _name, _prio, _level));
      }

      virtual bool content(int, const String&) override
      {
        return false;
      }

      virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String& name) override
      {
        if(name == "Patch")
          return std::make_shared<PatchParser>(_patches);

        return nullptr;
      }
    }; // class PartitionParser

    template<typename RootMesh_>
    class MeshNodeParser :
      public Xml::MarkupParser
    {
    public:
      typedef RootMesh_ RootMeshType;
      typedef RootMeshNode<RootMesh_> RootMeshNodeType;
      typedef MeshAtlas<RootMesh_> MeshAtlasType;

    protected:
      /// root mesh node
      RootMeshNodeType& _root_node;
      /// mesh atlas
      MeshAtlasType& _mesh_atlas;
      /// partition set
      PartitionSet* _part_set;
      /// mesh node linker
      MeshNodeLinker<RootMesh_>& _linker;
      /// mesh declarator
      String _mesh_decl;

    public:
      explicit MeshNodeParser(RootMeshNodeType& root_node, MeshAtlasType& mesh_atlas, PartitionSet* part_set,  MeshNodeLinker<RootMesh_>& linker) :
        _root_node(root_node),
        _mesh_atlas(mesh_atlas),
        _part_set(part_set),
        _linker(linker),
        _mesh_decl()
      {
      }

      virtual bool attribs(std::map<String,bool>& attrs) const override
      {
        attrs.emplace("version", true); // mandatory
        attrs.emplace("mesh", false); // optional
        return true;
      }

      virtual void create(int iline, const String& sline, const String& name, const std::map<String, String>& attrs, bool) override
      {
        // verify marker name
        if(name != "FeatMeshFile")
          throw Xml::GrammarError(iline, sline, "Invalid root markup; expected '<FeatMeshFile ...>'");

        // verify file version
        int version(0);
        if(!attrs.find("version")->second.parse(version))
          throw Xml::ContentError(iline, sline, "Cannot parse file version");
        if(version != 1)
          throw Xml::ContentError(iline, sline, "Invalid file version; expected 1");

        // have a mesh declaration?
        auto it = attrs.find("mesh");
        if(it != attrs.end())
        {
          // okay, store the declarator
          _mesh_decl = it->second;
        }
      }

      virtual void close(int, const String&) override
      {
      }

      virtual bool content(int, const String&) override
      {
        return false; // no content allowed
      }

      virtual std::shared_ptr<MarkupParser> markup(int iline, const String& sline, const String& name) override
      {
        // check for valid markups
        if(name == "Info")
        {
          // valid info markup
          return std::make_shared<Xml::DummyParser>();
        }
        if(name == "Chart")
        {
          // create a chart parser
          return std::make_shared<ChartParser<RootMeshType>>(_mesh_atlas);
        }
        if(name == "Mesh")
        {
          if(_root_node.get_mesh() != nullptr)
          {
            // that's the second one!
            throw Xml::GrammarError(iline, sline, "Invalid second Mesh in file");
          }

          // valid mesh markup
          return std::make_shared<MeshParser<RootMeshType>>(_root_node);
        }
        if(name == "MeshPart")
        {
          // valid meshpart markup
          return std::make_shared<MeshPartParser<RootMeshType>>(_root_node, _linker);
        }
        if(name == "Partition")
        {
          // valid partition markup
          if(_part_set != nullptr)
            return std::make_shared<PartitionParser>(*_part_set);
          else
          {
            // the user is not interested in partitions, so create a dummy parser
            return std::make_shared<Xml::DummyParser>();
          }
        }

        // anything else is invalid
        return nullptr;
      }
    }; // class MeshNodeParser

    /// \endcond

    /**
     * \brief Mesh file reader class
     *
     * This class implements the reader which parses the XML-based FEAT mesh files
     * into the corresponding Geometry class objects of type MeshAtlas, RootMeshNode and
     * PartitionSet.
     *
     * The basic usage is as follows:
     * -# Create an object of the class MeshFileReader
     * -# Add all input file streams to the object by calling the add_stream() function.
     * -# Read the root markup by calling the read_root_markup() function.
     * -# Create MeshAtlas and RootMeshNode objects of the corresponding mesh- and shape-type.
     * -# Parse the file contents into the corresponding MeshAtlas, RootMeshNode and PartitionSet objects
     *    by calling the parse() function.
     *
     * For more details on meshes, see the related doxygen page \ref mesh_file_format.
     *
     * \author Peter Zajac
     */
    class MeshFileReader
    {
    public:
      /// mesh type enumeration
      enum class MeshType
      {
        /// unknown mesh type
        unknown = 0,
        /// conformal mesh type
        conformal
      };

      /// shape type enumeration
      enum class ShapeType
      {
        /// unknown shape type
        unknown = 0,
        /// simplex shape type
        simplex,
        /// hypercube shape type
        hypercube
      };

    protected:
      /// Our internally managed streams
      std::deque<std::shared_ptr<std::stringstream>> _streams;
      /// Our Xml scanner objects
      std::deque<std::shared_ptr<Xml::Scanner>> _scanners;
      /// Did we read the root markup yet?
      bool _have_root_markup;
      /// The parsed mesh type
      String _mesh_type_string;

      // Split up mesh type
      MeshType _mesh_type;
      // Split up shape type
      ShapeType _shape_type;
      // Mesh shape dimension
      int _shape_dim;
      // Mesh world dimension
      int _world_dim;

    public:
      /// default constuctor
      explicit MeshFileReader() :
        _have_root_markup(false),
        _mesh_type_string(),
        _mesh_type(MeshType::unknown),
        _shape_type(ShapeType::unknown),
        _shape_dim(0),
        _world_dim(0)
      {
      }

      /**
       * \brief Input-Stream constructor
       *
       * \param[in] is
       * The input stream that is to be parsed.
       */
      explicit MeshFileReader(std::istream& is) :
        _have_root_markup(false),
        _mesh_type_string(),
        _mesh_type(MeshType::unknown),
        _shape_type(ShapeType::unknown),
        _shape_dim(0),
        _world_dim(0)
      {
        // create a new scanner object
        _scanners.push_back(std::make_shared<Xml::Scanner>(is));
      }

      MeshFileReader(const MeshFileReader&) = delete;
      MeshFileReader& operator=(const MeshFileReader&) = delete;

      /// virtual destructor
      virtual ~MeshFileReader()
      {
        _scanners.clear();
        _streams.clear();
      }

      /**
       * \brief Adds an input stream to the list of streams to be parsed.
       *
       * \param[in] is
       * The input stream that is to be added to the list.
       */
      void add_stream(std::istream& is)
      {
        XASSERTM(!_have_root_markup, "cannot add new stream after reading root");

        // create a new scanner object
        _scanners.push_back(std::make_shared<Xml::Scanner>(is));
      }

      /**
       * \brief Adds a list of mesh files to the list of streams to be parsed.
       *
       * \param[in] comm
       * The communicator to be used for parallel file I/O.
       *
       * \param[in] filenames
       * A list of mesh filenames that are to be parsed.
       *
       * \param[in] dirpath
       * The path in which the mesh files are located.
       */
      void add_mesh_files(const Dist::Comm& comm, const std::deque<String>& filenames, String dirpath = "")
      {
        XASSERTM(!_have_root_markup, "cannot add new stream after reading root");

        if(filenames.empty())
          return;

        // ensure that the directory path ends with a path separator
        if(!dirpath.empty())
        {
          /// \platformswitch windows uses backslash as path separator
#ifdef _WIN32
          if((dirpath.back() != '/') && (dirpath.back() != '\\'))
            dirpath.push_back('\\');
#else
          if(dirpath.back() != '/')
            dirpath.push_back('/');
#endif
        }

        // read all files
        for(std::size_t i(0); i < filenames.size(); ++i)
        {
          // build file path
          String filepath = dirpath + filenames.at(i);

          // create a new stream
          auto stream = std::make_shared<std::stringstream>();

          // read the stream
          DistFileIO::read_common(*stream, filepath, comm);

          // add stream to deque
          _streams.push_back(stream);

          // add to mesh reader
          add_stream(*_streams.back());
        }
      }

      /**
       * \brief Adds a list of mesh files to the list of streams to be parsed.
       *
       * \param[in] filenames
       * A list of mesh filenames that are to be parsed.
       *
       * \param[in] dirpath
       * The path in which the mesh files are located.
       */
      void add_mesh_files(const std::deque<String>& filenames, String dirpath = "")
      {
        Dist::Comm comm = Dist::Comm::world();
        add_mesh_files(comm, filenames, dirpath);
      }

      /**
       * \brief Returns the mesh-type string from the root markup node.
       *
       * \note
       * This function can only be called after calling the #read_root_markup()
       * function as it will always return an empty string otherwise.
       *
       * \returns
       * The mesh-type string from the root markup node.
       */
      const String& get_meshtype_string() const
      {
        return _mesh_type_string;
      }

      /**
       * \brief Returns the mesh-type from the root markup node.
       *
       * \note
       * This function can only be called after calling the #read_root_markup()
       * function as it will always return MeshType::unknown otherwise.
       *
       * \returns
       * The mesh-type from the root markup node.
       */
      MeshType get_mesh_type() const
      {
        return _mesh_type;
      }

      /**
       * \brief Returns the shape-type from the root markup node.
       *
       * \note
       * This function can only be called after calling the #read_root_markup()
       * function as it will always return ShapeType::unknown otherwise.
       *
       * \returns
       * The shape-type from the root markup node.
       */
      ShapeType get_shape_type() const
      {
        return _shape_type;
      }

      /**
       * \brief Returns the shape-dimension from the root markup node.
       *
       * \note
       * This function can only be called after calling the #read_root_markup()
       * function as it will always return 0 otherwise.
       *
       * \returns
       * The shape-dimension from the root markup node.
       */
      int get_shape_dim() const
      {
        return _shape_dim;
      }

      /**
       * \brief Returns the world-dimension from the root markup node.
       *
       * \note
       * This function can only be called after calling the #read_root_markup()
       * function as it will always return 0 otherwise.
       *
       * \returns
       * The world-dimension from the root markup node.
       */
      int get_world_dim() const
      {
        return _world_dim;
      }

      /**
       * \brief Reads the root markup from the input stream and analyses it.
       */
      void read_root_markup()
      {
        XASSERTM(!_have_root_markup, "root has already been read");

        // loop over all scanners
        for(auto it = _scanners.begin(); it != _scanners.end(); ++it)
        {
          // get the scanner
          Xml::Scanner& scanner = **it;

          // try to read root node
          scanner.read_root();

          // verify root marker
          if(scanner.get_cur_name() != "FeatMeshFile")
            scanner.throw_grammar("Invalid mesh file");

          // get attributes
          const auto& attr = scanner.get_cur_attribs();

          // check version attribute
          auto it_v = attr.find("version");
          if(it_v == attr.end())
            scanner.throw_grammar("Mesh file is missing mandatory attribute 'version'");
          int iversion(0);
          if(!(it_v->second.parse(iversion)))
            scanner.throw_grammar("Failed to parse mesh file version attribute");
          if(iversion != 1)
            scanner.throw_grammar("Invalid mesh file version");

          // check mesh attribute
          auto it_m = attr.find("mesh");
          if(it_m != attr.end())
          {
            // did we read a mesh type already?
            if(!_mesh_type_string.empty())
            {
              // make sure that the mesh types are identical
              XASSERTM(it_m->second == _mesh_type_string, "conflicting mesh types");
            }
            else
            {
              // save the mesh type string
              _mesh_type_string = it_m->second;

              // split up mesh type string
              std::deque<String> mts = _mesh_type_string.split_by_string(":");

              if(mts.size() != std::size_t(4))
                scanner.throw_grammar("Mesh file root markup has invalid 'mesh' attribute");

              // check mesh type
              if(mts.at(0) == "conformal")
                _mesh_type = MeshType::conformal;
              else
                scanner.throw_grammar("Invalid 'mesh' attribute mesh type");

              // check shape type
              if(mts.at(1) == "simplex")
                _shape_type = ShapeType::simplex;
              else if(mts.at(1) == "hypercube")
                _shape_type = ShapeType::hypercube;
              else
                scanner.throw_grammar("Invalid 'mesh' attribute shape type");

              // parse shape dimension
              if(!mts.at(2).parse(_shape_dim) || (_shape_dim <= 0))
                scanner.throw_content("Invalid 'mesh' attribute shape dimension");

              // parse world dimension
              if(!mts.at(3).parse(_world_dim) || (_world_dim <= 0))
                scanner.throw_content("Invalid 'mesh' attribute world dimension");
            }
          }

          // continue with next stream
        }

        // alright, root markup read
        _have_root_markup = true;
      }

      /**
       * \brief Parses the mesh file into a mesh node and a mesh atlas.
       *
       * \param[in,out] linker
       * A linker that the post-parse tasks are to be added to.
       *
       * \param[in,out] root_mesh_node
       * The root mesh node into which the mesh and the mesh parts are to be added.
       *
       * \param[in,out] mesh_atlas
       * The mesh atlas into which charts are to be added. Is also used to search
       * for charts for mesh parts.
       *
       * \param[in,out] part_set
       * A pointer to the partition set that partitions are added to.
       * May be \p nullptr, if the partitions are to be ignored.
       */
      template<typename RootMesh_>
      void parse(
        MeshNodeLinker<RootMesh_>& linker,
        RootMeshNode<RootMesh_>& root_mesh_node,
        MeshAtlas<RootMesh_>& mesh_atlas,
        PartitionSet* part_set = nullptr)
      {
        // read root markup unless it has already been read
        if(!_have_root_markup)
          read_root_markup();

        // loop over all scanners
        for(auto it = _scanners.begin(); it != _scanners.end(); ++it)
        {
          // create a corresponding parser and scan
          (*it)->set_root_parser(std::make_shared<MeshNodeParser<RootMesh_>>(root_mesh_node, mesh_atlas, part_set, linker));
          (*it)->scan();
        }
      }

      /**
       * \brief Parses the mesh file into a mesh node and a mesh atlas.
       *
       * \param[in,out] root_mesh_node
       * The root mesh node into which the mesh and the mesh parts are to be added.
       *
       * \param[in,out] mesh_atlas
       * The mesh atlas into which charts are to be added. Is also used to search
       * for charts for mesh parts.
       *
       * \param[in,out] part_set
       * A pointer to the partition set that partitions are added to.
       * May be \p nullptr, if the partitions are to be ignored.
       */
      template<typename RootMesh_>
      void parse(
        RootMeshNode<RootMesh_>& root_mesh_node,
        MeshAtlas<RootMesh_>& mesh_atlas,
        PartitionSet* part_set = nullptr)
      {
        // create linker
        MeshNodeLinker<RootMesh_> linker(root_mesh_node, mesh_atlas);

        // parse
        parse(linker, root_mesh_node, mesh_atlas, part_set);

        // execute linker
        linker.execute();
      }

      /**
       * \brief Parses the mesh file into a mesh node and a mesh atlas.
       *
       * \param[in,out] mesh_atlas
       * The mesh atlas into which charts are to be added. Is also used to search
       * for charts for mesh parts.
       *
       * \param[in,out] part_set
       * A pointer to the partition set that partitions are added to.
       * May be \p nullptr, if the partitions are to be ignored.
       *
       * \returns
       * A shared pointer to the root mesh node containing the mesh and the mesh parts.
       * This mesh node is automatically linked to the given mesh atlas.
       */
      template<typename RootMesh_>
      std::shared_ptr<RootMeshNode<RootMesh_>> parse(
        MeshAtlas<RootMesh_>& mesh_atlas,
        PartitionSet* part_set = nullptr)
      {
        // create new root mesh node
        std::shared_ptr<RootMeshNode<RootMesh_>> root_mesh_node =
          std::make_shared<RootMeshNode<RootMesh_>>(nullptr, &mesh_atlas);

        // parse
        this->parse(*root_mesh_node, mesh_atlas, part_set);

        // return root mesh node
        return root_mesh_node;
      }
    }; // class MeshFileReader

    /**
     * \brief World dimension dependent helper class for parsing charts
     *
     * \tparam RootMesh_
     * Type of the mesh the chart is supposed to work with.
     *
     * \tparam world_dim
     * World dimension of the RootMesh_. This is separate because it is used for explicit specialization.
     *
     * This class filters which chart(-parser) classes will be instantiated for each world_dim value, because (all?)
     * charts are tied to a certain dimension. This is the generic implementation for values of world_dim which do
     * not have any charts, and serves as an interface documentation.
     *
     */
    template<typename RootMesh_, int world_dim>
    struct DimensionalChartHelper
    {
      static_assert(RootMesh_::world_dim == world_dim, "Nonmatching world_dim.");
      /**
       * \brief Creates a parser for a chart, its type identified by a name String
       *
       * \param[in] name
       * Type name of the chart to parse
       *
       * \param[out] chart
       * Where to parse the chart to
       *
       */
      static std::shared_ptr<Xml::MarkupParser> markup(const String& DOXY(name),
      Atlas::ChartBase<RootMesh_>*& DOXY(chart))
      {
        return nullptr;
      }
    };

    /// \cond internal
    template<typename RootMesh_>
    struct DimensionalChartHelper<RootMesh_,2>
    {
      static std::shared_ptr<Xml::MarkupParser> markup(const String& name,
      Atlas::ChartBase<RootMesh_>*& chart)
      {
        if(name == "Circle")
          return std::make_shared<Atlas::CircleChartParser<RootMesh_>>(chart);
        if(name == "Bezier")
          return std::make_shared<Atlas::BezierChartParser<RootMesh_>>(chart);

        return nullptr;
      }
    };

    template<typename RootMesh_>
    struct DimensionalChartHelper<RootMesh_,3>
    {
      static std::shared_ptr<Xml::MarkupParser> markup(const String& name,
      Atlas::ChartBase<RootMesh_>*& chart)
      {
        if(name == "Sphere")
          return std::make_shared<Atlas::SphereChartParser<RootMesh_>>(chart);
        if(name == "SurfaceMesh")
          return std::make_shared<Atlas::SurfaceMeshChartParser<RootMesh_>>(chart);
        if(name == "Extrude")
          return std::make_shared<Atlas::ExtrudeChartParser<RootMesh_>>(chart);

        return nullptr;
      }
    };
    /// \endcond

#ifdef FEAT_EICKT
    extern template class MeshNodeLinker<ConformalMesh<Shape::Simplex<2>, 2, Real>>;
    extern template class MeshNodeLinker<ConformalMesh<Shape::Simplex<3>, 3, Real>>;
    extern template class MeshNodeLinker<ConformalMesh<Shape::Hypercube<2>, 2, Real>>;
    extern template class MeshNodeLinker<ConformalMesh<Shape::Hypercube<3>, 3, Real>>;

    extern template void MeshFileReader::parse<ConformalMesh<Shape::Simplex<2>, 2, Real>>(
      MeshNodeLinker<ConformalMesh<Shape::Simplex<2>, 2, Real>>&,
      RootMeshNode<ConformalMesh<Shape::Simplex<2>, 2, Real>>&,
      MeshAtlas<ConformalMesh<Shape::Simplex<2>, 2, Real>>&,
      PartitionSet*);
    extern template void MeshFileReader::parse<ConformalMesh<Shape::Simplex<3>, 3, Real>>(
      MeshNodeLinker<ConformalMesh<Shape::Simplex<3>, 3, Real>>&,
      RootMeshNode<ConformalMesh<Shape::Simplex<3>, 3, Real>>&,
      MeshAtlas<ConformalMesh<Shape::Simplex<3>, 3, Real>>&,
      PartitionSet*);
    extern template void MeshFileReader::parse<ConformalMesh<Shape::Hypercube<2>, 2, Real>>(
      MeshNodeLinker<ConformalMesh<Shape::Hypercube<2>, 2, Real>>&,
      RootMeshNode<ConformalMesh<Shape::Hypercube<2>, 2, Real>>&,
      MeshAtlas<ConformalMesh<Shape::Hypercube<2>, 2, Real>>&,
      PartitionSet*);
    extern template void MeshFileReader::parse<ConformalMesh<Shape::Hypercube<3>, 3, Real>>(
      MeshNodeLinker<ConformalMesh<Shape::Hypercube<3>, 3, Real>>&,
      RootMeshNode<ConformalMesh<Shape::Hypercube<3>, 3, Real>>&,
      MeshAtlas<ConformalMesh<Shape::Hypercube<3>, 3, Real>>&,
      PartitionSet*);
#endif // FEAT_EICKT
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_MESH_FILE_READER_HPP
