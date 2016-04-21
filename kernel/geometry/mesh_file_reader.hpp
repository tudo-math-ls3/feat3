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
#include <kernel/geometry/atlas/polyline.hpp>
#include <kernel/geometry/atlas/surface_mesh.hpp>
#include <kernel/geometry/atlas/sphere.hpp>
#include <kernel/geometry/atlas/tube.hpp>
#include <kernel/util/xml_scanner.hpp>

namespace FEAST
{
  namespace Geometry
  {
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
        if(_atlas.find_mesh_chart(_chart_name) != nullptr)
          throw InternalError(String("Chart '") + _chart_name + "' already exists in atlas");

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
        // What have we here?
        if(name == "Bezier")       return std::make_shared<Atlas::BezierChartParser<RootMesh_>>(_chart);
        if(name == "Circle")       return std::make_shared<Atlas::CircleChartParser<RootMesh_>>(_chart);
        if(name == "Polyline")     return std::make_shared<Atlas::PolylineChartParser<RootMesh_>>(_chart);
        // \todo: Implement SphereChartParser
        //if(name == "Sphere")       return std::make_shared<Atlas::SphereChartParser<RootMesh_>>(_chart);
        if(name == "SurfaceMesh")  return std::make_shared<Atlas::SurfaceMeshChartParser<RootMesh_>>(_chart);
        // \todo: Implement TubeChartParser
        //if(name == "Tube")       return std::make_shared<Atlas::TubeChartParser<RootMesh_>>(_chart);
        if(name == "Extrude")      return std::make_shared<Atlas::ExtrudeChartParser<RootMesh_>>(_chart);

        return nullptr;
      }
    }; // class ChartParser<...>

    template<int num_coords_, int stride_, typename Coord_>
    class VerticesParser :
      public Xml::MarkupParser
    {
    public:
      typedef VertexSet<num_coords_, stride_, Coord_> VertexSetType;

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
        std::deque<String> scoords;
        sline.split_by_charset(scoords);

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

    /// \cond internal
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
          throw InternalError(__func__, __FILE__, __LINE__, "Thou shall not arrive here");
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
    /// \endcond

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

        // initialise the corresponding topology:
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
        std::deque<String> sidx;
        sline.split_by_charset(sidx);

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

        // initialise the corresponding mapping:
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
      typedef MeshAttribute<DataType_> AttribType;

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
        // Because adding the MeshAttribute to a MeshPart passes ownership, we just need to nullify the pointer here
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
        _attrib = new AttribType(_count, int(_my_dim), 0);
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
        std::deque<String> svals;
        sline.split_by_charset(svals);

        // check size
        if(svals.size() != std::size_t(_my_dim))
          throw Xml::ContentError(iline, sline, "Invalid number of values");

        // try to parse all coords
        for(int i(0); i < int(_my_dim); ++i)
        {
          if(!svals.at(std::size_t(i)).parse((*_attrib)[_read][i]))
            throw Xml::ContentError(iline, sline, "Failed to parse value");
        }

        // okay, another one processed
        ++_read;

        return true;
      }
    };

    template<typename Mesh_>
    class MeshParser;

    template<typename Shape_, int num_coords_, int stride_, typename Coord_>
    class MeshParser<ConformalMesh<Shape_, num_coords_, stride_, Coord_>> :
      public Xml::MarkupParser
    {
    public:
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;

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
        std::deque<String> stype;
        attrs.find("type")->second.split_by_charset(stype, ":");

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
        std::deque<String> ssize;
        attrs.find("size")->second.split_by_charset(ssize);
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

        // initialise
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
          return std::make_shared<VerticesParser<num_coords_, stride_, Coord_>>(_mesh->get_vertex_set());
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

    template<typename Shape_, int num_coords_, int stride_, typename Coord_>
    class MeshPartParser<ConformalMesh<Shape_, num_coords_, stride_, Coord_>> :
      public Xml::MarkupParser
    {
    public:
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
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
      MeshAtlas<MeshType>& _atlas;
      MeshPartType* _mesh_part;
      ChartType* _chart;
      TopoType _topo_type;
      String _name, _parent, _chart_name;
      std::deque<int> _have_topology;
      std::deque<int> _have_mapping;
      std::vector<Index> _sizes;

    public:
      explicit MeshPartParser(RootMeshNode<MeshType>& root_node, MeshAtlas<MeshType>& atlas) :
        _root_node(root_node),
        _atlas(atlas),
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
            _chart = _atlas.find_mesh_chart(_chart_name);
            if(_chart == nullptr)
              throw Xml::ContentError(iline, sline, "Unknown chart name: " + _chart_name);
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
          _topo_type = TopoType::parent;
        else
          throw Xml::ContentError(iline, sline, "Invalid topology attribute '" + stopo_type + "'");

        // split sizes
        std::deque<String> ssize;
        attrs.find("size")->second.split_by_charset(ssize);
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

        // initialise
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
          // Create the MeshPart's topology from the parent mesh's topology if told so
          _mesh_part->deduct_topology(*(_root_node.get_mesh()->get_topology()));
        }

        // Finally, insert mesh part into root node
        _root_node.add_mesh_part(_name, _mesh_part, _chart_name, _chart);

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
      /// mesh declarator
      String _mesh_decl;

    public:
      explicit MeshNodeParser(RootMeshNodeType& root_node, MeshAtlasType& mesh_atlas) :
        _root_node(root_node),
        _mesh_atlas(mesh_atlas),
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
          return std::make_shared<MeshPartParser<RootMeshType>>(_root_node, _mesh_atlas);
        }

        // anything else is invalid
        return nullptr;
      }
    }; // class MeshFileParser

    /**
     * \brief Mesh file reader class
     *
     * This class implements the reader which parses the XML-based Feat mesh files.
     *
     * \author Peter Zajac
     */
    class MeshFileReader
    {
    public:
      enum class MeshType
      {
        unknown = 0,
        conformal
      };

      enum class ShapeType
      {
        unknown = 0,
        simplex,
        hypercube
      };

    protected:
      /// The input stream to read from
      std::istream& _istream;
      /// Our Xml scanner object
      Xml::Scanner _scanner;
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
      explicit MeshFileReader(std::istream& is) :
        _istream(is),
        _scanner(_istream),
        _have_root_markup(false),
        _mesh_type_string(),
        _mesh_type(MeshType::unknown),
        _shape_type(ShapeType::unknown),
        _shape_dim(0),
        _world_dim(0)
      {

      }

      MeshFileReader(const MeshFileReader&) = delete;
      MeshFileReader& operator=(const MeshFileReader&) = delete;

      virtual ~MeshFileReader()
      {
      }

      const String& get_meshtype_string() const
      {
        return _mesh_type_string;
      }

      MeshType get_mesh_type() const
      {
        return _mesh_type;
      }

      ShapeType get_shape_type() const
      {
        return _shape_type;
      }

      int get_shape_dim() const
      {
        return _shape_dim;
      }

      int get_world_dim() const
      {
        return _world_dim;
      }

      /**
       * \brief Reads the root markup from the input stream and analyses it.
       */
      void read_root_markup()
      {
        // try to read root node
        _scanner.read_root();

        // verify root marker
        if(_scanner.get_cur_name() != "FeatMeshFile")
          _scanner.throw_grammar("Invalid mesh file");

        // get attributes
        const auto& attr = _scanner.get_cur_attribs();

        // check version attribute
        auto it_v = attr.find("version");
        if(it_v == attr.end())
          _scanner.throw_grammar("Mesh file is missing mandatory attribute 'version'");
        int iversion(0);
        if(!(it_v->second.parse(iversion)))
          _scanner.throw_grammar("Failed to parse mesh file version attribute");
        if(iversion != 1)
          _scanner.throw_grammar("Invalid mesh file version");

        // check mesh attribute
        auto it_m = attr.find("mesh");
        if(it_m != attr.end())
        {
          // save the mesh type string
          _mesh_type_string = it_m->second;

          // split up mesh type string
          std::deque<String> mts;
          _mesh_type_string.split_by_charset(mts, ":");

          if(mts.size() != std::size_t(4))
            _scanner.throw_grammar("Mesh file root markup has invalid 'mesh' attribute");

          // check mesh type
          if(mts.at(0) == "conformal")
            _mesh_type = MeshType::conformal;
          else
            _scanner.throw_grammar("Invalid 'mesh' attribute mesh type");

          // check shape type
          if(mts.at(1) == "simplex")
            _shape_type = ShapeType::simplex;
          else if(mts.at(1) == "hypercube")
            _shape_type = ShapeType::hypercube;
          else
            _scanner.throw_grammar("Invalid 'mesh' attribute shape type");

          // parse shape dimension
          if(!mts.at(2).parse(_shape_dim) || (_shape_dim <= 0))
            _scanner.throw_content("Invalid 'mesh' attribute shape dimension");

          // parse world dimension
          if(!mts.at(3).parse(_world_dim) || (_world_dim <= 0))
            _scanner.throw_content("Invalid 'mesh' attribute world dimension");
        }

        // alright, root markup read
        _have_root_markup = true;
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
       */
      template<typename RootMesh_>
      void parse(RootMeshNode<RootMesh_>& root_mesh_node, MeshAtlas<RootMesh_>& mesh_atlas)
      {
        // read root markup unless it has already been read
        if(!_have_root_markup)
          read_root_markup();

        // create a corresponding parser
        _scanner.set_root_parser(std::make_shared<MeshNodeParser<RootMesh_>>(root_mesh_node, mesh_atlas));
        _scanner.scan();
      }
    }; // class MeshFileReader
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_FILE_READER_HPP
