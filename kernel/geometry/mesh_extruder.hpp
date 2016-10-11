#pragma once
#ifndef KERNEL_GEOMETRY_MESH_EXTRUDER_HPP
#define KERNEL_GEOMETRY_MESH_EXTRUDER_HPP 1

#include <kernel/geometry/atlas/extrude.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/index_calculator.hpp>
#include <kernel/geometry/partition_set.hpp>
#include <kernel/util/xml_scanner.hpp>

#include <sstream>

namespace FEAT
{
  namespace Geometry
  {
    // Forward declarations of generic extruder class templates
    // Note: The following templates are only implemented for
    //       SourceMesh_ = ConformalMesh<Hypercube<2>, 2, 2, ...>

    template<typename SourceMesh_>
    class MeshExtruder;

    template<typename SourceMesh_>
    class MeshExtruderFactory;

    template<typename SourceMesh_>
    class MeshPartExtruderFactory;

    template<typename SourceMesh_>
    class MeshPartSliceExtruderFactory;

    /**
     * \brief Quadrilateral-to-Hexahedral mesh extruder
     *
     * \author Peter Zajac
     */
    template<typename Coord_>
    class MeshExtruder<ConformalMesh<Shape::Hypercube<2>, 2, 2, Coord_>>
    {
    public:
      typedef Coord_ CoordType;

      typedef ConformalMesh<Shape::Hypercube<2>, 2, 2, Coord_> QuadMesh;
      typedef ConformalMesh<Shape::Hypercube<3>, 3, 3, Coord_> HexaMesh;

      typedef MeshAtlas<QuadMesh> QuadAtlas;
      typedef MeshAtlas<HexaMesh> HexaAtlas;

      typedef Atlas::ChartBase<QuadMesh> QuadChart;
      typedef Atlas::ChartBase<HexaMesh> HexaChart;

      typedef RootMeshNode<QuadMesh> QuadRootNode;
      typedef RootMeshNode<HexaMesh> HexaRootNode;

      typedef typename QuadMesh::VertexSetType QuadVertexSet;
      typedef typename HexaMesh::VertexSetType HexaVertexSet;

      typedef typename QuadMesh::IndexSetHolderType QuadTopology;
      typedef typename HexaMesh::IndexSetHolderType HexaTopology;

      typedef MeshPart<QuadMesh> QuadPart;
      typedef MeshPart<HexaMesh> HexaPart;

      typedef typename QuadPart::TargetSetHolderType QuadTrgSetHolder;
      typedef typename HexaPart::TargetSetHolderType HexaTrgSetHolder;

      typedef typename QuadPart::MeshAttributeType QuadAttrib;
      typedef typename HexaPart::MeshAttributeType HexaAttrib;

    protected:
      /// the number of slices in z-direction
      const Index _slices;
      /// minimal and maximal z-coordinate
      const CoordType _z_min, _z_max;
      /// mesh-part names for z-min/max boundary regions
      String _zmin_part_name, _zmax_part_name;
      /// transformation origin
      CoordType _origin_x, _origin_y;
      /// transformation angles
      CoordType _angle_y, _angle_p, _angle_r;
      /// transformation offset
      CoordType _offset_x, _offset_y, _offset_z;

    public:
      /**
       * \brief Creates the mesh extruder object
       *
       * \param[in] slices
       * The number of mesh slices in Z-direction. Must be > 0.
       *
       * \param[in] z_min, z_max
       * The minimal and maximal Z-coordinate for the extruded 3D mesh.
       *
       * \param[in] zmin_part_name, zmax_part_name
       * The names of the mesh-parts which are to be generated for the boundary
       * regions at \p z_min and \p z_max. If set to an empty string,
       * no mesh-part is generated for the corresponding boundary region.
       */
      explicit MeshExtruder(Index slices, CoordType z_min, CoordType z_max, String zmin_part_name, String zmax_part_name) :
        _slices(slices),
        _z_min(z_min),
        _z_max(z_max),
        _zmin_part_name(zmin_part_name),
        _zmax_part_name(zmax_part_name),
        _origin_x(CoordType(0)),
        _origin_y(CoordType(0)),
        _angle_y(CoordType(0)),
        _angle_p(CoordType(0)),
        _angle_r(CoordType(0)),
        _offset_x(CoordType(0)),
        _offset_y(CoordType(0)),
        _offset_z(CoordType(0))
      {
        XASSERT(slices > Index(0));
        XASSERT(z_min < z_max);
      }

      /**
       * \brief Sets the transformation origin
       *
       * \param[in] x,y
       * The coordinates of the 2D origin.
       */
      void set_origin(CoordType x, CoordType y)
      {
        _origin_x = x;
        _origin_y = y;
      }

      /**
       * \brief Sets the transformation angles
       *
       * \param[in] yaw, pitch, roll
       * The yaw-pitch-roll angles for the transformation.
       */
      void set_angles(CoordType yaw, CoordType pitch, CoordType roll)
      {
        _angle_y = yaw;
        _angle_p = pitch;
        _angle_r = roll;
      }

      /**
       * \brief Sets the transformation offset
       *
       * \param[in] x,y,z
       * The coordinates of the 3D offset.
       */
      void set_offset(CoordType x, CoordType y, CoordType z)
      {
        _offset_x = x;
        _offset_y = y;
        _offset_z = z;
      }

      /**
       * \brief Computes the extruded entity counts.
       *
       * \param[out] hexa_num_entities
       * The array of length 4 that receives the extruded entity counts.
       *
       * \param[in] quad_num_entities
       * The array of length 3 to be extruded.
       */
      void extrude_num_entities(Index hexa_num_entities[], const Index quad_num_entities[]) const
      {
        XASSERT(hexa_num_entities != nullptr);
        XASSERT(quad_num_entities != nullptr);
        hexa_num_entities[0] = quad_num_entities[0] * (_slices + 1);
        hexa_num_entities[1] = quad_num_entities[0] * _slices + quad_num_entities[1] * (_slices + 1);
        hexa_num_entities[2] = quad_num_entities[1] * _slices + quad_num_entities[2] * (_slices + 1);
        hexa_num_entities[3] = quad_num_entities[2] * _slices;
      }

      /**
       * \brief Tries to extrudes a quadrilateral mesh chart.
       *
       * \tparam SubChart_
       * The type of the quadrilateral sub-chart that is to be tested.
       *
       * \param[out] hexa_chart
       * A pointer to the hexa chart.
       *
       * \param[in] quad_chart
       * The chart that is to be extruded.
       *
       * \returns
       * \c true, if the chart was extruded, otherwise \c false.
       */
      template<typename SubChart_>
      bool try_extrude_chart(HexaChart*& hexa_chart, const QuadChart& quad_chart) const
      {
        // try down-cast to sub-chart
        const SubChart_* sub_chart = dynamic_cast<const SubChart_*>(&quad_chart);
        if(sub_chart == nullptr)
          return false;

        // success; create the extured hexa chart
        auto* the_chart = new Atlas::Extrude<HexaMesh, SubChart_>(new SubChart_(*sub_chart));
        hexa_chart = the_chart;

        // set our transformation
        the_chart->set_origin(_origin_x, _origin_y);
        the_chart->set_offset(_offset_x, _offset_y, _offset_z);
        the_chart->set_angles(_angle_y, _angle_p, _angle_r);
        return true;
      }

      /**
       * \brief Extrudes and returns a quadrilateral mesh chart.
       *
       * \param[in] quad_chart
       * The chart that is to be extruded.
       *
       * \returns
       * The extruded hexahedral chart.
       *
       * \attention
       * The returned object must be deleted by the caller or inserted into an atlas.
       */
      HexaChart* extrude_chart(const QuadChart& quad_chart) const
      {
        // The following ugly piece of code is a trial-&-error sequence using
        // dynamic-cast to figure out the actual type of the 2D chart.
        // That's not an elegant solution, but I don't have any other ideas...

        HexaChart* hexa_chart = nullptr;

        // Is it a Circle chart?
        if(this->template try_extrude_chart<Atlas::Circle<QuadMesh>>(hexa_chart, quad_chart))
          return hexa_chart;

        // Is it a Spline chart?
        if(this->template try_extrude_chart<Atlas::Spline<QuadMesh>>(hexa_chart, quad_chart))
          return hexa_chart;

        // If we come out here then we failed to extrude the 2D chart for some reason...
        throw InternalError(__func__, __FILE__, __LINE__, "Could not extrude 2D chart to 3D");
      }

      /**
       * \brief Extrudes an quadrilateral mesh atlas.
       *
       * This function extrudes all charts in the input quadrilateral mesh atlas
       * and inserts the extruded charts into the hexahedral mesh atlas.
       *
       * \param[in,out] hexa_atlas
       * The hexahedral mesh atlas that receives the extruded charts.
       *
       * \param[in] quad_atlas
       * The quadrilateral mesh atlas whose charts are to be extruded.
       */
      void extrude_atlas(HexaAtlas& hexa_atlas, const QuadAtlas& quad_atlas) const
      {
       // get chart names
        auto names = quad_atlas.get_chart_names();
        for(auto& name : names)
        {
          // get input chart
          const QuadChart* quad_chart = quad_atlas.find_mesh_chart(name);
          XASSERT(quad_chart != nullptr);

          // extrude chart
          HexaChart* hexa_chart = extrude_chart(*quad_chart);

          // add to hexa atlas
          hexa_atlas.add_mesh_chart(name, hexa_chart);
        }
      }

      /**
       * \brief Extrudes a quadrilateral vertex set.
       *
       * \param[in,out] hexa_vtx
       * The hexahedral mesh vertex set that receives the extruded vertex set.
       * Is assumed to be allocated to the correct size.
       *
       * \param[in] quad_vtx
       * The quadrilateral mesh vertex set that is to be extruded.
       */
      void extrude_vertex_set(HexaVertexSet& hexa_vtx, const QuadVertexSet& quad_vtx) const
      {
        // get number of source vertices
        const Index quad_num_verts = quad_vtx.get_num_vertices();

        // validate vertex count in hexa vertex set
        XASSERT(hexa_vtx.get_num_vertices() == (quad_num_verts * (_slices+1)));

        // compute rigid rotation matrix
        Tiny::Matrix<CoordType, 3, 3> rotation;
        rotation.set_rotation_3d(_angle_y, _angle_p, _angle_r);
        Tiny::Vector<CoordType, 3> tmp1, tmp2;

        // loop over all slices
        for(Index j(0); j <= _slices; ++j)
        {
          // compute slice index offset
          const Index vo = j * quad_num_verts;

          // compute slice z-coord
          const CoordType z = _z_min + (CoordType(j) / CoordType(_slices))*(_z_max - _z_min);

          // loop over all vertices
          for(Index i(0); i < quad_num_verts; ++i)
          {
            // apply our rigid transformation
            tmp1[0] = quad_vtx[i][0] - _origin_x;
            tmp1[1] = quad_vtx[i][1] - _origin_y;
            tmp1[2] = z;
            tmp2.set_mat_vec_mult(rotation, tmp1);
            hexa_vtx[vo+i][0] = tmp2[0] + _offset_x;
            hexa_vtx[vo+i][1] = tmp2[1] + _offset_y;
            hexa_vtx[vo+i][2] = tmp2[2] + _offset_z;
          }
        }
      }

      /**
       * \brief Extrudes a quadrilateral topology (aka "index set holder").
       *
       * \param[in,out] hexa_topo
       * The hexahedral topology to be computed. Must be allocated to the correct dimensions.
       *
       * \param[in] quad_topo
       * The quadrilateral topology to be extruded.
       */
      void extrude_topology(HexaTopology& hexa_topo, const QuadTopology& quad_topo) const
      {
        // fetch index sets
        const auto& quad_edge = quad_topo.template get_index_set<1,0>();
        const auto& quad_quad = quad_topo.template get_index_set<2,0>();
        auto& hexa_edge = hexa_topo.template get_index_set<1,0>();
        auto& hexa_quad = hexa_topo.template get_index_set<2,0>();
        auto& hexa_hexa = hexa_topo.template get_index_set<3,0>();

        const Index quad_num_entities[] =
        {
          quad_edge.get_index_bound(),
          quad_edge.get_num_entities(),
          quad_quad.get_num_entities()
        };

        // validate dimensions
        XASSERT(hexa_edge.get_index_bound()  == (quad_num_entities[0] * (_slices+1)));
        XASSERT(hexa_edge.get_num_entities() == (quad_num_entities[0] * _slices + quad_num_entities[1] * (_slices+1)));
        XASSERT(hexa_quad.get_num_entities() == (quad_num_entities[1] * _slices + quad_num_entities[2] * (_slices+1)));
        XASSERT(hexa_hexa.get_num_entities() == (quad_num_entities[2] * _slices));

        // loop over all slices (edges in x/y-direction)
        for(Index j(0); j <= _slices; ++j)
        {
          // compute slice vertex index offset
          const Index vo = j * quad_num_entities[0];

          // compute slice edge index offset
          const Index eo = j * quad_num_entities[1];

          // loop over all edges
          for(Index i(0); i < quad_num_entities[1]; ++i)
          {
            hexa_edge[eo+i][0] = vo + quad_edge[i][0];
            hexa_edge[eo+i][1] = vo + quad_edge[i][1];
          }
        }

        // loop over all slice intervals (edges in z-direction)
        for(Index j(0); j < _slices; ++j)
        {
          // compute slice vertex index offset
          const Index vo = j * quad_num_entities[0];

          // compute slice edge index offset
          const Index eo = j * quad_num_entities[0] + (_slices+1) * quad_num_entities[1];

          // loop over all vertices
          for(Index i(0); i < quad_num_entities[0]; ++i)
          {
            hexa_edge[eo+i][0] = vo + i;
            hexa_edge[eo+i][1] = vo + i + quad_num_entities[0];
          }
        }

        // loop over all slices (quads in x/y-direction)
        for(Index j(0); j <= _slices; ++j)
        {
          // compute slice vertex index offset
          const Index vo = j * quad_num_entities[0];

          // compute slice quad index offset
          const Index qo = j * quad_num_entities[2];

          // loop over all quads
          for(Index i(0); i < quad_num_entities[2]; ++i)
          {
            hexa_quad[qo+i][0] = vo + quad_quad[i][0];
            hexa_quad[qo+i][1] = vo + quad_quad[i][1];
            hexa_quad[qo+i][2] = vo + quad_quad[i][2];
            hexa_quad[qo+i][3] = vo + quad_quad[i][3];
          }
        }

        // loop over all slice intervals (quads in z-direction)
        for(Index j(0); j < _slices; ++j)
        {
          // compute slice vertex index offset
          const Index vo = j * quad_num_entities[0];

          // compute slice quad index offset
          const Index qo = j * quad_num_entities[1] + (_slices+1) * quad_num_entities[2];

          // loop over all edges
          for(Index i(0); i < quad_num_entities[1]; ++i)
          {
            hexa_quad[qo+i][0] = vo + quad_edge[i][0];
            hexa_quad[qo+i][1] = vo + quad_edge[i][1];
            hexa_quad[qo+i][2] = vo + quad_edge[i][0] + quad_num_entities[0];
            hexa_quad[qo+i][3] = vo + quad_edge[i][1] + quad_num_entities[0];
          }
        }

        // loop over all slices (quads in x/y-direction)
        for(Index j(0); j < _slices; ++j)
        {
          // compute slice vertex index offset
          const Index vo = j * quad_num_entities[0];

          // compute slice quad index offset
          const Index ho = j * quad_num_entities[2];

          // loop over all quads
          for(Index i(0); i < quad_num_entities[2]; ++i)
          {
            hexa_hexa[ho+i][0] = vo + quad_quad[i][0];
            hexa_hexa[ho+i][1] = vo + quad_quad[i][1];
            hexa_hexa[ho+i][2] = vo + quad_quad[i][2];
            hexa_hexa[ho+i][3] = vo + quad_quad[i][3];
            hexa_hexa[ho+i][4] = vo + quad_quad[i][0] + quad_num_entities[0];
            hexa_hexa[ho+i][5] = vo + quad_quad[i][1] + quad_num_entities[0];
            hexa_hexa[ho+i][6] = vo + quad_quad[i][2] + quad_num_entities[0];
            hexa_hexa[ho+i][7] = vo + quad_quad[i][3] + quad_num_entities[0];
          }
        }

        // build redundant index sets
        RedundantIndexSetBuilder<Shape::Hexahedron>::compute(hexa_topo);
      }

      /**
       * \brief Extrudes a quadrilateral mapping (aka "target set holder")
       *
       * \param[in,out] hexa_mapp
       * The hexahedral mapping to be computed. Is assumed to be allocated to the correct dimensions.
       *
       * \param[in] quad_mapp
       * The quadrilateral mapping to be extruded.
       *
       * \param[in] quad_parent
       * The quadrilateral parent mesh.
       */
      void extrude_mapping(HexaTrgSetHolder& hexa_mapp, const QuadTrgSetHolder& quad_mapp, const QuadMesh& quad_parent) const
      {
        // fetch target sets
        const auto& qv = quad_mapp.template get_target_set<0>();
        const auto& qe = quad_mapp.template get_target_set<1>();
        const auto& qq = quad_mapp.template get_target_set<2>();
        auto& hv = hexa_mapp.template get_target_set<0>();
        auto& he = hexa_mapp.template get_target_set<1>();
        auto& hq = hexa_mapp.template get_target_set<2>();
        auto& hh = hexa_mapp.template get_target_set<3>();

        const Index quad_num_entities[] =
        {
          quad_parent.get_num_entities(0),
          quad_parent.get_num_entities(1),
          quad_parent.get_num_entities(2)
        };

        const Index quad_part_num_entities[] =
        {
          qv.get_num_entities(),
          qe.get_num_entities(),
          qq.get_num_entities()
        };

        // compute vertices from vertices
        for(Index j(0); j <= _slices; ++j)
        {
          const Index vom = j * quad_num_entities[0];
          const Index vop = j * quad_part_num_entities[0];

          for(Index i(0); i < quad_part_num_entities[0]; ++i)
            hv[vop+i] = vom + qv[i];
        }

        // compute edges from edges
        for(Index j(0); j <= _slices; ++j)
        {
          const Index eom = j * quad_num_entities[1];
          const Index eop = j * quad_part_num_entities[1];

          for(Index i(0); i < quad_part_num_entities[1]; ++i)
            he[eop+i] = eom + qe[i];
        }

        // compute quads from quads
        for(Index j(0); j <= _slices; ++j)
        {
          const Index qom = j * quad_num_entities[2];
          const Index qop = j * quad_part_num_entities[2];

          for(Index i(0); i < quad_part_num_entities[2]; ++i)
            hq[qop+i] = qom + qq[i];
        }

        // compute edges from extruded vertices
        for(Index j(0); j < _slices; ++j)
        {
          const Index eom = j * quad_num_entities[0] + (_slices+1)*quad_num_entities[1];
          const Index eop = j * quad_part_num_entities[0] + (_slices+1)*quad_part_num_entities[1];

          for(Index i(0); i < quad_part_num_entities[0]; ++i)
            he[eop+i] = eom + qv[i];
        }

        // compute quads from extruded edges
        for(Index j(0); j < _slices; ++j)
        {
          const Index qom = j * quad_num_entities[1] + (_slices+1)*quad_num_entities[2];
          const Index qop = j * quad_part_num_entities[1] + (_slices+1)*quad_part_num_entities[2];

          for(Index i(0); i < quad_part_num_entities[1]; ++i)
            hq[qop+i] = qom + qe[i];
        }

        // compute hexas from extruded quads
        for(Index j(0); j < _slices; ++j)
        {
          const Index hom = j * quad_num_entities[2];
          const Index hop = j * quad_part_num_entities[2];

          for(Index i(0); i < quad_part_num_entities[2]; ++i)
            hh[hop+i] = hom + qq[i];
        }
      }

      void extrude_slice_mapping(HexaTrgSetHolder& hexa_mapp, const QuadMesh& quad_parent, const Index slice) const
      {
        XASSERT(slice <= _slices);

        const Index nv = quad_parent.get_num_entities(0);
        const Index ne = quad_parent.get_num_entities(1);
        const Index nq = quad_parent.get_num_entities(2);

        // compute offsets
        const Index ov = slice * nv;
        const Index oe = slice * ne;
        const Index oq = slice * nq;

        // fill vertex mapping
        auto& hv = hexa_mapp.template get_target_set<0>();
        for(Index i(0); i < nv; ++i)
          hv[i] = ov + i;

        // fill vertex mapping
        auto& he = hexa_mapp.template get_target_set<1>();
        for(Index i(0); i < ne; ++i)
          he[i] = oe + i;

        // fill vertex mapping
        auto& hq = hexa_mapp.template get_target_set<2>();
        for(Index i(0); i < nq; ++i)
          hq[i] = oq + i;
      }

      /**
       * \brief Extrudes a quadrilateral mesh attribute.
       *
       * \param[out] hexa_attrib
       * The hexahedral mesh attribute to be computed.
       *
       * \param[in] quad_attrib
       * The quadrilateral mesh attribute to be extruded.
       */
      void extrude_attribute(HexaAttrib*& hexa_attrib, const QuadAttrib& quad_attrib) const
      {
        const Index quad_num_verts = quad_attrib.get_num_vertices();
        const int quad_num_coords = quad_attrib.get_num_coords();

        // create attribute
        XASSERT(hexa_attrib == nullptr);
        hexa_attrib = new HexaAttrib(quad_num_verts * (_slices+1), quad_num_coords+1, 0);

        // loop over all slices
        for(Index j(0); j <= _slices; ++j)
        {
          // compute slice index offset
          const Index vo = j * quad_num_verts;

          // compute slice z-coord
          const CoordType z = _z_min + (CoordType(j) / CoordType(_slices))*(_z_max - _z_min);

          // loop over all vertices
          for(Index i(0); i < quad_num_verts; ++i)
          {
            for(int k(0); k < quad_num_coords; ++k)
            {
              (*hexa_attrib)[vo+i][k] = quad_attrib[i][k];
            }
            (*hexa_attrib)[vo+i][quad_num_coords] = z;
          }
        }
      }

      void extrude_partition(Partition& hexa_parti, const Partition& quad_parti) const
      {
        const Adjacency::Graph& quad_graph = quad_parti.get_patches();

        const Index num_patches = quad_parti.get_num_patches();
        const Index quad_num_elems = quad_parti.get_num_elements();

        // allocate adjacency graph
        Adjacency::Graph hexa_graph(quad_graph.get_num_nodes_domain(),
          _slices * quad_graph.get_num_nodes_image(), _slices * quad_graph.get_num_indices());

        // get arrays
        const Index* q_ptr = quad_graph.get_domain_ptr();
        const Index* q_idx = quad_graph.get_image_idx();
        Index* h_ptr = hexa_graph.get_domain_ptr();
        Index* h_idx = hexa_graph.get_image_idx();

        // loop over all patches
        h_ptr[0] = Index(0);
        for(Index i(0); i < num_patches; ++i)
        {
          Index hp = h_ptr[i];
          // loop over all elements in patch
          for(Index j(q_ptr[i]); j < q_ptr[i+1]; ++j)
          {
            // add the element for each slice
            for(Index k(0); k < _slices; ++k, ++hp)
              h_idx[hp] = (k*quad_num_elems) + q_idx[j];
          }
          h_ptr[i+1] = hp;
        }

        // create the partition
        hexa_parti = Partition(std::move(hexa_graph), quad_parti.get_name(),
          quad_parti.get_priority(), quad_parti.get_level());
      }

      void extrude_partition_set(PartitionSet& hexa_part_set, const PartitionSet& quad_part_set) const
      {
        const auto& quad_parts = quad_part_set.get_partitions();
        for(auto it = quad_parts.begin(); it != quad_parts.end(); ++it)
        {
          Partition hexa_parti;
          extrude_partition(hexa_parti, *it);
          hexa_part_set.add_partition(std::move(hexa_parti));
        }
      }

      /**
       * \brief Extrudes a quadrilateral root mesh node
       *
       * This function extrudes the root mesh as well as all mesh parts of the root node.
       *
       * \attention
       * This function does \b not extrude the atlas of the root node.
       *
       * \param[in,out] hexa_root_node
       * The hexahedral root mesh node that is to be created.
       *
       * \param[in] quad_root_node
       * The quadrilateral root mesh node that is to be extruded.
       *
       * \param[in] hexa_atlas
       * The hexahedral atlas that has already been extruded from the corresponding quadrilateral atlas.
       * May be \c nullptr if there are no charts associated with the quadrilateral root mesh node.
       */
      void extrude_root_node(HexaRootNode& hexa_root_node, const QuadRootNode& quad_root_node, const HexaAtlas* hexa_atlas) const
      {
        // get root mesh
        const QuadMesh* quad_mesh = quad_root_node.get_mesh();
        XASSERT(quad_mesh != nullptr);

        // extrude root mesh
        {
          MeshExtruderFactory<QuadMesh> factory(*this, *quad_mesh);
          hexa_root_node.set_mesh(new HexaMesh(factory));
        }

        // extrude mesh parts
        auto part_names = quad_root_node.get_mesh_part_names();
        for(const auto& name : part_names)
        {
          // get the quad part
          const QuadPart* quad_part = quad_root_node.find_mesh_part(name);

          if(quad_part == nullptr)
            continue;

          // get the chart name
          String chart_name = quad_root_node.find_mesh_part_chart_name(name);
          const HexaChart* hexa_chart(nullptr);

          // try to find the chart
          if(!chart_name.empty() && (hexa_atlas != nullptr))
            hexa_chart = hexa_atlas->find_mesh_chart(chart_name);

          // create mesh part extruder factory
          MeshPartExtruderFactory<QuadMesh> factory(*this, *quad_part, *quad_mesh);

          // create hexa mesh part and insert into root node
          hexa_root_node.add_mesh_part(name, new HexaPart(factory), chart_name, hexa_chart);
        }

        // add zmin/zmax mesh-parts
        if(!_zmin_part_name.empty())
        {
          MeshPartSliceExtruderFactory<QuadMesh> factory(*this, *quad_mesh, _zmin_part_name, Index(0));
          hexa_root_node.add_mesh_part(_zmin_part_name, new HexaPart(factory));
        }
        if(!_zmax_part_name.empty())
        {
          MeshPartSliceExtruderFactory<QuadMesh> factory(*this, *quad_mesh, _zmax_part_name, _slices);
          hexa_root_node.add_mesh_part(_zmax_part_name, new HexaPart(factory));
        }
      }
    }; // class MeshExtruder<ConformalMesh<Hypercube<2>,...>>

    /**
     * \brief Mesh Extruder factory
     *
     * This class implements the Factory interface using a MeshExtruder object to
     * extrude a given 2D quadrilateral mesh to a 3D hexahedral mesh.
     *
     * \author Peter Zajac
     */
    template<typename Coord_>
    class MeshExtruderFactory<ConformalMesh<Shape::Hypercube<2>, 2, 2, Coord_>> :
      public Factory<ConformalMesh<Shape::Hypercube<3>, 3, 3, Coord_>>
    {
    public:
      typedef ConformalMesh<Shape::Hypercube<2>, 2, 2, Coord_> QuadMesh;
      typedef ConformalMesh<Shape::Hypercube<3>, 3, 3, Coord_> HexaMesh;
      typedef MeshExtruder<QuadMesh> MeshExtruderType;

      typedef typename HexaMesh::VertexSetType VertexSetType;
      typedef typename HexaMesh::IndexSetHolderType IndexSetHolderType;

    protected:
      const MeshExtruderType& _mesh_extruder;
      const QuadMesh& _quad_mesh;
      Index _hexa_num_entities[4];

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mesh_extruder
       * The object that defines the mesh extrusion.
       *
       * \param[in] quad_mesh
       * The quadrilateral mesh that is to be extruded.
       */
      explicit MeshExtruderFactory(const MeshExtruderType& mesh_extruder, const QuadMesh& quad_mesh) :
        _mesh_extruder(mesh_extruder),
        _quad_mesh(quad_mesh)
      {
        const Index quad_num_entities[3] =
        {
          _quad_mesh.get_num_entities(0),
          _quad_mesh.get_num_entities(1),
          _quad_mesh.get_num_entities(2)
        };
        _mesh_extruder.extrude_num_entities(_hexa_num_entities, quad_num_entities);
      }

      virtual Index get_num_entities(int dim) override
      {
        XASSERT(dim < 4);
        return _hexa_num_entities[dim];
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        _mesh_extruder.extrude_vertex_set(vertex_set, _quad_mesh.get_vertex_set());
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
      {
        _mesh_extruder.extrude_topology(index_set_holder, _quad_mesh.get_index_set_holder());
      }
    }; // class MeshExtruderFactory<ConformalMesh<Shape::Hypercube<2>,...>>

    /**
     * \brief MeshPart Extruder factory
     *
     * This class implements the Factory interface using a MeshExtruder object to
     * extrude a given 2D quadrilateral mesh-part to a 3D hexahedral mesh-part.
     *
     * \author Peter Zajac
     */
    template<typename Coord_>
    class MeshPartExtruderFactory<ConformalMesh<Shape::Hypercube<2>, 2, 2, Coord_>> :
      public Factory<MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, 3, Coord_>>>
    {
    public:
      typedef ConformalMesh<Shape::Hypercube<2>, 2, 2, Coord_> QuadMesh;
      typedef ConformalMesh<Shape::Hypercube<3>, 3, 3, Coord_> HexaMesh;

      typedef MeshPart<QuadMesh> QuadPart;
      typedef MeshPart<HexaMesh> HexaPart;

      typedef MeshExtruder<QuadMesh> MeshExtruderType;

      typedef typename HexaPart::IndexSetHolderType IndexSetHolderType;
      typedef typename HexaPart::TargetSetHolderType TargetSetHolderType;
      typedef typename HexaPart::MeshAttributeContainer MeshAttributeContainer;
      typedef typename HexaPart::MeshAttributeType HexaAttribute;

    protected:
      const MeshExtruderType& _mesh_extruder;
      const QuadPart& _quad_part;
      const QuadMesh& _quad_parent;
      Index _hexa_num_entities[4];

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mesh_extruder
       * The object that defines the mesh extrusion.
       *
       * \param[in] quad_part
       * The quadrilateral mesh-part that is to be extruded.
       *
       * \param[in] quad_parent
       * The quadrilateral parent mesh of the to-be-extruded mesh part.
       */
      explicit MeshPartExtruderFactory(const MeshExtruderType& mesh_extruder, const QuadPart& quad_part, const QuadMesh& quad_parent) :
        _mesh_extruder(mesh_extruder),
        _quad_part(quad_part),
        _quad_parent(quad_parent)
      {
        const Index quad_num_entities[3] =
        {
          _quad_part.get_num_entities(0),
          _quad_part.get_num_entities(1),
          _quad_part.get_num_entities(2)
        };
        _mesh_extruder.extrude_num_entities(_hexa_num_entities, quad_num_entities);
      }

      virtual Index get_num_entities(int dim) override
      {
        XASSERT(dim < 4);
        return _hexa_num_entities[dim];
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _mesh_extruder.extrude_mapping(target_set_holder, _quad_part.get_target_set_holder(), _quad_parent);
      }

      virtual void fill_index_sets(IndexSetHolderType*& index_set_holder) override
      {
        if(_quad_part.has_topology())
        {
          index_set_holder = new IndexSetHolderType(_hexa_num_entities);
          _mesh_extruder.extrude_topology(*index_set_holder, *_quad_part.get_topology());
        }
      }

      virtual void fill_attribute_sets(MeshAttributeContainer& attribute_container) override
      {
        // extrude attributes
        for(const auto& quad_attrib : _quad_part.get_mesh_attributes())
        {
          HexaAttribute* hexa_attrib = nullptr;
          _mesh_extruder.extrude_attribute(hexa_attrib, *(quad_attrib.second));
          attribute_container.insert(std::make_pair(quad_attrib.first, hexa_attrib));
        }
      }
    }; // class MeshPartExtruderFactory<ConformalMesh<Shape::Hypercube<2>,...>>

    /**
     * \brief MeshPart Extruder factory
     *
     * This class implements the Factory interface using a MeshExtruder object to
     * extrude a given 2D quadrilateral mesh-part to a 3D hexahedral mesh-part.
     *
     * \author Peter Zajac
     */
    template<typename Coord_>
    class MeshPartSliceExtruderFactory<ConformalMesh<Shape::Hypercube<2>, 2, 2, Coord_>> :
      public Factory<MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, 3, Coord_>>>
    {
    public:
      typedef ConformalMesh<Shape::Hypercube<2>, 2, 2, Coord_> QuadMesh;
      typedef ConformalMesh<Shape::Hypercube<3>, 3, 3, Coord_> HexaMesh;

      typedef MeshPart<QuadMesh> QuadPart;
      typedef MeshPart<HexaMesh> HexaPart;

      typedef MeshExtruder<QuadMesh> MeshExtruderType;

      typedef typename HexaPart::IndexSetHolderType IndexSetHolderType;
      typedef typename HexaPart::TargetSetHolderType TargetSetHolderType;
      typedef typename HexaPart::MeshAttributeContainer MeshAttributeContainer;
      typedef typename HexaPart::MeshAttributeType HexaAttribute;

    protected:
      const MeshExtruderType& _mesh_extruder;
      const QuadMesh& _quad_parent;
      Index _slice;
      String _name;

    public:
      explicit MeshPartSliceExtruderFactory(const MeshExtruderType& mesh_extruder, const QuadMesh& quad_parent, const String& name, const Index slice) :
        _mesh_extruder(mesh_extruder),
        _quad_parent(quad_parent),
        _slice(slice),
        _name(name)
      {
      }

      virtual Index get_num_entities(int dim) override
      {
        XASSERT(dim < 4);
        return (dim < 3) ? _quad_parent.get_num_entities(dim) : Index(0);
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _mesh_extruder.extrude_slice_mapping(target_set_holder, _quad_parent, _slice);
      }

      virtual void fill_index_sets(IndexSetHolderType*& index_set_holder) override
      {
        index_set_holder = nullptr;
      }

      virtual void fill_attribute_sets(MeshAttributeContainer& /*attribute_set_holder*/) override
      {
      }
    }; // class MeshPartSliceExtruderFactory<ConformalMesh<Shape::Hypercube<2>,...>>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_MESH_EXTRUDER_HPP
