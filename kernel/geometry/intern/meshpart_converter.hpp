#pragma once
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/atlas/bezier.hpp>
#include <kernel/adjacency/graph.hpp>


namespace FEAT::Geometry::Intern
{

  template<typename Shape_>
  class MeshPartToChartExporter
  {
  public:
    // typedef Mesh_ MeshType;
    typedef ConformalMesh<Shape_> MeshType;
    typedef MeshPart<MeshType> MeshPartType;
    static constexpr int dim = MeshType::shape_dim;
    typedef typename MeshType::CoordType DataType;


    static std::unique_ptr<Atlas::ChartBase<MeshType>> meshpart_to_bezier(const MeshType& mesh, const MeshPartType& meshpart)
    {
      // only works in 2d
      if constexpr(dim == 2)
      {
        const auto& vtx_set = mesh.get_vertex_set();
        const auto& edge_to_vert = mesh.template get_index_set<1, 0>();
        const auto& face_to_edge = mesh.template get_index_set<2, 1>();
        const auto& face_to_vert = mesh.template get_index_set<2, 0>();
        Adjacency::Graph vert_to_edge(Adjacency::RenderType::injectify_transpose, edge_to_vert);
        Adjacency::Graph edge_to_face(Adjacency::RenderType::injectify_transpose, face_to_edge);

        const auto& target_set = meshpart.template get_target_set<1>();

        std::vector<Index> is_target(mesh.get_num_entities(1));
        std::for_each(is_target.begin(), is_target.end(), [](Index& a){a = Index(0);});
        for(Index k = 0; k < target_set.get_num_entities(); ++k)
        {
          is_target[target_set[k]] = Index(1);
        }

        XASSERTM(target_set.get_num_entities() > 0u, "Trying to create bezier chart from non 1d meshpart");

        std::vector<Tiny::Vector<DataType, dim>> vtx_line;
        vtx_line.reserve(target_set.get_num_entities()+1);

        bool not_closed = false;
        Index cur_edge = target_set[0];
        Index prev_vert = ~Index(0);
        const auto* dom_ptr = vert_to_edge.get_domain_ptr();
        const auto* img_ptr = vert_to_edge.get_image_idx();
        // first of all, see if we start at a beginning or somewhere inside our mehspart
        for(int k = 0; k < edge_to_vert.get_num_indices(); ++k)
        {
          Index cur_vert = edge_to_vert[target_set[0]][k];
          bool is_boundary = true;
          for(Index l = dom_ptr[cur_vert]; l < dom_ptr[cur_vert+1]; ++l)
          {
            if(img_ptr[l] != cur_edge && is_target[img_ptr[l]] == Index(1))
              is_boundary = false;
          }
          if(is_boundary)
          {
            prev_vert = cur_vert;
            not_closed = true;
          }
        }

        if(prev_vert == ~Index(0))
        {
          prev_vert = edge_to_vert[target_set[0]][0];
        }

        const Index start_edge = target_set[0];
        vtx_line.push_back(vtx_set[prev_vert]);

        // track while running through our edges which orientation we are having
        Index orientation = Index(2);

        while(true)
        {
          Index cur_vert = Index(0);
          // first of all get the vertex not yet handled
          for(int k = 0; k < edge_to_vert.get_num_indices(); ++k)
          {
            if(edge_to_vert[cur_edge][k] != prev_vert)
            {
              cur_vert = edge_to_vert[cur_edge][k];
              break;
            }
          }

          vtx_line.push_back(vtx_set[cur_vert]);

          // TODO: we could do this smarter by comparing orientation of the cell and the referance egde, but this is more striaghtforward
          if(orientation != ~Index(0))
          {
            // compute mid vertex of our line segment
            auto mid_e = DataType(0.5) * (vtx_set[cur_vert] + vtx_set[prev_vert]);
            // compute mean value of our cell
            auto mid_c = Tiny::Vector<DataType, dim>();
            mid_c.format();
            if(edge_to_face.degree(cur_edge) != Index(1))
            {
              orientation = ~Index(0);
            }
            if(orientation != ~Index(0))
            {
              // get cell of our edge
              const auto& index_tuple = face_to_vert[edge_to_face.get_image_idx()[edge_to_face.get_domain_ptr()[cur_edge]]];
              for(int k = 0; k < index_tuple.num_indices; ++k)
              {
                mid_c.axpy(DataType(1)/DataType(index_tuple.num_indices), vtx_set[index_tuple[k]]);
              }

              const auto outer_normal = mid_e - mid_c;
              const auto line_segment = vtx_set[cur_vert] - vtx_set[prev_vert];

              // rotate line_segment by -90 degrees (i.e. to the right)
              const auto right_normal = Tiny::Vector<DataType, dim>{line_segment[1], -line_segment[0]};

              // if vectors are on the same halfplane, we have positive, i.e. counterclockwise orientation
              const Index new_orient = Tiny::dot(right_normal, outer_normal) > 0 ? Index(1) : Index(0);
              if(orientation == Index(2)) orientation = new_orient;
              XASSERTM(orientation == new_orient, "Non orientable");
            }
          }


          Index next_edge = ~Index(0);
          // find next edge
          for(Index l = dom_ptr[cur_vert]; l < dom_ptr[cur_vert+1]; ++l)
          {
            if(img_ptr[l] != cur_edge && is_target[img_ptr[l]] == Index(1))
            {
              next_edge = img_ptr[l];
              break;
            }
          }
          if(next_edge == ~Index(0))
          {
            not_closed = true;
            break;
          }

          if(next_edge == start_edge)
          {
            break;
          }

          cur_edge = next_edge;
          prev_vert = cur_vert;
        }

        std::unique_ptr<Atlas::Bezier<MeshType>> bezier = std::make_unique<Atlas::Bezier<MeshType>>(!not_closed, orientation == Index(0) ? 1 : -1);
        for(Index k = 0; k < vtx_line.size(); ++k)
        {
          bezier->push_vertex(vtx_line[k]);
          bezier->push_param(Tiny::Vector<DataType, 1>{DataType(k)});
        }

        return bezier;

      }
      else
      {
        XABORTM("Bezier chart only valid in 2d");
      }

      return std::unique_ptr<Atlas::ChartBase<MeshType>>{nullptr};
    }

  };

}