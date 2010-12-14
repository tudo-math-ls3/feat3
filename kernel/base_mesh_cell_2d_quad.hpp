#pragma once
#ifndef KERNEL_BASE_MESH_CELL_2D_QUAD_HPP
#define KERNEL_BASE_MESH_CELL_2D_QUAD_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_vertex.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_cell_1d_edge.hpp>

namespace FEAST
{
  namespace BaseMesh
  {

    /**
    * \brief 2D base mesh cell of type quad
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    * \author Peter Zajac
    */
// COMMENT_HILMAR: Um Code-Redundanz zu vermeiden, koennte wir ueberlegen, eine weitere Klasse Cell2D einzufuehren,
// die von Cell<2, space_dim_, world_dim_> erbt, und von der dann wieder um Quad und Tri erben. Darin koennte
// man zum Beispiel die Funktion _edge_has_correct_orientation() implementieren.
    template<
      unsigned char space_dim_,
      unsigned char world_dim_>
    class Quad
      : public Cell<2, space_dim_, world_dim_>
    {
      /// shortcuts for type Vertex<world_dim_>
      typedef Vertex<world_dim_> Vertex_;

      /// shortcuts for type Edge<space_dim_, world_dim_>
      typedef Edge<space_dim_, world_dim_> Edge_;

    private:
      /// vertices of the quad
      Vertex_* _vertices[4];

      /// edges of the quad
      Edge_* _edges[4];


      /// returns true when edge with local index iedge has the same orientation as the quad
      inline bool _edge_has_correct_orientation(unsigned char iedge)
      {
        // the orientation of the edge is correct (i.e. the same as that of the quad), when its start vertex within
        // the quad (i.e. the vertex with the same local index) is local vertex 0 within the edge structure
        return (vertex(iedge) == edge(iedge)->vertex(0));
      }


    public:
      /// CTOR
      Quad(Vertex_* v0, Vertex_* v1, Vertex_* v2, Vertex_* v3, Edge_* e0, Edge_* e1, Edge_* e2, Edge_* e3)
      {
        _vertices[0] = v0;
        _vertices[1] = v1;
        _vertices[2] = v2;
        _vertices[3] = v3;
        _edges[0] = e0;
        _edges[1] = e1;
        _edges[2] = e2;
        _edges[3] = e3;

        unsigned char temp[2] = {4,4};
        this->_init_neighbours(2, temp);
// TODO: Eigentlich haette ich das lieber in die Konstruktoren-Liste gepackt, also sowas in der Art:
//    : CellData<2, space_dim_, world_dim_>({4,4})
// (was nicht kompiliert). Wie kann man denn on-the-fly ein Array anlegen und durchreichen?
      }


      /// returns number of vertices
      inline unsigned char num_vertices() const
      {
        return 4;
      }


      /// returns vertex at given index
      inline Vertex_* vertex(unsigned char const index) const
      {
        assert(index < num_vertices());
        return _vertices[index];
      }


      /// returns number of edges
      inline unsigned char num_edges() const
      {
        return 4;
      }


      /// returns edge at given index
      inline Edge_* edge(unsigned char const index) const
      {
        assert(index < num_edges());
        return _edges[index];
      }


      /// subdivision routine splitting a quad and storing parent/child information
// COMMENT_HILMAR: this is currently hard-wired to splitting the quad into four quads. Later, this is parameterised
// via the information in the SubdivisionData object.
      inline void subdivide(SubdivisionData<2, space_dim_, world_dim_>& subdiv_data)
      {
        // assure that this cell has not been divided yet
        if(!this->active())
        {
          std::cerr << "Quad " << this->index() << " is already subdivided! Aborting program.";
          exit(1);
        }

        // clear all vectors of created entities in the SubdivisionData object
        subdiv_data.clear_created();

        /// vertices that this action creates and/or reuses
        Vertex_* new_vertices[5];

        /// edges that this action creates and/or reuses
        Cell<1, space_dim_, world_dim_>* new_edges[12];

        // local numbering (old and new)
        //         k2                                       e5     e4
        //   w3---------w2          -----v2------         -------------
        //   |           |          |     |     |       e6|    e10    |e3
        //   |           |          |  q3 | q2  |         |     |     |
        // k3|           |k1  ---> v3-----v4----v1        --e11----e9--
        //   |           |          |  q0 | q1  |       e7|     |     |e2
        //   |           |          |     |     |         |     e8    |
        //   w0---------w1          -----v0------         -------------
        //         k0                                        e0    e1

        // store old active-mask of each edge, because it gets overwritten once edges are getting split below
        bool old_edge_active_mask[4];

        SubdivisionData<1, space_dim_, world_dim_> subdiv_data_edge;

        // loop over all edges and split them eventually, creating new vertices and edges on the way
        for(unsigned char iedge(0) ; iedge < 4 ; ++iedge)
        {
          // if edge has no children, create them
          if (edge(iedge)->active())
          {
            // store old active mask
            old_edge_active_mask[iedge] = true;
            // create new vertex

            // subdivide edge
            edge(iedge)->subdivide(subdiv_data_edge);

            // add the created vertices/edges to the vector of created vertices/edges
            subdiv_data.created_vertices.push_back(subdiv_data_edge.created_vertex);

            // COMMENT_HILMAR: Not sure whether the order plays a role here... to be on the safe side, order them
            // according to the array new_edges[].
            if(_edge_has_correct_orientation(iedge))
            {
              subdiv_data.created_edges.push_back(subdiv_data_edge.created_cells[0]);
              subdiv_data.created_edges.push_back(subdiv_data_edge.created_cells[1]);
            }
            else
            {
              subdiv_data.created_edges.push_back(subdiv_data_edge.created_cells[1]);
              subdiv_data.created_edges.push_back(subdiv_data_edge.created_cells[0]);
            }
          }
          else // edge has children, reuse them
          {
            // store old active mask
            old_edge_active_mask[iedge] = false;
          }

          // add new vertex to array of new vertices
          new_vertices[iedge] = edge(iedge)->child(0)->vertex(1);

          // add new edges to array of new edges, respect the orientation of the edge
          if(_edge_has_correct_orientation(iedge))
          {
            // if the edge has the orientation of the quad, then child 0 is the first edge
            new_edges[2*iedge]   = edge(iedge)->child(0);
            new_edges[2*iedge+1] = edge(iedge)->child(1);
          }
          else
          {
            // if the edge does not have the orientation of the quad, then child 1 is the first edge.
            new_edges[2*iedge]   = edge(iedge)->child(1);
            new_edges[2*iedge+1] = edge(iedge)->child(0);
          }
// COMMENT_HILMAR: Beachte, dass wir auch in dem Fall, dass eine Kante schon Kinder hatte, *nicht* annehmen koennen,
// dass sie vom Nachbarelement erstellt worden ist. Beispiel: 2 benachbarte Zellen C0 und C1. C0 zerteilt sich, hat
// somit die gemeinsame Kante e zerteilt, sodass die Reihenfolge der Kinder der lokalen Orientierung von C0 entspricht.
// Irgendwann später zerteilt sich C1 und nutzt aus, dass die Kinder der Kante e ja von C0 angelegt worden sein
// muessen und dreht deren Reihenfolge also einfach um. Hier stimmt die Annahme also noch.
// Jetzt vergroebert sich C0 wieder, Kante e bleibt bestehen, da sie ja noch von C1 benutzt wird. Dann irgendwann
// verfeinert sich C0 wieder: Nun ist es aber falsch anzunehmen, dass der Nachbar (C1) die Kante angelegt haben muss!
// In Wirklichkeit war C0 es selbst, hat das aber inzwischen "vergessen".
        } // for(unsigned char iedge(0) ; iedge < 4 ; ++iedge)

        // create new midpoint and its incident edges (these are always new, have no children and cannot be reused)
        double x1 = new_vertices[0]->coords(0);
        double y1 = new_vertices[0]->coords(1);
        double x2 = new_vertices[2]->coords(0);
        double y2 = new_vertices[2]->coords(1);
        double x3 = new_vertices[1]->coords(0);
        double y3 = new_vertices[1]->coords(1);
        double x4 = new_vertices[3]->coords(0);
        double y4 = new_vertices[3]->coords(1);
        double p[2];
        // TODO factor out common subexpressions
        p[0] = ( (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) );
        p[1] = ( (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) );
        new_vertices[4] = new Vertex<world_dim_>(p);

        subdiv_data.created_vertices.push_back(new_vertices[4]);

        for (unsigned char i(0) ; i < 4 ; ++i)
        {
          new_edges[i+8] = new Edge_(new_vertices[i], new_vertices[4]);
          subdiv_data.created_edges.push_back(new_edges[i+8]);
        }

        // set number of children to 4
        this->_set_num_children(4);

        // finally, create and add new quads
        // TODO: get rid of manual downcasts here!
        _set_child(0, new Quad(vertex(0), new_vertices[0], new_vertices[4], new_vertices[3],
                               (Edge_*)new_edges[0], (Edge_*)new_edges[8], (Edge_*)new_edges[11], (Edge_*)new_edges[7]));
        _set_child(1, new Quad(new_vertices[0], vertex(1), new_vertices[1], new_vertices[4],
                               (Edge_*)new_edges[1], (Edge_*)new_edges[2], (Edge_*)new_edges[9], (Edge_*)new_edges[8]));
        _set_child(2, new Quad(new_vertices[4], new_vertices[1], vertex(2), new_vertices[2],
                               (Edge_*)new_edges[9], (Edge_*)new_edges[3], (Edge_*)new_edges[4], (Edge_*)new_edges[10]));
        _set_child(3, new Quad(new_vertices[3], new_vertices[4], new_vertices[2], vertex(3),
                               (Edge_*)new_edges[11], (Edge_*)new_edges[10], (Edge_*)new_edges[5], (Edge_*)new_edges[6]));

        for (unsigned char i(0) ; i < 4 ; ++i)
        {
          this->child(i)->set_parent(this);
          subdiv_data.created_cells.push_back(this->child(i));
        }
      } // subdivide()


      /// print information about this quad
      inline void print(std::ostream& stream)
      {
        stream << "Quad";
        Item::print(stream);
        stream << ": [";

        for(int i(0) ; i < 4 ; ++i)
        {
          stream << "E" << _edges[i]->index();
          if(_edge_has_correct_orientation(i))
          {
            stream << "(+)";
          }
          else
          {
            stream << "(-)";
          }
          if(i<3)
          {
            stream << ", ";
          }
          else
          {
            stream << "]";
          }
        }
        Cell<2, space_dim_, world_dim_>::print_history(stream);
        CellData<2, space_dim_, world_dim_>::print(stream);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_CELL_2D_QUAD_HPP
