#pragma once
#ifndef KERNEL_BASE_MESH_2D_HPP
#define KERNEL_BASE_MESH_2D_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector
#include <typeinfo>  // for typeid()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_vertex.hpp>
#include <kernel/base_mesh_cell.hpp>
#include <kernel/base_mesh_cell_1d_edge.hpp>
#include <kernel/base_mesh_cell_2d_tri.hpp>
#include <kernel/base_mesh_cell_2d_quad.hpp>
#include <kernel/base_mesh_cell_3d_tetra.hpp>
#include <kernel/base_mesh_cell_3d_hexa.hpp>

namespace FEAST
{
  namespace BaseMesh
  {
    /**
    * \brief 3D base mesh
    *
    * \author Dominik Goeddeke
    * \author Hilmar Wobker
    */
// COMMENT_HILMAR: Um Code-Redundanz zu vermeiden, koennten wir ueberlegen, eine Klasse BaseMesh einzufuehren,
// die als Elternklasse fuer BaseMesh1D/2D/3D dient. Und/Oder (aehnlich wie bei Cell) dimension-abhaengige Interfaces
// definieren.
    template<unsigned char world_dim_>
    class BaseMesh3D
    {
      /// shortcuts various cell types to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;
      typedef Edge<3, world_dim_> Edge_;
      typedef Cell<1, 3, world_dim_> Cell_1D_;
      typedef Tri<3, world_dim_> Tri_;
      typedef Quad<3, world_dim_> Quad_;
      typedef Cell<2, 3, world_dim_> Cell_2D_;
      typedef Tetra<3, world_dim_> Tetra_;
      typedef Hexa<3, world_dim_> Hexa_;
      typedef Cell<3, 3, world_dim_> Cell_;

    private:

      /* *****************
      * member variables *
      *******************/
      /// array of vertices
      std::vector<Vertex_*> _vertices;

      /// array of edges
      std::vector<Cell_1D_*> _edges;

      /// array of faces
      std::vector<Cell_2D_*> _faces;

      /// array of cells
      std::vector<Cell_*> _cells;

      /* **********
      * functions *
      ************/
      /**
      * \brief templated function to remove vertices, edges and cells from the corresponding vector
      *
      * The item is swapped to the end of the list and then deleted and removed.
      * TODO: This code is replicated in all base_mesh_xD classes since there is no generic class to inherit from
      */
      template<typename T_>
      inline void remove(std::vector<T_>& v, T_ item)
      {
        assert(item->index() < v.size());
        v[item->index()] = v.back();
        v[item->index()]->set_index(item->index());
        v[v.size()-1] = item;
        item->set_index(v.size()-1);
        delete item;
        v.pop_back();
      }

    public:

      /* ***************************
      * constructors & destructors *
      *****************************/
      /// default CTOR, currently generates a test mesh
      BaseMesh3D()
      {
        // first, a bunch of vertices
        Vertex_* v = new Vertex_();
        v->set_coord(0, 0.0); v->set_coord(1, 0.0); v->set_coord(2, 0.0);
        add(v);

        v = new Vertex_();
        v->set_coord(0, 1.0); v->set_coord(1, 0.0); v->set_coord(2, 0.0);
        add(v);

        v = new Vertex_();
        v->set_coord(0, 1.0); v->set_coord(1, 0.0); v->set_coord(2, 1.0);
        add(v);

        v = new Vertex_();
        v->set_coord(0, 0.0); v->set_coord(1, 0.0); v->set_coord(2, 1.0);
        add(v);

        v = new Vertex_();
        v->set_coord(0, 0.0); v->set_coord(1, 1.0); v->set_coord(2, 0.0);
        add(v);

        v = new Vertex_();
        v->set_coord(0, 1.0); v->set_coord(1, 1.0); v->set_coord(2, 0.0);
        add(v);

        v = new Vertex_();
        v->set_coord(0, 1.0); v->set_coord(1, 1.0); v->set_coord(2, 1.0);
        add(v);

        v = new Vertex_();
        v->set_coord(0, 0.0); v->set_coord(1, 1.0); v->set_coord(2, 1.0);
        add(v);

        v = new Vertex_();
        v->set_coord(0, 0.5); v->set_coord(1, 2.0); v->set_coord(2, 0.5);
        add(v);

        // then, a bunch of edges
        // e0
        Edge_* e = new Edge_(_vertices[0],_vertices[1]);
        add(e);
        // e1
        e = new Edge_(_vertices[1],_vertices[2]);
        add(e);
        // e2
        e = new Edge_(_vertices[2],_vertices[3]);
        add(e);
        // e3
        e = new Edge_(_vertices[0],_vertices[3]); // deliberately permuted
        add(e);
        // e4
        e = new Edge_(_vertices[0],_vertices[4]);
        add(e);
        // e5
        e = new Edge_(_vertices[1],_vertices[5]);
        add(e);
        // e6
        e = new Edge_(_vertices[2],_vertices[6]);
        add(e);
        // e7
        e = new Edge_(_vertices[3],_vertices[7]);
        add(e);
        // e8
        e = new Edge_(_vertices[4],_vertices[5]);
        add(e);
        // e9
        e = new Edge_(_vertices[5],_vertices[6]);
        add(e);
        // e10
        e = new Edge_(_vertices[6],_vertices[7]);
        add(e);
        // e11
        e = new Edge_(_vertices[7],_vertices[4]);
        add(e);
        // e12
        e = new Edge_(_vertices[5],_vertices[7]);
        add(e);
        // e13
        e = new Edge_(_vertices[8],_vertices[4]);
        add(e);
        // e14
        e = new Edge_(_vertices[8],_vertices[5]);
        add(e);
        // e15
        e = new Edge_(_vertices[8],_vertices[6]);
        add(e);
        // e16
        e = new Edge_(_vertices[8],_vertices[7]);
        add(e);

        // now for faces (this is where the ordering nightmare starts)
        /*
        * Local ordering for tets:
        * Edge 0: 0,1
        * Edge 1: 1,2
        * Edge 2: 2,0
        * Edge 3: 0,3
        * Edge 4: 1,3
        * Edge 5: 2,3
        * Face 0: 0,1,2
        * Face 1: 0,3,1
        * Face 2: 1,3,2
        * Face 3: 2,3,0
        */
        // erste Haelfte des Hex-Deckels (nur in cell1, d.h. folgendes global2local-mapping:
        // vertex4 ist 0, vertex5 ist 1, vertex7 ist 2 und vertex8 ist 3)
        // im ersten Tet ist das face0
        Tri_* tri = new Tri_(_vertices[4], _vertices[5], _vertices[7], _edges[8], _edges[12], _edges[11]);
        add(tri);
        // Vorderseite des ersten Tets (lokal face 1)
        tri = new Tri_(_vertices[4], _vertices[8], _vertices[5], _edges[13], _edges[14], _edges[8]);
        add(tri);
        // face senkrecht auf Kante 12 (geteilt zwischen den beiden Tets)
        // im ersten tet lokal face2, im zweiten tet lokal face 3)
        tri = new Tri_(_vertices[5], _vertices[8], _vertices[7], _edges[14], _edges[16], _edges[12]);
        add(tri);
        // "linkes" face des ersten tets (lokal face 3)
        tri = new Tri_(_vertices[7], _vertices[8], _vertices[4], _edges[16], _edges[13], _edges[11]);
        add(tri);

        // zweite Haelfte des Hex-Deckels (nur in cell2, d.h. folgendes global2local-mapping:
        // vertex5 ist 0, vertex6 ist 1, vertex7 ist 2 und vertex8 ist 3)
        // im zweiten tet ist das lokal face0
        tri = new Tri_(_vertices[5], _vertices[6], _vertices[7], _edges[9], _edges[10], _edges[12]);
        add(tri);
        // "rechte Seite" des zweiten Tets (lokal face 1)
        tri = new Tri_(_vertices[5], _vertices[8], _vertices[6], _edges[14], _edges[15], _edges[9]);
        add(tri);
        // "Rueckseite" des zweiten tets (lokal face 2)
        tri = new Tri_(_vertices[6], _vertices[8], _vertices[7], _edges[15], _edges[16], _edges[10]);
        add(tri);

        /*
        * local ordering for hexas:
        * Face 0: 0,1,2,3
        * Face 1: 0,4,5,1
        * Face 2: 1,5,6,2
        * Face 3: 2,6,7,3
        * Face 4: 0,3,7,4
        * Face 5: 4,7,6,5
        *
        * Edge 0: 0,1
        * Edge 1: 1,2
        * Edge 2: 2,3
        * Edge 3: 3,0
        * Edge 4: 0,4
        * Edge 5: 1,5
        * Edge 6: 2,6
        * Edge 7: 3,7
        * Edge 8: 4,5
        * Edge 9: 5,6
        * Edge10: 6,7
        * Edge11: 7,4
        */
        Quad_* quad = new Quad_(_vertices[0], _vertices[1], _vertices[2], _vertices[3],
                                _edges[0], _edges[1], _edges[2], _edges[3]);
        add(quad);
        quad = new Quad_(_vertices[0], _vertices[4], _vertices[5], _vertices[1],
                         _edges[4], _edges[8], _edges[5], _edges[0]);
        add(quad);
        quad = new Quad_(_vertices[1], _vertices[5], _vertices[6], _vertices[2],
                         _edges[5], _edges[9], _edges[6], _edges[1]);
        add(quad);
        quad = new Quad_(_vertices[2], _vertices[6], _vertices[7], _vertices[3],
                         _edges[6], _edges[10], _edges[7], _edges[2]);
        add(quad);
        quad = new Quad_(_vertices[0], _vertices[3], _vertices[7], _vertices[4],
                         _edges[3], _edges[7], _edges[11], _edges[4]);
        add(quad);
        quad = new Quad_(_vertices[4], _vertices[7], _vertices[6], _vertices[5],
                         _edges[11], _edges[10], _edges[9], _edges[8]);
        add(quad);

        // finally, cells (let's hope I didn't fuck up the local numbering within the faces)
        Hexa_* hex = new Hexa_(_vertices[0], _vertices[1], _vertices[2], _vertices[3],
                               _vertices[4], _vertices[5], _vertices[6], _vertices[7],
                               _edges[0], _edges[1], _edges[2], _edges[3], _edges[4], _edges[5],
                               _edges[6], _edges[7], _edges[8], _edges[9], _edges[10], _edges[11],
                               _faces[7], _faces[8], _faces[9], _faces[10], _faces[11], _faces[12]);
        add(hex);

        Tetra_* tet = new Tetra_(_vertices[4], _vertices[5], _vertices[7], _vertices[8],
                                 _edges[8], _edges[12], _edges[11], _edges[13], _edges[14], _edges[16],
                                 _faces[0], _faces[1], _faces[2], _faces[3]);
        add(tet);
        tet = new Tetra_(_vertices[5], _vertices[6], _vertices[7], _vertices[8],
                         _edges[9], _edges[10], _edges[12], _edges[14], _edges[15], _edges[16],
                         _faces[4], _faces[5], _faces[6], _faces[2]);
        add(tet);

        // haha, not finished yet! Now that everything is there, add neighbourhood information
        // this is emulated file parser part 2
        // and this is a pain in the lower end of the back *because* local orderings are again causing headaches

//        // the hex has the two tets as neighbours at its "top" face, aka at its local face 5
//        _cells[0]->face_neighbours(5).push_back(_cells[1]);
//        _cells[0]->face_neighbours(5).push_back(_cells[2]);
//        // the first tet sees the hex at its bottom face (local face 0) and the other tet at its "right" face (local 2)
//        _cells[1]->face_neighbours(0).push_back(_cells[0]);
//        _cells[1]->face_neighbours(2).push_back(_cells[2]);
//        // and the second tet sees the hex at its bottom face (local 0) and the other tet at its "left" face (local 3)
//        _cells[2]->face_neighbours(0).push_back(_cells[0]);
//        _cells[2]->face_neighbours(3).push_back(_cells[1]);

      }

      /// default destructor
      ~BaseMesh3D()
      {
        // delete all cells and their associated information
        // (pop_back calls destructor of the element being removed, so do not use an iterator because fiddling about with
        // the std::vector invalidates it. Also, pop_back calls default DTOR of BMI*, so we have to manually call
        // delete here)
        while (!_cells.empty())
        {
          delete _cells.back();
          _cells.pop_back();
        }

        // delete all faces and their associated information
        while (!_faces.empty())
        {
          delete _faces.back();
          _faces.pop_back();
        }

        // delete all edges and their associated information
        while (!_edges.empty())
        {
          delete _edges.back();
          _edges.pop_back();
        }

        // delete all vertices and their associated information
        while (!_vertices.empty())
        {
          delete _vertices.back();
          _vertices.pop_back();
        }
      }


      /* *********************************
      * getters & setters & manipulators *
      ************************************/

      /// returns number of vertices in this mesh
      inline global_index_t num_vertices() const
      {
        return _vertices.size();
      }

      /// returns number of edges in this mesh
      inline global_index_t num_edges() const
      {
        // TODO: potentiell falsch, auch Kanten koennen inaktiv sein und duerfen dann beim Transfer zu den Rechenprozessen
        // nicht mitgezaehlt werden!
        return _edges.size();
      }

      /// returns number of faces in this mesh
      inline global_index_t num_faces() const
      {
        // TODO: potentiell falsch, auch Faces koennen inaktiv sein und duerfen dann beim Transfer zu den Rechenprozessen
        // nicht mitgezaehlt werden!
        return _faces.size();
      }

      /// returns number of cells in this mesh (this is potentially expensive)
      inline global_index_t num_cells() const
      {
        // TODO: nochmal implementieren, wenn inaktive immer schoen ans Ende geschoben werden und es einen index gibt,
        // der die Position des letzten aktiven merkt
        global_index_t counter = 0;
        for (global_index_t i(0) ; i < _cells.size() ; ++i)
        {
          if (_cells[i]->active())
          {
            counter++;
          }
        }
        return counter;
      }

      /// returns vertex at given index
      inline Vertex_* vertex(global_index_t const index)
      {
        assert(index < _vertices.size());
        return _vertices[index];
      }

      /// returns edge at given index
      inline Cell_1D_* edge(global_index_t const index)
      {
        assert(index < _edges.size());
        return _edges[index];
      }

      /// returns face at given index
      inline Cell_2D_* face(global_index_t const index)
      {
        assert(index < _faces.size());
        return _faces[index];
      }

      /// returns cell at given index
      inline Cell_* cell(global_index_t const index) const
      {
        assert(index < _cells.size());
        return _cells[index];
      }

      /// adds given vertex to base mesh and sets its index
      inline void add(Vertex_* v)
      {
        _vertices.push_back(v);
        v->set_index(_vertices.size()-1);
      }

      /// adds given edge to base mesh and sets its index
      inline void add(Cell_1D_* e)
      {
        _edges.push_back(e);
        e->set_index(_edges.size()-1);
      }

      /// adds given face to base mesh and sets its index
      inline void add(Cell_2D_* f)
      {
        _faces.push_back(f);
        f->set_index(_faces.size()-1);
      }

      /// adds given cell to base mesh and sets its index
      inline void add(Cell_* c)
      {
        _cells.push_back(c);
        c->set_index(_cells.size()-1);
      }

      /// deletes given vertex
      inline void remove(Vertex_* v)
      {
        remove<Vertex_*>(_vertices, v);
      }

      /// deletes given edge
      inline void remove(Cell_1D_* e)
      {
        remove<Cell_1D_*>(_edges, e);
      }

      /// deletes given face
      inline void remove(Cell_2D_* f)
      {
        remove<Cell_2D_*>(_faces, f);
      }

      /// deletes given cell
      inline void remove(Cell_* c)
      {
        remove<Cell_*>(_cells, c);
      }

      inline void add_created_items(SubdivisionData<2, 2, world_dim_>& subdiv_data)
      {
        for(unsigned int i(0) ; i < subdiv_data.created_vertices.size() ; ++i)
        {
          add(subdiv_data.created_vertices[i]);
        }
        for(unsigned int i(0) ; i < subdiv_data.created_edges.size() ; ++i)
        {
          add(subdiv_data.created_edges[i]);
        }
        for(unsigned int i(0) ; i < subdiv_data.created_faces.size() ; ++i)
        {
          add(subdiv_data.created_faces[i]);
        }
        for(unsigned int i(0) ; i < subdiv_data.created_cells.size() ; ++i)
        {
          add(subdiv_data.created_cells[i]);
        }
      }

      /**
      * \brief prints this BaseMesh to the given ostream.
      *
      * According to http://www.cplusplus.com/reference/iostream, ostream is the superclass for both stdout, stderr and
      * a (stream associated with) an arbitrary already opened ASCII file, so this routine can be used for logging
      * and and printf()-style debugging at the same time. Neat.
      *
      * \param[in] stream
      *            Stream to dump this BaseMesh into.
      */
      void print(std::ostream& stream)
      {
        stream << "---------------------------------------------------" << std::endl;
        stream << "|               DUMPING BASE MESH                  " << std::endl;
        stream << "---------------------------------------------------" << std::endl;
        stream << num_vertices() <<" Vertices:" << std::endl;
        for (unsigned int ivert(0) ; ivert < _vertices.size() ; ++ivert)
        {
          _vertices[ivert]->print(stream);
          stream << std::endl;
        }
        stream << num_edges() << " Edges" << std::endl;
        for (unsigned int iedge(0) ; iedge < _edges.size() ; ++iedge)
        {
          _edges[iedge]->print(stream);
          stream << std::endl;
        }
        stream << num_faces() << " Faces" << std::endl;
        for (unsigned int iface(0) ; iface < _faces.size() ; ++iface)
        {
          _faces[iface]->print(stream);
          stream << std::endl;
        }
        stream << num_cells() << " Cells" << std::endl;
        for (unsigned int icell(0) ; icell < _cells.size() ; ++icell)
        {
          _cells[icell]->print(stream);
          stream << std::endl;
        }
        stream << "---------------------------------------------------" << std::endl;
      }
    }; // class BaseMesh3D
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_2D_HPP
