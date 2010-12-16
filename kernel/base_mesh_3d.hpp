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
#include <kernel/base_mesh_cell_2d_quad.hpp>
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
      /// shortcut for type Vertex<world_dim_>
      typedef Vertex<world_dim_> Vertex_;

      /// shortcut for type Edge<2, world_dim_>
      typedef Edge<3, world_dim_> Edge_;

      /// shortcut for type Cell<1, 3, world_dim_>
      typedef Cell<1, 3, world_dim_> Cell_1D_;

      /// shortcut for type Quad<3, world_dim_>
      typedef Quad<3, world_dim_> Quad_;

      /// shortcut for type Cell<2, 3, world_dim_>
      typedef Cell<2, 3, world_dim_> Cell_2D_;

      /// shortcut for type Cell<3, 3, world_dim_>
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
        // Base mesh example consisting of one hexa, faces viewed from inside the hexa in direction (0, -1, 0)
        // bottom face (x3=0) [0,0,0]x[1,1,0]      top face (x3=1) [0,0,1]x[1,1,1]
        //  (0,1,0)    (1,1,0)                     (0,0,1)     (1,0,1)
        //    v3---e2---v2                            v4---e6---v7
        //    |          |                            |          |
        //   e3    f0   e1                           e7    f1   e5
        //    |          |                            |          |
        //    v0---e0---v1                            v5---e4---v6
        //  (0,0,0)    (1,0,0)                     (0,1,1)     (1,1,1)
        //
        // front face (x2=0) [0,0,0]x[1,0,1]      back face (x2=1) [0,1,0]x[1,1,1]
        //  (1,0,1)    (0,0,1)                     (0,1,1)     (1,1,1)
        //    v7---e6---v4                            v5---e4---v6
        //    |          |                            |          |
        //   e9    f2   e8                           e11   f3   e10
        //    |          |                            |          |
        //    v1---e0---v0                            v3---e2---v2
        //  (1,0,0)    (0,0,0)                     (0,1,0)     (1,1,0)
        //
        // right face (x1=0) [0,0,0]x[0,1,1]      left face (x2=1) [0,1,0]x[1,1,1]
        //  (1,0,1)    (0,0,1)                     (0,1,1)     (1,1,1)
        //    v4---e7---v5                            v6---e5---v7
        //    |          |                            |          |
        //   e8    f4   e11                          e10   f5   e9
        //    |          |                            |          |
        //    v0---e3---v3                            v2---e1---v1
        //  (1,0,0)    (0,0,0)                     (0,1,0)     (1,1,0)

        // create the 8 vertices
        // v0 = (0,0,0)
        Vertex_* v = new Vertex_();
        v->set_coord(0, 0.0);
        v->set_coord(1, 0.0);
        v->set_coord(2, 0.0);
        add(v);

        // v1 = (1,0,0)
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, 0.0);
        v->set_coord(2, 0.0);
        add(v);

        // v2 = (1,1,0)
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, 1.0);
        v->set_coord(2, 0.0);
        add(v);

        // v3 = (0,1,0)
        v = new Vertex_();
        v->set_coord(0, 0.0);
        v->set_coord(1, 1.0);
        v->set_coord(2, 0.0);
        add(v);

        // v4 = (0,0,1)
        v = new Vertex_();
        v->set_coord(0, 0.0);
        v->set_coord(1, 0.0);
        v->set_coord(2, 1.0);
        add(v);

        // v5 = (0,1,1)
        v = new Vertex_();
        v->set_coord(0, 0.0);
        v->set_coord(1, 1.0);
        v->set_coord(2, 1.0);
        add(v);

        // v6 = (1,1,1)
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, 1.0);
        v->set_coord(2, 1.0);
        add(v);

        // v7 = (1,0,1)
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, 0.0);
        v->set_coord(2, 1.0);
        add(v);

        // create the 12 edges
        // e0
        Edge_* e = new Edge_(_vertices[0], _vertices[1]);
        add(e);
        // e1
        e = new Edge_(_vertices[1], _vertices[2]);
        add(e);
        // e2
        e = new Edge_(_vertices[2], _vertices[3]);
        add(e);
        // e3
        e = new Edge_(_vertices[3], _vertices[0]);
        add(e);

        // e4
        e = new Edge_(_vertices[5], _vertices[6]);
        add(e);
        // e5
        e = new Edge_(_vertices[6], _vertices[7]);
        add(e);
        // e6
        e = new Edge_(_vertices[7], _vertices[4]);
        add(e);
        // e7
        e = new Edge_(_vertices[4], _vertices[5]);
        add(e);

        // e8
        e = new Edge_(_vertices[0], _vertices[4]);
        add(e);
        // e9
        e = new Edge_(_vertices[7], _vertices[1]);
        add(e);
        // e10
        e = new Edge_(_vertices[2], _vertices[6]);
        add(e);
        // e11
        e = new Edge_(_vertices[5], _vertices[3]);
        add(e);

        // quad face 0
        Quad_* f = new Quad_(_vertices[0], _vertices[1], _vertices[2], _vertices[3],
                             _edges[0], _edges[1], _edges[2], _edges[3]);
        add(f);
        // quad face 1
        f = new Quad_(_vertices[5], _vertices[6], _vertices[7], _vertices[4],
                      _edges[4], _edges[5], _edges[6], _edges[7]);
        add(f);
        // quad face 2
        f = new Quad_(_vertices[1], _vertices[0], _vertices[4], _vertices[7],
                      _edges[0], _edges[8], _edges[6], _edges[9]);
        add(f);
        // quad face 3
        f = new Quad_(_vertices[3], _vertices[2], _vertices[6], _vertices[5],
                      _edges[2], _edges[10], _edges[4], _edges[11]);
        add(f);
        // quad face 4
        f = new Quad_(_vertices[0], _vertices[3], _vertices[5], _vertices[4],
                      _edges[3], _edges[11], _edges[7], _edges[8]);
        add(f);
        // quad face 5
        f = new Quad_(_vertices[2], _vertices[1], _vertices[7], _vertices[6],
                      _edges[1], _edges[9], _edges[5], _edges[10]);
        add(f);

        // hex 0
        Hexa<3,3>* h = new Hexa<3,3>(_vertices[0], _vertices[1], _vertices[2], _vertices[3],
                           _vertices[4], _vertices[5], _vertices[6], _vertices[7],
                           _edges[0], _edges[1], _edges[2], _edges[3], _edges[4], _edges[5],
                           _edges[6], _edges[7], _edges[8], _edges[9], _edges[10], _edges[11],
                           _faces[0], _faces[1], _faces[2], _faces[3], _faces[4], _faces[5]);
        add(h);
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
