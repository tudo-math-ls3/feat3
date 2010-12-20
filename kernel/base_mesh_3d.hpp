#pragma once
#ifndef KERNEL_BASE_MESH_2D_HPP
#define KERNEL_BASE_MESH_2D_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector

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
      * \brief templated function to remove vertices, edges, faces and cells from the corresponding vector
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
        /*
        * Test mesh consisting of 4 hexas:
        * View direction (0,0,-1) (from the front), (x,y) \in [0,2]x[0,2]:
        *   plane z=0 (back face):        plane z=1 (center):           plane z=2 (front):
        *
        * (0,2,0)             (2,2,0)   (0,2,1)             (2,2,1)   (0,2,2)             (2,2,2)
        *                                 13---26---14---27---15        05---08---06---09---07
        *                                  |         |         |         |         |         |
        *                                 23   13   24   14   25        05   01   06   02   07      y
        *                                  |         |         |         |         |         |      /\
        *             18---35---19        10---21---11---22---12        02---03---03---04---04      |
        *              |         |                   |         |                   |         |      |
        *             33   17   34                  19   12   20                  01   00   02      |
        *              |         |                   |         |                   |         |       --------> x
        *             16---32---17                  08---18---09                  00---00---01
        * (0,0,0)             (2,0,0)   (0,0,1)             (2,0,1)   (0,0,2)             (2,0,2)
        *
        * View direction (-1,0,0) (from the right), (y,z) \in [0,2]x[2,0]:
        *   plane x=0 (right):          plane x=1 (center):         plane x=2 (left):
        *
        * (0,2,2)             (0,2,0)   (1,2,2)             (1,2,0)   (2,2,2)             (2,2,0)
        *   05---15---13                  06---16---14                  07---17---15
        *    |         |                   |         |                   |         |
        *   05   07   23                  06   08   24                  07   09   25                          y
        *    |         |                   |         |                   |         |                          /\
        *   02---12---10                  03---13---11---30---18        04---14---12---31---19                |
        *                                  |         |         |         |         |         |                |
        *                                 01   03   19   16   33        02   04   20   17   34                |
        *                                  |         |         |         |         |         |     z <--------
        *                                 00---19---08---28---16        01---11---09---29---17
        * (0,0,2)             (0,0,0)   (1,0,2)             (1,0,0)   (2,0,2)             (2,0,0)
        *
        * View direction (0,-1,0) (from the top), (x,z) \in [0,2]x[2,0]:
        *   plane y=0 (bottom):         plane y=1 (center):         plane y=2 (top):
        *
        * (0,0,0)             (2,0,0)   (0,1,0)             (2,1,0)   (0,2,0)             (2,2,0)
        *             16---32---17                  18---35---19
        *              |         |                   |         |
        *             28   15   29                  30   18   31                                     --------> x
        *              |         |                   |         |                                    |
        *             08---18---09        10---21---11---22---12        13---26---14---27---15      |
        *              |         |         |         |         |         |         |         |      |
        *             10   02   11        12   05   13   06   14        15   10   16   11   17      \/
        *              |         |         |         |         |         |         |         |      z
        *             00---00---01        02---03---03---04---04        05---08---06---09---07
        * (0,0,2)             (2,0,2)   (0,1,2)             (2,1,2)   (0,2,2)             (2,2,2)
        *
        */
        // first, a bunch of vertices
        // v0 = (1,0,2)
        double c[3] = {1,0,2};
        Vertex_* v = new Vertex_(c);
        add(v);

        // v1 = (2,0,2)
        c[0] = 2;
        v = new Vertex_(c);
        add(v);

        // v2 = (0,1,2)
        c[0] = 0;
        c[1] = 1;
        v = new Vertex_(c);
        add(v);

        // v3 = (1,1,2)
        c[0] = 1;
        v = new Vertex_(c);
        add(v);

        // v4 = (2,1,2)
        c[0] = 2;
        v = new Vertex_(c);
        add(v);

        // v5 = (0,2,2)
        c[0] = 0;
        c[1] = 2;
        v = new Vertex_(c);
        add(v);

        // v6 = (1,2,2)
        c[0] = 1;
        v = new Vertex_(c);
        add(v);

        // v7 = (2,2,2)
        c[0] = 2;
        v = new Vertex_(c);
        add(v);

        // v8 = (1,0,1)
        c[0] = 1;
        c[1] = 0;
        c[2] = 1;
        v = new Vertex_(c);
        add(v);

        // v9 = (2,0,1)
        c[0] = 2;
        v = new Vertex_(c);
        add(v);

        // v10 = (0,1,1)
        c[0] = 0;
        c[1] = 1;
        v = new Vertex_(c);
        add(v);

        // v11 = (1,1,1)
        c[0] = 1;
        v = new Vertex_(c);
        add(v);

        // v12 = (2,1,1)
        c[0] = 2;
        v = new Vertex_(c);
        add(v);

        // v13 = (0,2,1)
        c[0] = 0;
        c[1] = 2;
        v = new Vertex_(c);
        add(v);

        // v14 = (1,2,1)
        c[0] = 1;
        v = new Vertex_(c);
        add(v);

        // v15 = (2,2,1)
        c[0] = 2;
        v = new Vertex_(c);
        add(v);

        // v16 = (1,0,0)
        c[0] = 1;
        c[1] = 0;
        c[2] = 0;
        v = new Vertex_(c);
        add(v);

        // v17 = (2,0,0)
        c[0] = 2;
        v = new Vertex_(c);
        add(v);

        // v18 = (1,1,0)
        c[0] = 1;
        c[1] = 1;
        v = new Vertex_(c);
        add(v);

        // v19 = (2,1,0)
        c[0] = 2;
        v = new Vertex_(c);
        add(v);


        // then, a bunch of edges
        // e0
        Edge_* e = new Edge_(_vertices[0],_vertices[1]);
        add(e);
        // e1
        e = new Edge_(_vertices[0],_vertices[3]);
        add(e);
        // e2
        e = new Edge_(_vertices[1],_vertices[4]);
        add(e);
        // e3
        e = new Edge_(_vertices[2],_vertices[3]);
        add(e);
        // e4
        e = new Edge_(_vertices[3],_vertices[4]);
        add(e);
        // e5
        e = new Edge_(_vertices[2],_vertices[5]);
        add(e);
        // e6
        e = new Edge_(_vertices[3],_vertices[6]);
        add(e);
        // e7
        e = new Edge_(_vertices[4],_vertices[7]);
        add(e);
        // e8
        e = new Edge_(_vertices[5],_vertices[6]);
        add(e);
        // e9
        e = new Edge_(_vertices[6],_vertices[7]);
        add(e);
        // e10
        e = new Edge_(_vertices[0],_vertices[8]);
        add(e);
        // e11
        e = new Edge_(_vertices[1],_vertices[9]);
        add(e);
        // e12
        e = new Edge_(_vertices[2],_vertices[10]);
        add(e);
        // e13
        e = new Edge_(_vertices[3],_vertices[11]);
        add(e);
        // e14
        e = new Edge_(_vertices[4],_vertices[12]);
        add(e);
        // e15
        e = new Edge_(_vertices[5],_vertices[13]);
        add(e);
        // e16
        e = new Edge_(_vertices[6],_vertices[14]);
        add(e);
        // e17
        e = new Edge_(_vertices[7],_vertices[15]);
        add(e);
        // e18
        e = new Edge_(_vertices[8],_vertices[9]);
        add(e);
        // e19
        e = new Edge_(_vertices[8],_vertices[11]);
        add(e);
        // e20
        e = new Edge_(_vertices[9],_vertices[12]);
        add(e);
        // e21
        e = new Edge_(_vertices[10],_vertices[11]);
        add(e);
        // e22
        e = new Edge_(_vertices[11],_vertices[12]);
        add(e);
        // e23
        e = new Edge_(_vertices[10],_vertices[13]);
        add(e);
        // e24
        e = new Edge_(_vertices[11],_vertices[14]);
        add(e);
        // e25
        e = new Edge_(_vertices[12],_vertices[15]);
        add(e);
        // e26
        e = new Edge_(_vertices[13],_vertices[14]);
        add(e);
        // e27
        e = new Edge_(_vertices[14],_vertices[15]);
        add(e);
        // e28
        e = new Edge_(_vertices[8],_vertices[16]);
        add(e);
        // e29
        e = new Edge_(_vertices[9],_vertices[17]);
        add(e);
        // e30
        e = new Edge_(_vertices[11],_vertices[18]);
        add(e);
        // e31
        e = new Edge_(_vertices[12],_vertices[19]);
        add(e);
        // e32
        e = new Edge_(_vertices[16],_vertices[17]);
        add(e);
        // e33
        e = new Edge_(_vertices[16],_vertices[18]);
        add(e);
        // e34
        e = new Edge_(_vertices[17],_vertices[19]);
        add(e);
        // e35
        e = new Edge_(_vertices[18],_vertices[19]);
        add(e);

        // now faces
        // f0
        Quad_* quad = new Quad_(_vertices[0], _vertices[1], _vertices[3], _vertices[4],
                                _edges[0], _edges[4], _edges[1], _edges[2]);
        add(quad);
        // f1
        quad = new Quad_(_vertices[2], _vertices[3], _vertices[5], _vertices[6],
                         _edges[3], _edges[8], _edges[5], _edges[6]);
        add(quad);
        // f2
        quad = new Quad_(_vertices[3], _vertices[4], _vertices[6], _vertices[7],
                         _edges[4], _edges[9], _edges[6], _edges[7]);
        add(quad);
        // f3
        quad = new Quad_(_vertices[0], _vertices[8], _vertices[3], _vertices[11],
                         _edges[10], _edges[13], _edges[1], _edges[19]);
        add(quad);
        // f4
        quad = new Quad_(_vertices[1], _vertices[9], _vertices[4], _vertices[12],
                         _edges[11], _edges[14], _edges[2], _edges[20]);
        add(quad);
        // f5
        quad = new Quad_(_vertices[2], _vertices[3], _vertices[10], _vertices[11],
                         _edges[3], _edges[21], _edges[12], _edges[13]);
        add(quad);
        // f6
        quad = new Quad_(_vertices[3], _vertices[4], _vertices[11], _vertices[12],
                         _edges[4], _edges[22], _edges[13], _edges[14]);
        add(quad);
        // f7
        quad = new Quad_(_vertices[2], _vertices[10], _vertices[5], _vertices[13],
                         _edges[12], _edges[15], _edges[5], _edges[23]);
        add(quad);
        // f8
        quad = new Quad_(_vertices[3], _vertices[11], _vertices[6], _vertices[14],
                         _edges[13], _edges[16], _edges[6], _edges[24]);
        add(quad);
        // f9
        quad = new Quad_(_vertices[4], _vertices[12], _vertices[7], _vertices[15],
                         _edges[14], _edges[17], _edges[7], _edges[25]);
        add(quad);
        // f10
        quad = new Quad_(_vertices[5], _vertices[6], _vertices[13], _vertices[14],
                         _edges[8], _edges[26], _edges[15], _edges[16]);
        add(quad);
        // f11
        quad = new Quad_(_vertices[6], _vertices[7], _vertices[14], _vertices[15],
                         _edges[9], _edges[27], _edges[16], _edges[17]);
        add(quad);
        // f12
        quad = new Quad_(_vertices[8], _vertices[9], _vertices[11], _vertices[12],
                         _edges[18], _edges[22], _edges[19], _edges[20]);
        add(quad);
        // f13
        quad = new Quad_(_vertices[10], _vertices[11], _vertices[13], _vertices[14],
                         _edges[21], _edges[26], _edges[23], _edges[24]);
        add(quad);
        // f14
        quad = new Quad_(_vertices[11], _vertices[12], _vertices[14], _vertices[15],
                         _edges[22], _edges[27], _edges[24], _edges[25]);
        add(quad);
        // f15
        quad = new Quad_(_vertices[8], _vertices[9], _vertices[16], _vertices[17],
                         _edges[18], _edges[32], _edges[28], _edges[29]);
        add(quad);
        // f16
        quad = new Quad_(_vertices[8], _vertices[16], _vertices[11], _vertices[18],
                         _edges[28], _edges[30], _edges[19], _edges[33]);
        add(quad);
        // f17
        quad = new Quad_(_vertices[9], _vertices[17], _vertices[12], _vertices[19],
                         _edges[29], _edges[31], _edges[20], _edges[34]);
        add(quad);
        // f18
        quad = new Quad_(_vertices[11], _vertices[12], _vertices[18], _vertices[19],
                         _edges[22], _edges[35], _edges[30], _edges[31]);
        add(quad);
        // f19
        quad = new Quad_(_vertices[16], _vertices[17], _vertices[18], _vertices[19],
                         _edges[32], _edges[35], _edges[33], _edges[34]);
        add(quad);

        // finally, cells
        // c0
        Hexa_* hex = new Hexa_(_vertices[0], _vertices[1], _vertices[8], _vertices[9],
                               _vertices[3], _vertices[4], _vertices[11], _vertices[12],
                               _edges[0], _edges[18], _edges[4], _edges[22], _edges[10], _edges[11],
                               _edges[13], _edges[14], _edges[1], _edges[2], _edges[19], _edges[20],
                               _faces[2], _faces[6], _faces[0], _faces[12], _faces[3], _faces[4]);
        add(hex);
        // c1
        hex = new Hexa_(_vertices[2], _vertices[3], _vertices[10], _vertices[11],
                        _vertices[5], _vertices[6], _vertices[13], _vertices[14],
                        _edges[3], _edges[21], _edges[8], _edges[26], _edges[12], _edges[13],
                        _edges[15], _edges[16], _edges[5], _edges[6], _edges[23], _edges[24],
                        _faces[5], _faces[10], _faces[1], _faces[13], _faces[7], _faces[8]);
        add(hex);
        // c2
        hex = new Hexa_(_vertices[3], _vertices[4], _vertices[11], _vertices[12],
                        _vertices[6], _vertices[7], _vertices[14], _vertices[15],
                        _edges[4], _edges[22], _edges[9], _edges[27], _edges[13], _edges[14],
                        _edges[16], _edges[17], _edges[6], _edges[7], _edges[24], _edges[25],
                        _faces[6], _faces[11], _faces[2], _faces[14], _faces[8], _faces[9]);
        add(hex);
        // c3
        hex = new Hexa_(_vertices[8], _vertices[9], _vertices[16], _vertices[17],
                        _vertices[11], _vertices[12], _vertices[18], _vertices[19],
                        _edges[18], _edges[32], _edges[22], _edges[35], _edges[28], _edges[29],
                        _edges[30], _edges[31], _edges[19], _edges[20], _edges[33], _edges[34],
                        _faces[15], _faces[18], _faces[12], _faces[19], _faces[16], _faces[17]);
        add(hex);


        // neighbourhood
        // face neighbours
        _cells[0]->add_neighbour(2, 1, _cells[2]);
        _cells[0]->add_neighbour(2, 3, _cells[3]);
        _cells[1]->add_neighbour(2, 5, _cells[2]);
        _cells[2]->add_neighbour(2, 0, _cells[0]);
        _cells[2]->add_neighbour(2, 4, _cells[1]);
        _cells[3]->add_neighbour(2, 2, _cells[0]);

        // edge neighbours
        _cells[0]->add_neighbour(1, 6, _cells[1]);
        _cells[1]->add_neighbour(1, 5, _cells[0]);
        _cells[2]->add_neighbour(1, 1, _cells[3]);
        _cells[3]->add_neighbour(1, 2, _cells[2]);

        // vertex neighbours
        _cells[1]->add_neighbour(0, 3, _cells[3]);
        _cells[3]->add_neighbour(0, 4, _cells[1]);
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
