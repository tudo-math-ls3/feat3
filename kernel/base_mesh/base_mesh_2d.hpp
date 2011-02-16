#pragma once
#ifndef KERNEL_BASE_MESH_2D_HPP
#define KERNEL_BASE_MESH_2D_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh/vertex.hpp>
#include <kernel/base_mesh/cell.hpp>
#include <kernel/base_mesh/cell_1d_edge.hpp>
#include <kernel/base_mesh/cell_2d_quad.hpp>
#include <kernel/base_mesh/cell_2d_tri.hpp>
#include <kernel/graph.hpp>

namespace FEAST
{
  namespace BaseMesh
  {
    /**
    * \brief 2D base mesh
    *
    * \author Dominik Goeddeke
    * \author Hilmar Wobker
    */
// COMMENT_HILMAR: Um Code-Redundanz zu vermeiden, koennten wir ueberlegen, eine Klasse BaseMesh einzufuehren,
// die als Elternklasse fuer BaseMesh1D/2D/3D dient. Und/Oder (aehnlich wie bei Cell) dimension-abhaengige Interfaces
// definieren.
// COMMENT_HILMAR: Ich befuerchte, wir muessen auch das BaseMesh nach space_dimension_ templatisieren. Ansonsten müssen
// alle Klassen "oben drüber", die BaseMesh benutzen (z.B. Load Balancer), auch also ...1D, ...2D, ...3D implementiert
// werden.
    template<unsigned char world_dim_>
    class BaseMesh2D
    {
      /// shortcuts various cell types to save typing of template parameters
      typedef Vertex<world_dim_> Vertex_;
      typedef Edge<2, world_dim_> Edge_;
      typedef Cell<1, 2, world_dim_> Cell_1D_;
      typedef Tri<2, world_dim_> Tri_;
      typedef Quad<2, world_dim_> Quad_;
      typedef Cell<2, 2, world_dim_> Cell_;


    private:

      /* *****************
      * member variables *
      *******************/
      /// array of vertices
      std::vector<Vertex_*> _vertices;

      /// array of edges
      std::vector<Cell_1D_*> _edges;

      /// array of cells
      std::vector<Cell_*> _cells;

      /// graph describing the connectivity of the base mesh
// COMMENT_HILMAR: probably not really needed... just temporarily placed here, later we need a graph structure for
// the connectivity of process patches.
      Graph* _graph;

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
      inline void _remove(std::vector<T_>& v, T_ item)
      {
        assert(item->index() < v.size());
        v[item->index()] = v.back();
        v[item->index()]->set_index(item->index());
        v[v.size()-1] = item;
        item->set_index(v.size()-1);
        delete item;
        v.pop_back();
      }


      /// adds given vertex to base mesh and sets its index
      inline void _add(Vertex_* v)
      {
        _vertices.push_back(v);
        v->set_index(_vertices.size()-1);
      }


      /// adds given edge to base mesh and sets its index
      inline void _add(Cell_1D_* e)
      {
        _edges.push_back(e);
        e->set_index(_edges.size()-1);
      }


      /// adds given cell to base mesh and sets its index
      inline void _add(Cell_* c)
      {
        _cells.push_back(c);
        c->set_index(_cells.size()-1);
      }


      /// deletes given vertex
      inline void _remove(Vertex_* v)
      {
        _remove<Vertex_*>(_vertices, v);
      }


      /// deletes given edge
      inline void _remove(Cell_1D_* e)
      {
        _remove<Cell_1D_*>(_edges, e);
      }


      /// deletes given cell
      inline void _remove(Cell_* c)
      {
        _remove<Cell_*>(_cells, c);
      }


    public:

      // Parsing of the mesh files is outsourced to class FileParser2D. Make this class friend such that it has
      // access to private members of BaseMesh2D.
      friend class FileParser2D;

      /* ***************************
      * constructors & destructors *
      *****************************/

      /**
      * \brief  default CTOR for a 2D base mesh
      *
      * Ceates a base mesh basing on the provided mesh file.
      *
      * \param[in] file_name
      * name of the mesh file
      */
      BaseMesh2D()
        : _vertices(nullptr),
          _edges(nullptr),
          _cells(nullptr),
          _graph(nullptr)
      {
      }

/*
      /// test CTOR, generates a hard-wired test mesh
      BaseMesh2D()
      {
        // Base mesh example consisting of three quads and two tris:
        //    v0---e0---v1---e1---v2 \.
        //    |          |         |    \.
        //   e2    c0   e3   c1   e4  c2  \ e5
        //    |          |         |         \.
        //    v3---e6---v4---e7---v5----e8---v6
        //                      /  |         |
        //                  e9/ c3 e10  c4  e11
        //                /        |         |
        //              v7---e12--v8---e13---v9

        // create the ten vertices
        // v0
        Vertex_* v = new Vertex_();
        v->set_coord(0, 0.0);
        v->set_coord(1, 1.0);
        _add(v);

        // v1
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, 1.0);
        _add(v);

        // v2
        v = new Vertex_();
        v->set_coord(0, 2.0);
        v->set_coord(1, 1.0);
        _add(v);

        // v3
        v = new Vertex_();
        v->set_coord(0, 0.0);
        v->set_coord(1, 0.0);
        _add(v);

        // v4
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, 0.0);
        _add(v);

        // v5
        v = new Vertex_();
        v->set_coord(0, 2.0);
        v->set_coord(1, 0.0);
        _add(v);

        // v6
        v = new Vertex_();
        v->set_coord(0, 3.0);
        v->set_coord(1, 0.0);
        _add(v);

        // v7
        v = new Vertex_();
        v->set_coord(0, 1.0);
        v->set_coord(1, -1.0);
        _add(v);

        // v8
        v = new Vertex_();
        v->set_coord(0, 2.0);
        v->set_coord(1, -1.0);
        _add(v);

        // v9
        v = new Vertex_();
        v->set_coord(0, 3.0);
        v->set_coord(1, -1.0);
        _add(v);


        // create the 14 edges
        // (just to ease manual sanity checks, always use the vertex of smaller global index as start vertex)
        // e0
        Edge_* e = new Edge_(_vertices[0], _vertices[1], 0);
        _add(e);
        // e1
        e = new Edge_(_vertices[1], _vertices[2], 0);
        _add(e);
        // e2
        e = new Edge_(_vertices[0], _vertices[3], 0);
        _add(e);
        // e3
        e = new Edge_(_vertices[1], _vertices[4], 0);
        _add(e);
        // e4
        e = new Edge_(_vertices[2], _vertices[5], 0);
        _add(e);
        // e5
        e = new Edge_(_vertices[2], _vertices[6], 0);
        _add(e);
        // e6
        e = new Edge_(_vertices[3], _vertices[4], 0);
        _add(e);
        // e7
        e = new Edge_(_vertices[4], _vertices[5], 0);
        _add(e);
        // e8
        e = new Edge_(_vertices[5], _vertices[6], 0);
        _add(e);
        // e9
        e = new Edge_(_vertices[5], _vertices[7], 0);
        _add(e);
        // e10
        e = new Edge_(_vertices[5], _vertices[8], 0);
        _add(e);
        // e11
        e = new Edge_(_vertices[6], _vertices[9], 0);
        _add(e);
        // e12
        e = new Edge_(_vertices[7], _vertices[8], 0);
        _add(e);
        // e13
        e = new Edge_(_vertices[8], _vertices[9], 0);
        _add(e);

        // create quad cell c0
        Quad_* quad =
          new Quad_(_vertices[3], _vertices[4], _vertices[0], _vertices[1],
                    _edges[6], _edges[0], _edges[2], _edges[3], 0);
        _add(quad);

        // create quad cell c1
        quad = new Quad_(_vertices[4], _vertices[5], _vertices[1], _vertices[2],
                         _edges[7], _edges[1], _edges[3], _edges[4], 0);
        _add(quad);

        // create tri cell c2
        Tri_* tri = new Tri_(_vertices[5], _vertices[6], _vertices[2],
                             _edges[8], _edges[5], _edges[4], 0);
        _add(tri);
        // create tri cell c3
        tri = new Tri_(_vertices[7], _vertices[8], _vertices[5],
                       _edges[12], _edges[10], _edges[9], 0);
        _add(tri);

        // create quad cell c4
        quad = new Quad_(_vertices[8], _vertices[9], _vertices[5], _vertices[6],
                         _edges[13], _edges[8], _edges[10], _edges[11], 0);
        _add(quad);

        // set neighbourhood information (emulated file parser part 2)
        _cells[0]->add_neighbour(SDIM_EDGE, 3, _cells[1]);

        _cells[1]->add_neighbour(SDIM_EDGE, 2, _cells[0]);
        _cells[1]->add_neighbour(SDIM_EDGE, 3, _cells[2]);
        _cells[1]->add_neighbour(SDIM_VERTEX, 1, _cells[4]);
        _cells[1]->add_neighbour(SDIM_VERTEX, 1, _cells[3]);

        _cells[2]->add_neighbour(SDIM_EDGE, 2, _cells[1]);
        _cells[2]->add_neighbour(SDIM_EDGE, 0, _cells[4]);
        _cells[2]->add_neighbour(SDIM_VERTEX, 0, _cells[3]);

        _cells[3]->add_neighbour(SDIM_EDGE, 1, _cells[4]);
        _cells[3]->add_neighbour(SDIM_VERTEX, 2, _cells[2]);
        _cells[3]->add_neighbour(SDIM_VERTEX, 2, _cells[1]);

        _cells[4]->add_neighbour(SDIM_EDGE, 1, _cells[2]);
        _cells[4]->add_neighbour(SDIM_EDGE, 2, _cells[3]);
        _cells[4]->add_neighbour(SDIM_VERTEX, 2, _cells[1]);
      }
*/

      /// default destructor
      ~BaseMesh2D()
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


      /// returns number of edges in this mesh (including inactive ones)
      inline global_index_t num_edges() const
      {
        // TODO: potentiell falsch, auch Kanten koennen inaktiv sein und duerfen dann beim Transfer zu den
        // Rechenprozessen nicht mitgezaehlt werden!
        return _edges.size();
      }


      /// returns number of cells in this mesh (including inactive ones)
      inline global_index_t num_cells() const
      {
        return _cells.size();
      }


      /// returns number of active cells in this mesh (this is potentially expensive)
      inline global_index_t num_active_cells() const
      {
        // TODO: nochmal implementieren, wenn inaktive immer schoen ans Ende geschoben werden und es einen index gibt,
        // der die Position des letzten aktiven merkt
        global_index_t counter = 0;
        for(global_index_t i(0) ; i < num_cells() ; ++i)
        {
          if(cell(i)->active())
          {
            counter++;
          }
        }
        return counter;
      }


      /// returns vertex at given index
      inline Vertex_* vertex(global_index_t const index)
      {
        assert(index < num_vertices());
        return _vertices[index];
      }


      /// returns edge at given index
      inline Cell_1D_* edge(global_index_t const index)
      {
        assert(index < num_edges());
        return _edges[index];
      }


      /// returns cell at given index
      inline Cell_* cell(global_index_t const index) const
      {
        assert(index < num_cells());
        return _cells[index];
      }


      inline void add_created_items(SubdivisionData<2, 2, world_dim_>* subdiv_data)
      {
        for(unsigned int i(0) ; i < subdiv_data->created_vertices.size() ; ++i)
        {
          _add(subdiv_data->created_vertices[i]);
        }
        for(unsigned int i(0) ; i < subdiv_data->created_edges.size() ; ++i)
        {
          _add(subdiv_data->created_edges[i]);
        }
        for(unsigned int i(0) ; i < subdiv_data->created_cells.size() ; ++i)
        {
          _add(subdiv_data->created_cells[i]);
        }
      }


      /**
      * \brief sets numbers (not indices!) in all active cells (has to be called after subdivision of a cell etc.)
      *
      * \note Stupidly runs linearly through the vector of cells and overwrites existing numbers.
      * A different numbering strategy (try to keep existing numbers, fill gaps, ...) may be more clever... this will
      * be modified later.
      */
      inline void set_cell_numbers() const
      {
        global_index_t counter = 0;
        for(global_index_t i(0) ; i < num_cells() ; ++i)
        {
          if(cell(i)->active())
          {
            cell(i)->set_number(counter);
            counter++;
          }
          else
          {
            cell(i)->unset_number();
          }
        }
      }


      /// creates connectivity graph from information stored in this BaseMesh
// COMMENT_HILMAR: This is just an intermediate solution to artificially connect the base mesh to the load balancer.
// I.e., we assume here that each process receives exactly one BMC and that the connectivity graph relevant for the
// load balancer actually is the connectivity graph of the base mesh. Later, there will be the matrix patch layer and
// the process patch layer, which both have their own connectivity structure. The load balancer then actually needs the
// connectivity graph of the process patch layer. We also do not distinguish between edge and vertex neighbours here.
      void create_graph()
      {
        global_index_t n_active_cells = num_active_cells();
        // allocate index array
        unsigned int* index = new unsigned int[n_active_cells + 1];

        // graph data structure is filled by two sweeps through the cell list
        // first sweep: count neighbours of each cell, and maintain running total to fill index array
        // treat last index entry separately because cell array has one less entry than index array
        unsigned int num_neighbours_so_far = 0;
        for (global_index_t icell=0 ; icell < num_cells() ; ++icell)
// TODO: wir brauchen einen iterator fuer aktive Zellen!
        {
          global_index_t ipos(0);

          if(cell(icell)->active())
          {
            // set neighbours counted so far
            index[ipos] = num_neighbours_so_far;
//std::cout << "Setting index[" << ipos << "] = " << num_neighbours_so_far << std::endl;
            // count neighbours (here: edge neighbours and vertex neighbours)
            for (unsigned char sdim(0) ; sdim < 2 ; ++sdim)
            {
              num_neighbours_so_far += cell(icell)->num_neighbours_subdim((subdim)sdim);
            }
            ++ipos;
          }
        }
        index[n_active_cells] = num_neighbours_so_far;
//std::cout << "Setting index[" << n_active_cells << "] = " << num_neighbours_so_far << std::endl;

        // second sweep through data structure
        // second sweep adds actual neighbour cell numbers in the appropriate places into array neighbours
        // again, treat last loop instance separately
        unsigned int* neighbours = new unsigned int[index[n_active_cells]];
        num_neighbours_so_far = 0;
        for (global_index_t icell=0 ; icell < n_active_cells ; icell++)
// TODO: wir brauchen einen iterator fuer aktive Zellen!
        {
          Cell_* c = cell(icell);
          if (c->active())
          {
            for (unsigned char sdim(0) ; sdim < 2 ; ++sdim)
            {
              for(unsigned char item(0) ; item < c->num_subitems_per_subdim((subdim)sdim) ; ++item)
              {
                std::vector<Cell_*>& neigh_cells = c->neighbours_item((subdim)sdim, item);
                for(unsigned char k(0) ; k < neigh_cells.size() ; ++k)
                {
                  neighbours[num_neighbours_so_far] = neigh_cells[k]->number();
//std::cout << "neighbours[" << num_neighbours_so_far << "] = " << neighbours[num_neighbours_so_far] << std::endl;
                  ++num_neighbours_so_far;
                }
              }
            }
          }
        }

        // now, create graph object
        // temporarily, do not distinguish edge neighbours and diagonal neighbours
        if (_graph != nullptr)
        {
          delete _graph;
          _graph = nullptr;
        }
        _graph = new Graph(num_cells(), index, neighbours);
      }


      /// validates the base mesh and all of its cells
      void validate() const
      {
        std::cout << "Validating cells..." << std::endl;
        for(unsigned int icell(0) ; icell < _cells.size() ; ++icell)
        {
          cell(icell)->validate();
        }
        std::cout << "...done!" << std::endl;

        // COMMENT_HILMAR: add further validations...
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
        stream << _vertices.size() <<" vertices:" << std::endl;
        for(unsigned int ivert(0) ; ivert < _vertices.size() ; ++ivert)
        {
          _vertices[ivert]->print(stream);
          stream << std::endl;
        }
        stream << _edges.size() << " edges" << std::endl;
        for(unsigned int iedge(0) ; iedge < _edges.size() ; ++iedge)
        {
          _edges[iedge]->print(stream);
          stream << std::endl;
        }
        stream << _cells.size() << " cells" << std::endl;
        for(unsigned int icell(0) ; icell < _cells.size() ; ++icell)
        {
          _cells[icell]->print(stream);
          stream << std::endl;
        }
        stream << "---------------------------------------------------" << std::endl;
      }
    }; // class BaseMesh2D
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_2D_HPP
