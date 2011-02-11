#pragma once
#ifndef KERNEL_BASE_MESH_HPP
#define KERNEL_BASE_MESH_HPP 1

// includes, system
#include <iostream>
#include <vector>

// includes, Feast
#include <kernel/graph.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/file_reader_ascii.hpp>

namespace FEAST
{

  /**
  * \brief abstract class describing a base mesh cell
  *
  * Concrete extensions are available for 2D and 3D etc.
  *
  * COMMENT_DOMINIK: This is <em>very</em> preliminary, and only serves as a container to be able to implement
  *                  mesh I/O. Actual data structures need to be changed later on.
  *
  * \author Dominik Goeddeke
  */
  class BaseMeshCell
  {

  protected:
    /* *****************
    * member variables *
    *******************/
    /// global number of this cell
    unsigned int _number;

    /// dimension of this cell
    unsigned int _dimension;

    /// type of this cell, TODO: Add proper types later, when they become necessary
    unsigned int _type;

    /// parallel block (= rank within process group) that this cell is supposed to reside on
    unsigned int _parallel_block;

    /// matrix block this cell is supposed to contribute to
    unsigned int _matrix_block;


  public:
    /* ******
    * enums *
    ********/
    /// ID for 2D cells, Q1 elements and tensor product structure
    static const unsigned int TYPE_2D_Q1_TP = 0;


    /* ***************************
    * constructors & destructors *
    *****************************/
    /**
    * \brief Default CTOR
    *
    * \param[in] number
    * Number of this cell. Uniqueness must be ensured by BaseMesh class.
    *
    * \param[in] dimension
    * Dimension of this cell (1,2,3).
    */
    BaseMeshCell(unsigned int number, unsigned int dimension)
      : _number(number),
        _dimension(dimension),
        _type(0),
        _parallel_block(0),
        _matrix_block(0)
    {
    }

    /**
    * \brief Set-all CTOR
    */
    BaseMeshCell(unsigned int number,
                 unsigned int dimension,
                 unsigned int type,
                 unsigned int parallel_block,
                 unsigned int matrix_block)
      : _number(number),
        _dimension(dimension),
        _type(type),
        _parallel_block(parallel_block),
        _matrix_block(matrix_block)
    {
    }

    // TODO: CTOR WITH NEIGHBOURS //DG, Oct19,2010
    // TODO: (SHALLOW) COPY CTOR //DG, Oct19,2010

    /// default DTOR (must be virtual)
    virtual ~BaseMeshCell()
    {
    }

    /* ******************
    * getters & setters *
    ********************/
    inline unsigned int number() const
    {
      return _number;
    }
    inline unsigned int dimension() const
    {
      return _dimension;
    }
    inline unsigned int type() const
    {
      return _type;
    }
    inline void set_type(const unsigned int type)
    {
      _type = type;
    }
    inline unsigned int parallel_block() const
    {
      return _parallel_block;
    }
    inline void set_parallel_block(const unsigned int parallel_block)
    {
      _parallel_block = parallel_block;
    }
    inline unsigned int matrix_block() const
    {
      return _matrix_block;
    }
    inline void set_matrix_block(const unsigned int matrix_block)
    {
      _matrix_block = matrix_block;
    }


    /* *************************
    * virtual member functions *
    ***************************/
    /// prints this cell into a single line in std::cout (note: this routine must be extended in derived classes)
    virtual inline void print()
    {
      std::cout << _number << ": ";
      switch (_type)
      {
        case TYPE_2D_Q1_TP:
          std::cout << "2D-Q1-TP";
          break;
        default:
          std::cout << "unknown";
      }
      std::cout << ", ParallelBlock " << _parallel_block << ", MatrixBlock " << _matrix_block;
    }


    /* ******************************
    * pure virtual member functions *
    ********************************/
    /// returns number of neighbours of given dimension (point neighbours are dimension-0 neighbours etc.)
    virtual unsigned int num_neighbours(const unsigned int dimension) const = 0;

    /// adds given cell to list of given-dimension neighbours
    virtual void add_neighbour(const unsigned int dimension, BaseMeshCell * neighbour) = 0;

    /**
    * \brief Returns given neighbour of given dimension
    *
    * TODO: Ich wollte das hier mit nem Iterator machen, aber ich habe es nicht hinbekommen. So muss man jetzt
    * aussenrum ne Schleife ueber alle Nachbarn der gewuenschten Dimension machen (Ende der Schleife kennt man
    * ueber num_neighbours()) damit man in diese Routine den entsprechenden Index stopfen kann.
    */
    virtual BaseMeshCell * get_neighbour(const unsigned int dimension, const unsigned int index) const = 0;

  }; // class BaseMeshCell


  /**
  * \brief Base mesh cell in two dimensions
  *
  * COMMENT_DOMINIK: This is <em>very</em> preliminary, and only serves as a container to be able to implement
  *                  grid I/O. Actual data structures need to be changed later on.
  *
  * \author Dominik Goeddeke
  */
  class BaseMeshCell2D : public BaseMeshCell
  {

  private:
    /* *****************
    * member variables *
    *******************/
    /// edge neighbours (dimension-1 neighbours)
    std::vector<BaseMeshCell*> _edge_neighbours;

    /// point neighbours (dimension-0 neighbours)
    std::vector<BaseMeshCell*> _point_neighbours;


  public:
    /* *************************
    * constructor & destructor *
    ***************************/
    /// Default CTOR
    BaseMeshCell2D(unsigned int number)
      : BaseMeshCell(number, 2)
    {
    }

    ~BaseMeshCell2D()
    {
      // IMPORTANT: Do not delete neighbour arrays
    }

    /* *****************
    * member functions *
    *******************/
    /// returns number of neighbours of given dimension (point neighbours are dimension-0 neighbours etc.)
    inline unsigned int num_neighbours(const unsigned int dimension) const
    {
      switch (dimension)
      {
        case 0:
          return _point_neighbours.size();
        case 1:
          return _edge_neighbours.size();
        default:
          throw InternalError("2D cells do not have dimension-2 neighbours (ie, faces are not connected by faces to other faces");
      }
    }

    /// adds given cell to list of given-dimension neighbours
    inline void add_neighbour(const unsigned int dimension, BaseMeshCell * neighbour)
    {
      switch (dimension)
      {
        case 0:
          _point_neighbours.push_back(neighbour);
          break;
        case 1:
          _edge_neighbours.push_back(neighbour);
          break;
        default:
          throw InternalError("2D cells do not have dimension-2 neighbours");
      }
    }

    /**
    * \brief Returns given neighbour of given dimension
    *
    * TODO: Ich wollte das hier mit nem Iterator machen, aber ich habe es nicht hinbekommen. So muss man jetzt
    * aussenrum ne Schleife ueber alle Nachbarn der gewuenschten Dimension machen (Ende der Schleife kennt man
    * ueber num_neighbours()) damit man in diese Routine den entsprechenden Index stopfen kann.
    */
    inline BaseMeshCell * get_neighbour(const unsigned int dimension, const unsigned int index) const
    {
      switch (dimension)
      {
        case 0:
          assert(index < _point_neighbours.size());
          return _point_neighbours.at(index);
        case 1:
          assert(index < _edge_neighbours.size());
          return _edge_neighbours.at(index);
        default:
          throw InternalError("2D cells do not have dimension-2 neighbours");
      }
    }

    /// simple debugging routine: prints this cell to std::cout
    virtual void print()
    {
      BaseMeshCell::print();
      // TODO Mit Iteratoren habe ich das nicht hinbekommen
      std::cout << ". EdgeNeighbours: ";
      for (unsigned int i=0; i<_edge_neighbours.size(); ++i)
        std::cout << _edge_neighbours.at(i)->number() << " ";
      std::cout << ". PointNeighbours: ";
      for (unsigned int i=0; i<_point_neighbours.size(); ++i)
        std::cout << _point_neighbours.at(i)->number() << " ";
      std::cout << std::endl;
    }

  }; // class BaseMeshCell2D




  /**
  * \brief Class describing a base mesh
  *
  * \author Hilmar Wobker
  * \author Dominik Goeddeke
  */
  class BaseMesh
  {

  private:
    /* *****************
    * member variables *
    *******************/
    /// description of this base mesh
    std::string _description;

    /// number of base mesh cells
    unsigned int _num_cells;

    /// actual mesh, stored temporarily as a fixed-length array of pointers to mesh cells
    std::vector<BaseMeshCell*> _cells;

    /// number of partitions
    unsigned int _num_partitions;

    /// graph describing the connectivity of the base mesh (derived from _cells)
    Graph* _graph;


//    /* *****************
//    * member functions *
//    *******************/
//    /// reads (currently, basically skips) FEAST_GEOM and FEAST_FGEOM parts
//    void parse_geom(FileReaderASCII * meshfile)
//    {
//      for (unsigned int i = 0 ; i < 2 ; ++i)
//      {
//        // header
//        if (i==0)
//          meshfile->read("FEAST_GEOM");
//        else
//          meshfile->read("FEAST_FGEOM");
//        meshfile->read("2D");
//        int version_major;
//        meshfile->read(version_major);
//        int version_minor;
//        meshfile->read(version_minor);
//
//        // number of boundary components
//        unsigned int num_boundaries;
//        meshfile->read(num_boundaries);
//
//        // loop over boundary components
//        for (unsigned int ibc=0; ibc<num_boundaries; ++ibc)
//        {
//          unsigned int num_segments;
//          meshfile->read(num_segments);
//          for (unsigned int iseg = 0 ; iseg < num_segments ; ++iseg)
//          {
//            int seg_type;
//            meshfile->read(seg_type);
//            switch (seg_type)
//            {
//              case 0: // line: skip three lines
//                meshfile->skip_line();
//                meshfile->skip_line();
//                meshfile->skip_line();
//                break;
//              case 1: // circle: skip four lines
//                meshfile->skip_line();
//                meshfile->skip_line();
//                meshfile->skip_line();
//                meshfile->skip_line();
//                break;
//              default:
//                throw InternalError("Unknown boundary type found.");
//                break;
//            }
//          }
//        }
//      }
//    } // parse_geom()
//
//    /// reads MESH portion of FEAST files into this BaseMesh
//    void parse_mesh(FileReaderASCII * meshfile)
//    {
//      // header
//      meshfile->read("FEAST_MESH");
//      meshfile->read("2D");
//      int version_major;
//      meshfile->read(version_major);
//      int version_minor;
//      meshfile->read(version_minor);
//      int type;
//      meshfile->read(type);
//      // number of entities
//      unsigned int num_nodes;
//      meshfile->read(num_nodes);
//      meshfile->read(_num_cells);
//      unsigned int num_edges;
//      meshfile->read(num_edges);
//      // std::cout << num_nodes << " " << num_edges << " " << _num_cells << std::endl;
//
//      // nodes
//      for (unsigned int inode = 0 ; inode < num_nodes ; ++inode)
//      {
//        double x, y;
//        unsigned int type;
//        meshfile->read(x,y,type);
//        //std::cout << "skipping node " << x << " " << y << " " << type << std::endl;
//      }
//
//      // edges
//      for (unsigned int iedge = 0 ; iedge < num_edges ; ++iedge)
//      {
//        unsigned int n1,n2;
//        meshfile->read(n1,n2);
//        //std::cout << "skipping edge " << n1 << " " << n2 << std::endl;
//      }
//
//      // cells. Idea: First allocate everything ...
//      for (unsigned int icell = 0 ; icell < _num_cells ; ++icell)
//      {
//        BaseMeshCell * c = new BaseMeshCell2D(icell);
//        _cells.push_back(c);
//      }
//      // ... and then fill with actual data
//      for (unsigned int icell = 0 ; icell < _num_cells ; ++icell)
//      {
//        // type of this cell
//        unsigned int type;
//        meshfile->read(type);
//        switch(type)
//        {
//          case BaseMeshCell::TYPE_2D_Q1_TP:
//          {
//            _cells.at(icell)->set_type(type);
//            // skip node indices
//            unsigned int n1,n2,n3,n4;
//            meshfile->read(n1,n2,n3,n4);
//            //std::cout << "Nodes: " << n1 << " "<< n2 << " " << n3 << " " << n4 << std::endl;
//            // skip edge indices
//            unsigned int e1,e2,e3,e4;
//            meshfile->read(e1,e2,e3,e4);
//            // read edge neighbours (note the special case, only in preparation for 3D)
//            int en[4];
//            meshfile->read(en[0], en[1], en[2], en[3]);
//            //std::cout << "EN: " << en[0] << " "<< en[1] << " " << en[2] << " " << en[3] << std::endl;
//            for (unsigned int i = 0 ; i < 4 ; ++i)
//            {
//              if (en[i] > 0)
//              {
//                _cells.at(icell)->add_neighbour(1, _cells.at(en[i]-1));
//              }
//              else
//              {
//                int num_extra = abs(en[i]);
//                switch(num_extra)
//                {
//                  case 1:
//                  {
//                    int i1;
//                    meshfile->read(i1);
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i1-1));
//                    break;
//                  }
//                  case 2:
//                  {
//                    int i1,i2;
//                    meshfile->read(i1, i2);
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i1-1));
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i2-1));
//                    break;
//                  }
//                  case 3:
//                  {
//                    int i1,i2, i3;
//                    meshfile->read(i1, i2, i3);
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i1-1));
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i2-1));
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i3-1));
//                    break;
//                  }
//                  case 4:
//                  {
//                    int i1,i2,i3,i4;
//                    meshfile->read(i1, i2, i3, i4);
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i1-1));
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i2-1));
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i3-1));
//                    _cells.at(icell)->add_neighbour(1, _cells.at(i4-1));
//                    break;
//                  }
//                }
//              } // switch (num_extra)
//            } // for-loop over edge neighbours
//
//            // read point neighbours (note the special case)
//            int pn[4];
//            meshfile->read(pn[0], pn[1], pn[2], pn[3]);
//            for (unsigned int i = 0 ; i < 4 ; ++i)
//            {
//              if (pn[i] > 0)
//                _cells.at(icell)->add_neighbour(0, _cells.at(pn[i]-1));
//              else
//              {
//                int num_extra = abs(pn[i]);
//                switch(num_extra)
//                {
//                  case 1:
//                  {
//                    int i1;
//                    meshfile->read(i1);
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i1-1));
//                    break;
//                  }
//                  case 2:
//                  {
//                    int i1,i2;
//                    meshfile->read(i1, i2);
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i1-1));
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i2-1));
//                    break;
//                  }
//                  case 3:
//                  {
//                    int i1,i2, i3;
//                    meshfile->read(i1, i2, i3);
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i1-1));
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i2-1));
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i3-1));
//                    break;
//                  }
//                  case 4:
//                  {
//                    int i1,i2,i3,i4;
//                    meshfile->read(i1, i2, i3, i4);
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i1-1));
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i2-1));
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i3-1));
//                    _cells.at(icell)->add_neighbour(0, _cells.at(i4-1));
//                    break;
//                  }
//                }
//              }
//            } // loop over point neighbours
//            // finally, skip the line with refinement level etc.
//            int ref_level, ref_mode;
//            double ref_factor1, ref_factor2, ref_factor3;
//            meshfile->read(ref_level, ref_mode, ref_factor1, ref_factor2, ref_factor3);
//            // and continue with next cell
//            break;
//          }
//          default:
//            throw InternalError("Unknown cell type <" + StringUtils::stringify(type) + "> found.");
//        }
//      }
//    } // parse_mesh
//
//    /// parses partition part of FEAST files
//    void parse_partition(FileReaderASCII * meshfile)
//    {
//      // header
//      meshfile->read("FEAST_PART");
//      meshfile->read("2D");
//      int version_major;
//      meshfile->read(version_major);
//      int version_minor;
//      meshfile->read(version_minor);
//      // number of entities
//      unsigned int num_cells;
//      meshfile->read(num_cells);
//      if (num_cells != _num_cells)
//        throw InternalError("Number of BaseMeshCells from mesh file does not match that of partition file.");
//      meshfile->read(_num_partitions);
//      // TODO: THIS IS A HACK UNTIL HILMAR AND I DECIDE WHAT DTO DO WITH THE BASEMESH
//      if (_num_cells != _num_partitions)
//        throw InternalError("Sorry. Hilmar and Dominik currently only support one basemeshcell per process.");
//
//      // cells
//      for (unsigned int icell = 0 ; icell < _num_cells; ++icell)
//      {
//        unsigned int pb, mb;
//        meshfile->read(pb, mb);
//        _cells.at(icell)->set_parallel_block(pb-1);
//        _cells.at(icell)->set_matrix_block(mb-1);
//      }
//    }
//
//    /// creates connectivity graph from information stored in this BaseMesh
//    void create_graph()
//    {
//      // allocate index array
//      unsigned int* index = new unsigned int[_num_cells+1];
//
//      // graph data structure is filled by two sweeps through the cell list
//      // first sweep: count neighbours of each cell, and maintain running total to fill index array
//      // treat last index entry separately because cell array has one less entry than index array
//      unsigned int num_neighbours_so_far = 0;
//      for (unsigned int icell=0 ; icell < _num_cells ; ++icell)
//      {
//        // set neighbours counted so far to current cell
//        index[icell] = num_neighbours_so_far;
//        //std::cout << "Setting index[" << icell << "] = " << num_neighbours_so_far << std::endl;
//        // count neighbours (dimension-1 downto dimension=0 aka point-neighbours)
//        for (unsigned int dim = 0 ; dim < _cells.at(icell)->dimension() ; ++dim)
//        {
//          num_neighbours_so_far += (_cells.at(icell)->num_neighbours(dim));
//        }
//      }
//      index[_num_cells] = num_neighbours_so_far;
//      //std::cout << "Setting index[" << _num_cells << "] = " << num_neighbours_so_far << std::endl;
//
//      // second sweep through data structure
//      // second sweep adds actual neighbour cell numbers in the appropriate places into array neighbours
//      // again, treat last loop instance separately
//      unsigned int* neighbours = new unsigned int[index[_num_cells]];
//      num_neighbours_so_far = 0;
//      for (unsigned int icell=0 ; icell < _num_cells ; icell++)
//        for (unsigned int dim = 0 ; dim < _cells.at(icell)->dimension() ; ++dim)
//          for (unsigned int neigh = 0 ; neigh < _cells.at(icell)->num_neighbours(dim) ; ++neigh)
//          {
//            neighbours[num_neighbours_so_far] = _cells.at(icell)->get_neighbour(dim,neigh)->number();
//            //std::cout << "neighbours[" << num_neighbours_so_far << "] = " << neighbours[num_neighbours_so_far] << std::endl;
//            ++num_neighbours_so_far;
//          }
//
//      // now, create graph object
//      // temporarily, do not distinguish edge neighbours and diagonal neighbours
//      if (_graph != nullptr)
//      {
//        delete _graph;
//        _graph = nullptr;
//      }
//      _graph = new Graph(_num_cells, index, neighbours);
//    }


  public:
    /* *****************
    * member variables *
    *******************/

    /* *************************
    * constructor & destructor *
    ***************************/
    BaseMesh()
      : _num_cells(0), _graph(nullptr)
    {
    }

    ~BaseMesh()
    {
      // delete actual cells
      while (!_cells.empty())
      {
        // note that pop_back calls the removed element's destructor
        _cells.pop_back();
      }

      // delete graph derived from cells
      if (_graph != nullptr)
      {
        delete _graph;
        _graph = nullptr;
      }
    }


    /* ******************
    * getters & setters *
    ********************/
    inline unsigned int num_cells() const
    {
      return _num_cells;
    }

    /**
    * \brief getter for the graph
    *
    * \return Graph pointer #_graph
    */
    inline Graph* graph() const
    {
      return _graph;
    }


    /* *****************
    * member functions *
    *******************/
    /// prints this base mesh to std::cout
    void print() const
    {
      for (unsigned int i=0; i<_cells.size(); ++i)
        _cells.at(i)->print();
    }

//    /// Reads given FEAST-V3-inline file and parses its contents into this BaseMesh
//    void read_mesh(const std::string& filename)
//    {
//      FileReaderASCII * meshfile = new FileReaderASCII(filename, '#', true);
//
//      // file header
//      meshfile->read("FEAST");
//      meshfile->read("2D");
//      int version_major;
//      meshfile->read(version_major);
//      int version_minor;
//      meshfile->read(version_minor);
//      if (version_major != 3 && version_minor != 0)
//        throw InternalError("Only file format 3.0 (inline) is currently supported");
//      meshfile->read(_description);
//
//      // remaining parts of this file
//      parse_geom(meshfile);
//      parse_mesh(meshfile);
//      parse_partition(meshfile);
//      // HACK: do not read in FEAST_BC and FEAST_PROP, instead, just close the file and get on with life
//
//      // clean up file reader
//      delete meshfile;
//
//      // debug
//      print();
//
//      // and convert to Graph structure
//      create_graph();
//    }

    /**
    * \brief Hilmar's debug utility routine
    */
    void dummy_read_mesh()
    {
      // read mesh
      // ...

      // not implemented yet, so manually create a graph describing some base mesh
      // subdomain  num_neigh   neigh             num_diag    diag
      //    0          4        3,4,5,1              2        4,5
      //    1          6        0,3,4,5,6,2          3        3,4,6
      //    2          4        1,5,6,7              2        5,7
      //    3          6        8,9,4,5,1,0          3        9,5,1
      //    4          6        5,1,0,3,8,9          3        1,0,8
      //    5          6        6,2,1,0,3,4          3        2,0,3
      //    6          4        7,2,1,5              1        1
      //    7          3        12,2,6               1        2
      //    8          6        13,14,10,9,4,3       3        14,10,4
      //    9          6        4,3,8,13,14,10       3        3,13,14
      //   10          6        9,8,13,14,15,11      3        8,13,15
      //   11          4        10,14,15,12          1        14
      //   12          3        11,15,7              1        15
      //   13          4        14,10,9,8            2        10,9
      //   14          6        15,11,10,9,8,13      3        11,9,8
      //   15          4        12,11,10,14          2        12,10
      unsigned int const num_base_cells = 16;
      unsigned int* index = new unsigned int[num_base_cells+1];
      index[0]  =  0;
      index[1]  =  4;
      index[2]  = 10;
      index[3]  = 14;
      index[4]  = 20;
      index[5]  = 26;
      index[6]  = 32;
      index[7]  = 36;
      index[8]  = 39;
      index[9]  = 45;
      index[10] = 51;
      index[11] = 57;
      index[12] = 61;
      index[13] = 64;
      index[14] = 68;
      index[15] = 74;
      index[16] = 78;

      unsigned int* neigh = new unsigned int[index[num_base_cells]];
      // neighbours of subdomain 0
      neigh[0]  =  3;  neigh[1] =  4;  neigh[2] =  5;  neigh[3] =  1;
      // neighbours of subdomain 1
      neigh[4]  =  0;  neigh[5] =  3;  neigh[6] =  4;  neigh[7] =  5;  neigh[8] =  6;  neigh[9] =  2;
      // neighbours of subdomain 2
      neigh[10] =  1; neigh[11] =  5; neigh[12] =  6; neigh[13] =  7;
      // neighbours of subdomain 3
      neigh[14] =  8; neigh[15] =  9; neigh[16] =  4; neigh[17] =  5;  neigh[18] =  1; neigh[19] =  0;
      // neighbours of subdomain 4
      neigh[20] =  5; neigh[21] =  1; neigh[22] =  0; neigh[23] =  3;  neigh[24] =  8; neigh[25] =  9;
      // neighbours of subdomain 5
      neigh[26] =  6; neigh[27] =  2; neigh[28] =  1; neigh[29] =  0;  neigh[30] =  3; neigh[31] =  4;
      // neighbours of subdomain 6
      neigh[32] =  7; neigh[33] =  2; neigh[34] =  1; neigh[35] =  5;
      // neighbours of subdomain 7
      neigh[36] = 12; neigh[37] =  2; neigh[38] =  6;
      // neighbours of subdomain 8
      neigh[39] = 13; neigh[40] = 14; neigh[41] = 10; neigh[42] =  9;  neigh[43] =  4; neigh[44] =  3;
      // neighbours of subdomain 9
      neigh[45] =  4; neigh[46] =  3; neigh[47] =  8; neigh[48] = 13;  neigh[49] = 14; neigh[50] = 10;
      // neighbours of subdomain 10
      neigh[51] =  9; neigh[52] =  8; neigh[53] = 13; neigh[54] = 14;  neigh[55] = 15; neigh[56] = 11;
      // neighbours of subdomain 11
      neigh[57] = 10; neigh[58] = 14; neigh[59] = 15; neigh[60] = 12;
      // neighbours of subdomain 12
      neigh[61] = 11; neigh[62] = 15; neigh[63] =  7;
      // neighbours of subdomain 13
      neigh[64] = 14; neigh[65] = 10; neigh[66] =  9; neigh[67] =  8;
      // neighbours of subdomain 14
      neigh[68] = 15; neigh[69] = 11; neigh[70] = 10; neigh[71] =  9;  neigh[72] =  8; neigh[73] = 13;
      // neighbours of subdomain 15
      neigh[74] = 12; neigh[75] = 11; neigh[76] = 10; neigh[77] = 14;

      if (_graph != nullptr)
      {
        delete _graph;
        _graph = nullptr;
      }
      // create graph object
      // temporarily, do not distinguish edge neighbours and diagonal neighbours
      _graph = new Graph(num_base_cells, index, neigh);
      _graph->print();
    }
  }; // class BaseMesh
} // namespace FEAST

#endif // guard KERNEL_BASE_MESH_HPP
