#pragma once
#ifndef KERNEL_BASE_MESH_VERTEX_HPP
#define KERNEL_BASE_MESH_VERTEX_HPP 1

// includes, system
#include <iostream> // for std::ostream
//#include <iomanip>  // for std::setw

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/base_mesh/item.hpp>

/// FEAST namespace
namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
  namespace BaseMesh
  {
    /**
    * \brief Container for a base mesh vertex, templated with dimension of the entire base mesh.
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    */
    template<unsigned char world_dim_>
    class Vertex
      : public Item
    {
    private:

      /// coordinates of this vertex
      double _coords[world_dim_];


    public:

      /// default CTOR
      Vertex()
        : Item()
      {
        CONTEXT("BaseMesh::Vertex::Vertex()");
        for(int i(0) ; i < world_dim_ ; ++i)
        {
          _coords[i] = 0.0;
        }
      }


      /// CTOR
      Vertex(double coords[])
        : Item()
      {
        CONTEXT("BaseMesh::Vertex::Vertex()");
        for(int i(0) ; i < world_dim_ ; ++i)
        {
          _coords[i] = coords[i];
        }
      }


      /// DTOR
      ~Vertex()
      {
        CONTEXT("BaseMesh::Vertex::~Vertex()");
      }


      /// returns a specific coordinate
      inline double coord(unsigned char const index) const
      {
        CONTEXT("BaseMesh::Vertex::coord()");
        ASSERT(index < world_dim_, "");
        return _coords[index];
      }


      /// returns a pointer to the coordinate array
      inline double const* coords() const
      {
        CONTEXT("BaseMesh::Vertex::coords()");
        return _coords;
      }


      /// sets the coordinates of this vertex
      inline void set_coords(double const* coords)
      {
        CONTEXT("BaseMesh::Vertex::set_coords()");
        for (unsigned char index(0) ; index < world_dim_ ; ++index)
        {
          _coords[index] = coords[index];
        }
      }

      /// sets a  a specific coordinate
      inline void set_coord(unsigned char const index, double const coord)
      {
        CONTEXT("BaseMesh::Vertex::set_coord()");
        ASSERT(index < world_dim_, "");
        _coords[index] = coord;
      }


      /// prints this vertex to the given stream
      inline void print(std::ostream& stream)
      {
        CONTEXT("BaseMesh::Vertex::print()");
        stream << "V";
        print_index(stream);
//        stream.precision(3);
//        stream.setf(std::ios::fixed);
//        stream << ": (" << std::setw(7) << _coords[0];
        stream << ": (" << _coords[0];
        for(unsigned char i(1) ; i < world_dim_ ; ++i)
        {
//          stream  << ", " << std::setw(7) << _coords[i];
          stream  << ", " << _coords[i];
        }
        stream << ")";
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BASE_MESH_VERTEX_HPP
