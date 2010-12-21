#pragma once
#ifndef KERNEL_BASE_MESH_VERTEX_HPP
#define KERNEL_BASE_MESH_VERTEX_HPP 1

// includes, system
#include <iostream> // for std::ostream
//#include <iomanip>  // for std::setw
#include <cassert>  // for assert()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_item.hpp>

namespace FEAST
{
  namespace BaseMesh
  {
    /**
    * \brief Container for a base mesh vertex, templated with dimension of the entire base mesh.
    *
    * \author Hilmar Wobker
    * \author Dominik Goeddeke
    * \author Peter Zajac
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
        for(int i(0) ; i < world_dim_ ; ++i)
        {
          _coords[i] = 0.0;
        }
      }


      /// CTOR
      Vertex(double coords[])
        : Item()
      {
        for(int i(0) ; i < world_dim_ ; ++i)
        {
          _coords[i] = coords[i];
        }
      }


      /// DTOR
      ~Vertex()
      {
      }


      /// returns a specific coordinate
      inline double coord(unsigned char const index) const
      {
        assert(index < world_dim_);
        return _coords[index];
      }


      /// returns a pointer to the coordinate array
      inline double const* coords() const
      {
        return _coords;
      }


      /// sets the coordinates of this vertex
      inline void set_coords(double const* coords)
      {
        for (unsigned char index(0) ; index < world_dim_ ; ++index)
        {
          _coords[index] = coords[index];
        }
      }

      /// sets a  a specific coordinate
      inline void set_coord(unsigned char const index, double const coord)
      {
        assert(index < world_dim_);
        _coords[index] = coord;
      }


      /// prints this vertex to the given stream
      virtual void print(std::ostream& stream)
      {
        stream << "V";
        Item::print(stream);
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
