#pragma once
#ifndef KERNEL_GEOMETRY_STANDARD_VERTEX_REFINER_HPP
#define KERNEL_GEOMETRY_STANDARD_VERTEX_REFINER_HPP 1

// includes, FEAT
#include <kernel/geometry/intern/vertex_abacus.hpp>
#include <kernel/geometry/intern/standard_refinement_traits.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        typename VertexSet_>
      struct StandardVertexRefiner;

      /**
       * \brief Vertex refiner implementation for Vertex shape
       *
       * \author Peter Zajac
       */
      template<typename VertexSet_>
      struct StandardVertexRefiner<Shape::Vertex, VertexSet_>
      {
        typedef VertexSet_ VertexSetType;

        static Index refine(
          Index offset,
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in)
        {
          // get the number of coarse mesh vertices
          Index num_verts = vertex_set_in.get_num_vertices();
          ASSERT_(vertex_set_out.get_num_vertices() >= num_verts);

          // create a vertex-abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // loop over all vertices
          for(Index i(0); i < num_verts; ++i)
          {
            // copy source vertex
            abacus.copy(vertex_set_out[offset + i], vertex_set_in[i]);
          }

          // return number of created vertices
          return num_verts;
        }
      }; // class StandardVertexRefiner<Vertex,...>

      /**
       * \brief Vertex refiner implementation for Hypercube<...> shape
       *
       * \author Peter Zajac
       */
      template<
        int cell_dim_,
        typename VertexSet_>
      struct StandardVertexRefiner<Shape::Hypercube<cell_dim_>, VertexSet_>
      {
        typedef Shape::Hypercube<cell_dim_> ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef VertexSet_ VertexSetType;

        static Index refine(
          Index offset,
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const IndexSetType& index_set_in)
        {
          typedef typename IndexSetType::ConstIndexVectorReference ConstIndexVectorReference;
          typedef typename VertexSetType::VertexReference VertexReference;
          typedef typename VertexSetType::CoordType CoordType;

          // scaling factor
          static const CoordType scale = CoordType(1) / CoordType(IndexSetType::num_indices);

          // get number of cells
          Index num_cells = index_set_in.get_num_entities();

          // create a vertex abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // loop over all cells
          for(Index i(0); i < num_cells; ++i)
          {
            // get input index vector
            ConstIndexVectorReference idx_in = index_set_in[i];

            // get output vertex
            VertexReference vtx_out = vertex_set_out[offset + i];

            // clear output vertex
            abacus.clear(vtx_out);

            // add all other vertices onto it
            for(int k(0); k < IndexSetType::num_indices; ++k)
            {
              abacus.add(vtx_out, vertex_set_in[idx_in[k]]);
            }

            // scale the output vertex
            abacus.scale(vtx_out, scale);
          }

          // return number of created vertices
          return num_cells;
        }
      }; // struct StandardVertexRefiner<Hypercube<...>,...>

      /**
       * \brief Vertex refiner implementation for Simplex<1> shape
       *
       * \author Constantin Christof
       */
      template<
        typename VertexSet_>
      struct StandardVertexRefiner<Shape::Simplex<1>, VertexSet_>
      {
        typedef Shape::Simplex<1> ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef VertexSet_ VertexSetType;

        static Index refine(
          Index offset,
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const IndexSetType& index_set_in)
        {
          typedef typename IndexSetType::ConstIndexVectorReference ConstIndexVectorReference;
          typedef typename VertexSetType::VertexReference VertexReference;
          typedef typename VertexSetType::CoordType CoordType;

          // scaling factor
          static const CoordType scale = CoordType(1) / CoordType(IndexSetType::num_indices);

          // get number of cells
          Index num_cells = index_set_in.get_num_entities();

          // create a vertex abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // loop over all cells
          for(Index i(0); i < num_cells; ++i)
          {
            // get input index vector
            ConstIndexVectorReference idx_in = index_set_in[i];

            // get output vertex
            VertexReference vtx_out = vertex_set_out[offset + i];

            // clear output vertex
            abacus.clear(vtx_out);

            // add all other vertices onto it
            for(int k(0); k < IndexSetType::num_indices; ++k)
            {
              abacus.add(vtx_out, vertex_set_in[idx_in[k]]);
            }

            // scale the output vertex
            abacus.scale(vtx_out, scale);
          }

          // return number of created vertices
          return num_cells;
        }
      }; // struct VertexRefiner<Simplex<...>,...>

      /**
       * \brief Vertex refiner implementation for Simplex<2> shape
       *
       * \author Constantin Christof
       */
      template<
        typename VertexSet_>
      struct StandardVertexRefiner<Shape::Simplex<2>, VertexSet_>
      {
        typedef Shape::Simplex<2> ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef VertexSet_ VertexSetType;

        static Index refine(
          Index /*offset*/,
          VertexSetType& /*vertex_set_out*/,
          const VertexSetType& /*vertex_set_in*/,
          const IndexSetType& /*index_set_in*/)
        {
          // return number of created vertices
          return 0;
        }
      }; // struct StandardVertexRefiner<Simplex<...>,...>

      /**
       * \brief Vertex refiner implementation for Simplex<3> shape
       *
       * \author Constantin Christof
       */
      template<
        typename VertexSet_>
      struct StandardVertexRefiner<Shape::Simplex<3>, VertexSet_>
      {
        typedef Shape::Simplex<3> ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef VertexSet_ VertexSetType;

        static Index refine(
          Index offset,
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const IndexSetType& index_set_in)
        {
          typedef typename IndexSetType::ConstIndexVectorReference ConstIndexVectorReference;
          typedef typename VertexSetType::VertexReference VertexReference;
          typedef typename VertexSetType::CoordType CoordType;

          // scaling factor
          static const CoordType scale = CoordType(1) / CoordType(IndexSetType::num_indices);

          // get number of cells
          Index num_cells = index_set_in.get_num_entities();

          // create a vertex abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // loop over all cells
          for(Index i(0); i < num_cells; ++i)
          {
            // get input index vector
            ConstIndexVectorReference idx_in = index_set_in[i];

            // get output vertex
            VertexReference vtx_out = vertex_set_out[offset + i];

            // clear output vertex
            abacus.clear(vtx_out);

            // add all other vertices onto it
            for(int k(0); k < IndexSetType::num_indices; ++k)
            {
              abacus.add(vtx_out, vertex_set_in[idx_in[k]]);
            }

            // scale the output vertex
            abacus.scale(vtx_out, scale);
          }

          // return number of created vertices
          return num_cells;
        }
      }; // struct StandardVertexRefiner<Simplex<3>,...>

      /**
       * \brief Standard Vertex Refinement Wrapper class template
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        typename VertexSet_>
      struct StandardVertexRefineWrapper
      {
        typedef VertexSet_ VertexSetType;
        typedef IndexSetHolder<Shape_> IndexSetHolderType;
        typedef IndexSet<Shape::FaceTraits<Shape_, 0>::count> IndexSetType;

        static Index refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const IndexSetHolderType& index_set_holder_in)
        {
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension-1>::ShapeType FacetType;

          // recursive call of VertexRefineWrapper
          Index offset = StandardVertexRefineWrapper<FacetType, VertexSet_>
            ::refine(vertex_set_out, vertex_set_in, index_set_holder_in);

          // get index set of current shape
          const IndexSetType& index_set_in = index_set_holder_in.
            template get_index_set_wrapper<Shape_::dimension>().template get_index_set<0>();

          // call VertexRefiner
          Index num_verts = StandardVertexRefiner<Shape_, VertexSet_>
            ::refine(offset, vertex_set_out, vertex_set_in, index_set_in);

          // validate number of created vertices
          ASSERT(num_verts == (StandardRefinementTraits<Shape_, 0>::count * index_set_in.get_num_entities()),
            "VertexRefiner output does not match StdRefTraits prediction");

          // return new offset
          return offset + num_verts;
        }
      }; // struct StandardVertexRefineWrapper<...>

      template<typename VertexSet_>
      struct StandardVertexRefineWrapper<Shape::Vertex, VertexSet_>
      {
        typedef VertexSet_ VertexSetType;
        typedef IndexSetHolder<Shape::Vertex> IndexSetHolderType;

        static Index refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const IndexSetHolderType& /*index_set_holder_in*/)
        {
          // call VertexRefiner
          Index num_verts =  StandardVertexRefiner<Shape::Vertex,VertexSet_>
            ::refine(0, vertex_set_out, vertex_set_in);

          // validate number of created vertices
          ASSERT(num_verts == (StandardRefinementTraits<Shape::Vertex, 0>::count * vertex_set_in.get_num_vertices()),
            "VertexRefiner output does not match StdRefTraits prediction");

          // return new offset
          return num_verts;
        }
      }; // struct StandardVertexRefineWrapper<Vertex,...>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_STANDARD_VERTEX_REFINER_HPP
