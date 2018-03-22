#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STANDARD_ATTRIB_REFINER_HPP
#define KERNEL_GEOMETRY_INTERN_STANDARD_ATTRIB_REFINER_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>
#include <kernel/geometry/attribute_set.hpp>
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
        typename AttribSet_>
      struct StandardAttribRefiner;

      /**
       * \brief Attribute refiner implementation for Vertex shape
       *
       * \author Peter Zajac
       */
      template<typename AttribSet_>
      struct StandardAttribRefiner<Shape::Vertex, AttribSet_>
      {
        typedef AttribSet_ AttribSetType;

        static Index refine(
          const Index offset,
          AttribSetType& attrib_set_out,
          const AttribSetType& attrib_set_in)
        {
          XASSERT(attrib_set_out.get_dimension() == attrib_set_in.get_dimension());

          // get the number of coarse mesh attributes
          const Index num_values = attrib_set_in.get_num_values();
          const int dim_attrib = attrib_set_in.get_dimension();

          XASSERT(attrib_set_out.get_num_values() >= offset + num_values);

          // loop over all attributes
          for(Index i(0); i < num_values; ++i)
          {
            for(int j(0); j < dim_attrib; ++j)
            {
              attrib_set_out(offset+i, j) = attrib_set_in(i, j);
            }
          }

          // return number of created vertices
          return num_values;
        }
      }; // class StandardAttribRefiner<Vertex,...>

      /**
       * \brief Attribute refiner implementation for Hypercube<...> shape
       *
       * \author Peter Zajac
       */
      template<int cell_dim_, typename AttribSet_>
      struct StandardAttribRefiner<Shape::Hypercube<cell_dim_>, AttribSet_>
      {
        typedef Shape::Hypercube<cell_dim_> ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef AttribSet_ AttribSetType;

        static Index refine(
          Index offset,
          AttribSetType& attrib_set_out,
          const AttribSetType& attrib_set_in,
          const IndexSetType& index_set_in)
        {
          // dimensional assert
          XASSERT(attrib_set_out.get_dimension() == attrib_set_in.get_dimension());

          typedef typename AttribSetType::DataType DataType;

          // scaling factor
          static const DataType scale = DataType(1) / DataType(IndexSetType::num_indices);

          // get number of cells
          const Index num_cells = index_set_in.get_num_entities();
          const int dim_attrib = attrib_set_in.get_dimension();

          XASSERT(attrib_set_out.get_num_values() >= offset + num_cells);

          // loop over all cells
          for(Index i(0); i < num_cells; ++i)
          {
            // get input index tuple
            const auto& idx_in = index_set_in[i];

            // clear output attributes
            for(int j(0); j < dim_attrib; ++j)
              attrib_set_out(offset+i, j) = DataType(0);

            // sum up the input attributes
            for(int k(0); k < IndexSetType::num_indices; ++k)
            {
              for(int j(0); j < dim_attrib; ++j)
                attrib_set_out(offset+i, j) += attrib_set_in(idx_in[k], j);
            }

            // scale attribute
            for(int j(0); j < dim_attrib; ++j)
              attrib_set_out(offset+i, j) *= scale;
          }

          // return number of created vertices
          return num_cells;
        }
      }; // struct StandardAttribRefiner<Hypercube<...>,...>

      /**
       * \brief Attrib refiner implementation for Simplex<...> shape
       *
       * \author Peter Zajac
       */
      template<int cell_dim_, typename AttribSet_>
      struct StandardAttribRefiner<Shape::Simplex<cell_dim_>, AttribSet_>
      {
        typedef Shape::Simplex<cell_dim_> ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef AttribSet_ AttribSetType;

        static Index refine(
          Index offset,
          AttribSetType& attrib_set_out,
          const AttribSetType& attrib_set_in,
          const IndexSetType& index_set_in)
        {
          // dimensional assert
          XASSERT(attrib_set_out.get_dimension() == attrib_set_in.get_dimension());

          typedef typename AttribSetType::DataType DataType;

          // scaling factor
          static const DataType scale = DataType(1) / DataType(IndexSetType::num_indices);

          // get number of cells
          const Index num_cells = index_set_in.get_num_entities();
          const int dim_attrib = attrib_set_in.get_dimension();

          XASSERT(attrib_set_out.get_num_values() >= offset + num_cells);

          // loop over all cells
          for(Index i(0); i < num_cells; ++i)
          {
            // get input index tuple
            const auto& idx_in = index_set_in[i];

            // clear output attributes
            for(int j(0); j < dim_attrib; ++j)
              attrib_set_out(offset+i, j) = DataType(0);

            // sum up the input attributes
            for(int k(0); k < IndexSetType::num_indices; ++k)
            {
              for(int j(0); j < dim_attrib; ++j)
                attrib_set_out(offset+i, j) += attrib_set_in(idx_in[k], j);
            }

            // scale attribute
            for(int j(0); j < dim_attrib; ++j)
              attrib_set_out(offset+i, j) *= scale;
          }

          // return number of created vertices
          return num_cells;
        }
      }; // struct AttribRefiner<Simplex<...>,...>

      /**
       * \brief Attrib refiner implementation for Simplex<2> shape
       *
       * \author Peter Zajac
       */
      template<typename AttribSet_>
      struct StandardAttribRefiner<Shape::Simplex<2>, AttribSet_>
      {
        typedef Shape::Simplex<2> ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef AttribSet_ AttribSetType;

        static Index refine(
          Index /*offset*/,
          AttribSetType& /*attrib_set_out*/,
          const AttribSetType& /*attrib_set_in*/,
          const IndexSetType& /*index_set_in*/)
        {
          // return number of created attributes
          return 0;
        }
      }; // struct StandardAttribRefiner<Simplex<2>,...>

      /**
       * \brief Standard Vertex Refinement Wrapper class template
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        typename AttribSet_>
      struct StandardAttribRefineWrapper
      {
        typedef AttribSet_ AttribSetType;
        typedef IndexSetHolder<Shape_> IndexSetHolderType;
        typedef IndexSet<Shape::FaceTraits<Shape_, 0>::count> IndexSetType;

        static Index refine(
          AttribSetType& attrib_set_out,
          const AttribSetType& attrib_set_in,
          const IndexSetHolderType& index_set_holder_in)
        {
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension-1>::ShapeType FacetType;

          // recursive call of AttribRefineWrapper
          Index offset = StandardAttribRefineWrapper<FacetType, AttribSet_>
            ::refine(attrib_set_out, attrib_set_in, index_set_holder_in);

          // get index set of current shape
          const IndexSetType& index_set_in = index_set_holder_in.
            template get_index_set_wrapper<Shape_::dimension>().template get_index_set<0>();

          // call AttribRefiner
          Index num_values = StandardAttribRefiner<Shape_, AttribSet_>
            ::refine(offset, attrib_set_out, attrib_set_in, index_set_in);

          // validate number of created vertices
          XASSERTM(num_values == (StandardRefinementTraits<Shape_, 0>::count * index_set_in.get_num_entities()),
            "AttribRefiner output does not match StdRefTraits prediction");

          // return new offset
          return offset + num_values;
        }
      }; // struct StandardAttribRefineWrapper<...>

      template<typename AttribSet_>
      struct StandardAttribRefineWrapper<Shape::Vertex, AttribSet_>
      {
        typedef AttribSet_ AttribSetType;
        typedef IndexSetHolder<Shape::Vertex> IndexSetHolderType;

        static Index refine(
          AttribSetType& attrib_set_out,
          const AttribSetType& attrib_set_in,
          const IndexSetHolderType& /*index_set_holder_in*/)
        {
          // call AttribRefiner
          Index num_values =  StandardAttribRefiner<Shape::Vertex,AttribSet_>
            ::refine(0, attrib_set_out, attrib_set_in);

          // validate number of created values
          XASSERTM(num_values == (StandardRefinementTraits<Shape::Vertex, 0>::count * attrib_set_in.get_num_values()),
            "AttribRefiner output does not match StdRefTraits prediction");

          // return new offset
          return num_values;
        }
      }; // struct StandardAttribRefineWrapper<Vertex,...>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_STANDARD_ATTRIB_REFINER_HPP
