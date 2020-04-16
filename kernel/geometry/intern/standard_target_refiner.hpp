// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STANDARD_TARGET_REFINER_HPP
#define KERNEL_GEOMETRY_INTERN_STANDARD_TARGET_REFINER_HPP 1

// includes, FEAT
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/target_set.hpp>
#include <kernel/geometry/intern/target_index_mapping.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        int cell_dim_>
      struct StandardTargetRefiner;

      template<>
      struct StandardTargetRefiner<Shape::Vertex, 0>
      {
        typedef Shape::Vertex ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in)
        {
          const Index num_verts = target_set_holder_in.get_num_entities(0);
          const TargetSetType& target_set_in = target_set_holder_in.get_target_set<0>();

          for(Index i(0); i < num_verts; ++i)
          {
            target_set_out[i] = index_offsets[0] + target_set_in[offset + i];
          }

          return num_verts;
        }
      };

      /**
        * \brief StandardTargetRefiner implementation for Simplex<...> shape: Vertex indices
        *
        * \author Constantin Christof
        */
      template<int shape_dim_>
      struct StandardTargetRefiner<Shape::Simplex<shape_dim_>, 0>
      {
        typedef Shape::Simplex<shape_dim_> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& /*index_set_holder_src*/,
          const IndexSetHolderType& /*index_set_holder_trg*/)
        {
          // fetch indices
          const TargetSet& target_set_in = target_set_holder_in.template get_target_set<shape_dim_>();
          Index num_cells = target_set_in.get_num_entities();

          // set target indices
          for(Index i(0); i < num_cells; ++i)
          {
            target_set_out[offset + i] = index_offsets[shape_dim_] + target_set_in[i];
          }

          return num_cells;
        }
      }; // StandardTargetRefiner<Shape::Simplex<...>, 0>

      /**
        * \brief StandardTargetRefiner implementation for Simplex<1> shape: Edge indices
        *
        * \author Constantin Christof
        */
      template<>
      struct StandardTargetRefiner<Shape::Simplex<1>, 1>
      {
        typedef Shape::Simplex<1> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          typedef TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

          // fetch index offset
          const Index ioe = index_offsets[1];

          // typedef index set types
          typedef IndexSet<2> IndexSetTypeEV;

          // typedef index vector type
          typedef IndexSetTypeEV::IndexTupleType IndexTupleTypeEV;

          // fetch the edge-vertex index sets
          const IndexSetTypeEV& index_src_e_v = index_set_holder_src.get_index_set<1,0>();
          const IndexSetTypeEV& index_trg_e_v = index_set_holder_trg.get_index_set<1,0>();

          // fetch target indices
          const TargetSet& target_set_v = target_set_holder_in.get_target_set<0>();
          const TargetSet& target_set_e = target_set_holder_in.get_target_set<1>();
          Index num_cells = target_set_e.get_num_entities();

          // set target indices
          for(Index i(0); i < num_cells; ++i)
          {
            // fetch edge target index
            const Index trg_e = target_set_e[i];

            // fetch target edge-vert index vector
            const IndexTupleTypeEV& trg_e_v = index_trg_e_v[trg_e];

            // fetch source edge-vert index vector
            const IndexTupleTypeEV& src_e_v = index_src_e_v[i];

            // create a target-index mapping object
            TargetIndexMappingType tim(trg_e_v, src_e_v, target_set_v);

            // fetch fine mesh target indices
            Index& e_0 = target_set_out[offset + 2*i + 0];
            Index& e_1 = target_set_out[offset + 2*i + 1];

            // calculate fine edge target indices
            e_0 = ioe + 2*trg_e + tim.map(0);
            e_1 = ioe + 2*trg_e + tim.map(1);
          }
          return 2*num_cells;
        }
      }; // StandardTargetRefiner<Shape::Simplex<1>,1>

      /**
        * \brief StandardTargetRefiner implementation for Simplex<2> shape: Vertex indices
        *
        * \author Constantin Christof
        */
      template<>
      struct StandardTargetRefiner<Shape::Simplex<2>, 0>
      {
        typedef Shape::Simplex<2> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& /*target_set_out*/,
          const Index /*offset*/,
          const Index* /*index_offsets*/,
          const TargetSetHolderType& /*target_set_holder_in*/,
          const IndexSetHolderType& /*index_set_holder_src*/,
          const IndexSetHolderType& /*index_set_holder_trg*/)
        {
          return 0;
        }
      }; // StandardTargetRefiner<Shape::Simplex<2>, 0>

      /**
        * \brief StandardTargetRefiner implementation for Simplex<2> shape: Edge indices
        *
        * \author Constantin Christof
        */
      template<>
      struct StandardTargetRefiner<Shape::Simplex<2>, 1>
      {
        typedef Shape::Simplex<2> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          typedef TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

          // fetch index offset
          const Index ioe = index_offsets[2];

          // typedef index set types
          typedef IndexSet<3> IndexSetTypeTV;

          // typedef index vector type
          typedef IndexSetTypeTV::IndexTupleType IndexTupleTypeTV;

          // fetch the tria-vertex index sets
          const IndexSetTypeTV& index_src_t_v = index_set_holder_src.get_index_set<2,0>();
          const IndexSetTypeTV& index_trg_t_v = index_set_holder_trg.get_index_set<2,0>();

          // fetch target indices
          const TargetSet& target_set_v = target_set_holder_in.get_target_set<0>();
          const TargetSet& target_set_t = target_set_holder_in.get_target_set<2>();
          Index num_cells = target_set_t.get_num_entities();

          // set target indices
          for(Index i(0); i < num_cells; ++i)
          {
            // fetch tria target index
            const Index trg_t = target_set_t[i];

            // fetch target triangle-vert index vector
            const IndexTupleTypeTV& trg_t_v = index_trg_t_v[trg_t];

            // fetch source edge-vert index vector
            const IndexTupleTypeTV& src_t_v = index_src_t_v[i];

            // create a target-index mapping object
            TargetIndexMappingType tim(trg_t_v, src_t_v, target_set_v);

            // fetch fine mesh target indices
            Index& e_0 = target_set_out[offset + 3*i + 0];
            Index& e_1 = target_set_out[offset + 3*i + 1];
            Index& e_2 = target_set_out[offset + 3*i + 2];

            // calculate fine edge target indices
            e_0 = ioe + 3*trg_t + tim.map(0);
            e_1 = ioe + 3*trg_t + tim.map(1);
            e_2 = ioe + 3*trg_t + tim.map(2);
          }
          return 3*num_cells;
        }
      }; // StandardTargetRefiner<Shape::Simplex<2>,1>

      /**
        * \brief StandardTargetRefiner implementation for Simplex<2> shape: Tri indices
        *
        * \author Constantin Christof
        */
      template<>
      struct StandardTargetRefiner<Shape::Simplex<2>, 2>
      {
        typedef Shape::Simplex<2> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          typedef TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

          // fetch index offset
          const Index iot = index_offsets[2];

          // typedef index set types
          typedef IndexSet<3> IndexSetTypeTV;

          // typedef index vector type
          typedef IndexSetTypeTV::IndexTupleType IndexTupleTypeTV;

          // fetch the tria-vertex index sets
          const IndexSetTypeTV& index_src_t_v = index_set_holder_src.get_index_set<2,0>();
          const IndexSetTypeTV& index_trg_t_v = index_set_holder_trg.get_index_set<2,0>();

          // fetch target indices
          const TargetSet& target_set_v = target_set_holder_in.get_target_set<0>();
          const TargetSet& target_set_t = target_set_holder_in.get_target_set<2>();
          Index num_cells = target_set_t.get_num_entities();

          // set target indices
          for(Index i(0); i < num_cells; ++i)
          {
            // fetch tria target index
            const Index trg_t = target_set_t[i];

            // fetch target tria-vert index vector
            const IndexTupleTypeTV& trg_t_v = index_trg_t_v[trg_t];

            // fetch source edge-vert index vector
            const IndexTupleTypeTV& src_t_v = index_src_t_v[i];

            // create a target-index mapping object
            TargetIndexMappingType tim(trg_t_v, src_t_v, target_set_v);

            // fetch fine mesh target indices
            Index& t_0 = target_set_out[offset + 4*i + 0];
            Index& t_1 = target_set_out[offset + 4*i + 1];
            Index& t_2 = target_set_out[offset + 4*i + 2];
            Index& t_3 = target_set_out[offset + 4*i + 3];

            // calculate fine triangle target indices
            t_0 = iot + 4*trg_t + tim.map(0);
            t_1 = iot + 4*trg_t + tim.map(1);
            t_2 = iot + 4*trg_t + tim.map(2);
            t_3 = iot + 4*trg_t + 3;
          }

          return 4*num_cells;
        }
      }; // StandardTargetRefiner<Shape::Simplex<2>,2>

      /**
        * \brief StandardTargetRefiner implementation for Simplex<3> shape
        *
        * \author Constantin Christof
        */
      template<int cell_dim_>
      struct StandardTargetRefiner<Shape::Simplex<3>, cell_dim_>
      {
        typedef Shape::Simplex<3> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& /*target_set_out*/,
          const Index /*offset*/,
          const Index* /*index_offsets*/,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& /*index_set_holder_src*/,
          const IndexSetHolderType& /*index_set_holder_trg*/)
        {
          // We do not have an implementation of 3D tetrahedral submeshes yet due to the lack of a
          // tetrahedral congruency sampler (see Congruency::Sampler<...> class template) which is needed
          // to implement this function. However, this function will also be called when refining a
          // 2D submesh of a 3D tetrahedral mesh -- in this case the number of tetrahedra is zero,
          // so we'll check for that to ensure that we don't prohibit this valid scenario.
          Index num_cells = target_set_holder_in.get_target_set<3>().get_num_entities();
          XASSERTM(num_cells == 0, "TargetSet refinement not implemented for Tetrahedra");

          // okay, no tetrahedron to be refined...
          return 0;
        }
      }; // StandardTargetRefiner<Shape::Simplex<3>,...>

      /**
        * \brief StandardTargetRefiner implementation for Simplex<3> shape: Vertex indices
        *
        * \note This class is only necessary to resolve the ambiguity of the previous specialisations.
        *
        * \author Constantin Christof
        */
      template<>
      struct StandardTargetRefiner<Shape::Simplex<3>,0>
      {
        typedef Shape::Simplex<3> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& /*target_set_out*/,
          const Index /*offset*/,
          const Index* /*index_offsets*/,
          const TargetSetHolderType& /*target_set_holder_in*/,
          const IndexSetHolderType& /*index_set_holder_src*/,
          const IndexSetHolderType& /*index_set_holder_trg*/)
        {
          return 0;
        }
      }; // StandardTargetRefiner<Shape::Simplex<3>,...>
      /**
        * \brief StandardTargetRefiner implementation for Hypercube<...> shape: Vertex indices
        *
        * \author Peter Zajac
        */
      template<int shape_dim_>
      struct StandardTargetRefiner<Shape::Hypercube<shape_dim_>, 0>
      {
        typedef Shape::Hypercube<shape_dim_> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& /*index_set_holder_src*/,
          const IndexSetHolderType& /*index_set_holder_trg*/)
        {
          // fetch indices
          const TargetSet& target_set_in = target_set_holder_in.template get_target_set<shape_dim_>();
          Index num_cells = target_set_in.get_num_entities();

          // set target indices
          for(Index i(0); i < num_cells; ++i)
          {
            target_set_out[offset + i] = index_offsets[shape_dim_] + target_set_in[i];
          }

          return num_cells;
        }
      }; // StandardTargetRefiner<Shape::Hypercube<...>, 0>

      /**
        * \brief StandardTargetRefiner implementation for Hypercube<1> shape: Edge indices
        *
        * \author Peter Zajac
        */
      template<>
      struct StandardTargetRefiner<Shape::Hypercube<1>, 1>
      {
        typedef Shape::Hypercube<1> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          typedef TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

          // fetch index offset
          const Index ioe = index_offsets[1];

          // typedef index set types
          typedef IndexSet<2> IndexSetTypeEV;

          // typedef index vector type
          typedef IndexSetTypeEV::IndexTupleType IndexTupleTypeEV;

          // fetch the edge-vertex index sets
          const IndexSetTypeEV& index_src_e_v = index_set_holder_src.get_index_set<1,0>();
          const IndexSetTypeEV& index_trg_e_v = index_set_holder_trg.get_index_set<1,0>();

          // fetch target indices
          const TargetSet& target_set_v = target_set_holder_in.get_target_set<0>();
          const TargetSet& target_set_e = target_set_holder_in.get_target_set<1>();
          Index num_cells = target_set_e.get_num_entities();

          // set target indices
          for(Index i(0); i < num_cells; ++i)
          {
            // fetch edge target index
            const Index trg_e = target_set_e[i];

            // fetch target edge-vert index vector
            const IndexTupleTypeEV& trg_e_v = index_trg_e_v[trg_e];

            // fetch source edge-vert index vector
            const IndexTupleTypeEV& src_e_v = index_src_e_v[i];

            // create a target-index mapping object
            TargetIndexMappingType tim(trg_e_v, src_e_v, target_set_v);

            // fetch fine mesh target indices
            Index& e_0 = target_set_out[offset + 2*i + 0];
            Index& e_1 = target_set_out[offset + 2*i + 1];

            // calculate fine edge target indices
            e_0 = ioe + 2*trg_e + tim.map(0);
            e_1 = ioe + 2*trg_e + tim.map(1);
          }

          return 2*num_cells;
        }
      }; // StandardTargetRefiner<Shape::Hypercube<1>,1>


      /**
        * \brief StandardTargetRefiner implementation for Hypercube<2> shape: Edge indices
        *
        * \author Peter Zajac
        */
      template<>
      struct StandardTargetRefiner<Shape::Hypercube<2>, 1>
      {
        typedef Shape::Hypercube<2> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          typedef TargetIndexMapping<ShapeType, 1> TargetIndexMappingType;

          // fetch index offset
          const Index ioe = index_offsets[2];

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeQV;

          // typedef index vector type
          typedef IndexSetTypeQV::IndexTupleType IndexTupleTypeQV;

          // fetch the quad-vertex index sets
          const IndexSetTypeQV& index_src_q_v = index_set_holder_src.get_index_set<2,0>();
          const IndexSetTypeQV& index_trg_q_v = index_set_holder_trg.get_index_set<2,0>();

          // fetch target indices
          const TargetSet& target_set_v = target_set_holder_in.get_target_set<0>();
          const TargetSet& target_set_q = target_set_holder_in.get_target_set<2>();
          Index num_cells = target_set_q.get_num_entities();

          // set target indices
          for(Index i(0); i < num_cells; ++i)
          {
            // fetch quad target index
            const Index trg_q = target_set_q[i];

            // fetch target quad-vert index vector
            const IndexTupleTypeQV& trg_q_v = index_trg_q_v[trg_q];

            // fetch source edge-vert index vector
            const IndexTupleTypeQV& src_q_v = index_src_q_v[i];

            // create a target-index mapping object
            TargetIndexMappingType tim(trg_q_v, src_q_v, target_set_v);

            // fetch fine mesh target indices
            Index& e_0 = target_set_out[offset + 4*i + 0];
            Index& e_1 = target_set_out[offset + 4*i + 1];
            Index& e_2 = target_set_out[offset + 4*i + 2];
            Index& e_3 = target_set_out[offset + 4*i + 3];

            // calculate fine edge target indices
            e_0 = ioe + 4*trg_q + tim.map(0);
            e_1 = ioe + 4*trg_q + tim.map(1);
            e_2 = ioe + 4*trg_q + tim.map(2);
            e_3 = ioe + 4*trg_q + tim.map(3);
          }

          return 4*num_cells;
        }
      }; // StandardTargetRefiner<Shape::Hypercube<2>,1>

      /**
        * \brief StandardTargetRefiner implementation for Hypercube<2> shape: Quad indices
        *
        * \author Peter Zajac
        */
      template<>
      struct StandardTargetRefiner<Shape::Hypercube<2>, 2>
      {
        typedef Shape::Hypercube<2> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index* index_offsets,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          typedef TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

          // fetch index offset
          const Index ioe = index_offsets[2];

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeQV;

          // typedef index vector type
          typedef IndexSetTypeQV::IndexTupleType IndexTupleTypeQV;

          // fetch the quad-vertex index sets
          const IndexSetTypeQV& index_src_q_v = index_set_holder_src.get_index_set<2,0>();
          const IndexSetTypeQV& index_trg_q_v = index_set_holder_trg.get_index_set<2,0>();

          // fetch target indices
          const TargetSet& target_set_v = target_set_holder_in.get_target_set<0>();
          const TargetSet& target_set_q = target_set_holder_in.get_target_set<2>();
          Index num_cells = target_set_q.get_num_entities();

          // set target indices
          for(Index i(0); i < num_cells; ++i)
          {
            // fetch quad target index
            const Index trg_q = target_set_q[i];

            // fetch target quad-vert index vector
            const IndexTupleTypeQV& trg_q_v = index_trg_q_v[trg_q];

            // fetch source edge-vert index vector
            const IndexTupleTypeQV& src_q_v = index_src_q_v[i];

            // create a target-index mapping object
            TargetIndexMappingType tim(trg_q_v, src_q_v, target_set_v);

            // fetch fine mesh target indices
            Index& q_0 = target_set_out[offset + 4*i + 0];
            Index& q_1 = target_set_out[offset + 4*i + 1];
            Index& q_2 = target_set_out[offset + 4*i + 2];
            Index& q_3 = target_set_out[offset + 4*i + 3];

            // calculate fine quad target indices
            q_0 = ioe + 4*trg_q + tim.map(0);
            q_1 = ioe + 4*trg_q + tim.map(1);
            q_2 = ioe + 4*trg_q + tim.map(2);
            q_3 = ioe + 4*trg_q + tim.map(3);
          }

          return 4*num_cells;
        }
      }; // StandardTargetRefiner<Shape::Hypercube<2>,2>

      /**
        * \brief StandardTargetRefiner implementation for Hypercube<3> shape
        *
        * \author Peter Zajac
        */
      template<int cell_dim_>
      struct StandardTargetRefiner<Shape::Hypercube<3>, cell_dim_>
      {
        typedef Shape::Hypercube<3> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& /*target_set_out*/,
          const Index /*offset*/,
          const Index* /*index_offsets*/,
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& /*index_set_holder_src*/,
          const IndexSetHolderType& /*index_set_holder_trg*/)
        {
          // We do not have an implementation of 3D hexahedral submeshes yet due to the lack of a
          // hexahedral congruency sampler (see Congruency::Sampler<...> class template) which is needed
          // to implement this function. However, this function will also be called when refining a
          // 2D submesh of a 3D hexahedral mesh -- in this case the number of hexahedra is zero,
          // so we'll check for that to ensure that we don't prohibit this valid scenario.
          Index num_cells = target_set_holder_in.get_target_set<3>().get_num_entities();
          XASSERTM(num_cells == 0, "TargetSet refinement not implemented for Hexahedra");

          // okay, no hexahedra to be refined...
          return 0;
        }
      }; // StandardTargetRefiner<Shape::Hypercube<3>,...>

      /**
        * \brief StandardTargetRefiner implementation for Hypercube<3> shape: Vertex indices
        *
        * \note This class is only necessary to resolve the ambiguity of the previous specialisations.
        *
        * \author Peter Zajac
        */
      template<>
      struct StandardTargetRefiner<Shape::Hypercube<3>,0>
      {
        typedef Shape::Hypercube<3> ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          TargetSetType& /*target_set_out*/,
          const Index /*offset*/,
          const Index* /*index_offsets*/,
          const TargetSetHolderType& /*target_set_holder_in*/,
          const IndexSetHolderType& /*index_set_holder_src*/,
          const IndexSetHolderType& /*index_set_holder_trg*/)
        {
          return 0;
        }
      }; // StandardTargetRefiner<Shape::Hypercube<3>,...>

      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */

      template<
        typename Shape_,
        int cell_dim_,
        int shape_dim_ = Shape_::dimension>
      class TargetRefineShapeWrapper
      {
      public:
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "TargetRefineShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_)
            + "," + stringify(shape_dim_) + ">";
        }

        static Index refine(
          TargetSetType& target_set_out,
          const Index index_offsets[],
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          // recursive call of TargetRefineShapeWrapper
          typedef typename Shape::FaceTraits<ShapeType, shape_dim_ - 1>::ShapeType FacetType;
          Index offset = TargetRefineShapeWrapper<FacetType, cell_dim_>::refine(
            target_set_out,
            index_offsets,
            target_set_holder_in,
            index_set_holder_src,
            index_set_holder_trg);

          // call target refiner
          StandardTargetRefiner<ShapeType, cell_dim_>::refine(
            target_set_out,
            offset,
            index_offsets,
            target_set_holder_in,
            index_set_holder_src,
            index_set_holder_trg);

          // return new offset
          return offset + Intern::StandardRefinementTraits<ShapeType, cell_dim_>::count *
            target_set_holder_in.get_num_entities(shape_dim_);
        }
      };

      template<
        typename Shape_,
        int cell_dim_>
      class TargetRefineShapeWrapper<Shape_, cell_dim_, cell_dim_>
      {
      public:
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef TargetSetHolder<CellType> TargetSetHolderType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "TargetRefineShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_)
            + "," + stringify(cell_dim_) + ">";
        }

        static Index refine(
          TargetSetType& target_set_out,
          const Index index_offsets[],
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          // call target refiner
          StandardTargetRefiner<ShapeType, cell_dim_>::refine(
            target_set_out,
            0,
            index_offsets,
            target_set_holder_in,
            index_set_holder_src,
            index_set_holder_trg);

          // return new offset
          return Intern::StandardRefinementTraits<ShapeType, cell_dim_>::count
            * target_set_holder_in.get_num_entities(cell_dim_);
        }
      };

      template<int cell_dim_>
      class TargetRefineShapeWrapper<Shape::Vertex, cell_dim_, cell_dim_>
      {
      public:
        typedef Shape::Vertex ShapeType;
        typedef TargetSet TargetSetType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef TargetSetHolder<CellType> TargetSetHolderType;
        typedef IndexSetHolder<Shape::Vertex> IndexSetHolderType;

        static String name()
        {
          return "TargetRefineShapeWrapper<Vertex," + stringify(cell_dim_) + "," + stringify(cell_dim_) + ">";
        }

        static Index refine(
          TargetSetType& target_set_out,
          const Index index_offsets[],
          const TargetSetHolderType& target_set_holder_in,
          const IndexSetHolderType& /*index_set_holder_src*/,
          const IndexSetHolderType& /*index_set_holder_trg*/)
        {
          // call target refiner
          StandardTargetRefiner<ShapeType, cell_dim_>
            ::refine(target_set_out, 0, index_offsets, target_set_holder_in);

          // return new offset
          return Intern::StandardRefinementTraits<ShapeType, cell_dim_>::count
            * target_set_holder_in.get_num_entities(cell_dim_);
        }
      };

      /* ********************************************************************************************************* */

      template<
        typename Shape_,
        int cell_dim_ = Shape_::dimension>
      class TargetRefineWrapper
      {
      public:
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef TargetSetHolder<CellType> TargetSetHolderTypeOut;
        typedef TargetSetHolder<ShapeType> TargetSetHolderTypeIn;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "TargetRefineWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ">";
        }

        static void refine(
          TargetSetHolderTypeOut& target_set_holder_out,
          const Index num_entities_trg[],
          const TargetSetHolderTypeIn& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          // recursive call of TargetRefineWrapper
          TargetRefineWrapper<ShapeType, cell_dim_ - 1>::refine(
            target_set_holder_out,
            num_entities_trg,
            target_set_holder_in,
            index_set_holder_src,
            index_set_holder_trg);

          // calculate index offsets
          Index index_offsets[Shape_::dimension+1];
          Intern::EntityCounter<Intern::StandardRefinementTraits, ShapeType, cell_dim_>
            ::offset(index_offsets, num_entities_trg);

          // call shape wrapper
          TargetRefineShapeWrapper<ShapeType, cell_dim_>::refine(
            target_set_holder_out.template get_target_set<cell_dim_>(),
            index_offsets,
            target_set_holder_in,
            index_set_holder_src,
            index_set_holder_trg);
        }
      };

      template<typename Shape_>
      class TargetRefineWrapper<Shape_, 0>
      {
      public:
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<Shape::Vertex> TargetSetHolderTypeOut;
        typedef TargetSetHolder<ShapeType> TargetSetHolderTypeIn;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "TargetRefineWrapper<" + Shape_::name() + ",0>";
        }

        static void refine(
          TargetSetHolderTypeOut& target_set_holder_out,
          const Index num_entities_trg[],
          const TargetSetHolderTypeIn& target_set_holder_in,
          const IndexSetHolderType& index_set_holder_src,
          const IndexSetHolderType& index_set_holder_trg)
        {
          // calculate index offsets
          Index index_offsets[Shape_::dimension+1];
          Intern::EntityCounter<Intern::StandardRefinementTraits, ShapeType, 0>::offset(index_offsets, num_entities_trg);

          // call shape wrapper
          TargetRefineShapeWrapper<ShapeType, 0>::refine(
            target_set_holder_out.template get_target_set<0>(),
            index_offsets,
            target_set_holder_in,
            index_set_holder_src,
            index_set_holder_trg);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT


#endif // KERNEL_GEOMETRY_INTERN_STANDARD_TARGET_REFINER_HPP
