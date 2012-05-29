#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_HYPERCUBE_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_HYPERCUBE_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal/standard_refinement/target_refiner.hpp>
#include <kernel/geometry/congruency/target_index_mapping.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Conformal
    {
      /// \cond internal
      namespace StandardRefinement
      {
        template<int shape_dim_>
        struct TargetRefiner<Shape::Hypercube<shape_dim_>, 0>
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
        }; // TargetRefiner<Shape::Hypercube<...>, 0>

        template<>
        struct TargetRefiner<Shape::Hypercube<1>, 1>
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
            typedef Congruency::TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

            // fetch index offset
            const Index ioe = index_offsets[1];

            // typedef index set types
            typedef IndexSet<2> IndexSetTypeEV;

            // typedef index vector type
            typedef IndexSetTypeEV::IndexVectorType IndexVectorTypeEV;

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
              const IndexVectorTypeEV& trg_e_v = index_trg_e_v[trg_e];

              // fetch source edge-vert index vector
              const IndexVectorTypeEV& src_e_v = index_src_e_v[i];

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
        }; // TargetRefiner<Shape::Hypercube<1>,1>


        template<>
        struct TargetRefiner<Shape::Hypercube<2>, 1>
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
            typedef Congruency::TargetIndexMapping<ShapeType, 1> TargetIndexMappingType;

            // fetch index offset
            const Index ioe = index_offsets[2];

            // typedef index set types
            typedef IndexSet<4> IndexSetTypeQV;

            // typedef index vector type
            typedef IndexSetTypeQV::IndexVectorType IndexVectorTypeQV;

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
              const IndexVectorTypeQV& trg_q_v = index_trg_q_v[trg_q];

              // fetch source edge-vert index vector
              const IndexVectorTypeQV& src_q_v = index_src_q_v[i];

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
        }; // TargetRefiner<Shape::Hypercube<2>,1>

        template<>
        struct TargetRefiner<Shape::Hypercube<2>, 2>
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
            typedef Congruency::TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

            // fetch index offset
            const Index ioe = index_offsets[2];

            // typedef index set types
            typedef IndexSet<4> IndexSetTypeQV;

            // typedef index vector type
            typedef IndexSetTypeQV::IndexVectorType IndexVectorTypeQV;

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
              const IndexVectorTypeQV& trg_q_v = index_trg_q_v[trg_q];

              // fetch source edge-vert index vector
              const IndexVectorTypeQV& src_q_v = index_src_q_v[i];

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
        }; // TargetRefiner<Shape::Hypercube<2>,2>

        template<int cell_dim_>
        struct TargetRefiner<Shape::Hypercube<3>, cell_dim_>
        {
          typedef Shape::Hypercube<3> ShapeType;
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
            // We do not have an implementation of 3D hexahedral submeshes yet due to the lack of a
            // hexahedral congruency sampler (see Congruency::Sampler<...> class template) which is needed
            // to implement this function. However, this function will also be called when refining a
            // 2D submesh of a 3D hexahedral mesh -- in this case the number of hexahedra is zero,
            // so we'll check for that to ensure that we don't prohibit this valid scenario.
            Index num_cells = target_set_holder_in.get_target_set<3>().get_num_entities();
            if(num_cells != 0)
            {
              throw Exception("TargetSet refinement not implemented for Hexahedra");
            }

            // okay, no hexahedra to be refined...
            return 0;
          }
        }; // TargetRefiner<Shape::Hypercube<3>,...>
      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_HYPERCUBE_HPP
