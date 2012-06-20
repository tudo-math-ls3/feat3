#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_SIMPLEX_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_SIMPLEX_HPP 1

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
        struct TargetRefiner<Shape::Simplex<shape_dim_>, 0>
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
        }; // TargetRefiner<Shape::Simplex<...>, 0>

        template<>
        struct TargetRefiner<Shape::Simplex<2>, 0>
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
        }; // TargetRefiner<Shape::Simplex<...>, 0> for triangles

        template<>
        struct TargetRefiner<Shape::Simplex<1>, 1>
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
        }; // TargetRefiner<Shape::Simplex<1>,1>

        template<>
        struct TargetRefiner<Shape::Simplex<2>, 1>
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
            typedef Congruency::TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

            // fetch index offset
            const Index ioe = index_offsets[2];

            // typedef index set types
            typedef IndexSet<3> IndexSetTypeTV;

            // typedef index vector type
            typedef IndexSetTypeTV::IndexVectorType IndexVectorTypeTV;

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
              const IndexVectorTypeTV& trg_t_v = index_trg_t_v[trg_t];

              // fetch source edge-vert index vector
              const IndexVectorTypeTV& src_t_v = index_src_t_v[i];

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
        }; // TargetRefiner<Shape::Simplex<2>,1>

        template<>
        struct TargetRefiner<Shape::Simplex<2>, 2>
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
            typedef Congruency::TargetIndexMapping<ShapeType, 0> TargetIndexMappingType;

            // fetch index offset
            const Index iot = index_offsets[2];

            // typedef index set types
            typedef IndexSet<3> IndexSetTypeTV;

            // typedef index vector type
            typedef IndexSetTypeTV::IndexVectorType IndexVectorTypeTV;

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
              const IndexVectorTypeTV& trg_t_v = index_trg_t_v[trg_t];

              // fetch source edge-vert index vector
              const IndexVectorTypeTV& src_t_v = index_src_t_v[i];

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
        }; // TargetRefiner<Shape::Simplex<2>,2>

        template<int cell_dim_>
        struct TargetRefiner<Shape::Simplex<3>, cell_dim_>
        {
          typedef Shape::Simplex<3> ShapeType;
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
            // We do not have an implementation of 3D tetrahedron submeshes yet due to the lack of a
            // tetrahedral congruency sampler (see Congruency::Sampler<...> class template) which is needed
            // to implement this function. However, this function will also be called when refining a
            // 2D submesh of a 3D tetrahedral mesh -- in this case the number of tetrahedra is zero,
            // so we'll check for that to ensure that we don't prohibit this valid scenario.
            Index num_cells = target_set_holder_in.get_target_set<3>().get_num_entities();
            if(num_cells != 0)
            {
              throw InternalError("TargetSet refinement not implemented for Tetrahedra");
            }

            // okay, no tetrahedron to be refined...
            return 0;
          }
        }; // TargetRefiner<Shape::Simplex<3>,...>

        template<>
        struct TargetRefiner<Shape::Simplex<3>,0>
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
        }; // TargetRefiner<Shape::Simplex<3>,...>

      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_SIMPLEX_HPP
