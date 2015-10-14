#pragma once
#ifndef KERNEL_ASSEMBLY_DISCRETE_PROJECTOR_HPP
#define KERNEL_ASSEMBLY_DISCRETE_PROJECTOR_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Discrete vertex projector class
     *
     * This class interpolates a discrete finite-element function in the vertices of the space's underyling
     * mesh.
     *
     * \author Peter Zajac
     */
    class DiscreteVertexProjector
    {
    private:
      /// \cond internal
      struct SpaceConfig :
        public Space::ConfigBase
      {
        static constexpr bool need_value = true;
      };
      /// \endcond

    public:
      /**
       * \brief Projects a discrete function into the vertices.
       *
       * \param[out] vector
       * A reference to a vector object that shall receive the vertex interpolation of the discrete function.
       *
       * \param[in] coeff
       * A reference to the coefficient vector of the finite-element function.
       *
       * \param[in] space
       * A reference to the finite-element space.
       */
      template<
        typename VectorOut_,
        typename VectorIn_,
        typename Space_>
      static void project(
        VectorOut_& vector,
        const VectorIn_& coeff,
        const Space_& space)
      {
        typedef Space_ SpaceType;
        typedef typename SpaceType::TrafoType TrafoType;
        typedef typename TrafoType::MeshType MeshType;
        typedef typename MeshType::ShapeType ShapeType;

        static constexpr int shape_dim = ShapeType::dimension;
        static constexpr int nverts = Shape::FaceTraits<ShapeType, 0>::count;

        // define assembly traits
        typedef AsmTraits1<typename VectorOut_::DataType, Space_, Trafo::ConfigBase, SpaceConfig> AsmTraits;
        typedef typename AsmTraits::DataType DataType;

        // define the reference cell type
        typedef Shape::ReferenceCell<ShapeType> RefCell;

        // fetch the trafo and the mesh
        const TrafoType& trafo(space.get_trafo());
        const MeshType& mesh(space.get_mesh());

        // fetch the index set
        typedef typename MeshType::template IndexSet<shape_dim, 0>::Type IndexSetType;
        const IndexSetType& vert_idx(mesh.template get_index_set<shape_dim, 0>());

        // fetch the cell count
        const Index num_verts(mesh.get_num_entities(0));

        // create a clear output vector
        vector = VectorOut_(num_verts, DataType(0));

        // allocate an auxiliary count array
        int* aux = new int[num_verts];
        for(Index i(0); i < num_verts; ++i)
        {
          aux[i] = 0;
        }

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local vector data
        typename AsmTraits::LocalVectorType lvad;

        // create a vector gather-axpy
        typename VectorIn_::GatherAxpy gather_axpy(coeff);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // format local matrix
          lvad.format();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // fetch local vector
          gather_axpy(lvad, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all vertices of the cell
          for(int k(0); k < nverts; ++k)
          {
            typename AsmTraits::DomainPointType dom_point;

            // initialise domain point
            for(int i(0); i < shape_dim; ++i)
            {
              dom_point[i] = DataType(RefCell::coord(k, i));
            }

            // compute trafo data
            trafo_eval(trafo_data, dom_point);

            // compute basis function data
            space_eval(space_data, trafo_data);

            // compute function value
            DataType value(DataType(0));

            // basis function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate fe function
              value += lvad(i) * space_data.phi[i].value;
              // continue with next basis function
            }

            // fetch the vertex index
            Index vi = vert_idx(cell, Index(k));

            // add vertex contribution
            vector(vi, vector(vi) + value);

            // update contribution counter
            ++aux[vi];

            // continue with next vertex
          }

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // finally, scale the output vector
        for(Index i(0); i < num_verts; ++i)
        {
          if(aux[i] > 1)
          {
            vector(i, vector(i) / DataType(aux[i]));
          }
        }

        // release auxiliary array
        delete [] aux;
      }
    }; // class DiscreteVertexProjector<...>

    /**
     * \brief Discrete cell projector class
     *
     * This class interpolates a discrete finite-element function in the cells of the space's underyling
     * mesh by computing the integral mean over the cell using a cubature rule.
     *
     * \author Peter Zajac
     */
    class DiscreteCellProjector
    {
    private:
      /// \cond internal
      struct SpaceConfig :
        public Space::ConfigBase
      {
        static constexpr bool need_value = true;
      };
      /// \endcond

    public:
      /**
       * \brief Projects a discrete function into the cells.
       *
       * \param[out] vector
       * A reference to a vector object that shall receive the cell interpolation of the discrete function.
       *
       * \param[in] coeff
       * A reference to the coefficient vector of the finite-element function.
       *
       * \param[in] space
       * A reference to the finite-element space.
       *
       * \param[in] cubature_factory
       * The cubature factory that is to be used for integration.
       */
      template<
        typename VectorOut_,
        typename VectorIn_,
        typename Space_,
        typename CubatureFactory_>
      static void project(
        VectorOut_& vector,
        const VectorIn_& coeff,
        const Space_& space,
        const CubatureFactory_& cubature_factory)
      {
        typedef Space_ SpaceType;
        typedef typename SpaceType::TrafoType TrafoType;
        typedef typename TrafoType::MeshType MeshType;
        typedef typename MeshType::ShapeType ShapeType;

        // define assembly traits
        typedef AsmTraits1<typename VectorOut_::DataType, SpaceType, Trafo::ConfigBase, SpaceConfig> AsmTraits;
        typedef typename AsmTraits::DataType DataType;

        // define the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo and the mesh
        const TrafoType& trafo(space.get_trafo());
        const MeshType& mesh(space.get_mesh());

        // fetch the cell count
        const Index num_cells(mesh.get_num_entities(ShapeType::dimension));

        // create a clear output vector
        vector = VectorOut_(num_cells, DataType(0));

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local vector data
        typename AsmTraits::LocalVectorType lvad;

        // create a vector gather-axpy
        typename VectorIn_::GatherAxpy gather_axpy(coeff);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // format local matrix
          lvad.format();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // fetch local vector
          gather_axpy(lvad, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // compute function value
          DataType value(DataType(0));
          DataType area(DataType(0));

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            DataType v(DataType(0));

            // basis function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate functor and integrate
              v += lvad(i) * space_data.phi[i].value;
              // continue with next basis function
            }

            // compute weight
            DataType weight(trafo_data.jac_det * cubature_rule.get_weight(k));

            // update cell area
            area += weight;

            // update cell value
            value += v * weight;

            // continue with next vertex
          }

          // set contribution
          vector(cell, value / area);

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }
      }
    }; // class DiscreteCellProjector<...>
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_DISCRETE_PROJECTOR_HPP
