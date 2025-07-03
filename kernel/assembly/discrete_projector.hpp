// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/assembly/asm_traits.hpp>

#include <vector>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Discrete vertex projector class
     *
     * This class interpolates a discrete finite-element function in the vertices of the space's underlying
     * mesh.
     *
     * \author Peter Zajac
     */
    class DiscreteVertexProjector
    {
    public:
      /**
       * \brief Projects a discrete function into the vertices.
       *
       * \param[out] vector
       * A \transient reference to a vector object that shall receive the vertex interpolation of the discrete function.
       *
       * \param[in] coeff
       * A \transient reference to the coefficient vector of the finite-element function.
       *
       * \param[in] space
       * A \transient reference to the finite-element space.
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
        typedef AsmTraits1<typename VectorOut_::DataType, Space_, TrafoTags::none, SpaceTags::value> AsmTraits;
        typedef typename AsmTraits::DataType DataType;

        // get our value type
        typedef typename VectorOut_::ValueType ValueType;

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
        std::vector<int> aux(num_verts, 0);

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
        typename AsmTraits::template TLocalVector<ValueType> loc_vec;

        // create a vector gather-axpy
        typename VectorIn_::GatherAxpy gather_axpy(coeff);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // format local matrix
          loc_vec.format();

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch local vector
          gather_axpy(loc_vec, dof_mapping);

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

            // initialize domain point
            for(int i(0); i < shape_dim; ++i)
            {
              dom_point[i] = Shape::ReferenceCell<ShapeType>::template vertex<DataType>(k, i);
            }

            // compute trafo data
            trafo_eval(trafo_data, dom_point);

            // compute basis function data
            space_eval(space_data, trafo_data);

            // compute function value
            ValueType value(DataType(0));

            // basis function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate fe function
              Tiny::axpy(value, loc_vec[i], space_data.phi[i].value);

              // continue with next basis function
            }

            // fetch the vertex index
            Index vi = vert_idx(cell, k);

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
            vector(i, (DataType(1)  / DataType(aux[i])) * vector(i));
          }
        }
      }

      /**
       * \brief Projects the gradient of a scalar discrete function into the vertices.
       *
       * \param[out] vector
       * A \transient reference to a vector object that shall receive the vertex interpolation
       * of the discrete function's gradient. Must be an instance of LAFEM::DenseVectorBlocked.
       *
       * \param[in] coeff
       * A \transient reference to the coefficient vector of the finite-element function.
       * Must be an instance of LAFEM::DenseVector
       *
       * \param[in] space
       * A \transient reference to the finite-element space.
       */
      template<
        typename VectorOut_,
        typename VectorIn_,
        typename Space_>
      static void project_gradient(
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
        typedef AsmTraits1<typename VectorOut_::DataType, Space_, TrafoTags::none, SpaceTags::grad> AsmTraits;
        typedef typename AsmTraits::DataType DataType;

        // get our value type
        typedef typename VectorOut_::ValueType ValueOutType;
        typedef typename VectorIn_::ValueType ValueInType;

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
        std::vector<int> aux(num_verts, 0);

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
        typename AsmTraits::template TLocalVector<ValueInType> loc_vec;

        // create a vector gather-axpy
        typename VectorIn_::GatherAxpy gather_axpy(coeff);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // format local matrix
          loc_vec.format();

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch local vector
          gather_axpy(loc_vec, dof_mapping);

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

            // initialize domain point
            for(int i(0); i < shape_dim; ++i)
            {
              dom_point[i] = Shape::ReferenceCell<ShapeType>::template vertex<DataType>(k, i);
            }

            // compute trafo data
            trafo_eval(trafo_data, dom_point);

            // compute basis function data
            space_eval(space_data, trafo_data);

            // compute function value
            ValueOutType value(DataType(0));

            // basis function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate fe function
              Tiny::axpy(value, space_data.phi[i].grad, loc_vec[i]);

              // continue with next basis function
            }

            // fetch the vertex index
            Index vi = vert_idx(cell, k);

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

        ValueOutType* vv = vector.elements();

        // finally, scale the output vector
        for(Index i(0); i < num_verts; ++i)
        {
          if(aux[i] > 1)
          {
            vv[i] *= (DataType(1)  / DataType(aux[i]));
          }
        }
      }
    }; // class DiscreteVertexProjector<...>

    /**
     * \brief Discrete cell projector class
     *
     * This class interpolates a discrete finite-element function in the cells of the space's underlying
     * mesh by computing the integral mean over the cell using a cubature rule.
     *
     * \author Peter Zajac
     */
    class DiscreteCellProjector
    {
    public:
      /**
       * \brief Projects a discrete function into the cells using the barycentre cubature rule.
       *
       * \param[out] vector
       * A \transient reference to a vector object that shall receive the cell interpolation of the discrete function.
       *
       * \param[in] coeff
       * A \transient reference to the coefficient vector of the finite-element function.
       *
       * \param[in] space
       * A \transient reference to the finite-element space.
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
        Cubature::DynamicFactory cubature_factory("barycentre");
        project(vector, coeff, space, cubature_factory);
      }

      /**
       * \brief Projects a discrete function into the cells.
       *
       * \param[out] vector
       * A \transient reference to a vector object that shall receive the cell interpolation of the discrete function.
       *
       * \param[in] coeff
       * A \transient reference to the coefficient vector of the finite-element function.
       *
       * \param[in] space
       * A \transient reference to the finite-element space.
       *
       * \param[in] cubature_factory
       * The \transient cubature factory that is to be used for integration.
       */
      template<
        typename VectorOut_,
        typename VectorIn_,
        typename Space_>
      static void project(
        VectorOut_& vector,
        const VectorIn_& coeff,
        const Space_& space,
        const String& cubature_name)
      {
        Cubature::DynamicFactory cubature_factory(cubature_name);
        project(vector, coeff, space, cubature_factory);
      }

      /**
       * \brief Projects a discrete function into the cells.
       *
       * \param[out] vector
       * A \transient reference to a vector object that shall receive the cell interpolation of the discrete function.
       *
       * \param[in] coeff
       * A \transient reference to the coefficient vector of the finite-element function.
       *
       * \param[in] space
       * A \transient reference to the finite-element space.
       *
       * \param[in] cubature_factory
       * The \transient cubature factory that is to be used for integration.
       */
      template<typename VectorOut_, typename VectorIn_, typename Space_>
      static void project(
        VectorOut_& vector,
        const VectorIn_& coeff,
        const Space_& space,
        const Cubature::DynamicFactory& cubature_factory)
      {
        typedef Space_ SpaceType;
        typedef typename SpaceType::TrafoType TrafoType;
        typedef typename TrafoType::MeshType MeshType;
        typedef typename MeshType::ShapeType ShapeType;

        // define assembly traits
        typedef AsmTraits1<typename VectorOut_::DataType, SpaceType, TrafoTags::jac_det, SpaceTags::value> AsmTraits;
        typedef typename AsmTraits::DataType DataType;

        // get our value type
        typedef typename VectorOut_::ValueType ValueType;

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
        typename AsmTraits::template TLocalVector<ValueType> loc_vec;

        // create a vector gather-axpy
        typename VectorIn_::GatherAxpy gather_axpy(coeff);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // format local matrix
          loc_vec.format();

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch local vector
          gather_axpy(loc_vec, dof_mapping);

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

            ValueType val(DataType(0));

            // basis function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate functor and integrate
              Tiny::axpy(val, loc_vec[i], space_data.phi[i].value);
              // continue with next basis function
            }

            // compute weight
            DataType weight(trafo_data.jac_det * cubature_rule.get_weight(k));

            // update cell area
            area += weight;

            // update cell value
            value += weight * val;

            // continue with next vertex
          }

          // set contribution
          vector(cell, (DataType(1) / area) * value);

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }
      }

      /**
       * \brief Projects the gradient of a scalar discrete function into the cells.
       *
       * \param[out] vector
       * A \transient reference to a vector object that shall receive the cell interpolation
       * of the discrete function's gradient.
       *
       * \param[in] coeff
       * A \transient reference to the coefficient vector of the finite-element function.
       *
       * \param[in] space
       * A \transient reference to the finite-element space.
       *
       * \param[in] cubature_factory
       * The \transient cubature factory that is to be used for integration.
       */
      template<
        typename VectorOut_,
        typename VectorIn_,
        typename Space_>
      static void project_gradient(
        VectorOut_& vector,
        const VectorIn_& coeff,
        const Space_& space,
        const String& cubature_name)
      {
        Cubature::DynamicFactory cubature_factory(cubature_name);
        project_gradient(vector, coeff, space, cubature_factory);
      }

      /**
      * \brief Projects the gradient of a scalar discrete function into the cells.
      *
      * \param[out] vector
      * A \transient reference to a vector object that shall receive the cell interpolation of the discrete function.
      *
      * \param[in] coeff
      * A \transient reference to the coefficient vector of the finite-element function.
      *
      * \param[in] space
      * A \transient reference to the finite-element space.
      *
      * \param[in] cubature_factory
      * The \transient cubature factory that is to be used for integration.
      */
      template<typename VectorOut_, typename VectorIn_, typename Space_>
      static void project_gradient(
        VectorOut_& vector,
        const VectorIn_& coeff,
        const Space_& space,
        const Cubature::DynamicFactory& cubature_factory)
      {
        typedef Space_ SpaceType;
        typedef typename SpaceType::TrafoType TrafoType;
        typedef typename TrafoType::MeshType MeshType;
        typedef typename MeshType::ShapeType ShapeType;

        // define assembly traits
        typedef AsmTraits1<typename VectorOut_::DataType, SpaceType, TrafoTags::jac_det, SpaceTags::grad> AsmTraits;
        typedef typename AsmTraits::DataType DataType;

        // get our value type
        typedef typename VectorIn_::ValueType ValueInType;
        typedef typename VectorOut_::ValueType ValueOutType;

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
        typename AsmTraits::template TLocalVector<ValueInType> loc_vec;

        // create a vector gather-axpy
        typename VectorIn_::GatherAxpy gather_axpy(coeff);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // format local matrix
          loc_vec.format();

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch local vector
          gather_axpy(loc_vec, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // compute function value
          ValueOutType value(DataType(0));
          DataType area(DataType(0));

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            ValueOutType val(DataType(0));

            // basis function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate functor and integrate
              Tiny::axpy(val, space_data.phi[i].grad, loc_vec[i]);
              // continue with next basis function
            }

            // compute weight
            DataType weight(trafo_data.jac_det * cubature_rule.get_weight(k));

            // update cell area
            area += weight;

            // update cell value
            Tiny::axpy(value, val, weight);

            // continue with next vertex
          }

          // set contribution
          vector(cell, (DataType(1) / area) * value);

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }
      }

      /**
       * \brief Projects a discrete function into the cells of a once refined mesh.
       *
       * \param[out] vector
       * A \transient reference to a vector object that shall receive the refined cell interpolation of the discrete function.
       *
       * \param[in] coeff
       * A \transient reference to the coefficient vector of the finite-element function.
       *
       * \param[in] space
       * A \transient reference to the finite-element space.
       */
      template<typename VectorOut_, typename VectorIn_, typename Space_>
      static void project_refined(VectorOut_& vector, const VectorIn_& coeff, const Space_& space)
      {
        typedef Space_ SpaceType;
        typedef typename SpaceType::TrafoType TrafoType;
        typedef typename TrafoType::MeshType MeshType;
        typedef typename MeshType::ShapeType ShapeType;

        // define assembly traits
        typedef AsmTraits1<typename VectorOut_::DataType, SpaceType, TrafoTags::jac_det, SpaceTags::value> AsmTraits;
        typedef typename AsmTraits::DataType DataType;

        // get our value type
        typedef typename VectorOut_::ValueType ValueType;

        // define the cubature rule
        Cubature::DynamicFactory cubature_factory("refine:midpoint");
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo and the mesh
        const TrafoType& trafo(space.get_trafo());
        const MeshType& mesh(space.get_mesh());

        // fetch the cell count
        const Index num_cells(mesh.get_num_entities(ShapeType::dimension));
        const int num_points = cubature_rule.get_num_points();

        // create a clear output vector
        vector = VectorOut_(num_cells * Index(num_points), DataType(0));

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
        typename AsmTraits::template TLocalVector<ValueType> loc_vec;

        // create a vector gather-axpy
        typename VectorIn_::GatherAxpy gather_axpy(coeff);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // format local matrix
          loc_vec.format();

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch local vector
          gather_axpy(loc_vec, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            ValueType value(DataType(0));

            // basis function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate functor and integrate
              Tiny::axpy(value, loc_vec[i], space_data.phi[i].value);
              // continue with next basis function
            }

            // save value
            vector(cell*Index(num_points) + Index(k), value);

            // continue with next vertex
          }

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }
      }
    }; // class DiscreteCellProjector<...>
  } // namespace Assembly
} // namespace FEAT
