// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_ANALYTIC_PROJECTOR_HPP
#define KERNEL_ASSEMBLY_ANALYTIC_PROJECTOR_HPP 1

// includes, FEAT
#include <kernel/assembly/base.hpp>
#include <kernel/analytic/function.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/trafo/base.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Analytic vertex projector class
     *
     * This class interpolates an function in the vertices of the trafo's underlying mesh.
     *
     * \author Peter Zajac
     */
    class AnalyticVertexProjector
    {
    public:
      /**
       * \brief Projects an analytic function into the vertices.
       *
       * \param[out] vector
       * A reference to a vector object that shall receive the vertex interpolation of the analytic function.
       *
       * \param[in] function
       * A reference to the analytic function to be projected
       *
       * \param[in] trafo
       * A reference to the transformation defined on the mesh where to project the function.
       */
      template<typename VectorOut_, typename Function_, typename Trafo_>
      static void project(VectorOut_& vector, const Function_& function, const Trafo_& trafo)
      {
        // our mesh and shape types
        typedef Trafo_ TrafoType;
        typedef typename TrafoType::MeshType MeshType;

        // get our data and value types
        typedef typename VectorOut_::DataType DataType;

        // get the mesh
        const MeshType& mesh(trafo.get_mesh());

        // fetch the vertex set
        typedef typename MeshType::VertexSetType VertexSetType;
        const VertexSetType& vtx(mesh.get_vertex_set());

        // fetch the cell count
        const Index num_verts(mesh.get_num_entities(0));

        // create a clear output vector
        vector = VectorOut_(num_verts, DataType(0));

        // define the evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> EvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<EvalTraits> func_eval(function);

        // loop over all cells of the mesh
        for(Index i(0); i < num_verts; ++i)
        {
          // compute and store function value
          vector(i, func_eval.value(vtx[i]));

          // continue with next vertex
        }
      }
    }; // class AnalyticVertexProjector<...>

    /**
     * \brief Analytic cell projector class
     *
     * This class interpolates an analytic function in the cells of the trafo's underlying
     * mesh by computing the integral mean over the cell using a cubature rule.
     *
     * \author Peter Zajac
     */
    class AnalyticCellProjector
    {
    private:
      /// \cond internal
      static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
      /// \endcond

    public:
      /**
       * \brief Projects an analytic function into the cells using the barycentre cubature rule.
       *
       * \param[out] vector
       * A reference to a vector object that shall receive the vertex interpolation of the analytic function.
       *
       * \param[in] function
       * A reference to the analytic function to be projected
       *
       * \param[in] trafo
       * A reference to the transformation defined on the mesh where to project the function.
       */
      template<
        typename VectorOut_,
        typename Function_,
        typename Trafo_>
      static void project(
        VectorOut_& vector,
        const Function_& function,
        const Trafo_& trafo)
      {
        Cubature::DynamicFactory cubature_factory("barycentre");
        project(vector, function, trafo, cubature_factory);
      }

      /**
       * \brief Projects an analytic function into the cells.
       *
       * \param[out] vector
       * A reference to a vector object that shall receive the vertex interpolation of the analytic function.
       *
       * \param[in] function
       * A reference to the analytic function to be projected
       *
       * \param[in] trafo
       * A reference to the transformation defined on the mesh where to project the function.
       *
       * \param[in] cubature_factory
       * The cubature factory that is to be used for integration.
       */
      template<
        typename VectorOut_,
        typename Function_,
        typename Trafo_,
        typename CubatureFactory_>
      static void project(
        VectorOut_& vector,
        const Function_& function,
        const Trafo_& trafo,
        const CubatureFactory_& cubature_factory)
      {
        // our mesh and shape types
        typedef Trafo_ TrafoType;
        typedef typename TrafoType::MeshType MeshType;
        typedef typename MeshType::ShapeType ShapeType;

        // our shape dimension
        static constexpr int shape_dim = ShapeType::dimension;

        // get our data and value types
        typedef typename VectorOut_::DataType DataType;

        // create our cubature rule
        Cubature::Rule<ShapeType, DataType, DataType, Tiny::Vector<DataType, shape_dim>>
          cubature_rule(Cubature::ctor_factory, cubature_factory);

        // get our mesh
        const MeshType& mesh(trafo.get_mesh());

        // fetch the cell count
        const Index num_cells(mesh.get_num_entities(shape_dim));

        // create a clear output vector
        vector = VectorOut_(num_cells, DataType(0));

        // create a trafo evaluator
        typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;
        typedef typename TrafoEvaluator::template ConfigTraits<trafo_config> TrafoConfigTraits;
        typename TrafoConfigTraits::EvalDataType trafo_data;

        // create a trafo evaluator
        TrafoEvaluator trafo_eval(trafo);

        // define the evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> FuncEvalTraits;
        typedef typename FuncEvalTraits::ValueType ValueType;

        // create a function evaluator
        typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

        // loop over all cells of the mesh
        for(Index cell(0); cell < num_cells; ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          ValueType value(DataType(0));
          DataType area(DataType(0));

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute function value
            ValueType val = func_eval.value(trafo_data.img_point);

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

          // finish trafo evaluator
          trafo_eval.finish();

          // continue with next cell
        }
      }
    }; // class AnalyticCellProjector<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_ANALYTIC_PROJECTOR_HPP
