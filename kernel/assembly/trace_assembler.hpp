// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_TRACE_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_TRACE_ASSEMBLER_HPP 1

// includes, FEAT
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/geometry/intern/face_ref_trafo.hpp>
#include <kernel/geometry/intern/congruency_sampler.hpp>
#include <kernel/geometry/intern/congruency_trafo.hpp>
#include <kernel/geometry/intern/index_representative.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

#include <vector>

namespace FEAT
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Outer_, typename Shape_, int cell_dim_>
      class CompIndexMap
      {
      protected:
        const Outer_& _outer;
        int _idx;

      public:
        explicit CompIndexMap(const Outer_& outer, int idx) : _outer(outer), _idx(idx) {}

        Index operator[](int i) const
        {
          return _outer[Geometry::Intern::FaceIndexMapping<Shape_, cell_dim_, 0>::map(_idx,i)];
        }
      }; // class CompIndexMap

      template<int max_>
      class CommonDofMap
      {
      public:
        Index glob_dofs[max_];
        int loc1_dofs[max_];
        int loc2_dofs[max_];
        int num_dofs;

        explicit CommonDofMap() : num_dofs(0) {}

        void clear()
        {
          num_dofs = 0;
        }

        int push_1(Index dof, int loc)
        {
          for(int i(0); i < num_dofs; ++i)
          {
            if(dof == glob_dofs[i])
            {
              loc1_dofs[i] = loc;
              return i;
            }
          }
          ASSERT(num_dofs < max_);
          glob_dofs[num_dofs] = dof;
          loc1_dofs[num_dofs] = loc;
          loc2_dofs[num_dofs] = -1;
          return num_dofs++;
        }

        int push_2(Index dof, int loc)
        {
          for(int i(0); i < num_dofs; ++i)
          {
            if(dof == glob_dofs[i])
            {
              loc2_dofs[i] = loc;
              return i;
            }
          }
          ASSERT(num_dofs < max_);
          glob_dofs[num_dofs] = dof;
          loc1_dofs[num_dofs] = -1;
          loc2_dofs[num_dofs] = loc;
          return num_dofs++;
        }

        int loc_1(int k) const
        {
          return loc1_dofs[k];
        }

        int loc_2(int k) const
        {
          return loc2_dofs[k];
        }

        int get_num_local_dofs() const
        {
          return num_dofs;
        }

        Index get_index(int i) const
        {
          return glob_dofs[i];
        }
      }; // class CommonDofMap<...>

    } // namespace Intern
    /// \endcond

    /**
     * \brief Assembles Facet information for assembling operators
     *
     * This functionality similar to the OuterNormalComputer (which assembles only boundary facets, but can compute
     * weighted outer normals, too).
     *
     * \tparam Trafo_
     * The transformation
     *
     * \note After all MeshParts are added, compile_facets() needs to be called
     */
    template<typename Trafo_>
    class TraceAssembler
    {
    public:
      typedef Trafo_ TrafoType;
      typedef typename TrafoType::MeshType MeshType;
      typedef typename TrafoType::ShapeType ShapeType;

      static constexpr int shape_dim = ShapeType::dimension;
      static constexpr int facet_dim = shape_dim-1;

    protected:
      const TrafoType& _trafo;
      std::vector<int> _facet_mask, _cell_facet, _facet_ori;
      std::vector<Index> _cells, _facets;

    public:
      explicit TraceAssembler(const TrafoType& trafo) :
        _trafo(trafo),
        _facet_mask(trafo.get_mesh().get_num_entities(facet_dim), 0)
      {
      }

      void clear_facets()
      {
        _cell_facet.clear();
        _facet_ori.clear();
        _cells.clear();
        _facets.clear();
        for(auto& f : _facets)
          f = Index(0);
      }

      void add_facet(Index ifacet)
      {
        ASSERTM(ifacet < Index(_facets.size()), "invalid facet index");
        _facet_mask.at(ifacet) = 1;
      }

      void add_mesh_part(const Geometry::MeshPart<MeshType>& mesh_part)
      {
        const auto& trg = mesh_part.template get_target_set<facet_dim>();
        for(Index i(0); i < trg.get_num_entities(); ++i)
          _facet_mask.at(trg[i]) = 1;
      }

      void compile_facets(bool only_boundary = true)
      {
        _cell_facet.clear();
        _facet_ori.clear();
        _cells.clear();
        _facets.clear();

        // build elements-at-facet graph
        Adjacency::Graph elem_at_facet(Adjacency::RenderType::injectify_transpose,
          _trafo.get_mesh().template get_index_set<shape_dim, facet_dim>());

        // loop over all facets
        for(Index iface(0); iface < Index(_facet_mask.size()); ++iface)
        {
          // use this facet?
          if(_facet_mask[iface] == 0)
            continue;

          // ensure that this is a boundary facet if required
          if(only_boundary && (elem_at_facet.degree(iface) != Index(1)))
            throw InternalError(__func__, __FILE__, __LINE__, "facet is adjacent to more than 1 element");

          // add all elements
          for(auto it = elem_at_facet.image_begin(iface); it != elem_at_facet.image_end(iface); ++it)
          {
            Index icell = *it;

            // try to compute local facet index and orientation
            int loc_face(0), face_ori(0);
            if(!find_local_facet(iface, icell, loc_face, face_ori))
              throw InternalError(__func__, __FILE__, __LINE__, "failed to find local facet");

            // alright, add this facet to our list
            _facets.push_back(iface);
            _cells.push_back(icell);
            _cell_facet.push_back(loc_face);
            _facet_ori.push_back(face_ori);
          }
        }
      }

      void compile_all_facets(bool inner, bool outer)
      {
        _cell_facet.clear();
        _facet_ori.clear();
        _cells.clear();
        _facets.clear();

        // build elements-at-facet graph
        Adjacency::Graph elem_at_facet(Adjacency::RenderType::injectify_transpose,
          _trafo.get_mesh().template get_index_set<shape_dim, facet_dim>());

        // loop over all facets
        for(Index iface(0); iface < Index(_facet_mask.size()); ++iface)
        {
          const int degree = int(elem_at_facet.degree(iface));
          XASSERT(degree > 0);
          XASSERT(degree < 3);

          if((degree == 1) && !outer)
            continue;
          if((degree == 2) && !inner)
            continue;

          // add all elements
          for(auto it = elem_at_facet.image_begin(iface); it != elem_at_facet.image_end(iface); ++it)
          {
            Index icell = *it;

            // try to compute local facet index and orientation
            int loc_face(0), face_ori(0);
            if(!find_local_facet(iface, icell, loc_face, face_ori))
              throw InternalError(__func__, __FILE__, __LINE__, "failed to find local facet");

            // alright, add this facet to our list
            _facets.push_back(iface);
            _cells.push_back(icell);
            _cell_facet.push_back(loc_face);
            _facet_ori.push_back(face_ori);
          }
        }
      }

      template<
        typename Matrix_,
        typename Operator_,
        typename Space_,
        typename CubatureFactory_>
      void assemble_operator_matrix1(
        Matrix_& matrix,
        Operator_& operat,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1)) const
      {
        assemble_operator_matrix2(matrix, operat, space, space, cubature_factory, alpha);
      }

      template<
        typename Matrix_,
        typename Operator_,
        typename TestSpace_,
        typename TrialSpace_,
        typename CubatureFactory_>
      void assemble_operator_matrix2(
        Matrix_& matrix,
        Operator_& operat,
        const TestSpace_& test_space,
        const TrialSpace_& trial_space,
        const CubatureFactory_& cubature_factory,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1)) const
      {
        // validate matrix dimensions
        XASSERTM(matrix.rows() == test_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == trial_space.get_num_dofs(), "invalid matrix dimensions");

        // matrix type
        typedef Matrix_ MatrixType;
        // operator type
        typedef Operator_ OperatorType;
        // test-space type
        typedef TestSpace_ TestSpaceType;
        // trial-space type
        typedef TrialSpace_ TrialSpaceType;

        // assembly traits
        typedef AsmTraits2<
          typename MatrixType::DataType,
          TestSpaceType,
          TrialSpaceType,
          OperatorType::trafo_config,
          OperatorType::test_config,
          OperatorType::trial_config> AsmTraits;

        typedef typename AsmTraits::DataType DataType;

        // shape types
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

        // fetch the trafo
        const TrafoType& trafo = test_space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a trafo facet evaluator
        typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;
        TrafoFacetEvaluator trafo_facet_eval(trafo);

        typedef typename TrafoFacetEvaluator::template ConfigTraits<TrafoTags::jac_det>::EvalDataType TrafoFacetEvalData;

        // create space evaluators
        typename AsmTraits::TestEvaluator test_eval(test_space);
        typename AsmTraits::TrialEvaluator trial_eval(trial_space);

        // create dof-mappings
        typename AsmTraits::TestDofMapping test_dof_mapping(test_space);
        typename AsmTraits::TrialDofMapping trial_dof_mapping(trial_space);

        // create a operator evaluator
        typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        TrafoFacetEvalData trafo_facet_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;
        typename AsmTraits::TrialEvalData trial_data;

        // create local matrix data
        typename AsmTraits::LocalMatrixType lmd;

        // create cubature rule
        typename Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_axpy(matrix);

        // trafo matrices and vectors
        Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat;
        Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat;
        Tiny::Vector<DataType, shape_dim> face_vec;
        Tiny::Vector<DataType, facet_dim> ori_vec;

        face_mat.format();
        ori_mat.format();
        face_vec.format();
        ori_vec.format();

        // loop over all cells of the mesh
        for(Index f(0); f < Index(_facets.size()); ++f)
        {
          // get facet index
          const Index face = _facets[f];
          const Index cell = _cells[f];

          // compute facet trafos
          Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat, face_vec, _cell_facet[f]);

          // compute orientation trafos
          Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat, ori_vec, _facet_ori[f]);

          // prepare trafo evaluators
          trafo_facet_eval.prepare(face);
          trafo_eval.prepare(cell);

          // prepare space evaluators
          test_eval.prepare(trafo_eval);
          trial_eval.prepare(trafo_eval);

          // prepare operator evaluator
          oper_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_test_dofs = test_eval.get_num_local_dofs();
          int num_loc_trial_dofs = trial_eval.get_num_local_dofs();

          // format local matrix
          lmd.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // get cubature point
            auto cub_pt = cubature_rule.get_point(k);

            // transform to local facet
            auto cub_cf = (face_mat * ((ori_mat * cub_pt) + ori_vec)) + face_vec;

            // compute trafo data
            trafo_facet_eval(trafo_facet_data, cub_pt);
            trafo_eval(trafo_data, cub_cf);

            // compute basis function data
            test_eval(test_data, trafo_data);
            trial_eval(trial_data, trafo_data);

            // prepare bilinear operator
            oper_eval.set_point(trafo_data);

            // test function loop
            for(int i(0); i < num_loc_test_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_trial_dofs; ++j)
              {
                // evaluate operator and integrate
                lmd(i,j) += trafo_facet_data.jac_det * cubature_rule.get_weight(k) *
                  oper_eval(trial_data.phi[j], test_data.phi[i]);
                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }

          // finish operator evaluator
          oper_eval.finish();

          // finish evaluators
          trial_eval.finish();
          test_eval.finish();
          trafo_eval.finish();
          trafo_facet_eval.finish();

          // initialise dof-mappings
          test_dof_mapping.prepare(cell);
          trial_dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(lmd, test_dof_mapping, trial_dof_mapping, alpha);

          // finish dof mapping
          trial_dof_mapping.finish();
          test_dof_mapping.finish();

          // continue with next cell
        }

        // okay, that's it
      }

      template<
        typename Vector_,
        typename Functional_,
        typename CubatureFactory_,
        typename Space_>
      void assemble_functional_vector(
        Vector_& vector,
        const Functional_& functional,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        typename Vector_::DataType alpha = typename Vector_::DataType(1)) const
      {
        // validate vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");

        // vector type
        typedef Vector_ VectorType;
        // functional type
        typedef Functional_ FunctionalType;
        // space type
        typedef Space_ SpaceType;

        // assembly traits
        typedef AsmTraits1<
          typename VectorType::DataType,
          SpaceType,
          FunctionalType::trafo_config,
          FunctionalType::test_config> AsmTraits;

        typedef typename AsmTraits::DataType DataType;

        // shape types
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

        // fetch the trafo
        const TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a trafo facet evaluator
        typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;
        TrafoFacetEvaluator trafo_facet_eval(trafo);

        typedef typename TrafoFacetEvaluator::template ConfigTraits<TrafoTags::jac_det>::EvalDataType TrafoFacetEvalData;

        // create a space evaluator and evaluation data
        typename AsmTraits::TestEvaluator test_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create a functional evaluator
        typename FunctionalType::template Evaluator<AsmTraits> func_eval(functional);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        TrafoFacetEvalData trafo_facet_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;

        // create local vector data
        typename AsmTraits::LocalVectorType lvad;

        // create cubature rule
        typename Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename VectorType::ScatterAxpy scatter_axpy(vector);

        // trafo matrices and vectors
        Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat;
        Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat;
        Tiny::Vector<DataType, shape_dim> face_vec;
        Tiny::Vector<DataType, facet_dim> ori_vec;

        face_mat.format();
        ori_mat.format();
        face_vec.format();
        ori_vec.format();

        // loop over all cells of the mesh
        for(Index f(0); f < Index(_facets.size()); ++f)
        {
          // get facet index
          const Index face = _facets[f];
          const Index cell = _cells[f];

          // compute facet trafos
          Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat, face_vec, _cell_facet[f]);

          // compute orientation trafos
          Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat, ori_vec, _facet_ori[f]);

          // prepare trafo evaluators
          trafo_facet_eval.prepare(face);
          trafo_eval.prepare(cell);

          // prepare test-space evaluator
          test_eval.prepare(trafo_eval);

          // prepare functional evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = test_eval.get_num_local_dofs();

          // format local matrix
          lvad.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // get cubature point
            auto cub_pt = cubature_rule.get_point(k);

            // transform to local facet
            auto cub_cf = (face_mat * ((ori_mat * cub_pt) + ori_vec)) + face_vec;

            // compute trafo data
            trafo_facet_eval(trafo_facet_data, cub_pt);
            trafo_eval(trafo_data, cub_cf);

            // compute test basis function data
            test_eval(test_data, trafo_data);

            // prepare functional
            func_eval.set_point(trafo_data);

            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate functional and integrate
              lvad(i) += trafo_facet_data.jac_det * cubature_rule.get_weight(k) * func_eval(test_data.phi[i]);
              // continue with next trial function
            }
            // continue with next test function
          }

          // finish functional evaluator
          func_eval.finish();

          // finish evaluators
          test_eval.finish();
          trafo_eval.finish();
          trafo_facet_eval.finish();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(lvad, dof_mapping, alpha);

          // finish dof-mapping
          dof_mapping.finish();

          // continue with next cell
        }
      }

      /**
       * \brief Assembles a flow accumulator.
       *
       * This function assembles a so-called accumulator on a sub-dimensional
       * mesh region, which can be used to assemble body forces on a boundary.
       *
       * \attention
       * This assembly function is somewhat provisional - use at own risk!
       *
       * The accumulator that is assembled by this function has to provide
       * an overloaded "operator()" with the following function parameters:
       * - cubature weight (scalar)
       * - mapped image point (Tiny::Vector)
       * - jacobi matrix (Tiny::Matrix)
       * - velocity value (Tiny::Vector)
       * - velocity gradient (Tiny::Matrix)
       * - pressure value (scalar)
       *
       * \param[inout] accum
       * The accumulator to be assembled. The "operator()" of this object
       * is called for each cubature point on each facet.
       *
       * \param[in] vector_v
       * The velocity vector.
       *
       * \param[in] vector_p
       * The pressure vector.
       *
       * \param[in] space_v
       * The velocity space.
       *
       * \param[in] space_p
       * The pressure space.
       *
       * \param[in] cubature_factory
       * The cubature factory that is to be used for integration.
       */
      template<
        typename Accum_,
        typename DataType_,
        typename IndexType_,
        int dim_,
        typename SpaceV_,
        typename SpaceP_,
        typename CubatureFactory_>
      void assemble_flow_accum(
        Accum_& accum,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector_v,
        const LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>& vector_p,
        const SpaceV_& space_v,
        const SpaceP_& space_p,
        const CubatureFactory_& cubature_factory)
      {
        // validate vector dimensions
        XASSERTM(vector_v.size() == space_v.get_num_dofs(), "invalid velocity vector size");
        XASSERTM(vector_p.size() == space_p.get_num_dofs(), "invalid pressure vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VeloVector;
        typedef LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> PresVector;

        // assembly traits
        typedef Assembly::AsmTraits2<
          DataType_,
          SpaceP_,
          SpaceV_,
          TrafoTags::none,
          SpaceTags::value,
          SpaceTags::value|SpaceTags::grad> AsmTraits;

        typedef typename AsmTraits::DataType DataType;

        // shape types
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

        // fetch the trafo
        const TrafoType& trafo = space_v.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a trafo facet evaluator
        typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;
        TrafoFacetEvaluator trafo_facet_eval(trafo);

        /// \compilerhack PGI and Intel(18) do not understand complex template statements
#if defined(FEAT_COMPILER_PGI) || (defined(FEAT_COMPILER_INTEL) && FEAT_COMPILER_INTEL >= 1800 && FEAT_COMPILER_INTEL < 1900)
        static constexpr TrafoTags trafo_facet_eval_tags = TrafoTags::img_point|TrafoTags::jac_det|TrafoTags::jac_mat;
        typedef typename TrafoFacetEvaluator::template ConfigTraits <trafo_facet_eval_tags>::EvalDataType TrafoFacetEvalData;
#else
        typedef typename TrafoFacetEvaluator::template ConfigTraits
          <TrafoTags::img_point|TrafoTags::jac_det|TrafoTags::jac_mat>::EvalDataType TrafoFacetEvalData;
#endif

        // create a space evaluator and evaluation data
        typename AsmTraits::TrialEvaluator space_eval_v(space_v);
        typename AsmTraits::TestEvaluator space_eval_p(space_p);

        // create a dof-mapping
        typename AsmTraits::TrialDofMapping dof_mapping_v(space_v);
        typename AsmTraits::TestDofMapping dof_mapping_p(space_p);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        TrafoFacetEvalData trafo_facet_data;

        // create space evaluation data
        typename AsmTraits::TrialEvalData space_data_v;
        typename AsmTraits::TestEvalData space_data_p;

        // create cubature rule
        typename Assembly::Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename VeloVector::GatherAxpy gather_v(vector_v);
        typename PresVector::GatherAxpy gather_p(vector_p);

        // get maximum number of local dofs
        static constexpr int max_local_dofs_v = AsmTraits::max_local_trial_dofs;
        static constexpr int max_local_dofs_p = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs_v> LocalVectorTypeV;
        typedef Tiny::Vector<DataType, max_local_dofs_p> LocalVectorTypeP;
        LocalVectorTypeV local_vector_v;
        LocalVectorTypeP local_vector_p;

        // our local velocity gradient
        Tiny::Vector<DataType, dim_> loc_value_v;
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        // trafo matrices and vectors
        Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat;
        Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat;
        Tiny::Vector<DataType, shape_dim> face_vec;
        Tiny::Vector<DataType, facet_dim> ori_vec;

        face_mat.format();
        ori_mat.format();
        face_vec.format();
        ori_vec.format();

        // loop over all cells of the mesh
        for(Index f(0); f < Index(this->_facets.size()); ++f)
        {
          // get facet index
          const Index face = this->_facets[f];
          const Index cell = this->_cells[f];

          // compute facet trafos
          Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat, face_vec, this->_cell_facet[f]);

          // compute orientation trafos
          Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat, ori_vec, this->_facet_ori[f]);

          // prepare trafo evaluators
          trafo_facet_eval.prepare(face);
          trafo_eval.prepare(cell);

          // prepare space evaluators
          space_eval_v.prepare(trafo_eval);
          space_eval_p.prepare(trafo_eval);

          // initialise dof-mappings
          dof_mapping_v.prepare(cell);
          dof_mapping_p.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs_v = space_eval_v.get_num_local_dofs();
          const int num_loc_dofs_p = space_eval_p.get_num_local_dofs();

          // gather our local velocity dofs
          local_vector_v.format();
          local_vector_p.format();
          gather_v(local_vector_v, dof_mapping_v);
          gather_p(local_vector_p, dof_mapping_p);

          // finish dof-mapping
          dof_mapping_p.finish();
          dof_mapping_v.finish();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // get cubature point
            auto cub_pt = cubature_rule.get_point(k);

            // transform to local facet
            auto cub_cf = (face_mat * ((ori_mat * cub_pt) + ori_vec)) + face_vec;

            // compute trafo data
            trafo_facet_eval(trafo_facet_data, cub_pt);
            trafo_eval(trafo_data, cub_cf);

            // compute test basis function data
            space_eval_v(space_data_v, trafo_data);
            space_eval_p(space_data_p, trafo_data);

            // compute local velocity value
            loc_value_v.format();
            for(int i(0); i < num_loc_dofs_v; ++i)
              loc_value_v.axpy(space_data_v.phi[i].value, local_vector_v[i]);

            // compute local velocity gradient
            loc_grad_v.format();
            for(int i(0); i < num_loc_dofs_v; ++i)
              loc_grad_v.add_outer_product(local_vector_v[i], space_data_v.phi[i].grad);

            // compute local pressure value
            DataType loc_value_p = DataType(0);
            for(int i(0); i < num_loc_dofs_p; ++i)
              loc_value_p += local_vector_p[i] * space_data_p.phi[i].value;

            // call accumulator
            accum(
              trafo_facet_data.jac_det * cubature_rule.get_weight(k),
              trafo_facet_data.img_point,
              trafo_facet_data.jac_mat,
              loc_value_v,
              loc_grad_v,
              loc_value_p
            );

            // continue with next basis function
          }

          // finish evaluators
          space_eval_p.finish();
          space_eval_v.finish();
          trafo_eval.finish();
          trafo_facet_eval.finish();

          // continue with next cell
        }
      }

      /**
       * \brief Assembles the jump-stabilisation operator onto a matrix.
       *
       * This function assembles the jump stabilisation operator:
       *   \f[J(\varphi,\psi) = \gamma \sum_E (s\cdot J_E)^{p} \int_E [\nabla \varphi]\cdot[\nabla\psi]\f]
       *
       * \attention
       * The matrix must have an extended stencil, which must have been assembled by calling
       * Assembly::SymbolicAssembler::assemble_matrix_ext1() !
       *
       * \param[inout] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] space
       * The finite element space to be used.
       *
       * \param[in] cubature_factory
       * The cubature for integration. Note that this is a cubature rule on the facets.
       *
       * \param[in] gamma
       * The scaling factor gamma for the jump stabilisation operator.
       *
       * \param[in] jacdet_scal
       * The scaling factor \e s for the Jacobian determinant factor.
       *
       * \param[in] jacdet_expo
       * The exponent \e p for the Jacobian determinant factor.
       */
      template<
        typename Matrix_,
        typename Space_,
        typename CubatureFactory_>
      void assemble_jump_stabil_operator_matrix(
        Matrix_& matrix,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        typename Matrix_::DataType gamma = typename Matrix_::DataType(1),
        typename Matrix_::DataType jacdet_scal = typename Matrix_::DataType(2),
        typename Matrix_::DataType jacdet_expo = typename Matrix_::DataType(2)) const
      {
        // validate matrix dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");

        // matrix type
        typedef Matrix_ MatrixType;
        // test-space type
        typedef Space_ SpaceType;

        // assembly traits
        typedef AsmTraits1<
          typename MatrixType::DataType,
          SpaceType,
          TrafoTags::none,
          SpaceTags::grad> AsmTraits;

        typedef typename AsmTraits::DataType DataType;

        // shape types
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

        // fetch the trafo
        const TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval_1(trafo), trafo_eval_2(trafo);

        // create a trafo facet evaluator
        typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;
        typedef typename TrafoFacetEvaluator::template ConfigTraits<TrafoTags::jac_det>::EvalDataType TrafoFacetEvalData;
        TrafoFacetEvaluator trafo_facet_eval(trafo);

        // create space evaluators
        typename AsmTraits::SpaceEvaluator space_eval_1(space), space_eval_2(space);

        // create dof-mappings
        typename AsmTraits::DofMapping dof_mapping_1(space), dof_mapping_2(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data_1, trafo_data_2;
        TrafoFacetEvalData trafo_facet_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data_1, space_data_2;

        // create cubature rule
        typename Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_axpy(matrix);

        // common DOF mapping
        static constexpr int max_common_dofs = 2 * AsmTraits::max_local_test_dofs;
        Intern::CommonDofMap<max_common_dofs> common_map;

        // jump gradients
        typename AsmTraits::SpaceEvalTraits::BasisGradientType jump_grad[max_common_dofs];

        // local matrix data
        Tiny::Matrix<DataType, max_common_dofs, max_common_dofs> loc_mat;

        // trafo matrices and vectors
        Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat_1, face_mat_2;
        Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat_1, ori_mat_2;
        Tiny::Vector<DataType, shape_dim> face_vec_1, face_vec_2;
        Tiny::Vector<DataType, facet_dim> ori_vec_1, ori_vec_2;

        face_mat_1.format();
        face_mat_2.format();
        ori_mat_1.format();
        ori_mat_2.format();
        face_vec_1.format();
        face_vec_2.format();
        ori_vec_1.format();
        ori_vec_2.format();

        // loop over all cells of the mesh
        for(std::size_t f(0); f < _facets.size(); /*++f*/)
        {
          // get facet index
          const Index face = _facets[f];
          const Index cell_1 = _cells[f];
          const int cell_facet_1 = _cell_facet[f];
          const int facet_ori_1 = _facet_ori[f];
          ++f;

          // check for next cell
          const bool inner = ((f < _facets.size()) && (face == _facets[f]));
          const Index cell_2 = (inner ? _cells[f] : cell_1);
          const int cell_facet_2 = (inner ? _cell_facet[f] : cell_facet_1);
          const int facet_ori_2 = (inner ? _facet_ori[f] : facet_ori_1);
          if(inner)
            ++f;

          // prepare dof mappings
          dof_mapping_1.prepare(cell_1);
          dof_mapping_2.prepare(cell_2);

          // build local dof map
          common_map.clear();
          for(int i(0); i < dof_mapping_1.get_num_local_dofs(); ++i)
            common_map.push_1(dof_mapping_1.get_index(i), i);
          if(inner)
          {
            for(int i(0); i < dof_mapping_2.get_num_local_dofs(); ++i)
              common_map.push_2(dof_mapping_2.get_index(i), i);
          }

          // finish dof mapping
          dof_mapping_2.finish();
          dof_mapping_1.finish();

          // get number of common local dofs
          const int num_local_dofs = common_map.get_num_local_dofs();

          // compute facet trafos
          Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat_1, face_vec_1, cell_facet_1);
          Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat_2, face_vec_2, cell_facet_2);

          // compute orientation trafos
          Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat_1, ori_vec_1, facet_ori_1);
          Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat_2, ori_vec_2, facet_ori_2);

          // prepare trafo evaluators
          trafo_facet_eval.prepare(face);
          trafo_eval_1.prepare(cell_1);
          trafo_eval_2.prepare(cell_2);

          // prepare space evaluators
          space_eval_1.prepare(trafo_eval_1);
          space_eval_2.prepare(trafo_eval_2);

          // format local matrix
          loc_mat.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // get cubature point
            auto cub_pt = cubature_rule.get_point(k);

            // transform to local facets
            auto cub_cf_1 = (face_mat_1 * ((ori_mat_1 * cub_pt) + ori_vec_1)) + face_vec_1;
            auto cub_cf_2 = (face_mat_2 * ((ori_mat_2 * cub_pt) + ori_vec_2)) + face_vec_2;

            // compute trafo data
            trafo_facet_eval(trafo_facet_data, cub_pt);
            trafo_eval_1(trafo_data_1, cub_cf_1);
            trafo_eval_2(trafo_data_2, cub_cf_2);

            // compute basis function data
            space_eval_1(space_data_1, trafo_data_1);
            space_eval_2(space_data_2, trafo_data_2);

            // compute weight factor
            const DataType weight = trafo_facet_data.jac_det * cubature_rule.get_weight(k) *
              Math::pow(jacdet_scal * trafo_facet_data.jac_det, jacdet_expo);

            // compute jump gradients
            for(int i(0); i < num_local_dofs; ++i)
            {
              jump_grad[i].format();
              const int i_1 = common_map.loc_1(i);
              const int i_2 = common_map.loc_2(i);
              // note the different signs to compute the jump
              if(i_1 > -1) jump_grad[i] += space_data_1.phi[i_1].grad;
              if(i_2 > -1) jump_grad[i] -= space_data_2.phi[i_2].grad;
            }

            // assemble jump stabilisation operator
            for(int i(0); i < num_local_dofs; ++i)
            {
              for(int j(0); j < num_local_dofs; ++j)
              {
                loc_mat(i,j) += weight * Tiny::dot(jump_grad[i], jump_grad[j]);
              }
            }

            // continue with next cubature point
          }

          // finish evaluators
          space_eval_2.finish();
          space_eval_1.finish();
          trafo_eval_2.finish();
          trafo_eval_1.finish();
          trafo_facet_eval.finish();

          // incorporate local matrix
          scatter_axpy(loc_mat, common_map, common_map, gamma);

          // continue with next cell
        }

        // okay, that's it
      }

      /**
       * \brief Assembles the jump operator onto a matrix.
       *
       * This function assembles the jump operator:
       *   \f[J(\varphi,\psi) = \alpha \sum_E \int_E [\varphi]\cdot[\psi]\f]
       *
       * \attention
       * The matrix must have an extended stencil, which must have been assembled by calling
       * Assembly::SymbolicAssembler::assemble_matrix_ext1() !
       *
       * \param[inout] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] space
       * The finite element space to be used.
       *
       * \param[in] cubature_factory
       * The cubature for integration. Note that this is a cubature rule on the facets.
       *
       * \param[in] alpha
       * The scaling factor alpha for the jump operator.
       */
      template<
        typename Matrix_,
        typename Space_,
        typename CubatureFactory_>
      void assemble_jump_operator_matrix(
        Matrix_& matrix,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1)) const
      {
        // validate matrix dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");

        // matrix type
        typedef Matrix_ MatrixType;
        // test-space type
        typedef Space_ SpaceType;

        // assembly traits
        typedef AsmTraits1<
          typename MatrixType::DataType,
          SpaceType,
          TrafoTags::none,
          SpaceTags::value> AsmTraits;

        typedef typename AsmTraits::DataType DataType;

        // shape types
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

        // fetch the trafo
        const TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval_1(trafo), trafo_eval_2(trafo);

        // create a trafo facet evaluator
        typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;
        typedef typename TrafoFacetEvaluator::template ConfigTraits<TrafoTags::jac_det>::EvalDataType TrafoFacetEvalData;
        TrafoFacetEvaluator trafo_facet_eval(trafo);

        // create space evaluators
        typename AsmTraits::SpaceEvaluator space_eval_1(space), space_eval_2(space);

        // create dof-mappings
        typename AsmTraits::DofMapping dof_mapping_1(space), dof_mapping_2(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data_1, trafo_data_2;
        TrafoFacetEvalData trafo_facet_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data_1, space_data_2;

        // create cubature rule
        typename Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_axpy(matrix);

        // common DOF mapping
        static constexpr int max_common_dofs = 2 * AsmTraits::max_local_test_dofs;
        Intern::CommonDofMap<max_common_dofs> common_map;

        // jump gradients
        typename AsmTraits::SpaceEvalTraits::BasisValueType jump_value[max_common_dofs];

        // local matrix data
        Tiny::Matrix<DataType, max_common_dofs, max_common_dofs> loc_mat;

        // trafo matrices and vectors
        Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat_1, face_mat_2;
        Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat_1, ori_mat_2;
        Tiny::Vector<DataType, shape_dim> face_vec_1, face_vec_2;
        Tiny::Vector<DataType, facet_dim> ori_vec_1, ori_vec_2;

        face_mat_1.format();
        face_mat_2.format();
        ori_mat_1.format();
        ori_mat_2.format();
        face_vec_1.format();
        face_vec_2.format();
        ori_vec_1.format();
        ori_vec_2.format();

        // loop over all cells of the mesh
        for(std::size_t f(0); f < _facets.size(); /*++f*/)
        {
          // get facet index
          const Index face = _facets[f];
          const Index cell_1 = _cells[f];
          const int cell_facet_1 = _cell_facet[f];
          const int facet_ori_1 = _facet_ori[f];
          ++f;

          // check for next cell
          const bool inner = ((f < _facets.size()) && (face == _facets[f]));
          const Index cell_2 = (inner ? _cells[f] : cell_1);
          const int cell_facet_2 = (inner ? _cell_facet[f] : cell_facet_1);
          const int facet_ori_2 = (inner ? _facet_ori[f] : facet_ori_1);
          if(inner)
            ++f;

          // prepare dof mappings
          dof_mapping_1.prepare(cell_1);
          dof_mapping_2.prepare(cell_2);

          // build local dof map
          common_map.clear();
          for(int i(0); i < dof_mapping_1.get_num_local_dofs(); ++i)
            common_map.push_1(dof_mapping_1.get_index(i), i);
          if(inner)
          {
            for(int i(0); i < dof_mapping_2.get_num_local_dofs(); ++i)
              common_map.push_2(dof_mapping_2.get_index(i), i);
          }

          // finish dof mapping
          dof_mapping_2.finish();
          dof_mapping_1.finish();

          // get number of common local dofs
          const int num_local_dofs = common_map.get_num_local_dofs();

          // compute facet trafos
          Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat_1, face_vec_1, cell_facet_1);
          Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat_2, face_vec_2, cell_facet_2);

          // compute orientation trafos
          Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat_1, ori_vec_1, facet_ori_1);
          Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat_2, ori_vec_2, facet_ori_2);

          // prepare trafo evaluators
          trafo_facet_eval.prepare(face);
          trafo_eval_1.prepare(cell_1);
          trafo_eval_2.prepare(cell_2);

          // prepare space evaluators
          space_eval_1.prepare(trafo_eval_1);
          space_eval_2.prepare(trafo_eval_2);

          // format local matrix
          loc_mat.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // get cubature point
            auto cub_pt = cubature_rule.get_point(k);

            // transform to local facets
            auto cub_cf_1 = (face_mat_1 * ((ori_mat_1 * cub_pt) + ori_vec_1)) + face_vec_1;
            auto cub_cf_2 = (face_mat_2 * ((ori_mat_2 * cub_pt) + ori_vec_2)) + face_vec_2;

            // compute trafo data
            trafo_facet_eval(trafo_facet_data, cub_pt);
            trafo_eval_1(trafo_data_1, cub_cf_1);
            trafo_eval_2(trafo_data_2, cub_cf_2);

            // compute basis function data
            space_eval_1(space_data_1, trafo_data_1);
            space_eval_2(space_data_2, trafo_data_2);

            // compute weight factor
            const DataType weight = trafo_facet_data.jac_det * cubature_rule.get_weight(k);

            // compute jump gradients
            for(int i(0); i < num_local_dofs; ++i)
            {
              jump_value[i] = DataType(0);
              const int i_1 = common_map.loc_1(i);
              const int i_2 = common_map.loc_2(i);
              // note the different signs to compute the jump
              if(i_1 > -1) jump_value[i] += space_data_1.phi[i_1].value;
              if(i_2 > -1) jump_value[i] -= space_data_2.phi[i_2].value;
            }

            // assemble jump stabilisation operator
            for(int i(0); i < num_local_dofs; ++i)
            {
              for(int j(0); j < num_local_dofs; ++j)
              {
                loc_mat(i,j) += weight * jump_value[i] * jump_value[j];
              }
            }

            // continue with next cubature point
          }

          // finish evaluators
          space_eval_2.finish();
          space_eval_1.finish();
          trafo_eval_2.finish();
          trafo_eval_1.finish();
          trafo_facet_eval.finish();

          // incorporate local matrix
          scatter_axpy(loc_mat, common_map, common_map, alpha);

          // continue with next cell
        }

        // okay, that's it
      }

    protected:
      bool find_local_facet(Index face, Index cell, int& facet, int& ori)
      {
        typedef typename Shape::FaceTraits<ShapeType, facet_dim>::ShapeType FacetType;
        static constexpr int num_facets = Shape::FaceTraits<ShapeType, shape_dim-1>::count;
        static constexpr int num_vaf = Shape::FaceTraits<FacetType, 0>::count;

        typedef Geometry::Intern::IndexRepresentative<FacetType> FacetRepr;

        const auto& vert_at_elem = _trafo.get_mesh().template get_index_set<shape_dim, 0>();
        const auto& vert_at_face = _trafo.get_mesh().template get_index_set<facet_dim, 0>();
        //const auto& face_at_elem = _trafo.get_mesh().template get_index_set<shape_dim, facet_dim>();

        Index face_verts[num_vaf];
        FacetRepr::compute(face_verts, vert_at_face[face]);

        for(int li(0); li < num_facets; ++li)
        {
          Intern::CompIndexMap<decltype(vert_at_elem[0]), ShapeType, shape_dim-1> cim(vert_at_elem[cell], li);
          Index lf_verts[num_vaf];
          FacetRepr::compute(lf_verts, cim);
          bool match = true;
          for(int k(0); k < num_vaf; ++k)
          {
            if(lf_verts[k] != face_verts[k])
            {
              match = false;
              break;
            }
          }
          if(match)
          {
            ori = Geometry::Intern::CongruencySampler<FacetType>::compare(vert_at_face[face], cim);
            facet = li;
            return true;
          }
        }
        return false;
      }
    }; // class TraceAssembler<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_TRACE_ASSEMBLER_HPP
