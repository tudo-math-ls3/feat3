// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/assembly/trace_assembler.hpp>

namespace Turb
{
  template<typename Trafo>
  class TurbTraceAssembler : public Assembly::TraceAssembler<Trafo>
  {
  protected:
    typedef Trafo TrafoType;
    typedef typename TrafoType::MeshType MeshType;
    typedef typename TrafoType::ShapeType ShapeType;

    static constexpr int shape_dim = ShapeType::dimension;
    static constexpr int facet_dim = shape_dim-1;

    // this functions evaluates u_t, see equation (21)
    template<typename DT_>
    DT_ eval_u_t(DT_ k_val, Tiny::Vector<DT_, shape_dim> loc_v)
    {
      DT_ value1 = pow(_C_mu, 0.25) * sqrt(k_val);
      DT_ value2 = loc_v.norm_euclid() / _y_plus_star;
      return std::max(value1, value2);
    }

    // assembles the - u_t/y_*w_i term, see equation (23)
    template<typename DT_>
    DT_ eval_velocity_functional(DT_ u_t, DT_ test_eval)
    {
      return - (u_t / _y_plus_star) * test_eval;
    }

    // this function calculates local term in equation (27)
    template<typename DT_>
    DT_ eval_dissipation_functional(DT_ u_t, DT_ eps_val, DT_ test_val)
    {
      return (_kappa * u_t / _sigma_e) * eps_val * test_val;
    }

    using Assembly::TraceAssembler<Trafo>::_facets;
    using Assembly::TraceAssembler<Trafo>::_cells;
    using Assembly::TraceAssembler<Trafo>::_facet_mask;
    using Assembly::TraceAssembler<Trafo>::_facet_ori;
    using Assembly::TraceAssembler<Trafo>::_cell_facet;

  public:
    const TrafoType& _trafo;
    FEAT::Real _wall_distance;
    FEAT::Real _C_mu;
    FEAT::Real _y_plus_star;
    FEAT::Real _sigma_e;
    FEAT::Real _kappa;

    TurbTraceAssembler(
      const TrafoType& trafo,
      FEAT::Real wall_distance,
      FEAT::Real C_mu = 0.09,
      FEAT::Real y_plus_star = 11.06,
      FEAT::Real sigma_e = 1.3,
      FEAT::Real kappa = 0.41
    )
    : Assembly::TraceAssembler<TrafoType>(trafo),
      _trafo(trafo),
      _wall_distance(wall_distance),
      _C_mu(C_mu),
      _y_plus_star(y_plus_star),
      _sigma_e(sigma_e),
      _kappa(kappa)
    {}

    /**
      * \brief Assembles the surface integral for velocity
      *
      * \param[in,out] vec_out_velo
      * A \transient reference to the vector that represents the function to be integrated (i.e. a dual vector)
      *
      * \param[in] primal_turb_kinetic
      * A \transient reference to the primal (and synced 1) kinetic vector k
      *
      * \param[in] primal_dissipation_rate
      * A \transient reference to the primal (and synced 1) dissipation rate e
      *
      * \param[in] primal_conv_vec
      * A \transient reference to the primal (and synced 1) convection vector
      *
      * \param[in] space_turb
      * A \transient reference to the finite element space of the transportation terms
      *
      * \param[in] space_velo
      * A \transient reference to the finite element space of the convection
      *
      * \param[in] cubature_nane
      * The name of the cubature rule to be used for integration
      *
      * \param[in] alpha
      * The scaling factor for the linear functional.
      *
      * \returns The surface integral of the discrete function
      */
    template<
      typename Vector_,
      typename ConvVector_,
      typename Space_,
      typename ConvSpace_>
    void assemble_functional_vector_velocity(
      ConvVector_&  vec_out_velo,
      const Vector_& primal_turb_kinetic,
      const Vector_& primal_dissipation_rate,
      const ConvVector_& primal_conv_vec,
      const Space_& space_turb,
      const ConvSpace_& space_velo,
      const String& cubature_name,
      typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      Cubature::DynamicFactory cubature_factory(cubature_name);
      return assemble_functional_vector_velocity(vec_out_velo, primal_turb_kinetic, primal_dissipation_rate, primal_conv_vec, space_turb,
              space_velo, cubature_factory, alpha);
    }

    /**
      * \brief Assembles the surface integral for velocity
      *
      * \param[in,out] vec_out_velo
      * A \transient reference to the vector that represents the function to be integrated (i.e. a dual vector)
      *
      * \param[in] primal_turb_kinetic
      * A \transient reference to the primal (and synced 1) kinetic vector k
      *
      * \param[in] primal_dissipation_rate
      * A \transient reference to the primal (and synced 1) dissipation rate e
      *
      * \param[in] primal_conv_vec
      * A \transient reference to the primal (and synced 1) convection vector
      *
      * \param[in] space_turb
      * A \transient reference to the finite element space of the transportation terms
      *
      * \param[in] space_velo
      * A \transient reference to the finite element space of the convection
      *
      * \param[in] cubature_factory
      * A \transient reference to the cubature factory to be used for integration.
      *
      * \param[in] alpha
      * The scaling factor for the linear functional.
      *
      * \returns The surface integral of the discrete function
      */
    template<
      typename Vector_,
      typename ConvVector_,
      typename CubatureFactory_,
      typename Space_,
      typename ConvSpace_>
    void assemble_functional_vector_velocity(
      ConvVector_&  vec_out_velo,
      const Vector_& primal_turb_kinetic,
      const Vector_& primal_dissipation_rate,
      const ConvVector_& primal_conv_vec,
      const Space_& space_turb,
      const ConvSpace_& space_conv,
      const CubatureFactory_& cubature_factory,
      typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      // validate vector dimensions
      XASSERTM(vec_out_velo.size() == space_conv.get_num_dofs(), "invalid vector size");
      XASSERTM(primal_turb_kinetic.size() == space_turb.get_num_dofs(), "invalid vector size");
      XASSERTM(primal_dissipation_rate.size() == space_turb.get_num_dofs(), "invalid vector size");
      XASSERTM(primal_conv_vec.size() == space_conv.get_num_dofs(), "invalid vector size");

      // vector type
      typedef Vector_ VectorType;
      typedef ConvVector_ ConvVectorType;
      // space type
      typedef Space_ SpaceType;
      typedef ConvSpace_ ConvSpaceType;

      // assembly traits
      typedef Assembly::AsmTraits1<
        typename VectorType::DataType,
        SpaceType,
        TrafoTags::jac_det,
        SpaceTags::value
        > TurbAsmTraits;

      // assembly traits
      typedef Assembly::AsmTraits1<
        typename ConvVectorType::DataType,
        ConvSpaceType,
        TrafoTags::jac_det,
        SpaceTags::value
        > ConvAsmTraits;

      //typedef typename ConvAsmTraits::DataType DataType;

      // shape types
      typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

      // fetch the trafo
      const TrafoType& trafo = space_turb.get_trafo();

      // create a trafo evaluator
      typename TurbAsmTraits::TrafoEvaluator trafo_eval(trafo);

      // create a trafo facet evaluator
      typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;
      TrafoFacetEvaluator trafo_facet_eval(trafo);

      typedef typename TrafoFacetEvaluator::template ConfigTraits<TrafoTags::jac_det>::EvalDataType TrafoFacetEvalData;

      // create a space evaluator and evaluation data
      typename TurbAsmTraits::TestEvaluator test_eval_turb(space_turb);
      typename ConvAsmTraits::TestEvaluator test_eval_conv(space_conv);

      // create a dof-mapping
      typename TurbAsmTraits::DofMapping dof_mapping_turb(space_turb);
      typename ConvAsmTraits::DofMapping dof_mapping_conv(space_conv);

      // create trafo evaluation data
      typename TurbAsmTraits::TrafoEvalData trafo_data_turb;
      TrafoFacetEvalData trafo_facet_data;

      // create space evaluation data
      typename TurbAsmTraits::TestEvalData test_data_turb;
      typename ConvAsmTraits::TestEvalData test_data_conv;

      // the value type of our out vector is our value type
      typedef typename VectorType::ValueType ValueType;
      typedef typename ConvVectorType::ValueType ConvVectorValue;

      static constexpr int max_local_dofs_turb = TurbAsmTraits::max_local_test_dofs;
      static constexpr int max_local_dofs_conv = ConvAsmTraits::max_local_test_dofs;

      // define our local vector types
      typedef Tiny::Vector<ValueType, max_local_dofs_turb> LocalTurbVectorType;
      typedef Tiny::Vector<ConvVectorValue, max_local_dofs_conv> LocalConvVectorType;

      // create local vector data for gathering and scattering
      LocalConvVectorType loc_vec_out;
      LocalTurbVectorType loc_vec_kinetic_dofs;
      LocalTurbVectorType loc_vec_dissipation_dofs;
      LocalConvVectorType loc_conv_dofs;

      // create cubature rule
      typename Assembly::Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

      // create vector gather for our inputs
      typename VectorType::GatherAxpy gather_kinetic_axpy(primal_turb_kinetic);
      typename VectorType::GatherAxpy gather_dissipation_axpy(primal_dissipation_rate);
      typename ConvVectorType::GatherAxpy gather_conv_axpy(primal_conv_vec);

      // create scatter axpy for output
      typename ConvVectorType::ScatterAxpy scatter_conv(vec_out_velo);


      // trafo matrices and vectors
      Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat;
      Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat;
      Tiny::Vector<DataType, shape_dim> face_vec;
      Tiny::Vector<DataType, facet_dim> ori_vec;
      ConvVectorValue loc_v;

      face_mat.format();
      ori_mat.format();
      face_vec.format();
      ori_vec.format();
      loc_vec_out.format();
      loc_vec_kinetic_dofs.format();
      loc_vec_dissipation_dofs.format();
      loc_conv_dofs.format();
      loc_v.format();

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

        // compute orientation of actual cell facet
        const int cell_facet_ori = Geometry::Intern::CongruencySampler<FacetType>::orientation(_facet_ori[f])
          * Shape::ReferenceCell<ShapeType>::facet_orientation(_cell_facet[f]);

        // prepare trafo evaluators
        trafo_facet_eval.prepare(face);
        trafo_eval.prepare(cell);

        // prepare test-space evaluator
        test_eval_turb.prepare(trafo_eval);
        test_eval_conv.prepare(trafo_eval);

        // fetch number of local dofs
        const int num_loc_dofs_turb = test_eval_turb.get_num_local_dofs();
        const int num_loc_dofs_conv = test_eval_conv.get_num_local_dofs();

        // gather our local transport vectors
        dof_mapping_turb.prepare(cell);
        loc_vec_dissipation_dofs.format();
        gather_dissipation_axpy(loc_vec_dissipation_dofs, dof_mapping_turb);

        loc_vec_kinetic_dofs.format();
        gather_kinetic_axpy(loc_vec_kinetic_dofs, dof_mapping_turb);

        dof_mapping_turb.finish();

        //gather local conv dofs for our cell
        dof_mapping_conv.prepare(cell);

        loc_conv_dofs.format();
        gather_conv_axpy(loc_conv_dofs, dof_mapping_conv);

        // format local vector
        loc_vec_out.format();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // get cubature point
          auto cub_pt = cubature_rule.get_point(k);

          // transform to local facet
          auto cub_cf = (face_mat * ((ori_mat * cub_pt) + ori_vec)) + face_vec;

          // compute trafo data
          trafo_facet_eval(trafo_facet_data, cub_pt);
          trafo_eval(trafo_data_turb, cub_cf);

          // compute normal vector
          trafo_data_turb.normal = Tiny::orthogonal(trafo_facet_data.jac_mat).normalize();
          if(cell_facet_ori < 0)
            trafo_data_turb.normal.negate();

          // compute test basis function data
          test_eval_turb(test_data_turb, trafo_data_turb);
          test_eval_conv(test_data_conv, trafo_data_turb);

          // gather local velocity value
          loc_v.format();
          for(int i = 0; i < num_loc_dofs_conv; ++i)
          {
            loc_v.axpy(test_data_conv.phi[i].value, loc_conv_dofs[i]);
          }

          // gather local epsilon and kappa
          DataType eps_val = DataType(0);
          DataType k_val = DataType(0);
          for(int i = 0; i < num_loc_dofs_turb; ++i)
          {
            eps_val += test_data_turb.phi[i].value * loc_vec_dissipation_dofs[i];
            k_val += test_data_turb.phi[i].value * loc_vec_kinetic_dofs[i];
          }

          // calculate all relevant values
          DataType u_t = eval_u_t(k_val, loc_v);


          // test function loop
          for(int i(0); i < num_loc_dofs_conv; ++i)
          {
            const DataType value = trafo_facet_data.jac_det * cubature_rule.get_weight(k) * eval_velocity_functional(u_t, test_data_conv.phi[i].value);
            // evaluate functional and integrate
            loc_vec_out[i].axpy(value, loc_v);
            // continue with next trial function
          }
          // continue with next test function
        }

        // finish evaluators
        test_eval_conv.finish();
        test_eval_turb.finish();
        trafo_eval.finish();
        trafo_facet_eval.finish();

        // incorporate local matrix
        scatter_conv(loc_vec_out, dof_mapping_conv, alpha);

        // finish dof-mapping
        dof_mapping_conv.finish();

        // continue with next cell
      }
    }
    /**
      * \brief Assembles the surface integral for dissipation
      *
      * \param[in,out] vec_out_dissipation
      * A \transient reference to the vector that represents the function to be integrated (i.e. a dual vector)
      *
      * \param[in] primal_turb_kinetic
      * A \transient reference to the primal (and synced 1) kinetic vector k
      *
      * \param[in] primal_dissipation_rate
      * A \transient reference to the primal (and synced 1) dissipation rate e
      *
      * \param[in] primal_conv_vec
      * A \transient reference to the primal (and synced 1) convection vector
      *
      * \param[in] space_turb
      * A \transient reference to the finite element space of the transportation terms
      *
      * \param[in] space_velo
      * A \transient reference to the finite element space of the convection
      *
      * \param[in] cubature_nane
      * The name of the cubature rule to be used for integration
      *
      * \param[in] alpha
      * The scaling factor for the linear functional.
      *
      * \returns The surface integral of the discrete function
      */
    template<
      typename Vector_,
      typename ConvVector_,
      typename Space_,
      typename ConvSpace_>
    void assemble_functional_vector_dissipation(
      Vector_&  vec_out_dissipation,
      const Vector_& primal_turb_kinetic,
      const Vector_& primal_dissipation_rate,
      const ConvVector_& primal_conv_vec,
      const Space_& space_turb,
      const ConvSpace_& space_velo,
      const String& cubature_name,
      typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      Cubature::DynamicFactory cubature_factory(cubature_name);
      return assemble_functional_vector_dissipation(vec_out_dissipation, primal_turb_kinetic, primal_dissipation_rate, primal_conv_vec, space_turb,
              space_velo, cubature_factory, alpha);
    }

    /**
      * \brief Assembles the surface integral for dissipation
      *
      * \param[in,out] vec_out_dissipation
      * A \transient reference to the vector that represents the function to be integrated (i.e. a dual vector)
      *
      * \param[in] primal_turb_kinetic
      * A \transient reference to the primal (and synced 1) kinetic vector k
      *
      * \param[in] primal_dissipation_rate
      * A \transient reference to the primal (and synced 1) dissipation rate e
      *
      * \param[in] primal_conv_vec
      * A \transient reference to the primal (and synced 1) convection vector
      *
      * \param[in] space_turb
      * A \transient reference to the finite element space of the transportation terms
      *
      * \param[in] space_velo
      * A \transient reference to the finite element space of the convection
      *
      * \param[in] cubature_factory
      * A \transient reference to the cubature factory to be used for integration.
      *
      * \param[in] alpha
      * The scaling factor for the linear functional.
      *
      * \returns The surface integral of the discrete function
      */
    template<
      typename Vector_,
      typename ConvVector_,
      typename CubatureFactory_,
      typename Space_,
      typename ConvSpace_>
    void assemble_functional_vector_dissipation(
      Vector_&  vec_out_dissipation,
      const Vector_& primal_turb_kinetic,
      const Vector_& primal_dissipation_rate,
      const ConvVector_& primal_conv_vec,
      const Space_& space_turb,
      const ConvSpace_& space_conv,
      const CubatureFactory_& cubature_factory,
      typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      // validate vector dimensions
      XASSERTM(vec_out_dissipation.size() == space_turb.get_num_dofs(), "invalid vector size");
      XASSERTM(primal_turb_kinetic.size() == space_turb.get_num_dofs(), "invalid vector size");
      XASSERTM(primal_dissipation_rate.size() == space_turb.get_num_dofs(), "invalid vector size");
      XASSERTM(primal_conv_vec.size() == space_conv.get_num_dofs(), "invalid vector size");

      // vector type
      typedef Vector_ VectorType;
      typedef ConvVector_ ConvVectorType;
      // space type
      typedef Space_ SpaceType;
      typedef ConvSpace_ ConvSpaceType;

      // assembly traits
      typedef Assembly::AsmTraits1<
        typename VectorType::DataType,
        SpaceType,
        TrafoTags::jac_det,
        SpaceTags::value
        > TurbAsmTraits;

      // assembly traits
      typedef Assembly::AsmTraits1<
        typename ConvVectorType::DataType,
        ConvSpaceType,
        TrafoTags::jac_det,
        SpaceTags::value
        > ConvAsmTraits;

      //typedef typename TurbAsmTraits::DataType DataType;

      // shape types
      typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

      // fetch the trafo
      const TrafoType& trafo = space_turb.get_trafo();

      // create a trafo evaluator
      typename TurbAsmTraits::TrafoEvaluator trafo_eval(trafo);

      // create a trafo facet evaluator
      typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;
      TrafoFacetEvaluator trafo_facet_eval(trafo);

      typedef typename TrafoFacetEvaluator::template ConfigTraits<TrafoTags::jac_det>::EvalDataType TrafoFacetEvalData;

      // create a space evaluator and evaluation data
      typename TurbAsmTraits::TestEvaluator test_eval_turb(space_turb);
      typename ConvAsmTraits::TestEvaluator test_eval_conv(space_conv);

      // create a dof-mapping
      typename TurbAsmTraits::DofMapping dof_mapping_turb(space_turb);
      typename ConvAsmTraits::DofMapping dof_mapping_conv(space_conv);

      // create trafo evaluation data
      typename TurbAsmTraits::TrafoEvalData trafo_data_turb;
      TrafoFacetEvalData trafo_facet_data;

      // create space evaluation data
      typename TurbAsmTraits::TestEvalData test_data_turb;
      typename ConvAsmTraits::TestEvalData test_data_conv;

      // the value type of our out vector is our value type
      typedef typename VectorType::ValueType ValueType;
      typedef typename ConvVectorType::ValueType ConvVectorValue;

      static constexpr int max_local_dofs_turb = TurbAsmTraits::max_local_test_dofs;
      static constexpr int max_local_dofs_conv = ConvAsmTraits::max_local_test_dofs;

      // define our local vector types
      typedef Tiny::Vector<ValueType, max_local_dofs_turb> LocalTurbVectorType;
      typedef Tiny::Vector<ConvVectorValue, max_local_dofs_conv> LocalConvVectorType;

      // create local vector data for gathering and scattering
      LocalTurbVectorType loc_vec_out;
      LocalTurbVectorType loc_vec_kinetic_dofs;
      LocalTurbVectorType loc_vec_dissipation_dofs;
      LocalConvVectorType loc_conv_dofs;

      // create cubature rule
      typename Assembly::Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

      // create vector gather for our inputs
      typename VectorType::GatherAxpy gather_kinetic_axpy(primal_turb_kinetic);
      typename VectorType::GatherAxpy gather_dissipation_axpy(primal_dissipation_rate);
      typename ConvVectorType::GatherAxpy gather_conv_axpy(primal_conv_vec);

      // create scatter axpy for output
      typename VectorType::ScatterAxpy scatter_dissipation(vec_out_dissipation);

      // trafo matrices and vectors
      Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat;
      Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat;
      Tiny::Vector<DataType, shape_dim> face_vec;
      Tiny::Vector<DataType, facet_dim> ori_vec;
      ConvVectorValue loc_v;

      face_mat.format();
      ori_mat.format();
      face_vec.format();
      ori_vec.format();
      loc_vec_out.format();
      loc_vec_kinetic_dofs.format();
      loc_vec_dissipation_dofs.format();
      loc_conv_dofs.format();
      loc_v.format();

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

        // compute orientation of actual cell facet
        const int cell_facet_ori = Geometry::Intern::CongruencySampler<FacetType>::orientation(_facet_ori[f])
          * Shape::ReferenceCell<ShapeType>::facet_orientation(_cell_facet[f]);

        // prepare trafo evaluators
        trafo_facet_eval.prepare(face);
        trafo_eval.prepare(cell);

        // prepare test-space evaluator
        test_eval_turb.prepare(trafo_eval);
        test_eval_conv.prepare(trafo_eval);

        // fetch number of local dofs
        const int num_loc_dofs_turb = test_eval_turb.get_num_local_dofs();
        const int num_loc_dofs_conv = test_eval_conv.get_num_local_dofs();


        //gather local conv dofs for our cell
        dof_mapping_conv.prepare(cell);

        loc_conv_dofs.format();
        gather_conv_axpy(loc_conv_dofs, dof_mapping_conv);

        dof_mapping_conv.finish();

        // gather our local transport vectors
        dof_mapping_turb.prepare(cell);
        loc_vec_dissipation_dofs.format();
        gather_dissipation_axpy(loc_vec_dissipation_dofs, dof_mapping_turb);

        loc_vec_kinetic_dofs.format();
        gather_kinetic_axpy(loc_vec_kinetic_dofs, dof_mapping_turb);

        // format local vector
        loc_vec_out.format();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // get cubature point
          auto cub_pt = cubature_rule.get_point(k);

          // transform to local facet
          auto cub_cf = (face_mat * ((ori_mat * cub_pt) + ori_vec)) + face_vec;

          // compute trafo data
          trafo_facet_eval(trafo_facet_data, cub_pt);
          trafo_eval(trafo_data_turb, cub_cf);

          // compute normal vector
          trafo_data_turb.normal = Tiny::orthogonal(trafo_facet_data.jac_mat).normalize();
          if(cell_facet_ori < 0)
            trafo_data_turb.normal.negate();

          // compute test basis function data
          test_eval_turb(test_data_turb, trafo_data_turb);
          test_eval_conv(test_data_conv, trafo_data_turb);

          // gather local velocity value
          loc_v.format();
          for(int i = 0; i < num_loc_dofs_conv; ++i)
          {
            loc_v.axpy(test_data_conv.phi[i].value, loc_conv_dofs[i]);
          }

          // gather local epsilon and kappa
          DataType eps_val = DataType(0);
          DataType k_val = DataType(0);
          for(int i = 0; i < num_loc_dofs_turb; ++i)
          {
            eps_val += test_data_turb.phi[i].value * loc_vec_dissipation_dofs[i];
            k_val += test_data_turb.phi[i].value * loc_vec_kinetic_dofs[i];
          }

          // calculate all relevant values
          DataType u_t = eval_u_t(k_val, loc_v);


          // test function loop
          for(int i(0); i < num_loc_dofs_turb; ++i)
          {
            // evaluate functional and integrate
            Tiny::axpy(loc_vec_out(i), eval_dissipation_functional(u_t, eps_val, test_data_turb.phi[i].value),
              trafo_facet_data.jac_det * cubature_rule.get_weight(k));
            // continue with next trial function
          }
          // continue with next test function
        }

        // finish evaluators
        test_eval_conv.finish();
        test_eval_turb.finish();
        trafo_eval.finish();
        trafo_facet_eval.finish();

        // incorporate local matrix
        scatter_dissipation(loc_vec_out, dof_mapping_turb, alpha);

        // finish dof-mapping
        dof_mapping_turb.finish();

        // continue with next cell
      }
    }
  }; // class TurbTraceAssembler
} // namespace Turb
