#pragma once
#ifndef KERNEL_ASSEMBLY_ASM_TRAITS_HPP
#define KERNEL_ASSEMBLY_ASM_TRAITS_HPP 1

// includes, FEAT
#include <kernel/assembly/base.hpp>
#include <kernel/cubature/dynamic_factory.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Space_, typename DataType_>
      struct EvalPolicyFetcher
      {
        typedef typename Space_::TrafoType TrafoType;
        typedef typename TrafoType::ShapeType ShapeType;
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalPolicy EvalPolicy;
      };

      template<typename TrafoEvaluator_>
      class CubatureTraits
      {
      public:
        /// evaluation policy
        typedef typename TrafoEvaluator_::EvalPolicy EvalPolicy;

        /// data type
        typedef typename EvalPolicy::DataType DataType;

        /// shape type
        typedef typename EvalPolicy::ShapeType ShapeType;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;

        /// cubature rule type
        typedef Cubature::Rule<ShapeType, DataType, DataType, DomainPointType> RuleType;
      };
    }
    /// \endcond

    /**
      * \brief Common single-space assembly traits class template
      *
      * This class template takes care of defining the necessary classes for assembly with one single
      * finite element space.
      *
      * This class can e.g. be used as a base class for
      * - assembly of a linear functional, e.g. a right-hand-side vector
      * - assembly of a bilinear operator with identical test- and trial-spaces
      * - assembly of a trilinear operator with identical test-, trial- and multiplier-spaces
      * - post-processing of a primal vector, e.g. L2- or H1-error calculation
      *
      * \tparam DataType_
      * The data type that is used to be for the assembly.
      *
      * \tparam Space_
      * The finite element space that is to be used as a (test- and trial-) space.
      *
      * \tparam TrafoConfig_
      * A trafo config class defining additional trafo requirements, e.g. from a (bi)linear functor.
      *
      * \tparam SpaceConfig_
      * A space config class defining additional space requirements, e.g. from a (bi)linear functor.
      *
      * \author Peter Zajac
      */
    template<
      typename DataType_,
      typename Space_,
      TrafoTags trafo_config_,
      SpaceTags space_config_>
    class AsmTraits1 :
      public Intern::EvalPolicyFetcher<Space_, DataType_>::EvalPolicy
    {
    public:
      /// data type
      typedef DataType_ DataType;
      /// space type
      typedef Space_ SpaceType;
      /// test-space type
      typedef SpaceType TestSpaceType;
      /// trial-space type
      typedef SpaceType TrialSpaceType;
      /// multiplier space type
      typedef SpaceType MultSpaceType;

      /// trafo type
      typedef typename SpaceType::TrafoType TrafoType;
      /// shape type
      typedef typename TrafoType::ShapeType ShapeType;
      /// mesh type
      typedef typename TrafoType::MeshType MeshType;

      /// trafo evaluator type
      typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;

      /// trafo cell iterator type
      typedef typename TrafoEvaluator::CellIterator CellIterator;

      /// space evaluator types
      typedef typename SpaceType::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;
      typedef SpaceEvaluator TestEvaluator;
      typedef SpaceEvaluator TrialEvaluator;
      typedef SpaceEvaluator MultEvaluator;

      /// space evaluator traits
      typedef typename SpaceEvaluator::SpaceEvalTraits SpaceEvalTraits;
      typedef SpaceEvalTraits TestEvalTraits;
      typedef SpaceEvalTraits TrialEvalTraits;
      typedef SpaceEvalTraits MultEvalTraits;

      // define test- and trial-space configs
      static constexpr SpaceTags space_config = space_config_;
      static constexpr SpaceTags test_config = space_config_;
      static constexpr SpaceTags trial_config = space_config_;
      static constexpr SpaceTags mult_config = space_config_;

      // now fetch the trafo config from the space
      typedef typename SpaceEvaluator::template ConfigTraits<space_config> SpaceConfigTraits;
      static constexpr TrafoTags space_trafo_config = SpaceConfigTraits::trafo_config;
      static constexpr TrafoTags test_trafo_config = SpaceConfigTraits::trafo_config;
      static constexpr TrafoTags trial_trafo_config = SpaceConfigTraits::trafo_config;
      static constexpr TrafoTags mult_trafo_config = SpaceConfigTraits::trafo_config;

      /// assembly trafo config: derive from user-defined trafo config
      static constexpr TrafoTags trafo_config = trafo_config_ | space_trafo_config | TrafoTags::jac_det;

      /// trafo evaluation data type
      typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;
      typedef TrafoEvalData TrafoData;

      /// space evaluation data types
      typedef typename SpaceEvaluator::template ConfigTraits<space_config>::EvalDataType SpaceEvalData;
      typedef SpaceEvalData TestEvalData;
      typedef SpaceEvalData TrialEvalData;
      typedef SpaceEvalData MultEvalData;

      /// basis function data types
      typedef typename SpaceEvalData::BasisDataType SpaceBasisData;
      typedef SpaceBasisData BasisData;
      typedef SpaceBasisData TestBasisData;
      typedef SpaceBasisData TrialBasisData;
      typedef SpaceBasisData MultBasisData;

      /// dof-mapping types
      typedef typename SpaceType::DofMappingType DofMapping;
      typedef DofMapping TestDofMapping;
      typedef DofMapping TrialDofMapping;
      typedef DofMapping MultDofMapping;

      /// trafo domain dimension
      static constexpr int domain_dim = TrafoEvaluator::domain_dim;
      /// trafo image dimension
      static constexpr int image_dim = TrafoEvaluator::image_dim;

      /// maximum local dofs
      static constexpr int max_local_test_dofs  = TestEvaluator::max_local_dofs;
      static constexpr int max_local_trial_dofs = TrialEvaluator::max_local_dofs;
      static constexpr int max_local_mult_dofs  = MultEvaluator::max_local_dofs;

      /// local vector type
      typedef Tiny::Vector<DataType, SpaceEvaluator::max_local_dofs> LocalVectorType;
      typedef LocalVectorType LocalTestVectorType;
      typedef LocalVectorType LocalTrialVectorType;
      typedef LocalVectorType LocalMultVectorType;

      /// local matrix type
      typedef Tiny::Matrix<DataType, SpaceEvaluator::max_local_dofs, SpaceEvaluator::max_local_dofs> LocalMatrixType;

      /// cubature rule type
      typedef typename Intern::CubatureTraits<TrafoEvaluator>::RuleType CubatureRuleType;
    }; // class AsmTraits1

    /**
      * \brief Common test-/trial-space assembly traits class template
      *
      * This class template takes care of defining the necessary classes for assembly with a combination of different
      * test- and trial-spaces using the same transformation.
      *
      * This class can e.g. be used as a base class for
      * - assembly of a bilinear operator with different test- and trial-spaces
      * - assembly of a trilinear operator with different test- and trial-spaces but identical trial- and
      *   multiplier-spaces
      *
      * \tparam DataType_
      * The data type that is used to be for the assembly.
      *
      * \tparam TestSpace_
      * The finite element space that is to be used as the test-space.
      *
      * \tparam TrialSpace_
      * The finite element space that is to be used as the trial-space. Must be defined on the same  trafo object as
      * \p TestSpace_.
      *
      * \tparam TrafoConfig_
      * A trafo config class defining additional trafo requirements, e.g. from a (bi)linear functor.
      *
      * \tparam TestConfig_, TrialConfig_
      * Two space config classes defining additional test- and trial-space requirements, e.g. from a (bi)linear functor.
      *
      * \author Peter Zajac
      */
    template<
      typename DataType_,
      typename TestSpace_,
      typename TrialSpace,
      TrafoTags trafo_config_,
      SpaceTags test_config_,
      SpaceTags trial_config_>
    class AsmTraits2 :
      public Intern::EvalPolicyFetcher<TestSpace_, DataType_>::EvalPolicy
    {
    public:
      /// data type
      typedef DataType_ DataType;
      /// test-space type
      typedef TestSpace_ TestSpaceType;
      /// trial-space type
      typedef TrialSpace TrialSpaceType;
      /// mult-space type
      typedef TrialSpaceType MultSpaceType;

      /// trafo type
      typedef typename TestSpaceType::TrafoType TrafoType;
      /// shape type
      typedef typename TrafoType::ShapeType ShapeType;
      /// mesh type
      typedef typename TrafoType::MeshType MeshType;

      /// trafo evaluator type
      typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;

      /// trafo cell iterator type
      typedef typename TrafoEvaluator::CellIterator CellIterator;

      /// space evaluator types
      typedef typename TestSpaceType::template Evaluator<TrafoEvaluator>::Type TestEvaluator;
      typedef typename TrialSpaceType::template Evaluator<TrafoEvaluator>::Type TrialEvaluator;
      typedef TrialEvaluator MultEvaluator;

      /// space evaluator traits
      typedef typename TestEvaluator::SpaceEvalTraits TestEvalTraits;
      typedef typename TrialEvaluator::SpaceEvalTraits TrialEvalTraits;
      typedef TrialEvalTraits MultEvalTraits;

      // define test- and trial-space configs
      static constexpr SpaceTags test_config = test_config_;
      static constexpr SpaceTags trial_config = trial_config_;
      static constexpr SpaceTags mult_config = trial_config_;

      // now fetch the trafo configs from the spaces
      typedef typename TestEvaluator::template ConfigTraits<test_config> TestConfigTraits;
      typedef typename TrialEvaluator::template ConfigTraits<trial_config> TrialConfigTraits;
      static constexpr TrafoTags test_trafo_config = TestConfigTraits::trafo_config;
      static constexpr TrafoTags trial_trafo_config = TrialConfigTraits::trafo_config;
      static constexpr TrafoTags mult_trafo_config = TrialConfigTraits::trafo_config;

      /// trafo config: combine space and assembly trafo configs
      static constexpr TrafoTags trafo_config = trafo_config_ | test_trafo_config | trial_trafo_config | TrafoTags::jac_det;

      /// trafo evaluation data type
      typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;
      typedef TrafoEvalData TrafoData;

      /// space evaluation data types
      typedef typename TestEvaluator::template ConfigTraits<test_config>::EvalDataType TestEvalData;
      typedef typename TrialEvaluator::template ConfigTraits<trial_config>::EvalDataType TrialEvalData;
      typedef TrialEvalData MultEvalData;

      /// basis function data types
      typedef typename TestEvalData::BasisDataType TestBasisData;
      typedef typename TrialEvalData::BasisDataType TrialBasisData;
      typedef TrialBasisData MultBasisData;

      /// dof-mapping types
      typedef typename TestSpaceType::DofMappingType TestDofMapping;
      typedef typename TrialSpaceType::DofMappingType TrialDofMapping;
      typedef TrialDofMapping MultDofMapping;

      /// trafo domain dimension
      static constexpr int domain_dim = TrafoEvaluator::domain_dim;
      /// trafo image dimension
      static constexpr int image_dim = TrafoEvaluator::image_dim;

      /// maximum local dofs
      static constexpr int max_local_test_dofs  = TestEvaluator::max_local_dofs;
      static constexpr int max_local_trial_dofs = TrialEvaluator::max_local_dofs;
      static constexpr int max_local_mult_dofs  = MultEvaluator::max_local_dofs;

      /// local vector type
      typedef Tiny::Vector<DataType, TestEvaluator::max_local_dofs> LocalTestVectorType;
      typedef Tiny::Vector<DataType, TrialEvaluator::max_local_dofs> LocalTrialVectorType;
      typedef LocalTrialVectorType LocalMultVectorType;

      /// local matrix type
      typedef Tiny::Matrix<DataType, TestEvaluator::max_local_dofs, TrialEvaluator::max_local_dofs> LocalMatrixType;

      /// cubature rule type
      typedef typename Intern::CubatureTraits<TrafoEvaluator>::RuleType CubatureRuleType;
    }; // class AsmTraits2

    /**
      * \brief Common single-space block matrix assembly traits class template
      *
      * This class template takes care of defining the necessary classes for the assembly with a vector valued operator
      * with one single scalar finite element space. An example for this is the Du : Dv operator.
      *
      * This class can e.g. be used as a base class for
      * - assembly of a linear functional, e.g. a right-hand-side vector
      * - assembly of a bilinear operator with identical test- and trial-spaces
      * - assembly of a trilinear operator with identical test-, trial- and multiplier-spaces
      * - post-processing of a primal vector, e.g. L2- or H1-error calculation
      *
      * \tparam DataType_
      * The data type that is used to be for the assembly.
      *
      * \tparam Space_
      * The finite element space that is to be used as a (test- and trial-) space.
      *
      * \tparam BlockHeight_
      * Height of the blocks, i.e. the dimension of the space the operator maps to.
      *
      * \tparam BlockWidth_
      * Width of the blocks, i.e. the dimension of the space the operator maps from.
      *
      * \tparam TrafoConfig_
      * A trafo config class defining additional trafo requirements, e.g. from a (bi)linear functor.
      *
      * \tparam SpaceConfig_
      * A space config class defining additional space requirements, e.g. from a (bi)linear functor.
      *
      * \author Jordi Paul
      */
    template<
      typename DataType_,
      typename OperatorValueType_,
      typename Space_,
      TrafoTags trafo_config_,
      SpaceTags space_config_>
    class AsmTraits1Blocked :
      public Intern::EvalPolicyFetcher<Space_, DataType_>::EvalPolicy
    {
    public:
      /// data type
      typedef DataType_ DataType;
      /// type an operator returns
      typedef OperatorValueType_ OperatorValueType;
      static constexpr int BlockHeight = OperatorValueType::m;
      static constexpr int BlockWidth = OperatorValueType::n;
      /// type a functional returns
      typedef FEAT::Tiny::Vector<DataType,BlockWidth> FunctionalDataType;
      /// space type
      typedef Space_ SpaceType;
      /// test-space type
      typedef SpaceType TestSpaceType;
      /// trial-space type
      typedef SpaceType TrialSpaceType;
      /// multiplier space type
      typedef SpaceType MultSpaceType;

      /// trafo type
      typedef typename SpaceType::TrafoType TrafoType;
      /// shape type
      typedef typename TrafoType::ShapeType ShapeType;
      /// mesh type
      typedef typename TrafoType::MeshType MeshType;

      /// trafo evaluator type
      typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;

      /// trafo cell iterator type
      typedef typename TrafoEvaluator::CellIterator CellIterator;

      /// space evaluator types
      typedef typename SpaceType::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;
      typedef SpaceEvaluator TestEvaluator;
      typedef SpaceEvaluator TrialEvaluator;
      typedef SpaceEvaluator MultEvaluator;

      /// space evaluator traits
      typedef typename SpaceEvaluator::SpaceEvalTraits SpaceEvalTraits;
      typedef SpaceEvalTraits TestEvalTraits;
      typedef SpaceEvalTraits TrialEvalTraits;
      typedef SpaceEvalTraits MultEvalTraits;

      // define test- and trial-space configs
      static constexpr SpaceTags space_config = space_config_;
      static constexpr SpaceTags test_config = space_config_;
      static constexpr SpaceTags trial_config = space_config_;
      static constexpr SpaceTags mult_config = space_config_;

      // now fetch the trafo config from the space
      typedef typename SpaceEvaluator::template ConfigTraits<space_config> SpaceConfigTraits;
      static constexpr TrafoTags space_trafo_config = SpaceConfigTraits::trafo_config;
      static constexpr TrafoTags test_trafo_config = SpaceConfigTraits::trafo_config;
      static constexpr TrafoTags trial_trafo_config = SpaceConfigTraits::trafo_config;
      static constexpr TrafoTags mult_trafo_config = SpaceConfigTraits::trafo_config;

      /// assembly trafo config: derive from user-defined trafo config
      static constexpr TrafoTags trafo_config = trafo_config_ | space_trafo_config | TrafoTags::jac_det;

      /// trafo evaluation data type
      typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;
      typedef TrafoEvalData TrafoData;

      /// space evaluation data types
      typedef typename SpaceEvaluator::template ConfigTraits<space_config>::EvalDataType SpaceEvalData;
      typedef SpaceEvalData TestEvalData;
      typedef SpaceEvalData TrialEvalData;
      typedef SpaceEvalData MultEvalData;

      /// basis function data types
      typedef typename SpaceEvalData::BasisDataType SpaceBasisData;
      typedef SpaceBasisData BasisData;
      typedef SpaceBasisData TestBasisData;
      typedef SpaceBasisData TrialBasisData;
      typedef SpaceBasisData MultBasisData;

      /// dof-mapping types
      typedef typename SpaceType::DofMappingType DofMapping;
      typedef DofMapping TestDofMapping;
      typedef DofMapping TrialDofMapping;
      typedef DofMapping MultDofMapping;

      /// trafo domain dimension
      static constexpr int domain_dim = TrafoEvaluator::domain_dim;
      /// trafo image dimension
      static constexpr int image_dim = TrafoEvaluator::image_dim;

      /// maximum local dofs
      static constexpr int max_local_test_dofs  = TestEvaluator::max_local_dofs;
      static constexpr int max_local_trial_dofs = TrialEvaluator::max_local_dofs;
      static constexpr int max_local_mult_dofs  = MultEvaluator::max_local_dofs;

      /// local vector type
      typedef Tiny::Vector<DataType, SpaceEvaluator::max_local_dofs> FunctionalValueType;
      typedef Tiny::Vector<FunctionalValueType, SpaceEvaluator::max_local_dofs> LocalVectorType;
      typedef Tiny::Vector<FunctionalValueType, BlockHeight> LocalTestVectorType;
      typedef Tiny::Vector<FunctionalValueType, BlockWidth>  LocalTrialVectorType;
      typedef Tiny::Vector<FunctionalValueType, BlockWidth> LocalMultVectorType;

      /// local matrix type
      typedef Tiny::Matrix<OperatorValueType, SpaceEvaluator::max_local_dofs, SpaceEvaluator::max_local_dofs>
        LocalMatrixType;

      /// cubature rule type
      typedef typename Intern::CubatureTraits<TrafoEvaluator>::RuleType CubatureRuleType;
    }; // class AsmTraits1Blocked
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_ASM_TRAITS_HPP
