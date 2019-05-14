// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_DOMAIN_ASSEMBLER_HELPERS_HPP
#define KERNEL_ASSEMBLY_DOMAIN_ASSEMBLER_HELPERS_HPP 1

// includes, FEAT
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/basic_assembly_jobs.hpp>
#include <kernel/assembly/function_integral_jobs.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Assembles a bilinear operator into a matrix with identical test- and trial-spaces.
     *
     * \param[inout] matrix
     * The matrix that is to be assembled.
     *
     * \param[in] bilinear_operator
     * A reference to the operator implementing the BilinearOperator interface to be assembled.
     *
     * \param[in] space
     * A reference to the finite-element to be used as the test- and trial-space.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Matrix_, typename BilOp_, typename Space_>
    void assemble_bilinear_operator_matrix_1(DomainAssembler<Trafo_>& dom_asm, Matrix_& matrix,
      const BilOp_& bilinear_operator, const Space_& space, const String& cubature,
      const typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix column count");
      XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix row count");

      BilinearOperatorMatrixAssemblyJob1<BilOp_, Matrix_, Space_> job(
        bilinear_operator, matrix, space, cubature, alpha);
      dom_asm.assemble(job);
    }

    /**
     * \brief Assembles a bilinear operator into a matrix with different test- and trial-spaces.
     *
     * \param[inout] matrix
     * The matrix that is to be assembled.
     *
     * \param[in] bilinear_operator
     * A reference to the operator implementing the BilinearOperator interface to be assembled.
     *
     * \param[in] test_space
     * A reference to the finite-element test-space to be used.
     *
     * \param[in] trial_space
     * A reference to the finite-element trial-space to be used.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Matrix_, typename BilOp_, typename TestSpace_, typename TrialSpace_>
    void assemble_bilinear_operator_matrix_2(DomainAssembler<Trafo_>& dom_asm,
      Matrix_& matrix, const BilOp_& bilinear_operator, const TestSpace_& test_space,
      const TrialSpace_& trial_space, const String& cubature,
      const typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == test_space.get_trafo(), "domain assembler and test space have different trafos");
      XASSERTM(dom_asm.get_trafo() == trial_space.get_trafo(), "domain assembler and trial space have different trafos");
      XASSERTM(matrix.columns() == trial_space.get_num_dofs(), "invalid matrix column count");
      XASSERTM(matrix.rows() == test_space.get_num_dofs(), "invalid matrix row count");

      BilinearOperatorMatrixAssemblyJob2<BilOp_, Matrix_, TestSpace_, TrialSpace_> job(
        bilinear_operator, matrix, test_space, trial_space, cubature, alpha);
      dom_asm.assemble(job);
    }

    /**
     * \brief Assembles a linear functional into a vector
     *
     * \param[inout] vector
     * The vector that is to be assembled.
     *
     * \param[in] linear_functional
     * A reference to the functional implementing the Assembly::LinearFunctional interface to be assembled.
     *
     * \param[in] space
     * A reference to the finite-element to be used as the test-space.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Vector_, typename LinFunc_, typename Space_>
    void assemble_linear_functional_vector(DomainAssembler<Trafo_>& dom_asm,
      Vector_& vector, const LinFunc_& linear_functional, const Space_& space,
      const String& cubature, const typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector length");

      LinearFunctionalAssemblyJob<LinFunc_, Vector_, Space_> job(
        linear_functional, vector, space, cubature, alpha);
      dom_asm.assemble(job);
    }

    /**
     * \brief Assembles a force function into a vector
     *
     * \param[inout] vector
     * The vector that is to be assembled.
     *
     * \param[in] function
     * A reference to the force function implementing the Analytic::Function interface to be assembled.
     *
     * \param[in] space
     * A reference to the finite-element to be used as the test-space.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Vector_, typename Function_, typename Space_>
    void assemble_force_function_vector(DomainAssembler<Trafo_>& dom_asm,
      Vector_& vector, const Function_& function, const Space_& space,
      const String& cubature, const typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector length");

      ForceFunctionalAssemblyJob<Function_, Vector_, Space_> job(
        function, vector, space, cubature, alpha);
      dom_asm.assemble(job);
    }

    /**
     * \brief Assembles the integral of an analytic function
     *
     * \tparam DataType_
     * The (scalar) datatype in which the assembly is to be performed. Must always be given.
     *
     * \param[in] function
     * A reference to the analytic function that is to be integrated.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \returns
     * The integrals of the input function stored in an Assembly::FunctionIntegralInfo instance,
     * whose exact type is determined by the Assembly::AnalyticFunctionIntegral helper class.
     */
    template<int max_der_, typename DataType_, typename Trafo_, typename Function_>
    typename AnalyticFunctionIntegral<DataType_, Function_>::Type integrate_analytic_function(
      DomainAssembler<Trafo_>& dom_asm, const Function_& function, const String& cubature)
    {
      AnalyticFunctionIntegralJob<DataType_, Function_, Trafo_, max_der_> job(function, dom_asm.get_trafo(), cubature);
      dom_asm.assemble(job);
      return job.result();
    }

    /**
     * \brief Assembles the integral of a discrete finite element function
     *
     * \param[in] vector
     * A reference to the coefficient vector of the finite element function that is to be integrated.
     *
     * \param[in] space
     * A reference to the finite element space that the coefficient vector belongs to.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \returns
     * The integrals of the input function stored in an Assembly::FunctionIntegralInfo instance,
     * whose exact type is determined by the Assembly::DiscreteFunctionIntegral helper class.
     */
    template<int max_der_, typename Vector_, typename Trafo_, typename Space_>
    typename DiscreteFunctionIntegral<Vector_, Space_>::Type integrate_discrete_function(
       DomainAssembler<Trafo_>& dom_asm, const Vector_& vector, const Space_& space, const String& cubature)
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid coefficient vector length");
      DiscreteFunctionIntegralJob<Vector_, Space_, max_der_> job(vector, space, cubature);
      dom_asm.assemble(job);
      return job.result();
    }

    /**
     * \brief Assembles the integral of an (analytic - discrete) error function
     *
     * \param[in] function
     * A reference to the analytic function that is to be integrated.
     *
     * \param[in] vector
     * A reference to the coefficient vector of the finite element function that is to be integrated.
     *
     * \param[in] space
     * A reference to the finite element space that the coefficient vector belongs to.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \returns
     * The integrals of the error function stored in an Assembly::FunctionIntegralInfo instance,
     * whose exact type is determined by the Assembly::DiscreteFunctionIntegral helper class.
     */
    template<int max_der_, typename Function_, typename Vector_, typename Trafo_, typename Space_>
    typename DiscreteFunctionIntegral<Vector_, Space_>::Type integrate_error_function(
       DomainAssembler<Trafo_>& dom_asm, const Function_& function,
       const Vector_& vector, const Space_& space, const String& cubature)
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid coefficient vector length");
      ErrorFunctionIntegralJob<Function_, Vector_, Space_, max_der_> job(function, vector, space, cubature);
      dom_asm.assemble(job);
      return job.result();
    }
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_DOMAIN_ASSEMBLER_HELPERS_HPP
