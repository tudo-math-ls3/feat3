#pragma once
#ifndef KERNEL_ASSEMBLY_BILINEAR_FUNCTOR_BASE_HPP
#define KERNEL_ASSEMBLY_BILINEAR_FUNCTOR_BASE_HPP 1

// includes, FEAST
#include <kernel/assembly/base.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Base class for Bilinear functors
     *
     * This class acts as a base class and interface documentation for bilinear operator functors which
     * are used by the BilinearOperator assembly class template.
     *
     * \author Peter Zajac
     */
    class BilinearFunctorBase
    {
    public:
      /// use 'base' trafo configuration
      typedef Trafo::ConfigBase TrafoConfig;
      /// use 'base' test space configuration
      typedef Space::ConfigBase TestConfig;
      /// use 'base' trial space configuration
      typedef Space::ConfigBase TrialConfig;

      /**
       * \brief Bilinear Functor Evaluator class template
       *
       * \tparam AsmTraits_
       * The assembly traits class.
       *
       * \author Peter Zajac
       */
      template<typename AsmTraits_>
      class Evaluator
      {
      public:
        /// trafo evaluator type
        typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;
        /// data type
        typedef typename AsmTraits_::DataType DataType;
        /// trafo data type
        typedef typename AsmTraits_::TrafoData TrafoData;
        /// test function data type
        typedef typename AsmTraits_::TestBasisData TestBasisData;
        /// trial function data type
        typedef typename AsmTraits_::TrialBasisData TrialBasisData;

      public:
        /**
         * \brief Prepares the evaluator for a given cell
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& /*trafo_eval*/)
        {
          // do nothing
        }

        /**
         * \brief Releases the evaluator from the current cell.
         */
        void finish()
        {
          // do nothing
        }

#ifdef DOXYGEN
        /**
         * \brief Evaluation operator
         *
         * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
         * a single point.
         *
         * \param[in] tau
         * The transformation data in the current evaluation point. \see TrafoEvalData
         *
         * \param[in] phi
         * The trial function data in the current evaluation point. \see SpaceEvalData
         *
         * \param[in] psi
         * The test function data in the current evaluation point. \see SpaceEvalData
         *
         * \returns
         * The value of the bilinear functor.
         */
        DataType operator()(const TrafoData& tau, const TrialBasisData& phi, const TestBasisData& psi) const;
#endif // DOXYGEN
      }; // class BilinearFunctorBase::Evaluator<...>
    }; // class BilinearFunctorBase
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_BILINEAR_FUNCTOR_BASE_HPP
