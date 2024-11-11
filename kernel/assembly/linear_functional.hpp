// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/assembly/base.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Base class for Linear Functionals
     *
     * This class acts as a base class and interface documentation for linear functionals which are used by
     * the LinearFunctionalAssembler assembly class template.
     *
     * \author Peter Zajac
     */
    class LinearFunctional
    {
    public:
#ifdef DOXYGEN
      static constexpr TrafoTags trafo_config = ...;
      static constexpr SpaceTags test_config = ...;
#endif // DOXYGEN

      /**
       * \brief Linear Functional Evaluator class template
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
#ifdef DOXYGEN
        /// value type (either equal to DataType or Tiny::Vector/Matrix of DataType's)
        typedef ... ValueType;
#endif // DOXYGEN

      public:
        /**
         * \brief Prepares the evaluator for a given cell
         *
         * \param[in] trafo_eval
         * A \transient reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& DOXY(trafo_eval))
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

        /**
         * \brief Point initialization function
         *
         * This function is called to initialize the evaluator for a new evaluation point.
         *
         * \param[in] tau
         * The \transient transformation data in the current evaluation point. \see Trafo::EvalData
         */
        void set_point(const TrafoData& DOXY(tau))
        {
          // do nothing
        }

#ifdef DOXYGEN
        /**
         * \brief Evaluation function
         *
         * This function evaluates the linear functional for a given test function in a single cubature point.
         *
         * \param[in] psi
         * The \transient test function data in the current evaluation point. \see Space::EvalData
         *
         * \returns
         * The value of the linear functional.
         */
        ValueType eval(const TestBasisData& psi) const;
#endif // DOXYGEN
      }; // class LinearFunctional::Evaluator<...>
    }; // class LinearFunctional
  } // namespace Assembly
} // namespace FEAT
