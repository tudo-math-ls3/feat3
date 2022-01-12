// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_MEAN_FILTER_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_MEAN_FILTER_ASSEMBLER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/analytic/common.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Mean Filter assembler class
     *
     * This class assembles an integral mean filter for a given finite element space.
     *
     * \author Peter Zajac
     */
    class MeanFilterAssembler
    {
    public:
      /**
       * \brief Assembles an integral mean filter
       *
       * \param[out] vec_prim, vec_dual
       * The \transient references to the primal and dual vectors for the mean filter.
       *
       * \param[in] space
       * The \transient finite element space for which the filter is to be assembled.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to use for integration
       */
      template<typename MemType_, typename DataType_, typename IndexType_, typename Space_>
      static void assemble(
        LAFEM::DenseVector<MemType_, DataType_, IndexType_>& vec_prim,
        LAFEM::DenseVector<MemType_, DataType_, IndexType_>& vec_dual,
        const Space_& space, const String& cubature_name)
      {
        Cubature::DynamicFactory cubature_factory(cubature_name);

        assemble(vec_prim, vec_dual, space, cubature_factory);
      }

      /**
       * \brief Assembles an integral mean filter
       *
       * \param[out] vec_prim, vec_dual
       * The \transient references to the primal and dual vectors for the mean filter.
       *
       * \param[in] space
       * The \transient finite element space for which the filter is to be assembled.
       *
       * \param[in] cubature_factory
       * A cubature factory for integration.
       */
      template<typename MemType_, typename DataType_, typename IndexType_, typename Space_>
      static void assemble(
        LAFEM::DenseVector<MemType_, DataType_, IndexType_>& vec_prim,
        LAFEM::DenseVector<MemType_, DataType_, IndexType_>& vec_dual,
        const Space_& space, const Cubature::DynamicFactory& cubature_factory)
      {
        // allocate primal and dual vectors
        LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> vec_v(space.get_num_dofs(), DataType_(0));
        LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> vec_w(space.get_num_dofs(), DataType_(0));

        // create a constant 1-function and its corresponding force functional
        Analytic::Common::ConstantFunction<Space_::world_dim, DataType_> one_func(DataType_(1));
        Assembly::Common::ForceFunctional<decltype(one_func)> one_force(one_func);

        // interpolate 1-function into vector v
        Assembly::Interpolator::project(vec_v, one_func, space);

        // assemble 1-function force into vector w
        Assembly::LinearFunctionalAssembler::assemble_vector(vec_w, one_force, space, cubature_factory);

        // convert mem types
        vec_prim.convert(vec_v);
        vec_dual.convert(vec_w);
      }

      /**
       * \brief Assembles an integral mean filter
       *
       * \param[out] filter
       * A \transient reference to the filter to be assembled. The filter is automatically allocated
       * by this function, so it does not need to be allocated beforehand.
       *
       * \param[in] space
       * A \transient reference to the finite element space for which the filter is to be assembled.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to use for integration
       *
       * \param[in] sol_mean
       * The desired integral mean of the solution vector. Defaults to 0.
       *
       * \warning This function does not work for global mean filters!
       */
      template<typename MemType_, typename DataType_, typename IndexType_, typename Space_>
      static void assemble(
        LAFEM::MeanFilter<MemType_, DataType_, IndexType_>& filter,
        const Space_& space, const String& cubature_name,
        const DataType_ sol_mean = DataType_(0))
      {
        Cubature::DynamicFactory cubature_factory(cubature_name);
        assemble(filter, space, cubature_factory, sol_mean);
      }

      /**
       * \brief Assembles an integral mean filter
       *
       * \param[out] filter
       * A \transient reference to the filter to be assembled. The filter is automatically allocated
       * by this function, so it does not need to be allocated beforehand.
       *
       * \param[in] space
       * A \transient reference to the finite element space for which the filter is to be assembled.
       *
       * \param[in] cubature_factory
       * A cubature factory for integration.
       *
       * \param[in] sol_mean
       * The desired integral mean of the solution vector. Defaults to 0.
       *
       * \warning This function does not work for global mean filters!
       */
      template<typename MemType_, typename DataType_, typename IndexType_, typename Space_>
      static void assemble(
        LAFEM::MeanFilter<MemType_, DataType_, IndexType_>& filter,
        const Space_& space, const Cubature::DynamicFactory& cubature_factory,
        const DataType_ sol_mean = DataType_(0))
      {
        // allocate primal and dual vectors
        LAFEM::DenseVector<MemType_, DataType_, IndexType_> vec_prim, vec_dual;

        // assemble vectors
        assemble(vec_prim, vec_dual, space, cubature_factory);

        // create the filter
        filter = LAFEM::MeanFilter<MemType_, DataType_, IndexType_>
          (std::move(vec_prim), std::move(vec_dual), sol_mean);
      }
    }; // class MeanFilterAssembler
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_MEAN_FILTER_ASSEMBLER_HPP
