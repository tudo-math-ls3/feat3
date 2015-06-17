#pragma once
#ifndef KERNEL_ASSEMBLY_MEAN_FILTER_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_MEAN_FILTER_ASSEMBLER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/common_functionals.hpp>

namespace FEAST
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
       * \brief Assembles an intergral mean filter
       *
       * \param[out] filter
       * The filter to be assembled.
       *
       * \param[in] space
       * The finite element space for which the filter is to be assembled.
       *
       * \param[in] cubature_factory
       * A cubature factory for integration.
       *
       * \todo This functions needs to synchronise for parallel vectors
       */
      template<typename MemType_, typename DataType_, typename IndexType_, typename Space_, typename CubatureFactory_>
      static void assemble(
        LAFEM::MeanFilter<MemType_, DataType_, IndexType_>& filter,
        const Space_& space, const CubatureFactory_& cubature_factory)
      {
        // allocate primal and dual vectors
        LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> vec_v(space.get_num_dofs(), DataType_(0));
        LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> vec_w(space.get_num_dofs(), DataType_(0));

        // create a constant 1-function and its corresponding force functional
        Assembly::Common::ConstantFunction one_func(Real(1));
        Assembly::Common::ForceFunctional<Assembly::Common::ConstantFunction> one_force(one_func);

        // interpolate 1-function into vector v
        Assembly::Interpolator::project(vec_v, one_func, space);

        // assemble 1-function force into vector w
        Assembly::LinearFunctionalAssembler::assemble_vector(vec_w, one_force, space, cubature_factory);

        // compute the dot-product to obtain domain volume
        DataType_ vol = vec_w.dot(vec_v);

        // create the filter
        LAFEM::MeanFilter<Mem::Main, DataType_, IndexType_> mean_filter(std::move(vec_v), std::move(vec_w), vol);

        // convert to correct memory architecture
        filter.convert(mean_filter);
      }
    }; // class MeanFilterAssembler
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_MEAN_FILTER_ASSEMBLER_HPP
