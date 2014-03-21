#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_HPP
#define KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_HPP 1

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief base class template for Rumpf functionals
     *
     * \tparam MemoryType_
     * Memory architecture
     *
     * \tparam DataType_
     * Our data type
     *
     **/
    template<typename MemoryType_, typename DataType_>
    class RumpfFunctionalBase
    {
      public:
        /// Our data type
        typedef DataType_ DataType;
      protected:
        /// Factor for the Frobenius norm part in the mesh quality functional
        DataType _fac_norm;
        /// Factor for the det part in the mesh quality functional
        DataType _fac_det;
        /// Factor for the 1/det part in the mesh quality functional
        DataType _fac_det2;
        /// Factor for the cof part in the mesh quality functional
        DataType _fac_cof;
        /// Regularisation parameter for the 1/det part in the mesh quality functional
        DataType _fac_reg;


        /** \brief Constructor
         *
         * \details
         *
         * \param[in] fac_norm_
         * Factor for the Rumpf functional
         *
         * \param[in] fac_det_
         * Factor for the Rumpf functional
         *
         * \param[in] fac_det2_
         * Factor for the Rumpf functional
         *
         * \param[in] fac_cof_
         * Factor for the Rumpf functional
         *
         * \param[in] fac_reg_
         * Regularisation factor for the Rumpf functional
         *
         */
        explicit RumpfFunctionalBase(
            const DataType_ fac_norm_,
            const DataType_ fac_det_,
            const DataType_ fac_det2_,
            const DataType_ fac_cof_,
            const DataType_ fac_reg_) :
          _fac_norm(fac_norm_),
          _fac_det(fac_det_),
          _fac_det2(fac_det2_),
          _fac_cof(fac_cof_),
          _fac_reg(fac_reg_)
      {
      }

        virtual ~RumpfFunctionalBase()
        {
        }

        /**
         * \brief Prepares the functional for evaluation
         *
         * This could be something like evaluating the levelset function on the new mesh.
         *
         **/
        virtual void prepare()
        {
        }
    };

    /**
     * \brief Base class template for Rumpf functionals
     *
     * \tparam MemoryType_
     * Memory architecture
     *
     * \tparam DataType_
     * Our data type
     *
     * \tparam TrafoType_
     * The underlying transformation
     *
     * The actual implementation has to be supplied by the implementation of specialisations in TrafoType_.
     *
     **/
    template<typename MemoryType_, typename DataType_, typename TrafoType_>
    class RumpfFunctional;

    /**
     * \brief Base class template for Rumpf functionals with levelsets
     *
     * \tparam MemoryType_
     * Memory architecture
     *
     * \tparam DataType_
     * Our data type
     *
     * \tparam TrafoType_
     * The underlying transformation
     *
     * The actual implementation has to be supplied by the implementation of specialisations in TrafoType_.
     *
     **/
    template<typename MemoryType_, typename DataType_, typename TrafoType_>
    class RumpfFunctionalLevelset;
  }
}

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_HPP
