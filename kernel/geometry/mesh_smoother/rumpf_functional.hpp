#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_HPP
#define KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_HPP 1

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Base class template for Rumpf functionals
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

        /// Destructor
        virtual ~RumpfFunctionalBase()
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
     * \tparam ShapeType_
     * Shape type of the underlying transformation
     *
     * The actual implementation has to be supplied by the specialisations in ShapeType_.
     *
     **/
#ifndef DOXYGEN
    template<typename MemoryType_, typename DataType_, typename ShapeType_>
    class RumpfFunctional;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation.
#else
    template<typename MemoryType_, typename DataType_, typename ShapeType_>
    class RumpfFunctional
    {
      /**
       * \brief Constructor
       *
       * \details
       *
       * \param[in] fac_norm_
       * Factor for the Rumpf functional
       *
       * \param[in] fac_det_
       * Factor for the Rumpf functional
       *
       * \param[in] fac_cof_
       * Factor for the Rumpf functional
       *
       * \param[in] fac_reg_
       * Regularisation factor for the Rumpf functional
       *
       * Because we want the functional to have a minimum if the mesh consists of correctly scaled versions of the
       * optimal element (that is, every local transformation is the identity modulo translation), fac_det and
       * fac_det2 are now coupled.
       *
       * \f$ f(d) = c_s d^p + \frac{1}{d + \sqrt{ \mathrm{fac_reg}^2 + d^2}} \f$
       * We want to have \f$ f'(1) = 0 \Leftrightarrow c_s = \frac{1}{d \sqrt{ \mathrm{fac_reg}^2 + d^2} +
       * \mathrm{fac_reg}^2 + d^2} \f$ for \f$ p = 1 \f$ or
       * \f$ c_s = \frac{1}{2d (d \sqrt{ \mathrm{fac_reg}^2 + d^2} + \mathrm{fac_reg}^2 + d^2)} \f$ for \f$ d = 1 \f$ for p = 2.
       *
       **/
      RumpfFunctional(
        const DataType_ fac_norm_,
        const DataType_ fac_det_,
        const DataType_ fac_cof_,
        const DataType_ fac_reg_);

      /**
       * \brief Computes value the Rumpf functional on one element.
       *
       * \tparam Tx_
       * Type for the cell coordinates x
       *
       * \tparam Th_
       * Type for the local mesh size
       *
       * \param[in] x
       * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$.
       *
       * \param[in] h
       * (world_dim) - vector that holds the dimensions of the target reference cell
       *
       * \returns
       * Local contribution to functional value
       *
       * Computes the local functional value for one cell defined by the coordinates x with regard
       * to the reference cell which depends on the shape type. The functional value is a function of these quantities,
       * in this case
       * \f$ f( det(A), \| A \|^2_F) = fac_norm (\| A \|^2_F - 2)^2 + fac_det det(A) + fac_det2 \frac{1}{det(A) + \sqrt(\mathrm{fac_reg}^2 + det^2(A))} \f$
       *
       **/
      template<typename Tx_, typename Th_>
      DataType_ compute_local_functional(const Tx_& x, const Th_& h);

      /**
       * \copydoc compute_local_functional()
       *
       * \tparam Tx_
       * Type for the cell coordinates x
       *
       * \tparam Th_
       * Type for the local mesh size
       *
       * \param[in] func_norm
       * Contribution from the Frobenius norm
       *
       * \param[in] func_det
       * Contribution from the det term
       *
       * \param[in] func_det2
       * Contribution from the 1/det term
       *
       * Debug variant that also returns the contributions of the different terms seperately.
       *
       **/
      template<typename Tx_, typename Th_>
      DataType_ compute_local_functional(const Tx_& x, const Th_& h,
      DataType_& func_norm,
      DataType_& func_det,
      DataType_& func_det2);

      /**
       * \brief Computes the det term on one element
       *
       * \tparam Tx_
       * Type for the cell coordinates x
       *
       * \tparam Th_
       * Type for the local mesh size
       *
       * \param[in] x
       * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$.
       *
       * \param[in] h
       * (world_dim) - vector that holds the dimensions of the target reference cell
       *
       * \returns
       * The local contribution to the 1/det term
       *
       * Computes \f$ \int_K (det(A^T A))^{\frac{p}{2}} dx \f$ for the cell defined by the coordinates x with regard to the reference
       * cell \f$ K^* \f$.
       *
       **/
      template<typename Tx_, typename Th_ >
      DataType_ compute_det_A( const Tx_& x, const Th_& h);

      /**
       * \brief Computes the 1/det term on one element
       *
       * \tparam Tx_
       * Type for the cell coordinates x
       *
       * \tparam Th_
       * Type for the local mesh size
       *
       * \param[in] x
       * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$.
       *
       * \param[in] h
       * (world_dim) - vector that holds the dimensions of the target reference cell
       *
       * \returns
       * The local contribution to the 1/det term
       *
       * Computes \f$ \int_K \frac{1}{\sqrt{det(A^T A)} + \sqrt{det(A^T A) + fac_reg^2}} dx \f$ for the cell defined by the coordinates x with regard to the reference
       * cell \f$ K^* \f$.
       *
       **/
      template<typename Tx_, typename Th_ >
      DataType_ compute_det2_A( const Tx_& x, const Th_& h);

      /**
       * \brief Computes the Frobenius norm term for one cell
       *
       * Computes \f$ \int_K \| A \|^2_F  dx \f$ for the cell defined by the coordinates x with regard to the reference
       * cell \f$ K^* \f$.
       *
       * \tparam Tx_
       * Type for the cell coordinates x
       *
       * \tparam Th_
       * Type for the local mesh size h
       *
       * \param[in] x
       * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$.
       *
       * \param[in] h
       * (world_dim) - vector that holds the dimensions of the target reference cell
       *
       * \returns
       * The local contribution to the Frobenius norm term
       *
       **/
      template<typename Tx_, typename Th_ >
      DataType_ compute_norm_A( const Tx_& x, const Th_& h);

      /**
       * \brief Computes the functional gradient for one cell
       *
       * Computes the functional gradient for the cell defined by the coordinates x with regard to the reference
       * cell \f$ K^* \f$.
       *
       * \tparam Tx_
       * Type for the cell coordinates x
       *
       * \tparam Th_
       * Type for the local mesh size h
       *
       * \tparam Tgrad_
       * Type for the local gradient matrix grad
       *
       * \param[in] x
       * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$.
       *
       * \param[in] h
       * (world_dim) - vector that holds the dimensions of the target reference cell
       *
       * \param[in] grad
       * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local contribution to the
       * global functional gradient
       *
       * \returns
       * The local contribution grad to global functional gradient
       *
       **/
      template<typename Tx_, typename Th_, typename Tgrad_>
      void compute_local_grad( const Tx_& x, const Th_& h, Tgrad_& grad);
    };
#endif

    /**
     * \brief Base class template for Rumpf functionals with levelsets
     *
     * This inherits from RumpfFunctional.
     *
     * \tparam MemoryType_
     * Memory architecture
     *
     * \tparam DataType_
     * Our data type
     *
     * \tparam ShapeType_
     * Shape type of the underlying transformation
     *
     * The actual implementation has to be supplied by the implementation of specialisations in ShapeType_.
     *
     **/
#ifndef DOXYGEN
    template<typename MemoryType_, typename DataType_, typename ShapeType_>
    class RumpfFunctionalLevelset;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation.
    // Only components that are not present in the base class are documented here.
#else
    template<typename MemoryType_, typename DataType_, typename ShapeType_>
    class RumpfFunctionalLevelset
    {
      /**
       * \brief Computes the additional levelset penalty term for the Rumpf functional
       *
       * \tparam Tl_
       * Type for the object containing the local levelset values, i.e. FEAST::Tiny::Vector
       *
       * \param[in] lvlset_vals
       * The values of the levelset function in the grid vertices for the current element
       *
       * \returns
       * The local contribution to the global levelset penalty term
       *
       **/
      template<typename Tl_>
      DataType_ compute_lvlset_penalty(const Tl_& lvlset_vals);

      /**
       * \brief Adds the gradient of the additional levelset penalty term
       *
       * \tparam Tl_
       * Type for the object containing the local levelset values, i.e. FEAST::Tiny::Vector
       *
       * \tparam Tlg_
       * Type for the object containing the local levelset gradient values, i.e. FEAST::Tiny::Matrix
       *
       * \tparam Tgrad_
       * Type for the object containing the functional gradient, i.e. FEAST::Tiny::Matrix
       *
       * \param[in] lvlset_vals
       * The values of the levelset function in the grid vertices for the current element
       *
       * \param[in] lvlset_grad_vals
       * The values of the levelset function's gradient in the grid vertices for the current element
       *
       * \param[in] grad
       * Local functional gradient to which the gradient of the levelset penalty term is added
       *
       * \param[in] lvlset_constraint_last
       * Last computed levelset constraint
       *
       **/
      template<typename Tl_, typename Tlg_, typename Tgrad_>
      void add_lvlset_penalty_grad(const Tl_& lvlset_vals, const Tlg_& lvlset_grad_vals, Tgrad_& grad, DataType_& lvlset_constraint_last);
    };
#endif
    /**
     * \brief Base class template for Rumpf functionals with levelsets that use the Q1 Hack
     *
     * Each Hypercube<d> is subdivided into d! simplices and the P1 functional is evaluated for each of the possible
     * ways of performing that subdivision.
     *
     * \tparam MemoryType_
     * Memory architecture
     *
     * \tparam DataType_
     * Our data type
     *
     * \tparam ShapeType_
     * Shape type of the underlying transformation
     *
     * The actual implementation has to be supplied by the implementation of specialisations in ShapeType_.
     *
     **/
    template<typename MemoryType_, typename DataType_, typename ShapeType_>
    class RumpfFunctionalLevelsetQ1Hack;

    template<typename MemoryType_, typename DataType_, typename ShapeType_>
    class RumpfFunctional_D2;

    template<typename MemoryType_, typename DataType_, typename ShapeType_>
    class RumpfFunctionalLevelset_D2;
    }
  }

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_HPP
