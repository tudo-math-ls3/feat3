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
     * \tparam DataType_
     * Our data type
     *
     **/
    template<typename DataType_>
    class RumpfFunctionalBase
    {
      public:
        /// Our data type
        typedef DataType_ DataType;
      public:
        /// Factor for the Frobenius norm part in the mesh quality functional
        DataType _fac_norm;
        /// Factor for the det part in the mesh quality functional
        DataType _fac_det;
        /// Factor for the 1/det part in the mesh quality functional
        DataType _fac_rec_det;
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
         * \param[in] fac_rec_det_
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
          const DataType_ fac_rec_det_,
          const DataType_ fac_cof_,
          const DataType_ fac_reg_) :
          _fac_norm(fac_norm_),
          _fac_det(fac_det_),
          _fac_rec_det(fac_rec_det_),
          _fac_cof(fac_cof_),
          _fac_reg(fac_reg_)
          {
          }

        /// Destructor
        virtual ~RumpfFunctionalBase()
        {
        }

        virtual void print()
        {
          std::cout << "RumpfFunctionalBase characteristics: " << std::endl;
          std::cout << "fac_norm = " << scientify(_fac_norm) << ", fac_det = " << scientify(_fac_det) << ", fac_rec_det = " << scientify(_fac_rec_det) << std::endl;
          std::cout << "fac_cof =  " << scientify(_fac_cof) << ", fac_reg = " << scientify(_fac_reg) << std::endl;
        }

    };

    /**
     * \brief Base class template for Rumpf functionals
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
    template<typename DataType_, typename ShapeType_>
    class RumpfFunctional;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation.
#else
    template<typename DataType_, typename ShapeType_>
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
       * fac_rec_det are now coupled.
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
       * \f$ f( det(A), \| A \|^2_F) = fac_norm (\| A \|^2_F - 2)^2 + fac_det det(A) + fac_rec_det \frac{1}{det(A) + \sqrt(\mathrm{fac_reg}^2 + det^2(A))} \f$
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
       * \param[in] func_rec_det
       * Contribution from the 1/det term
       *
       * Debug variant that also returns the contributions of the different terms seperately.
       *
       **/
      template<typename Tx_, typename Th_>
      DataType_ compute_local_functional(const Tx_& x, const Th_& h,
      DataType_& func_norm,
      DataType_& func_det,
      DataType_& func_rec_det);

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
      DataType_ compute_rec_det_A( const Tx_& x, const Th_& h);

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
     * \brief Base class template for the levelset functionality of Rumpf functionals with levelsets
     *
     * This is completely seperate from RumpfFunctional.
     *
     * \tparam DataType_
     * Our data type
     *
     * \tparam ShapeType_
     * Shape type of the underlying transformation
     *
     * The actual implementation has to be supplied by the implementation. Currently, the generic implementation fits
     * all shape types.
     *
     **/
#ifndef DOXYGEN
    template<typename DataType_, typename ShapeType_>
    class RumpfFunctionalLevelset;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation.
    // Only components that are not present in the base class are documented here.
#else
    template<typename DataType_, typename ShapeType_>
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
     * \brief Class template for Rumpf functionals using the Q1 Hack
     *
     * Each Hypercube<d> is subdivided into d! simplices and the P1 functional is evaluated for each of the possible
     * ways of performing that subdivision.
     *
     * \tparam DataType_
     * Our data type
     *
     * \tparam ShapeType_
     * Shape type of the underlying transformation
     *
     * \tparam BaseClass_
     * The class of funtional that is to be used on each simplex
     *
     * The actual implementation has to be supplied by the implementation of specialisations in ShapeType_. Currently,
     * this is basically a wrapper that makes specialisations in Hypercube<d> use the specialisation of BaseClass in
     * Simplex<d>, to keep the template interface clean and the parameters consistent.
     *
     **/
    template<typename DataType_, typename ShapeType_, template<typename, typename> class BaseClass_>
    class RumpfFunctionalQ1Hack;

    /**
     * \brief Class template for Rumpf functionals with squared det and 1/det terms
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
    template<typename DataType_, typename ShapeType_>
    class RumpfFunctional_D2;

    /**
     * \brief Class template for Rumpf functionals aware of a concentration function
     *
     * The concentration function prescribes the mesh density by prescribing the local optimal scales.
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
    template<typename DataType_, typename ShapeType_>
    class RumpfFunctionalConc;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation.
#else
    template<typename DataType_, typename ShapeType_>
    class RumpfFunctionalConc{
      /**
       * \brief Adds the gradient h part to the local gradient
       *
       * Due to the chain rule, the derivatives of h wrt. the vertex coordinates x have been computed and
       * contribute to the local gradient.
       *
       * \tparam Tgrad_
       * Type of the local gradient, i.e. Tiny::Matrix
       *
       * \tparam Tx_
       * Type of the local vertex coordinates, i.e. Tiny::Matrix
       *
       * \tparam Th_
       * Type of the local optimal scales, i.e. Tiny::Vector
       *
       * \tparam Tgrad_
       * Type of the local gradient of h, i.e. Tiny::Vector of length dimension*(number of local vertices).
       * This is due to the structure of DenseVectorBlocked.
       *
       * \param grad
       * The local gradient of the functional
       *
       * \param x
       * The local vertex coordinates
       *
       * \param h
       * The local optimal scales
       *
       * \param grad_h
       * The gradient of h wrt. the local vertex coordinates
       *
       **/
      template<typename Tgrad_, typename Tx_, typename Th_, typename Tgradh_>
      void add_grad_h_part(Tgrad_& grad, const Tx_& x, const Th_& h, const Tgradh_& grad_h);
    };
#endif

    template<typename DataType_, typename ShapeType_>
    class RumpfFunctionalConc_D2;

    template<typename DataType_, typename ShapeType_>
    class RumpfFunctionalLevelsetMonitor;

  }

}
#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_HPP
