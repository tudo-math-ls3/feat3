#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONAL_HPP
#define KERNEL_MESHOPT_RUMPF_FUNCTIONAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/mpi_cout.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <iostream>

namespace FEAT
{
  namespace Meshopt
  {
    /**
     * \brief Base class template for Rumpf functionals
     *
     * For a detailed description, \see RumpfFunctional
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
        DataType _fac_frobenius;
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
         * \param[in] fac_frobenius_
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
          const DataType_ fac_frobenius_,
          const DataType_ fac_det_,
          const DataType_ fac_rec_det_,
          const DataType_ fac_cof_,
          const DataType_ fac_reg_) :
          _fac_frobenius(fac_frobenius_),
          _fac_det(fac_det_),
          _fac_rec_det(fac_rec_det_),
          _fac_cof(fac_cof_),
          _fac_reg(fac_reg_)
          {
          }

        /// \brief Virtual destructor
        virtual ~RumpfFunctionalBase()
        {
        }

        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        static String name()
        {
          return "RumpfFunctionalBase";
        }

        /// \brief Print basic information
        void print()
        {
          Util::mpi_cout_pad_line("fac_frobenius",stringify_fp_sci(_fac_frobenius));
          Util::mpi_cout_pad_line("fac_det",stringify_fp_sci(_fac_det));
          Util::mpi_cout_pad_line("fac_rec_det",stringify_fp_sci(_fac_rec_det));
          Util::mpi_cout_pad_line("fac_cof",stringify_fp_sci(_fac_cof));
          Util::mpi_cout_pad_line("fac_reg",stringify_fp_sci(_fac_reg));
        }

    };

    /**
     * \brief Functionals for measuring and optimising mesh quality
     *
     * \tparam DataType_
     * Our data type
     *
     * \tparam TrafoType_
     * The transformation mapping reference cells to the mesh
     *
     * The actual implementation has to be supplied by the specialisations in TrafoType_.
     *
     * A RumpfFunctional computes the cell local contribution to the local functional value and its gradient.
     * The basic archetype is of the form
     *
     * \f[
     *   L( \nabla R_K ) = \int_K fac\_norm (\| \nabla R_K \|^2_F - d)^2
     *                          + fac\_cof (\| \mathrm{Cof}(\nabla R_K) \|^2_F - d)^2 +
     *                          +fac\_det (\det \nabla R_K)^p
     *   + \frac{fac\_rec\_det}{\det \nabla R_K + \sqrt(\mathrm{fac\_reg}^2 + (\det \nabla R_K)^2)}
     * \f]
     *
     * For the computation of the gradient of \f$ L \f$, we need to calculate some directional derivatives wrt.
     * matrices.
     *
     * Let \f$ M: \mathbb{R} \to \mathbb{R}^{d \times d} be a matrix valued mapping. Then we have the following
     * identities:
     * \f{align*}
     *    \frac{d}{dt}M(t)^{-1} = - M^{-1}(t) \frac{d M}{d t}(t) M^{-1}(t) \\
     *    \frac{d}{dt} \| M(t) \|_F^2 =  M(t) : \frac{d M}{d t}(t)\\
     *    \frac{d}{dt} \mathrm{det}( M(t) ) =  \mathrm{det}(M(t)) \left( M^{-1}(t) : \frac{d M}{d t}(t) \right)\\
     *    \frac{d}{dt} \mathrm{cof}( M(t) ) =  M(t)^{-T} : \frac{d M}{d t}(t) \mathrm{cof}(M(t))
     *    - M^{-T} \left(\frac{d M}{d t}(t) \right)^T \mathrm{cof}(M(t))
     * \f]
     *
     */
#ifndef DOXYGEN
    template<typename DataType_, typename TrafoType_>
    class RumpfFunctional;
    // Note:
    // The following block serves as an element interface documentation and is therefore only
    // visible to doxygen. The actual functionality has to be supplied by the implementation.
#else
    template<typename DataType_, typename TrafoType_>
    class RumpfFunctional
    {
      /**
       * \brief Constructor
       *
       * \details
       *
       * \param[in] fac_frobenius_
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
       * \f$ c_s = \frac{1}{2d (d \sqrt{ \mathrm{fac\_reg}^2 + d^2} + \mathrm{fac\_reg}^2 + d^2)} \f$ for \f$ d = 1 \f$ for p = 2.
       *
       **/
      RumpfFunctional(
        const DataType_ fac_frobenius_,
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
       * (number of vertices per cell) \f$ \times \f$ (world_dim) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$.
       *
       * \param[in] h
       * (world_dim) - vector that holds the dimensions of the target reference cell
       *
       * \returns
       * Local contribution to functional value
       *
       * Computes the local functional value for one cell defined by the coordinates x with regard
       * to the reference cell which depends on the shape type.
       *
       */
      template<typename Tx_, typename Th_>
      DataType_ eval_fval_grad(const Tx_& x, const Th_& h);

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
      DataType_ eval_fval_cellwise(const Tx_& x, const Th_& h,
      DataType_& func_norm,
      DataType_& func_det,
      DataType_& func_rec_det);

      /**
       * \brief Adds the gradient h part to the local gradient
       *
       * Due to the chain rule, the derivatives of h wrt. the vertex coordinates x have been computed and
       * contribute to the local gradient.
       *
       * The optimal scales \f$ h \f$ depend on the concentration \f$ c \f$ by the relation
       * \f[
       *   \forall K_k \in \mathcal{T}_h: h_k = \left( \frac{c (K_k)}{\sum_{l=1}^N c(K_l) \sum_{l=1}^N
       *   \mathrm{det} \nabla R_{T,l}(\phi)} \right)^{\frac{1}{d}}
       * \f]
       *
       * The concentration function in turn indirectly depends on the vertex locations
       * \f$ \{ x_i : i=1, \dots, n \} \f$ via
       * \f[ c(K_k) = (\alpha + \sum_{x_j \in K_k} \varphi(x_j))^\beta, \f]
       * so that we have to take this dependency into account for the full gradient. Define
       * \f[
       *   s_d := \sum_{l=1}^N \mathrm{det} \nabla R_{T,l}(\Phi), s_c :=  \frac{c (K_k)}{\sum_{l=1}^N c(K_l)}.
       * \f]
       * So for each \f$ x_j \f$ we arrive at
       * \f{align*}
       *   \frac{\partial h(K_k)}{\partial x_j} & = \frac{1}{d} \left( \frac{c(K_k)}{s_d s_c} \right)^
       *   {\frac{1}{d}-1} \left[ c(K_k) ( \frac{\partial s_d}{\partial x_j} s_c + s_d
       *   \frac{\partial s_c}{\partial x_j} ) + \frac{\partial c(K_k)}{\partial x_j} s_d s_c \right]
       *   (s_d s_c)^{-2}
       * \f}
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

    /**
     * \brief Class template for Rumpf functionals splitting hypercubes into simplices
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
     */
    template<typename DataType_, typename ShapeType_, template<typename, typename> class BaseClass_>
    class RumpfFunctionalQ1Split;

    /// \cond internal
    template<typename DataType_, typename ShapeType_>
    class RumpfFunctionalUnrolled;
    /// \endcond
  } // namespace Meshopt

} // namespace FEAT
#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONAL_HPP
