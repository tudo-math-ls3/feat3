#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONAL_HPP
#define KERNEL_MESHOPT_RUMPF_FUNCTIONAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/math.hpp>
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
     */
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
            XASSERTM(_fac_frobenius >= DataType_(0), "_fac_frobenius must be >= 0!\n");
            XASSERTM(_fac_det >= DataType_(0), "_fac_det must be >= 0!\n");
            XASSERTM(_fac_rec_det >= DataType_(0), "_fac_rec_det must be >= 0!\n");
            XASSERTM(_fac_cof >= DataType_(0), "_fac_cof must be >= 0!\n");
            XASSERTM(_fac_reg >= DataType_(0), "_fac_reg must be >= 0!\n");
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
          Index pad_width(30);
          Dist::Comm comm_world(Dist::Comm::world());

          String msg;

          msg = String("fac_frobenius").pad_back(pad_width, '.') + String(": ") + stringify_fp_sci(_fac_frobenius);
          comm_world.print(msg);

          msg = String("fac_cof").pad_back(pad_width, '.') + String(": ") + stringify_fp_sci(_fac_cof);
          comm_world.print(msg);

          msg = String("fac_det").pad_back(pad_width, '.') + String(": ") + stringify_fp_sci(_fac_det);
          comm_world.print(msg);

          msg = String("fac_rec_det").pad_back(pad_width, '.') + String(": ") + stringify_fp_sci(_fac_rec_det);
          comm_world.print(msg);

          msg = String("fac_reg").pad_back(pad_width, '.') + String(": ") + stringify_fp_sci(_fac_reg);
          comm_world.print(msg);

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
     *                          + fac\_cof (\| \mathrm{Cof}(\nabla R_K) \|^2_F - d)^2
     *                          + fac\_det (\det \nabla R_K)^{p_d}
     *   + \frac{fac\_rec\_det}{A\det \nabla R_K + \sqrt{(\mathrm{fac\_reg}^2 + (\det \nabla R_K)^2)}}
     * \f]
     *
     * For the computation of the gradient of \f$ L \f$, we need to calculate some directional derivatives wrt.
     * matrices.
     *
     * Let \f$ M: \mathbb{R} \to \mathbb{R}^{d \times d} \f$ be a matrix valued mapping. Then we have the following
     * identities:
     * \f{align*}
     *    \frac{d}{dt} M(t)^{-1} &= - M^{-1}(t) \frac{d M}{d t}(t) M^{-1}(t) \\
     *    \frac{d}{dt} \| M(t) \|_F^2 &=  M(t) : \frac{d M}{d t}(t)\\
     *    \frac{d}{dt} \mathrm{cof}( M(t) ) &=  M(t)^{-T} : \frac{d M}{d t}(t) \mathrm{cof}(M(t))
     *    - M^{-T} \left(\frac{d M}{d t}(t) \right)^T \mathrm{cof}(M(t)) \\
     *    \frac{d}{dt} \mathrm{det}( M(t) ) &=  \mathrm{det}(M(t)) \left( M^{-1}(t) : \frac{d M}{d t}(t) \right)
     * \f}
     * With this, we can compute directional derivatives of a function \f$ f: \mathbb{R}^{d \times d} \to V \f$ with
     * e.g. \f$ V = \mathbb{R} \f$ or \f$ V = \mathbb{R}^{d \times d} \f$ with regard to matrices according to
     * \f[
     *   M(t) = A + tB, \quad f'(A)B = \left. \frac{d f(M(t))}{dt}\right|_{t = 0}.
     * \f]
     * Here, \f$ M(0) = A \f$ and \f$ \frac{d M}{d t}(t) = B \f$. For computing the directional derivatives of this
     * RumpfFunctional, set \f$ A = \nabla R_K \f$ (the current FE function) and \f$ B = \nabla \eta \f$
     * (test function). Because the form is nonlinear, it cannot be represented as a matrix where the columns are
     * associated with the trial functions.
     *
     * So we get
     * \f{align*}
     *   \frac{d}{dt} \left( (\| \nabla R_K \|_F^2 - d)^2\right)' \eta &= 4(\| \nabla R_K \|_F^2 - d)
     *   \nabla R_K : \nabla \eta \\
     *   \frac{d}{dt} \left( (\| \mathrm{cof}(\nabla R_K) \|_F^2 - d)^2\right)' \eta & =
     *     4(\| \mathrm{cof}(\nabla R_K) \|_F^2 - d) ~ \mathrm{cof}(\nabla R_K) :
     *      \left( \left( \left( (\nabla R_K)^{-1} :\nabla \eta ~ I_d - (\nabla R_K)^{-1} \nabla \eta \right) \right) \mathrm{cof}(\nabla R_K) \right) \\
     *   \frac{d}{dt} \left( \mathrm{det}(\nabla R_K)^p \right)' \eta & = p \mathrm{det}(\nabla R_K )^p
     *   \nabla R_K)^{-1} : \nabla \eta.
     * \f}
     *
     * Note that the test functiona \f$ \eta \f$ might not be the standard FE basis functions as defined by the FE
     * space used (call those \f$ \hat{\eta} \f$). They might be with regard to a different reference cell, so a
     * transformation is involved, which is realised by passing the [c] material tensor [/c] \f$ M \f$ (which can be
     * identified with a matrix in this case) to corresponding functions and using the relation
     * \f[
     *   \nabla \eta = M^T \nabla \hat{\eta},
     * \f]
     * which is nothing but the chain rule applied to \f$ \eta := \hat{\eta} \circ M\f$.
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
      public:
        /// Our baseclass
        typedef RumpfFunctionalBase<DataType_> BaseClass;

        /// Our data type
        typedef DataType_ DataType;
        /// Our shape dimension
        static constexpr int shape_dim = shape_dim_;
        /// Our world dimension - only world_dim == shape_dim is supported
        static constexpr int world_dim = shape_dim;
        /// Shape type of the underlying transformation
        typedef Shape::Hypercube<shape_dim_> ShapeType;
        /// The transformation this functional works on
        typedef Trafo::Standard::Mapping<Geometry::ConformalMesh<ShapeType, world_dim, world_dim, DataType_>> TrafoType;
        /// The FE space associated with the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space TrafoSpace;

        /// Type for a pack of local vertex coordinates
        typedef Tiny::Matrix<DataType_, Shape::FaceTraits<ShapeType,0>::count, world_dim> Tx;
        /// Type for the gradient of the local cell sizes
        typedef Tiny::Vector<DataType_, Shape::FaceTraits<ShapeType,0>::count*world_dim> Tgradh;

        /// Type for the gradient of the local reference mapping
        typedef Tiny::Matrix<DataType_, shape_dim, world_dim> TgradR;

        /// Type for evaluating the transformation
        typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;
        /// Type for evaluating the FE functions
        typedef typename TrafoSpace::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;

        /// We need the Trafo to evaluate the image point, the Jacobian and its determinant
        static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_mat | TrafoTags::jac_det;
        /// For the FE space, we only need the gradients on the reference cell
        static constexpr SpaceTags space_config = SpaceTags::ref_grad;

        /// Data from evaluating the transformation
        typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;
        /// Data from evaluating FE spaces
        typedef typename SpaceEvaluator::template ConfigTraits<space_config>::EvalDataType SpaceEvalData;

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
       * \param[in] exponent_det_
       * The exponent for the det term, called \f$ p_d \f$ above.
       *
       * Because we want the functional to have a minimum if the mesh consists of correctly scaled versions of the
       * optimal element (that is, every local transformation is the identity modulo translation), fac_det and
       * fac_rec_det are now coupled.
       *
       * \f$ f(d) = c_s d^p + \frac{1}{d + \sqrt{ \mathrm{fac\_reg}^2 + d^2}} \f$
       * We want to have \f$ f'(1) = 0 \Leftrightarrow c_s = \frac{1}{d \sqrt{ \mathrm{fac\_reg}^2 + d^2} +
       * \mathrm{fac\_reg}^2 + d^2} \f$ for \f$ p = 1 \f$ or
       * \f$ c_s = \frac{1}{2d (d \sqrt{ \mathrm{fac\_reg}^2 + d^2} + \mathrm{fac\_reg}^2 + d^2)} \f$ for \f$ d = 1 \f$ for p = 2.
       *
       **/
      RumpfFunctional(
        const DataType_ fac_frobenius_,
        const DataType_ fac_det_,
        const DataType_ fac_cof_,
        const DataType_ fac_reg_,
        const int exponent_det_);

      /**
       * \brief Computes the functional's value and gradient on one cell
       *
       * \param[out] fval
       * The functional value on this cell.
       *
       * \param[out] grad
       * The contribution of this cell to the global gradient.
       *
       * \param[in] mat_tensor
       * The material function.
       *
       * \param[in] trafo_eval
       * Evaluator for the transformation.
       *
       * \param[in] space_eval
       * Evaluator for the FE spaces.
       *
       * \param[in] x
       * (number of vertices per cell) \f$ \times \f$ (world_dim) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$. These are the coefficients for the FE function representing the transformation.
       *
       * \param[in] h
       * (world_dim) - vector that holds the dimensions of the target reference cell
       *
       * Computes the local functional value for one cell defined by the coordinates x with regard
       * to the reference cell which depends on the shape type.
       *
       */
      void eval_fval_grad( DataType& DOXY(fval), Tx& DOXY(grad), const TgradR& DOXY(mat_tensor),
      const TrafoEvaluator& DOXY(trafo_eval), const SpaceEvaluator& DOXY(space_eval),
      const Tx& DOXY(x), const DataType& DOXY(h));

      /**
       * \brief Computes the different contributions to the functional's value on one cell
       *
       * \param[out] fval
       * The functional value on this cell.
       *
       * \param[in] mat_tensor
       * The material function.
       *
       * \param[in] trafo_eval
       * Evaluator for the transformation.
       *
       * \param[in] space_eval
       * Evaluator for the FE spaces.
       *
       * \param[in] x
       * (number of vertices per cell) \f$ \times \f$ (world_dim) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$. These are the coefficients for the FE function representing the transformation.
       *
       * \param[in] h
       * The local optimal scale.
       *
       * \param[out] fval_frobenius
       * Contribution of \f$ \| \nabla R_K \|_F^2 \f$.
       *
       * \param[out] fval_cof
       * Contribution of \f$ \| \mathbb{cof}(\nabla R_K) \|_F^2 \f$.
       *
       * \param[out] fval_det
       * Contribution of \f$ \mathbb{det}(\nabla R_K)\f$.
       *
       * \note: This is for debugging and visualisation purposes and does not compute the local gradient.
       *
       */
      void eval_fval_cellwise( DataType& DOXY(fval), const TgradR& DOXY(mat_tensor),
      const TrafoEvaluator& DOXY(trafo_eval), const SpaceEvaluator& DOXY(space_eval),
      const Tx& DOXY(x), const DataType& DOXY(h),
      DataType& DOXY(fval_frobenius), DataType& DOXY(fval_cof), DataType& DOXY(fval_det));

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
       * \param grad
       * The local gradient of the functional to which the contribution is added.
       *
       * \param[in] mat_tensor
       * The material function.
       *
       * \param[in] trafo_eval
       * Evaluator for the transformation.
       *
       * \param[in] space_eval
       * Evaluator for the FE spaces.
       *
       * \param[in] x
       * (number of vertices per cell) \f$ \times \f$ (world_dim) - matrix that holds the local coordinates that
       * define the cell \f$ K \f$. These are the coefficients for the FE function representing the transformation.
       *
       * \param[in] h
       * The local optimal scale
       *
       * \param grad_h
       * The gradient of h wrt. the local vertex coordinates
       *
       */
        void add_grad_h_part(Tx& grad, const TgradR& mat_tensor, const TrafoEvaluator& trafo_eval, const SpaceEvaluator& space_eval, const Tx& DOXY(x), const DataType& h, const Tgradh& grad_h)

    };
#endif

    /// \cond internal
    // These functionals are kept only for debugging purposes
    template<typename DataType_, typename TrafoType_>
    class RumpfFunctionalUnrolled;
    /// \endcond
  } // namespace Meshopt

} // namespace FEAT
#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONAL_HPP
