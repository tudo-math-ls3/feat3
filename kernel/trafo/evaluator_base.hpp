#pragma once
#ifndef KERNEL_TRAFO_EVALUATOR_BASE_HPP
#define KERNEL_TRAFO_EVALUATOR_BASE_HPP 1

// includes, FEAST
#include <kernel/trafo/eval_data.hpp>

namespace FEAST
{
  namespace Trafo
  {
    /**
     * \brief Trafo Evaluator CRTP base-class template
     *
     * This class template is a CRTP base class used by various transformation evaluators to outsource
     * common wrapper code which is independent of the actual transformation in use.
     *
     * \tparam Evaluator_
     * The evaluator class that derives from this class template.
     *
     * \tparam EvalPolicy_
     * The evaluation policy of the trafo evaluator.
     *
     * \author Peter Zajac
     */
    template<
      typename Trafo_,
      typename Evaluator_,
      typename EvalPolicy_>
    class EvaluatorBase
    {
    public:
      /// trafo type
      typedef Trafo_ TrafoType;
      /// shape type
      typedef typename TrafoType::ShapeType ShapeType;

      /// evaluation policy
      typedef EvalPolicy_ EvalPolicy;

      /// evaluation data type
      typedef typename EvalPolicy::DataType DataType;

      /// trafo coefficient type
      typedef typename EvalPolicy::TrafoCoeffType TrafoCoeffType;

      /// domain coordinate type
      typedef typename EvalPolicy::DomainCoordType DomainCoordType;
      /// domain point type
      typedef typename EvalPolicy::DomainPointType DomainPointType;
      /// domain point reference
      typedef typename EvalPolicy::DomainPointRef DomainPointRef;
      /// domain point const-reference
      typedef typename EvalPolicy::DomainPointConstRef DomainPointConstRef;

      /// image coordinate type
      typedef typename EvalPolicy::ImageCoordType ImageCoordType;
      /// image point type
      typedef typename EvalPolicy::ImagePointType ImagePointType;
      /// image point reference
      typedef typename EvalPolicy::ImagePointRef ImagePointRef;
      /// image point const-reference
      typedef typename EvalPolicy::ImagePointConstRef ImagePointConstRef;

      /// jacobian matrix coefficient type
      typedef typename EvalPolicy::JacMatCoeff JacMatCoeff;
      /// jacobian matrix type
      typedef typename EvalPolicy::JacMatType JacMatType;
      /// jacobian matrix reference
      typedef typename EvalPolicy::JacMatRef JacMatRef;
      /// jacobian matrix const-reference
      typedef typename EvalPolicy::JacMatConstRef JacMatConstRef;

      /// jacobian determinant type
      typedef typename EvalPolicy::JacDetType JacDetType;

      /// dummy enumeration
      enum
      {
        /// domain dimension
        domain_dim = EvalPolicy::domain_dim,
        /// image dimension
        image_dim = EvalPolicy::image_dim
      };

      // Note:
      // The following block serves as an element interface documentation and is therefore only
      // visible to doxygen. The actual functionality has to be supplied by the implementation.
#ifdef DOXYGEN
      /**
       * \brief Capability enumeration
       *
       * This enumeration specifies the capabilites of the evaluator, i.e. its values specifiy whether
       * the evaluator is capable of computing image point, jacobian matrices, etc.
       */
      enum EvaluatorCapabilities
      {
        /**
         * \brief Domain-Point capability
         *
         * This entry specifies whether the evaluator is capable of computing domain point coordinates.\n
         * This value is always non-zero.
         */
        can_dom_point = ...,

        /**
         * \brief Image-Point capability
         *
         * This entry specifies whether the evaluator is capable of computing image point coordinates.\n
         * If this value is non-zero, the evaluator implements the #map_point member function.\n
         * See #map_point for details.
         */
        can_img_point = ...,

        /**
         * \brief Jacobian-Matrix capability
         *
         * This entry specifies whether the evaluator is capable of computing jacobian matrices.\n
         * If this value is non-zero, the evaluator implements the #calc_jac_mat member function.\n
         * See #calc_jac_mat for details.
         */
        can_jac_mat = ...,

        /**
         * \brief Jacobian-Inverse-Matrix capability
         *
         * This entry specifies whether the evaluator is capable of computing jacobian inverse matrices.\n
         * If this value is non-zero, the evaluator implements the #calc_jac_inv member function.\n
         * See #calc_jac_inv for details.
         */
        can_jac_inv = ...,

        /**
         * \brief Jacobian-Determinant capability
         *
         * This entry specifies whether the evaluator is capable of computing jacobian determinants.\n
         * If this value is non-zero, the evaluator implements the #calc_jac_det member function.\n
         * See #calc_jac_det for details.
         */
        can_jac_det = ...,
      };
#endif // DOXYGEN

    protected:
      /// \cond internal
      Evaluator_& me()
      {
        return static_cast<Evaluator_&>(*this);
      }

      const Evaluator_& me() const
      {
        return static_cast<const Evaluator_&>(*this);
      }
      /// \endcond

      /// trafo reference
      const TrafoType& _trafo;

      /// currently active cell index
      Index _cell_index;

      /// constructor
      explicit EvaluatorBase(const TrafoType& trafo) :
        _trafo(trafo),
        _cell_index(~Index(0))
      {
      }

    public:
      /**
       * \brief Returns a reference to the trafo object.
       */
      const TrafoType& get_trafo() const
      {
        return _trafo;
      }

      /**
       * \brief Returns the number of cells in the mesh.
       *
       * \returns
       * The number of cells in the underlying mesh.
       */
      Index get_num_cells() const
      {
        return _trafo.get_mesh().get_num_entities(ShapeType::dimension);
      }

      /**
       * \brief Returns the index of the currently active cell.
       */
      Index get_cell_index() const
      {
        return _cell_index;
      }

      /**
       * \brief Prepares the evaluator for a given cell.
       *
       * \param[in] cell_index
       * The index of the cell for which the evaluator is to be prepared.
       */
      void prepare(Index cell_index)
      {
        // store cell index
        _cell_index = cell_index;
      }

      /**
       * \brief Finishes the evaluator for the currently active cell.
       */
      void finish()
      {
        // reset cell index
        _cell_index = ~Index(0);
      }

      // Note:
      // The following block serves as an element interface documentation and is therefore only
      // visible to doxygen. The actual functionality has to be supplied by the implementation.
#ifdef DOXYGEN
      /**
       * \brief Maps a domain point from the reference cell to the currently active cell.
       *
       * \param[out] img_point
       * A reference to the image point on the active cell that is to be computed.
       *
       * \param[in] dom_point
       * A reference to the domain point on the reference cell that is to be mapped.
       */
      void map_point(ImagePointRef img_point, DomainPointConstRef dom_point) const;

      /**
       * \brief Calculates the jacobian matrix for a given domain point.
       *
       * \param[out] jac_mat
       * A reference to the jacobian matrix that is to be computed.
       *
       * \param[in] dom_point
       * A reference to the domain point on the reference cell for which the jacobian matrix is to be computed.
       */
      void calc_jac_mat(JacMatRef jac_mat, DomainPointConstRef dom_point) const;
#endif // DOXYGEN

      /**
       * \brief Computes the jacobian determinant for a domain point.
       *
       * \param[in] dom_point
       * The domain point at which the jacobian determinant is to be computed.
       *
       * \returns
       * The jacobian determinant at the domain point.
       */
      JacDetType calc_jac_det(DomainPointConstRef dom_point) const
      {
        // compute jacobian matrix
        JacMatType jac_mat;
        me().calc_jac_mat(jac_mat, dom_point);

        // return its volume
        return jac_mat.vol();
      }

      /**
       * \brief Computes the inverse jacobian matrix for a domain point.
       *
       * \param[out] jac_inv
       * A reference to the inverse jacobian matrix that is to be computed.
       *
       * \param[in] dom_point
       * The domain point at which the jacobian matrix is to be computed.
       */
      void calc_jac_inv(JacMatRef jac_inv, DomainPointConstRef dom_point) const
      {
        // compute jacobian matrix
        JacMatType jac_mat;
        me().calc_jac_mat(jac_mat, dom_point);

        // invert jacobian matrix
        jac_inv.set_inverse(jac_mat);
      }
    }; // class EvaluatorBase<...>
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_EVALUATOR_BASE_HPP
