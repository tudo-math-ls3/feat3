// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONALS_2D_P1_UNROLLED_HPP
#define KERNEL_MESHOPT_RUMPF_FUNCTIONALS_2D_P1_UNROLLED_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>
#include <kernel/meshopt/base.hpp>
#include <kernel/meshopt/rumpf_functional.hpp>

namespace FEAT
{
  namespace Meshopt
  {
    /// \cond internal

    /**
     * \brief Class template for Rumpf functionals, standard trafo on conformal Simplex<2> meshes.
     *
     * This is the variant where the local functional's contributions are computed directly with Maple-generated
     * code, which means there is excessive manual loop unrolling.
     */
    template<typename DataType_>
    class RumpfFunctionalUnrolled<DataType_,
    Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, DataType_>>> :
      public RumpfFunctionalBase<DataType_>
      {
        public:
          /// Our baseclass
          typedef RumpfFunctionalBase<DataType_> BaseClass;

          /// Our data type
          typedef DataType_ DataType;
          /// Shape type of the underlying transformation
          typedef Shape::Simplex<2> ShapeType;
          /// The transformation this functional works on
          typedef Trafo::Standard::Mapping<Geometry::ConformalMesh<ShapeType, 2, DataType_>> TrafoType;
          /// Our world dimension
          static constexpr int world_dim = TrafoType::world_dim;
          /// Our shape dimension
          static constexpr int shape_dim = ShapeType::dimension;
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

        private:
          /// Exponent for the det based terms
          const int _exponent_det;
          /// Do we need to compute \f$ \| \nabla R_K \|_F \f$?
          const bool _compute_frobenius;
          /// Do we need to compute \f$ \| \mathrm{Cof}(\nabla R_K) \|_F \f$?
          const bool _compute_cof;
          /// Do we need to compute \f$ det(\nabla R_K) \f$?
          const bool _compute_det;

        public:
          /**
           * \brief Constructor
           */
          explicit RumpfFunctionalUnrolled(
            const DataType fac_frobenius_,
            const DataType fac_det_,
            const DataType fac_cof_,
            const DataType fac_reg_,
            const int exponent_det_) :
            BaseClass( fac_frobenius_,
            fac_det_,
            fac_det_*(Math::sqrt( Math::sqr(fac_reg_) + DataType(1) )*Math::pow( DataType(1) + Math::sqrt(Math::sqr(fac_reg_) + DataType(1)), DataType(exponent_det_))),
            fac_cof_,
            fac_reg_),
            _exponent_det(exponent_det_),
            _compute_frobenius( (fac_frobenius_ > DataType(0)) ),
            _compute_cof( false ),
            _compute_det( (fac_det_ > 0) )
            {
              XASSERTM(exponent_det_ == 1 || exponent_det_ == 2,"exponent_det must be 1 or 2!");
              XASSERTM(fac_cof_ == DataType(0), "In 2d, the cofactor term is redundant, so set fac_cof == 0.");
            }

          /**
           * \brief The class name
           *
           * \returns String with the class name
           */
          static String name()
          {
            return "RumpfFunctionalUnrolled<"+ShapeType::name()+">";
          }

          /**
           * \brief Prints object parameters
           */
          String info() const
          {
            const Index pad_width(30);
            return name() + ":" + BaseClass::info() + String("\nexponent_det").pad_back(pad_width, '.')
              + String(": ") + stringify(_exponent_det);
          }

          /**
           * \brief Computes the functional gradient for one cell
           */
          void eval_fval_grad(DataType& fval, Tx& grad, const TgradR& DOXY(mat_tensor), const TrafoEvaluator& DOXY(trafo_eval), const SpaceEvaluator& DOXY(space_eval), const Tx& x, const DataType& h)
          {

            fval = DataType(0);
            grad.format(DataType(0));

            DataType fval_frobenius(0);
            DataType fval_det(0);
            DataType fval_rec_det(0);

            if(_compute_frobenius)
            {
              fval_frobenius = this->_fac_frobenius*compute_frobenius_part(x,h);
              add_grad_frobenius_part(grad, x, h);
            }

            if(_compute_det)
            {
              if(_exponent_det == 1)
              {
                fval_det = this->_fac_det*compute_det_1_part(x,h);
                fval_rec_det = this->_fac_rec_det*compute_rec_det_1_part(x,h);
                add_grad_det_1_part(grad, x, h);
                add_grad_rec_det_1_part(grad, x, h);
              }
              else
              {
                fval_det = this->_fac_det*compute_det_2_part(x,h);
                fval_rec_det = this->_fac_rec_det*compute_rec_det_2_part(x,h);
                add_grad_det_2_part(grad, x, h);
                add_grad_rec_det_2_part(grad, x, h);
              }
            }

            fval = fval_frobenius + fval_det + fval_rec_det;

          }

          void eval_fval_cellwise(DataType& fval, const TgradR& DOXY(mat_tensor), const TrafoEvaluator& DOXY(trafo_eval), const SpaceEvaluator& DOXY(space_eval), const Tx& x, const DataType& h, DataType& fval_frobenius, DataType& fval_cof, DataType& fval_det)
          {
            fval = DataType(0);
            fval_frobenius = DataType(0);
            fval_cof = DataType(0);
            fval_det = DataType(0);

            if(_compute_frobenius)
            {
              fval_frobenius = this->_fac_frobenius*compute_frobenius_part(x,h);
            }

            if(_compute_det)
            {
              if(_exponent_det == 1)
              {
                fval_det = this->_fac_det*compute_det_1_part(x,h);
                fval_det += this->_fac_rec_det*compute_rec_det_1_part(x,h);
              }
              else
              {
                fval_det = this->_fac_det*compute_det_2_part(x,h);
                fval_det += this->_fac_rec_det*compute_rec_det_2_part(x,h);
              }
            }

            fval = fval_frobenius + fval_det;

          }

          /**
           * \brief Adds the part coming from the chain rule and the derivative wrt. h
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          void add_grad_h_part(Tx& grad, const TgradR& DOXY(mat_tensor), const TrafoEvaluator& DOXY(trafo_eval), const SpaceEvaluator& DOXY(space_eval), const Tx& x, const DataType& h, const Tgradh& grad_h)
#else
            void NOINLINE add_grad_h_part(Tx& grad, const TgradR& DOXY(mat_tensor), const TrafoEvaluator& DOXY(trafo_eval), const SpaceEvaluator& DOXY(space_eval), const Tx& x, const DataType& h, const Tgradh& grad_h)
#endif
            {

              DataType frobenius_der_h(0);
              DataType det_der_h(0);
              DataType rec_det_der_h(0);

              if(_compute_frobenius)
              {
                frobenius_der_h = this->_fac_frobenius*compute_frobenius_der_h_part(x,h);
              }

              if(_compute_det)
              {
                if(_exponent_det == 1)
                {
                  det_der_h = this->_fac_det * compute_det_1_der_h_part(x,h);
                  rec_det_der_h = this->_fac_rec_det * compute_rec_det_1_der_h_part(x,h);
                }
                else
                {
                  det_der_h = this->_fac_det * compute_det_2_der_h_part(x,h);
                  rec_det_der_h = this->_fac_rec_det * compute_rec_det_2_der_h_part(x,h);
                }

              }

              DataType der_h(frobenius_der_h + det_der_h + rec_det_der_h);

              for(int i(0); i < Tx::m; ++i)
              {
                for(int d(0); d < Tx::n; ++d)
                {
                  grad(i,d) += der_h*grad_h(i*Tx::n + d);
                }
              }
            } // add_grad_h_part


          /**
           * \brief Computes the Frobenius norm term for one cell
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_frobenius_part(const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_frobenius_part(const Tx& x, const DataType& h)
#endif
            {
              DataType fval_frobenius_part(0);
              fval_frobenius_part = DataType(4)/DataType(9)*Math::pow(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)),DataType(2))/(h*h*h*h);
              return fval_frobenius_part;
            }

          /**
           * \brief Computes the det term on one element
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_det_1_part( const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_det_1_part( const Tx& x, const DataType& h)
#endif
            {
              DataType fval_det_1_part(0);
              fval_det_1_part = DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h);
              return fval_det_1_part;
            }
          /**
           * \brief Computes the det term on one element
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_det_2_part( const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_det_2_part( const Tx& x, const DataType& h)
#endif
            {
              DataType fval_det_2_part(0);
              fval_det_2_part = DataType(4)/DataType(3)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h);
              return fval_det_2_part;
            }

          /**
           * \brief Computes the 1/det term on one element
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_rec_det_1_part( const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_rec_det_1_part( const Tx& x, const DataType& h)
#endif
            {
              DataType fval_rec_det_1_part(0);
              fval_rec_det_1_part = DataType(1)/(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3));
              return fval_rec_det_1_part;
            }

          /**
           * \brief Computes the 1/det term on one element
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_rec_det_2_part( const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_rec_det_2_part( const Tx& x, const DataType& h)
#endif
            {
              DataType fval_rec_det_2_part(0);
              fval_rec_det_2_part = Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(2));
              return fval_rec_det_2_part;
            }
          /**
           * \brief Computes the Frobenius norm term for one cell
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_frobenius_der_h_part(const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_frobenius_der_h_part(const Tx& x, const DataType& h)
#endif
            {
              DataType frobenius_der_h_part(0);
              frobenius_der_h_part = DataType(16)/DataType(3)*(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)))/(h*h*h)-DataType(16)/DataType(9)*Math::pow(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)),DataType(2))/(h*h*h*h*h);
              return frobenius_der_h_part;
            }

          /**
           * \brief Computes the det term on one element
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_det_1_der_h_part( const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_det_1_der_h_part( const Tx& x, const DataType& h)
#endif
            {
              DataType det_1_der_h_part(0);
              det_1_der_h_part = -DataType(4)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h);
              return det_1_der_h_part;
            }
          /**
           * \brief Computes the det term on one element
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_det_2_der_h_part( const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_det_2_der_h_part( const Tx& x, const DataType& h)
#endif
            {
              DataType det_2_der_h_part(0);
              det_2_der_h_part = -DataType(16)/DataType(3)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h*h);
              return det_2_der_h_part;
            }

          /**
           * \brief Computes the 1/det term on one element
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_rec_det_1_der_h_part( const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_rec_det_1_der_h_part( const Tx& x, const DataType& h)
#endif
            {
              DataType rec_det_1_der_h_part(0);
              rec_det_1_der_h_part = -Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(2))*(-DataType(4)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h)-DataType(8)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h*h));
              return rec_det_1_der_h_part;
            }

          /**
           * \brief Computes the 1/det term on one element
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          DataType compute_rec_det_2_der_h_part( const Tx& x, const DataType& h)
#else
            DataType NOINLINE compute_rec_det_2_der_h_part( const Tx& x, const DataType& h)
#endif
            {
              DataType rec_det_2_der_h_part(0);
              rec_det_2_der_h_part = -DataType(2)*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(3))*(-DataType(4)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h)-DataType(8)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h*h));
              return rec_det_2_der_h_part;
            }

          /**
           * \brief Computes the gradient of the Frobenius norm term
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          void add_grad_frobenius_part(Tx& grad_frobenius_part, const Tx& x, const DataType& h)
#else
            void NOINLINE add_grad_frobenius_part(Tx& grad_frobenius_part, const Tx& x, const DataType& h)
#endif
            {
              grad_frobenius_part(0,0) += DataType(8)/DataType(9)*this->_fac_frobenius*(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)))/(h*h*h*h)*(-DataType(4)*x(0,0)+DataType(2)*x(1,0)+DataType(2)*x(2,0));

              grad_frobenius_part(0,1) += DataType(8)/DataType(9)*this->_fac_frobenius*(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)))/(h*h*h*h)*(-DataType(4)*x(0,1)+DataType(2)*x(1,1)+DataType(2)*x(2,1));

              grad_frobenius_part(1,0) += DataType(8)/DataType(9)*this->_fac_frobenius*(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)))/(h*h*h*h)*(DataType(2)*x(0,0)-DataType(4)*x(1,0)+DataType(2)*x(2,0));

              grad_frobenius_part(1,1) += DataType(8)/DataType(9)*this->_fac_frobenius*(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)))/(h*h*h*h)*(DataType(2)*x(0,1)-DataType(4)*x(1,1)+DataType(2)*x(2,1));

              grad_frobenius_part(2,0) += DataType(8)/DataType(9)*this->_fac_frobenius*(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)))/(h*h*h*h)*(DataType(2)*x(0,0)+DataType(2)*x(1,0)-DataType(4)*x(2,0));

              grad_frobenius_part(2,1) += DataType(8)/DataType(9)*this->_fac_frobenius*(DataType(3)*h*h-DataType(2)*Math::sqr(x(0,0))+DataType(2)*x(0,0)*x(1,0)+DataType(2)*x(0,0)*x(2,0)-DataType(2)*Math::sqr(x(0,1))+DataType(2)*x(0,1)*x(1,1)+DataType(2)*x(0,1)*x(2,1)-DataType(2)*Math::sqr(x(1,0))+DataType(2)*x(1,0)*x(2,0)-DataType(2)*Math::sqr(x(1,1))+DataType(2)*x(1,1)*x(2,1)-DataType(2)*Math::sqr(x(2,0))-DataType(2)*Math::sqr(x(2,1)))/(h*h*h*h)*(DataType(2)*x(0,1)+DataType(2)*x(1,1)-DataType(4)*x(2,1));

            }


          /**
           * \brief Computes the gradient of the Frobenius norm term
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          void add_grad_det_1_part(Tx& grad_det_1_part, const Tx& x, const DataType& h)
#else
            void NOINLINE add_grad_det_1_part(Tx& grad_det_1_part, const Tx& x, const DataType& h)
#endif
            {
              grad_det_1_part(0,0) += DataType(2)/DataType(3)*this->_fac_det*Math::sqrt(DataType(3))*(x(1,1)-x(2,1))/Math::sqr(h);

              grad_det_1_part(0,1) += DataType(2)/DataType(3)*this->_fac_det*Math::sqrt(DataType(3))*(-x(1,0)+x(2,0))/Math::sqr(h);

              grad_det_1_part(1,0) += DataType(2)/DataType(3)*this->_fac_det*Math::sqrt(DataType(3))*(-x(0,1)+x(2,1))/Math::sqr(h);

              grad_det_1_part(1,1) += DataType(2)/DataType(3)*this->_fac_det*Math::sqrt(DataType(3))*(x(0,0)-x(2,0))/Math::sqr(h);

              grad_det_1_part(2,0) += DataType(2)/DataType(3)*this->_fac_det*Math::sqrt(DataType(3))*(x(0,1)-x(1,1))/Math::sqr(h);

              grad_det_1_part(2,1) += DataType(2)/DataType(3)*this->_fac_det*Math::sqrt(DataType(3))*(-x(0,0)+x(1,0))/Math::sqr(h);
            }

          /**
           * \brief Computes the gradient of the Frobenius norm term
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          void add_grad_det_2_part(Tx& grad_det_2_part, const Tx& x, const DataType& h)
#else
            void NOINLINE add_grad_det_2_part(Tx& grad_det_2_part, const Tx& x, const DataType& h)
#endif
            {
              grad_det_2_part(0,0) += DataType(8)/DataType(3)*this->_fac_det*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(1,1)-x(2,1));

              grad_det_2_part(0,1) += DataType(8)/DataType(3)*this->_fac_det*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(1,0)+x(2,0));

              grad_det_2_part(1,0) += DataType(8)/DataType(3)*this->_fac_det*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(0,1)+x(2,1));

              grad_det_2_part(1,1) += DataType(8)/DataType(3)*this->_fac_det*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(0,0)-x(2,0));

              grad_det_2_part(2,0) += DataType(8)/DataType(3)*this->_fac_det*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(0,1)-x(1,1));

              grad_det_2_part(2,1) += DataType(8)/DataType(3)*this->_fac_det*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(0,0)+x(1,0));
            }

          /**
           * \brief Computes the gradient of the Frobenius norm term
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          void add_grad_rec_det_1_part(Tx& grad_rec_det_1_part, const Tx& x, const DataType& h)
#else
            void NOINLINE add_grad_rec_det_1_part(Tx& grad_rec_det_1_part, const Tx& x, const DataType& h)
#endif
            {
              grad_rec_det_1_part(0,0) += -this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(2))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(1,1)-x(2,1))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(1,1)-x(2,1)));

              grad_rec_det_1_part(0,1) += -this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(2))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(-x(1,0)+x(2,0))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(1,0)+x(2,0)));

              grad_rec_det_1_part(1,0) += -this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(2))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(-x(0,1)+x(2,1))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(0,1)+x(2,1)));

              grad_rec_det_1_part(1,1) += -this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(2))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)-x(2,0))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(0,0)-x(2,0)));

              grad_rec_det_1_part(2,0) += -this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(2))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,1)-x(1,1))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(0,1)-x(1,1)));

              grad_rec_det_1_part(2,1) += -this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(2))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(-x(0,0)+x(1,0))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(0,0)+x(1,0)));
            }

          /**
           * \brief Computes the gradient of the Frobenius norm term
           */
          /// \compilerhack icc < 16.0 gets confused by NOINLINE
#if defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600)
          void add_grad_rec_det_2_part(Tx& grad_rec_det_2_part, const Tx& x, const DataType& h)
#else
            void NOINLINE add_grad_rec_det_2_part(Tx& grad_rec_det_2_part, const Tx& x, const DataType& h)
#endif
            {
              grad_rec_det_2_part(0,0) += -DataType(2)*this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(3))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(1,1)-x(2,1))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(1,1)-x(2,1)));

              grad_rec_det_2_part(0,1) += -DataType(2)*this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(3))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(-x(1,0)+x(2,0))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(1,0)+x(2,0)));

              grad_rec_det_2_part(1,0) += -DataType(2)*this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(3))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(-x(0,1)+x(2,1))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(0,1)+x(2,1)));

              grad_rec_det_2_part(1,1) += -DataType(2)*this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(3))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)-x(2,0))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(0,0)-x(2,0)));

              grad_rec_det_2_part(2,0) += -DataType(2)*this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(3))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,1)-x(1,1))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(x(0,1)-x(1,1)));

              grad_rec_det_2_part(2,1) += -DataType(2)*this->_fac_rec_det*Math::pow(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/Math::sqr(h)+Math::sqrt(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h))/DataType(3),-DataType(3))*(DataType(2)/DataType(3)*Math::sqrt(DataType(3))*(-x(0,0)+x(1,0))/Math::sqr(h)+DataType(4)*Math::pow(DataType(9)*this->_fac_reg*this->_fac_reg+DataType(12)*Math::pow(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0),DataType(2))/(h*h*h*h),-DataType(1)/DataType(2))*(x(0,0)*x(1,1)-x(0,0)*x(2,1)-x(1,0)*x(0,1)+x(2,0)*x(0,1)+x(1,0)*x(2,1)-x(1,1)*x(2,0))/(h*h*h*h)*(-x(0,0)+x(1,0)));
            }


      }; // class RumpfFunctionalUnrolled


    //#ifdef FEAT_EICKT
    //    extern template class RumpfFunctionalUnrolled<double,
    //    Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, double>>>;
    //#endif
    /// \endcond
  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONALS_2D_P1_UNROLLED_HPP
