#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONALS_P1_HPP
#define KERNEL_MESHOPT_RUMPF_FUNCTIONALS_P1_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/shape.hpp>
#include <kernel/meshopt/base.hpp>
#include <kernel/meshopt/rumpf_functional.hpp>

#include <kernel/eval_tags.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/cubature/rule.hpp>

namespace FEAT
{
  namespace Meshopt
  {
    /// \cond internal

    /**
     * \brief Class template for Rumpf functionals, standard P1 trafo on a conformal simplex mesh
     */
    template<typename DataType_, int shape_dim_>
    class RumpfFunctional<DataType_,
    Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<shape_dim_>, shape_dim_, shape_dim_, DataType_>>> :
      public RumpfFunctionalBase<DataType_>
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
          typedef Shape::Simplex<shape_dim_> ShapeType;
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

          // Data from evaluating the transformation
          typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;
          // Data from evaluating FE spaces
          typedef typename SpaceEvaluator::template ConfigTraits<space_config>::EvalDataType SpaceEvalData;


        private:
          /// The gradient of the local reference mapping of the normalised reference cell to the current cell
          TgradR grad_R;
          /// The inverse matrix of the local reference mapping of the normalised reference cell to the current cell
          TgradR inv_grad_R;
          /// The cofactor matrix of the local reference mapping of the normalised reference cell to the current cell
          TgradR cof_grad_R;

          /// Factory for creating the cubature rule
          Cubature::DynamicFactory _cubature_factory;
          /// The cubature rule for integrating
          Cubature::Rule<ShapeType, DataType, DataType, Tiny::Vector<DataType, world_dim>> _cubature_rule;

          /// \f$ \| \nabla R_K \|_F \f$
          DataType _frobenius_grad_R;
          /// \f$ \| \mathrm{Cof}(\nabla R_K) \|_F \f$
          DataType _frobenius_cof_grad_R;
          /// \f$ det(\nabla R_K) \f$
          DataType _det_grad_R;
          /// The volume of the normalised reference cell
          DataType _normalised_ref_cell_vol;

          /// Exponent of the det terms
          const int _exponent_det;

          /// Do we need to compute \f$ \| \nabla R_K \|_F \f$?
          const bool _compute_frobenius;
          /// Do we need to compute \f$ \| \mathrm{Cof}(\nabla R_K) \|_F \f$?
          const bool _compute_cof;
          /// Do we need to compute \f$ det(\nabla R_K) \f$?
          const bool _compute_det;
          /// Do we need to compute \f$ (\nabla R_K)^{-1} \f$?
          const bool _compute_inverse;

        public:
          /**
           * \brief Constructor
           */
          explicit RumpfFunctional(
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
            grad_R(0),
            inv_grad_R(0),
            _cubature_factory("auto-degree:0"),
            _cubature_rule(Cubature::ctor_factory,_cubature_factory),
            _frobenius_grad_R(0),
            _frobenius_cof_grad_R(0),
            _det_grad_R(0),
            _normalised_ref_cell_vol(1),
            _exponent_det(exponent_det_),
            _compute_frobenius( (fac_frobenius_ > DataType(0)) ),
            _compute_cof( (fac_cof_ > 0) ),
            _compute_det( (fac_det_ > 0) ),
            _compute_inverse( _compute_det || _compute_cof)
            {
              XASSERTM(exponent_det_ == 1 || exponent_det_ == 2,"exponent_det must be 1 or 2!");
              XASSERTM(world_dim == 3 || fac_cof_ == DataType(0),"In 2d, the cofactor and frobenius norm term are redundant, so set fac_cof = 0.");

              // vol(Simplex<d>) = d!
              for(int d(1); d < world_dim; ++d)
              {
                _normalised_ref_cell_vol /= DataType(d+1);
              }

            }

          /**
           * \brief The class name
           *
           * \returns String with the class name
           */
          static String name()
          {
            return "RumpfFunctional<"+ShapeType::name()+">";
          }

          /**
           * \brief Prints object parameters
           */
          void print()
          {
            Index pad_width(30);
            Dist::Comm comm_world(Dist::Comm::world());

            String msg;

            msg = name()+":";
            comm_world.print(msg);

            BaseClass::print();

            msg = String("exponent_det").pad_back(pad_width, '.') + String(": ") + stringify(_exponent_det);
            comm_world.print(msg);

            msg = String("cubature_rule").pad_back(pad_width, '.') += String(": ") += _cubature_rule.get_name();
            comm_world.print(msg);
          }

          void set_point(const TrafoEvalData& trafo_data, const TgradR& mat_tensor)
          {

            // We abuse inv_grad_R as temporary storage
            inv_grad_R = trafo_data.jac_mat*mat_tensor;

            // grad_R = mat_tensor^T * grad(trafo) because we need to differentiate in the scaled coordinate system
            grad_R.set_transpose(inv_grad_R);

            // Set the inverse matrix needed for everything related to the derivative of det(grad_R)
            if(_compute_inverse)
            {
              inv_grad_R.set_inverse(grad_R);
            }

            if(_compute_frobenius)
            {
              _frobenius_grad_R = grad_R.norm_frobenius();
            }

            if(_compute_det)
            {
              _det_grad_R = grad_R.det();
            }

            if(_compute_cof)
            {
              cof_grad_R.set_cofactor(grad_R);
              _frobenius_cof_grad_R = cof_grad_R.norm_frobenius();
            }

          }

          /**
           * \brief Computes the Frobenius norm term for one cell
           */
          DataType compute_frobenius_part()
          {
            return Math::sqr(Math::sqr(_frobenius_grad_R) - DataType(TgradR::n))/_normalised_ref_cell_vol;
          }

          /**
           * \brief Computes the cof term on one element
           */
          DataType compute_cof_part()
          {
            return Math::sqr(Math::sqr(_frobenius_cof_grad_R) - DataType(TgradR::n))/_normalised_ref_cell_vol;
          }

          /**
           * \brief Computes the det term on one element
           */
          DataType compute_det_part()
          {
            if(_exponent_det == 1)
            {
              return _det_grad_R/_normalised_ref_cell_vol;
            }
            else
            {
              return Math::sqr(_det_grad_R)/_normalised_ref_cell_vol;
            }
          }

          /**
           * \brief Computes the det term on one element
           */
          DataType compute_rec_det_part()
          {
            DataType fac;
            if(_det_grad_R <= DataType(0))
            {
              fac = DataType(1)/this->_fac_reg;
            }
            else
            {
              fac = DataType(1)/( _det_grad_R + Math::sqrt(Math::sqr(this->_fac_reg) + Math::sqr(_det_grad_R)) );
            }

            if(_exponent_det == 1)
            {
              return fac/_normalised_ref_cell_vol;
            }
            else
            {
              return Math::sqr(fac)/_normalised_ref_cell_vol;
            }
          }

          /**
           * \brief Computes the functional gradient for one cell
           */
          void eval_fval_grad(DataType& fval, Tx& grad, const TgradR& mat_tensor, const TrafoEvaluator& trafo_eval, const SpaceEvaluator& space_eval, const Tx& DOXY(x), const DataType& DOXY(h))
          {
            TrafoEvalData trafo_data;
            SpaceEvalData space_data;

            grad.format(DataType(0));
            fval = DataType(0);

            DataType fval_frobenius(0);
            DataType fval_cof(0);
            DataType fval_det(0);
            DataType fval_rec_det(0);

            for(int k(0); k < _cubature_rule.get_num_points(); ++k)
            {
              // Evaluate trafo and FE space
              trafo_eval(trafo_data, _cubature_rule.get_point(k));
              space_eval(space_data, trafo_data);
              // Compute quantities for nonlinear form
              set_point(trafo_data, mat_tensor);

              DataType weight(_cubature_rule.get_weight(k));

              if(_compute_frobenius)
              {
                fval_frobenius += this->_fac_frobenius*weight*compute_frobenius_part();
                // The factor jac_det cancels out with the chain rule for the derivative of the test function
                add_grad_frobenius(grad, space_data, mat_tensor, this->_fac_frobenius*weight);
              }

              if(_compute_cof)
              {
                fval_cof += this->_fac_cof*weight*compute_cof_part();
                // The factor jac_det cancels out with the chain rule for the derivative of the test function
                add_grad_cof(grad, space_data, mat_tensor, this->_fac_cof*weight);
              }

              if(_compute_det)
              {
                fval_det += this->_fac_det*weight*compute_det_part();
                add_grad_det(grad, space_data, mat_tensor, this->_fac_det*weight);

                fval_rec_det += this->_fac_rec_det*weight*compute_rec_det_part();
                add_grad_rec_det(grad, space_data, mat_tensor, this->_fac_rec_det*weight);
              }

            }

            fval = fval_frobenius + fval_cof + fval_det + fval_rec_det;

          }

          /**
           * \brief Computes the functional gradient for one cell
           */
          void eval_fval_cellwise(DataType& fval, const TgradR& mat_tensor, const TrafoEvaluator& trafo_eval, const SpaceEvaluator& space_eval, const Tx& DOXY(x), const DataType& DOXY(h), DataType& fval_frobenius, DataType& fval_cof, DataType& fval_det)
          {

            TrafoEvalData trafo_data;
            SpaceEvalData space_data;

            fval = DataType(0);
            fval_frobenius = DataType(0);
            fval_cof = DataType(0);
            fval_det = DataType(0);

            TgradR grad_rec_det(0);

            for(int k(0); k < _cubature_rule.get_num_points(); ++k)
            {
              // Evaluate trafo and FE space
              trafo_eval(trafo_data, _cubature_rule.get_point(k));
              space_eval(space_data, trafo_data);

              // Compute quantities for nonlinear form
              set_point(trafo_data, mat_tensor);

              DataType weight(_cubature_rule.get_weight(k));

              if(_compute_frobenius)
              {
                fval_frobenius += weight*this->_fac_frobenius*compute_frobenius_part();
              }

              if(_compute_cof)
              {
                fval_cof += weight*this->_fac_cof*compute_frobenius_part();
              }

              if(_compute_det)
              {
                fval_det += weight*this->_fac_det*compute_det_part();
                fval_det += weight*this->_fac_rec_det*compute_rec_det_part();
              }

            }

            fval = fval_frobenius + fval_cof + fval_det;

          }

          /**
           * \brief Computes the gradient of the Frobenius norm term
           */
          void add_grad_frobenius(Tx& grad, const SpaceEvalData& space_data, const TgradR& mat_tensor, const DataType fac)
          {
            DataType my_fac(fac*DataType(4)*(Math::sqr(_frobenius_grad_R) - DataType(TgradR::n)));
            my_fac /= _normalised_ref_cell_vol;

            for(int i(0); i < SpaceEvalData::max_local_dofs; ++i)
            {
              for(int d(0); d < world_dim; ++d)
              {
                grad[i] += my_fac*(mat_tensor[d]*space_data.phi[i].ref_grad[d])*grad_R;
              }
            }
          }

          /**
           * \brief Computes the gradient of the cofactor matrix frobenius norm term
           */
          void add_grad_cof(Tx& grad, const SpaceEvalData& space_data, const TgradR& mat_tensor, const DataType fac)
          {
            add_grad_cof1(grad, space_data, mat_tensor, fac);
            add_grad_cof2(grad, space_data, mat_tensor, fac);
          }

          /**
           * \brief Computes the gradient of the cofactor matrix frobenius norm term
           */
          void add_grad_cof1(Tx& grad, const SpaceEvalData& space_data, const TgradR& mat_tensor, const DataType fac)
          {
            DataType my_fac(fac*DataType(4)*(Math::sqr(_frobenius_cof_grad_R) - DataType(TgradR::n)));
            my_fac /= _normalised_ref_cell_vol;

            TgradR cof_grad_R_transpose(0);
            cof_grad_R_transpose.set_transpose(cof_grad_R);


            for(int i(0); i < SpaceEvalData::max_local_dofs; ++i)
            {
              auto transformed_phi_ref_grad = space_data.phi[i].ref_grad*mat_tensor;

              for(int d(0); d < world_dim; ++d)
              {
                grad[i][d] += my_fac*(
                  Tiny::dot( inv_grad_R[d], transformed_phi_ref_grad)*Math::sqr(_frobenius_cof_grad_R) );
                //  - Tiny::dot(cof_grad_R[d], cof_grad_R_transpose*(inv_grad_R*transformed_phi_ref_grad)));
              }
            }
          }

          /**
           * \brief Computes the gradient of the cofactor matrix frobenius norm term
           */
          void add_grad_cof2(Tx& grad, const SpaceEvalData& space_data, const TgradR& mat_tensor, const DataType fac)
          {
            DataType my_fac(fac*DataType(4)*(Math::sqr(_frobenius_cof_grad_R) - DataType(TgradR::n)));
            my_fac /= _normalised_ref_cell_vol;

            TgradR cof_grad_R_transpose(0);
            cof_grad_R_transpose.set_transpose(cof_grad_R);


            for(int i(0); i < SpaceEvalData::max_local_dofs; ++i)
            {
              auto transformed_phi_ref_grad = space_data.phi[i].ref_grad*mat_tensor;

              for(int d(0); d < world_dim; ++d)
              {
                grad[i][d] += my_fac*(
                  //Tiny::dot( inv_grad_R[d], transformed_phi_ref_grad)*Math::sqr(_frobenius_cof_grad_R) );
                  - Tiny::dot(cof_grad_R[d], cof_grad_R_transpose*(inv_grad_R*transformed_phi_ref_grad)));
              }
            }
          }

          /**
           * \brief Computes the gradient of the Frobenius norm term
           */
          void add_grad_det(Tx& grad, const SpaceEvalData& space_data, const TgradR& mat_tensor, const DataType fac)
          {
            DataType my_fac(fac);
            if(_exponent_det == 1)
            {
              my_fac *= _det_grad_R/_normalised_ref_cell_vol;
            }
            else
            {
              my_fac *= DataType(2)*Math::sqr(_det_grad_R)/_normalised_ref_cell_vol;
            }

            for(int i(0); i < SpaceEvalData::max_local_dofs; ++i)
            {
              for(int d(0); d < world_dim; ++d)
              {
                grad[i] += my_fac*inv_grad_R*(mat_tensor[d]*space_data.phi[i].ref_grad[d]);
              }
            }
          }

          /**
           * \brief Computes the gradient of the 1/det term
           */
          void add_grad_rec_det(Tx& grad, const SpaceEvalData& space_data, const TgradR& mat_tensor, const DataType fac)
          {
            DataType my_fac(-fac/_normalised_ref_cell_vol);
            if(_exponent_det == 1)
            {
              my_fac *= _det_grad_R / Math::sqrt(Math::sqr(this->_fac_reg) + Math::sqr(_det_grad_R));
              if(_det_grad_R > DataType(0))
              {
                my_fac /= (_det_grad_R + Math::sqrt(Math::sqr(this->_fac_reg) + Math::sqr(_det_grad_R)));
              }
              else
              {
                my_fac /= this->_fac_reg;
              }
            }
            else
            {
              my_fac *= DataType(2)*Math::sqr(_det_grad_R) / (Math::sqr(this->_fac_reg) + Math::sqr(_det_grad_R));
              if(_det_grad_R > DataType(0))
              {
                my_fac /= Math::sqr(_det_grad_R + Math::sqrt(Math::sqr(this->_fac_reg) + Math::sqr(_det_grad_R)));
              }
              else
              {
                my_fac /= Math::sqr(this->_fac_reg);
              }
            }

            for(int i(0); i < SpaceEvalData::max_local_dofs; ++i)
            {
              for(int d(0); d < world_dim; ++d)
              {
                grad[i] += my_fac*inv_grad_R*(mat_tensor[d]*space_data.phi[i].ref_grad[d]);
              }
            }
          }

          /**
           * \brief Adds the part coming from the chain rule involving h to the local gradient
           */
          void add_grad_h_part(Tx& grad, const TgradR& mat_tensor, const TrafoEvaluator& trafo_eval, const SpaceEvaluator& space_eval, const Tx& DOXY(x), const DataType& h, const Tgradh& grad_h)
          {

            TrafoEvalData trafo_data;
            SpaceEvalData space_data;

            DataType frobenius_der_h(0);
            DataType det_der_h(0);
            DataType rec_det_der_h(0);

            for(int k(0); k < _cubature_rule.get_num_points(); ++k)
            {
              // Evaluate trafo and FE space
              trafo_eval(trafo_data, _cubature_rule.get_point(k));
              space_eval(space_data, trafo_data);

              // Compute quantities for nonlinear form
              set_point(trafo_data, mat_tensor);

              DataType weight(_cubature_rule.get_weight(k));

              // Add outer parts of the chain rule derivatives, which are different in every integration point
              frobenius_der_h += weight*DataType(4)*(Math::sqr(_frobenius_grad_R) - DataType(world_dim))*
                Math::sqr(_frobenius_grad_R)/_normalised_ref_cell_vol;

              det_der_h += weight*_det_grad_R/_normalised_ref_cell_vol;

              DataType fac;
              if(_det_grad_R >= DataType(0))
              {
                fac = Math::sqrt(Math::sqr(this->_fac_reg) + Math::sqr(_det_grad_R));
              }
              else
              {
                fac = this->_fac_reg;
              }
              rec_det_der_h += weight*_det_grad_R/(fac*(_det_grad_R + fac))/_normalised_ref_cell_vol;

            }

            // The inner part of the chain rule derivative is constant wrt. to space, so we can just scale at this point
            frobenius_der_h *= this->_fac_frobenius*DataType(-1)/h;
            det_der_h *= this->_fac_det*DataType(-2)/h;
            rec_det_der_h *= this->_fac_rec_det*DataType(2)/h;

            DataType der_h(frobenius_der_h + det_der_h + rec_det_der_h);

            // Add the local contributions to the global gradient
            for(int i(0); i < Tx::m; ++i)
            {
              for(int d(0); d < Tx::n; ++d)
              {
                grad(i,d) += der_h*grad_h(i*Tx::n + d);
              }
            }
          } // add_grad_h_part

      }; // class RumpfFunctional
#ifdef FEAT_EICKT
    extern template class RumpfFunctional<double,
    Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, double>>>;

    extern template class RumpfFunctional<double,
    Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, double>>>;
#endif // FEAT_EICKT
    /// \endcond
  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONALS_P1_HPP
