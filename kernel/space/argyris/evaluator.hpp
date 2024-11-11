// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/space/evaluator_base.hpp>
#include <kernel/geometry/intern/sub_index_mapping.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Argyris
    {
      /**
       * \brief Argyris Element Evaluator class template declaration.
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_,
        typename Shape_ = typename Space_::ShapeType>
      class Evaluator DOXY({});

      /**
       * \brief Argyris Element Evaluator implementation for Simplex<2> shape
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Simplex<2> > :
        public EvaluatorBase<
          Evaluator<
            Space_,
            TrafoEvaluator_,
            SpaceEvalTraits_,
            Shape::Simplex<2> >,
          TrafoEvaluator_,
          SpaceEvalTraits_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Evaluator, TrafoEvaluator_, SpaceEvalTraits_> BaseClass;

        /// space type
        typedef Space_ SpaceType;

        /// space evaluation traits
        typedef SpaceEvalTraits_ SpaceEvalTraits;

        /// trafo evaluator type
        typedef TrafoEvaluator_ TrafoEvaluator;

        /// trafo evaluator traits
        typedef typename TrafoEvaluator::EvalTraits TrafoEvalTraits;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;

        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

        static constexpr SpaceTags eval_caps = (SpaceTags::value | SpaceTags::grad | SpaceTags::hess);

        template<SpaceTags cfg_>
        struct ConfigTraits
        {
          /// evaluation data configuration
          static constexpr SpaceTags config = cfg_;

          /// trafo configuration
          static constexpr TrafoTags trafo_config = TrafoTags::img_point;

          /// evaluation data typedef
          typedef Space::EvalData<SpaceEvalTraits, config> EvalDataType;
        };

      protected:
        // coefficient matrix
        Tiny::Matrix<DataType, 21, 21> _coeff;
        // triangle barycentre
        ImagePointType _barycentre;

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const SpaceType& DOXY(space))
        {
        }

        /**
         * \brief Returns the number of local DOFs.
         *
         * \returns
         * The number of local dofs.
         */
        int get_num_local_dofs() const
        {
          return 21;
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& trafo_eval)
        {
          // domain point
          DomainPointType dom_point;
          ImagePointType img_point;

          // create sub-index mapping
          typedef Geometry::Intern::SubIndexMapping<Shape::Simplex<2>, 1, 0> SimType;

          // compute edge orientations
          SimType sim(
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,0>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,1>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<1,0>());
          DataType edge_ori[3] =
          {
            DataType(1) - DataType(2*sim.map(0, 0)),
            DataType(1) - DataType(2*sim.map(1, 0)),
            DataType(1) - DataType(2*sim.map(2, 0))
          };

          // initialize edge vectors
          ImagePointType edge_vec[3];
          edge_vec[0].format();
          edge_vec[1].format();
          edge_vec[2].format();

          // compute barycentre
          dom_point[0] = dom_point[1] = DataType(1) / DataType(3);
          trafo_eval.map_point(_barycentre, dom_point);

          // nodal matrix
          Tiny::Matrix<DataType, 21, 21> node_mat;
          node_mat.format();

          // loop over all vertices
          for(int vi(0); vi < 3; ++vi)
          {
            // compute vertex coordinate
            dom_point[0] = (vi == 1) ? DataType(1) : DataType(0);
            dom_point[1] = (vi == 2) ? DataType(1) : DataType(0);
            trafo_eval.map_point(img_point, dom_point);
            img_point -= _barycentre;

            // update edge vectors
            edge_vec[(vi+1) % 3] += img_point;
            edge_vec[(vi+2) % 3] -= img_point;

            // initialize monomial powers
            Tiny::Vector<DataType, 6> vx, vy;
            vx(0) = vy(0) = DataType(1);
            for(int l(0); l < 5; ++l)
            {
              vx(l+1) = vx(l) * img_point[0];
              vy(l+1) = vy(l) * img_point[1];
            }

            // set monomial values
            int k(0);
            for(int i(0); i < 6; ++i)
            {
              for(int j(0); i+j < 6; ++j, ++k)
              {
                // monomial value
                node_mat(k, 6*vi+0) = vx(i) * vy(j);
                // dx-derivative
                if(i > 0)
                  node_mat(k, 6*vi+1) = DataType(i) * vx(i-1) * vy(j);
                // dy-derivative
                if(j > 0)
                  node_mat(k, 6*vi+2) = DataType(j) * vy(j-1) * vx(i);
                // dxx-derivative
                if(i > 1)
                  node_mat(k, 6*vi+3) = DataType(i) * DataType(i-1) * vx(i-2) * vy(j);
                // dyy-derivative
                if(j > 1)
                  node_mat(k, 6*vi+4) = DataType(j) * DataType(j-1) * vy(j-2) * vx(i);
                // dxy-derivative
                if(i*j > 0)
                  node_mat(k, 6*vi+5) = DataType(i) * vx(i-1) * DataType(j) * vy(j-1);
              }
            }
          }

          for(int ei(0); ei < 3; ++ei)
          {
            // compute edge midpoint
            dom_point[0] = (ei != int(1) ? DataType(0.5) : DataType(0));
            dom_point[1] = (ei != int(2) ? DataType(0.5) : DataType(0));
            trafo_eval.map_point(img_point, dom_point);
            img_point -= _barycentre;

            // compute edge normal and apply orientation
            DataType dnorm = edge_vec[ei].norm_euclid();
            DataType dnx = +(edge_vec[ei](1) * edge_ori[ei]) / dnorm;
            DataType dny = -(edge_vec[ei](0) * edge_ori[ei]) / dnorm;

            // initialize monomial powers
            Tiny::Vector<DataType, 6> vx, vy;
            vx(0) = vy(0) = DataType(1);
            for(int l(0); l < 5; ++l)
            {
              vx(l+1) = vx(l) * img_point[0];
              vy(l+1) = vy(l) * img_point[1];
            }

            // set monomial values
            int k(0);
            for(int i(0); i < 6; ++i)
            {
              for(int j(0); i+j < 6; ++j, ++k)
              {
                // edge normal derivative
                if(i > 0)
                  node_mat(k, 18+ei) += DataType(i) * dnx * vx(i-1) * vy(j);
                if(j > 0)
                  node_mat(k, 18+ei) += DataType(j) * dny * vy(j-1) * vx(i);
              }
            }
          }

          // invert nodal matrix to get the coefficient matrix
          _coeff.set_inverse(node_mat);
        }

        template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
        void eval_values(
          EvalData<SpaceEvalTraits, space_cfg_>& phi,
          const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& tau) const
        {
          // initialize monomial powers
          Tiny::Vector<DataType, 6> vx, vy;
          vx(0) = vy(0) = DataType(1);
          for(int l(0); l < 5; ++l)
          {
            vx(l+1) = vx(l) * (tau.img_point[0] - _barycentre[0]);
            vy(l+1) = vy(l) * (tau.img_point[1] - _barycentre[1]);
          }

          // compute basis function values
          for(int l(0); l < 21; ++l)
          {
            int k(0);
            DataType v(DataType(0));
            for(int i(0); i < 6; ++i)
              for(int j(0); i+j < 6; ++j, ++k)
                v += _coeff(l,k) * vx(i) * vy(j);
            phi.phi[l].value = v;
          }
        }

        template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
        void eval_gradients(
          EvalData<SpaceEvalTraits, space_cfg_>& phi,
          const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& tau) const
        {
          // initialize monomial powers
          Tiny::Vector<DataType, 6> vx, vy;
          vx(0) = vy(0) = DataType(1);
          for(int l(0); l < 5; ++l)
          {
            vx(l+1) = vx(l) * (tau.img_point[0] - _barycentre[0]);
            vy(l+1) = vy(l) * (tau.img_point[1] - _barycentre[1]);
          }

          // compute basis function gradients
          for(int l(0); l < 21; ++l)
          {
            int k(0);
            DataType dx(DataType(0)), dy(DataType(0));
            for(int i(0); i < 6; ++i)
            {
              for(int j(0); i+j < 6; ++j, ++k)
              {
                if(i > 0)
                  dx += _coeff(l,k) * DataType(i) * vx(i-1) * vy(j);
                if(j > 0)
                  dy += _coeff(l,k) * DataType(j) * vy(j-1) * vx(i);
              }
            }
            phi.phi[l].grad[0] = dx;
            phi.phi[l].grad[1] = dy;
          }
        }

        template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
        void eval_hessians(
          EvalData<SpaceEvalTraits, space_cfg_>& phi,
          const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& tau) const
        {
          // initialize monomial powers
          Tiny::Vector<DataType, 6> vx, vy;
          vx(0) = vy(0) = DataType(1);
          for(int l(0); l < 5; ++l)
          {
            vx(l+1) = vx(l) * (tau.img_point[0] - _barycentre[0]);
            vy(l+1) = vy(l) * (tau.img_point[1] - _barycentre[1]);
          }

          // compute basis function hessians
          for(int l(0); l < 21; ++l)
          {
            int k(0);
            DataType dxx(DataType(0)), dyy(DataType(0)), dxy(DataType(0));
            for(int i(0); i < 6; ++i)
            {
              for(int j(0); i+j < 6; ++j, ++k)
              {
                if(i > 1)
                  dxx += _coeff(l,k) * DataType(i) * DataType(i-1) * vx(i-2) * vy(j);
                if(j > 1)
                  dyy += _coeff(l,k) * DataType(j) * DataType(j-1) * vy(j-2) * vx(i);
                if(i*j > 0)
                  dxy += _coeff(l,k) * DataType(i) * vx(i-1) * DataType(j) * vy(j-1);
              }
            }
            phi.phi[l].hess(0,0) = dxx;
            phi.phi[l].hess(1,1) = dyy;
            phi.phi[l].hess(0,1) = phi.phi[l].hess(1,0) = dxy;
          }
        }
      }; // class Evaluator<...,Simplex<2>>
    } // namespace Argyris
  } // namespace Space
} // namespace FEAT
