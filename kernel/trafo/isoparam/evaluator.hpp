// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_ISOPARAM_EVALUATOR_HPP
#define KERNEL_TRAFO_ISOPARAM_EVALUATOR_HPP 1

#include <kernel/trafo/standard/evaluator.hpp>
#include <kernel/geometry/atlas/chart.hpp>

namespace FEAT
{
  namespace Trafo
  {
    namespace Isoparam
    {
      template<
        typename Trafo_,
        typename EvalPolicy_,
        int degree_,
        typename Shape_ = typename EvalPolicy_::ShapeType>
      class Evaluator DOXY({});


      /// \cond internal
      namespace Intern
      {
        template<int n_>
        struct Basis
        {
          // interpolation points
          template<typename T_>
          static inline void points(T_ v[])
          {
            static const T_ s = T_(2) / T_(n_);
            for(int i(0); i <= n_; ++i)
              v[i] *= (T_(i)*s - T_(1));
          }

          // basis function values for generic P_n
          template<typename T_>
          static inline void val(T_ v[], const T_ x)
          {
            static const T_ s = T_(2) / T_(n_);

            for(int j(0); j <= n_; ++j)
            {
              v[j] = T_(1);
              for(int i(0); i <= n_; ++i)
                v[j] *= (i != j ? (T_(1) + x - T_(i)*s) / (T_(j)*s - T_(i)*s) : T_(1));
            }
          }

          // first order derivatives
          template<typename T_>
          static inline void d1(T_ v[], const T_ x)
          {
            static const T_ s = T_(2) / T_(n_);

            for(int j(0); j <= n_; ++j)
            {
              v[j] = T_(0);
              for(int k(0); k <= n_; ++k)
              {
                if(k == j)
                  continue;
                T_ t = T_(1);
                for(int i(0); i <= n_; ++i)
                  t *= ((i != k) && (i != j) ? (T_(1) + x - T_(i)*s) / (T_(j)*s - T_(i)*s) : T_(1));
                v[j] += t / (T_(j)*s - T_(k)*s);
              }
            }
          }

          // second order derivatives
          template<typename T_>
          static inline void d2(T_ v[], const T_ x)
          {
            static const T_ s = T_(2) / T_(n_);

            for(int j(0); j <= n_; ++j)
            {
              v[j] = T_(0);
              for(int k(0); k <= n_; ++k)
              {
                if(k == j)
                  continue;
                T_ t = T_(0);
                for(int i(0); i <= n_; ++i)
                {
                  if((i == k) || (i == j))
                    continue;
                  T_ r = T_(1);
                  for(int l(0); l <= n_; ++l)
                    r *= ((l != k) && (l != j) && (l != i) ? (T_(1) + x - T_(l)*s) / (T_(j)*s - T_(l)*s) : T_(1));
                  t += r / (T_(j)*s - T_(i)*s);
                }
                v[j] += t / (T_(j)*s - T_(k)*s);
              }
            }
          }
        }; // struct Basis<n_>

        template<>
        struct Basis<1>
        {
          // interpolation points for P1
          template<typename T_>
          static inline void points(T_ v[])
          {
            v[0] = -T_(1);
            v[1] = +T_(1);
          }

          // basis function values for P1
          template<typename T_>
          static inline void val(T_ v[], const T_ x)
          {
            v[0] = T_(0.5) * (T_(1) - x);
            v[1] = T_(0.5) * (T_(1) + x);
          }

          // first order derivatives
          template<typename T_>
          static inline void d1(T_ v[], const T_)
          {
            v[0] = -T_(0.5);
            v[1] = +T_(0.5);
          }

          // second order derivatives
          template<typename T_>
          static inline void d2(T_ v[], const T_)
          {
            v[0] = T_(0);
            v[1] = T_(0);
          }
        }; // Basis<1>

        template<>
        struct Basis<2>
        {
          // interpolation points for P2
          template<typename T_>
          static inline void points(T_ v[])
          {
            v[0] = -T_(1);
            v[1] =  T_(0);
            v[2] = +T_(1);
          }

          // basis function values for P2
          template<typename T_>
          static inline void val(T_ v[], const T_ x)
          {
            v[0] = T_(0.5) * x * (x - T_(1));
            v[1] = (T_(1) - x) * (T_(1) + x);
            v[2] = T_(0.5) * x * (x + T_(1));
          }

          // first order derivatives
          template<typename T_>
          static inline void d1(T_ v[], const T_ x)
          {
            v[0] = x - T_(0.5);
            v[1] = -T_(2) * x;
            v[2] = x + T_(0.5);
          }

          // second order derivatives
          template<typename T_>
          static inline void d2(T_ v[], const T_)
          {
            v[0] = T_(1);
            v[1] = -T_(2);
            v[2] = T_(1);
          }
        }; // Basis<2>

        template<>
        struct Basis<3>
        {
          // interpolation points for P3
          template<typename T_>
          static inline void points(T_ v[])
          {
            v[0] = -T_(1);
            v[1] = -T_(1)/T_(3);
            v[2] = +T_(1)/T_(3);
            v[3] = +T_(1);
          }

          // basis function values for P3
          template<typename T_>
          static inline void val(T_ v[], const T_ x)
          {
            v[0] = T_(0.0625) * (-T_(1) + x*( T_(1) + T_(9)*x*(T_(1) - x)));
            v[1] = T_(0.5625) * (T_(1) + x*(-T_(3) + x*(-T_(1) + T_(3)*x)));
            v[2] = T_(0.5625) * (T_(1) + x*( T_(3) + x*(-T_(1) - T_(3)*x)));
            v[3] = T_(0.0625) * (-T_(1) + x*(-T_(1) + T_(9)*x*(T_(1) + x)));
          }

          // first order derivatives
          template<typename T_>
          static inline void d1(T_ v[], const T_ x)
          {
            v[0] = T_(0.0625) * ( T_(1) + x*(T_(18) - T_(27)*x));
            v[1] = T_(0.0625) * (-T_(27) + x*(-T_(18) + T_(81)*x));
            v[2] = T_(0.0625) * ( T_(27) + x*(-T_(18) - T_(81)*x));
            v[3] = T_(0.0625) * (-T_(1) + x*(T_(18) + T_(27)*x));
          }

          // second order derivatives
          template<typename T_>
          static inline void d2(T_ v[], const T_ x)
          {
            v[0] = T_(1.125) * (T_(1) - T_(3)*x);
            v[1] = T_(1.125) * (-T_(1) + T_(9)*x);
            v[2] = T_(1.125) * (-T_(1) - T_(9)*x);
            v[3] = T_(1.125) * (T_(1) + T_(3)*x);
          }
        }; // Basis<3>

        /*
        // experimental: use Gauss-Quadrature points
        template<>
        struct Basis<3>
        {
          // interpolation points for P3
          template<typename T_>
          static inline void points(T_ v[])
          {
            v[0] = -T_(1);
            v[1] = -(v[2] = Math::sqrt(T_(1)/T_(3)));
            v[3] = +T_(1);
          }

          // basis function values for P3
          template<typename T_>
          static inline void val(T_ v[], const T_ x)
          {
            static const T_ S3 = Math::sqrt(T_(3));
            v[0] = T_(0.25)*(-T_(1) + x*( T_(1) + T_(3)*x*(T_(1) - x)));
            v[1] = T_(0.75)*( T_(1) + x*(-S3 + x*(-T_(1) + S3*x)));
            v[2] = T_(0.75)*( T_(1) + x*( S3 + x*(-T_(1) - S3*x)));
            v[3] = T_(0.25)*(-T_(1) + x*(-T_(1) + T_(3)*x*(T_(1) + x)));
          }

          // first order derivatives
          template<typename T_>
          static inline void d1(T_ v[], const T_ x)
          {
            static const T_ S3 = Math::sqrt(T_(3));
            static const T_ S27 = Math::sqrt(T_(27));
            v[0] = T_(0.25)*( T_(1) + x*(T_(6) - T_(9)*x));
            v[1] = T_(0.75)*(-S3 + x*(-T_(2) + S27*x));
            v[2] = T_(0.75)*( S3 + x*(-T_(2) - S27*x));
            v[3] = T_(0.25)*(-T_(1) + x*(T_(6) + T_(9)*x));
          }

          // second order derivatives
          template<typename T_>
          static inline void d2(T_ v[], const T_ x)
          {
            static const T_ S27 = Math::sqrt(T_(27));
            v[0] = T_(1.5)*( T_(1) - T_(3)*x);
            v[1] = T_(1.5)*(-T_(1) + sqrt(27)*x);
            v[2] = T_(1.5)*(-T_(1) - sqrt(27)*x);
            v[3] = T_(1.5)*(-T_(1) + T_(3)*x);
          }
        }; // Basis<3>
        */
      } // namespace Intern
      /// \endcond

      /* ************************************************************************************* */

      /**
       * \brief Specialization of iso-parametric trafo evaluator for Vertex shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_,
        int degree_>
      class Evaluator<Trafo_, EvalPolicy_, degree_, Shape::Vertex > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, degree_, Shape::Vertex >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Vertex ShapeType;
        /// trafo type using this evaluator
        typedef Trafo_ TrafoType;
        /// trafo evaluation traits
        typedef EvalPolicy_ EvalPolicy;

        /// type of the underlying mesh
        typedef typename TrafoType::MeshType MeshType;

        /// evaluation data type
        typedef typename EvalPolicy::DataType DataType;
        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;
        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;
        /// jacobian matrix type
        typedef typename EvalPolicy::JacobianMatrixType JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef typename EvalPolicy::JacobianInverseType JacobianInverseType;
        /// jacobian determinant type
        typedef typename EvalPolicy::JacobianDeterminantType JacobianDeterminantType;
        /// hessian tensor type
        typedef typename EvalPolicy::HessianTensorType HessianTensorType;
        /// hessian inverse tensor type
        typedef typename EvalPolicy::HessianInverseType HessianInverseType;

        /// domain dimension
        static constexpr int domain_dim = EvalPolicy::domain_dim;
        /// image dimension
        static constexpr int image_dim = EvalPolicy::image_dim;

        /// we can compute domain and image points
        static constexpr TrafoTags eval_caps = TrafoTags::dom_point | TrafoTags::img_point;

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim];

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] trafo
         * A reference to the trafo using this evaluator.
         */
        explicit Evaluator(const TrafoType& trafo) :
          BaseClass(trafo)
        {
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] cell_index
         * The index of the cell for which the evaluator is to be prepared.
         */
        void prepare(Index cell_index)
        {
          // prepare base-class
          BaseClass::prepare(cell_index);

          // fetch the mesh from the trafo
          const MeshType& mesh = this->_trafo.get_mesh();

          // fetch the vertex set from the mesh
          typedef typename MeshType::VertexSetType VertexSetType;
          const VertexSetType& vertex_set = mesh.get_vertex_set();

          // fetch the vertex
          typedef typename VertexSetType::VertexType VertexType;
          const VertexType& vtx = vertex_set[cell_index];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
          {
            _coeff[i] = DataType(vtx[i]);
          }
        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         */
        void map_point(ImagePointType& img_point, const DomainPointType& DOXY(dom_point)) const
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = _coeff[i];
          }
        }
      }; // class Evaluator<Vertex,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialization of iso-parametric trafo evaluator for Hypercube<1> shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_,
        int degree_>
      class Evaluator<Trafo_, EvalPolicy_, degree_, Shape::Hypercube<1> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, degree_, Shape::Hypercube<1> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Hypercube<1> ShapeType;
        /// trafo type using this evaluator
        typedef Trafo_ TrafoType;
        /// trafo evaluation traits
        typedef EvalPolicy_ EvalPolicy;

        /// type of the underlying mesh
        typedef typename TrafoType::MeshType MeshType;
        /// type of the chart
        typedef Geometry::Atlas::ChartBase<MeshType> ChartType;

        /// evaluation data type
        typedef typename EvalPolicy::DataType DataType;
        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;
        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;
        /// jacobian matrix type
        typedef typename EvalPolicy::JacobianMatrixType JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef typename EvalPolicy::JacobianInverseType JacobianInverseType;
        /// jacobian determinant type
        typedef typename EvalPolicy::JacobianDeterminantType JacobianDeterminantType;
        /// hessian tensor type
        typedef typename EvalPolicy::HessianTensorType HessianTensorType;
        /// hessian inverse tensor type
        typedef typename EvalPolicy::HessianInverseType HessianInverseType;

        /// domain dimension
        static constexpr int domain_dim = EvalPolicy::domain_dim;
        /// image dimension
        static constexpr int image_dim = EvalPolicy::image_dim;

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

      protected:
        /// the chart pointer array
        const std::vector<const ChartType*>& _charts;
        /// the coefficients of the trafo
        Tiny::Vector<DataType, image_dim> _iso_coeff[degree_+1];

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] trafo
         * A reference to the trafo using this evaluator.
         */
        explicit Evaluator(const TrafoType& trafo) :
          BaseClass(trafo),
          _charts(trafo.get_charts_vector(1))
        {
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] cell_index
         * The index of the cell for which the evaluator is to be prepared.
         */
        void prepare(Index cell_index)
        {
          // prepare base-class
          BaseClass::prepare(cell_index);

          // fetch the mesh from the trafo
          const MeshType& mesh = this->_trafo.get_mesh();

          // fetch the vertex set from the mesh
          typedef typename MeshType::VertexSetType VertexSetType;
          const VertexSetType& vertex_set = mesh.get_vertex_set();

          // fetch the index set
          typedef typename MeshType::template IndexSet<domain_dim, 0>::Type IndexSetType;
          const IndexSetType& index_set = mesh.template get_index_set<domain_dim, 0>();

          // for shorter indexing
          static constexpr int n = degree_;

          // store vertex points
          _iso_coeff[0] = vertex_set[index_set(cell_index, 0)];
          _iso_coeff[n] = vertex_set[index_set(cell_index, 1)];

          // interpolate inner edge points linearly
          //DataType v[degree_+1];
          //Intern::Basis<degree_>::points(v);
          for(int i(1); i < n; ++i)
            _iso_coeff[i].set_convex(DataType(i) / DataType(n), _iso_coeff[0], _iso_coeff[n]);
            //_iso_coeff[i].set_convex((v[i]+DataType(1))*DataType(0.5), _iso_coeff[0], _iso_coeff[n]);

          // get chart for this edge
          const ChartType* chart = _charts.at(cell_index);
          if(chart != nullptr)
          {
            // loop over all inner edge points and project them
            for(int i(1); i < degree_; ++i)
            {
              _iso_coeff[i] = chart->project(_iso_coeff[i]);
            }
          }
        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         */
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          DataType v[degree_+1];
          Intern::Basis<degree_>::val(v, dom_point[0]);

          img_point.format();
          for(int k(0); k <= degree_; ++k)
            img_point.axpy(v[k], _iso_coeff[k]);
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         */
        void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& dom_point) const
        {
          DataType v[degree_+1];

          Intern::Basis<degree_>::d1(v, dom_point[0]);

          jac_mat.format();
          for(int k(0); k < image_dim; ++k)
          {
            for(int i(0); i <= degree_; ++i)
              jac_mat(k,0) += v[i] * _iso_coeff[i][k];
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         */
        void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& dom_point) const
        {
          DataType v[degree_+1];

          Intern::Basis<degree_>::d2(v, dom_point[0]);

          hess_ten.format();
          for(int i(0); i < image_dim; ++i)
          {
            for(int k(0); k <= degree_; ++k)
              hess_ten(i,0,0) += v[k] * _iso_coeff[k][i];
          }
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          XASSERTM(false, "volume computation not available for isoparametric trafo");
          return DataType(0);
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        DataType width_directed(const ImagePointType&) const
        {
          XASSERTM(false, "cell width computation not available for isoparametric trafo");
          return volume();
        }
      }; // class Evaluator<Hypercube<1>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialization of iso-parametric trafo evaluator for Hypercube<2> shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_,
        int degree_>
      class Evaluator<Trafo_, EvalPolicy_, degree_, Shape::Hypercube<2> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, degree_, Shape::Hypercube<2> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Hypercube<2> ShapeType;
        /// trafo type using this evaluator
        typedef Trafo_ TrafoType;
        /// trafo evaluation traits
        typedef EvalPolicy_ EvalPolicy;

        /// type of the underlying mesh
        typedef typename TrafoType::MeshType MeshType;

        /// type of the chart
        typedef Geometry::Atlas::ChartBase<MeshType> ChartType;

        /// evaluation data type
        typedef typename EvalPolicy::DataType DataType;
        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;
        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;
        /// jacobian matrix type
        typedef typename EvalPolicy::JacobianMatrixType JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef typename EvalPolicy::JacobianInverseType JacobianInverseType;
        /// jacobian determinant type
        typedef typename EvalPolicy::JacobianDeterminantType JacobianDeterminantType;
        /// hessian tensor type
        typedef typename EvalPolicy::HessianTensorType HessianTensorType;
        /// hessian inverse tensor type
        typedef typename EvalPolicy::HessianInverseType HessianInverseType;

        /// domain dimension
        static constexpr int domain_dim = EvalPolicy::domain_dim;
        /// image dimension
        static constexpr int image_dim = EvalPolicy::image_dim;

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

      protected:
        /// the chart pointer arrays
        const std::vector<const ChartType*>& _charts_1;
        const std::vector<const ChartType*>& _charts_2;
        /// the coefficients of the trafo
        Tiny::Vector<DataType, image_dim> _iso_coeff[degree_+1][degree_+1];

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] trafo
         * A reference to the trafo using this evaluator.
         */
        explicit Evaluator(const TrafoType& trafo) :
          BaseClass(trafo),
          _charts_1(trafo.get_charts_vector(1)),
          _charts_2(trafo.get_charts_vector(2))
        {
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] cell_index
         * The index of the cell for which the evaluator is to be prepared.
         */
        void prepare(Index cell_index)
        {
          // prepare base-class
          BaseClass::prepare(cell_index);

          // fetch the mesh from the trafo
          const MeshType& mesh = this->_trafo.get_mesh();

          // fetch the vertex set from the mesh
          typedef typename MeshType::VertexSetType VertexSetType;
          const VertexSetType& vertex_set = mesh.get_vertex_set();

          // fetch the vertex-index set
          typedef typename MeshType::template IndexSet<2, 0>::Type IndexSetType;
          const IndexSetType& index_set = mesh.template get_index_set<2, 0>();

          // fetch the edge-index set
          typedef typename MeshType::template IndexSet<2, 1>::Type EdgeIndexSetType;
          const EdgeIndexSetType& edges_at_quad = mesh.template get_index_set<2, 1>();

          // for shorter indexing
          static constexpr int n = degree_;

          // store vertex points
          _iso_coeff[0][0] = vertex_set[index_set(cell_index, 0)];
          _iso_coeff[0][n] = vertex_set[index_set(cell_index, 1)];
          _iso_coeff[n][0] = vertex_set[index_set(cell_index, 2)];
          _iso_coeff[n][n] = vertex_set[index_set(cell_index, 3)];

          // vN0---vNN
          //  |     |    ^ y
          //  |     |    |
          // v00---v0N   +--> x

          // interpolate inner edge points linearly
          //DataType v[degree_+1];
          //Intern::Basis<degree_>::points(v);
          for(int i(1); i < degree_; ++i)
          {
            const DataType alpha = DataType(i) / DataType(degree_);
            //const DataType alpha = (v[i] + DataType(1)) * DataType(0.5);

            // bottom edge
            _iso_coeff[0][i].set_convex(alpha, _iso_coeff[0][0], _iso_coeff[0][n]);
            // top edge
            _iso_coeff[n][i].set_convex(alpha, _iso_coeff[n][0], _iso_coeff[n][n]);
            // left edge
            _iso_coeff[i][0].set_convex(alpha, _iso_coeff[0][0], _iso_coeff[n][0]);
            // right edge
            _iso_coeff[i][n].set_convex(alpha, _iso_coeff[0][n], _iso_coeff[n][n]);
          }

          // iso-coeff array index offsets and increments
          const int ox[4] = {1, 1, 0, n};
          const int oy[4] = {0, n, 1, 1};
          const int ix[4] = {1, 1, 0, 0};
          const int iy[4] = {0, 0, 1, 1};

          // loop over all four edges of the element
          for(int e = 0; e < 4; ++e)
          {
            // get global edge index
            const Index edge_idx = edges_at_quad(cell_index, e);

            // get chart for that edge
            const ChartType* chart = _charts_1.at(edge_idx);

            // no chart?
            if(chart == nullptr)
              continue;

            // initialize indices for this edge
            int y = oy[e];
            int x = ox[e];

            // loop over all edge points and project them
            for(int k(1); k < degree_; ++k, y += iy[e], x += ix[e])
            {
              _iso_coeff[y][x] = chart->project(_iso_coeff[y][x]);
            }
          }

          Tiny::Vector<DataType, image_dim> p_horz, p_vert;

          // interpolate inner quad points bilinearly
          for(int i(1); i < degree_; ++i)
          {
            const DataType alpha_v = DataType(i) / DataType(degree_);
            //const DataType alpha_v = (v[i] + DataType(1)) * DataType(0.5);

            for(int j(1); j < degree_; ++j)
            {
              const DataType alpha_h = DataType(j) / DataType(degree_);
              //const DataType alpha_h = (v[j] + DataType(1)) * DataType(0.5);

              // interpolate horizontally (between left and right edge)
              p_horz.set_convex(alpha_h, _iso_coeff[i][0], _iso_coeff[i][n]);

              // interpolate vertically between bottom and top edge
              p_vert.set_convex(alpha_v, _iso_coeff[0][j], _iso_coeff[n][j]);

              // combine both to obtain iso point
              _iso_coeff[i][j].set_convex(DataType(0.5), p_horz, p_vert);
            }
          }

          // get chart for this quad
          const ChartType* chart = _charts_2.at(cell_index);
          if(chart != nullptr)
          {
            // loop over all inner quad points and project them
            for(int i(1); i < degree_; ++i)
            {
              for(int j(1); j < degree_; ++j)
              {
                _iso_coeff[i][j] = chart->project(_iso_coeff[i][j]);
              }
            }
          }
        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         */
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          DataType vx[degree_+1], vy[degree_+1];
          Intern::Basis<degree_>::val(vx, dom_point[0]);
          Intern::Basis<degree_>::val(vy, dom_point[1]);

          img_point.format();
          for(int i(0); i <= degree_; ++i)
          {
            for(int j(0); j <= degree_; ++j)
              img_point.axpy(vy[i]*vx[j], _iso_coeff[i][j]);
          }
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         */
        void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& dom_point) const
        {
          DataType vx[degree_+1], vy[degree_+1];
          DataType dx[degree_+1], dy[degree_+1];

          Intern::Basis<degree_>::val(vx, dom_point[0]);
          Intern::Basis<degree_>::val(vy, dom_point[1]);
          Intern::Basis<degree_>::d1(dx, dom_point[0]);
          Intern::Basis<degree_>::d1(dy, dom_point[1]);

          jac_mat.format();
          for(int k(0); k < image_dim; ++k)
          {
            for(int i(0); i <= degree_; ++i)
            {
              for(int j(0); j <= degree_; ++j)
              {
                jac_mat(k,0) += vy[i] * dx[j] * _iso_coeff[i][j][k];
                jac_mat(k,1) += dy[i] * vx[j] * _iso_coeff[i][j][k];
              }
            }
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         */
        void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& dom_point) const
        {
          DataType vx[degree_+1], vy[degree_+1];
          DataType dx[degree_+1], dy[degree_+1];
          DataType qx[degree_+1], qy[degree_+1];

          Intern::Basis<degree_>::val(vx, dom_point[0]);
          Intern::Basis<degree_>::val(vy, dom_point[1]);
          Intern::Basis<degree_>::d1(dx, dom_point[0]);
          Intern::Basis<degree_>::d1(dy, dom_point[1]);
          Intern::Basis<degree_>::d2(qx, dom_point[0]);
          Intern::Basis<degree_>::d2(qy, dom_point[1]);

          hess_ten.format();
          for(int k(0); k < image_dim; ++k)
          {
            for(int i(0); i <= degree_; ++i)
            {
              for(int j(0); j <= degree_; ++j)
              {
                hess_ten(k,0,0) += vy[i] * qx[j] * _iso_coeff[i][j][k];
                hess_ten(k,0,1) += dy[i] * dx[j] * _iso_coeff[i][j][k];
                hess_ten(k,1,1) += qy[i] * vx[j] * _iso_coeff[i][j][k];
              }
            }
            hess_ten(k,1,0) = hess_ten(k,0,1);
          }
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          XASSERTM(false, "volume computation not available for isoparametric trafo");
          return DataType(0);
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        DataType width_directed(const ImagePointType&) const
        {
          XASSERTM(false, "cell width computation not available for isoparametric trafo");
          return DataType(0);
        }
      }; // class Evaluator<Hypercube<2>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialization of iso-parametric trafo evaluator for Hypercube<3> shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_,
        int degree_>
      class Evaluator<Trafo_, EvalPolicy_, degree_, Shape::Hypercube<3> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, degree_, Shape::Hypercube<3> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Hypercube<3> ShapeType;
        /// trafo type using this evaluator
        typedef Trafo_ TrafoType;
        /// trafo evaluation traits
        typedef EvalPolicy_ EvalPolicy;

        /// type of the underlying mesh
        typedef typename TrafoType::MeshType MeshType;

        /// type of the chart
        typedef Geometry::Atlas::ChartBase<MeshType> ChartType;

        /// evaluation data type
        typedef typename EvalPolicy::DataType DataType;
        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;
        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;
        /// jacobian matrix type
        typedef typename EvalPolicy::JacobianMatrixType JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef typename EvalPolicy::JacobianInverseType JacobianInverseType;
        /// jacobian determinant type
        typedef typename EvalPolicy::JacobianDeterminantType JacobianDeterminantType;
        /// hessian tensor type
        typedef typename EvalPolicy::HessianTensorType HessianTensorType;
        /// hessian inverse tensor type
        typedef typename EvalPolicy::HessianInverseType HessianInverseType;

        /// domain dimension
        static constexpr int domain_dim = EvalPolicy::domain_dim;
        /// image dimension
        static constexpr int image_dim = EvalPolicy::image_dim;

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

      protected:
        /// the chart pointer arrays
        const std::vector<const ChartType*>& _charts_1;
        const std::vector<const ChartType*>& _charts_2;
        const std::vector<const ChartType*>& _charts_3;
        /// the coefficients of the trafo
        Tiny::Vector<DataType, image_dim> _iso_coeff[degree_+1][degree_+1][degree_+1];

        // auxiliary helper functions for convex combinations:

        // q <- u1 + x*(u2-u1)
        template<typename T_, typename C_>
        static void _c1(C_& q, T_ x, const C_& u1, const C_& u2)
        {
          for(int i(0); i < C_::n; ++i)
            q[i] = u1[i] + x * (u2[i] - u1[i]);
        }

        // q <- 1/2 * (u1 + v1 + x*(u2-u1) + y*(v2-v1))
        template<typename T_, typename C_>
        static void _c2(C_& q, T_ x, const C_& u1, const C_& u2, T_ y, const C_& v1, const C_& v2)
        {
          for(int i(0); i < C_::n; ++i)
            q[i] = T_(0.5) * (u1[i] + v1[i] + x * (u2[i] - u1[i]) + y * (v2[i] - v1[i]));
        }

        // q <- 1/3 * (u1 + v1 + w1 + x*(u2-u1) + y*(v2-v1) + z*(w2-w1))
        template<typename T_, typename C_>
        static void _c3(C_& q, T_ x, const C_& u1, const C_& u2, T_ y, const C_& v1, const C_& v2, T_ z, const C_& w1, const C_& w2)
        {
          for(int i(0); i < C_::n; ++i)
            q[i] = (u1[i] + v1[i] + w1[i] + x * (u2[i] - u1[i]) + y * (v2[i] - v1[i]) + z * (w2[i] - w1[i])) / T_(3);
        }

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] trafo
         * A reference to the trafo using this evaluator.
         */
        explicit Evaluator(const TrafoType& trafo) :
          BaseClass(trafo),
          _charts_1(trafo.get_charts_vector(1)),
          _charts_2(trafo.get_charts_vector(2)),
          _charts_3(trafo.get_charts_vector(3))
        {
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] cell_index
         * The index of the cell for which the evaluator is to be prepared.
         */
        void prepare(Index cell_index)
        {
          // prepare base-class
          BaseClass::prepare(cell_index);

          // fetch the mesh from the trafo
          const MeshType& mesh = this->_trafo.get_mesh();

          // fetch the vertex set from the mesh
          typedef typename MeshType::VertexSetType VertexSetType;
          const VertexSetType& vertex_set = mesh.get_vertex_set();

          // fetch the vertex-index set
          typedef typename MeshType::template IndexSet<3, 0>::Type IndexSetType;
          const IndexSetType& index_set = mesh.template get_index_set<3, 0>();

          // fetch the edge-index set
          typedef typename MeshType::template IndexSet<3, 1>::Type EdgeIndexSetType;
          const EdgeIndexSetType& edges_at_hexa = mesh.template get_index_set<3, 1>();

          // fetch the quad-index set
          typedef typename MeshType::template IndexSet<3, 2>::Type QuadIndexSetType;
          const QuadIndexSetType& quads_at_hexa = mesh.template get_index_set<3, 2>();

          // for shorter indexing
          static constexpr int n = degree_;

          // store vertex points
          _iso_coeff[0][0][0] = vertex_set[index_set(cell_index, 0)];
          _iso_coeff[0][0][n] = vertex_set[index_set(cell_index, 1)];
          _iso_coeff[0][n][0] = vertex_set[index_set(cell_index, 2)];
          _iso_coeff[0][n][n] = vertex_set[index_set(cell_index, 3)];
          _iso_coeff[n][0][0] = vertex_set[index_set(cell_index, 4)];
          _iso_coeff[n][0][n] = vertex_set[index_set(cell_index, 5)];
          _iso_coeff[n][n][0] = vertex_set[index_set(cell_index, 6)];
          _iso_coeff[n][n][n] = vertex_set[index_set(cell_index, 7)];

          //     vNN0------vNNN
          //    ./'|      ./'|
          // vN00------vN0N  |
          //   |   |     |   |    z
          //   | v0N0----|-v0NN   ^   y
          //   |./'      |./'     |./'
          // v000------v00N       +--> x

          // interpolate inner edge points linearly
          for(int i(1); i < degree_; ++i)
          {
            const DataType alpha = DataType(i) / DataType(degree_);

            // edges parallel to X axis
            // 00x edge
            _c1(_iso_coeff[0][0][i], alpha, _iso_coeff[0][0][0], _iso_coeff[0][0][n]);
            // 0Nx edge
            _c1(_iso_coeff[0][n][i], alpha, _iso_coeff[0][n][0], _iso_coeff[0][n][n]);
            // N0x edge
            _c1(_iso_coeff[n][0][i], alpha, _iso_coeff[n][0][0], _iso_coeff[n][0][n]);
            // NNx edge
            _c1(_iso_coeff[n][n][i], alpha, _iso_coeff[n][n][0], _iso_coeff[n][n][n]);

            // edges parallel to Y axis
            // 0y0 edge
            _c1(_iso_coeff[0][i][0], alpha, _iso_coeff[0][0][0], _iso_coeff[0][n][0]);
            // 0yN edge
            _c1(_iso_coeff[0][i][n], alpha, _iso_coeff[0][0][n], _iso_coeff[0][n][n]);
            // Ny0 edge
            _c1(_iso_coeff[n][i][0], alpha, _iso_coeff[n][0][0], _iso_coeff[n][n][0]);
            // NyN edge
            _c1(_iso_coeff[n][i][n], alpha, _iso_coeff[n][0][n], _iso_coeff[n][n][n]);

            // edges parallel to Z axis
            // z00 edge
            _c1(_iso_coeff[i][0][0], alpha, _iso_coeff[0][0][0], _iso_coeff[n][0][0]);
            // z0N edge
            _c1(_iso_coeff[i][0][n], alpha, _iso_coeff[0][0][n], _iso_coeff[n][0][n]);
            // zN0 edge
            _c1(_iso_coeff[i][n][0], alpha, _iso_coeff[0][n][0], _iso_coeff[n][n][0]);
            // zNN edge
            _c1(_iso_coeff[i][n][n], alpha, _iso_coeff[0][n][n], _iso_coeff[n][n][n]);
          }

          // iso-coeff array index offsets and increments
          const int eox[12] = {1, 1, 1, 1, 0, n, 0, n, 0, n, 0, n};
          const int eoy[12] = {0, n, 0, n, 1, 1, 1, 1, 0, 0, n, n};
          const int eoz[12] = {0, 0, n, n, 0, 0, n, n, 1, 1, 1, 1};
          const int eix[12] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
          const int eiy[12] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0};
          const int eiz[12] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1};

          // loop over all twelve edges of the element
          for(int e = 0; e < 12; ++e)
          {
            // get global edge index
            const Index edge_idx = edges_at_hexa(cell_index, e);

            // get chart for that edge
            const ChartType* chart = _charts_1.at(edge_idx);

            // no chart?
            if(chart == nullptr)
              continue;

            // initialize indices for this edge
            int z = eoz[e];
            int y = eoy[e];
            int x = eox[e];

            // loop over all edge points and project them
            for(int i(1); i < degree_; ++i, z += eiz[e], y += eiy[e], x += eix[e])
            {
              _iso_coeff[z][y][x] = chart->project(_iso_coeff[z][y][x]);
            }
          }

          // interpolate inner quad points bilinearly
          for(int i(1); i < degree_; ++i)
          {
            const DataType alpha_i = DataType(i) / DataType(degree_);

            for(int j(1); j < degree_; ++j)
            {
              const DataType alpha_j = DataType(j) / DataType(degree_);

              // quads parallel to XY plane
              _c2(_iso_coeff[0][i][j],
                alpha_i, _iso_coeff[0][0][j], _iso_coeff[0][n][j],
                alpha_j, _iso_coeff[0][i][0], _iso_coeff[0][i][n]);
              _c2(_iso_coeff[n][i][j],
                alpha_i, _iso_coeff[n][0][j], _iso_coeff[n][n][j],
                alpha_j, _iso_coeff[n][i][0], _iso_coeff[n][i][n]);

              // quads parallel to XZ plane
              _c2(_iso_coeff[i][0][j],
                alpha_i, _iso_coeff[0][0][j], _iso_coeff[n][0][j],
                alpha_j, _iso_coeff[i][0][0], _iso_coeff[i][0][n]);
              _c2(_iso_coeff[i][n][j],
                alpha_i, _iso_coeff[0][n][j], _iso_coeff[n][n][j],
                alpha_j, _iso_coeff[i][n][0], _iso_coeff[i][n][n]);

              // quads parallel to YZ plane
              _c2(_iso_coeff[i][j][0],
                alpha_i, _iso_coeff[0][j][0], _iso_coeff[n][j][0],
                alpha_j, _iso_coeff[i][0][0], _iso_coeff[i][n][0]);
              _c2(_iso_coeff[i][j][n],
                alpha_i, _iso_coeff[0][j][n], _iso_coeff[n][j][n],
                alpha_j, _iso_coeff[i][0][n], _iso_coeff[i][n][n]);
            }
          }

          // iso-coeff array index offsets and i/j increments
          const int qox[6] = {1, 1, 1, 1, 0, n};
          const int qoy[6] = {1, 1, 0, n, 1, 1};
          const int qoz[6] = {0, n, 1, 1, 1, 1};
          const int qix[6] = {0, 0, 0, 0, 0, 0};
          const int qjx[6] = {1, 1, 1, 1, 0, 0};
          const int qiy[6] = {1, 1, 0, 0, 0, 0};
          const int qjy[6] = {0, 0, 0, 0, 1, 1};
          const int qiz[6] = {0, 0, 1, 1, 1, 1};
          const int qjz[6] = {0, 0, 0, 0, 0, 0};

          // loop over all six quads of the element
          for(int q = 0; q < 6; ++q)
          {
            // get global quad index
            const Index quad_idx = quads_at_hexa(cell_index, q);

            // get chart for that quad
            const ChartType* chart = _charts_2.at(quad_idx);

            // no chart?
            if(chart == nullptr)
              continue;

            // initialize indices for this quad
            int z = qoz[q];
            int y = qoy[q];
            int x = qox[q];

            // loop over all inner quad points and project them
            for(int i(1); i < degree_; ++i, z += qiz[q], y += qiy[q], x += qix[q])
            {
              for(int j(1); j < degree_; ++j, z += qjz[q], y += qjy[q], x += qjx[q])
              {
                _iso_coeff[z][y][x] = chart->project(_iso_coeff[z][y][x]);
              }
            }
          }

          // interpolate inner hexa points linearly
          for(int i(1); i < degree_; ++i)
          {
            const DataType alpha_i = DataType(i) / DataType(degree_);
            for(int j(1); j < degree_; ++j)
            {
              const DataType alpha_j = DataType(j) / DataType(degree_);
              for(int k(1); k < degree_; ++k)
              {
                const DataType alpha_k = DataType(k) / DataType(degree_);
                _c3(_iso_coeff[i][j][k],
                  alpha_i, _iso_coeff[0][j][k], _iso_coeff[n][j][k],
                  alpha_j, _iso_coeff[i][0][k], _iso_coeff[i][n][k],
                  alpha_k, _iso_coeff[i][j][0], _iso_coeff[i][j][n]);
              }
            }
          }

          // get chart for this hexa
          const ChartType* chart = _charts_3.at(cell_index);
          if(chart != nullptr)
          {
            // loop over all inner hexa points and project them
            for(int i(1); i < degree_; ++i)
            {
              for(int j(1); j < degree_; ++j)
              {
                for(int k(1); k < degree_; ++k)
                {
                  _iso_coeff[i][j][k] = chart->project(_iso_coeff[i][j][k]);
                }
              }
            }
          }
        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         */
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          DataType vx[degree_+1], vy[degree_+1], vz[degree_+1];
          Intern::Basis<degree_>::val(vx, dom_point[0]);
          Intern::Basis<degree_>::val(vy, dom_point[1]);
          Intern::Basis<degree_>::val(vz, dom_point[2]);

          img_point.format();
          for(int i(0); i <= degree_; ++i)
          {
            for(int j(0); j <= degree_; ++j)
            {
              for(int k(0); k <= degree_; ++k)
              {
                img_point.axpy(vz[i]*vy[j]*vx[k], _iso_coeff[i][j][k]);
              }
            }
          }
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         */
        void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& dom_point) const
        {
          DataType vx[degree_+1], vy[degree_+1], vz[degree_+1];
          DataType dx[degree_+1], dy[degree_+1], dz[degree_+1];

          Intern::Basis<degree_>::val(vx, dom_point[0]);
          Intern::Basis<degree_>::val(vy, dom_point[1]);
          Intern::Basis<degree_>::val(vz, dom_point[2]);
          Intern::Basis<degree_>::d1(dx, dom_point[0]);
          Intern::Basis<degree_>::d1(dy, dom_point[1]);
          Intern::Basis<degree_>::d1(dz, dom_point[2]);

          jac_mat.format();
          for(int l(0); l < image_dim; ++l)
          {
            for(int i(0); i <= degree_; ++i)
            {
              for(int j(0); j <= degree_; ++j)
              {
                for(int k(0); k <= degree_; ++k)
                {
                  jac_mat(l,0) += vz[i] * vy[j] * dx[k] * _iso_coeff[i][j][k][l];
                  jac_mat(l,1) += vz[i] * dy[j] * vx[k] * _iso_coeff[i][j][k][l];
                  jac_mat(l,2) += dz[i] * vy[j] * vx[k] * _iso_coeff[i][j][k][l];
                }
              }
            }
          }
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         */
        void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& dom_point) const
        {
          DataType vx[degree_+1], vy[degree_+1], vz[degree_+1];
          DataType dx[degree_+1], dy[degree_+1], dz[degree_+1];
          DataType qx[degree_+1], qy[degree_+1], qz[degree_+1];

          Intern::Basis<degree_>::val(vx, dom_point[0]);
          Intern::Basis<degree_>::val(vy, dom_point[1]);
          Intern::Basis<degree_>::val(vz, dom_point[2]);
          Intern::Basis<degree_>::d1(dx, dom_point[0]);
          Intern::Basis<degree_>::d1(dy, dom_point[1]);
          Intern::Basis<degree_>::d1(dz, dom_point[2]);
          Intern::Basis<degree_>::d2(qx, dom_point[0]);
          Intern::Basis<degree_>::d2(qy, dom_point[1]);
          Intern::Basis<degree_>::d2(qz, dom_point[2]);

          hess_ten.format();
          for(int l(0); l < image_dim; ++l)
          {
            for(int i(0); i <= degree_; ++i)
            {
              for(int j(0); j <= degree_; ++j)
              {
                for(int k(0); k <= degree_; ++k)
                {
                  hess_ten(l,0,0) += vz[i] * vy[j] * qx[k] * _iso_coeff[i][j][k][l];
                  hess_ten(l,0,1) += vz[i] * dy[j] * dx[k] * _iso_coeff[i][j][k][l];
                  hess_ten(l,0,2) += dz[i] * vy[j] * dx[k] * _iso_coeff[i][j][k][l];
                  hess_ten(l,1,1) += vz[i] * qy[j] * vx[k] * _iso_coeff[i][j][k][l];
                  hess_ten(l,1,2) += dz[i] * dy[j] * vx[k] * _iso_coeff[i][j][k][l];
                  hess_ten(l,2,2) += qz[i] * vy[j] * vx[k] * _iso_coeff[i][j][k][l];
                }
              }
            }
            hess_ten(l,1,0) = hess_ten(l,0,1);
            hess_ten(l,2,0) = hess_ten(l,0,2);
            hess_ten(l,2,1) = hess_ten(l,1,2);
          }
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          XASSERTM(false, "volume computation not available for isoparametric trafo");
          return DataType(0);
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        DataType width_directed(const ImagePointType&) const
        {
          XASSERTM(false, "cell width computation not available for isoparametric trafo");
          return DataType(0);
        }
      }; // class Evaluator<Hypercube<3>,...>
    } // namespace Isoparam
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_ISOPARAM_EVALUATOR_HPP
