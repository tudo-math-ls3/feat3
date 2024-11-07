// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/analytic/static_wrapper.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/cgal.hpp>

namespace FEAT
{
  namespace Analytic
  {
    /**
     * \brief Analytic Distance namespace
     *
     * This namespace encapsulated commonly used functions,
     * that represent distance functions.
     */
    namespace Distance
    {
      /**
       * \brief Analytic distance function
       *
       * This class implements the AnalyticFunction interface representing the distance function
       * \f$ f(x) =\| x - x_0 \|_2 \f$
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \warning Because the function is differentiable everywhere except for x_0, Bad Things (TM) might happen if
       * someone wants to compute the gradient or hessian there, as the functions return 0.
       *
       * \tparam ImgPointType_
       * The type of \f$ x_0 \f$.
       *
       * \author Jordi Paul
       */
      template<int dim_, typename DataType_>
      class DistanceFunction :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// our origin
          const PointType _origin;

        public:
          /// Constructor
          explicit Evaluator(const DistanceFunction& function) :
            _origin(function._origin)
          {
          }

          ValueType value(const PointType& point) const
          {
            return (point - _origin).norm_euclid();
          }

          GradientType gradient(const PointType& point) const
          {
            DataType norm (this->value(point));

            if(norm < Math::eps<DataType>())
              return GradientType::null();

            return (DataType(1) / norm) * (point - _origin);
          }

          HessianType hessian(const PointType& point) const
          {
            DataType norm (this->value(point));
            if(norm < Math::eps<DataType>())
              return HessianType::null();

            HessianType hess (DataType(0));
            norm = DataType(1)/norm;
            DataType denom = Math::sqr(norm)*norm;

            for(int i(0); i < dim_; ++i)
            {
              hess[i][i] = norm;
              for(int j(0); j < dim_; ++j)
              {
                hess[i][j] -= ((point[i] - _origin[i]) * (point[j] - _origin[j]))*denom;
              }
            }
            return hess;
          }
        }; // class DistanceFunction::Evaluator<...>

      public:
        /// Point to calculate the distance from
        PointType _origin;

      public:
        /// default constructor
        DistanceFunction()
        {
          _origin.format();
        }

        /// Constructor
        explicit DistanceFunction(const PointType& origin) :
          _origin(origin)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType& origin)
        {
          _origin = origin;
        }
      }; // class DistanceFunction

      /**
       * \brief Scaled and displaced analytic distance function
       *
       * This class implements the AnalyticFunction interface representing the scaled and displaced distance function
       * \f$ f(x) = a + b \| x - x_0 \|_2 \f$
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \warning Because the function is differentiable everywhere except for x_0, Bad Things (TM) might happen if
       * someone wants to compute the gradient or hessian there, as the functions return 0.
       *
       * \tparam ImgPointType_
       * The type of \f$ x_0 \f$.
       *
       * \author Jordi Paul
       */
      template<int dim_, typename DataType_>
      class DistanceFunctionSD :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;

        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// our origin
          const PointType _origin;
          /// displacement and scaling factor
          const DataType _a, _b;

        public:
          /// Constructor
          explicit Evaluator(const DistanceFunctionSD& function) :
            _origin(function._origin),
            _a(function._a),
            _b(function._b)
          {
          }

          ValueType value(const PointType& point) const
          {
            return _a + _b*(point - _origin).norm_euclid();
          }

          GradientType gradient(const PointType& point) const
          {
            DataType norm ((point - _origin).norm_euclid());
            if(norm < Math::eps<DataType>())
              return GradientType::null();

            return (_b / norm) * (point - _origin);
          }

          HessianType hessian(const PointType& point) const
          {
            DataType norm ((point - _origin).norm_euclid());
            if(norm <= Math::eps<DataType>())
              return HessianType::null();

            HessianType hess(DataType(0));

            norm = DataType(1)/norm;
            DataType denom = _b*Math::sqr(norm)*norm;

            for(int i(0); i < dim_; ++i)
            {
              hess[i][i] = _b*norm;
              for(int j(0); j < dim_; ++j)
              {
                hess[i][j] -= ((point[i] - _origin[i]) * (point[j] - _origin[j]))*denom;
              }
            }

            return hess;
          }
        }; // class DistanceFunctionSD::Evaluator<...>

      private:
        /// The point to which the distance to is calculated
        PointType _origin;
        /// Displacement of the function
        DataType_ _a;
        /// Scaling factor
        DataType_ _b;

      public:
        /// Constructor
        explicit DistanceFunctionSD(const PointType& origin, const DataType_ a_, const DataType_ b_) :
          _origin(origin),
          _a(a_),
          _b(b_)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType& origin)
        {
          _origin = origin;
        }
      }; // class DistanceFunctionSD

      /**
       * \brief Scaled and displaced analytic distance from i-th coordinate axis/plane function
       *
       * This class implements the AnalyticFunction interface representing the scaled and displaced distance function
       * \f$ f(x) = b | (x - x_0)_i | \f$
       *
       * This class supports function values and gradients for all dimensions.
       *
       * \tparam component_
       * Index of the coordinate axis/plane.
       *
       * \tparam ImgPointType_
       * The type of \f$ x_0 \f$.
       *
       * \author Jordi Paul
       */
      template<int component_, int dim_, typename DataType_>
      class PlaneDistanceFunctionSD :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;

        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// our origin
          const PointType _origin;
          /// our scaling
          const DataType_ _b;

        public:
          /// Constructor
          explicit Evaluator(const PlaneDistanceFunctionSD& function) :
            _origin(function._origin),
            _b(function._b)
          {
          }

          ValueType value(const PointType& point) const
          {
            return _b * Math::abs(point[component_] - _origin[component_]);
          }

          GradientType gradient(const PointType& point) const
          {
            DataType norm = this->value(point);
            if(norm < Math::eps<DataType>())
              return GradientType::null();

            GradientType grad;
            grad.format();
            grad[component_] = _b * Math::signum(point[component_] - _origin[component_]);
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
          }
        }; // class PlaneDistanceFunctionSD::Evaluator<...>

      private:
        /// The point to which the distance to is calculated
        PointType _origin;
        /// Scaling factor
        DataType_ _b;

      public:
        /// Constructor
        explicit PlaneDistanceFunctionSD(const PointType& origin, const DataType_ b_) :
          _origin(origin),
          _b(b_)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType& origin)
        {
          _origin = origin;
        }
      }; // class PlaneDistanceFunctionSD

      /**
       * \brief Analytic inverse distance function
       *
       * This class implements the AnalyticFunction interface representing the inverse distance function
       * \f$ f(x) = \frac{1}{\| x - x_0 \|_2} \f$
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \warning Because the function is differentiable everywhere except for x_0, Bad Things (TM) might happen if
       * someone wants to compute the gradient or hessian there, as the functions return 0.
       *
       * \tparam ImgPointType_
       * The type of \f$ x_0 \f$.
       *
       * \author Peter Zajac
       */
      template<int dim_, typename DataType_>
      class InverseDistanceFunction :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// our origin
          const PointType _origin;

        public:
          /// Constructor
          explicit Evaluator(const InverseDistanceFunction& function) :
            _origin(function._origin)
          {
          }

          ValueType value(const PointType& point) const
          {
            return DataType(1) / (point - _origin).norm_euclid();
          }

          GradientType gradient(const PointType& point) const
          {
            const DataType norm_sqr = (point - _origin).norm_euclid_sqr();

            if(norm_sqr < Math::eps<DataType>())
              return GradientType::null();

            return (DataType(1) / Math::pow(norm_sqr, DataType(1.5))) * (_origin - point);
          }

          HessianType hessian(const PointType& point) const
          {
            const DataType norm_sqr = (point - _origin).norm_euclid_sqr();
            if(norm_sqr < Math::eps<DataType>())
              return HessianType::null();

            const DataType denom = DataType(1) / Math::pow(norm_sqr, DataType(2.5));
            HessianType hess (DataType(0));
            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                if(i == j)
                {
                  for(int k(0); k < dim_; ++k)
                    hess[i][j] += DataType(k == i ? 2 : -1) * Math::sqr(point[k] - _origin[k]) * denom;
                }
                else
                {
                  hess[i][j] = DataType(3) * (point[i] - _origin[i]) * (point[j] - _origin[j]) * denom;
                }
              }
            }
            return hess;
          }
        }; // class InverseDistanceFunction::Evaluator<...>

      public:
        /// Point to calculate the distance from
        PointType _origin;

      public:
        /// default constructor
        InverseDistanceFunction()
        {
          _origin.format();
        }

        /// Constructor
        explicit InverseDistanceFunction(const PointType& origin) :
          _origin(origin)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType& origin)
        {
          _origin = origin;
        }
      }; // class InverseDistanceFunction

#ifdef FEAT_HAVE_CGAL
      /**
       * \brief CGAL based analytic signed distance function
       *
       * This function can be used as an interface to an anayltic function of a cgal 3D mesh distance evaluation,
       * for example to evaluate an interpolation.
       * Thereby the definition is as follows:
       * dist(point) >= 0 if point is inside or on the boundary of the mesh.
       * dist(point) < 0 if point is outside of the geometry.
       *
       * \tparam DT_ The underlying datatype of the evaluation.
       * \note The actual cgal call always works with double precision.
       *
       * \author Maximilian Esser
       */
      template<typename DT_>
      class CGALDistanceFunction : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Scalar ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 3;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = false;
        /// We can't compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

      private:
        /// Wrapper to CGAL library
        Geometry::CGALWrapper<DT_> _cgal_wrapper;
      public:

        /**
         * \brief Constructor
         *
         * \param[in] filestream_
         *            A stringstream containing the file data.
         *
         * \param[in] filemode_
         *            The file mode of the data.
         */
        explicit CGALDistanceFunction(std::istream& filestream_, Geometry::CGALFileMode filemode_) :
          _cgal_wrapper(filestream_, filemode_)
        {
        }

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// cgal datatype
          typedef typename CGALDistanceFunction::DataType IDT_;
          const Geometry::CGALWrapper<IDT_>* _cgal_wrapper;

        public:
          explicit Evaluator(const CGALDistanceFunction& function) :
            _cgal_wrapper(&function._cgal_wrapper)
          {
          }

          ValueType value(const PointType& point)
          {
            return ((_cgal_wrapper->point_inside(IDT_(point[0]), IDT_(point[1]), IDT_(point[2])) ? ValueType(1) : ValueType(-1)) *
                    ValueType(Math::sqrt(_cgal_wrapper->squared_distance(IDT_(point[0]), IDT_(point[1]), IDT_(point[2])))));
          }
        }; // class CGALDistFunc::Evaluator<...>
      }; // class CGALDistFunc<...>
#endif // FEAT_HAVE_CGAL
    } // namespace Distance
  } // namespace Analytic
} // namespace FEAT
