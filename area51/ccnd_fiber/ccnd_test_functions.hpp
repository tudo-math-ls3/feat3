#pragma once
#ifndef AREA51_CCND_TEST_FUNCTION
#define AREA51_CCND_TEST_FUNCTION 1

#include <area51/ccnd_fiber/ccnd_fiber_common.hpp>

namespace CCND_FIBER
{
  using namespace FEAT;

  //template for pow operation... inlined...
  template<uint exp_>
  inline DataType pow(const DataType base)
  {
    DataType ret{DataType(1)};
    for(uint i(0u); i < exp_; ++i)
    {
      ret *= base;
    }
    return ret;
  }

  template<>
  inline DataType pow<1>(DataType base)
  {
    return base;
  }

  template<>
  inline DataType pow<2>(DataType base)
  {
    return base * base;
  }

  template<>
  inline DataType pow<3>(DataType base)
  {
    return pow<2>(base) * base;
  }

  template<>
  inline DataType pow<4>(DataType base)
  {
    return pow<2>(base) * pow<2>(base);
  }

  template<>
  inline DataType pow<5>(DataType base)
  {
    return pow<3>(base) * pow<2>(base);
  }

  template<>
  inline DataType pow<6>(DataType base)
  {
    return pow<3>(base) * pow<3>(base);
  }

  template<>
  inline DataType pow<7>(DataType base)
  {
    return pow<3>(base) * pow<4>(base);
  }

  template<>
  inline DataType pow<8>(DataType base)
  {
    return pow<4>(base) * pow<4>(base);
  }

  template<int dim_>
  class SimpleTensor;

  template<>
  class SimpleTensor<2>
  {
  public:
    class Velo;
    class Pres;
    class RHS;
    class Robin;
    class Orient;
    class Fourth_Moment;

    DataType _vmax;
    DataType _mu;
    DataType _rho;
    DataType _N_s;
    DataType _N_p;
    DataType _period;
    DataType _epsilon;
    DataType _alpha;
    bool _konv;

    explicit SimpleTensor(DataType vmax, DataType mu, DataType rho, DataType N_s, DataType N_p, DataType epsilon = DataType(1), DataType period = DataType(2) * DataType(Math::pi<DataType>()), DataType alpha = DataType(1.), bool konv = false) :
    _vmax(vmax),
    _mu(mu),
    _rho(rho),
    _N_s(N_s),
    _N_p(N_p),
    _period(period),
    _epsilon(epsilon),
    _alpha(alpha),
    _konv(konv)
    {}

    class Velo :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;
      static constexpr bool can_grad = true;

    protected:
      const DataType _vmax;
      const DataType _period;
      const DataType _alpha;

    private:
      DataType _t;


    public:
      explicit Velo(const SimpleTensor& func) :
      _vmax(func._vmax),
      _period(func._period),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;

        const DataType _vmax, _period, _alpha, _t;

      public:
        explicit Evaluator(const Velo& function) :
        _vmax(function._vmax),
        _period(function._period),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType exp = Math::exp(_alpha * _t);

          ValueType val;
          val[0] = exp * _vmax * Math::cos(x) * Math::sin(y);
          val[1] = DataType(-1) * exp * _vmax * Math::sin(x) * Math::cos(y);
          return val;
        }

        GradientType gradient(const PointType& point)
        {
          //with gradient we mean the Jacobian Matrix here
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType exp = Math::exp(_alpha * _t);

          GradientType grad;
          grad[0][0] = - exp * _vmax * _period * Math::sin(x) * Math::sin(y);
          grad[0][1] = exp * _vmax * _period * Math::cos(x) * Math::cos(y);
          grad[1][0] = - exp * _vmax * _period * Math::cos(x) * Math::cos(y);
          grad[1][1] = exp * _vmax * _period * Math::sin(x) * Math::sin(y);

          return grad;
        }

      };
    }; //class Velo


    class Pres :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Scalar ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _period;
      DataType _alpha;

    private:
      DataType _t;

    public:
      explicit Pres(const SimpleTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _period(func._period),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _period, _alpha, _t;

      public:
        explicit Evaluator(const Pres& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _period(function._period),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType exp = Math::exp(_alpha * _t);

          ValueType val;
          val = -DataType(2) * exp * _mu * _period * _vmax * Math::sin(x) * Math::sin(y);
          return val;
        }

      };
    }; // class Pres

    /**
     * \brief The assoicated right side for an orientation matrix
     */
    class RHS :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _period;
      DataType _alpha;
      bool _konv;

    private:
      DataType _t;

    public:
      explicit RHS(const SimpleTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _period(func._period),
      _alpha(func._alpha),
      _konv(func._konv),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho,_period, _alpha, _t;
        const bool _konv;

      public:
        explicit Evaluator(const RHS& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _period(function._period),
        _alpha(function._alpha),
        _t(function._t),
        _konv(function._konv)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType exp = Math::exp(_alpha * _t);

          ValueType val;
          {
            const DataType tau_grad = - DataType(2) * exp * _mu * _period * _period * _vmax * Math::cos(x) * Math::sin(y);
            const DataType convection = _konv ? (- DataType(0.5) * exp * exp * _period * _rho * _vmax * _vmax * Math::sin(DataType(2) * x)) : DataType(0);
            const DataType p_grad = - DataType(2) * exp * _mu * _period * _period * _vmax * Math::cos(x) * Math::sin(y);
            const DataType react = _alpha * exp * _vmax * Math::cos(x) * Math::sin(y);

            val[0] = (react + convection + p_grad - tau_grad)/_rho;
          }

          {
            const DataType tau_grad = DataType(2) * exp * _mu * _period * _period * _vmax * Math::cos(y) * Math::sin(x);
            const DataType convection = _konv ? (- DataType(0.5) * exp * exp * _period * _rho * _vmax * _vmax * Math::sin(DataType(2) * y)) : DataType(0);
            const DataType p_grad = - DataType(2) * exp * _mu * _period * _period * _vmax * Math::cos(y) * Math::sin(x);
            const DataType react = DataType(-1) * _alpha * exp * _vmax * Math::sin(x) * Math::cos(y);

            val[1] = (react + convection + p_grad - tau_grad)/_rho;
          }

          return val;
        }
      };
    }; //class RHS

    class Robin :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    private:
      DataType _t;

    public:
      explicit Robin(const SimpleTensor&) :
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

      public:
        explicit Evaluator(const Robin&)
        {
        }

        ValueType value(const PointType&)
        {
          ValueType val;
          val[0] = val[1] = DataType(0);
          return val;
        }
      };
    }; //class Robin

    class Orient :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon;
    public:
      explicit Orient(const SimpleTensor& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      void set_time(const DataType)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon;
      public:
        explicit Evaluator(const Orient& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType&)
        {
          ValueType val;
          val[0] = DataType(0);
          val[2] = DataType(0);
          //           val[2] = val[1];
          val[1] = DataType(0);

          return val;
        }
      };
    }; // class Orient

    class Fourth_Moment :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<5> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon;
    public:
      explicit Fourth_Moment(const SimpleTensor& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      void set_time(const DataType)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon;
      public:
        explicit Evaluator(const Fourth_Moment& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType&)
        {
          ValueType val;
          val[0] = DataType(0);
          val[1] = DataType(0);
          val[2] = DataType(0);
          val[3] = DataType(0);
          val[4] = DataType(0);

          return val;
        }
      };
    }; // class Fourth_Moment

  }; //SimpleTensor


  template<int dim_>
  class FullTensor;


   /////////////////////////////////////////////////
  // Bubble function constructed to meet zero Neumann boundary on right side for space dependent 4th order orientation, Christophs pyhsical formulation...
  ////////////////////////////////////////////////////
  /**
   * \brief Divergence-free 2-d bubble function class
   *
   * This class initialises all necessary parameters and provides an interface to access generalized test functions.
   */
  template<>
  class FullTensor<2>
  {
  public:
    class Velo;
    class Pres;
    class RHS;
    class RHS_noreg;
    class Robin;
    class Orient;
    class Fourth_Moment;

    DataType _vmax;
    DataType _mu;
    DataType _rho;
    DataType _N_s;
    DataType _N_p;
    DataType _period;
    DataType _epsilon;
    DataType _alpha;
    bool _konv;

    explicit FullTensor(DataType vmax, DataType mu, DataType rho, DataType N_s, DataType N_p, DataType epsilon = DataType(1), DataType period = DataType(2) * DataType(Math::pi<DataType>()), DataType alpha = DataType(0.), bool konv = false) :
      _vmax(vmax),
      _mu(mu),
      _rho(rho),
      _N_s(N_s),
      _N_p(N_p),
      _period(period),
      _epsilon(epsilon),
      _alpha(alpha),
      _konv(konv)
      {}

    class Velo :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;
      static constexpr bool can_grad = true;

    protected:
      const DataType _vmax;
      const DataType _period;


    public:
      explicit Velo(const FullTensor& func) :
      _vmax(func._vmax),
      _period(func._period)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;

        const DataType _vmax, _period;

      public:
        explicit Evaluator(const Velo& function) :
        _vmax(function._vmax),
        _period(function._period)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];

          ValueType val;
          val[0] = _vmax * Math::cos(x) * Math::sin(y);
          val[1] = DataType(-1) * _vmax * Math::sin(x) * Math::cos(y);
          return val;
        }

        GradientType gradient(const PointType& point)
        {
          //with gradient we mean the Jacobian Matrix here
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];

          GradientType grad;
          grad[0][0] = - _vmax * _period * Math::sin(x) * Math::sin(y);
          grad[0][1] = _vmax * _period * Math::cos(x) * Math::cos(y);
          grad[1][0] = -_vmax * _period * Math::cos(x) * Math::cos(y);
          grad[1][1] = _vmax * _period * Math::sin(x) * Math::sin(y);

          return grad;
        }

      };
    }; //class Velo

    class Pres :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Scalar ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;

    public:
      explicit Pres(const FullTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_p, _period, _epsilon;

      public:
        explicit Evaluator(const Pres& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];

          ValueType val;
          val = -DataType(2) * _mu * _period * _vmax * Math::sin(x) * Math::sin(y);
          return val;
        }

      };
    }; // class Pres

    /**
     * \brief The assoicated right side for an orientation matrix
     */
    class RHS :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;
      DataType _alpha;
      bool _konv;

    public:
      explicit RHS(const FullTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _konv(func._konv)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _epsilon, _alpha;
        const bool _konv;

      public:
        explicit Evaluator(const RHS& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _konv(function._konv)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];

          ValueType val;
          {
            //first val component
            const DataType taugrad = - DataType(2) * _mu *(_N_p *((_period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(2) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) *
                                      pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                      pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                      pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * _vmax * Math::cos(x) * Math::sin(y) *(_epsilon + pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /
                                      (DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax * Math::cos(x) *
                                      Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * Math::sin(x) * Math::sin(y) *(_epsilon + pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(4) * _period *
                                      pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                      pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) +
                                      DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y))) *(DataType(4) * _period *
                                      pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                      pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) +
                                      DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) + pow<2>(_period) * _vmax * Math::cos(x) * Math::sin(y) -(DataType(4) * _N_s * pow<2>(_period) * pow<3>(_vmax) * Math::cos(x) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)))
                                      /(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) +(DataType(2) * _N_s * pow<2>(_period) * _vmax * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::sin(y)))) /(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) -(DataType(2) * _N_s * _period * _vmax * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * Math::sin(x) - DataType(2) * _period * pow<2>(_vmax) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)))) - DataType(2) * _N_p * _mu *((pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) +
                                      pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                      pow<3>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(DataType(2) *
                                      pow<2>(_period) * pow<5>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) *(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) +
                                      pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * pow<5>(_vmax) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<2>(Math::sin(y)) *
                                      (DataType(4) * _period * pow<4>(_vmax) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                      pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) *
                                      pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * pow<5>(_vmax) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(4) * _period * pow<4>(_vmax) *
                                      pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) +
                                      DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                      pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))));

            const DataType pgrad = DataType(-2) * _mu * _period * _period * _vmax * Math::cos(x) * Math::sin(y);

            const DataType konvektion = _konv ? (DataType(-0.5) * _period * _rho * _vmax * _vmax * Math::sin(DataType(2) * x)) : DataType(0);

            const DataType react = _alpha * _vmax * Math::cos(x) * Math::sin(y);

            val[0] = (react + konvektion + pgrad - taugrad)/_rho;
          }

          {
            //second component
            const DataType taugrad = DataType(2) * _mu *(_N_p *((_period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(2) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) - DataType(2) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                        pow<2>(Math::sin(x)) *(Math::sin(y)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                        pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(y)) * pow<5>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) *
                                        pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * _vmax * Math::cos(y) * Math::sin(x) *(_epsilon + pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(8) *
                                        _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax * Math::cos(y) *
                                        Math::sin(x) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                        pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * Math::sin(x) * Math::sin(y) *(_epsilon + pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(4) * _period * pow<4>(_vmax) *
                                        pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) +
                                        DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                        pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x))) *(DataType(4) * _period * pow<4>(_vmax) *
                                        pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) +
                                        DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                        pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) + pow<2>(_period) * _vmax * Math::cos(y) * Math::sin(x) -(DataType(4) * _N_s * pow<2>(_period) * pow<3>(_vmax) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y))) /
                                        (DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) +(DataType(2) * _N_s * pow<2>(_period) * _vmax * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))))
                                        /(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) -(DataType(2) * _N_s * _period * _vmax * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(Math::cos(y)) *
                                        pow<2>(Math::sin(x))) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * Math::cos(y) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                        pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)))) + DataType(2) * _N_p * _mu *((pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(y)) * pow<5>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) +
                                        pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                        pow<2>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) *
                                        pow<5>(_vmax) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                        pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) *
                                        pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * pow<5>(_vmax) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(4) *
                                        _period * pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                        pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) +
                                        DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * pow<5>(_vmax) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) *
                                        pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                        DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) *
                                        pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))));

            const DataType pgrad = DataType(-2) * _mu * _period * _period * _vmax * Math::cos(y) * Math::sin(x);

            const DataType konvektion = _konv ? (DataType(-0.5)* _period * _rho * _vmax * _vmax * Math::sin(DataType(2)*y)) : DataType(0);

            const DataType react = -_alpha * _vmax * Math::sin(x) * Math::cos(y);

            val[1] = (react + konvektion + pgrad - taugrad)/_rho;
          }

          return val;
        }
      };
    }; // class RHS

    class RHS_noreg :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _alpha;
      bool _konv;

    public:
      explicit RHS_noreg(const FullTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _alpha(func._alpha),
      _konv(func._konv)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _alpha;
        const bool _konv;

      public:
        explicit Evaluator(const RHS_noreg& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _alpha(function._alpha),
        _konv(function._konv)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];

          ValueType val;
          {
            //first val component
            const DataType taugrad = DataType(-2) * _mu *(_N_p *((pow<2>(_period) * pow<5>(_vmax) * pow<5>(Math::cos(x)) * pow<5>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                        pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) *
                                        pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * pow<5>(_vmax) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y)) *(DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) -
                                        DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) *
                                        pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))
                                        -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                        pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<3>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) *
                                        pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * pow<5>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(4) * _period *
                                        pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                        pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) *
                                        pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) + pow<2>(_period) * _vmax * Math::cos(x) * Math::sin(y) +(DataType(2) * _N_s * pow<2>(_period) * pow<3>(_vmax) * pow<3>(Math::cos(x)) * pow<3>(Math::sin(y))) /(pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                        pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) -(DataType(4) * _N_s * pow<2>(_period) * pow<3>(_vmax) * Math::cos(x) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) *
                                        pow<2>(Math::sin(x))) -(DataType(2) * _N_s * _period * pow<3>(_vmax) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * Math::sin(x) - DataType(2) * _period * pow<2>(_vmax) * Math::cos(x) *
                                        Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)))) - DataType(2) * _N_p * _mu *((pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y)))
                                        /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) *
                                        pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                        pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) *(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) *
                                        pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<3>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) *
                                        pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * pow<5>(_vmax) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(4) *
                                        _period * pow<4>(_vmax) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                        pow<3>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                        pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * pow<5>(_vmax) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(4) * _period * pow<4>(_vmax) * pow<4>(Math::cos(x)) * Math::cos(y) *
                                        pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) + DataType(4) * _period *
                                        pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                        pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))));

            const DataType pgrad = DataType(-2) * _mu * _period * _period * _vmax * Math::cos(x) * Math::sin(y);

            const DataType konvektion = _konv ? (DataType(-0.5) * _period * _rho * _vmax * _vmax * Math::sin(DataType(2) * x)) : DataType(0);

            const DataType react = _alpha * _vmax * Math::cos(x) * Math::sin(y);

            val[0] = (react + konvektion + pgrad - taugrad)/_rho;
          }

          {
            //second component
            const DataType taugrad = DataType(2) * _mu *(_N_p *((pow<2>(_period) * pow<5>(_vmax) * pow<5>(Math::cos(y)) * pow<5>(Math::sin(x))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                          pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(y)) * pow<5>(Math::sin(x)) * pow<2>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                          pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * pow<5>(_vmax) * pow<4>(Math::cos(y)) * pow<5>(Math::sin(x)) * Math::sin(y) *(DataType(4) * _period * pow<4>(_vmax) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) *
                                          _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                          pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -
                                          (DataType(3) * pow<2>(_period) * pow<5>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                          pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) *
                                          pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * pow<5>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(4) * _period * pow<4>(_vmax) *
                                          pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) +
                                          DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                          pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) + pow<2>(_period) * _vmax * Math::cos(y) * Math::sin(x) +(DataType(2) * _N_s * pow<2>(_period) * pow<3>(_vmax) * pow<3>(Math::cos(y)) * pow<3>(Math::sin(x))) /(pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) +
                                          pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) -(DataType(4) * _N_s * pow<2>(_period) * pow<3>(_vmax) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y))) /(pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)))
                                          -(DataType(2) * _N_s * _period * pow<3>(_vmax) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * Math::cos(y) *
                                          pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)))) + DataType(2) * _N_p * _mu *((pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(y)) * pow<5>(Math::sin(x)) * pow<2>(Math::sin(y))) /
                                          (pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<2>(Math::cos(x)) *
                                          pow<3>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))
                                          +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                          pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y))) /(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) *
                                          pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * pow<5>(_vmax) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(4) * _period *
                                          pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                          pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) *
                                          pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * pow<5>(_vmax) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) *
                                          pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) *
                                          pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                          pow<2>(Math::sin(y))));

            const DataType pgrad = DataType(-2) * _mu * _period * _period * _vmax * Math::cos(y) * Math::sin(x);

            const DataType konvektion = _konv ? (DataType(-0.5)* _period * _rho * _vmax * _vmax * Math::sin(DataType(2)*y)) : DataType(0);

            const DataType react = -_alpha * _vmax * Math::sin(x) * Math::cos(y);

            val[1] = (react + konvektion + pgrad - taugrad)/_rho;
          }

          return val;
        }
      };
    }; // class RHS_noreg

    class Robin :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;

    public:
      explicit Robin(const FullTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _epsilon;

      public:
        explicit Evaluator(const Robin& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          ValueType val;
          val[0] = -(DataType(2) * _mu * _period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(4) * _N_p * pow<2>(_epsilon) + DataType(16) * _N_s * pow<2>(_epsilon) + _N_p * pow<6>(_vmax) * pow<6>(Math::cos(x)) * pow<6>(Math::sin(y)) + DataType(2) * _N_s * pow<6>(_vmax) * pow<6>(Math::cos(x)) *
                    pow<6>(Math::sin(y)) + DataType(2) * _N_p * _epsilon * pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + DataType(2) * _N_p * _epsilon * pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + DataType(16) * _N_s * _epsilon * pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) +
                    DataType(2) * _N_s * _epsilon * pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y))) + DataType(2) * _mu * _period * _vmax * pow<4>(Math::cos(y)) * pow<5>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _N_s * _epsilon * pow<4>(_vmax) - _N_p * pow<6>(_vmax) *
                    pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + DataType(2) * _N_s * pow<6>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) + DataType(2) * _mu * _period * _vmax * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _N_p * _epsilon * pow<2>(_vmax) +
                    DataType(4) * _N_s * pow<6>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) - DataType(2) * _N_p * _epsilon * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + DataType(4) * _N_s * _epsilon * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /((DataType(2) * _epsilon +
                    pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) *(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                    pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))));
          val[1] = (DataType(2) * _N_p * _mu * _period * pow<5>(_vmax) * Math::cos(x) * Math::cos(y) *(pow<2>(Math::cos(x)) - DataType(1)) *(pow<2>(Math::cos(y)) - DataType(1)) *(pow<2>(Math::cos(x)) - pow<2>(Math::cos(y)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) +
                    DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) - DataType(4) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<4>(Math::cos(y)) - DataType(4) * pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(4) * pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)));
          return val;
        }
      };
    }; // class Robin

    class Orient :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon;
    public:
      explicit Orient(const FullTensor& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon;
      public:
        explicit Evaluator(const Orient& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType& point)
        {
           const DataType x = _period * point[0];
           const DataType y = _period * point[1];
           const DataType temp_norm = DataType(2) * _epsilon + pow<2>(_vmax) * (pow<2>(Math::cos(x) * Math::sin(y)) + pow<2>(Math::cos(y) * Math::sin(x)));

          ValueType val;
          val[0] = (_vmax * Math::cos(x) * Math::sin(y) * _vmax * Math::cos(x) * Math::sin(y) + _epsilon) / temp_norm;
          val[2] = _vmax * Math::cos(x) * Math::sin(y) * DataType(-1) * _vmax * Math::sin(x) * Math::cos(y) / temp_norm;
//           val[2] = val[1];
          val[1] = (_vmax * Math::sin(x) * Math::cos(y) * _vmax * Math::sin(x) * Math::cos(y) + _epsilon) / temp_norm;

          return val;
        }
      };
    }; // class Orient

    class Fourth_Moment :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<5> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon;
    public:
      explicit Fourth_Moment(const FullTensor& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon;
      public:
        explicit Evaluator(const Fourth_Moment& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType temp_norm = DataType(8) * _epsilon + pow<4>(_vmax) * (pow<4>(Math::cos(x) * Math::sin(y)) + pow<4>(Math::cos(y) * Math::sin(x))
                                     + DataType(2)*(pow<2>(Math::cos(x) * Math::sin(y) * Math::sin(x) * Math::cos(y))));


          //we want to build the tensor A_(klmn} = v_k * v_l * v_m * v_n + \delta_{klmn}*eps
          //we represent this by a vector of the form:
          // A_{klmn} = val(d^3*k + d^2*l + d*m + n)
          //we have 2 different combs for all same, which are in val[0] and val[15]
          // 2 different combs for one different
          ValueType val;
//           val[0] = (pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + DataType(3) * _epsilon) / temp_norm;
//           val[1] = (pow<4>(_vmax) * pow<3>(Math::cos(x) * Math::sin(y)) * pow<1>(-Math::sin(x) * Math::cos(y))) / temp_norm;
//           val[2] = val[1];
//           val[3] = (pow<4>(_vmax) * pow<2>(Math::cos(x) * Math::sin(y)) * pow<2>(-Math::sin(x) * Math::cos(y)) + _epsilon) / temp_norm;
//           val[4] = val[1];
//           val[5] = val[3];
//           val[6] = val[3];
//           val[7] = (pow<4>(_vmax) * pow<1>(Math::cos(x) * Math::sin(y)) * pow<3>(-Math::sin(x) * Math::cos(y))) / temp_norm;
//           val[8] = val[1];
//           val[9] = val[3];
//           val[10] = val[3];
//           val[11] = val[7];
//           val[12] = val[3];
//           val[13] = val[7];
//           val[14] = val[7];
//           val[15] = (pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(3) * _epsilon) / temp_norm;

          //for this we use the definition of the fourth_order_contraction operator, see tiny_algebra.hpp
          val[0] = (pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + DataType(3) * _epsilon) / temp_norm;
          val[1] = (pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(3) * _epsilon) / temp_norm;
          val[2] = (pow<4>(_vmax) * pow<3>(Math::cos(x) * Math::sin(y)) * pow<1>(-Math::sin(x) * Math::cos(y))) / temp_norm;
          val[3] = (pow<4>(_vmax) * pow<1>(Math::cos(x) * Math::sin(y)) * pow<3>(-Math::sin(x) * Math::cos(y))) / temp_norm;
          val[4] = (pow<4>(_vmax) * pow<2>(Math::cos(x) * Math::sin(y)) * pow<2>(-Math::sin(x) * Math::cos(y)) + _epsilon) / temp_norm;

          return val;
        }
      };
    }; // class Fourth_Moment

  };//class FullTensor 2D spezialization

    /////////////////////////////////////////////////
  // A non-trivial divergence free velocity field and all associated functions with pyhsical, regularized tensors
  ////////////////////////////////////////////////////
  /**
   * \brief Divergence-free 3-D bubble function class
   *
   * This class initialises all necessary parameters and provides an interface to access generalized test functions.
   */
  template<>
  class FullTensor<3>
  {
  public:
    class Velo;
    class Pres;
    class RHS;
    class Robin;
    class Orient;
    class Fourth_Moment;

    DataType _vmax;
    DataType _mu;
    DataType _rho;
    DataType _N_s;
    DataType _N_p;
    DataType _period;
    DataType _epsilon;
    DataType _alpha;
    bool _konv;

    explicit FullTensor(DataType vmax, DataType mu, DataType rho, DataType N_s, DataType N_p, DataType epsilon = DataType(1), DataType period = DataType(2) * DataType(Math::pi<DataType>()), DataType alpha = DataType(0.), bool konv = false) :
      _vmax(vmax),
      _mu(mu),
      _rho(rho),
      _N_s(N_s),
      _N_p(N_p),
      _period(period),
      _epsilon(epsilon),
      _alpha(alpha),
      _konv(konv)
      {}

    class Velo :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;
      static constexpr bool can_grad = true;

    protected:
      const DataType _vmax;
      const DataType _period;


    public:
      explicit Velo(const FullTensor& func) :
      _vmax(func._vmax),
      _period(func._period)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;

        const DataType _vmax, _period;

      public:
        explicit Evaluator(const Velo& function) :
        _vmax(function._vmax),
        _period(function._period)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];

          ValueType val;
          val[0] = _vmax * Math::sin(x) * Math::sin(y) * z;
          val[1] = DataType(-1) * _vmax * Math::cos(x) * Math::cos(y)*z;
          val[2] = DataType(-1) * _vmax * _period * Math::cos(x) * Math::sin(y) * z * z;
          return val;
        }

        GradientType gradient(const PointType& point)
        {
          //with gradient we mean the Jacobian Matrix here
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];

          GradientType grad;
          grad[0][0] =  _vmax * _period * Math::cos(x) * Math::sin(y) * z;
          grad[0][1] = _vmax * _period * Math::sin(x) * Math::cos(y) * z;
          grad[0][2] = _vmax * Math::sin(x) * Math::sin(y);
          grad[1][0] = _vmax * _period * Math::sin(x) * Math::cos(y) * z;
          grad[1][1] = _vmax * _period * Math::cos(x) * Math::sin(y) * z;
          grad[1][2] = - _vmax * Math::cos(x) * Math::cos(y);
          grad[2][0] = _vmax * _period * _period * Math::sin(x) * Math::sin(y) * z * z;
          grad[2][1] = - _vmax * _period * _period * Math::cos(x) * Math::cos(y) * z * z;
          grad[2][2] = DataType(-2) *  _vmax * _period * Math::cos(x) * Math::sin(y) * z;

          return grad;
        }

      };
    }; //class Velo

    class Pres :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Scalar ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _period;

    public:
      explicit Pres(const FullTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _period(func._period)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _period;

      public:
        explicit Evaluator(const Pres& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _period(function._period)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];

          ValueType val;
          val = DataType(-2) * _mu * _period * _vmax * z * Math::cos(x) * Math::sin(y);
          return val;
        }

      };
    }; // class Pres

    /**
     * \brief The assoicated right side for an orientation matrix of typedef
     * A = (lambda1 0; 0 lambda2)
     */
    /**
     * \brief The associated pressure function to SteadyBubble2d
     */
    class RHS :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;
      DataType _alpha;
      bool _konv;

    public:
      explicit RHS(const FullTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _konv(func._konv)
      {
      }

      void print_variables() const
      {
        std::cout << "RHS print:\n";
        std::cout << "Vmax : " << _vmax << "\n";
        std::cout << "Mu : " << _mu << "\n";
        std::cout << "Rho : " << _rho << "\n";
        std::cout << "Ns : " << _N_s << "\n";
        std::cout << "Np : " << _N_p << "\n";
        std::cout << "period : " << _period << "\n";
        std::cout << "epsilon : " << _epsilon << "\n";
        std::cout << "Alpha : " << _alpha << "\n";
        std::cout << "Konv : " << _konv << "\n";
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _epsilon, _alpha;
        const bool _konv;

      public:
        explicit Evaluator(const RHS& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _konv(function._konv)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];

          //yes this is highly optimisable... if you want, have fun
          ValueType val;
          //RHS for first component
          {
            DataType grtau_d = - DataType(2) * _mu * pow<2>(_period) * _vmax * z * Math::sin(x) *(Math::sin(y));
            DataType grtau_np_first_term = DataType(2) * _N_p * _mu *((DataType(2) * pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax *
                                                pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) *
                                                z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(16) * pow<4>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(12) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *
                                                ((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                                Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(12) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<4>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(12) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) *
                                                _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *
                                                pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) *
                                                pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) *
                                                pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_second_term = - DataType(2) * _N_p * _mu *((DataType(2) * _period * _vmax * z * Math::cos(y) * Math::sin(x) *(DataType(2) * _period * pow<4>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) - DataType(2) * _period * pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *
                                                (_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)) *
                                                Math::sin(x) *(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(2) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) *
                                                pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(8) * pow<4>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) *
                                                _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) *
                                                Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * _vmax * z * Math::cos(y) * Math::sin(x) *(_epsilon +
                                                pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax *
                                                pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *
                                                (Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *
                                                pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<5>(_vmax) * pow<5>(z) *
                                                pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<5>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) +
                                                DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<3>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) *
                                                pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                (Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *
                                                ((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) *
                                                pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_third_term = - DataType(2) * _N_p * _mu *((_period * _vmax * z * Math::cos(x) * Math::sin(y) *(DataType(2) * _period * pow<4>(_vmax) * pow<4>(z) * Math::cos(x) *
                                                pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<4>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) *
                                                pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y)) *
                                                ((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(z) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y)))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * _vmax * z *
                                                Math::cos(x) * Math::sin(y) *(DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) - DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)))) /(DataType(15) *
                                                _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) *
                                                pow<3>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                                (DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(y)) * pow<5>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) *
                                                _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(6) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                                DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(8) * pow<2>(_period) *
                                                pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +
                                                (pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *
                                                (Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *
                                                (_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) *
                                                Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * Math::cos(x) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) *
                                                _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *
                                                (Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                                (DataType(2) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(z) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) -
                                                DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(_vmax) * pow<5>(z) * Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) /
                                                DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<5>(_vmax) * pow<5>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) *
                                                pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_vmax *
                                                Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period *
                                                pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_ns_first_term = DataType(2) * _N_s * _mu *((DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +
                                                (pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) *
                                                z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -((_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) *
                                                Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax *
                                                z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(5) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                                Math::sin(x) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -
                                                (pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) *
                                                pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_vmax) * z * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *
                                                (Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_vmax) *
                                                pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<2>(_vmax) *
                                                z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) *
                                                Math::sin(y) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) *
                                                _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_second_term = - DataType(2) * _N_s * _mu *((pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *
                                                  (Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                                  (pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                  pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax *
                                                  pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                  pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                  pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(6) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) *(Math::sin(y))) /
                                                  (DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * z *
                                                  Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                  pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                  pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *
                                                  (DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) *
                                                  pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                  pow<2>(Math::sin(y))) -(_period * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) /
                                                  DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(_vmax) * pow<3>(z) *
                                                  Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) *
                                                  pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * Math::cos(x) * Math::cos(y) *
                                                  Math::sin(x) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                  pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_vmax * Math::sin(x) *
                                                  (Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) *
                                                  pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                  pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<2>(Math::sin(y)) *
                                                  (DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                  Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                                  (_period * pow<2>(_vmax) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) *
                                                  _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                  Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_third_term = - DataType(2) * _N_s * _mu *((DataType(2) * pow<2>(_period) * _vmax * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                  pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *
                                                  (Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *
                                                  (Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -
                                                  (DataType(4) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                  pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                  pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) *(Math::sin(y))) /(DataType(3) *
                                                  _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * _vmax * z * Math::cos(x) *
                                                  Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                                  Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                  pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(_vmax) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_period * _vmax * Math::cos(x) *
                                                  (Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(_vmax) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax *
                                                  pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) +
                                                  DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<3>(_vmax) * pow<3>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * Math::sin(y) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) *
                                                  Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) *
                                                  _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType alp_term = _alpha * _vmax * z * Math::sin(x) * Math::sin(y);
            DataType konv_term = _konv ? (DataType(-0.5) * (_period * _rho * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(y)) * Math::sin(DataType(2) * x))) : DataType(0);
            DataType grap = DataType(2) * _mu * _period * _period * _vmax * z * Math::sin(x) * Math::sin(y);
            //calculate value[0] through (alpha*v + konv + grap - grtau)/rho, where grap is the gradient of the pressure, and grtau = grtau_d + grtau_np + grtau_ns_first_term
            val[0] = (alp_term + konv_term + grap - (grtau_d + grtau_np_first_term + grtau_np_second_term + grtau_np_third_term + grtau_ns_first_term + grtau_ns_second_term + grtau_ns_third_term))/_rho;
          }
          //second component
          {
            DataType grtau_d = DataType(2) * _mu * _period *_period * _vmax * z * Math::cos(x) * Math::cos(y);

            DataType grtau_np_first_term = DataType(2) * _N_p * _mu *((DataType(2) * pow<2>(_period) * _vmax * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * _vmax * z * Math::cos(y) * Math::sin(x) *
                                              (DataType(2) * _period * pow<4>(_vmax) * pow<4>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<4>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) /(DataType(15) *
                                              _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<5>(Math::cos(x)) *
                                              Math::cos(y) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                              (pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * Math::cos(x) * Math::cos(y) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)) *((_vmax *
                                              Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(8) * pow<4>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<3>(Math::cos(x)) * Math::cos(y) *
                                              pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                              (DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) /
                                              DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +
                                              (pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) *
                                              Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                              (DataType(4) * _period * _vmax * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) *
                                              Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *
                                              (pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * Math::cos(x) * Math::cos(y) * pow<3>(Math::sin(x)) *
                                              pow<3>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax *
                                              Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                              Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *
                                              (pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<5>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) *
                                              pow<2>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period *
                                              pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<3>(_period) *
                                              pow<5>(_vmax) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                              Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax *
                                              pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) +
                                              DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                              (DataType(4) * _period * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                              DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) *
                                              pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /
                                              pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_second_term = - DataType(2) * _N_p * _mu * ((DataType(16) * pow<4>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<5>(Math::cos(x)) * Math::cos(y) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) *
                                              pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                              pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) *
                                              pow<2>(_period) * _vmax * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(12) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *
                                              (Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                              pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(12) *
                                              pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<4>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<5>(Math::cos(x)) * Math::cos(y) * pow<4>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) *
                                              pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                              (DataType(12) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                              DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) *
                                              pow<5>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) *
                                              pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /
                                              pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) *
                                              pow<6>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) *
                                              pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /
                                              pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) *
                                              pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(_vmax) * z *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_third_term = DataType(2) * _N_p * _mu *((pow<2>(_period) * _vmax * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                              (_period * _vmax * z * Math::cos(x) * Math::sin(y) *(DataType(2) * _period * pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) - DataType(2) * _period * pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /
                                              (DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(pow<2>(_period) * _vmax * z * Math::cos(x) *
                                              Math::cos(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *(DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) -
                                              DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) *(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<5>(Math::cos(y)) * pow<2>(Math::sin(x))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * _vmax * z * Math::cos(x) * Math::cos(y) *
                                              (_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(6) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax *
                                              pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) +(DataType(8) * pow<2>(_period) * pow<5>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(y) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) /
                                              DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +
                                              (pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(_vmax) *
                                              pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                              (Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *
                                              (_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                              pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y))) *(DataType(2) * _period * pow<2>(_vmax) *
                                              pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                              (Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) *
                                              Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) * Math::sin(y) *
                                              ((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<5>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<4>(Math::cos(y)) * pow<2>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _period *
                                              pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                              (Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                              Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                              Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax *
                                              Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_ns_first_term = DataType(2) * _N_s * _mu *((DataType(2) * pow<2>(_period) * _vmax * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) /
                                              DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *
                                              (Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) *
                                              pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * _vmax * z * Math::cos(x) *
                                              Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) /
                                              DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<3>(_vmax) * pow<3>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                              Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              Math::cos(y) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                              Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_second_term = DataType(2) * _N_s * _mu *((pow<2>(_period) * _vmax * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *
                                              (Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) *
                                              pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -
                                              (DataType(2) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(6) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /
                                              (DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * z * Math::cos(y) *
                                              Math::sin(x) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                              Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *
                                              (DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                              pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                              (_period * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /
                                              (DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(_period * pow<2>(_vmax) * pow<3>(z) * Math::cos(x) *
                                              Math::sin(x) * pow<2>(Math::sin(y)) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(_vmax) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *
                                              ((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) *
                                              _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_vmax *
                                              Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period *
                                              pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *
                                              ((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<2>(Math::sin(y)) *(DataType(2) * pow<3>(_period) *
                                              pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *
                                              (Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_third_term =  - DataType(2) * _N_s * _mu *((DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +
                                              (pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax *
                                              z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) *
                                              pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_vmax) * z * Math::cos(x) * Math::cos(y) *
                                              Math::sin(x) * Math::sin(y) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(5) * pow<2>(_period) * pow<3>(_vmax) * pow<3>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /
                                              (DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) *
                                              pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) *
                                              pow<3>(_vmax) * pow<4>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) *
                                              pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y))) -(pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(_vmax) * z *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) *
                                              _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType alp_term = - _alpha * _vmax * z * Math::cos(x) * Math::cos(y);
            DataType konv_term = _konv ? (DataType(0.5) * (_period * _rho * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * Math::sin(DataType(2) * y))) : DataType(0);
            DataType grap = DataType(-2) * _mu * _period * _period * _vmax * z * Math::cos(x) * Math::cos(y);
            //calculate value[0] through (grap - grtau)/rho, where grap is the gradient of the pressure, and grtau = grtau_d + grtau_np + grtau_ns_first_term
            val[1] = (alp_term + konv_term + grap - (grtau_d + grtau_np_first_term + grtau_np_second_term + grtau_np_third_term + grtau_ns_first_term + grtau_ns_second_term + grtau_ns_third_term))/_rho;
          }

          //calculate third component
          {
            DataType grtau_d = DataType(4) * _mu *((_period * _vmax * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2)) - DataType(4) * _mu * _period * _vmax * Math::cos(x) *(Math::sin(y));

            DataType grtau_np_first_term = DataType(2) * _N_p * _mu *((DataType(2) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * Math::cos(y) *
                                                pow<3>(Math::sin(y)) - DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) *(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) /
                                                DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) -
                                                DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(2) * pow<5>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<4>(Math::cos(y)) *(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(3) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(8) * pow<5>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) *
                                                pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                                DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<3>(_period) * pow<5>(_vmax) *
                                                pow<6>(z) * pow<3>(Math::cos(x)) * pow<4>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(10) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *
                                                ((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x)
                                                *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<4>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<5>(Math::cos(x)) * Math::cos(y) * pow<4>(Math::sin(y)) *(DataType(2) * _period *
                                                pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                (Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) -
                                                DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                                (DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +
                                                (pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_second_term = - DataType(2) * _N_p * _mu *((DataType(2) * _period * _vmax * Math::cos(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_period) * pow<4>(_vmax) * pow<8>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period * _vmax * Math::cos(x) * Math::sin(y) *(_epsilon +
                                                pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period * _vmax * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(16) * pow<5>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<5>(Math::sin(y))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) *
                                                pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(2) * pow<5>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<5>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) *
                                                _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(14) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) *
                                                Math::sin(x) * pow<4>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(14) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) *
                                                ((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(14) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *
                                                (_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period *
                                                _vmax * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) *
                                                z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) -(DataType(4) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_period) * pow<4>(_vmax) * pow<8>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<3>(_period) * pow<4>(_vmax) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax *
                                                pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<3>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                                pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<3>(_period) * pow<4>(_vmax) * pow<7>(z) * pow<3>(Math::cos(x)) * Math::sin(x) *
                                                pow<4>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) *
                                                z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))));

            DataType grtau_np_third_term = - DataType(2) * _N_p * _mu *((DataType(2) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) *
                                                Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) - DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *((_period * _vmax *
                                                Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon +
                                                pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<5>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * Math::cos(x) * pow<4>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon +
                                                pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(3) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) *
                                                pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(8) * pow<5>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) *((_vmax *
                                                Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(10) * pow<3>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) *
                                                Math::sin(x) * pow<3>(Math::sin(y)) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y)) *
                                                (DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                                pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) *
                                                _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(6) * pow<3>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) *
                                                Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<4>(_period) * pow<5>(_vmax) * pow<8>(z) * pow<4>(Math::cos(x)) * Math::sin(x) *
                                                pow<5>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) *
                                                pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /
                                                pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) *
                                                pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) -
                                                DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *
                                                (Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) *
                                                pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_ns_first_term = DataType(2) * _N_s * _mu *(((_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /
                                                (DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +((_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y))) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *
                                                (Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) *
                                                Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +((_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +
                                                (pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) *
                                                pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * Math::cos(x) * pow<2>(Math::sin(x)) *
                                                pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(_period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) *
                                                ((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::sin(x) *
                                                pow<3>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period *
                                                pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *
                                                (Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_vmax) *
                                                pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(3) *
                                                _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * Math::sin(y) *
                                                (DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) *
                                                Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) +(pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) /
                                                DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) *
                                                Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))));

            DataType grtau_ns_second_term = - DataType(2) * _N_s * _mu *((DataType(4) * _period * _vmax * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(16) * pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::sin(y))) /(DataType(3) *
                                                _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                                (DataType(2) * pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * Math::cos(x) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(6) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax *
                                                pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) +(DataType(6) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                                DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(4) * _period * _vmax *
                                                z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) *
                                                Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period *
                                                pow<2>(_vmax) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(_vmax) * z *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(_vmax) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_third_term =  DataType(2) * _N_s * _mu *(((_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +((_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *((_period * _vmax * Math::cos(x) *(Math::sin(y))) /
                                                DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) +(((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) *
                                                _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                                (pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * Math::cos(x) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_vmax *
                                                Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon +
                                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<3>(_vmax) * pow<4>(z) * Math::cos(x) *
                                                pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) -(_period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * Math::sin(x) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /
                                                (DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(_vmax) * pow<2>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *
                                                ((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((_period * _vmax * Math::cos(y) *(Math::sin(x))) /
                                                DataType(2) +(pow<3>(_period) * _vmax * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) *
                                                Math::sin(y) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(_vmax) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType alp_term = - _alpha * _period * _vmax * pow<2>(z) * Math::cos(x) * Math::sin(y);
            DataType konv_term = _konv ? (pow<2>(_period) * _rho * pow<2>(_vmax) * pow<3>(z) * (pow<2>(Math::cos(x)) + pow<2>(Math::sin(y)))) : DataType(0);
            DataType grap = DataType(-2) * _mu * _period * _vmax * Math::cos(x) * Math::sin(y);
                        //calculate value[0] through (grap - grtau)/rho, where grap is the gradient of the pressure, and grtau = grtau_d + grtau_np + grtau_ns_first_term
            val[2] = (alp_term + konv_term + grap - (grtau_d + grtau_np_first_term + grtau_np_second_term + grtau_np_third_term + grtau_ns_first_term + grtau_ns_second_term + grtau_ns_third_term))/_rho;
          }
          //all done return val
          return val;
        }
      };
    }; // class RHS

    class Robin :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;

    public:
      explicit Robin(const FullTensor& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _epsilon;

      public:
        explicit Evaluator(const Robin& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];

          ValueType val;
          //Again build the different components one after the other
          {
            DataType h_d = DataType(4) * _mu * _period * _vmax * z * Math::cos(x) * Math::sin(y); //because of stupid choice of p it 2 times the normal h in this component...... which is zero
            DataType h_np = - DataType(2) * _N_p * _mu *((DataType(2) * _period * _vmax * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period * _vmax * z * Math::cos(x) * Math::sin(y) *(DataType(3) *
                              _epsilon + pow<4>(_vmax) * pow<4>(z) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) *
                              pow<2>(Math::sin(y)))) -(_period * _vmax * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<5>(_vmax) * pow<5>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *
                              pow<4>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                              (DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) *
                              _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) *
                              Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType h_ns = - DataType(2) * _N_s * _mu *((DataType(2) * _period * pow<2>(_vmax) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                              DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period *
                              _vmax * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<3>(_vmax) * pow<3>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            val[0] = (h_d + h_np + h_ns)/_rho;
          }

          {
            DataType h_d = DataType(2) * _mu * _period * _vmax * z * Math::cos(y) * Math::sin(x);

            DataType h_np = DataType(2) * _N_p * _mu *((DataType(2) * _period * _vmax * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<4>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) *
                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period * pow<5>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) *
                                pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period *
                                pow<5>(_vmax) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<3>(_period) * pow<5>(_vmax) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *
                                ((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) *
                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(_vmax) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_vmax * Math::sin(x) *
                                (Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType h_ns = DataType(2) * _N_s * _mu *((_period * _vmax * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /
                                (DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(_vmax) * pow<3>(z) *
                                Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *
                                ((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * pow<3>(_vmax) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) *
                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            val[1] = (h_d + h_np + h_ns)/_rho;
          }

          {
            DataType h_d = DataType(2) * _mu * ((_vmax * Math::sin(x) * Math::sin(y))/DataType(2) + (_period * _period * _vmax * z * z * Math::sin(x) * Math::sin(y))/DataType(2));

            DataType h_np = DataType(2) * _N_p * _mu *((DataType(2) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) *
                              pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) *
                              pow<5>(_vmax) * pow<8>(z) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) *
                              pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(_vmax) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                              pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) *
                              pow<2>(_period) * pow<4>(_vmax) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /
                              (DataType(15) * _epsilon + pow<2>(pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType h_ns = DataType(2) * _N_s * _mu *((((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /(DataType(3) *
                              _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +((_epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) *
                              pow<2>(Math::sin(y))) *((_vmax * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_vmax) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) *
                              Math::sin(y) *((_vmax * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * _vmax * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) *
                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) *(Math::sin(y))) /(DataType(3) * _epsilon +
                              pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(_vmax) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            val[2] = (h_d + h_np + h_ns)/_rho;
          }

          return val;
        }
      };
    }; // class Robin

    class Orient :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<6> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon;
    public:
      explicit Orient(const FullTensor& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon;
      public:
        explicit Evaluator(const Orient& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];

          const DataType v1 = _vmax * Math::sin(x) * Math::sin(y) * z;
          const DataType v2 = DataType(-1) * _vmax * Math::cos(x) * Math::cos(y)*z;
          const DataType v3 = DataType(-1) * _vmax * _period * Math::cos(x) * Math::sin(y) * z * z;
          const DataType temp_norm = DataType(3) * _epsilon + v1 * v1 + v2 * v2 + v3 * v3;




          //the matrix will be represented by a vector by A_ij = f(j*dim + i)
          ValueType val;
//           val[0] = (v1 * v1 + _epsilon) / temp_norm;
//           val[1] = (v1 * v2) / temp_norm;
//           val[2] = (v1 * v3) / temp_norm;
//           val[3] = (v1 * v2) / temp_norm;
//           val[4] = (v2 * v2 + _epsilon) /temp_norm;
//           val[5] = (v2 * v3) / temp_norm;
//           val[6] = (v1 * v3) / temp_norm;
//           val[7] = (v2 * v3) / temp_norm;
//           val[8] = (v3 * v3 + _epsilon) / temp_norm;

          val[0] = (v1 * v1 + _epsilon) / temp_norm;
          val[1] = (v2 * v2 + _epsilon) /temp_norm;
          val[2] = (v3 * v3 + _epsilon) / temp_norm;
          val[3] = (v1 * v2) / temp_norm;
          val[4] = (v1 * v3) / temp_norm;
          val[5] = (v2 * v3) / temp_norm;

          return val;
        }
      };
    }; // class Orient

    class Fourth_Moment :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<15> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon;
    public:
      explicit Fourth_Moment(const FullTensor& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon;
      public:
        explicit Evaluator(const Fourth_Moment& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];

          const DataType v_0 = _vmax * Math::sin(x) * Math::sin(y) * z;
          const DataType v_1 = DataType(-1) * _vmax * Math::cos(x) * Math::cos(y)*z;
          const DataType v_2 = DataType(-1) * _vmax * _period * Math::cos(x) * Math::sin(y) * z * z;

          const DataType temp_norm = DataType(15) * _epsilon + pow<2>(v_0 * v_0 + v_1 * v_1 + v_2 * v_2);

          ValueType val;

          //symmetric version
          val[0] = (v_0 * v_0 * v_0 * v_0 + DataType(3) * _epsilon) / temp_norm;
          val[1] = (v_1 * v_1 * v_1 * v_1 + DataType(3) * _epsilon) / temp_norm;
          val[2] = (v_2 * v_2 * v_2 * v_2 + DataType(3) * _epsilon) / temp_norm;
          val[3] = (v_0 * v_0 * v_0 * v_1) / temp_norm;
          val[4] = (v_0 * v_0 * v_0 * v_2) / temp_norm;
          val[5] = (v_1 * v_1 * v_1 * v_0) / temp_norm;
          val[6] = (v_1 * v_1 * v_1 * v_2) / temp_norm;
          val[7] = (v_2 * v_2 * v_2 * v_0) / temp_norm;
          val[8] = (v_2 * v_2 * v_2 * v_1) / temp_norm;
          val[9] = (v_0 * v_0 * v_1 * v_1 + _epsilon) / temp_norm;
          val[10] = (v_0 * v_0 * v_2 * v_2 + _epsilon) / temp_norm;
          val[11] = (v_1 * v_1 * v_2 * v_2 + _epsilon) / temp_norm;
          val[12] = (v_0 * v_0 * v_1 * v_2) / temp_norm;
          val[13] = (v_1 * v_1 * v_0 * v_2) / temp_norm;
          val[14] = (v_2 * v_2 * v_0 * v_1) / temp_norm;

          return val;

        }
      };
    }; // class Fourth_Moment


  }; //class FullTensor 3D spezialization






  template<int dim_>
  class FullTensorTimeExpo;


   /////////////////////////////////////////////////
  // Bubble function constructed to meet zero Neumann boundary on right side for space dependent 4th order orientation, Christophs pyhsical formulation...
  ////////////////////////////////////////////////////
  /**
   * \brief Divergence-free 2-d bubble function class
   *
   * This class initialises all necessary parameters and provides an interface to access generalized test functions.
   */
  template<>
  class FullTensorTimeExpo<2>
  {
  public:
    class Velo;
    class Pres;
    class RHS;
    class Robin;
    class Orient;
    class Fourth_Moment;

    DataType _vmax;
    DataType _mu;
    DataType _rho;
    DataType _N_s;
    DataType _N_p;
    DataType _period;
    DataType _epsilon;
    DataType _alpha;
    bool _konv;

    explicit FullTensorTimeExpo(DataType vmax, DataType mu, DataType rho, DataType N_s, DataType N_p, DataType epsilon = DataType(1), DataType period = DataType(2) * DataType(Math::pi<DataType>()), DataType alpha = DataType(0.), bool konv = false) :
      _vmax(vmax),
      _mu(mu),
      _rho(rho),
      _N_s(N_s),
      _N_p(N_p),
      _period(period),
      _epsilon(epsilon),
      _alpha(alpha),
      _konv(konv)
      {}

    class Velo :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;
      static constexpr bool can_grad = true;

    protected:
      const DataType _vmax;
      const DataType _period;
      const DataType _alpha;

    private:
      DataType _t;


    public:
      explicit Velo(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _period(func._period),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;

        const DataType _vmax, _period, _alpha, _t;

      public:
        explicit Evaluator(const Velo& function) :
        _vmax(function._vmax),
        _period(function._period),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType exp = Math::exp(_alpha * _t);

          ValueType val;
          val[0] = exp * _vmax * Math::cos(x) * Math::sin(y);
          val[1] = DataType(-1) * exp * _vmax * Math::sin(x) * Math::cos(y);
          return val;
        }

        GradientType gradient(const PointType& point)
        {
          //with gradient we mean the Jacobian Matrix here
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType exp = Math::exp(_alpha * _t);

          GradientType grad;
          grad[0][0] = - exp * _vmax * _period * Math::sin(x) * Math::sin(y);
          grad[0][1] = exp * _vmax * _period * Math::cos(x) * Math::cos(y);
          grad[1][0] = - exp * _vmax * _period * Math::cos(x) * Math::cos(y);
          grad[1][1] = exp * _vmax * _period * Math::sin(x) * Math::sin(y);

          return grad;
        }

      };
    }; //class Velo

    class Pres :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Scalar ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;
      DataType _alpha;

    private:
      DataType _t;

    public:
      explicit Pres(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_p, _period, _epsilon, _alpha, _t;

      public:
        explicit Evaluator(const Pres& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];

          ValueType val;
          val = -DataType(2) * _mu * _period * _vmax * Math::sin(x) * Math::sin(y);
          return val;
        }

      };
    }; // class Pres

    /**
     * \brief The assoicated right side for an orientation matrix
     */
    class RHS :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;
      DataType _alpha;
      bool _konv;

    private:
      DataType _t;

    public:
      explicit RHS(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _konv(func._konv),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      void print_variables() const
      {
        std::cout << "RHS print:\n";
        std::cout << "Vmax : " << _vmax << "\n";
        std::cout << "Mu : " << _mu << "\n";
        std::cout << "Rho : " << _rho << "\n";
        std::cout << "Ns : " << _N_s << "\n";
        std::cout << "Np : " << _N_p << "\n";
        std::cout << "period : " << _period << "\n";
        std::cout << "epsilon : " << _epsilon << "\n";
        std::cout << "Alpha : " << _alpha << "\n";
        std::cout << "Konv : " << _konv << "\n";
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _epsilon, _alpha, _t;
        const bool _konv;

      public:
        explicit Evaluator(const RHS& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t),
        _konv(function._konv)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType exp = Math::exp(_alpha * _t);
//           const DataType vmax_exp = _vmax * exp;

          ValueType val;
          {
            //first val component
            const DataType taugrad = - DataType(2) * _mu *(_N_p *((_period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(2) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) *
                                      pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                      pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                      pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * _vmax * Math::cos(x) * Math::sin(y) *(_epsilon + pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /
                                      (DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax * Math::cos(x) *
                                      Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * Math::sin(x) * Math::sin(y) *(_epsilon + pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(4) * _period *
                                      pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                      pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) +
                                      DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y))) *(DataType(4) * _period *
                                      pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                      pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) +
                                      DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) + pow<2>(_period) * _vmax * Math::cos(x) * Math::sin(y) -(DataType(4) * _N_s * pow<2>(_period) * pow<3>(_vmax) * Math::cos(x) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)))
                                      /(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) +(DataType(2) * _N_s * pow<2>(_period) * _vmax * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::sin(y)))) /(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) -(DataType(2) * _N_s * _period * _vmax * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * Math::sin(x) - DataType(2) * _period * pow<2>(_vmax) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)))) - DataType(2) * _N_p * _mu *((pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) +
                                      pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                      pow<3>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(DataType(2) *
                                      pow<2>(_period) * pow<5>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) *(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                      pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) +
                                      pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * pow<5>(_vmax) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<2>(Math::sin(y)) *
                                      (DataType(4) * _period * pow<4>(_vmax) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                      pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) *
                                      pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * pow<5>(_vmax) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(4) * _period * pow<4>(_vmax) *
                                      pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) +
                                      DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                      pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))));

            const DataType pgrad = DataType(-2) * _mu * _period * _period * _vmax * Math::cos(x) * Math::sin(y);

            const DataType konvektion = _konv ? (DataType(-0.5) * _period * _rho * _vmax * _vmax * Math::sin(DataType(2) * x)) : DataType(0);

            const DataType react = _vmax * Math::cos(x) * Math::sin(y);

            val[0] = (_alpha * exp * react + exp * exp * konvektion + pgrad - exp * taugrad)/_rho;
          }

          {
            //second component
            const DataType taugrad = DataType(2) * _mu *(_N_p *((_period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(2) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) - DataType(2) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                        pow<2>(Math::sin(x)) *(Math::sin(y)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *
                                        pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(y)) * pow<5>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) *
                                        pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * _vmax * Math::cos(y) * Math::sin(x) *(_epsilon + pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(8) *
                                        _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * _vmax * Math::cos(y) *
                                        Math::sin(x) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                        pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * _vmax * Math::sin(x) * Math::sin(y) *(_epsilon + pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(4) * _period * pow<4>(_vmax) *
                                        pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) +
                                        DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                        pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * _vmax * Math::sin(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x))) *(DataType(4) * _period * pow<4>(_vmax) *
                                        pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * Math::sin(y) - DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) +
                                        DataType(4) * _period * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) *
                                        pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) + pow<2>(_period) * _vmax * Math::cos(y) * Math::sin(x) -(DataType(4) * _N_s * pow<2>(_period) * pow<3>(_vmax) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y))) /
                                        (DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) +(DataType(2) * _N_s * pow<2>(_period) * _vmax * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))))
                                        /(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) -(DataType(2) * _N_s * _period * _vmax * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_vmax) * pow<2>(Math::cos(y)) *
                                        pow<2>(Math::sin(x))) *(DataType(2) * _period * pow<2>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(_vmax) * Math::cos(y) * pow<2>(Math::sin(x)) *(Math::sin(y)))) / pow<2>(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::cos(x)) *
                                        pow<2>(Math::sin(y)) + pow<2>(_vmax) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)))) + DataType(2) * _N_p * _mu *((pow<2>(_period) * pow<5>(_vmax) * pow<3>(Math::cos(y)) * pow<5>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) +
                                        pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<5>(_vmax) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                        pow<2>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) *
                                        pow<5>(_vmax) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) *
                                        pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(DataType(3) * pow<2>(_period) * pow<5>(_vmax) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) *
                                        pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_period * pow<5>(_vmax) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(4) *
                                        _period * pow<4>(_vmax) * Math::cos(x) * pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                                        pow<2>(Math::sin(y)) + DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) +
                                        DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) -(_period * pow<5>(_vmax) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) *
                                        pow<4>(Math::cos(y)) * pow<3>(Math::sin(x)) - DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) - DataType(4) * _period * pow<4>(_vmax) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                        DataType(4) * _period * pow<4>(_vmax) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) / pow<2>(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(_vmax) *
                                        pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))));

            const DataType pgrad = DataType(-2) * _mu * _period * _period * _vmax * Math::cos(y) * Math::sin(x);

            const DataType konvektion = _konv ? (DataType(-0.5)* _period * _rho * _vmax * _vmax * Math::sin(DataType(2)*y)) : DataType(0);

            const DataType react = - _vmax * Math::sin(x) * Math::cos(y);

            val[1] = (_alpha * exp * react + exp * exp * konvektion + pgrad - exp * taugrad)/_rho;
          }

          return val;
        }
      };
    }; // class RHS


    class Robin :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;
      DataType _alpha;

    private:
      DataType _t;

    public:
      explicit Robin(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _epsilon, _alpha, _t;

      public:
        explicit Evaluator(const Robin& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType exp = Math::exp(_alpha * _t);
//           const DataType vmax_exp = _vmax * exp;

          ValueType val;
          /*val[0] = -(DataType(2) * _mu * _period * vmax_exp * Math::sin(x) * Math::sin(y) *(DataType(4) * _N_p * pow<2>(_epsilon) + DataType(16) * _N_s * pow<2>(_epsilon) + _N_p * pow<6>(vmax_exp) * pow<6>(Math::cos(x)) * pow<6>(Math::sin(y)) + DataType(2) * _N_s * pow<6>(vmax_exp) * pow<6>(Math::cos(x)) *
                    pow<6>(Math::sin(y)) + DataType(2) * _N_p * _epsilon * pow<2>(vmax_exp) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + DataType(2) * _N_p * _epsilon * pow<4>(vmax_exp) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + DataType(16) * _N_s * _epsilon * pow<2>(vmax_exp) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) +
                    DataType(2) * _N_s * _epsilon * pow<4>(vmax_exp) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y))) + DataType(2) * _mu * _period * vmax_exp * pow<4>(Math::cos(y)) * pow<5>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _N_s * _epsilon * pow<4>(vmax_exp) - _N_p * pow<6>(vmax_exp) *
                    pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + DataType(2) * _N_s * pow<6>(vmax_exp) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) + DataType(2) * _mu * _period * vmax_exp * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _N_p * _epsilon * pow<2>(vmax_exp) +
                    DataType(4) * _N_s * pow<6>(vmax_exp) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) - DataType(2) * _N_p * _epsilon * pow<4>(vmax_exp) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + DataType(4) * _N_s * _epsilon * pow<4>(vmax_exp) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /((DataType(2) * _epsilon +
                    pow<2>(vmax_exp) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) + pow<2>(vmax_exp) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x))) *(DataType(8) * _epsilon + pow<4>(vmax_exp) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + pow<4>(vmax_exp) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(2) * pow<4>(vmax_exp) * pow<2>(Math::cos(x)) *
                    pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))));
          val[1] = (DataType(2) * _N_p * _mu * _period * pow<5>(vmax_exp) * Math::cos(x) * Math::cos(y) *(pow<2>(Math::cos(x)) - DataType(1)) *(pow<2>(Math::cos(y)) - DataType(1)) *(pow<2>(Math::cos(x)) - pow<2>(Math::cos(y)))) /(DataType(8) * _epsilon + pow<4>(vmax_exp) * pow<4>(Math::cos(x)) + pow<4>(vmax_exp) * pow<4>(Math::cos(y)) +
                    DataType(2) * pow<4>(vmax_exp) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) - DataType(4) * pow<4>(vmax_exp) * pow<2>(Math::cos(x)) * pow<4>(Math::cos(y)) - DataType(4) * pow<4>(vmax_exp) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(4) * pow<4>(vmax_exp) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)));
          */
          val[0] = DataType(2) * _mu * _period * _vmax * Math::sin(x) * Math::sin(y) - DataType(2) * _mu *(exp * _period * _vmax * Math::sin(x) * Math::sin(y) +(DataType(2) * _N_s * exp * _period * _vmax * Math::sin(x) * Math::sin(y) *
                    (_epsilon - pow<2>(_vmax) * pow<2>(Math::sin(y)) *(pow<2>(Math::sin(x)) - DataType(1)))) /(DataType(2) * _epsilon + pow<2>(_vmax) * pow<2>(Math::sin(x)) + pow<2>(_vmax) * pow<2>(Math::sin(y)) - DataType(2) * pow<2>(_vmax) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) +(_N_p * exp * _period * _vmax *
                    Math::sin(x) * Math::sin(y) *(DataType(2) * _epsilon + pow<4>(_vmax) * pow<4>(Math::sin(y)) - pow<4>(_vmax) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) - pow<4>(_vmax) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) + pow<4>(_vmax) * pow<4>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(8) *
                    _epsilon + pow<4>(_vmax) * pow<4>(Math::sin(x)) + pow<4>(_vmax) * pow<4>(Math::sin(y)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) - DataType(4) * pow<4>(_vmax) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) - DataType(4) * pow<4>(_vmax) * pow<4>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                    DataType(4) * pow<4>(_vmax) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y))));

          val[1] = (DataType(2) * _N_p * exp * _mu * _period * pow<5>(_vmax) * Math::cos(x) * Math::cos(y) *(pow<2>(Math::cos(x)) - DataType(1)) *(pow<2>(Math::cos(y)) - DataType(1)) *(pow<2>(Math::cos(x)) - pow<2>(Math::cos(y)))) /(DataType(8) * _epsilon + pow<4>(_vmax) * pow<4>(Math::cos(x)) + pow<4>(_vmax) *
                    pow<4>(Math::cos(y)) + DataType(2) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) - DataType(4) * pow<4>(_vmax) * pow<2>(Math::cos(x)) * pow<4>(Math::cos(y)) - DataType(4) * pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(4) * pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)));

          return val;
        }
      };
    }; // class Robin

    class Orient :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon, _alpha;

    private:
      DataType _t;

    public:
      explicit Orient(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon, _alpha, _t;
      public:
        explicit Evaluator(const Orient& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
           const DataType x = _period * point[0];
           const DataType y = _period * point[1];
//            const DataType exp = Math::exp(_alpha * _t);
//            const DataType vmax_exp = _vmax * exp;
           const DataType temp_norm = DataType(2) * _epsilon + pow<2>(_vmax) * (pow<2>(Math::cos(x) * Math::sin(y)) + pow<2>(Math::cos(y) * Math::sin(x)));

          ValueType val;
          val[0] = (_vmax * Math::cos(x) * Math::sin(y) * _vmax * Math::cos(x) * Math::sin(y) + _epsilon) / temp_norm;
          val[2] = _vmax * Math::cos(x) * Math::sin(y) * DataType(-1) * _vmax * Math::sin(x) * Math::cos(y) / temp_norm;
//           val[2] = val[1];
          val[1] = (_vmax * Math::sin(x) * Math::cos(y) * _vmax * Math::sin(x) * Math::cos(y) + _epsilon) / temp_norm;

          return val;
        }
      };
    }; // class Orient

    class Fourth_Moment :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<5> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon, _alpha;
    private:
      DataType _t;
    public:
      explicit Fourth_Moment(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon, _alpha, _t;
      public:
        explicit Evaluator(const Fourth_Moment& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
//           const DataType exp = Math::exp(_alpha * _t);
//           const DataType vmax_exp = _vmax * exp;
          const DataType temp_norm = DataType(8) * _epsilon + pow<4>(_vmax) * (pow<4>(Math::cos(x) * Math::sin(y)) + pow<4>(Math::cos(y) * Math::sin(x))
                                     + DataType(2)*(pow<2>(Math::cos(x) * Math::sin(y) * Math::sin(x) * Math::cos(y))));


          //we want to build the tensor A_(klmn} = v_k * v_l * v_m * v_n + \delta_{klmn}*eps
          //we represent this by a vector of the form:
          // A_{klmn} = val(d^3*k + d^2*l + d*m + n)
          //we have 2 different combs for all same, which are in val[0] and val[15]
          // 2 different combs for one different
          ValueType val;
//           val[0] = (pow<4>(vmax_exp) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + DataType(3) * _epsilon) / temp_norm;
//           val[1] = (pow<4>(vmax_exp) * pow<3>(Math::cos(x) * Math::sin(y)) * pow<1>(-Math::sin(x) * Math::cos(y))) / temp_norm;
//           val[2] = val[1];
//           val[3] = (pow<4>(vmax_exp) * pow<2>(Math::cos(x) * Math::sin(y)) * pow<2>(-Math::sin(x) * Math::cos(y)) + _epsilon) / temp_norm;
//           val[4] = val[1];
//           val[5] = val[3];
//           val[6] = val[3];
//           val[7] = (pow<4>(vmax_exp) * pow<1>(Math::cos(x) * Math::sin(y)) * pow<3>(-Math::sin(x) * Math::cos(y))) / temp_norm;
//           val[8] = val[1];
//           val[9] = val[3];
//           val[10] = val[3];
//           val[11] = val[7];
//           val[12] = val[3];
//           val[13] = val[7];
//           val[14] = val[7];
//           val[15] = (pow<4>(vmax_exp) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(3) * _epsilon) / temp_norm;

          //for this we use the definition of the fourth_order_contraction operator, see tiny_algebra.hpp
          val[0] = (pow<4>(_vmax) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)) + DataType(3) * _epsilon) / temp_norm;
          val[1] = (pow<4>(_vmax) * pow<4>(Math::cos(y)) * pow<4>(Math::sin(x)) + DataType(3) * _epsilon) / temp_norm;
          val[2] = (pow<4>(_vmax) * pow<3>(Math::cos(x) * Math::sin(y)) * pow<1>(-Math::sin(x) * Math::cos(y))) / temp_norm;
          val[3] = (pow<4>(_vmax) * pow<1>(Math::cos(x) * Math::sin(y)) * pow<3>(-Math::sin(x) * Math::cos(y))) / temp_norm;
          val[4] = (pow<4>(_vmax) * pow<2>(Math::cos(x) * Math::sin(y)) * pow<2>(-Math::sin(x) * Math::cos(y)) + _epsilon) / temp_norm;

          return val;
        }
      };
    }; // class Fourth_Moment

  };//class FullTensorTimeExpo 2D spezialization









      /////////////////////////////////////////////////
  // A non-trivial divergence free velocity field and all associated functions with pyhsical, regularized tensors
  ////////////////////////////////////////////////////
  /**
   * \brief Divergence-free 3-D bubble function class
   *
   * This class initialises all necessary parameters and provides an interface to access generalized test functions.
   */
  template<>
  class FullTensorTimeExpo<3>
  {
  public:
    class Velo;
    class Pres;
    class RHS;
    class Robin;
    class Orient;
    class Fourth_Moment;

    DataType _vmax;
    DataType _mu;
    DataType _rho;
    DataType _N_s;
    DataType _N_p;
    DataType _period;
    DataType _epsilon;
    DataType _alpha;
    bool _konv;

    explicit FullTensorTimeExpo(DataType vmax, DataType mu, DataType rho, DataType N_s, DataType N_p, DataType epsilon = DataType(1), DataType period = DataType(2) * DataType(Math::pi<DataType>()), DataType alpha = DataType(2), bool konv = false) :
      _vmax(vmax),
      _mu(mu),
      _rho(rho),
      _N_s(N_s),
      _N_p(N_p),
      _period(period),
      _epsilon(epsilon),
      _alpha(alpha),
      _konv(konv)
      {}

    class Velo :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;
      static constexpr bool can_grad = true;

    private:
      //the time
      DataType _t;

    protected:
      const DataType _vmax;
      const DataType _period;
      const DataType _alpha;


    public:
      explicit Velo(const FullTensorTimeExpo& func) :
      _t(DataType(0.)),
      _vmax(func._vmax),
      _period(func._period),
      _alpha(func._alpha)
      {
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;

        const DataType _vmax, _period, _alpha, _t;

      public:
        explicit Evaluator(const Velo& function) :
        _vmax(function._vmax),
        _period(function._period),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];
          const DataType exp = Math::exp(_alpha * _t);

          ValueType val;
          val[0] = exp * _vmax * Math::sin(x) * Math::sin(y) * z;
          val[1] = exp * DataType(-1) * _vmax * Math::cos(x) * Math::cos(y)*z;
          val[2] = exp * DataType(-1) * _vmax * _period * Math::cos(x) * Math::sin(y) * z * z;
          return val;
        }

        GradientType gradient(const PointType& point)
        {
          //with gradient we mean the Jacobian Matrix here
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];
          const DataType exp = Math::exp(_alpha * _t);

          GradientType grad;
          grad[0][0] =  exp * _vmax * _period * Math::cos(x) * Math::sin(y) * z;
          grad[0][1] = exp * _vmax * _period * Math::sin(x) * Math::cos(y) * z;
          grad[0][2] = exp * _vmax * Math::sin(x) * Math::sin(y);
          grad[1][0] = exp * _vmax * _period * Math::sin(x) * Math::cos(y) * z;
          grad[1][1] = exp * _vmax * _period * Math::cos(x) * Math::sin(y) * z;
          grad[1][2] = - exp * _vmax * Math::cos(x) * Math::cos(y);
          grad[2][0] = exp * _vmax * _period * _period * Math::sin(x) * Math::sin(y) * z * z;
          grad[2][1] = - exp * _vmax * _period * _period * Math::cos(x) * Math::cos(y) * z * z;
          grad[2][2] = DataType(-2) * exp *  _vmax * _period * Math::cos(x) * Math::sin(y) * z;

          return grad;
        }


      };

      void set_time(const DataType t)
      {
        _t = t;
      }
    }; //class Velo

    class Pres :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Scalar ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _period;
      DataType _alpha;

    private:
      DataType _t;

    public:
      explicit Pres(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _period(func._period),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _period, _alpha, _t;

      public:
        explicit Evaluator(const Pres& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _period(function._period),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];
          const DataType exp = Math::exp(_alpha * _t);

          ValueType val;
          val = DataType(-2) * exp * _mu * _period * _vmax * z * Math::cos(x) * Math::sin(y);
          return val;
        }

      };
    }; // class Pres

    /**
     * \brief The assoicated right side for an orientation matrix of typedef
     * A = (lambda1 0; 0 lambda2)
     */
    /**
     * \brief The associated pressure function to SteadyBubble2d
     */
    class RHS :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;
      DataType _alpha;
      bool _konv;

    private:
      DataType _t;

    public:
      explicit RHS(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _konv(func._konv),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      void print_variables() const
      {
        std::cout << "RHS print:\n";
        std::cout << "Vmax : " << _vmax << "\n";
        std::cout << "Mu : " << _mu << "\n";
        std::cout << "Rho : " << _rho << "\n";
        std::cout << "Ns : " << _N_s << "\n";
        std::cout << "Np : " << _N_p << "\n";
        std::cout << "period : " << _period << "\n";
        std::cout << "epsilon : " << _epsilon << "\n";
        std::cout << "Alpha : " << _alpha << "\n";
        std::cout << "Konv : " << _konv << "\n";
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _epsilon, _alpha, _t;
        const bool _konv;

      public:
        explicit Evaluator(const RHS& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t),
        _konv(function._konv)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];
          const DataType exp = Math::exp(_alpha * _t);
          const DataType vmax_exp = _vmax * exp;

          //yes this is highly optimisable... if you want, have fun
          ValueType val;
          //RHS for first component
          {
            DataType grtau_d = - DataType(2) * _mu * pow<2>(_period) * vmax_exp * z * Math::sin(x) *(Math::sin(y));
            DataType grtau_np_first_term = DataType(2) * _N_p * _mu *((DataType(2) * pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp *
                                                pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) *
                                                z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(16) * pow<4>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(12) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *
                                                ((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                                Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(12) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<4>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(12) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) *
                                                vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *
                                                pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) *
                                                pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) *
                                                pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_second_term = - DataType(2) * _N_p * _mu *((DataType(2) * _period * vmax_exp * z * Math::cos(y) * Math::sin(x) *(DataType(2) * _period * pow<4>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) - DataType(2) * _period * pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *
                                                (_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)) *
                                                Math::sin(x) *(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(2) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(3) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) *
                                                pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(8) * pow<4>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) *
                                                _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) *
                                                Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * vmax_exp * z * Math::cos(y) * Math::sin(x) *(_epsilon +
                                                pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp *
                                                pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *
                                                (Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *
                                                pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<5>(vmax_exp) * pow<5>(z) *
                                                pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<5>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) +
                                                DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<3>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y)) *(DataType(2) * _period * pow<2>(vmax_exp) *
                                                pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                (Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *
                                                ((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) *
                                                pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_third_term = - DataType(2) * _N_p * _mu *((_period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(DataType(2) * _period * pow<4>(vmax_exp) * pow<4>(z) * Math::cos(x) *
                                                pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<4>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) *
                                                pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y)) *
                                                ((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y)))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * vmax_exp * z *
                                                Math::cos(x) * Math::sin(y) *(DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) - DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)))) /(DataType(15) *
                                                _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) *
                                                pow<3>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                                (DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(y)) * pow<5>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) *
                                                vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(6) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                                DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(8) * pow<2>(_period) *
                                                pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +
                                                (pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *
                                                (Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *
                                                (_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) *
                                                Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * Math::cos(x) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) *
                                                vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *
                                                (Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                                (DataType(2) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) -
                                                DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(vmax_exp) * pow<5>(z) * Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) /
                                                DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<5>(vmax_exp) * pow<5>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) *
                                                pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((vmax_exp *
                                                Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period *
                                                pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_ns_first_term = DataType(2) * _N_s * _mu *((DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +
                                                (pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) *
                                                z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -((_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) *
                                                Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * vmax_exp *
                                                z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(5) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                                Math::sin(x) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -
                                                (pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) *
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(vmax_exp) * z * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *
                                                (Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(vmax_exp) *
                                                pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<2>(vmax_exp) *
                                                z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) *
                                                Math::sin(y) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) *
                                                _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_second_term = - DataType(2) * _N_s * _mu *((pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *
                                                  (Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                                  (pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                  pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp *
                                                  pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                  pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                  pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(6) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) *(Math::sin(y))) /
                                                  (DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * vmax_exp * z *
                                                  Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                  pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                  pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * vmax_exp * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *
                                                  (DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) *
                                                  pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                  pow<2>(Math::sin(y))) -(_period * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) /
                                                  DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(vmax_exp) * pow<3>(z) *
                                                  Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) *
                                                  pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::cos(y) *
                                                  Math::sin(x) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                  pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((vmax_exp * Math::sin(x) *
                                                  (Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) *
                                                  pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                  pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<2>(Math::sin(y)) *
                                                  (DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                  Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                                  (_period * pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) *
                                                  _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                  Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_third_term = - DataType(2) * _N_s * _mu *((DataType(2) * pow<2>(_period) * vmax_exp * z * Math::sin(x) * Math::sin(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                  pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *
                                                  (Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *
                                                  (Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -
                                                  (DataType(4) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                  pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                  pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) *(Math::sin(y))) /(DataType(3) *
                                                  _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * vmax_exp * z * Math::cos(x) *
                                                  Math::sin(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                                  Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                  pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((_period * vmax_exp * Math::cos(x) *
                                                  (Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp *
                                                  pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) +
                                                  DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                  pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<3>(vmax_exp) * pow<3>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * Math::sin(y) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) *
                                                  Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) *
                                                  _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType alp_term = exp * _alpha * _vmax * z * Math::sin(x) * Math::sin(y);
            DataType konv_term = _konv ? (DataType(-0.5) * exp * (_period * _rho * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::cos(y)) * Math::sin(DataType(2) * x))) : DataType(0);
            DataType grap = DataType(2) * exp * exp *  _mu * _period * _period * _vmax * z * Math::sin(x) * Math::sin(y);
            //calculate value[0] through (alpha*v + konv + grap - grtau)/rho, where grap is the gradient of the pressure, and grtau = grtau_d + grtau_np + grtau_ns_first_term
            val[0] = (alp_term + konv_term + grap - (grtau_d + grtau_np_first_term + grtau_np_second_term + grtau_np_third_term + grtau_ns_first_term + grtau_ns_second_term + grtau_ns_third_term))/_rho;
          }
          //second component
          {
            DataType grtau_d = DataType(2) * _mu * _period *_period * vmax_exp * z * Math::cos(x) * Math::cos(y);

            DataType grtau_np_first_term = DataType(2) * _N_p * _mu *((DataType(2) * pow<2>(_period) * vmax_exp * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * vmax_exp * z * Math::cos(y) * Math::sin(x) *
                                              (DataType(2) * _period * pow<4>(vmax_exp) * pow<4>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<4>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)))) /(DataType(15) *
                                              _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<5>(Math::cos(x)) *
                                              Math::cos(y) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                              (pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * Math::cos(x) * Math::cos(y) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)) *((vmax_exp *
                                              Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(3) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(8) * pow<4>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<3>(Math::cos(x)) * Math::cos(y) *
                                              pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                              (DataType(4) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) /
                                              DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +
                                              (pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) *
                                              Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                              (DataType(4) * _period * vmax_exp * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) *
                                              Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *
                                              (pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * Math::cos(x) * Math::cos(y) * pow<3>(Math::sin(x)) *
                                              pow<3>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp *
                                              Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                              Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *
                                              (pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<5>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) *
                                              pow<2>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period *
                                              pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<3>(_period) *
                                              pow<5>(vmax_exp) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                              Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp *
                                              pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) +
                                              DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                              (DataType(4) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                              DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) *
                                              pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /
                                              pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_second_term = - DataType(2) * _N_p * _mu * ((DataType(16) * pow<4>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<5>(Math::cos(x)) * Math::cos(y) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) *
                                              pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                              pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) *
                                              pow<2>(_period) * vmax_exp * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(12) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *
                                              (Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                              pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(12) *
                                              pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<4>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<5>(Math::cos(x)) * Math::cos(y) * pow<4>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) *
                                              pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                              (DataType(12) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                              DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) *
                                              pow<5>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) *
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /
                                              pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<5>(vmax_exp) *
                                              pow<6>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /
                                              pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) *
                                              pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(vmax_exp) * z *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_third_term = DataType(2) * _N_p * _mu *((pow<2>(_period) * vmax_exp * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                              (_period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(DataType(2) * _period * pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) - DataType(2) * _period * pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /
                                              (DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(pow<2>(_period) * vmax_exp * z * Math::cos(x) *
                                              Math::cos(y) *(DataType(3) * _epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) -
                                              DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) *(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<5>(Math::cos(y)) * pow<2>(Math::sin(x))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * vmax_exp * z * Math::cos(x) * Math::cos(y) *
                                              (_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(6) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp *
                                              pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) +(DataType(8) * pow<2>(_period) * pow<5>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(y) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) /
                                              DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +
                                              (pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(vmax_exp) *
                                              pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                              (Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *
                                              (_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                              pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::cos(y))) *(DataType(2) * _period * pow<2>(vmax_exp) *
                                              pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                              (Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * Math::cos(y) *
                                              Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) * Math::sin(y) *
                                              ((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<5>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<4>(Math::cos(y)) * pow<2>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _period *
                                              pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                              (Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                              Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                              Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp *
                                              Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_ns_first_term = DataType(2) * _N_s * _mu *((DataType(2) * pow<2>(_period) * vmax_exp * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) /
                                              DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *
                                              (Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) *
                                              pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(4) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * vmax_exp * z * Math::cos(x) *
                                              Math::sin(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) /
                                              DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<3>(vmax_exp) * pow<3>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * Math::sin(y) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                              Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              Math::cos(y) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                              Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_second_term = DataType(2) * _N_s * _mu *((pow<2>(_period) * vmax_exp * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *
                                              (Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) *
                                              pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * vmax_exp * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -
                                              (DataType(2) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * Math::cos(x) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(x))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(6) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /
                                              (DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * vmax_exp * z * Math::cos(y) *
                                              Math::sin(x) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                              Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * vmax_exp * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *
                                              (DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                              pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                              (_period * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /
                                              (DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(_period * pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) *
                                              Math::sin(x) * pow<2>(Math::sin(y)) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *
                                              ((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) *
                                              _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((vmax_exp *
                                              Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period *
                                              pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *
                                              ((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<2>(Math::sin(y)) *(DataType(2) * pow<3>(_period) *
                                              pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *
                                              (Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_third_term =  - DataType(2) * _N_s * _mu *((DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +
                                              (pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z *
                                              pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) *
                                              pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * vmax_exp *
                                              z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * vmax_exp * z * Math::cos(x) * Math::cos(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(4) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) *
                                              pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<2>(vmax_exp) * z * Math::cos(x) * Math::cos(y) *
                                              Math::sin(x) * Math::sin(y) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(5) * pow<2>(_period) * pow<3>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) /
                                              (DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) *
                                              pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(_period) *
                                              pow<3>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) *
                                              pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                              pow<2>(Math::sin(y))) -(pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(vmax_exp) * z *
                                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) *
                                              _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType alp_term = - exp * _alpha * _vmax * z * Math::cos(x) * Math::cos(y);
            DataType konv_term = _konv ? (DataType(0.5) * exp * exp * (_period * _rho * pow<2>(_vmax) * pow<2>(z) * pow<2>(Math::sin(x)) * Math::sin(DataType(2) * y))) : DataType(0);
            DataType grap = DataType(-2) * exp * _mu * _period * _period * _vmax * z * Math::cos(x) * Math::cos(y);
            //calculate value[0] through (grap - grtau)/rho, where grap is the gradient of the pressure, and grtau = grtau_d + grtau_np + grtau_ns_first_term
            val[1] = (alp_term + konv_term + grap - (grtau_d + grtau_np_first_term + grtau_np_second_term + grtau_np_third_term + grtau_ns_first_term + grtau_ns_second_term + grtau_ns_third_term))/_rho;
          }

          //calculate third component
          {
            DataType grtau_d = DataType(4) * _mu *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2)) - DataType(4) * _mu * _period * vmax_exp * Math::cos(x) *(Math::sin(y));

            DataType grtau_np_first_term = DataType(2) * _N_p * _mu *((DataType(2) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * Math::cos(y) *
                                                pow<3>(Math::sin(y)) - DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) *(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) /
                                                DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) -
                                                DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(2) * pow<5>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<4>(Math::cos(y)) *(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(3) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(8) * pow<5>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) *
                                                pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                                DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<3>(_period) * pow<5>(vmax_exp) *
                                                pow<6>(z) * pow<3>(Math::cos(x)) * pow<4>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(10) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *
                                                ((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x)
                                                *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<3>(Math::cos(y)) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<4>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<5>(Math::cos(x)) * Math::cos(y) * pow<4>(Math::sin(y)) *(DataType(2) * _period *
                                                pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                (Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::cos(y)) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) -
                                                DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -
                                                (DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +
                                                (pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *
                                                Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_np_second_term = - DataType(2) * _N_p * _mu *((DataType(2) * _period * vmax_exp * Math::cos(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_period) * pow<4>(vmax_exp) * pow<8>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period * vmax_exp * Math::cos(x) * Math::sin(y) *(_epsilon +
                                                pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period * vmax_exp * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(16) * pow<5>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<5>(Math::sin(y))) /
                                                (DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) *
                                                pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) +(DataType(2) * pow<5>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(6) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<5>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) *
                                                _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(14) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) *
                                                Math::sin(x) * pow<4>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(14) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) *
                                                ((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(14) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *
                                                (_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period *
                                                vmax_exp * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) *
                                                z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) -(DataType(4) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(DataType(3) * _epsilon + pow<4>(_period) * pow<4>(vmax_exp) * pow<8>(z) * pow<4>(Math::cos(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<3>(_period) * pow<4>(vmax_exp) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp *
                                                pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<3>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) *
                                                pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<3>(_period) * pow<4>(vmax_exp) * pow<7>(z) * pow<3>(Math::cos(x)) * Math::sin(x) *
                                                pow<4>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) *
                                                z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))));

            DataType grtau_np_third_term = - DataType(2) * _N_p * _mu *((DataType(2) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) *
                                                Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) - DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::sin(x) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *((_period * vmax_exp *
                                                Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon +
                                                pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<5>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<5>(Math::cos(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * Math::cos(x) * pow<4>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<5>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon +
                                                pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(3) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) *
                                                pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                                                (DataType(8) * pow<5>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(4) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<4>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(2) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(y)) *((vmax_exp *
                                                Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(10) * pow<3>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) *
                                                Math::sin(x) * pow<3>(Math::sin(y)) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y)) *
                                                (DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                                pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) *
                                                _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(6) * pow<3>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) *
                                                Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<4>(_period) * pow<5>(vmax_exp) * pow<8>(z) * pow<4>(Math::cos(x)) * Math::sin(x) *
                                                pow<5>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) *
                                                pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /
                                                pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<5>(vmax_exp) *
                                                pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) * pow<3>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) -
                                                DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(DataType(4) * pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *
                                                (Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) *
                                                pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x))) *(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType grtau_ns_first_term = DataType(2) * _N_s * _mu *(((_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /
                                                (DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +((_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y))) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *
                                                (Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) *
                                                Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +((_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y))) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +
                                                (pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) *
                                                pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * Math::cos(x) * pow<2>(Math::sin(x)) *
                                                pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(_period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) *
                                                ((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::sin(x) *
                                                pow<3>(Math::sin(y)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period *
                                                pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *
                                                (Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(vmax_exp) *
                                                pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(3) *
                                                _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * Math::sin(y) *
                                                (DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) *
                                                Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) +(pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) /
                                                DataType(2)) *(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) + DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) *
                                                Math::cos(x) * pow<2>(Math::cos(y)) *(Math::sin(x)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))));

            DataType grtau_ns_second_term = - DataType(2) * _N_s * _mu *((DataType(4) * _period * vmax_exp * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(16) * pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::sin(y))) /(DataType(3) *
                                                _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                                (DataType(2) * pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * Math::cos(x) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(6) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp *
                                                pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) +(DataType(6) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                                                DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(4) * _period * vmax_exp *
                                                z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) *
                                                Math::cos(x) *(Math::cos(y))) / DataType(2)) *(DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period *
                                                pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * pow<2>(vmax_exp) * z *
                                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + DataType(2) * pow<2>(vmax_exp) * z * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + DataType(4) * pow<2>(_period) * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType grtau_ns_third_term =  DataType(2) * _N_s * _mu *(((_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) / DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +((_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *((_period * vmax_exp * Math::cos(x) *(Math::sin(y))) /
                                                DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) *
                                                pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) +(((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2)) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y))) *(DataType(2) *
                                                _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) * Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +
                                                (pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * Math::cos(x) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *((vmax_exp *
                                                Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon +
                                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<3>(vmax_exp) * pow<4>(z) * Math::cos(x) *
                                                pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::sin(y))) -(_period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * Math::sin(x) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /
                                                (DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) *
                                                Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *
                                                ((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) * Math::sin(y) *((_period * vmax_exp * Math::cos(y) *(Math::sin(x))) /
                                                DataType(2) +(pow<3>(_period) * vmax_exp * pow<2>(z) * Math::cos(y) *(Math::sin(x))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                                pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) * Math::cos(x) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) *
                                                Math::sin(y) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * Math::cos(y) * pow<2>(Math::sin(x)) *
                                                Math::sin(y) + DataType(2) * pow<3>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) - DataType(2) * _period * pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * Math::cos(y) *(Math::sin(y)))) / pow<2>(DataType(3) * _epsilon + pow<2>(vmax_exp) *
                                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            DataType alp_term = - _alpha * exp * _period * _vmax * pow<2>(z) * Math::cos(x) * Math::sin(y);
            DataType konv_term = _konv ? (pow<2>(_period) * exp * exp * _rho * pow<2>(_vmax) * pow<3>(z) * (pow<2>(Math::cos(x)) + pow<2>(Math::sin(y)))) : DataType(0);
            DataType grap = DataType(-2) * _mu * exp * _period * _vmax * Math::cos(x) * Math::sin(y);
                        //calculate value[0] through (grap - grtau)/rho, where grap is the gradient of the pressure, and grtau = grtau_d + grtau_np + grtau_ns_first_term
            val[2] = (alp_term + konv_term + grap - (grtau_d + grtau_np_first_term + grtau_np_second_term + grtau_np_third_term + grtau_ns_first_term + grtau_ns_second_term + grtau_ns_third_term))/_rho;
          }
          //all done return val
          return val;
        }
      };
    }; // class RHS

    class Robin :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _mu;
      DataType _rho;
      DataType _N_s;
      DataType _N_p;
      DataType _period;
      DataType _epsilon;
      DataType _alpha;

    private:
      DataType _t;

    public:
      explicit Robin(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _mu(func._mu),
      _rho(func._rho),
      _N_s(func._N_s),
      _N_p(func._N_p),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _t(DataType(0.))
      {
      }


      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _mu, _rho, _N_s, _N_p, _period, _epsilon, _alpha, _t;

      public:
        explicit Evaluator(const Robin& function) :
        _vmax(function._vmax),
        _mu(function._mu),
        _rho(function._rho),
        _N_s(function._N_s),
        _N_p(function._N_p),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];
          const DataType exp = Math::exp(_alpha * _t);
          const DataType vmax_exp = _vmax * exp;

          ValueType val;
          //Again build the different components one after the other
          {
            DataType h_d = DataType(4) * _mu * _period * vmax_exp * z * Math::cos(x) * Math::sin(y); //because of stupid choice of p it 2 times the normal h in this component...
            DataType h_np = - DataType(2) * _N_p * _mu *((DataType(2) * _period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(DataType(3) *
                              _epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<4>(Math::sin(x)) * pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) *
                              pow<2>(Math::sin(y)))) -(_period * vmax_exp * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<5>(vmax_exp) * pow<5>(z) * Math::cos(x) * pow<2>(Math::cos(y)) *
                              pow<4>(Math::sin(x)) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +
                              (DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * Math::cos(x) * pow<3>(Math::sin(x)) * pow<4>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) *
                              _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) *
                              Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType h_ns = - DataType(2) * _N_s * _mu *((DataType(2) * _period * pow<2>(vmax_exp) * pow<3>(z) * Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) /
                              DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period *
                              vmax_exp * z * Math::cos(x) * Math::sin(y) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(DataType(2) * _period * pow<3>(vmax_exp) * pow<3>(z) * Math::cos(x) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) *(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            val[0] = (h_d + h_np + h_ns)/_rho;
          }

          {
            DataType h_d = DataType(2) * _mu * _period * vmax_exp * z * Math::cos(y) * Math::sin(x);

            DataType h_np = DataType(2) * _N_p * _mu *((DataType(2) * _period * vmax_exp * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<4>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) *
                                pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period * pow<5>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<3>(Math::sin(x)) *
                                pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(_period *
                                pow<5>(vmax_exp) * pow<5>(z) * pow<4>(Math::cos(x)) * pow<3>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) *
                                pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<3>(_period) * pow<5>(vmax_exp) * pow<7>(z) * pow<4>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<4>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) +
                                pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<3>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<2>(Math::sin(y)) *
                                ((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) *
                                pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * _period * pow<4>(vmax_exp) * pow<5>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * pow<2>(Math::sin(x)) * pow<3>(Math::sin(y)) *((vmax_exp * Math::sin(x) *
                                (Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) +
                                pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType h_ns = DataType(2) * _N_s * _mu *((_period * vmax_exp * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                                pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * vmax_exp * z * Math::cos(y) * Math::sin(x) *(_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)))) /
                                (DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(vmax_exp) * pow<3>(z) *
                                Math::cos(x) * Math::sin(x) * pow<2>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                                pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(_period * pow<2>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(y) *
                                ((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                                pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) -(DataType(2) * _period * pow<3>(vmax_exp) * pow<3>(z) * pow<2>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<2>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) *
                                pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            val[1] = (h_d + h_np + h_ns)/_rho;
          }

          {
            DataType h_d = DataType(2) * _mu * ((vmax_exp * Math::sin(x) * Math::sin(y))/DataType(2) + (_period * _period * vmax_exp * z * z * Math::sin(x) * Math::sin(y))/DataType(2));

            DataType h_np = DataType(2) * _N_p * _mu *((DataType(2) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(x)) *
                              pow<4>(Math::sin(y)))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<4>(_period) *
                              pow<5>(vmax_exp) * pow<8>(z) * pow<4>(Math::cos(x)) * Math::sin(x) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<3>(Math::sin(x)) * pow<5>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                              pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) -(pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<4>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) *
                              pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) * pow<2>(_period) * pow<5>(vmax_exp) * pow<6>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * pow<3>(Math::sin(x)) *
                              pow<3>(Math::sin(y))) /(DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) +(DataType(2) *
                              pow<2>(_period) * pow<4>(vmax_exp) * pow<6>(z) * pow<3>(Math::cos(x)) * Math::cos(y) * Math::sin(x) * pow<3>(Math::sin(y)) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /
                              (DataType(15) * _epsilon + pow<2>(pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))));

            DataType h_ns = DataType(2) * _N_s * _mu *((((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2)) *(_epsilon + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y)))) /(DataType(3) *
                              _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +((_epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) *
                              pow<2>(Math::sin(y))) *((vmax_exp * Math::sin(x) *(Math::sin(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::sin(x) *(Math::sin(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * Math::sin(x) * pow<3>(Math::sin(y))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) *
                              pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(vmax_exp) * pow<2>(z) * Math::cos(x) * Math::cos(y) * Math::sin(x) *
                              Math::sin(y) *((vmax_exp * Math::cos(x) *(Math::cos(y))) / DataType(2) +(pow<2>(_period) * vmax_exp * pow<2>(z) * Math::cos(x) *(Math::cos(y))) / DataType(2))) /(DataType(3) * _epsilon + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) *
                              pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))) +(pow<2>(_period) * pow<3>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) * Math::sin(x) *(Math::sin(y))) /(DataType(3) * _epsilon +
                              pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::cos(y)) + pow<2>(vmax_exp) * pow<2>(z) * pow<2>(Math::sin(x)) * pow<2>(Math::sin(y)) + pow<2>(_period) * pow<2>(vmax_exp) * pow<4>(z) * pow<2>(Math::cos(x)) * pow<2>(Math::sin(y))));

            val[2] =  (h_d + h_np + h_ns)/_rho;
          }

          return val;
        }
      };
    }; // class Robin

    class Orient :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<6> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon, _alpha;

    private:
      DataType _t;
    public:
      explicit Orient(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon, _alpha, _t;
      public:
        explicit Evaluator(const Orient& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];
          const DataType exp = Math::exp(_alpha * _t);

          const DataType v1 = exp * _vmax * Math::sin(x) * Math::sin(y) * z;
          const DataType v2 = DataType(-1) * exp * _vmax * Math::cos(x) * Math::cos(y)*z;
          const DataType v3 = DataType(-1) * exp * _vmax * _period * Math::cos(x) * Math::sin(y) * z * z;
          const DataType temp_norm = DataType(3) * _epsilon + v1 * v1 + v2 * v2 + v3 * v3;




          //the matrix will be represented by a vector by A_ij = f(j*dim + i)
          ValueType val;
//           val[0] = (v1 * v1 + _epsilon) / temp_norm;
//           val[1] = (v1 * v2) / temp_norm;
//           val[2] = (v1 * v3) / temp_norm;
//           val[3] = (v1 * v2) / temp_norm;
//           val[4] = (v2 * v2 + _epsilon) /temp_norm;
//           val[5] = (v2 * v3) / temp_norm;
//           val[6] = (v1 * v3) / temp_norm;
//           val[7] = (v2 * v3) / temp_norm;
//           val[8] = (v3 * v3 + _epsilon) / temp_norm;

          val[0] = (v1 * v1 + _epsilon) / temp_norm;
          val[1] = (v2 * v2 + _epsilon) /temp_norm;
          val[2] = (v3 * v3 + _epsilon) / temp_norm;
          val[3] = (v1 * v2) / temp_norm;
          val[4] = (v1 * v3) / temp_norm;
          val[5] = (v2 * v3) / temp_norm;

          return val;
        }
      };
    }; // class Orient

    class Fourth_Moment :
    public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<15> ImageType;
      static constexpr bool can_value = true;

    protected:
      const DataType _vmax, _period, _epsilon, _alpha;

    private:
      DataType _t;
    public:
      explicit Fourth_Moment(const FullTensorTimeExpo& func) :
      _vmax(func._vmax),
      _period(func._period),
      _epsilon(func._epsilon),
      _alpha(func._alpha),
      _t(DataType(0))
      {
      }

      void set_time(const DataType t)
      {
        _t = t;
      }

      template<typename Traits_>
      class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _period, _epsilon, _alpha, _t;
      public:
        explicit Evaluator(const Fourth_Moment& function) :
        _vmax(function._vmax),
        _period(function._period),
        _epsilon(function._epsilon),
        _alpha(function._alpha),
        _t(function._t)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType x = _period * point[0];
          const DataType y = _period * point[1];
          const DataType z = point[2];
          const DataType exp = Math::exp(_alpha * _t);

          const DataType v_0 = exp * _vmax * Math::sin(x) * Math::sin(y) * z;
          const DataType v_1 = DataType(-1) * exp * _vmax * Math::cos(x) * Math::cos(y)*z;
          const DataType v_2 = DataType(-1) * exp * _vmax * _period * Math::cos(x) * Math::sin(y) * z * z;

          const DataType temp_norm = DataType(15) * _epsilon + pow<2>(v_0 * v_0 + v_1 * v_1 + v_2 * v_2);

          ValueType val;

          //symmetric version
          val[0] = (v_0 * v_0 * v_0 * v_0 + DataType(3) * _epsilon) / temp_norm;
          val[1] = (v_1 * v_1 * v_1 * v_1 + DataType(3) * _epsilon) / temp_norm;
          val[2] = (v_2 * v_2 * v_2 * v_2 + DataType(3) * _epsilon) / temp_norm;
          val[3] = (v_0 * v_0 * v_0 * v_1) / temp_norm;
          val[4] = (v_0 * v_0 * v_0 * v_2) / temp_norm;
          val[5] = (v_1 * v_1 * v_1 * v_0) / temp_norm;
          val[6] = (v_1 * v_1 * v_1 * v_2) / temp_norm;
          val[7] = (v_2 * v_2 * v_2 * v_0) / temp_norm;
          val[8] = (v_2 * v_2 * v_2 * v_1) / temp_norm;
          val[9] = (v_0 * v_0 * v_1 * v_1 + _epsilon) / temp_norm;
          val[10] = (v_0 * v_0 * v_2 * v_2 + _epsilon) / temp_norm;
          val[11] = (v_1 * v_1 * v_2 * v_2 + _epsilon) / temp_norm;
          val[12] = (v_0 * v_0 * v_1 * v_2) / temp_norm;
          val[13] = (v_1 * v_1 * v_0 * v_2) / temp_norm;
          val[14] = (v_2 * v_2 * v_0 * v_1) / temp_norm;

          return val;

        }
      };
    }; // class Fourth_Moment


  }; //class FullTensorTimeExpo 3D spezialization



}
#endif