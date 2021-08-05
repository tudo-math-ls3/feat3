#pragma once
#ifndef AREA51_CCND_FUNCTION_FIBER
#define AREA51_CCND_FUNCTION_FIBER 1

#include <area51/ccnd_fiber/ccnd_fiber_common.hpp>
#include <kernel/analytic/function.hpp>


namespace CCND_FIBER
{
  using namespace FEAT;
  // 3D steady inflow function used in bench1 and bench2

  //Generell template for ContractInflowFunction, should be specialised for the different dimensions
  template<int dim_>
  class ContractInflowFunction;

  template<>
  class ContractInflowFunction<3> :
  public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 3;
    typedef Analytic::Image::Vector<3> ImageType;
    static constexpr bool can_value = true;

  public:
    template<typename Traits_>
    class Evaluator :
    public Analytic::Function::Evaluator<Traits_>
    {
    protected:
      typedef typename Traits_::DataType DataType;
      typedef typename Traits_::PointType PointType;
      typedef typename Traits_::ValueType ValueType;


    public:
      explicit Evaluator(const ContractInflowFunction&)
      {
      }

      ValueType value(const PointType& point)
      {
        const DataType x = point[0];
        const DataType y = point[1];
//         const DataType z = point[2];

        ValueType val;
        val[0] = val[1] = DataType(0);
        //val[2] = -(1.0 - ((x*x + y*y) / (0.0225*0.0225)));
        val[2] = (1.0 - Math::sqr(x/0.0225) - Math::sqr(y/0.0225)) / 45.0;
        return val;
      }
    };
  }; // class ContractInflowFunction for dim = 3

  template<>
  class ContractInflowFunction<2> :
  public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 2;
    typedef Analytic::Image::Vector<2> ImageType;
    static constexpr bool can_value = true;

  public:
    template<typename Traits_>
    class Evaluator :
    public Analytic::Function::Evaluator<Traits_>
    {
    protected:
      typedef typename Traits_::DataType DataType;
      typedef typename Traits_::PointType PointType;
      typedef typename Traits_::ValueType ValueType;


    public:
      explicit Evaluator(const ContractInflowFunction&)
      {
      }

      ValueType value(const PointType& point)
      {
        //important, our contraction goes from -4 to 4... so we should scale the input point to [0,1]
        const DataType y = DataType(0.125) *(point[1] + DataType(4));
        constexpr DataType _vmax = DataType(0.3);
        constexpr DataType _d = DataType(1.);
        constexpr DataType _den = _d * _d;


        ValueType val;
        val[0] = (_vmax * DataType(4) * y * (_d - y)) / _den;
        val[1] = DataType(0);
        return val;
      }
    };
  }; // class ContractInflowFunction for dim = 2
}

#endif