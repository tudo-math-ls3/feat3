#pragma once
#ifndef KERNEL_ASSEMBLY_PARSED_FUNCTION_HPP
#define KERNEL_ASSEMBLY_PARSED_FUNCTION_HPP 1

#include <kernel/assembly/analytic_function.hpp>

// The contents of this file require the 'fparser' third-party library.
#if defined(FEAST_HAVE_FPARSER) || defined(DOXYGEN)

#include <vector>
#include <fparser.hh>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Parsed function implementation
     *
     * This class provides an implementatuon of the AnalyticFunction interface which
     * can parse and evaluate a formula in the up to three variables 'x', 'y' and 'z'
     * given as a string at runtime.
     *
     * \attention
     * This class is only available if FEAST is configured and build with support
     * for the <c>fparser</c> third-party library, which is used for the actual
     * parsing and evaluation.
     *
     * \note
     * This class does not only offer the possibility to evaluate the function values of
     * the parsed formula, but also offers the numerical evaluation of gradients
     * (first-order derivatives) and hessians (second-order derivatives) via
     * Richardson extrapolation applied onto second-order central difference quotients!
     * The initial 'h' for the difference quotient as well as the maximum number of
     * extrapolation steps can be adjusted by using the #set_grad_extrapol() and
     * #set_hess_extrapol() functions, respectively.
     *
     * For a full documentation and a list of supported expressions, see the documentation
     * of the undrlying parser at http://warp.povusers.org/FunctionParser/fparser.html
     *
     * The (up to) three variables of the function to be parsed are named \c x, \c y and \c z.
     *
     * There are two possibilities how to create an instance of this class:
     * If you have a 'simple' function without any symbolic constants, you can use the
     * constructor which takes the function string as well as the dimension as parameters.
     *
     * If your function requires additional (problem specific) constants (e.g. Reynolds number),
     * you first need to create an instance of this class by using the standard constructor,
     * then add all your required constants by using the #add_constant function and finally
     * parse your function string by calling the #parse function.
     *
     * This class already offers the following pre-defined constants:
     * - <c>pi</c> = 3.14159...
     * - <c>eps</c> = ~1E-16
     *
     * \tparam dim_
     * The dimension of the function, i.e. the number of variables. Must be 1 <= dim_ <= 3.
     *
     * \author Peter Zajac
     */
    template<int dim_>
    class ParsedFunction :
      public AnalyticFunction
    {
      static_assert((dim_ >= 1) && (dim_ <= 3), "unsupported function dimension");

    private:
      /// the actual function parser object
      ::FunctionParser _parser;
      /// maximum number of gradient extrapolation steps
      int _max_grad_steps;
      /// maximum number of hessian extrapolation steps
      int _max_hess_steps;
      /// initial h for gradient extrapolation
      double _init_grad_h;
      /// initial h for hessian extrapolation
      double _init_hess_h;

    public:
      /**
       * \brief Standard constructor
       *
       * This constructor adds the following constants to the underlying parser:
       * - <c>pi</c> = 3.14159...
       * - <c>eps</c> = ~1E-16
       *
       * Moreover, this constructor configures the Richardson extrapolation steps
       * for both the gradient and hessian evaluation to:
       * - initial h = 1E-2
       * - maximum steps = 10
       * which should be sufficient in most (non-exotic) cases.
       */
      explicit ParsedFunction() :
        _parser(),
        _max_grad_steps(10),
        _max_hess_steps(10),
        _init_grad_h(1E-2),
        _init_hess_h(1E-2)
      {
        // add constants to our parser
        _parser.AddConstant("pi", Math::pi<double>());
        _parser.AddConstant("eps", Math::eps<double>());
      }

      /**
       * \brief Constructor
       *
       * This constructor creates a parsed function from a String.
       *
       * \param[in] function
       * The expession that defines the function in the variables \c x, \c y, and \c z.
       */
      explicit ParsedFunction(const String& function) :
        ParsedFunction()
      {
        parse(function);
      }

      /**
       * \brief Adds a constant to the parser.
       *
       * \param[in] name
       * The name of the constant.
       *
       * \param[in] value
       * The value of the constant.
       */
      void add_constant(const String& name, double value)
      {
        _parser.AddConstant(name.c_str(), value);
      }

      /**
       * \brief Configures the Gradient evaluation.
       *
       * This function configures the Richardson extrapolation scheme for the
       * evaluation of function gradients (first-order derivatives).
       *
       * \param[in] initial_h
       * The initial h for the difference quotent. Must be > 0; default = 1E-2.
       *
       * \param[in] max_steps
       * The maximum number of Richardson extrapolation steps to be performed.
       * Must be > 0; default = 10.
       */
      void set_grad_extrapol(double initial_h, int max_steps)
      {
        if(initial_h < Math::eps<double>())
          throw InternalError("Initial h is too small or non-positive!");
        if(max_steps <= 0)
          throw InternalError("Invalid maximum extrapolation steps!");
        _init_grad_h = initial_h;
        _max_grad_steps = max_steps;
      }

      /**
       * \brief Configures the Hessian evaluation.
       *
       * This function configures the Richardson extrapolation scheme for the
       * evaluation of function hessians (second-order derivatives).
       *
       * \param[in] initial_h
       * The initial h for the difference quotent. Must be > 0; default = 1E-2.
       *
       * \param[in] max_steps
       * The maximum number of Richardson extrapolation steps to be performed.
       * Must be > 0; default = 10.
       */
      void set_hess_extrapol(double initial_h, int max_steps)
      {
        if(initial_h < Math::eps<double>())
          throw InternalError("Initial h is too small or non-positive!");
        if(max_steps <= 0)
          throw InternalError("Invalid maximum extrapolation steps!");
         _init_hess_h = initial_h;
        _max_hess_steps = max_steps;
      }

      /**
       * \brief Parses a function.
       *
       * \param[in] function
       * The expession that defines the function in the variables \c x, \c y, and \c z.
       */
      void parse(const String& function)
      {
          // add variables to our parser
        String vars("x");
        if(dim_ > 1) vars += ",y";
        if(dim_ > 2) vars += ",z";

        // try to parse the function
        if(_parser.Parse(function.c_str(), vars.c_str()) >= 0)
          throw InternalError("Failed to parse function");

        // optimise the parsed function
        _parser.Optimize();
      }

      /// we can compute function values
      static constexpr bool can_value = true;

      /// we can compute function gradients
      static constexpr bool can_grad = true;

      /// we can compute function gradients
      static constexpr bool can_hess = true;

      template<typename Config_>
      struct ConfigTraits
      {
        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          /// we always need image point coordinates
          static constexpr bool need_img_point = true;
        };
      };

      template<typename EvalTraits_>
      class Evaluator :
        public AnalyticFunction::Evaluator<EvalTraits_>
      {
      public:
        /// trafo evaluator data
        typedef typename EvalTraits_::TrafoEvaluator TrafoEvaluator;
        /// trafo data type
        typedef typename EvalTraits_::TrafoData TrafoData;
        /// coefficient data type
        typedef typename EvalTraits_::DataType DataType;
        /// value type
        typedef typename EvalTraits_::ValueType ValueType;
        /// gradient type
        typedef typename EvalTraits_::GradientType GradientType;
        /// hessian type
        typedef typename EvalTraits_::HessianType HessianType;

        /// verify our dimension
        static_assert(TrafoEvaluator::image_dim == dim_, "function evaluator dimension mismatch!");

      private:
        // Note: We need to create a copy of the parser rather than just a reference to it,
        // as the 'Eval' function of the parser is not thread-safe, whereas our Evaluator's
        // are always meant to be so!
        /// our own private copy of the function parser
        ::FunctionParser _parser;
        /// our gradient extrapolation table
        std::vector<GradientType> _grad;
        /// our hessian extrapolation table
        std::vector<HessianType> _hess;
        /// initial h for gradient extrapolation
        const double _init_grad_h;
        /// initial h for hessian extrapolation
        const double _init_hess_h;

      public:
        explicit Evaluator(const ParsedFunction& function) :
          _parser(function._parser),
          _grad(std::size_t(function._max_grad_steps)),
          _hess(std::size_t(function._max_hess_steps)),
          _init_grad_h(function._init_grad_h),
          _init_hess_h(function._init_hess_h)
        {
        }

        ValueType value(const TrafoData& tau)
        {
          // copy our image point into a Tiny::Vector<double>
          Tiny::Vector<double, dim_> v(tau.img_point);

          // and evaluate the parsed function
          return ValueType(_parser.Eval(v.v));
        }

        GradientType gradient(const TrafoData& tau)
        {
          // first, copy our image point into a Tiny::Vector<double>
          Tiny::Vector<double, dim_> v(tau.img_point);

          // next, choose the inital h
          double h(_init_grad_h);

          // Note: the '_grad' vector was already allocated to
          //       the correct size by our constructor

          // evaluate gradient
          _eval_grad(_grad[0], v.v, h);

          // Now comes the Richardson extrapolation loop
          const std::size_t n(_grad.size()-1);
          DataType def(1E+10);
          for(std::size_t i(0); i < n; ++i)
          {
            // evaluate next gradient
            _eval_grad(_grad[i+1], v.v, h *= 0.5);

            // initialise scaling fator
            DataType q = DataType(1);

            // perform extrapolation steps except for the last one
            for(std::size_t k(i); k > std::size_t(0); --k)
            {
              q *= DataType(4);
              (_grad[k] -= q*_grad[k+1]) *= DataType(1) / (DataType(1) - q);
            }

            // compute and check our defect
            DataType d = (_grad[1] - _grad[0]).norm_euclid_sqr();
            if(def <= d)
            {
              // The defect has increased, so we return the previous extrapolation result
              return _grad.front();
            }

            // remember current defect
            def = d;

            // perform last extrapolation step
            q *= DataType(4);
            (_grad[0] -= q*_grad[1]) *= DataType(1) / (DataType(1) - q);
          }

          // return our extrapolated gradient
          return _grad.front();
        }

        HessianType hessian(const TrafoData& tau)
        {
          // first, copy our image point into a Tiny::Vector<double>
          Tiny::Vector<double, dim_> v(tau.img_point);

          // next, choose the inital h
          double h(_init_hess_h);

          // Note: the '_hess' vector was already allocated to
          //       the correct size by our constructor

          // evaluate hessian
          _eval_hess(_hess[0], v.v, h);

          // Now comes the Richardson extrapolation loop
          const std::size_t n(_hess.size()-1);
          DataType def(1E+10);
          for(std::size_t i(0); i < n; ++i)
          {
            // evaluate next hessian
            _eval_hess(_hess[i+1], v.v, h *= 0.5);

            // initialise scaling fator
            DataType q = DataType(1);

            // perform extrapolation steps except for the last one
            for(std::size_t k(i); k > std::size_t(0); --k)
            {
              q *= DataType(4);
              (_hess[k] -= q*_hess[k+1]) *= DataType(1) / (DataType(1) - q);
            }

            // compute and check our defect
            DataType d = (_hess[1] - _hess[0]).hessian_sqr_norm();
            if(def <= d)
            {
              // The defect has increased, so we return the previous extrapolation result
              return _hess.front();
            }

            // remember current defect
            def = d;

            // perform last extrapolation step
            q *= DataType(4);
            (_hess[0] -= q*_hess[1]) *= DataType(1) / (DataType(1) - q);
          }

          // return our extrapolated hessian
          return _hess.front();
        }

      private:
        /// evaluates the first-order difference quotients
        void _eval_grad(GradientType& x, double* v, const double h)
        {
          // loop over all dimensions
          for(int i(0); i < dim_; ++i)
          {
            // backup coordinate i
            const double vi(v[i]);

            // evaluate f(v+h)
            v[i] = vi + h;
            DataType vr = DataType(_parser.Eval(v));

            // evaluate f(v-h)
            v[i] = vi - h;
            DataType vl = DataType(_parser.Eval(v));

            // compute difference quotient
            x[i] = (vr - vl) / DataType(2.0*h);

            // restore coord
            v[i] = vi;
          }
        }

        /// evaluates the second-order difference quotients
        void _eval_hess(HessianType& x, double* v, const double h)
        {
          // evaluate at center
          const DataType vc = DataType(_parser.Eval(v));

          // loop over all dimensions
          for(int i(0); i < dim_; ++i)
          {
            // backup coord
            const double vi(v[i]);

            // eval f(x+h)
            v[i] = vi + h;
            DataType vr = DataType(_parser.Eval(v));

            // eval f(x-h)
            v[i] = vi - h;
            DataType vl = DataType(_parser.Eval(v));

            // compute difference quotient
            x[i][i] = DataType(vr + vl - 2.0*vc) / DataType(h*h);

            // now the mixed derivatives
            for(int j(i+1); j < dim_; ++j)
            {
              // backup coord
              const double vj(v[j]);

              // we need four points here:
              // north-east
              v[i] = vi + h;
              v[j] = vj + h;
              DataType vne = DataType(_parser.Eval(v));

              // north-west
              v[i] = vi - h;
              v[j] = vj + h;
              DataType vnw = DataType(_parser.Eval(v));

              // south-east
              v[i] = vi + h;
              v[j] = vj - h;
              DataType vse = DataType(_parser.Eval(v));

              // north-west
              v[i] = vi - h;
              v[j] = vj - h;
              DataType vsw = DataType(_parser.Eval(v));

              // combinte into difference quotient
              x[i][j] = x[j][i] = ((vne + vsw) - (vnw + vse)) / DataType(4.0*h*h);

              // restore coord
              v[j] = vj;
            }

            // restore coord
            v[i] = vi;
          }
        }
      }; // class ParsedFunction::Evaluator<...>
    }; // class ParsedFunction
  } // namespace Assembly
} // namespace FEAST

#endif // FEAST_HAVE_FPARSER
#endif // KERNEL_ASSEMBLY_PARSED_FUNCTION_HPP
