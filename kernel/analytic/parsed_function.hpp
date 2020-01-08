// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ANALYTIC_PARSED_FUNCTION_HPP
#define KERNEL_ANALYTIC_PARSED_FUNCTION_HPP 1

#include <kernel/analytic/function.hpp>
#include <kernel/util/exception.hpp>

// The contents of this file require the 'fparser' third-party library.
#if defined(FEAT_HAVE_FPARSER) || defined(DOXYGEN)

#include <vector>
#include <fparser.hh>

namespace FEAT
{
  namespace Analytic
  {
    /**
     * \brief Parsed Function parse error
     *
     * An instance of this exception may is thrown if the parsing process
     * of a ParsedFunction formula fails for some reason.
     *
     * The error message contains a description of the error cause and
     * can be queried by calling the inherited what() method.
     *
     * \author Peter Zajac
     */
    class ParsedFunctionParseError : public ParseError
    {
    public:
      /// constructor
      explicit ParsedFunctionParseError(const String& msg) : ParseError(msg) {}
    };

    /**
     * \brief Parsed Function evaluation error
     *
     * An instance of this exception may is thrown if the evaluation of
     * a ParsedFunction evaluator fails for some reason, e.g. division by zero.
     *
     * \author Peter Zajac
     */
    class ParsedFunctionEvalError : public Exception
    {
    public:
      /// constructor
      explicit ParsedFunctionEvalError(const String& msg) : Exception(msg) {}
    };

    /**
     * \brief Parsed function implementation
     *
     * This class provides an implementation of the Analytic::Function interface which
     * can parse and evaluate a formula in the up to three variables 'x', 'y' and 'z'
     * given as a string at runtime.
     *
     * \attention
     * This class is only available if FEAT is configured and build with support
     * for the <c>fparser</c> third-party library, which is used for the actual
     * parsing and evaluation.
     *
     * For a full documentation and a list of supported expressions, see the documentation
     * of the underlying parser at http://warp.povusers.org/FunctionParser/fparser.html
     *
     * There are two possibilities how to create an instance of this class:
     * If you have a 'simple' function without any symbolic constants, you can use the
     * constructor which takes the function string as a single parameter.
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
     * <b>Note to Implementers:</b>\n
     * The <c>FunctionParser</c> class, which is defined by the 'fparser' library, is actually
     * just a typedef for the template class instance <c>FunctionParserBase<double></c>. This
     * may lead the reader to the (erroneous) conclusion, that one might easily templatise this
     * ParsedFunction class in the datatype, thus adding support for other interesting floating
     * point types like the <c>__float128</c> type of GCC's libquadmath. Unfortunately, most of
     * the auxiliary function templates implemented in the depths of the 'fparser' library are
     * only specialised for the common build-in types without offering any generic implementation.
     * Therefore, trying to use the <c>FunctionParserBase</c> class template with a somewhat
     * interesting type like <c>__float128</c> is unfortunately doomed to end in linker errors :(
     *
     * \author Peter Zajac
     */
    template<int dim_>
    class ParsedFunction :
      public Analytic::Function
    {
    public:
      /// validate our dimension
      static_assert((dim_ >= 1) && (dim_ <= 3), "unsupported function dimension");

      /// specify our domain dimension
      static constexpr int domain_dim = dim_;

      /// This function is always scalar
      typedef Analytic::Image::Scalar ImageType;

      /// we can compute function values
      static constexpr bool can_value = true;

      /// we cannot compute gradients
      static constexpr bool can_grad = false;

      /// we cannot compute hessians
      static constexpr bool can_hess = false;

    private:
      /// the actual function parser object
      ::FunctionParser _parser;

    public:
      /**
       * \brief Standard constructor
       *
       * This constructor adds the following constants to the underlying parser:
       * - <c>pi</c> = 3.14159...
       * - <c>eps</c> = ~1E-16
       */
      explicit ParsedFunction() :
        _parser()
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
       * The expression that defines the function in the variables \c x, \c y, and \c z.
       *
       * \throws ParsedFunctionParseError
       * An instance of the ParsedFunctionParseError exception is thrown if the
       * fparser library fails to parse the formula. The message of the exception
       * contains more information on the cause of the error and should be presented
       * the user in an appropriate way.
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
       * \brief Parses a function.
       *
       * \param[in] function
       * The expression that defines the function in the variables \c x, \c y, and \c z.
       *
       * \throws ParsedFunctionParseError
       * An instance of the ParsedFunctionParseError exception is thrown if the
       * fparser library fails to parse the formula. The message of the exception
       * contains more information on the cause of the error and should be presented
       * the user in an appropriate way.
       */
      void parse(const String& function)
      {
          // add variables to our parser
        String vars("x");
        if(dim_ > 1) vars += ",y";
        if(dim_ > 2) vars += ",z";

        // try to parse the function
        const int ret = _parser.Parse(function.c_str(), vars.c_str());
        if(ret >= 0)
        {
          String msg(_parser.ErrorMsg());
          msg.append("\n>>> '");
          msg.append(function);
          msg.append("'");
          if(ret < int(function.size()))
          {
            // ret contains the index of the first invalid character in the input string
            // append an additional line to mark the faulty character
            msg.append("\n>>>");
            msg.append(String(std::size_t(ret+2), '-'));
            msg.append("^");
          }
          throw ParsedFunctionParseError(msg);
        }

        // optimise the parsed function
        _parser.Optimize();
      }

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      public:
        /// coefficient data type
        typedef typename Traits_::DataType DataType;
        /// evaluation point type
        typedef typename Traits_::PointType PointType;
        /// value type
        typedef typename Traits_::ValueType ValueType;

      private:
        // Note: We need to create a copy of the parser rather than just a reference to it,
        // as the 'Eval' function of the parser is not thread-safe, whereas our Evaluators
        // are always meant to be so!
        /// our own private copy of the function parser
        ::FunctionParser _parser;

      public:
        explicit Evaluator(const ParsedFunction& function) :
          _parser(function._parser)
        {
          // force deep copy to avoid race conditions in case of multi-threading
          _parser.ForceDeepCopy();
        }

        ValueType value(const PointType& point)
        {
          // copy our image point into a Tiny::Vector<double>
          Tiny::Vector<double, dim_> v(point);

          // evaluate the parser
          const double val = _parser.Eval(v.v);

          // check for errors
          switch(_parser.EvalError())
          {
          case 0: // no error
            break;

          case 1: // division by zero
            throw ParsedFunctionEvalError("Error in ParsedFunction evaluation: division by zero");

          case 2: // sqrt error
          case 3: // log error
          case 4: // trigonometric error
            throw ParsedFunctionEvalError("Error in ParsedFunction evaluation: illegal input value");

          case 5: // recursion error
            throw ParsedFunctionEvalError("Error in ParsedFunction evaluation: maximum recursion depth reached");

          default: // ???
            throw ParsedFunctionEvalError("Error in ParsedFunction evaluation: unknown error");
          }

          // return value
          return ValueType(val);
        }
      }; // class ParsedFunction::Evaluator<...>
    }; // class ParsedFunction
  } // namespace Analytic
} // namespace FEAT

#endif // FEAT_HAVE_FPARSER
#endif // KERNEL_ANALYTIC_PARSED_FUNCTION_HPP
