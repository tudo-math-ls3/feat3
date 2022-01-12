// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ANALYTIC_PARSED_FUNCTION_HPP
#define KERNEL_ANALYTIC_PARSED_FUNCTION_HPP 1

#include <kernel/analytic/function.hpp>
#include <kernel/util/exception.hpp>

// The contents of this file require the 'fparser' third-party library.
#if defined(FEAT_HAVE_FPARSER) || defined(DOXYGEN)

#include <fparser.hh>

#include <array>
#include <vector>
#include <map>

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
     * \brief Parsed scalar function implementation
     *
     * This class provides an implementation of the Analytic::Function interface which
     * can parse and evaluate a formula in the up to three variables 'x', 'y' and 'z'
     * as well as user-defined extra variables given as a string at runtime.
     *
     * \attention
     * This class is only available if FEAT is configured and build with support
     * for the <c>fparser</c> third-party library, which is used for the actual
     * parsing and evaluation.
     *
     * For a full documentation and a list of supported expressions, see the documentation
     * of the underlying parser at the page \ref parsed_function.
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
     * may lead the reader to the (erroneous) conclusion, that one might easily templatize this
     * ParsedFunction class in the datatype, thus adding support for other interesting floating
     * point types like the <c>__float128</c> type of GCC's libquadmath. Unfortunately, most of
     * the auxiliary function templates implemented in the depths of the 'fparser' library are
     * only specialized for the common build-in types without offering any generic implementation.
     * Therefore, trying to use the <c>FunctionParserBase</c> class template with a somewhat
     * interesting type like <c>__float128</c> is unfortunately doomed to end in linker errors :(
     *
     * \author Peter Zajac
     */
    template<int dim_>
    class ParsedScalarFunction :
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
      /// our extra variable names
      std::map<String, double> _extra_vars;

    public:
      /**
       * \brief Standard constructor
       *
       * This constructor adds the following constants to the underlying parser:
       * - <c>pi</c> = 3.14159...
       * - <c>eps</c> = ~1E-16
       */
      explicit ParsedScalarFunction() :
        _parser(),
        _extra_vars()
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
      explicit ParsedScalarFunction(const String& function) :
        ParsedScalarFunction()
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
       * \brief Adds an extra variable to the parser.
       *
       * The value of the extra variable can be changed later on by using the #set_variable function.
       *
       * \param[in] name
       * The name of the variable. Must not be 'x', 'y' or 'z'.
       *
       * \param[in] value
       * The value for the extra variable.
       */
      void add_variable(const String& name, double value = 0.0)
      {
        XASSERTM(_extra_vars.find(name) == _extra_vars.end(), "extra variable '" + name + "' already added");
        XASSERTM(name.compare_no_case("x") != 0, "variable name 'x' is reserved and cannot be used as extra variable");
        XASSERTM(name.compare_no_case("y") != 0, "variable name 'y' is reserved and cannot be used as extra variable");
        XASSERTM(name.compare_no_case("z") != 0, "variable name 'z' is reserved and cannot be used as extra variable");
        _extra_vars.emplace(name, value);
      }

      /**
       * \brief Sets the value of an extra variable.
       *
       * \param[in] name
       * The name of the variable. Must have been added by using the #add_variable function.
       *
       * \param[in] value
       * The value for the extra variable.
       */
      void set_variable(const String& name, double value)
      {
        auto it = _extra_vars.find(name);
        XASSERTM(it != _extra_vars.end(), "no variable named '" + name + "' found");
        it->second = value;
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

        // add extra variables
        for(auto& v : _extra_vars)
          (vars += ",") += v.first;

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

        // optimize the parsed function
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
        /// a vector for the variable values
        std::vector<double> _vars;

      public:
        explicit Evaluator(const ParsedScalarFunction& function) :
          _parser(function._parser),
          _vars(std::size_t(dim_), 0.0)
        {
          // force deep copy to avoid race conditions in case of multi-threading
          _parser.ForceDeepCopy();
          // set extra vars
          for(const auto& v : function._extra_vars)
            _vars.push_back(v.second);
        }

        ValueType value(const PointType& point)
        {
          // copy coordinates
          for(int i(0); i < dim_; ++i)
            _vars[std::size_t(i)] = point[i];

          // evaluate the parser
          const double val = _parser.Eval(_vars.data());

          // check for errors
          switch(_parser.EvalError())
          {
          case 0: // no error
            break;

          case 1: // division by zero
            throw ParsedFunctionEvalError("Error in ParsedScalarFunction evaluation: division by zero");

          case 2: // sqrt error
          case 3: // log error
          case 4: // trigonometric error
            throw ParsedFunctionEvalError("Error in ParsedScalarFunction evaluation: illegal input value");

          case 5: // recursion error
            throw ParsedFunctionEvalError("Error in ParsedScalarFunction evaluation: maximum recursion depth reached");

          default: // ???
            throw ParsedFunctionEvalError("Error in ParsedScalarFunction evaluation: unknown error");
          }

          // return value
          return ValueType(val);
        }
      }; // class ParsedScalarFunction::Evaluator<...>
    }; // class ParsedScalarFunction

    /// for the sake of downwards compatibility
    template<int dim_>
    using ParsedFunction = ParsedScalarFunction<dim_>;


    /**
     * \brief Parsed vector function implementation
     *
     * This class provides an implementation of the Analytic::Function interface for vector
     * functions which can parse and evaluate a formula tuple in the up to three variables
     * 'x', 'y' and 'z' as well as user-defined extra variables given as a string at runtime.
     *
     * \tparam dom_dim_
     * The domain dimension of the vector function. Must be 1 <= dom_dim_ <= 3.
     *
     * \tparam img_dim_
     * The image dimension of the vector function. Must be >= 1.
     *
     * See the documentation of ParsedScalarFunction for more details.
     *
     * \author Peter Zajac
     */
    template<int dom_dim_, int img_dim_ = dom_dim_>
    class ParsedVectorFunction :
      public Analytic::Function
    {
    public:
      /// validate our dimensions
      static_assert((dom_dim_ >= 1) && (dom_dim_ <= 3), "unsupported domain dimension");
      static_assert(img_dim_ >= 1, "unsupported image dimension");

      /// specify our domain dimension
      static constexpr int domain_dim = dom_dim_;

      /// This function is always a vector field
      typedef Analytic::Image::Vector<img_dim_> ImageType;

      /// we can compute function values
      static constexpr bool can_value = true;

      /// we cannot compute gradients
      static constexpr bool can_grad = false;

      /// we cannot compute hessians
      static constexpr bool can_hess = false;

    private:
      /// our extra variable names
      std::map<String, double> _extra_vars;

      /// the actual function parser objects
      std::array<::FunctionParser, img_dim_> _parsers;

    public:
      /**
       * \brief Standard constructor
       *
       * This constructor adds the following constants to the underlying parser:
       * - <c>pi</c> = 3.14159...
       * - <c>eps</c> = ~1E-16
       */
      explicit ParsedVectorFunction()
      {
        // add constants to our parsers
        add_constant("pi", Math::pi<double>());
        add_constant("eps", Math::eps<double>());
      }

      /**
       * \brief Constructor
       *
       * This constructor creates a parsed function from a String.
       *
       * \param[in] function
       * The expression that defines the function in the variables \c x, \c y, and \c z.
       * The individual parts of the vector function components are separated by a single quotation mark <c>'</c>.
       *
       * \throws ParsedFunctionParseError
       * An instance of the ParsedFunctionParseError exception is thrown if the
       * fparser library fails to parse the formula. The message of the exception
       * contains more information on the cause of the error and should be presented
       * the user in an appropriate way.
       */
      explicit ParsedVectorFunction(const String& function) :
        ParsedVectorFunction()
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
        for(auto& p : _parsers)
          p.AddConstant(name.c_str(), value);
      }

      /**
       * \brief Adds an extra variable to the parser.
       *
       * The value of the extra variable can be changed later on by using the #set_variable function.
       *
       * \param[in] name
       * The name of the variable. Must not be 'x', 'y' or 'z'.
       *
       * \param[in] value
       * The value for the extra variable.
       */
      void add_variable(const String& name, double value = 0.0)
      {
        XASSERTM(_extra_vars.find(name) == _extra_vars.end(), "extra variable '" + name + "' already added");
        XASSERTM(name.compare_no_case("x") != 0, "variable name 'x' is reserved and cannot be used as extra variable");
        XASSERTM(name.compare_no_case("y") != 0, "variable name 'y' is reserved and cannot be used as extra variable");
        XASSERTM(name.compare_no_case("z") != 0, "variable name 'z' is reserved and cannot be used as extra variable");
        _extra_vars.emplace(name, value);
      }

      /**
       * \brief Sets the value of an extra variable.
       *
       * \param[in] name
       * The name of the variable. Must have been added by using the #add_variable function.
       *
       * \param[in] value
       * The value for the extra variable.
       */
      void set_variable(const String& name, double value)
      {
        auto it = _extra_vars.find(name);
        XASSERTM(it != _extra_vars.end(), "no variable named '" + name + "' found");
        it->second = value;
      }

      /**
       * \brief Parses a function.
       *
       * \param[in] function
       * The expression that defines the function in the variables \c x, \c y, and \c z.
       * The individual parts of the vector function components are separated by a single quotation mark <c>'</c>.
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
        if(dom_dim_ > 1) vars += ",y";
        if(dom_dim_ > 2) vars += ",z";

        // add extra variables
        for(auto& v : _extra_vars)
          (vars += ",") += v.first;

        String sfunc(function);
        if(sfunc.starts_with('['))
          sfunc.pop_front();
        if(sfunc.ends_with(']'))
          sfunc.pop_back();
        sfunc.trim_me();

        // split function string
        std::deque<String> sfuncs = sfunc.split_by_charset("'");
        if(sfuncs.size() != _parsers.size())
        {
          String msg("Invalid number of formulae: expected ");
          msg.append(stringify(_parsers.size()));
          msg.append(" but got ");
          msg.append(stringify(sfuncs.size()));
          throw ParsedFunctionParseError(msg);
        }

        // parse all functions
        for(std::size_t i(0); i < _parsers.size(); ++i)
        {
          // try to parse the function
          const int ret = _parsers.at(i).Parse(sfuncs.at(i).c_str(), vars.c_str());
          if(ret >= 0)
          {
            String msg(_parsers.at(i).ErrorMsg());
            msg.append("\n>>> '");
            msg.append(sfuncs.at(i));
            msg.append("'");
            if(ret < int(sfuncs.at(i).size()))
            {
              // ret contains the index of the first invalid character in the input string
              // append an additional line to mark the faulty character
              msg.append("\n>>>");
              msg.append(String(std::size_t(ret+2), '-'));
              msg.append("^");
            }
            throw ParsedFunctionParseError(msg);
          }

          // optimize the parsed function
          _parsers.at(i).Optimize();
        }
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
        /// our own private copy of the function parsers
        std::array< ::FunctionParser, img_dim_> _parsers;
        /// a vector for the variable values
        std::vector<double> _vars;

      public:
        explicit Evaluator(const ParsedVectorFunction& function) :
          _parsers(function._parsers),
          _vars(std::size_t(dom_dim_), 0.0)
        {
          // force deep copy to avoid race conditions in case of multi-threading
          for(auto& p : _parsers)
            p.ForceDeepCopy();

          // set extra vars
          for(const auto& v : function._extra_vars)
            _vars.push_back(v.second);
        }

        ValueType value(const PointType& point)
        {
          // copy coordinates
          for(int i(0); i < dom_dim_; ++i)
            _vars[std::size_t(i)] = point[i];

          // evaluate the parser
          ValueType val;
        for(std::size_t i(0); i < _parsers.size(); ++i)
          {
            // evaluate the parser
            val[int(i)] = DataType(_parsers[i].Eval(_vars.data()));

            // check for errors
            switch(_parsers[i] .EvalError())
            {
            case 0: // no error
              break;

            case 1: // division by zero
              throw ParsedFunctionEvalError("Error in ParsedScalarFunction evaluation: division by zero");

            case 2: // sqrt error
            case 3: // log error
            case 4: // trigonometric error
              throw ParsedFunctionEvalError("Error in ParsedScalarFunction evaluation: illegal input value");

            case 5: // recursion error
              throw ParsedFunctionEvalError("Error in ParsedScalarFunction evaluation: maximum recursion depth reached");

            default: // ???
              throw ParsedFunctionEvalError("Error in ParsedScalarFunction evaluation: unknown error");
            }
          }

          return val;
        }
      }; // class ParsedVectorFunction::Evaluator<...>
    }; // class ParsedVectorFunction
  } // namespace Analytic
} // namespace FEAT

#endif // FEAT_HAVE_FPARSER
#endif // KERNEL_ANALYTIC_PARSED_FUNCTION_HPP
