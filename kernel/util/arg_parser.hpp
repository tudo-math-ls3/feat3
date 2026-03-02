// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/runtime.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/string.hpp>

#include <algorithm>
#include <any>
#include <deque>
#include <functional>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace FEAT
{
  // Architecture Notes. Read class documentation of ArgParser first.
  // The core of this argument parser implementation consists of four classes
  // - ArgParser is the actual parser
  // - Parameter is the user-facing interface of a parameter
  // - ParameterCore is the inner interface of a parameter used by the ArgParser
  // - ParameterBuilder is a class implementing the builder pattern to define parameters
  //
  // The split into Parameter and ParameterCore is a matter of type-erasure.
  // We want to both present a neat and well-typed interface to the user and
  // automatically handle the parameters in the parser.
  // But the parser can not handle the templated Parameters directly, because
  // they form a non-homogenous collection. To solve this the Parameter class
  // is just a thin wrapper around a type-erased ParameterCore.
  // The parser can handle the ParameterCores, and the Parameter presents a
  // well-typed interface.
  //
  // To achieve this type-erasure the ParameterCore stores its parsed values
  // in a std::any (essentially a void* and a type tag). Any type dependent
  // logic required by the ParameterCore (parsing and validation, for now)
  // is handled by supplying functions that accept std::any to the ParameterCore.
  // These functions are constructed during the parameter definition
  // (when type information is still available), by wrapping the properly
  // typed functions into std::any_casts.
  //
  // Thus the ArgParser holds a list of shared pointers to ParameterCores and has
  // the Parameters as members. The parameters hold a shared pointer to their
  // respective ParameterCore as well. The ArgParser then uses its list of
  // parameter cores to match arguments to parameters during parsing.
  //
  // Start reading at ArgParser::parameter to see how parameters are created and defined.
  // Start reading at ArgParser::parse to see how parsing is implemented.

  namespace Intern
  {
    /// Check if type is some instantiation of std::deque.
    template<typename>
    struct is_std_deque : std::false_type
    {
    };

    /// Check if type is some instantiation of std::deque.
    template<typename T, typename A>
    struct is_std_deque<std::deque<T, A>> : std::true_type
    {
    };

    /// Check if type is some instantiation of std::deque.
    template<typename T>
    constexpr bool is_std_deque_v = is_std_deque<T>::value;

    /// Check if type has an operator>> implementation, i.e. can be parsed by String::parse.
    template<typename T, typename = void>
    struct has_input_operator : std::false_type
    {
    };

    /// Check if type has an operator>> implementation, i.e. can be parsed by String::parse.
    template<typename T>
    struct has_input_operator<T, std::void_t<decltype((std::declval<std::istringstream&>() >> std::declval<T&>()))>>
      : std::true_type
    {
    };

    /// Check if we can create a default parser for some type.
    template<typename T>
    constexpr bool can_default_parse_v = has_input_operator<T>::value;

    /// Check if we can create a default parser for some type.
    template<typename T, typename A>
    constexpr bool can_default_parse_v<std::deque<T, A>> = has_input_operator<T>::value;

    /// Parser for collections. Splits by whitespace, then parses values
    template<typename T_>
    static void parse_collection(std::any& any, const String& s)
    {
      std::deque<T_> collection;
      for(const String& value_str : s.split_by_whitespaces())
      {
        T_ tmp{};
        if(!value_str.parse(tmp))
        {
          throw std::invalid_argument("Parsing failed for string '" + s + "'");
        }
        collection.push_back(tmp);
      }

      any = collection;
    }

    /// Parser for single values
    template<typename T_>
    static void parse_single(std::any& any, const String& s)
    {
      if(!s.parse(std::any_cast<T_&>(any)))
      {
        throw std::invalid_argument("Parsing failed for string '" + s + "'");
      }
    }

    /**
     * \brief Create a type-erased parsing function for a given type
     *
     * \tparam T_ Type to parse to
     *
     * \returns A pointer to a type-erased parsing function
     *
     * If \c T_ is some kind of std::deque, then the parser will split the
     * given string by whitespace, parse each segment individually,
     * construct a deque of the results, and store that deque in the std::any.
     * Otherwise the parser will try to create a value of type \c T_ from
     * the String and write it to the std:.any&.
     *
     * The actual parsing is done by calling String::parse with the appropriate type.
     */
    template<typename T_>
    static void (*make_parser())(std::any&, const String&)
    {
      // If T_ has no input operator then instantiating a parser function for T_ will fail,
      // because String::parse depends on that operator.
      // But we want to give the user a chance to set a custom parser,
      // so we guard the instantiation of the default parser behind
      // can_default_parse_v

      if constexpr(!can_default_parse_v<T_>)
      {
        return nullptr;
      }
      else
      {
        if constexpr(is_std_deque_v<T_>)
        {
          return parse_collection<typename T_::value_type>;
        }
        else
        {
          return parse_single<T_>;
        }
      }
    }

    /**
     * \brief Create a type-erased printing function for a given type
     *
     * \tparam T_ Type to print
     *
     * \returns A pointer to a type-erased printing function
     *
     * The returned function writes the given value to the given
     * ostream. T_ must have an operator<<.
     */
    template<typename T_>
    static void (*make_printer())(std::ostream&, const std::any&)
    {
      if constexpr(is_std_deque_v<T_>)
      {
        return [](std::ostream& s, const std::any& v)
        {
          s << "[";
          for(const auto& i : std::any_cast<const T_&>(v))
          {
            s << i << ", ";
          }
          s << "]";
        };
      }
      else
      {
        return [](std::ostream& s, const std::any& v) { s << std::any_cast<const T_&>(v); };
      }
    }

    /**
     * \brief Enum for different command line parameters
     *
     * - flags are boolean parameters that expect zero arguments
     * - options are parameters that expect exactly one argument
     * - Collection_options are parameters that expect any number of arguments
     */
    enum class ParameterType : std::uint8_t
    {
      flag,
      option,
      collection_option,
    };

    /**
     * \brief Type-erased core of a Parameter
     *
     * This class enables the argument parser to treat all parameters as a homogenous collection.
     * It stores
     * - information for parsing (flags, environment variables, property paths)
     * - information for generating help texts
     * - metadata, like priority and whether the parameter is required to be set by the user
     * - the actual parameter values
     * - function pointers to type-erased function for type specific functionality.
     *
     * The type erased function pointers allow the argument parser to automatically
     * parse the parameters value, validate its correctness if needed, and print the parameter,
     * without needing to know the parameters actual type.
     */
    class ParameterCore
    {
    public:
      /// Type of type-erased parser function
      using ParserType = std::function<void(std::any&, const String&)>;
      /// Type of type-erased validator function
      using ValidatorType = std::function<void(const std::any&)>;
      /// Type of type-erased printer function
      using PrinterType = std::function<void(std::ostream&, const std::any&)>;

      /// Type of "A needs B" declaration with condition
      using Need = std::pair<std::weak_ptr<ParameterCore>, std::function<bool(const std::any&)>>;

      // Parsing information
      ParameterType type = ParameterType::flag;
      String short_flag{""};
      String long_flag{""};
      String environment_variable{""};
      String property_path{""};

      // Strings for usage documentation
      String custom_name;
      String placeholder;
      String help_string;

      // Metadata
      int priority = 0;
      bool required = false;
      bool set_by_user = false;
      std::deque<Need> needs;

      // Raw and parsed values
      std::deque<String> arguments;
      std::any default_value;
      std::any value;

      // Type-erased functions
      ParserType parser;
      ValidatorType validator;
      PrinterType printer;

      /// Constructor
      template<typename T_>
      explicit ParameterCore(T_ def_value) :
        default_value(def_value),
        value(default_value),
        printer(Intern::make_printer<T_>())
      {
        if constexpr(std::is_same_v<T_, bool>)
        {
          // Boolean flags do not consume any arguments beyond the flag itself
          type = ParameterType::flag;
        }
        else if constexpr(Intern::is_std_deque_v<T_>)
        {
          // Collection options append all arguments until the next flag
          type = ParameterType::collection_option;
        }
        else
        {
          type = ParameterType::option;
        }
      }

      /**
       * \brief Tries to create a user-readable name for this parameter
       *
       * Tries to derive a name from the long flag, the short flag, the property path or the environment variable.
       * If all four are empty the name is '<anonymous>'.
       */
      String name() const;

      /// Add an argument to this parameter
      void add_argument(const String& s, int incoming_priority);

      /**
       * \brief Validate this parameter
       *
       * Checks for requiredness, runs the validator if it is set, checks for dependencies between parameters
       */
      void validate(std::deque<String>& errors) const;

      /**
       * \brief Parse this parameter
       *
       * Convert the collected arguments into actual data types
       */
      void parse(std::deque<String>& errors);

      /// Reset this parameter to default
      void reset();
    };

  } // namespace Intern

  template<typename T_>
  class ParameterBuilder;

  /**
   * \brief Parameter of an ArgParser
   *
   * This class acts like a "smart-pointer" giving access to the parsed value for an parameter.
   */
  template<typename T_>
  class Parameter
  {
  private:
    std::shared_ptr<Intern::ParameterCore> _core;

  public:
    template<typename U_>
    friend class ParameterBuilder;

    explicit Parameter(std::shared_ptr<Intern::ParameterCore>&& core) : _core(std::move(core))
    {
    }

    explicit Parameter(std::shared_ptr<Intern::ParameterCore> core) : _core(std::move(core))
    {
    }

    /**
     * \brief Check if argument was set by user
     *
     * \return True, if the argument was explicitly set by the user, false otherwise.
     */
    operator bool() // NOLINT
    {
      return _core->set_by_user;
    }

    /// Accessor for parameter value
    const T_& value() const
    {
      return std::any_cast<const T_&>(_core->value);
    }

    /// Accessor for parameter value
    const T_& operator*() const
    {
      return std::any_cast<const T_&>(_core->value);
    }

    /// Accessor for parameter value
    const T_* operator->() const
    {
      return std::any_cast<const T_>(&(_core->value));
    }
  };

  /// Builder-pattern for defining parameters. See methods for details.
  template<typename T_>
  class ParameterBuilder
  {
  private:
    std::shared_ptr<Intern::ParameterCore> _core;

  public:
    explicit ParameterBuilder(std::shared_ptr<Intern::ParameterCore> core) : _core(std::move(core))
    {
    }

    /**
     * \brief Set the short flag for a parameter
     *
     * \param[in] s Short flag, must be of the form '-x' where x is some letter
     */
    ParameterBuilder& short_flag(String&& s)
    {
      XASSERTM(s.size() == 2, "Short flag must have length two");
      XASSERTM(s[0] == '-', "Short must start with a '-'");
      XASSERTM(s[1] != '-', "Short must not end with a '-'");

      _core->short_flag = std::move(s);
      return *this;
    }

    /**
     * \brief Set the long flag for a parameter
     *
     * \param[in] s Long flag, must be of the form '--foo'. Additional dashes within the flag are allowed, e.g.
     * '--foo-bar'
     */
    ParameterBuilder& long_flag(String&& s)
    {
      _core->long_flag = std::move(s);
      return *this;
    }

    /**
     * \brief Set the environment variable for a parameter
     *
     * \param[in] s Environment variable
     *
     * The parameters value will be taken from this environment variable during parsing,
     * if the environment variable is set.
     */
    ParameterBuilder& env(String&& s)
    {
      _core->environment_variable = std::move(s);
      return *this;
    }

    /**
     * \brief Set the property map path for a parameter
     *
     * \param[in] s PropertyMap path
     *
     * The parameter will be filled with the (parsed) value refered to by this path.
     * The path is always interpreted as relative to the root of the property map given
     * to the ArgParser.
     */
    ParameterBuilder& property(String&& s)
    {
      _core->property_path = std::move(s);
      return *this;
    }

    /**
     * \brief Set help string for this parameter
     *
     * \param[in] s Help string
     *
     * Sets the description of this parameter for the help text
     * shown by ArgParser::usage.
     *
     * The parameter definition
     * \code
     * parameter("mesh").long_flag("--mesh").placeholder("mesh").help_string("Mesh to operate on")
     * \endcode
     * will result in the following snippet in the help text
     * \verbatim
     * --mesh <mesh>
     *    Mesh to operate on
     * \endverbatim
     */
    ParameterBuilder& help_string(String&& s)
    {
      s.trim();
      _core->help_string = std::move(s);
      return *this;
    }

    /**
     * \brief Set a placeholder text for the argument of this parameter
     *
     * \param[in] s Placeholder
     *
     * The placeholder text will be shown in place of the parameters arguments
     * in the help text generated by ArgParser::usage.
     *
     * The parameter definition
     * \code
     * parameter("mesh").long_flag("--mesh").placeholder("mesh").help_string("Mesh to operate on")
     * \endcode
     * will result in the following snippet in the help text
     * \verbatim
     * --mesh <mesh>
     *    Mesh to operate on
     * \endverbatim
     */
    ParameterBuilder& placeholder(String&& s)
    {
      _core->placeholder = std::move(s);
      return *this;
    }

    /**
     * \brief Set a name for this parameter
     *
     * \param[in] s Name for parameter
     *
     * The name will be used in error messages and the ArgParser::display() method
     * to identify the parameter. If no name is set, the name will be derived from either
     * a flag, the property path, or the environment variable.
     * If none of these are set the name is '<anonymous>'.
     */
    ParameterBuilder& name(String&& s)
    {
      _core->custom_name = std::move(s);
      return *this;
    }

    /**
     * \brief Set a custom parser for this parameter
     *
     * \param[in] parser_fn Pointer to parser function
     *
     * Allows setting a custom ad-hoc parser for this parameter,
     * rather than relying on the usual String::parse method.
     *
     * The parsing function is expected to throw a std::invalid_argument
     * exception if any errors occur. The message given in the exception
     * will be shown to the user.
     */
    ParameterBuilder& parser(T_ (*parser_fn)(const String&))
    {
      // SAFETY:
      // Functions exist for the lifetime of the program.
      // Thus a valid function pointer will never dangle
      // and calling this lambda later is always safe
      const auto lambda = [parser_fn](std::any& any, const String& s) { any = parser_fn(s); };
      _core->parser = lambda;

      return *this;
    }

    /**
     * \brief Set a validator for this paramater
     *
     * \param[in] validator_fn Pointer to validator function
     *
     * Validators are run after argument parsing to ensure
     * that parameters fulfill additional properties beyond
     * just their type. This allows you to ensure that, for
     * example, a integer parameter is within a specific
     * of number, or that a file at a given file path actually exists.
     *
     * This specific overload of this method expects a validator function
     * that throws a std::invalid_argument exception if the given value
     * does not pass validation. The message given to the
     * std::invalid_argument exception is forwarded to the user.
     */
    ParameterBuilder& validator(void (*validator_fn)(const T_&))
    {
      XASSERTM(validator_fn != nullptr, "Validator must be valid function pointer!");

      // NOTE:
      // We accept a function pointer here, because it
      // forbids the user from passing us a lambda-expression
      // that captures something from its environment

      // SAFETY:
      // Functions exist for the lifetime of the program.
      // Thus a valid function pointer will never dangle
      // and calling this lambda later is always safe
      const auto lambda = [validator_fn](const std::any& any) { validator_fn(std::any_cast<const T_&>(any)); };
      _core->validator = lambda;

      return *this;
    }

    /**
     * \brief Set a validator for this paramater
     *
     * \param[in] validator_fn Pointer to validator function
     *
     * Validators are run after argument parsing to ensure
     * that parameters fulfill additional properties beyond
     * just their type. This allows you to ensure that, for
     * example, a integer parameter is within a specific
     * of number, or that a file at a given file path actually exists.
     *
     * This specific overload of this method expects a validator function
     * that returns a boolean indicating whether the given value passes
     * validation. A failed validation produces an error with a generic
     * "Failed boolean test" message. If you want to supply a more
     * useful error message to your users, see the non-boolean overload
     * of this function.
     */
    ParameterBuilder& validator(bool (*validator_fn)(const T_&))
    {
      XASSERTM(validator_fn != nullptr, "Validator must be valid function pointer!");

      // NOTE:
      // We accept a function pointer here, because it
      // forbids the user from passing us a lambda-expression
      // that captures something from its environment

      // SAFETY:
      // Functions exist for the lifetime of the program.
      // Thus a valid function pointer will never dangle
      // and calling this lambda later is always safe
      const auto lambda = [validator_fn](const std::any& any)
      {
        if(!validator_fn(std::any_cast<const T_&>(any)))
        {
          throw std::invalid_argument("Failed boolean test");
        }
      };
      _core->validator = lambda;

      return *this;
    }

    /**
     * \brief Mark this parameter as required
     *
     * A required parameter must be set by the user. It is
     * not enough for it to be default initialized.
     */
    ParameterBuilder& required()
    {
      _core->required = true;

      return *this;
    }

    /**
     * \brief Mark another parameter as needed by this one
     *
     * \param[in] parameter The parameter needed by this parameter
     *
     * Indicates that if this parameter is set by the user, then the given
     * parameter must be set by the user as well.
     *
     * Use this option if one parameter is not meaningful without another parameter.
     * A filename for example might be meaningless without setting an output directory.
     *
     * Alternatively consider combining the data into a single type with a custom
     * parser and combining the dependent parameters into a single parameter.
     *
     * \warning The parameter given to this function must already be fully defined.
     * In other words, it must come before the current parameter in the parser struct.
     *
     * \sa ParameterBuilder::needs_if
     */
    template<typename U_>
    ParameterBuilder& needs(const Parameter<U_>& parameter)
    {
      XASSERTM(parameter._core != nullptr, "Needed parameter must be defined already!");

      _core->needs.emplace_back(parameter._core, [](const std::any& /*unused*/) { return true; });

      return *this;
    }

    /**
     * \brief Conditionally mark another parameter as needed by this one
     *
     * \param[in] parameter The parameter needed by this parameter
     * \param[in] condition Function pointer to predicate
     *
     * Indicates that if this parameter is set by the user _and_ has a certain value,
     * as indicated by the condition function, then the given parameter must be set by
     * the user as well.
     *
     * Use this option if one parameter is sometimes not meaningful without another parameter.
     * For example, if this parameter chooses a solver to use,
     * then the different solver settings are only required if that solver is actually chosen.
     *
     * \code
     * parameter("Umfpack").long_flag("--solver").needs_if(max-iters, [](const auto& self) { self == "cg"; })
     * \endcode
     *
     * \warning The parameter given to this function must already be fully defined.
     * In other words, it must come before the current parameter in the parser struct.
     */
    template<typename U_>
    ParameterBuilder& needs_if(const Parameter<U_>& parameter, bool (*condition)(const T_&))
    {
      XASSERTM(parameter._core != nullptr, "Needed parameter must be defined already!");

      _core->needs.emplace_back(
        parameter._core,
        [condition](const std::any& value) { return condition(std::any_cast<const T_&>(value)); });

      return *this;
    }

    // NOTE: This conversion operator exists to convert the ParameterBuilder produced by ArgParser::parameter
    // into a parameter without forcing the user to call an additional function for every parameter.
    // For that reason it is deliberately implicit
    /// Conversion operator
    operator Parameter<T_>() // NOLINT
    {
      if(_core->placeholder.empty())
      {
        if(!_core->long_flag.empty())
        {
          _core->placeholder = _core->long_flag.substr(2);
        }
        else if(!_core->short_flag.empty())
        {
          _core->placeholder = _core->short_flag.substr(1);
        }
      }

      if(!_core->parser)
      {
        _core->parser = Intern::make_parser<T_>();
      }
      return Parameter<T_>(_core);
    }
  };

  /**
   * \brief Argument parser base class
   *
   * This class aims to be a more convenient alternative to the SimpleArgParser.
   * It allows defining fully integrated argument parsers that can:
   * - parse arguments from the command line, property maps, and environment variables
   * - store parsed arguments directly into a convenient struct
   * - generate a help text
   * - handle required and optional parameters
   * - handle custom parsing
   * - handle validation of arguments
   * - model dependencies between parameters.
   *
   * **Example**
   * \code
    // Define parser
    struct ExampleParser : public ArgParser
    {
      Parameter<String> mesh_path = parameter(String("mesh.xml")).short_flag("-m").long_flag("--mesh").required();

      Parameter<bool> write_vtk = parameter(false).long_flag("--vtk");
      Parameter<String> output_file = parameter(String("output.vtu")).long_flag("--output-file").needs(write_vtk);

      Parameter<std::uint32_t> max_iters = parameter(std::uint32_t(0)).long_flag("--max-iters").property("max-iters");
      Parameter<String> solver = parameter(String("Umfpack"))
                                   .long_flag("--solver")
                                   .property("solver")
                                   .validator([](const auto& s) { return s == "Umfpack" || s == "CG"; })
                                   .needs_if(max_iters, [](const auto& s) { return s == "CG"; });
    }

    // Use like this
    ExampleParser args;
    bool success = args.parse(argc, argv, property_map);

    if(!success)
    {
      for(const String& error : parser.errors())
      {
        std::cout << error << "\n";
      }
      Runtime::abort();
    }

    if(*args.write_vtk)
    {
      // Export VTK
    }
   * \endcode
   *
   * **Creating a parser**
   *
   * You create a parser by defining a class or struct that inherits from this base class.
   * Inheriting from this class makes a protected method ArgParser::parameter available to you
   * that allows you to define Parameters via a ParameterBuilder.
   * See the various methods of the ParameterBuilder for details on the different settings
   * you can use to define your parser.
   *
   * **Using a parser**
   *
   * After you have defined your parser, you can use it by simply creating an instance of it
   * and calling the parse method. The parse method returns a boolean indicating whether the
   * parameter parsing was successful. If parsing was unsuccessful the parser will contain a
   * list of errors that occured during parsing.
   * Note that for a better user experience the parser is optimistic and will try to produce
   * as many errors as possible rather than quitting after the first warning.
   *
   * After successfull parsing you can access the parsed values via the Parameter::value() method
   * or the overloads of operator* and operator-> on the Parameter class.
   *
   * The parser supports being reused. If you call parse again with different parameters,
   * then the parsed result will be identical to calling parse on a fresh parser.
   *
   * **Parsing priority**
   *
   * The parser considers, if provided, command line arguments, property map entries, environment
   * variables, and the supplied default values of parameters. During parsing the priority for
   * sources is as follows:
   *
   * Default Values < environment variables < property maps < command line arguments
   *
   * Roughly, the more explicit and interactive the method, to higher the prioritiy.
   * This means that setting a parameter via the command line will overide a value
   * for that parameter set in a property map.
   *
   * **Allowed types**
   *
   * The ArgParser supports all types that can be parsed via the String::parse method.
   * To add support for your own custom type, either supply an operator>> overload for your type,
   * which will be used by String::parse or supply a custom parser in the parameter definition.
   *
   * Boolean parameters are handled specially. Environment variables and properties are parsed as normal,
   * but setting the command line flag sets the parameters value to the inverse of the default value.
   * This allows you to define "--no-foo" flags that are on by default and get turned off by setting the flag.
   *
   * If a type T can be parsed, then the Parser also supports Parameters of type std::deque<T>.
   * These parameters are parsed by either splitting the value of a environment variable or property map
   * index at whitespace and parsing each value individually, or by parsing all following command line
   * arguments until the next flag occurs.
   *
   * **Help Text**
   *
   * You can create a help text with a list of all parameters with their flags, property map paths,
   * and environment variables via the ArgParser::help method.
   * The ParameterBuilder provides various settings to customize this help text for each parameter.
   *
   * **Display Text**
   *
   * You can create a table containing the values of all parameters and their respective sources
   * via the ArgParser::display method.
   *
   * **Errors**
   *
   * In the interest of correctness the parser treats most unexpected things as errors.
   * This includes:
   * - unused command line arguments
   * - failed calls to String::parse
   * - unset required parameters
   * - failed validators
   * - unfulfilled dependencies between parameters
   */
  class ArgParser
  {
  private:
    std::deque<std::shared_ptr<Intern::ParameterCore>> _parameters;
    std::deque<String> _errors;

    String _program{""};
    String _description{""};

    static constexpr int priority_default = 0;
    static constexpr int priority_env = 1;
    static constexpr int priority_property = 2;
    static constexpr int priority_cli = 3;

  public:
    const std::deque<String>& errors() const
    {
      return _errors;
    }

    /**
     * \brief Set program description
     *
     * \param[in] description Description
     *
     * This text will be shown before the usage information in the help text.
     */
    void set_description(String&& description)
    {
      _description = std::move(description);
    }

    /**
     * \brief Generate help text for this parser
     *
     * \returns A string containing the help text for this parser
     *
     * The help text gives information about all parameters supported by this parser
     */
    [[nodiscard]] String help_text()
    {
      // Sort parameters alphabetically
      std::sort(
        _parameters.begin(),
        _parameters.end(),
        [](auto& a, auto& b) { return a->name().compare_no_case(b->name()) <= 0; });

      std::stringstream stream;

      if(!_description.empty())
      {
        stream << _description << "\n\n";
      }
      stream << "Usage: " << _program << " [OPTIONS]\n\n";
      stream << "Options:\n";

      for(const auto& core : _parameters)
      {
        parameter_help(*core, stream);
      }

      return stream.str();
    }

    /**
     * \brief Produce formatted output of all values in this parser
     *
     * \returns A string containing the values of all parameters
     *
     * \note This function shows the arguments given by the user.
     * Parser might change these values however they like.
     */
    [[nodiscard]] String display()
    {
      // Sort parameters alphabetically
      std::sort(
        _parameters.begin(),
        _parameters.end(),
        [](auto& a, auto& b) { return a->name().compare_no_case(b->name()) <= 0; });

      std::stringstream stream;

      // Print header
      stream << String("Parameter").pad_back(40, ' ');
      stream << "|";
      stream << String("Source").pad_back(10, ' ');
      stream << "|";
      stream << "Value\n";
      stream << String("").pad_back(75, '-') << "\n";

      for(const auto& core : _parameters)
      {
        XASSERT(core->value.has_value());

        stream << core->name().pad_back(40, '.');
        stream << "|";

        switch(core->priority)
        {
        case priority_default:  stream << "default   |"; break;
        case priority_env:      stream << "env       |"; break;
        case priority_property: stream << "property  |"; break;
        case priority_cli:      stream << "cli       |"; break;
        default:                stream << "          |";
        }

        if(core->type == Intern::ParameterType::flag)
        {
          stream << (std::any_cast<const bool&>(core->value) ? "true" : "false");
        }
        else
        {
          XASSERT(bool(core->printer));
          core->printer(stream, core->value);
        }

        stream << "\n";
      }

      return stream.str();
    }

    /**
     * \brief Populate this parser
     *
     * \param[in] argc Number of command line arguments, as passed to main
     * \param[in] argv Command line arguments
     * \param[in] property_map Pointer to property map root
     *
     * Collect all parameters from the given sources, parse values, validate parser.
     * Note that the parser expects the first command line parameter to be the path
     * to the program binary, as is usual.
     */
    [[nodiscard]] bool parse(int argc, const char* const* argv, const PropertyMap* property_map = nullptr)
    {
      // Reset errors
      _errors.clear();

      if(!_check_for_unique_parameters())
      {
        return false;
      }

      // Set default values, reset priorities, clear arguments
      for(auto& core : _parameters)
      {
        core->reset();
      }

      // Collect values in environment variables
      _collect_env();

      // Collect values in property map
      if(property_map != nullptr)
      {
        _collect_property_map(property_map);
      }

      // Collect command line args
      _collect_args(argc, argv);

      // Parse collected data types into proper data types
      for(auto& core : _parameters)
      {
        core->parse(_errors);
      }

      _validate();

      return _errors.empty();
    }

    /**
     * \brief Populate this parser
     *
     * \param[in] property_map Reference to property map root
     *
     * Collect all parameters from the environment and given property map.
     */
    [[nodiscard]] bool parse(const PropertyMap& property_map)
    {
      return parse(0, nullptr, &property_map);
    }

    /**
     * \brief Populate this parser
     *
     * \param[in] argc Number of command line arguments, as passed to main
     * \param[in] argv Command line arguments
     * \param[in] property_map_path Path to file containing property map
     *
     * Collect all parameters from the given sources, parse values, validate parser.
     * Note that the parser expects the first command line parameter to be the path
     * to the program binary, as is usual.
     */
    [[nodiscard]] bool parse(int argc, const char* const* argv, const String& property_map_path)
    {
      PropertyMap pmap;
      pmap.read(property_map_path);
      return parse(argc, argv, &pmap);
    }

    /**
     * \brief Populate this parser
     *
     * \param[in] property_map_path Path to file containing property map
     *
     * Collect all parameters from the environment and given property map.
     */
    [[nodiscard]] bool parse(const String& property_map_path)
    {
      PropertyMap pmap;
      pmap.read(property_map_path);
      return parse(0, nullptr, &pmap);
    }

  protected:
    /**
     * \brief Add a Parameter to this parser
     *
     * \param[in] default_value Default value for this parameter
     *
     * \returns A ParameterBuilder for the new Parameter
     */
    template<typename T_>
    ParameterBuilder<T_> parameter(T_&& default_value = T_{})
    {
      // Create type-erased parameter core
      std::shared_ptr<Intern::ParameterCore> core = std::make_shared<Intern::ParameterCore>(std::forward<T_>(default_value));

      // Store type-erased core for parsing purposes.
      // This allows the parser to treat all parameters the same,
      // irrespective of the types of their values.
      _parameters.push_back(core);

      // Return builder to user for defining the parameter fully
      return ParameterBuilder<T_>(std::move(core));
    }

  private:
    /// Search for short flags, long flags, environment variables, or property paths that get used by multiple
    /// parameters
    bool _check_for_unique_parameters()
    {
      std::set<String> short_flags;
      std::set<String> long_flags;
      std::set<String> environment_variables;
      std::set<String> property_paths;

      bool result = true;

      auto handle_flag = [&](std::set<String>& set, const String& s)
      {
        if(s.empty())
        {
          return;
        }

        if(set.find(s) != set.end())
        {
          _errors.push_back("Error: Flag " + s + " is used by multiple parameters");
          result = false;
        }
        else
        {
          set.insert(s);
        }
      };

      for(auto& core : _parameters)
      {
        handle_flag(short_flags, core->short_flag);
        handle_flag(long_flags, core->long_flag);
        handle_flag(environment_variables, core->environment_variable);
        handle_flag(property_paths, core->property_path);
      }

      return result;
    }

    void _collect_env()
    {
      for(auto& core : _parameters)
      {
        const char* value = std::getenv(core->environment_variable.c_str());
        if(value != nullptr)
        {
          String s(value);
          core->add_argument(s, priority_env);
        }
      }
    }

    void _collect_property_map(const PropertyMap* property_map)
    {
      for(auto& core : _parameters)
      {
        auto [string, success] = property_map->query(core->property_path);
        if(success)
        {
          core->add_argument(string, priority_property);
        }
      }
    }

    void _collect_args(int argc, const char* const* argv)
    {
      if(argc == 0)
      {
        return;
      }

      _program = String(argv[0]);

      std::deque<String> args;

      for(int i(1); i < argc; i++)
      {
        args.emplace_back(argv[i]);
      }

      Index idx(0);
      while(idx < args.size())
      {
        const String& flag = args[idx++];

        if(!is_flag(flag))
        {
          XABORTM("Error: Expected flag in position " + stringify(idx));
        }

        bool handled = false;
        for(auto& core : _parameters)
        {
          if(flag.compare_no_case(core->short_flag) != 0 && flag.compare_no_case(core->long_flag) != 0)
          {
            continue;
          }

          handled = true;

          if(core->type == Intern::ParameterType::flag)
          {
            // Flag was given. Invert its state.
            // any_cast is safe because flags must be bool parameters
            core->value = !std::any_cast<const bool&>(core->default_value);
            core->set_by_user = true;
          }
          else if(core->type == Intern::ParameterType::option)
          {
            // Normal option. Consume a single argument
            core->add_argument(args[idx++], priority_cli);
          }
          else if(core->type == Intern::ParameterType::collection_option)
          {
            // Collect arguments until next flag
            while(idx < args.size() && !is_flag(args[idx]))
            {
              core->add_argument(args[idx++], priority_cli);
            }
          }
        }

        if(!handled)
        {
          // We found no parameter for this flag. Generate error
          _errors.push_back("Error: Unknown flag " + flag);

          // Then advance to next flag
          while(idx < args.size() && !is_flag(args[idx]))
          {
            idx++;
          }
        }
      }
    }

    /// Validate parameters after parsing
    void _validate()
    {
      for(auto& core : _parameters)
      {
        core->validate(_errors);
      }
    }

    /// Write flags, property path, environment variable and help text for given parameter to stringstream
    static void parameter_help(const Intern::ParameterCore& core, std::stringstream& s)
    {
      if(
        core.short_flag.empty() && core.long_flag.empty() && core.property_path.empty() &&
        core.environment_variable.empty())
      {
        return;
      }

      s << "\t";

      auto format_flag = [&s](Intern::ParameterType type, const String& flag, const String& placeholder)
      {
        if(type == Intern::ParameterType::flag)
        {
          s << flag;
        }
        if(type == Intern::ParameterType::option)
        {
          s << flag << " <" << placeholder << ">";
        }
        if(type == Intern::ParameterType::collection_option)
        {
          s << flag << " <" << placeholder << "s...>";
        }
      };

      auto format_other = [&s](const std::string_view prefix, const String& value)
      { s << "[" << prefix << ":" << value << "]"; };

      // -a,  --arg, [property:arg], [env:FEAT_ARG]

      bool has_prev = false;
      if(!core.short_flag.empty())
      {
        format_flag(core.type, core.short_flag, core.placeholder);
        has_prev = true;
      }

      if(!core.long_flag.empty())
      {
        if(has_prev)
        {
          s << ", ";
        }

        format_flag(core.type, core.long_flag, core.placeholder);

        has_prev = true;
      }

      if(!core.property_path.empty())
      {
        if(has_prev)
        {
          s << ", ";
        }

        format_other("property", core.property_path);

        has_prev = true;
      }

      if(!core.environment_variable.empty())
      {
        if(has_prev)
        {
          s << ", ";
        }

        format_other("env", core.environment_variable);

        has_prev = true;
      }

      if(has_prev)
      {
        s << "\n";
      }

      if(core.required)
      {
        s << "\t\tRequired argument.\n";
      }

      for(const String& line : core.help_string.split_by_charset("\n"))
      {
        s << "\t\t" << line.trim() << "\n";
      }

      s << "\n";
    }

    static bool is_flag(const String& s)
    {
      return !s.empty() && s[0] == '-';
    }
  };
} // namespace FEAT
