// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/runtime.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/string.hpp>

#include <any>
#include <deque>
#include <functional>
#include <iostream>
#include <memory>
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

    /// Check if type is some instantiation of std::map
    template<typename>
    struct is_std_map : std::false_type
    {
    };

    /// Check if type is some instantiation of std::map
    template<typename K, typename V, typename C, typename A>
    struct is_std_map<std::map<K, V, C, A>> : std::true_type
    {
    };

    /// Check if type is some instantiation of std::map
    template<typename T>
    constexpr bool is_std_map_v = is_std_map<T>::value;

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

    /// Check if we can create a default parser for some type.
    template<typename K, typename V, typename C, typename A>
    constexpr bool can_default_parse_v<std::map<K, V, C, A>> =
      has_input_operator<K>::value && has_input_operator<V>::value;

    /// Parser for collections. Splits by whitespace, then parses values
    template<typename T_>
    static void parse_deque(std::any& any, const String& s)
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

    /// Parser for maps. Splits by whitespace, treats adjacent strings as key-value pairs.
    template<typename K_, typename V_>
    static void parse_map(std::any& any, const String& s)
    {
      std::map<K_, V_> map;
      auto strs = s.split_by_whitespaces();

      if(!(strs.size() % 2 == 1))
      {
        throw std::invalid_argument("Parsing failed for string '" + s + "': Odd number of strings in key-value pairs.");
      }

      for(Index i(0); i < strs.size(); i += 2)
      {
        K_ key{};
        V_ value{};

        if(!strs[i].parse(key))
        {
          throw std::invalid_argument("Parsing failed for string '" + s + "': Could not parse key '" + strs[i] + "'");
        }

        if(!strs[i + 1].parse(value))
        {
          throw std::invalid_argument("Parsing failed for string '" + s + "': Could not parse value '" + strs[i] + "'");
        }

        map.insert({key, value});
      }

      any = map;
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
     * If \c T_ is some kind of std::map, then the parser will split the given
     * string by whitespace, parse keys and values in alternating order, place
     * the resulting key-value pairs in a std::map, and store that map in the
     * std::any.
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
          return parse_deque<typename T_::value_type>;
        }
        else if constexpr(is_std_map_v<T_>)
        {
          return parse_map<typename T_::key_type, typename T_::mapped_type>;
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
      else if constexpr(is_std_map_v<T_>)
      {
        return [](std::ostream& s, const std::any& v)
        {
          s << "{";
          for(const auto& [key, value] : std::any_cast<const T_&>(v))
          {
            s << key << ": " << value << ", ";
          }
          s << "}";
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
      /// Type of type-erased deriver function
      using DeriverType = std::function<void(std::any&, const std::deque<std::reference_wrapper<std::any>>&)>;
      /// Type of type-erased predicate function
      using PredicateType = std::function<bool(const std::any&)>;

      /// Type of "A needs B" declaration with condition
      using Need = std::pair<std::weak_ptr<ParameterCore>, PredicateType>;

      /// Type of "Hide A if B" declaration
      using HideIf = std::pair<std::weak_ptr<ParameterCore>, PredicateType>;

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
      bool bootstrap = false;
      bool set_by_user = false;
      bool hidden = false;
      std::deque<Need> needs;
      std::deque<HideIf> hide_if;
      std::deque<std::weak_ptr<ParameterCore>> derives_from;

      // Raw and parsed values
      std::deque<String> arguments;
      std::any default_value;
      std::any value;

      // Type-erased functions
      ParserType parser;
      ValidatorType validator;
      PrinterType printer;
      DeriverType deriver;

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
        else if constexpr(Intern::is_std_deque_v<T_> || Intern::is_std_map_v<T_>)
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

      /**
       * \brief Sets this parameters value based on other parameters
       *
       * Uses the parameters and function given by the user in
       * ParameterBuilder::derives_from() to set this parameters value. Does
       * nothing if no derives_from has been set.
       */
      void derive_self();

      /**
       * \brief Hides other parameters based on this parameters hide_if declarations
       */
      void hide();

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
    bool set_by_user() const
    {
      return _core->set_by_user;
    }

    /// Accessor for parameter value
    const T_& value() const
    {
      XASSERTM(
        !_core->hidden,
        "Trying to access hidden parameter " + _core->name() + ". Check your parameter definitions.");
      return std::any_cast<const T_&>(_core->value);
    }

    /// Accessor for parameter value
    const T_& operator*() const
    {
      XASSERTM(
        !_core->hidden,
        "Trying to access hidden parameter " + _core->name() + ". Check your parameter definitions.");
      return std::any_cast<const T_&>(_core->value);
    }

    /// Accessor for parameter value
    const T_* operator->() const
    {
      XASSERTM(
        !_core->hidden,
        "Trying to access hidden parameter " + _core->name() + ". Check your parameter definitions.");
      return std::any_cast<const T_>(&(_core->value));
    }
  };

  /// Builder-pattern for defining parameters. See methods for details.
  template<typename T_>
  class ParameterBuilder
  {
  private:
    std::shared_ptr<Intern::ParameterCore> _core;

    String _env_prefix;
    String _property_path_prefix;
    String _flag_prefix;

  public:
    explicit ParameterBuilder(
      std::shared_ptr<Intern::ParameterCore> core,
      String env_prefix,
      String property_path_prefix,
      String flag_prefix) :
      _core(std::move(core)),
      _env_prefix(std::move(env_prefix)),
      _property_path_prefix(std::move(property_path_prefix)),
      _flag_prefix(std::move(flag_prefix))
    {
    }

    /**
     * \brief Set the short flag for a parameter
     *
     * \param[in] s Short flag, can be of form "x" or "-x" where x is some letter.
     */
    ParameterBuilder& short_flag(String&& s)
    {
      // Prepend dash if required
      if(!s.starts_with("-"))
      {
        s = "-" + s;
      }

      XASSERTM(s.size() == 2, "Short flag must have length two");
      XASSERTM(s[0] == '-', "Short must start with a '-'");
      XASSERTM(s[1] != '-', "Short must not end with a '-'");

      _core->short_flag = std::move(s);
      return *this;
    }

    /**
     * \brief Set the long flag for a parameter
     *
     * \param[in] s long flag, can be of form "foo" or "--foo". Additional
     * dashes within the flag are allowed, e.g. '--foo-bar'
     */
    ParameterBuilder& long_flag(String s)
    {
      XASSERTM(!s.empty(), "Long flag must not be empty");

      // Remove leading dashes
      if(s.starts_with("--"))
      {
        s = s.substr(2);
      }

      // Append flag prefix if applicable
      if(!_flag_prefix.empty())
      {
        s = _flag_prefix + "-" + s;
      }

      _core->long_flag = "--" + s;
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
    ParameterBuilder& env(const String& s)
    {
      if(!_env_prefix.empty())
      {
        _core->environment_variable = _env_prefix + "_" + s;
      }
      else
      {
        _core->environment_variable = std::move(s);
      }
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
     *
     * \sa Path specification in PropertyMap
     */
    ParameterBuilder& property(String&& s)
    {
      if(!_property_path_prefix.empty())
      {
        _core->property_path = _property_path_prefix + "/" + s;
      }
      else
      {
        _core->property_path = std::move(s);
      }
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

    /**
     * \brief Conditionally mark another parameter as hidden
     *
     * \param[in] parameter The parameter to be conditionally hidden
     * \param[in] condition Function pointer to predicate
     *
     * Indicates that if this parameter is set by the user _and_ has a certain
     * value, as indicated by the condition function, then the given parameter
     * should be hidden. A parameter is hidden as soon as any of the hide_if
     * conditions regarding the parameter are fulfilled.
     *
     * To ensure that no parameter is accidentally hidden, accessing the value
     * of hidden parameter is an error and aborts the program.
     *
     * Use this option if one parameter is unused, depending on another parameters value.
     *
     * \warning The parameter given to this function must already be fully defined.
     * In other words, it must come before the current parameter in the parser struct.
     */
    template<typename U_>
    ParameterBuilder& hide_if(const Parameter<U_>& parameter, bool (*condition)(const T_&))
    {
      XASSERTM(parameter._core != nullptr, "Needed parameter must be defined already!");

      _core->hide_if.emplace_back(
        parameter._core,
        [condition](const std::any& value) { return condition(std::any_cast<const T_&>(value)); });

      return *this;
    }

    /**
     * \brief Mark this parameter as a bootstrap parameter
     *
     * A common FEAT idiom is parsing a path (or paths) to a property map from
     * the command line, loading that property map, and then parsing remaining
     * options from that map. Naively, this would required defining at least two
     * parsers, where the first parser just determines the paths of the property
     * maps and the second parser does the full parsing run.
     * To avoid this parameters can be marked as bootstrap parameters.
     * The ArgParser::bootstrap method will parse just those parameters, with no
     * further validation, and assemble a merged property map from the given
     * paths. With that merged property map the full parsing can then be run.
     *
     * \warning Bootstrap parameters must be of type std::deque<String>.
     */
    template<typename V_ = T_, typename = std::enable_if_t<std::is_same_v<std::deque<String>, V_>>>
    ParameterBuilder& bootstrap()
    {
      _core->bootstrap = true;

      return *this;
    }

    /**
     * \brief Derive this parameter from another parameter if it has not been set by the user
     *
     * \param[in] param Parameter this parameter derives from
     * \param[in] fn Function which calculates value from other parameter.
     *
     * Indicates that this parameter should have its value set based on the
     * value of one other parameter, if this parameter has not been set by the
     * user. If this parameter has been set by the user, the value set by the
     * user will be kept. A parameter that has no flags, properties, or
     * environment variables associated with it will always be derived.
     */
    template<typename U_>
    ParameterBuilder& derives_from(Parameter<U_>& param, T_ (*fn)(const U_&))
    {
      XASSERTM(param._core != nullptr, "Derivee needs to be initialized");

      // NOTE(mmuegge): It is totally possible to abstract over derives_from_1
      // and derives_from_2: The easiest way would probably be to convert the
      // std::deque passed into the deriver into a tuple of the same length as
      // the number of parameters this parameter is derived from, then cast the
      // tuple to the proper types via a fold-expression, then use std::apply to
      // call fn with the arguments from the tuple. All of that code is
      // possible, but just plain ugly and hard to understand. So I choose the
      // easy route for now.

      _core->derives_from.clear();
      _core->derives_from.emplace_back(param._core);
      _core->deriver = [fn](std::any& out, const std::deque<std::reference_wrapper<std::any>>& args)
      {
        XASSERT(args.size() >= 1);
        out = fn(std::any_cast<const U_&>(args[0]));
      };

      return *this;
    }

    /// ArgParser::derive_from wrapper for casting lambdas to function pointers
    template<typename U_, typename Fn_>
    ParameterBuilder& derives_from(Parameter<U_>& param, Fn_&& fn)
    {
      return derives_from(param, static_cast<T_ (*)(const U_&)>(std::forward<Fn_>(fn)));
    }

    /**
     * \brief Derive this parameter from two other parameters if it has not been set by the user
     *
     * \param[in] param_one First parameter this parameter derives from
     * \param[in] param_two Second parameter this parameter derives from
     * \param[in] fn Function which calculates value from other parameter.
     *
     * Indicates that this parameter should have its value set based on the
     * value of two other parameters, if this parameter has not been set by the
     * user. If this parameter has been set by the user, the value set by the
     * user will be kept. A parameter that has no flags, properties, or
     * environment variables associated with it will always be derived.
     */
    template<typename U_, typename V_>
    ParameterBuilder& derives_from(Parameter<U_>& param_one, Parameter<V_>& param_two, T_ (*fn)(const U_&, const V_&))
    {
      XASSERTM(param_one._core != nullptr, "Derivee needs to be initialized");
      XASSERTM(param_two._core != nullptr, "Derivee needs to be initialized");

      _core->derives_from.clear();
      _core->derives_from.emplace_back(param_one._core);
      _core->derives_from.emplace_back(param_two._core);
      _core->deriver = [fn](std::any& out, const std::deque<std::reference_wrapper<std::any>>& args)
      {
        XASSERT(args.size() >= 2);
        out = fn(std::any_cast<const U_&>(args[0]), std::any_cast<const V_&>(args[1]));
      };

      return *this;
    }

    /// ArgParser::derive_from wrapper for casting lambdas to function pointers
    template<typename U_, typename V_, typename Fn_>
    ParameterBuilder& derives_from(Parameter<U_>& param_one, Parameter<V_>& param_two, Fn_&& fn)
    {
      return derives_from(param_one, param_two, static_cast<T_ (*)(const U_&, const V_&)>(std::forward<Fn_>(fn)));
    }

    /**
     * \brief Derive this parameter from three other parameters if it has not been set by the user
     *
     * \param[in] param_one First parameter this parameter derives from
     * \param[in] param_two Second parameter this parameter derives from
     * \param[in] param_three Third parameter this parameter derives from
     * \param[in] fn Function which calculates value from other parameter.
     *
     * Indicates that this parameter should have its value set based on the
     * value of three other parameters, if this parameter has not been set by the
     * user. If this parameter has been set by the user, the value set by the
     * user will be kept. A parameter that has no flags, properties, or
     * environment variables associated with it will always be derived.
     */
    template<typename U_, typename V_, typename W_>
    ParameterBuilder& derives_from(
      const Parameter<U_>& param_one,
      const Parameter<V_>& param_two,
      const Parameter<W_>& param_three,
      T_ (*fn)(const U_&, const V_&))
    {
      XASSERTM(param_one._core != nullptr, "Derivee needs to be initialized");
      XASSERTM(param_two._core != nullptr, "Derivee needs to be initialized");
      XASSERTM(param_three._core != nullptr, "Derivee needs to be initialized");

      _core->derives_from.clear();
      _core->derives_from.emplace_back(param_one._core);
      _core->derives_from.emplace_back(param_two._core);
      _core->derives_from.emplace_back(param_three._core);
      _core->deriver = [fn](std::any& out, const std::deque<std::reference_wrapper<std::any>>& args)
      {
        XASSERT(args.size() >= 3);
        out =
          fn(std::any_cast<const U_&>(args[0]), std::any_cast<const V_&>(args[1]), std::any_cast<const W_&>(args[2]));
      };

      return *this;
    }

    /// ArgParser::derive_from wrapper for casting lambdas to function pointers
    template<typename U_, typename V_, typename W_, typename Fn_>
    ParameterBuilder&
    derives_from(Parameter<U_>& param_one, Parameter<V_>& param_two, Parameter<W_>& param_three, Fn_&& fn)
    {
      return derives_from(
        param_one,
        param_two,
        param_three,
        static_cast<T_ (*)(const U_&, const V_&, const W_&)>(std::forward<Fn_>(fn)));
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
   * If types K and V can be parsed, then the parser also supports Parameters of
   * type std::map<K, V>. These parameters are parsed similarly to std::deques,
   * except that pairs of arguments are treated as key-value pairs for the
   * std::map.
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
   *
   * **Composing parsers**
   *
   * You might have a parser A with several parameters and want to create a
   * parser B with the same parameters as A, plus some additional ones. There
   * are two mechanisms for composing parsers like that. Inheritance and
   * sub-parsers.
   *
   * Inheritance allows you to compose parser in a "flat" manner. You can have
   * parser B inherit from A, then B will naturally have all the same parameters
   * as A, plus any additional ones you define in the class definition of B. All
   * these parameters will be available as members of B.
   *
   * Sub-parsers allow you to compose parsers in a tree-like manner. You can add
   * A as a sub-parser of B via the ArgParser::parser member function, similar
   * to defining a parameter:
   * \code
    class B : public ArgParser
    {
      A a = parser<A>();
      //... additional parameters
    }
   * \endcode
   * A will automatically participate in the parsing of B, including tracking of
   * required parameters, dependencies between parameters, etc. No further
   * action by the user is required. You can access the parameters of A via the
   * member a of parser B, just like you would expect.
   *
   * The ArgParser::parser member function also allows you to set prefixes for a
   * parsers environment variables, property maps paths, and flags. For example,
   * assume that parser A has a parameter with a long flag of "--foo". If we
   * define a parser C via:
   * \code
    class C : public ArgParser
    {
      A a = parser<A>("A", "a", "a");
      //... additional parameters
    }
   * \endcode
   *
   * The the parameter associated with the "--foo" flag of A will be treated as
   * being defined with the flag "--a-foo". You can set separate prefixes for
   * environment variables, property map paths, and long flags.
   *
   * This prefix mechanism is especially useful when dealing with PropertyMaps.
   * It allows you to define a parser for each expected subsection of the
   * PropertyMap in isolation and then compose those parsers, setting the
   * subsection name as prefix for each parser.
   */
  class ArgParser
  {
  private:
    std::deque<std::shared_ptr<Intern::ParameterCore>> _parameters;
    std::deque<String> _errors;

    /// Path to program, i.e. argv[0]
    String _program{""};
    /// Description of this parser. Displayed in help text.
    String _description{""};

    /// Prefix for environment variables. Gets preprended to this parsers environment variables.
    String _env_prefix{""};
    /// Prefix for property map paths. Gets preprended to this parsers property map paths.
    String _property_prefix{""};
    /// Prefix for long flags. Gets preprended to this parsers long flags.
    String _flag_prefix{""};

    /// If true, unknown CLI flags do not produce errors
    bool _allow_unknown_cli_flags = false;

    static constexpr int priority_default = 0;
    static constexpr int priority_env = 1;
    static constexpr int priority_property = 2;
    static constexpr int priority_cli = 3;

  public:
    ArgParser() = default;

    /// Prefix constructor
    ArgParser(String env_prefix, String property_prefix, String flag_prefix) :
      _env_prefix(std::move(env_prefix)),
      _property_prefix(std::move(property_prefix)),
      _flag_prefix(std::move(flag_prefix))
    {
    }

    /**
     * \brief Allow parser to ignore unknown CLI flags
     *
     * The default behavior is to generate an error upon encountering a CLI flag
     * for which no parameter exists and declare the parse as failed. If this
     * option is set, then unknown flags will instead be silently ignored.
     *
     * \warning Erroring upon encountering an unknown flag is mean to protect
     * the user from errors due to misspelling a flag. Use this option with
     * care.
     */
    void allow_unknown_cli_flags(bool allow = true)
    {
      _allow_unknown_cli_flags = allow;
    }

    /// Error accessor
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
    [[nodiscard]] String help_text() const;

    /**
     * \brief Produce formatted output of all values in this parser
     *
     * \returns A string containing the values of all parameters
     *
     * \note This function shows the arguments given by the user.
     * Parser might change these values however they like.
     */
    [[nodiscard]] String display() const;

    /**
     * \brief Load PropertyMaps based on command line parameters or environment variables.
     *
     * \param[in,out] pmap Property map to merge additional property maps into.
     *
     * This method parses parameters marked as bootstrap parameters to load an
     * initial set of PropertyMaps for a full parse. Loaded property maps are merged
     * into the given PropertyMap. Entries are replaced by later maps.
     *
     * Errors that occur during the scanning of command line arguments and
     * during the parsing of bootstrap parameters are recorded. No further
     * validation is done.
     *
     * \note No parsing results are recorded into any parameters, including the
     * bootstrap parameters. Run ArgParser::parse() after bootstrapping your
     * PropertyMap to do the actual parsing.
     */
    [[nodiscard]] bool bootstrap(int argc, const char* const* argv, PropertyMap& pmap);

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
    [[nodiscard]] bool parse(int argc, const char* const* argv, const PropertyMap* property_map = nullptr);

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

    /**
     * \brief Adds the default values of all property-parameters to the given PropertyMap
     *
     * \param[in,out] pmap PropertyMap to fill
     *
     * This functions adds an entry to the given PropertyMap for each parameter with an associated PropertyMap path.
     * Each added entry contains that parameters default value.
     */
    void add_defaults_to_property_map(PropertyMap& pmap);

    /**
     * \brief Adds the values of all property-parameters to the given PropertyMap
     *
     * \param[in,out] pmap PropertyMap to fill
     *
     * This functions adds an entry to the given PropertyMap for each parameter with an associated PropertyMap path.
     * Each added entry contains that parameters parsed value.
     */
    void add_values_to_property_map(PropertyMap& pmap);

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
      std::shared_ptr<Intern::ParameterCore> core =
        std::make_shared<Intern::ParameterCore>(std::forward<T_>(default_value));

      // Store type-erased core for parsing purposes.
      // This allows the parser to treat all parameters the same,
      // irrespective of the types of their values.
      _parameters.push_back(core);

      // Return builder to user for defining the parameter fully
      return ParameterBuilder<T_>(std::move(core), _env_prefix, _property_prefix, _flag_prefix);
    }

    /**
     * \brief Add a sub-parser to this parser
     *
     * \param[in] env_prefix
     * All environment variables of the sub-parser will be prefixed with this string, separated by an underscore.
     *
     * \param[in] property_path_prefix
     * All PropertyMap paths of the sub-parser will be prefixed with this string, separated by a slash.
     *
     * \param[in] flag_prefix
     * All flags of the sub-parser will be prefixed with this string, separated by a dash.
     *
     * Prefixes are recursive. If B is a sub-parser of A with a flag prefix "b",
     * and C is a sub-parser of B with a flag prefix "b", then a flag "--foo" of
     * parser C will be treated as "--b-c-foo". Analogously for environment
     * variables and PropertyMap paths.
     *
     * The sub-parsers automatically take part in the parsing of the root
     * parser. No extra action by the user is required.
     */
    template<typename Parser_>
    Parser_ parser(String env_prefix = "", String property_path_prefix = "", String flag_prefix = "")
    {
      // NOTE: On the implementation of this function. Rather than parsing fully
      // recursively, and having to deal with tracking of used arguments and
      // error handling across parsers, we just add all parameters of the
      // sub-parser (with their prefix) to this parser. The parse method of this
      // parser will then parse the parameters the same way it parses its own
      // parameters. This might be terribly inconvenient later, but is very easy
      // to implement for now.
      const auto append = [](const String& a, const String& b, char delimiter) -> String
      {
        if(a.empty())
        {
          return b;
        }
        else if(!b.empty())
        {
          return a + delimiter + b;
        }
        else
        {
          return "";
        }
      };

      // Create sub-parser
      Parser_ parser(
        append(_env_prefix, env_prefix, '_'),
        append(_property_prefix, property_path_prefix, '/'),
        append(_flag_prefix, flag_prefix, '-'));

      // Add sub-parsers parameters to our own
      _parameters.insert(_parameters.end(), parser._parameters.begin(), parser._parameters.end());

      // Return parser to store it in current parser
      return parser;
    }

  private:
    /// Search for short flags, long flags, environment variables, or property paths that get used by multiple
    /// parameters
    bool _check_for_unique_parameters();

    void _collect_env();
    void _collect_property_map(const PropertyMap* property_map);
    void _collect_args(int argc, const char* const* argv);

    /// Validate parameters after parsing
    void _validate()
    {
      for(auto& core : _parameters)
      {
        core->validate(_errors);
      }
    }

    /// Derive parameters from other parameters
    void _derive()
    {
      for(auto& core : _parameters)
      {
        core->derive_self();
      }
    }

    /// Hide parameters if requested
    void _hide()
    {
      for(auto& core : _parameters)
      {
        core->hide();
      }
    }

    /// Create a PropertyMap from all property-parameters using the given extractor
    void _property_map(PropertyMap& pmap, const std::function<String(const Intern::ParameterCore&)>& extractor);

    /// Write flags, property path, environment variable and help text for given parameter to stringstream
    static void parameter_help(const Intern::ParameterCore& core, std::stringstream& s);

    static bool is_flag(const String& s)
    {
      return !s.empty() && s[0] == '-';
    }
  };
} // namespace FEAT
