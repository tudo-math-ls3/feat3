// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/arg_parser.hpp>

#include <numeric>
#include <sstream>
#include <set>

namespace FEAT::Intern
{
  /**
   * \brief Tries to create a user-readable name for this parameter
   *
   * Tries to derive a name from the long flag, the short flag, the property path or the environment variable.
   * If all four are empty the name is '<anonymous>'.
   */
  String ParameterCore::name() const
  {
    if(!custom_name.empty())
    {
      return custom_name;
    }
    if(!long_flag.empty())
    {
      return long_flag.substr(2);
    }

    if(!short_flag.empty())
    {
      return short_flag.substr(1);
    }

    if(!property_path.empty())
    {
      return property_path;
    }

    if(!environment_variable.empty())
    {
      return environment_variable;
    }

    return "<anonymous>";
  }

  /// Add an argument to this parameter
  void ParameterCore::add_argument(const String& s, int incoming_priority)
  {
    if(priority < incoming_priority)
    {
      // New argument is of higher priority than all previous arguments
      // Clear the list and set new priority
      arguments.clear();
      arguments.push_back(s.trim());
      priority = incoming_priority;
    }
    else if(priority == incoming_priority)
    {
      // New argument is of same priority
      // Store it for now
      arguments.push_back(s.trim());
    }
  }

  void ParameterCore::validate(std::deque<String>& errors) const
  {
    if(required && !set_by_user)
    {
      errors.push_back("Error: Parameter " + name() + " not given, but is required!");
      return;
    }

    if(validator)
    {
      try
      {
        validator(value);
      }
      catch(std::bad_any_cast& e)
      {
        XABORTM(e.what());
      }
      catch(std::invalid_argument& e)
      {
        errors.push_back("Error: Validation of parameter " + name() + " failed with message: " + e.what());
      }
    }

    for(const auto& [core, condition] : needs)
    {
      if(core.expired())
      {
        XABORTM("Expired weak ptr to parameter core of " + name());
      }

      auto ptr = core.lock();

      if(set_by_user && condition(value))
      {
        if(!ptr->set_by_user)
        {
          errors.push_back(
            "Error: Parameter " + name() + " needs parameter " + ptr->name() + ", but parameter " + ptr->name() +
            " was not set by user");
        }
      }
    }
  }

  void ParameterCore::parse(std::deque<String>& errors)
  {
    if(arguments.empty())
    {
      // Nothing to parse
      return;
    }

    set_by_user = true;

    // Merge arguments of collection options
    // NOTE: This is slightly inefficient, because
    // the collection parser is just going to split the string
    // by whitespace anyway. But CLI arguments come pre-split
    // and this way we can handle all argument sources with the same parser
    if(type == ParameterType::collection_option)
    {
      String merged = std::accumulate(
        arguments.begin(),
        arguments.end(),
        String(""),
        [](const String& a, const String& b) { return a + " " + b; });

      // Replace arguments with merged version
      arguments.clear();
      arguments.push_back(merged.trim());
    }

    if(arguments.size() > 1)
    {
      errors.push_back("Error: Parameter " + name() + " was set multiple times. Only the first value is considered!");
    }

    if(!parser)
    {
      errors.push_back(
        "Error: No parser for parameter " + name() +
        ". Ensure operator>> is available for your type or set a custom parser!");
      return;
    }

    try
    {
      parser(value, arguments.front());
    }
    catch(std::bad_any_cast& e)
    {
      XABORTM(e.what());
    }
    catch(std::invalid_argument& e)
    {
      errors.push_back("Error: Parsing for parameter " + name() + " failed with message: " + e.what());
    }
  }

  void ParameterCore::derive_self()
  {
    if(!set_by_user && deriver)
    {
      // Build deque of values
      std::deque<std::reference_wrapper<std::any>> args;

      for(const auto& core : derives_from)
      {
        if(core.expired())
        {
          XABORTM("Expired weak ptr to parameter core of " + name());
        }

        auto ptr = core.lock();

        args.emplace_back(ptr->value);
      }

      // Call deriver
      deriver(value, args);
    }
  }

  void ParameterCore::hide()
  {
    for(auto& [core, condition] : hide_if)
    {
      if(condition(value))
      {
        if(core.expired())
        {
          XABORTM("Expired weak ptr to parameter core of " + name());
        }

        auto ptr = core.lock();

        if(condition(value))
        {
          ptr->hidden = true;
        }
      }
    }
  }

  /// Reset this parameter to default
  void ParameterCore::reset()
  {
    value = default_value;
    priority = 0;
    arguments.clear();
    set_by_user = false;
  }
} // namespace FEAT::Intern

namespace FEAT
{
  [[nodiscard]] String ArgParser::help_text() const
  {
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

  [[nodiscard]] String ArgParser::display() const
  {
    std::stringstream stream;

    // Name padding is size of longest parameter name plus 3 chars for aesthetics
    String::size_type name_padding(0);
    for(const auto& core : _parameters)
    {
      name_padding = std::max(name_padding, core->name().size());
    }
    name_padding += 3;

    const String::size_type source_padding(10);
    const String::size_type header_width(75);

    // Print header
    stream << String("Parameter").pad_back(name_padding, ' ');
    stream << "|";
    stream << String("Source").pad_back(source_padding, ' ');
    stream << "|";
    stream << "Value\n";
    stream << String("").pad_back(header_width, '-') << "\n";

    for(const auto& core : _parameters)
    {
      // Skip hidden parameters
      if(core->hidden)
      {
        continue;
      }

      XASSERT(core->value.has_value());

      stream << core->name().pad_back(name_padding, '.');
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

  [[nodiscard]] bool ArgParser::bootstrap(int argc, const char* const* argv, PropertyMap& pmap)
  {
    // Reset errors
    _errors.clear();

    // Set default values, reset priorities, clear arguments
    for(auto& core : _parameters)
    {
      core->reset();
    }

    // Collect values in environment variables
    _collect_env();

    // Collect command line args
    _collect_args(argc, argv);

    // Parse bootstrap parameters and load property maps
    PropertyMap result;
    for(auto& core : _parameters)
    {
      if(core->bootstrap)
      {
        core->parse(_errors);
        // SAFETY: Cast to std::deque<String> is ok, because only Parameter<std::deque<String>> may be marked as bootstrap
        for(const String& path : std::any_cast<std::deque<String>>(core->value))
        {
          PropertyMap tmp_map;
          tmp_map.read(path);
          pmap.merge(tmp_map, true);
        }
      }
    }

    // Reset again, to discourage users from using this as an alternative to the parse function
    for(auto& core : _parameters)
    {
      core->reset();
    }

    return _errors.empty();
  }

  [[nodiscard]] bool ArgParser::parse(int argc, const char* const* argv, const PropertyMap* property_map)
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

    // Run validators
    _validate();

    // Derive values
    _derive();

    // Hide parameters
    _hide();

    // Call post_parse hook on successful parse
    if(_errors.empty())
    {
      post_parse();
    }

    // Sort parameters alphabetically for help and display output
    std::sort(
      _parameters.begin(),
      _parameters.end(),
      [](auto& a, auto& b) { return a->name().compare_no_case(b->name()) < 0; });

    return _errors.empty();
  }

  void ArgParser::add_defaults_to_property_map(PropertyMap& pmap)
  {
    _property_map(
      pmap,
      [](const Intern::ParameterCore& core)
      {
        std::stringstream ss;
        core.printer(ss, core.default_value);
        return ss.str();
      });
  }

  void ArgParser::add_values_to_property_map(PropertyMap& pmap)
  {
    _property_map(
      pmap,
      [](const Intern::ParameterCore& core)
      {
        std::stringstream ss;
        core.printer(ss, core.value);
        return ss.str();
      });
  }

  bool ArgParser::_check_for_unique_parameters()
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

  void ArgParser::_collect_env()
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

  void ArgParser::_collect_property_map(const PropertyMap* property_map)
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

  void ArgParser::_collect_args(int argc, const char* const* argv)
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
          if(std::any_cast<const bool&>(core->default_value))
          {
            core->add_argument("false", priority_cli);
          }
          else
          {
            core->add_argument("true", priority_cli);
          }
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
        if(!_allow_unknown_cli_flags)
        {
          // We found no parameter for this flag. Generate error
          _errors.push_back("Error: Unknown flag " + flag);
        }

        // Then advance to next flag
        while(idx < args.size() && !is_flag(args[idx]))
        {
          idx++;
        }
      }
    }
  }

  void ArgParser::_property_map(PropertyMap& pmap, const std::function<String(const Intern::ParameterCore&)>& extractor)
  {
    for(auto& core : _parameters)
    {
      if(core->property_path.empty())
      {
        // Not a property-parameter. Skip it.
        continue;
      }

      const auto path_components = core->property_path.split_by_charset("/");
      PropertyMap* section = &pmap;
      // Iterate over sections. Last entry of path_components is entry.
      for(Index i(0); i < path_components.size() - 1; i++)
      {
        PropertyMap* subsection = section->get_sub_section(path_components[i]);
        if(subsection == nullptr)
        {
          section = section->add_section(path_components[i]);
        }
        else
        {
          section = subsection;
        }
      }

      // Write entry
      section->add_entry(path_components.back(), extractor(*core));
    }
  }

   void ArgParser::parameter_help(const Intern::ParameterCore& core, std::stringstream& s)
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
} // namespace FEAT
