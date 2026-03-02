// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/arg_parser.hpp>

#include <numeric>

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

  /// Reset this parameter to default
  void ParameterCore::reset()
  {
    value = default_value;
    priority = 0;
    arguments.clear();
    set_by_user = false;
  }
} // namespace FEAT
