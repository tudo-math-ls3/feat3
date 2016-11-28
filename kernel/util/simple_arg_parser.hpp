#pragma once
#ifndef KERNEL_UTIL_SIMPLE_ARG_PARSER_HPP
#define KERNEL_UTIL_SIMPLE_ARG_PARSER_HPP 1

// includes, FEAT
#include <kernel/util/string.hpp>
#include <kernel/util/string_mapped.hpp>

// includes, system
#include <deque>
#include <map>
#include <utility>

namespace FEAT
{
  /**
   * \brief Simple argument parser implementation
   *
   * This class implements a simple parser for command line arguments.
   *
   * The parser interprets the command line arguments according to the following rules:
   * - Any command line argument starting with <c>\--</c> (double hyphen) is interpreted as
   *   an \e option.
   * - Any command line argument not starting with <c>\--</c> is interpreted as a \e parameter
   *   for the last option preceding the parameter.
   * - Any command line argument preceding the very first option is ignored by the parser.
   *
   * <b>Example:</b>\n
   * The command line call
        \verbatim
        my_app foobar --quiet --level 4 --precon spai jacobi
        \endverbatim
   * has a total of 8 arguments:
   * - 0: <c>"my_app"</c>: The name of the application's binary.
   * - 1: <c>"foobar"</c>: An argument ignored by the parser (as there is no option preceding it).
   * - 2: <c>"--quiet"</c>: The first option named <c>quiet</c> without any parameters.
   * - 3: <c>"--level"</c>: The second option named <c>level</c> with one parameter following.
   * - 4: <c>"4"</c>: The first and only parameter for the preceding option <c>level</c>.
   * - 5: <c>"--precon"</c>: The second option <c>precon</c> with two parameters following.
   * - 6: <c>"spai"</c>: The first parameter for the preceding option <c>precon</c>.
   * - 7: <c>"jacobi"</c>: The second parameter for the preceding option <c>precon</c>.
   *
   * Assume that the application's main function creates an SimpleArgParser object by
   * \code{.cpp}
     int main(int argc, char* argv[])
     {
       SimpleArgParser args(argc, argv);
       ...
     }
     \endcode
   *
   * The first step, which is highly recommended, is to add all option names, which are supported by
   * your application, to the argument parser by using the #support() function:
   * \code{.cpp}
       args.support("level");
       args.support("precon");
       args.support("quiet");
     \endcode
   * Afterwards, you can call the #query_unsupported() function to get a list of all options, which were
   * given in the command line arguments but were not marked as supported. Note that for the parser, it is
   * not an error if unsupported options were given - you have to take care of the error handling by
   * yourself. An appropriate error handling code could look like this:
   * \code{.cpp}
     std::deque<std::pair<int,String> > unsupported = args.query_unsupported();
     if( !unsupported.empty() )
     {
       // print all unsupported options to cerr
       for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
         std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;

       // abort program execution here, if you want to
     }
     \endcode
   *
   * Moreover, one may use the #check() member function to check whether an option was given or not.
   * - The call <code>args.check("solver")</code> will return -1, indicating the the option
   *   <c>solver</c> was not given in the command line.
   * - The call <code>args.check("quiet")</code> will return 0, indicating that the option
   *   <c>quiet</c> was given with no parameters.
   * - The call <code>args.check("level")</code> will return 1, incicating that the option
   *   <c>level</c> was given with one additional parameter.
   * - The call <code>args.check("precon")</code> will return 2, incicating that the option
   *   <c>precon</c> was given with two additional parameters.
   *
   * Furthermore, the #parse() member function can be used to parse the parameters of a particular option.
   * For instance, the <c>level</c> parameter of the example above can be queried by
     \code{.cpp}
     int level(0);
     args.parse("level", level);
     \endcode
   * The #parse() function will return the value 1 indicating that the option <c>level</c> exists and that
   * 1 parameter was parsed successfully.
   *
   * For the above example, the <c>parse</c> call in
     \code{.cpp}
     int quietness(0);
     args.parse("quiet", quietness);
     \endcode
   * will return 0 indicating that the <c>quiet</c> option does not have any parseable parameters.
   *
   * Moreover, the <c>parse</c> call in
     \code{.cpp}
     int precon(0);
     args.parse("precon", precon);
     \endcode
   * will return -6, indicating that the 6th command line argument <c>"spai"</c> could not be parsed as
   * an <c>int</c>.
   *
   * Furthermore, the #parse() function contains special handling for StringMapped parameters, i.e. in the
   * example above, one may use the following code to parse the <c>precon</c> parameters to corresponding
   * enumeration values.
     \code{.cpp}
     enum class Precon
     {
       none,
       spai,
       jacobi
     };
     std::map<String, Precon> precon_map;
     precon_map.insert(std::make_pair("none", Precon::none));
     precon_map.insert(std::make_pair("spai", Precon::spai));
     precon_map.insert(std::make_pair("jacobi", Precon::jacobi));

     Precon precon1(Precon::none), precon2(Precon::none);
     args.parse("precon", string_mapped(precon1, precon_map), string_mapped(precon2, precon_map));
     \endcode
   * The above <c>parse</c> call will return the value 2 indicating that both parameters were parsed
   * successfully and it will hold that <c>precon1 = Precon::spai</c> and <c>precon2 = Precon::jacobi</c>.
   *
   * Finally, this class contains a member function named #query(), which retrieves the argument index
   * a particular option and the deque of all parameters for this option. This member function can be used
   * for parsing custom and more complex parameter sets, which can not be covered by the #parse() function.
   *
   * \author Peter Zajac
   */
  class SimpleArgParser
  {
  private:
    /// number of skipped arguments
    int _num_skipped_args;
    /// command line arguments
    std::deque<String> _args;
    /// option-parameter map
    std::map<String, std::pair<int, std::deque<String> > > _opts;
    /// supported option set
    std::map<String, String> _supported;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] argc
     * The number of arguments passed to the applicaton.
     *
     * \param[in] argv
     * The array of all command line arguments passed to the application.
     */
    explicit SimpleArgParser(int argc, const char* const * const argv) :
      _num_skipped_args(0),
      _args(),
      _opts()
    {
      for(int i(0); i < argc; ++i)
        _args.push_back(String(argv[i]));
      _process();
    }

    /**
     * \brief Returns the total number of arguments passed in the constructor.
     *
     * \returns The number of arguments passed in the constructor, including the command name (first argrument).
     */
    int num_args() const
    {
      return int(_args.size());
    }

    /**
     * \brief Returns the number of arguments skipped by the parser.
     *
     * \returns The number of skipped arguments.
     */
    int num_skipped_args() const
    {
      return _num_skipped_args;
    }

    /**
     * \brief Returns an argument.
     *
     * \param[in] i
     * The index of the argument that is to be returned.
     *
     * \returns
     * The \p i-th argument passed to the constructor or an empty string if \p i is out-of-bounds.
     */
    String get_arg(int i) const
    {
      if((i < 0) || (i >= int(_args.size())))
        return String();
      return _args.at(std::size_t(i));
    }

    /**
     * \brief Returns a deque of all arguments skipped by the parser.
     */
    std::deque<String> skipped_args() const
    {
      std::deque<String> dsa;
      for(std::size_t i(0); i < std::size_t(_num_skipped_args); ++i)
        dsa.push_back(_args.at(i));
      return dsa;
    }

    /**
     * \brief Adds an option to the set of supported options.
     *
     * This function adds an option name to the set of supported options, which
     * is used by the #query_unsupported() function.
     *
     * \param[in] option
     * The name of the option that is to be added.
     *
     * \param[in] description
     * A descriptive string that describes the meaning of the option and its parameters.
     * This is used for the generation of a commented option list; see the documentation
     * of the #get_supported_help() function for details and remarks about this help string.
     *
     * \note
     * This function has no effect on the #check(), #query() and #parse() functions, i.e.
     * you can use these functions to query options independently of whether they have
     * been added to the set of supported options or not.
     */
    void support(const String& option, const String& description = String())
    {
      _supported.emplace(option, description);
    }

    /**
     * \brief Returns a string of all supported options and their help strings.
     *
     * This function returns a new-line separated string of all options and their
     * descriptions, which have been added to the parser using the #support() function.
     *
     * Each option is listed in the form
       \verbatim
       --option description
       \endverbatim
     *
     * Recommendations:
     * - If your option has no parameters, your description string should begin with
     *   a new-line character, so that the description starts in the next line.
     * - Parameters should be enclosed in angle-brackets and their names should be
     *   referenced in the description, e.g. choosing <c>option="level"</c> and
     *   <c>description="\<n\>\nSets the refinement level to \<n\>\n"</c>
     *   leads to the output
         \verbatim
         --level <n>
         Sets the refinement level to <n>.
         \endverbatim
     */
    String get_supported_help() const
    {
      String s;
      for(auto it = _supported.begin(); it != _supported.end(); ++it)
      {
        s += "--" + it->first + " " + it->second + "\n";
      }
      return s;
    }

    /**
     * \brief Queries all unsupported options passed in the command line.
     *
     * This function loops over all options, which were given in the command line, and
     * checks whether they have been added to the list of supported options by using
     * the #support() function.
     *
     * \returns
     * A deque of int-String pairs representing all unsupported options.
     * The first (int) component is the index of the argument in the command line,
     * whereas the second (String) component is the name of the unsupported option.
     */
    std::deque<std::pair<int,String> > query_unsupported() const
    {
      std::deque<std::pair<int,String> > result;

      // loop over all given options
      for(auto it = _opts.begin(); it != _opts.end(); ++it)
      {
        // get the option name and its index
        String name = it->first;
        int idx = it->second.first;

        // check whether we support that option
        if(_supported.find(name) == _supported.end())
        {
          // unsupported option; so push it
          result.push_back(std::make_pair(idx, name));
        }
      }

      // return our result
      return result;
    }

    /**
     * \brief Checks whether an option was given.
     *
     * \param[in] option
     * The name of the option (without the leading <c>\--</c>) that is to be checked for.
     *
     * \returns
     * A value \e n >= 0 specifying the number of parameters of the option if the
     * option was given, or -1 if the option was not given at all.
     */
    int check(const String& option) const
    {
      auto it =_opts.find(option);
      if(it == _opts.end())
        return -1;
      return int(it->second.second.size());
    }

    /**
     * \brief Query the parameters of an option.
     *
     * This function checks whether an option was given and returns a pointer to a pair
     * whose first member is the index of the command line argument identifying this
     * option and the second member is a deque of Strings representing the parameters
     * of this option.
     *
     * \param[in] option
     * The name of the option (without the leading <c>\--</c>) that is to be checked for.
     *
     * \returns
     * A const pointer to the pair of argument index and parameter deque or \c nullptr,
     * if the option was not given.
     */
    const std::pair<int, std::deque<String> >* query(const String& option) const
    {
      auto it = _opts.find(option);
      if(it == _opts.end())
        return nullptr;
      return &(it->second);
    }

    /**
     * \brief Parses the parameters of an option.
     *
     * This function uses the String::parse() member function to parse the parameter values. This function
     * also contains special handling for StringMapped objects, which allows to parse enumeration values
     * directly.
     *
     * \param[in] option
     * The name of the option (without the leading <c>\--</c>) that is to be queried.
     *
     * \param[in,out] prms
     * A list of references to the objects that shall be parsed as parameters.
     *
     * \returns
     * - 0, if the option was not given at all or if it was given, but without any parameters.
     * - \e n > 0, if the option was given and \e n parameters were parsed successfully.
     * - \e n < 0, if the option was given, but the \e n-th command line argument could not be parsed
     *   as a parameter. The faulty argument can be accessed by <c>get_arg(n)</c> for error handling.
     */
    template<typename... Prms_>
    int parse(const String& option, Prms_&&... prms) const
    {
      // try to find the option
      auto it = _opts.find(option);
      if(it == _opts.end())
        return 0;

      // try to parse arguments
      return _parse(it->second, std::size_t(0), std::forward<Prms_>(prms)...);
    }

  private:
    /// \cond internal
    void _process()
    {
      // get an option iterator
      auto opt = _opts.end();

      // we definitely skip the first argument
      _num_skipped_args = 1;

      // loop over all arguments
      for(int i(1); i < int(_args.size()); ++i)
      {
        // get the trimmed argument
        String arg = _args.at(std::size_t(i)).trim();

        // does it start with '--' ?
        if((arg.size() > std::size_t(1)) && (arg[0] == '-') && (arg[1] == '-'))
        {
          // it's an option
          String optname = arg.substr(2).trim();

          // let's insert it
          auto ik = _opts.insert(std::make_pair(optname, std::make_pair(i, std::deque<String>())));
          opt = ik.first;

          // clear existing deque
          if(!ik.second)
          {
            // option was already given, so overwrite it
            opt->second.first = i;
            opt->second.second.clear();
          }
        }
        else if(opt != _opts.end())
        {
          // push this argument as a parameter for the last processed option
          opt->second.second.push_back(arg);
        }
        else
        {
          // We did not process any options yet, so we skip this argument
          _num_skipped_args = i + 1;
        }
      }
    }

    template<typename Prm_>
    int _parse_prm(const std::pair<int,std::deque<String> >& isd, std::size_t offset, Prm_&& prm) const
    {
      // check whether we have another parameter
      if(offset >= isd.second.size())
        return 0;

      // try to parse the argument
      if(isd.second.at(offset).parse(std::forward<Prm_>(prm)))
        return +1; // success
      else
        return -isd.first - int(offset) - 1; // fail; return negative argument index
    }

    // overload for StringMapped
    template<typename ValueType_>
    int _parse_prm(const std::pair<int,std::deque<String> >& isd, std::size_t offset, StringMapped<ValueType_>&& prm) const
    {
      // check whether we have another parameter
      if(offset >= isd.second.size())
        return 0;

      // try to parse the argument
      if(prm.lookup(isd.second.at(offset)))
        return +1; // success
      else
        return -isd.first - int(offset) - 1; // fail; return negative argument index
    }

    template<typename Prm1_>
    int _parse(const std::pair<int,std::deque<String> >& isd, std::size_t offset, Prm1_&& prm1) const
    {
      int rtn = _parse_prm(isd, offset, std::forward<Prm1_>(prm1));
      if(rtn > 0)
      {
        // parse successful; no more arguments to be parsed
        return int(offset) + 1;
      }
      else if(rtn == 0)
      {
        // parameter missing
        return int(offset);
      }
      else
      {
        // failed to parse
        return rtn;
      }
    }

    template<typename Prm1_, typename... Prms_>
    int _parse(const std::pair<int,std::deque<String> >& isd, std::size_t offset, Prm1_&& prm1, Prms_&&... prms) const
    {
      int rtn = _parse_prm(isd, offset, std::forward<Prm1_>(prm1));
      if(rtn > 0)
      {
        // parse successful; continue parsing remaining parameters
        return _parse(isd, ++offset, std::forward<Prms_>(prms)...);
      }
      else if(rtn == 0)
      {
        // no more parameters, quit here
        return int(offset);
      }
      else
      {
        // failed to parse
        return rtn;
      }
    }
    /// \endcond
  }; // class SimpleArgParser
} // namespace FEAT

#endif // KERNEL_UTIL_SIMPLE_ARG_PARSER_HPP
