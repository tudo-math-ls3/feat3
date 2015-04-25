#pragma once
#ifndef KERNEL_UTIL_SIMPLE_ARG_PARSER_HPP
#define KERNEL_UTIL_SIMPLE_ARG_PARSER_HPP 1

// includes, FEAST
#include <kernel/util/string.hpp>
#include <kernel/util/string_mapped.hpp>

// includes, system
#include <deque>
#include <map>
#include <utility>

namespace FEAST
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
   *   for the last option preceeding the parameter.
   * - Any command line argument preceeding the very first option is ignored by the parser.
   *
   * <b>Example:</b>\n
   * The command line call
        \verbatim
        my_app foobar --quiet --level 4 --precon spai jacobi
        \endverbatim
   * has a total of 8 arguments:
   * - 0: <c>"my_app"</c>: The name of the application's binary.
   * - 1: <c>"foobar"</c>: An argument ignored by the parser (as there is no option preceeding it).
   * - 2: <c>"--quiet"</c>: The first option named <c>quiet</c> without any parameters.
   * - 3: <c>"--level"</c>: The second option named <c>level</c> with one parameter following.
   * - 4: <c>"4"</c>: The first and only parameter for the preceeding option <c>level</c>.
   * - 5: <c>"--precon"</c>: The second option <c>precon</c> with two parameters following.
   * - 6: <c>"spai"</c>: The first parameter for the preceeding option <c>precon</c>.
   * - 7: <c>"jacobi"</c>: The second parameter for the preceeding option <c>precon</c>.
   *
   * Assume that the application's main function creates an SimpleArgParser object by
   * \code{.cpp}
     int main(int argc, char* argv[])
     {
       SimpleArgParser args(argc, argv);
       ...
     }
     \endcode
   * Then one may use the #check() member function to check whether an option was given or not.
   * - The call <code>args.check("solver")</code> will return -1, indicating the the option
   *   <c>solver</c> was not given in the command line.
   * - The call <code>args.check("quiet")</code> will return 0, indicating that the option
   *   <c>quiet</c> was given with no parameters.
   * - The call <code>args.check("level")</code> will return 1, incicating that the option
   *   <c>level</c> was given with one additional parameter.
   * - The call <code>args.check("precon")</code> will return 2, incicating that the option
   *   <c>precon</c> was given with two additional parameters.
   *
   * Furthermore, the #query() member function can be used to parse the parameters of a particular option.
   * For instance, the <c>level</c> parameter of the example above can be queried by
     \code{.cpp}
     int level(0);
     args.query("level", level);
     \endcode
   * The #query() function will return the value 1 indicating that the option <c>level</c> exists and that
   * 1 parameter was parsed successfully.
   *
   * For the above example, the <c>query</c> call in
     \code{.cpp}
     int quietness(0);
     args.query("quiet", quietness);
     \endcode
   * will return 0 indicating that the <c>quiet</c> option does not have any parseable parameters.
   *
   * Moreover, the <c>query</c> call in
     \code{.cpp}
     int precon(0);
     args.query("precon", precon);
     \endcode
   * will return -6, indicating that the 6th command line argument <c>"spai"</c> could not be parsed as
   * an <c>int</c>.
   *
   * Finally, the #query() function contains special handling for StringMapped parameters, i.e. in the
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
     args.query("precon", string_mapped(precon1, precon_map), string_mapped(precon2, precon_map));
     \endcode
   * The above <c>query</c> call will return the value 2 indicating that both parameters were parsed
   * successfully and it will hold that <c>precon1 = Precon::spai</c> and <c>precon2 = Precon::jacobi</c>.
   * \author Peter Zajac
   */
  class SimpleArgParser
  {
  private:
    typedef std::deque<std::pair<int,String>> IntStringDeque;
    /// number of skipped arguments
    int _num_skipped_args;
    /// command line arguments
    std::deque<String> _args;
    /// option-parameter map
    std::map<String, IntStringDeque> _opts;

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
      return int(it->second.size());
    }

    /**
     * \brief Queries an option and parses its parameters.
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
    int query(const String& option, Prms_&&... prms) const
    {
      // try to find the option
      auto it = _opts.find(option);
      if(it == _opts.end())
        return 0;

      // try to parse arguments
      return _query(it->second, std::size_t(0), std::forward<Prms_>(prms)...);
    }

  private:
    /// \cond internal
    void _process()
    {
      // get an option iterator
      auto opt = _opts.end();

      // we definately skip the first argument
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
          auto ik = _opts.insert(std::make_pair(optname, IntStringDeque()));
          opt = ik.first;

          // clear existing deque
          if(!ik.second)
            opt->second.clear();
        }
        else if(opt != _opts.end())
        {
          // push this argument as a parameter for the last processed option
          opt->second.push_back(std::make_pair(i, arg));
        }
        else
        {
          // We did not process any options yet, so we skip this argument
          _num_skipped_args = i + 1;
        }
      }
    }

    template<typename Prm_>
    int _parse(const IntStringDeque& isd, std::size_t offset, Prm_&& prm) const
    {
      // check whether we have another parameter
      if(offset >= isd.size())
        return 0;

      // try to parse the argument
      const auto& x = isd.at(offset);
      if(x.second.parse(std::forward<Prm_>(prm)))
        return +1; // success
      else
        return -x.first; // fail; return negative argument index
    }

    // overload for StringMapped
    template<typename ValueType_>
    int _parse(const IntStringDeque& isd, std::size_t offset, StringMapped<ValueType_>&& prm) const
    {
      // check whether we have another parameter
      if(offset >= isd.size())
        return 0;

      // try to parse the argument
      const auto& x = isd.at(offset);
      if(prm.lookup(x.second))
        return +1; // success
      else
        return -x.first; // fail; return negative argument index
    }

    template<typename Prm1_>
    int _query(const IntStringDeque& isd, std::size_t offset, Prm1_&& prm1) const
    {
      int rtn = _parse(isd, offset, std::forward<Prm1_>(prm1));
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
    int _query(const IntStringDeque& isd, std::size_t offset, Prm1_&& prm1, Prms_&&... prms) const
    {
      int rtn = _parse(isd, offset, std::forward<Prm1_>(prm1));
      if(rtn > 0)
      {
        // parse successful; continue parsing remaining parameters
        return _query(isd, ++offset, std::forward<Prms_>(prms)...);
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
} // namespace FEAST

#endif // KERNEL_UTIL_SIMPLE_ARG_PARSER_HPP
