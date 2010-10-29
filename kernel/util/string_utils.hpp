#pragma once
#ifndef STRING_UTILS_HPP
#define STRING_UTILS_HPP

// includes, system
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <cassert>

namespace FEAST
{

  /**
  * \brief collection of various string utilities
  *
  * \author Dirk Ribbrock
  * \author Dominik Goeddeke
  */
  class StringUtils
  {
    public:
    /* *****************
    * member functions *
    *******************/

    /**
    * \brief converting an item to a string
    *
    * \param[in] item
    * The item to stringify
    */
    template <typename T_>
    static inline std::string stringify(const T_ & item)
    {
      std::ostringstream s;
      s << item;
      return s.str();
    }


    /**
    * \brief converting an item to a string (overload for std::string)
    *
    * \param[in] item
    * The item to stringify
    */
    static inline std::string stringify(const std::string & item)
    {
      return item;
    }


    /**
    * \brief converting an item to a string (overload for char)
    *
    * \param[in] item
    * The item to stringify
    */
    static inline std::string stringify(const char & item)
    {
      return std::string(1, item);
    }


    /**
    * \brief converting an item to a string (overload for unsigned char)
    *
    * \param[in] item
    * The item to stringify
    */
    static inline std::string stringify(const unsigned char & item)
    {
      return std::string(1, item);
    }


    /**
    * \brief converting an item to a string (overload for bool)
    *
    * \param[in] item
    * The item to stringify
    */
    static inline std::string stringify(const bool & item)
    {
      return item ? "true" : "false";
    }


    /**
    * \brief converting an item to a string (overload for char *, which isn't a screwup like other pointers)
    *
    * \param[in] item
    * The item to stringify
    */
    static inline std::string stringify(const char * const item)
    {
      return std::string(item);
    }


    /**
    * \brief Concatenate strings between begin and end iterator
    *
    * \todo also implement for std::vector iterator
    * COMMENT_HILMAR: is there a clever way to provide this functionality for all STL containers at once?
    *
    * \author Dirk Ribbrock
    */
    static std::string join(
      std::list<std::string>::const_iterator begin,
      std::list<std::string>::const_iterator end,
      const std::string& delimiter)
    {
      std::string result;

      if (begin != end)
      {
        while (true)
        {
          result += *begin;
          if (++begin == end)
          {
            break;
          }
          result += delimiter;
        }
      }
      result +="\n";
      return result;
    }
  }; // class StringUtils



  /**
  * \brief class for semi-automatically pretty printing contiguous blocks of strings
  *
  * The class can be used to produce pretty printed output blocks like the following:
  *
  * \verbatim
      PRE #######################################
      PRE #          SOLVER STATISTICS          #
      PRE #######################################
      PRE # it: 42                              #
      PRE # conv rate: 0.23                     #
      PRE # MFLOP/sec: 666                      #
      PRE #######################################
      PRE # some additional information whose length can not be estimated a priori
      PRE # and which should thus not use a closing '#' at the end of the line
      PRE #######################################
  * \endverbatim
  *
  * The width of the block, the delimiter symbol, and the prefix are set by the user. The strings to be written can be
  * centered or left aligned. For too long strings, the right delimiter can be omitted (is omitted automatically, resp.).
  *
  * Some construction rules: A blank is appended to the prefix. One blank is inserted between left delimiter and the
  * string. Between the string and the right delimiter there must be at least one blank, otherwise the right delimiter is
  * omitted.
  *
  * \todo add functionalities for producing tabular output like this:
  * \verbatim
      #########################################
      #   it |     conv |    MFLOP |     MFR  #
      # ------------------------------------- #
      #    1 |     0.23 |     42.1 |     666  #
      #    2 |     0.22 |     42.0 |     667  #
      #    3 |     0.21 |     42.0 |     668  #
      #    4 |     0.24 |     42.1 |     665  #
      #########################################
  * \endverbatim
  * \todo maybe use two different chars for separator lines and delimiter
  * \author Hilmar Wobker
  */
  class PrettyPrinter
  {

  private:

    /// string holding the pretty printed block of messages
    std::string _block;

    ///width of the pretty printed block (from left to right delimiter, excluding the prefix)
    unsigned int _width;

    /// character to be used to enclose lines and to build separator lines
    char _delim;

    /// string each line is prefixed with (useful in connection with grep commands)
    std::string _prefix;


  public:

    /**
    * \brief constructor using prefix
    *
    * \param[in] width
    * width of the lines in the pretty printed block
    *
    * \param[in] delim
    * character to be used to enclose lines and to build separator lines
    *
    * \param[in] prefix
    * prefix which can be useful in connection with grep commands
    */
    PrettyPrinter(
      unsigned int width,
      char delim,
      const std::string& prefix)
      : _block(""),
        _width(width),
        _delim(delim),
        _prefix(prefix)
    {
    }


    /**
    * \brief constructor using an empty prefix
    *
    * \param[in] width
    * width of the lines in the pretty printed block
    *
    * \param[in] delim
    * character to be used to enclose lines and to build separator lines
    */
    PrettyPrinter(
      unsigned int width,
      char delim)
      : _block(""),
        _width(width),
        _delim(delim),
        _prefix("")
    {
    }


    /// destructor
    ~PrettyPrinter()
    {
    }


    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the pretty printed block
    *
    * Return a reference in order to avoid making a copy. Make this return reference constant so that the user
    * cannot change the string.
    *
    * \return reference to the pretty printed block #_block
    */
    inline const std::string& block() const
    {
      return _block;
    }

    /**
    * \brief setter for the pretty printed block
    */
    inline void set_block(std::string block)
    {
      _block = block;
    }

    /**
    * \brief getter for the width of the pretty printed block
    *
    * \return width of the pretty printed block #_width
    */
    inline unsigned int width() const
    {
      return _width;
    }

    /**
    * \brief setter for the width of the pretty printed block
    */
    inline void set_width(unsigned int width)
    {
      _width = width;
    }

    /**
    * \brief getter for the prefix
    *
    * Return a reference in order to avoid making a copy. Make this return reference constant so that the user
    * cannot change the string.
    *
    * \return reference to the prefix #_prefix
    */
    inline const std::string& prefix() const
    {
      return _prefix;
    }

    /**
    * \brief setter for the prefix
    */
    inline void set_prefix(std::string prefix)
    {
      _prefix = prefix;
    }

    /**
    * \brief getter for the delimiter character
    *
    * \return delimiter #_delim
    */
    inline char delim() const
    {
      return _delim;
    }

    /**
    * \brief setter for the delimiter character
    */
    inline void set_delim(char delim)
    {
      _delim = delim;
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief pretty printing a separator line of length #_width consisting of the delimiter char #_delim
    *
    * This function can be used to automatically produce the 1., 3., 7. and last line of the example output in the class
    * description.
    */
    void add_line_sep()
    {
      _block += _prefix + " " + std::string(_width, _delim) + "\n";
    }


    /**
    * \brief adding a line of length #_width, with the given string \a s centered and #_delim as left and right delimiter
    *
    * This function can be used to automatically produce the 2. line of the example output in the class description.
    * If the number of blanks is odd, then on the left side one blank less is used. If the string is too long, it is
    * displayed left aligned without right delimiter.
    *
    * \param[in] s
    * the string to be centered in the line
    */
    void add_line_centered(const std::string& s)
    {
      unsigned int num_blanks_total((_width - s.size() - 2));
      if(num_blanks_total < 2)
      {
        // if the string is too long, produce a left aligned line without right delimiter
        add_line_no_right_delim(s);
      }
      else
      {
        // calculate blanks on the left and the right side
        unsigned int num_blanks_left(num_blanks_total / 2);
        unsigned int num_blanks_right(num_blanks_total%2 == 0 ? num_blanks_left : num_blanks_left + 1);
        // add the line to the pretty printed block
        _block +=   _prefix + " " + StringUtils::stringify(_delim) + std::string(num_blanks_left, ' ') + s
                  + std::string(num_blanks_right, ' ') + StringUtils::stringify(_delim) + "\n";
      }
    }


    /**
    * \brief adding a line of length #_width, with the given string \a s left aligned and #_delim as left and right
    *        delimiter
    *
    * This function can be used to automatically produce the 4., 5. and 6. line of the example output in the class
    * description. If the string is too long, it is displayed left aligned without right delimiter.
    *
    * \param[in] s
    * the string to be written
    */
    void add_line(const std::string& s)
    {
      int num_blanks((_width - s.size() - 2));
      if(num_blanks < 2)
      {
        // if the string is too long, produce a left aligned line without right delimiter
        add_line_no_right_delim(s);
      }
      else
      {
        // add the line to the pretty printed block
        _block +=   _prefix + " " + StringUtils::stringify(_delim) + " " + s + std::string(num_blanks - 1, ' ')
                  + StringUtils::stringify(_delim) + "\n";
      }
    }

    /**
    * \brief adding a line of length #_width, with the given string \a s left aligned and #_delim as left delimiter
    *
    * This function can be used to automatically produce the 8. and 9. line of the example output in the class
    * description.
    *
    * \param[in] s
    * the string to be written
    */
    void add_line_no_right_delim(const std::string& s)
    {
      _block += _prefix + " " + StringUtils::stringify(_delim) + " " + s + "\n";
    }

    /**
    * \brief reset the pretty printed block to an emtpy string (equivalent to PrettyPrinter::set_block(""))
    */
    inline void reset_block()
    {
      _block = "";
    }

    /**
    * \brief print the pretty printed block to the given stream
    */
    void print(std::ostream& stream)
    {
      stream << _block;
    }

  };



  /**
  * \brief class for storing and maintainig a prefix string which can be used for logging
  *
  * The prefix string can be created incrementally by "pushing" new substrings to it (see Prefix::push) and "popping"
  * them again (see Prefix::pop). This is, e.g., useful in recursive subroutine calls, e.g.
  *
  *   Prefix prefix = new Prefix();
  *   prefix.push("Proc" + stringify(Process::rank) + ":");
  *     // --> example screen output:
  *     $ Proc42: some log output
  *     $ Proc42: starting solver...
  *
  *   prefix.push("BiCG:")
  *     $ Proc42:BiCG: some log output
  *     $ Proc42:BiCG: some log output
  *     $ Proc42:BiCG: start MG preconditioner
  *
  *   prefix.push("MG:")
  *     $ Proc42:BiCG:MG: some log output
  *     $ Proc42:BiCG:MG: start smoother
  *
  *   prefix.push("Jacobi:")
  *     $ Proc42:BiCG:MG:Jacobi: some log output
  *     $ Proc42:BiCG:MG:Jacobi: smoother finished
  *
  *   prefix.pop();
  *     $ Proc42:BiCG:MG: some log output
  *     $ Proc42:BiCG:MG: MG finished
  *
  *   prefix.pop();
  *     $ Proc42:BiCG: some log output
  *     $ Proc42:BiCG: BiCG finished
  *
  *   prefix.pop();
  *     $ Proc42: solving finished
  *
  * \todo This class has not been tested yet!
  *
  * COMMENT_HILMAR:
  * How should we use this feature? Should there be some sort of process-global prefix which is known to all
  * objects and can be modified from everywhere? (If yes: where should it reside?) Or is the prefix object
  * passed through all the routines? (Which means adapting all the interfaces...)
  */
  class Prefix
  {
  private:
    /// string holding the prefix built by pushing and popping substrings
    std::string _s;

    /// start positions of the substrings inside the string #_s
    std::vector<unsigned int> _start_pos;

  public:
    /// constructor
    Prefix()
      : _s("")
    {
    }

    /// destructor
    ~Prefix()
    {
      _start_pos.clear();
    }

    /**
    * \brief appends a string to the prefix string
    *
    * This function receives a string and appends it to the string #_s.
    *
    * \param[in] string_to_append
    * string to be appended to the string #_s
    *
    * \author Hilmar Wobker
    */
    void push(std::string string_to_append)
    {
      // forbid pushing empty strings
      assert(string_to_append.size() > 0);
      // store the length of the current prefix string (which is at the same time the start position of the substring s
      // within the prefix string)
      _start_pos.push_back(_s.size());
      _s += string_to_append;
    }

    /**
    * \brief removes from the prefix string the substring appended last
    *
    * This function removes from the prefix string #_s the substring which has been appended last.
    *
    * \author Hilmar Wobker
    */
    void pop()
    {
      // assert that there is a substring to be removed
      assert(_start_pos.size() > 0);
      // replace the prefix string by its leading part resulting from removing the substring appended last
      _s = _s.substr(0, _start_pos.back());
      _start_pos.pop_back();
    }

    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the prefix string
    *
    * Return a reference in order to avoid making a copy. Make this return reference constant so that the user
    * cannot change the string.
    *
    * \return reference to prefix #_s
    */
    inline const std::string& s() const
    {
      return _s;
    }
  };

} // namespace FEAST

#endif //  #ifndef STRING_UTILS_HPP
