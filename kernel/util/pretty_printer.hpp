/* GENERAL_REMARK_BY_HILMAR:
 * When you decide on how the logger should work, this might also be related to this class (in the sense that this
 * pretty-printer-functionality is especially used for logging).
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_UTIL_PRETTY_PRINTER_HPP
#define KERNEL_UTIL_PRETTY_PRINTER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/exception.hpp>

// includes, system
#include <string>

namespace FEAST
{
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
  * Some construction rules: No blank is appended to the prefix (i.e., in the example above the prefix is "PRE ", and
  * not "PRE"). One blank is inserted between left delimiter and the string. Between the string and the right delimiter
  * there must be at least one blank, otherwise the right delimiter is omitted.
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
    String _block;

    /// width of the pretty printed block (from left to right delimiter, excluding the prefix)
    size_t _width;

    /// character to be used to enclose lines and to build separator lines
    char _delim;

    /// string each line is prefixed with (useful in connection with grep commands)
    String _prefix;


  public:

    /**
    * \brief CTOR using prefix
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
      size_t width,
      char delim,
      const String& prefix)
      : _block(""),
        _width(width),
        _delim(delim),
        _prefix(prefix)
    {
      CONTEXT("PrettyPrinter::PrettyPrinter()");
    }


    /**
    * \brief CTOR using an empty prefix
    *
    * \param[in] width
    * width of the lines in the pretty printed block
    *
    * \param[in] delim
    * character to be used to enclose lines and to build separator lines
    */
    PrettyPrinter(
      size_t width,
      char delim)
      : _block(""),
        _width(width),
        _delim(delim),
        _prefix("")
    {
      CONTEXT("PrettyPrinter::PrettyPrinter()");
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
    inline const String& block() const
    {
      CONTEXT("PrettyPrinter::block()");
      return _block;
    }

    /**
    * \brief setter for the pretty printed block
    */
    inline void set_block(String block)
    {
      CONTEXT("PrettyPrinter::set_block()");
      _block = block;
    }

    /**
    * \brief getter for the width of the pretty printed block
    *
    * \return width of the pretty printed block #_width
    */
    inline size_t width() const
    {
      CONTEXT("PrettyPrinter::width()");
      return _width;
    }

    /**
    * \brief setter for the width of the pretty printed block
    */
    inline void set_width(unsigned int width)
    {
      CONTEXT("PrettyPrinter::set_width()");
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
    inline const String& prefix() const
    {
      CONTEXT("PrettyPrinter::prefix()");
      return _prefix;
    }

    /**
    * \brief setter for the prefix
    */
    inline void set_prefix(String prefix)
    {
      CONTEXT("PrettyPrinter::set_prefix()");
      _prefix = prefix;
    }

    /**
    * \brief getter for the delimiter character
    *
    * \return delimiter #_delim
    */
    inline char delim() const
    {
      CONTEXT("PrettyPrinter::delim()");
      return _delim;
    }

    /**
    * \brief setter for the delimiter character
    */
    inline void set_delim(char delim)
    {
      CONTEXT("PrettyPrinter::set_delim()");
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
      CONTEXT("PrettyPrinter::add_line_sep()");
      _block += _prefix + String(_width, _delim) + "\n";
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
    void add_line_centered(const String& s)
    {
      CONTEXT("PrettyPrinter::add_line_centered()");
      if(_width < s.size() + 4)
      {
        // if the string is too long, produce a left aligned line without right delimiter
        add_line_no_right_delim(s);
      }
      else
      {
        // calculate blanks on the left and the right side
        size_t num_blanks_total = _width - s.size() - 2;
        size_t num_blanks_left(num_blanks_total / 2);
        size_t num_blanks_right(num_blanks_total % 2 == 0 ? num_blanks_left : num_blanks_left + 1);
        // add the line to the pretty printed block
        _block += _prefix + stringify(_delim) + String(num_blanks_left, ' ') + s
                  + String(num_blanks_right, ' ') + stringify(_delim) + "\n";
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
    void add_line(const String& s)
    {
      CONTEXT("PrettyPrinter::add_line()");
      if(_width < s.size() + 4)
      {
        // if the string is too long, produce a left aligned line without right delimiter
        add_line_no_right_delim(s);
      }
      else
      {
        // add the line to the pretty printed block
        _block += _prefix + stringify(_delim) + " " + s + String(_width - s.size() - 3, ' ')
                  + stringify(_delim) + "\n";
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
    void add_line_no_right_delim(const String& s)
    {
      CONTEXT("PrettyPrinter::add_line_no_right_delim()");
      _block += _prefix + stringify(_delim) + " " + s + "\n";
    }

    /**
    * \brief reset the pretty printed block to an emtpy string (equivalent to PrettyPrinter::set_block(""))
    */
    inline void reset_block()
    {
      CONTEXT("PrettyPrinter::reset_block()");
      _block = "";
    }

    /**
    * \brief print the pretty printed block to the given stream
    */
    void print(std::ostream& stream)
    {
      CONTEXT("PrettyPrinter::print()");
      stream << _block;
    }
  }; // class PrettyPrinter
} // namespace FEAST

#endif // KERNEL_UTIL_PRETTY_PRINTER_HPP
