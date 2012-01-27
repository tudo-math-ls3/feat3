#pragma once
#ifndef KERNEL_UTIL_PARAM_SECTION_HPP
#define KERNEL_UTIL_PARAM_SECTION_HPP 1

// includes, FEAST
#include <kernel/util/string.hpp>
#include <kernel/util/file_error.hpp>

// includes, system
#include <iostream>
#include <map>

//////////////////////////////////////
/////////// Section Struct ///////////
//////////////////////////////////////

namespace FEAST
{

  /**
   * \brief A class for storing data read from an INI-file
   *
   * \details
   * This class allows to parse and store data, which is read from an INI-file. The data is saved
   * within an tree-like structure.
   *
   * \internal
   *
   * Each section consists of a map called _values (map <String, String, NoCaseLess>), which contains
   * its key-value pairs, and a map named _sections (map<String, ParamSection*, NoCaseLess>), which
   * contains the names of the subordinated sections and pointers to these structs.
   *
   * \endinternal
   *
   * The supported INI-file format is as follows:<br>
   *  - '[ name ]' marks the beginning of a new section called 'name'.
   *  - '{' marks the beginning of a section's content.
   *  - '}' marks the end of a section.
   *  - '=' marks a key-value pair, that is defined by the lefthand side and the righthand side of the equals sign.
   *  - 'include\@PATH' marks, that an external file given by 'PATH' should be read at this point.
   *  - '#' marks the beginning of a comment.
   *  - '&' at the end of a line marks a line break.
   *
   *   \attention
   *  - a line must not contain more than one statement (e.g., 'a = 1  b = 4' in the same line would provoke an error).
   *  - each '[ name ]' must be followed by a '{' in the next line.
   *  - a comment cannot be continued by using '&', each comment line has to start with '#'.
   *  - a value or an 'include@', that is not assigned to any section (e.g., value1, value10, value11 in the example
   *    below), is associated with the ParamSection, the 'parse' function is called from.
   *  - a section's name or a key must not be an empty string.
   *
   * Example:
   *
   * \verbatim
     # This is a comment and the next line is a simple key-value pair in the root section.
     key = value

     # The next two lines open a new section named 'MySection'.
     [MySection]
     {
       # Two key-value pairs in section MySection.
       date = 1985-05-08
       message = Hello World! # This is a comment after a key-value pair

       # A subsection named 'MySubSection'
       [MySubSection]
       {
         # A value spanning over multiple lines
         pi = 3.1415926535 &
                8979323846 &
                2643383279...
       }
       # The previous line closed 'MySubSection', the next one will close 'MySection'
     }

     # Include another data file.
     include@ ./mydata.txt
     \endverbatim
   *
   * \author Constantin Christof
   */

  class ParamSection
  {
  public:
    /**
     * \brief Syntax Error exception class
     *
     * This class derives from FEAST::Exception and is thrown by the ParamSection parser when a syntax error
     * is detected.
     */
    class SyntaxError :
      public Exception
    {
    protected:
      /// name of the file containing the syntax error
      String _filename;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] message
       * A description of the syntax error.
       *
       * \param[in] filename
       * The name of the file in which the syntax error has been detected.
       */
      explicit SyntaxError(
        String message,
        String filename = ""):
        Exception(message + (filename.empty() ? "": " in file " + filename)),
        _filename(filename)
      {
      }

      /// virtual destructor
      virtual ~SyntaxError() throw()
      {
      }

      /// returns the filename
      String get_filename() const
      {
        return _filename;
      }
    }; // class ParamSection::SyntaxError

  private:
    /// value-map type
    typedef std::map<String, String, String::NoCaseLess> ValueMap;
    /// section-map type
    typedef std::map<String, ParamSection*, String::NoCaseLess> SectionMap;

    /// a map storing the key-value-pairs
    ValueMap _values;

    /// a map storing the sub-sections
    SectionMap _sections;

  public:
    /// Default Constructor
    ParamSection();

    /// Virtual Destructor
    virtual ~ParamSection();

    /**
     * \brief Adds a new key-value pair to this section.
     *
     * This function adds a new key-value pair to this section or, if desired, overwrites an existing one.
     *
     * \param[in] key
     * A String containing the key of the entry.
     *
     * \param[in] value
     * A String containing the value of the entry.
     *
     * \param[in] replace
     * Specifies the behaviour when an entry with the same key already exists in the section:
     *  - If \p replace = \c true, then the old value of the entry is replaced by the \p value parameter passed
     *    to this function.
     *  - If \p replace = \c false, then the old value of the entry is kept and this function returns \c false.
     *
     * \returns
     * \c true, if the key-value-pair has been stored or \c false, if an entry with the key already exists and
     * \p replace was set to \c false.
     */
    bool add_entry(String key, String value, bool replace = true);

    /**
     * \brief Adds a new sub-section to this section.
     *
     * This function adds a new sub-section and returns a pointer to it. If a sub-section with the corresponding
     * name already exists, a pointer to that section is returned instead of allocating a new one.
     *
     * \param[in] name
     * The name of the section to be added.
     *
     * \return
     * A pointer to a the \c ParamSection associated to the name.
     */
    ParamSection* add_section(String name);

    /**
     * \brief Erases a sub-section
     *
     * \param[in] name
     * The name of the sub-section to be erased.
     */
    void erase_section(String name);

    /**
     * \brief Erases a key-value pair.
     *
     * \param[in] key
     * The key of the entry to be erased.
     */
    void erase_entry(String key);

    /**
     * \brief Parses a file in INI-format.
     *
     * \param[in] filename
     * The name of the file to be parsed.
     */
    void parse(String filename);

    /**
     * \brief Parses an input stream.
     *
     * \param[in] ifs
     * A reference to an input stream to be parsed.
     */
    void parse(std::istream& ifs);

    /**
     * \brief Merges another ParamSection into \c this.
     *
     * \param[in] section
     * The section to be merged into \c this.
     *
     * \param[in] replace
     * Specifies the behaviour when conflicting key-value pairs are encountered. See add_entry() for more details.
     */
    void merge(const ParamSection& section, bool replace = true);

    /**
     * \brief Retrieves a value for a given key.
     *
     * This function returns the value string of a key-value pair.
     *
     * \param[in] key
     * The key of the entry whose value is to be returned.
     *
     * \returns
     * A pair<String, bool>, where the second component marks, whether an entry with the \c key has been found or not.
     * If the \c bool-component is \c true, then the String-component contains the value associated with the key,
     * otherwise the String-component is empty.
     */
    std::pair<String, bool> get_entry(String key) const;

    /**
     * \brief Returns a sub-section.
     *
     * \param[in] name
     * The name of the sub-section which is to be returned.
     *
     * \returns
     * A pointer to the ParamSection associated with \p name or \c nullptr if no section with that name exists.
     */
    ParamSection* get_section(String name);

    /** \copydoc get_section() */
    const ParamSection* get_section(String name) const;

    /**
     * \brief Dumps the section tree into an output stream.
     *
     * This function writes the whole section tree into an output stream, which can be parsed.
     *
     * \param[in,out] os
     * A refernce to an output stream to which to write to. It is silently assumed that the output stream is
     * open.
     *
     * \param[in] indent
     * Specifies the number of indenting whitespaces to be inserted at the beginning of each line.
     *
     * \note
     * The \p indent parameter is used internally for output formatting and does not need to be specified by the
     * caller.
     */
    void dump(std::ostream& os, String::size_type indent = 0) const;

    /**
     * \brief Dumps the section tree into a file.
     *
     * This function writes the whole section tree into a file, which can be parsed by the parse() function.
     *
     * \param[in] filename
     * The name of the file into which to dump to.
     *
     * \see dump(std::ostream&, String::size_type)
     */
    void dump(String filename) const;

  }; // class ParamSection
} // namespace FEAST

#endif // KERNEL_UTIL_PARAM_SECTION_HPP
