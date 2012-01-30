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
   * This class organises a tree-like structure of key-value pairs, which is parsed from one or multiple data files
   * stored in a custom extension of the INI data file format; see e.g. http://en.wikipedia.org/wiki/INI_file.
   *
   * The INI data file format is a simple text file which stores data in so-called <em>key-value pairs</em>:
   * \verbatim key = value\endverbatim
   * Both the \e key and \e value parts are parsed verbatim and stored as simple strings, where the leading and
   * trailing whitespaces are trimmed by the parser.
   *
   * As an extension of the standard INI format, multiple consecutive lines may be joined by appending a
   * <c>&</c> character as the last non-whitespace character in a line:
   * \verbatim key = value1 &
      value2\endverbatim
   * is equivalent to
   * \verbatim key = value1 value2\endverbatim
   *
   * Several key-value pairs can be grouped together in <em>sections</em>. A section begins with a section name
   * enclosed in a <c>[ ]</c> bracket pair, followed by a single line containing a <c>{</c> curly brace which
   * marks the beginning of a section's body. The end of a section's body is marked by a single <c>}</c>
   * curly brace in a separate line. Please note that a section may not only contain key-value pairs but also
   * sub-sections:
     \verbatim
     [Section]
     {
       key = value
       another_key = another value
       [SubSection]
       {
         key = subsection value
       }
     }
     \endverbatim
   *
   * The hash sign <c>#</c> denotes the beginning of a comment spanning to the end of the line. A comment may
   * appear anywhere within a file -- even behind a key-value pair.
     \verbatim
     # This is a comment in a separate line
     [Section] # A comment after a section marker
     {
       key = value # a comment behind a key-value pair
     } # end of [Section]
     \endverbatim
   * \note Comments can not be 'split' across several lines by the \c & character.
   *
   * In analogy to the \c #include keyword of the C preprocessor, an INI file can include another INI file
   * by using the \c \@include keyword followed by the filename. Please note that the file path must be either
   * absolute or relative to the application's working directory. If an \c \@include command is encountered
   * within a section body, the included file's contents parsed into that section.
     \verbatim
     # include 'myfile.txt' in the unnamed root section
     @include myfile.txt
     [Section]
     {
       # include the contents of 'another_file.txt' into [Section]
       @include another_file.txt
     }
     \endverbatim
   *
   *
   * A more complex example:
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
         pi = 3.1415926535&
                8979323846&
                2643383279...
       }
       # The previous line closed 'MySubSection', the next one will close 'MySection'
     }

     # Include another data file.
     include@ /home/feast/my_data.txt
     \endverbatim
   *
   * \todo further doxygenisation; describe special cases
   *
   * \author Constantin Christof
   * \author Peter Zajac
   */
  class ParamSection
  {
  public:
    /**
     * \brief Syntax Error exception class
     *
     * This class derives from FEAST::Exception and is thrown by the ParamSection parser when a syntax error
     * is detected.
     *
     * \author Constantin Christof
     * \author Peter Zajac
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

  protected:
    /// value-map type
    typedef std::map<String, String, String::NoCaseLess> ValueMap;
    /// section-map type
    typedef std::map<String, ParamSection*, String::NoCaseLess> SectionMap;

  public:
    typedef ValueMap::iterator EntryIterator;
    typedef ValueMap::const_iterator ConstEntryIterator;
    typedef SectionMap::iterator SectionIterator;
    typedef SectionMap::const_iterator ConstSectionIterator;

  protected:
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
     *
     * \returns
     * \c true if the section was erased or \c false if no section with that name was found.
     */
    bool erase_section(String name);

    /**
     * \brief Erases a key-value pair.
     *
     * \param[in] key
     * The key of the entry to be erased.
     *
     * \returns
     * \c true if the entry was erased or \c false if no entry with that key was found.
     */
    bool erase_entry(String key);

    /// Returns the first entry iterator.
    EntryIterator begin_entry()
    {
      CONTEXT("ParamSection::begin_entry()");
      return _values.begin();
    }

    /** \copydoc begin_entry() */
    ConstEntryIterator begin_entry() const
    {
      CONTEXT("ParamSection::begin_entry() [const]");
      return _values.begin();
    }

    /// Returns the last entry iterator.
    EntryIterator end_entry()
    {
      CONTEXT("ParamSection::end_entry()");
      return _values.end();
    }

    /** \copydoc end_entry() */
    ConstEntryIterator end_entry() const
    {
      CONTEXT("ParamSection::end_entry() const");
      return _values.end();
    }

    /// Returns the first section iterator.
    SectionIterator begin_section()
    {
      CONTEXT("ParamSection::begin_section()");
      return _sections.begin();
    }

    /** \copydoc begin_section() */
    ConstSectionIterator begin_section() const
    {
      CONTEXT("ParamSection::begin_section() [const]");
      return _sections.begin();
    }

    /// Returns the last section iterator
    SectionIterator end_section()
    {
      CONTEXT("ParamSection::end_section()");
      return _sections.end();
    }

    /** \copydoc end_section() */
    ConstSectionIterator end_section() const
    {
      CONTEXT("ParamSection::end_section() [const]");
      return _sections.end();
    }

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
