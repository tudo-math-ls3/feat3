#pragma once
#ifndef KERNEL_UTIL_PARAM_SECTION
#define KERNEL_UTIL_PARAM_SECTION 1

// include FEAST
#include <kernel/util/string.hpp>
#include <kernel/util/file_error.hpp>

// include system
#include <fstream>
#include <iostream>
#include <stack>
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
   * This class allows to store and handle data, which are read from an INI-file. The data are saved
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
   *  - 'include@PATH' marks, that an external file given by 'PATH' should be read at this point.
   *  - '#' marks the beginning of a comment.
   *  - '&' at the end of a line marks a line break.
   *
   *   \attention
   *  - a line must not contain more than one statement (e.g., 'a = 1  b = 4' in the same line would provoke an error).
   *  - each '[ name ]' must be followed by a '{' in the next line.
   *  - a comment cannot be continued by using '&', each comment line has to start with '#'.
   *  - an array is saved as an ordinary key-value pair, whose value is divided into several parts by a seperator.
   *  - the default seperator is ",".
   *  - a value or an 'include@', that is not assigned to any section (e.g., value1, value10, value11 in the example  below), is associated with the ParamSection, the 'parse' function is called from.
   *  - a section's name or a key must not be an empty string.
   *
   * Example:
   *
   * \verbatim

     value1 = 1

     [A]
     {
       array = 00/11/22/33/44/55/66/77/&     #Comment
               88/99/100/101/102&            #Comment
               /111/222

       value2 = 2
       value3 = 3

       [B]
       {
         value4 = 4
         [C]
         {
           value5 = 5
           value6 = 6
         }
       }

       [D]
       {
         include@/home/user/.../test
         value7 = string1
         value8 = string2
       }
     }

     [E]
     {
       value9 = 9
     }

     value10 = 42
     value11 = 23

     \endverbatim
   *
   * \author Constantin Christof
   */

  class ParamSection
  {

  public:
    class SyntaxError :
      public Exception
    {
    protected:
      String _filename;
    public:
      explicit SyntaxError(
        String message,
        String filename = ""):
        Exception(message + (filename.empty() ? "": " in file " + filename)),
        _filename(filename)
      {
      }
      virtual ~SyntaxError() throw()
      {
      }
      String get_filename() const
      {
        return _filename;
      }
    };

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


    //input//

  /**
   * \brief Adds a new key-value pair to _values.
   *
   * \param[in] key contains the pair's key.
   * \param[in] value contains the pair's value.
   * \param[in] replace marks, if an already existing pair with the name \p key should be replaced by
                the new one (replace is set to \c true by default)
   * \returns
   * \c true, if the key-value-pair has been stored or \c false, if an entry with the key already exists and
   * \p replace was set to \c false.
   * \sa add_section()
   */
    bool add_entry(String key, String value, bool replace = true);


  /**
   * \brief Creates a new section with the name '\p name' in _sections and returns a pointer to this entry
   *
   * \details If the section '\p name' already exists, the function returns a pointer to this entry.
   *
   * \param[in] name contains the new section's name
   * \return a ParamSection-Pointer to the new entry in _sections
   *
   * \sa add_entry()
   */
    ParamSection* add_section(String name);

  /**
   * \brief Parses the INI-file given by \p filename and inserts the data into the maps _sections and _values
   *
   * \param[in] filename contains the path to the INI-file
   */
    void parse(String filename);


  /**
   * \brief Parses the stream given by \p &f and inserts the data into the maps _sections and _values
   *
   * \param[in] f is the reference to the istream object, the parsed data is read from
   */
    void parse(std::istream& f);

  /**
   * \brief Merges 'this' ParamSection with \p section
   *
   * \param[in] section is the ParamSection, which should be included
   * \param[in] replace marks, if already existing key-value pairs should be overwritten
   */
    void merge(ParamSection* section, bool replace = true);

  /**
   * \brief Searches the map _values for an entry with the name '\p key' and returns the related value
   *
   * \param[in] key contains the searched key
   * \return a pair <String, bool>, where the second component marks, if the key has been found (\c false = not found),
   *  and the first entry contains the desired value (or is empty if there is no entry 'key')
   * \sa get_section()
   */
    std::pair<String, bool> get_entry(String key) const;


  /**
   * \brief Searches the map _values for an array named '\p key' and returns the array entry \p pos, if the array was   found.
   *
   * \param[in] key contains the searched key / the name of the array
   * \param[in] pos marks the position of the required entry within the array \p key
   * \param[in] sep defines the seperator that is used to structure the array (sep is set to ',' by default)
   * \return a pair <String, bool>, where the second component marks, if the access has been successful (\c false =    error occurred), and the first entry contains the desired value (or is empty if something went wrong)
   * \sa get_section()
   */
    std::pair<String, bool> get_entry(String key, String::size_type pos, String sep = ",") const;


  /**
   * \brief Searches the map _sections for an entry named '\p name' and returns a pointer to this entry (or nullptr)
   *
   * \param[in] name contains the searched section's name
   * \return a ParamSection- pointer to the element with the searched name or nullptr, if '\p name' has not been found
   * \sa get_entry()
   */
    ParamSection* get_section(String name);


  /**
   * \brief Searches the map _sections for an entry with the name '\p name' and returns a const pointer to this entry.
   *
   * \details If the desired entry could not be found, the function returns nullptr.
   *
   * \param[in] name contains the searched section's name
   * \return a const ParamSection- pointer to the entry with the searched name or nullptr, if '\p name' has not been    found
   * \sa get_entry()
   */
    const ParamSection* get_section(String name) const;


  /**
   * \brief Displays the content of the map _values via std::cout.
   *
   * \sa display_sections()
   */
    void display_values() const;

  /**
   * \brief Displays the content of the map _sections via std::cout.
   *
   * \sa display_values()
  */
    void display_sections() const;


  /**
   * \brief Dumps the section tree into an output stream.
   *
   * \details Passes the data structure to the ostream object \p os, adds a number of whitespaces in front of each
   * line, which is given by \p indent
   *
   * \param[in,out] os is the reference to the ostream object
   * \param[in] indent marks, how many whitespaces should be added
   */
    void dump(std::ostream& os, String::size_type indent = 0) const;

  /**
   * \brief Writes the data structure in a parsable INI-file.
   *
   * \param filename contains the path to the file, the data should be written in
   */
    void dump(String filename) const;


    //misc//

  /**
   * \brief Deletes the section given by \p name (and all of its content)
   *
   * \param[in] name contains the section's name
   * \sa erase_entry()
   */
    void erase_section(String name);

  /**
   * \brief Deletes the key-value pair given by 'key'
   *
   * \param[in] key contains the pair's key
   * \sa erase_section()
   */
    void erase_entry(String key);

  }; // class ParamSection
} // namespace FEAST
#endif
