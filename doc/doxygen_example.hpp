/**
 * \file doxygen_example.hpp
 * \brief demonstrates usage of doxygen
 *
 * This file demonstrates how to use doxygen documentation within the FEAST code.
 * \note This is work in progress! If someone finds further useful doxygen features or other improvements, please
 *       add them!
 */

#include <iostream>
#include <stdlib.h>
#include <string>

/**
 * \brief example class
 *
 * We use the JavaDoc comment style, i.e. with opening \c /** to identify the doxygen comment:
   \verbatim
   /**
    * \brief short description in one line, beginning with lower case, no full sentence, no period at the end
    *
    * Then the detailed description follows after a blank line. While the brief description uses the command
    * '\brief', the identifier '\detailed' for the detailed description is omitted. The detailed description is
    * optional (but highly recommended for important structures), and has to be written in full sentences.
    *
    * Further blank lines produce skips.
    *
    * Ignore the white space in the closing '* /' of this example comment.
    * I can't write it because that would end the C++/doxygen comment I'm currently writing...
    * (Don't know how to escape this.)
    * /
   \endverbatim
 * For comments without detailed description spanning only one line we allow the \c /// syntax in order not to
 * unnecessarily blow up the number of code lines:
   \verbatim
    /// description of member variable
    int some_member_variable
   \endverbatim
 * This is equivalent to
   \verbatim
   /**
    * \brief description of member variable
    * /
   int some_member_variable
   \endverbatim
 * See the doxygen manual for other special doxygen comment identifiers and be sure not to use them when you don't
 * want the specific comment to appear in the documentation.
 *
 * The doxygen comments are positioned directly above the item (class, struc, variable, function, ...) to be documented.
 *
 * Doxygen provides some special commands. The \c \\author command is mandatory in class comments and eventually in
 * function comments if the author of the function is not the author of the enclosing class. Other commands we could
 * use are:
 * \li \c \\bug for description of known bugs in a class or function (mentioning the issue ID of the described bug)
 * \li \c \\note, \c \\warning and \c \\todo for general hints and warnings like <em>'this does not really work yet',
 *      'that has not been tested in serial mode yet'</em> or <em>'this feature still has to be implemented'</em>
 * \li \c \\throw for the description of exceptions a class/function can throw
 *
 * These commands should be written at the end of the detailed description. We might add other useful commands
 * (I didn't read the entire doxygen manual yet...). Here are some examples for such commands.
 *
 * \author Arthur Dent
 * \bug This class does not work on tuesday.
 * \warning Always expect the unexpected, when using this class!
 * \note This is a note.
 * \todo Comments have to be improved.
 */
class DoxygenExample
{
  /* *********************************************************************************************
   * Eye-catching comments that are *not* to be included by doxygen should be written like this. *
   * Note the white space in the opening line '/* ***...' !!                                     *
   ***********************************************************************************************/

  /* other non-doxygen comments can be written like this */

  // or like this

  /*
   * or like this
   */

  /* *****************
   * private members *
   *******************/

  private:

    /* ************
     * attributes *
     **************/

    /// some private attribute
    int _some_private_attribute;

    /**
     * \brief specifier for some array size
     *
     * I really can't think of a useful detailed description right now.
     */
    int _size;

    /**
     * \brief some useful integer array
     *
     * This is an array of length #_size. (Note the linking feature.)
     * The array dimension has to be provided in the last line of the description. Examples:
     *
     * \li <em>Dimension:</em> [size]     (1D array with known size)
     * \li <em>Dimension:</em> [xsize][ysize]   (2D array with known size)
     * \li <em>Dimension:</em> [][]  (2D array with unknown size)
     * \li <em>Dimension:</em> [xsize][]   (2D array with partially known size)
     *
     * And now, as the last line, the actual dimension of this array.
     *
     * <em>Dimension:</em> [#_size]
     */
    int* _my_array;

    /* ******************
     * member functions *
     ********************/

    /**
     * \brief computes incredibly complicated stuff
     *
     * This is a long description of what the function does.
     * @param[in] useless_string input parameter to be ignored
     * @return a useful result
     * \sa answer_to_everything()
     */
    int _do_extremely_complicated_computation(std::string& useless_string)
    {
      return 42;
    }

  /* ****************
   * public members *
   ******************/

  public:

    /* ************
     * attributes *
     **************/

    /// for storing some useful information
    int i_am_public;

    /**
    * \brief for storing even more useful information
    *
    * This variable stores information that is even more useful than that of #i_am_public (note the linking feature).
    */
    int i_am_public_too;

    /* **************
     * constructors *
     ****************/

    /**
     * \brief simple constructor
     *
     * When using this constructor, then something happens.
     * @param[in] i some input parameter
     * @param[in] f another input parameter
     */
    DoxygenExample(int i, float f)
    {
      std::cout << "Constructor called with " << i << " and " << f << "." << std::endl;
    }

    /* ************
     * destructor *
     **************/

    /**
    * \brief simple destructor
    *
    * This is a simple destructor that annoys the user with some screen output.
    */
    ~DoxygenExample()
    {
      std::cout << "Object destroyed!" << std::endl;
    }

    /* ******************
     * members functions*
     ********************/

    /**
     * \brief gives the answer to everything
     *
     * This is a long description of what the function does.
     * @param[in] argc argument count passed to the main() method
     * @param[in] argv arguments passed to the main() method
     * @param[in] question arbitrary question
     * @return answer to the given question
     * \sa _do_extremely_complicated_computation()
     * \todo The detailed description needs to be improved.
     * \warning The computation takes roughly 10,000 years!
     */
    int answer_to_everything(int& argc, char* argv[], std::string& question)
    {
      return _do_extremely_complicated_computation(question);
    }
};
