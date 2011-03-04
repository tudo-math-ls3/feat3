#pragma once
#ifndef PREFIX_HPP
#define PREFIX_HPP

// includes, system
#include <string>
#include <vector>

// includes, FEAST
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>

/// FEAST namespace
namespace FEAST
{
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
  * objects and can be modified from everywhere (realised as singleton?)? Or is the prefix object
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
      CONTEXT("Prefix::Prefix()");
    }

    /// destructor
    ~Prefix()
    {
      CONTEXT("Prefix::~Prefix()");
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
      CONTEXT("Prefix::push()");
      // forbid pushing empty strings
      ASSERT(string_to_append.size() > 0, "String must not be of zero size.");
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
      CONTEXT("Prefix::pop()");
      // assert that there is a substring to be removed
      ASSERT(_start_pos.size() > 0, "There is no string to remove!");
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
      CONTEXT("Prefix::s()");
      return _s;
    }
  }; // class Prefix
} // namespace FEAST

#endif //  #ifndef PREFIX_HPP
