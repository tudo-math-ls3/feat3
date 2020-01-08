// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_PROPERTY_MAP_HPP
#define KERNEL_UTIL_PROPERTY_MAP_HPP 1

// includes, FEAT
#include <kernel/util/exception.hpp>
#include <kernel/util/dist.hpp>

// includes, system
#include <iostream>
#include <map>

namespace FEAT
{
  /**
   * \brief A class organising a tree of key-value pairs
   *
   * \see \ref ini_format
   *
   * <b>Path Specification:</b>\n
   * Some functions of this class, namely the #query() and
   * #query_section() functions, use \e paths for their search. Paths allow to search for
   * sections and key-value pairs (either relatively or absolutely) within the whole
   * property map tree instead of only within the current property map represented by
   * the ProperyMap object itself.
   *
   * For this, three characters are reserved as special characters in the path:
   * - The forward slash \c / represents a path separator (just as in file-system paths).
   * - The exclamation mark \c ! represents the root section.
   * - The tilde \c ~ represents the parent section (as the double-dot <c>..</c> in file-system paths).
   *
   * Examples:
   * - <c>SecName/KeyName</c>: searches for a key named \c KeyName inside a sub-section named \c SecName
   *   within the current section.
   * - <c>!/SecName/KeyName</c>: searches for key named \c KeyName inside a sub-section named \c SecName
   *   within the root section; independent of which section is represented by the current section.
   * - <c>~/KeyName</c>: searches for a key named \c KeyName inside the parent section of the
   *   current section.
   *
   * \author Constantin Christof
   * \author Peter Zajac
   */
  class PropertyMap
  {
  public:
    /// entry-map type
    typedef std::map<String, String, String::NoCaseLess> EntryMap;
    /// section-map type
    typedef std::map<String, PropertyMap*, String::NoCaseLess> SectionMap;

    /// entry iterator type
    typedef EntryMap::iterator EntryIterator;
    /// const entry iterator type
    typedef EntryMap::const_iterator ConstEntryIterator;
    /// section iterator type
    typedef SectionMap::iterator SectionIterator;
    /// const section iterator type
    typedef SectionMap::const_iterator ConstSectionIterator;

  protected:
    /// pointer to the parent node
    PropertyMap* _parent;

    /// a map storing the key-value-pairs
    EntryMap _values;

    /// a map storing the sub-sections
    SectionMap _sections;

  public:
    /// Default Constructor
    explicit PropertyMap(PropertyMap* parent = nullptr);
    /// Delete copy constructor
    PropertyMap(const PropertyMap&) = delete;
    /// Delete copy assignment operator
    PropertyMap& operator=(const PropertyMap&) = delete;

    /// Virtual Destructor
    virtual ~PropertyMap();

    /**
     * \brief Adds a new key-value pair to this map.
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
     * \brief Adds a new sub-section to this map.
     *
     * This function adds a new sub-section and returns a pointer to it. If a sub-section with the corresponding
     * name already exists, a pointer to that section is returned instead of allocating a new one.
     *
     * \param[in] name
     * The name of the section to be added.
     *
     * \return
     * A pointer to a the \c PropertyMap associated to the name.
     */
    PropertyMap* add_section(String name);

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
     * \brief Queries a value by its key path.
     *
     * \note
     * See the documentation of this class for details about paths.
     *
     * \param[in] key_path
     * A path to the key whose value is to be returned.
     *
     * \returns
     * A <c>pair<String,bool></c>, where the second component marks, whether an entry with the key path
     * has been found or not. If the <c>bool</c>-component is \c true, then the <c>String</c>-component
     * contains the value associated with the key, otherwise the <c>String</c>-component is empty.
     */
    std::pair<String, bool> query(String key_path) const;

    /**
     * \brief Queries a value by its key path.
     *
     * \note
     * See the documentation of this class for details about paths.
     *
     * \param[in] key_path
     * A path to the key whose value is to be returned.
     *
     * \param[in] default_value
     * A string that is to be returned in the case that the key was not found.
     *
     * \returns
     * A String containing the value corresponding to the key-path, or \p default_value if no such key was found.
     */
    String query(String key_path, String default_value) const;

    /**
     * \brief Queries a section by its section path.
     *
     * \note
     * See the documentation of this class for details about paths.
     *
     * \param[in] sec_path
     * A path to the section that is to be found.
     *
     * \returns
     * A (const) pointer to the section represented by \p path or \c nullptr if no such section was found.
     */
    PropertyMap* query_section(String sec_path);

    /** \copydoc query_section() */
    const PropertyMap* query_section(String sec_path) const;

    /**
     * \brief Retrieves a value by its key.
     *
     * \note
     * This function only searches the current section represented by \c this for the key.
     * If you want to search whole paths, use the #query() function instead.
     *
     * \param[in] key
     * The key of the entry whose value is to be returned.
     *
     * \returns
     * A <c>pair<String, bool></c>, where the second component marks, whether an entry with the \c key has been
     * found or not. If the <c>bool</c>-component is \c true, then the <c>String</c>-component contains the value
     * associated with the key, otherwise the <c>String</c>-component is empty.
     */
    std::pair<String, bool> get_entry(String key) const;

    /**
     * \brief Retrieves a sub-section by its name.
     *
     * \note
     * This function only search the current section represented by \c this for the sub-section.
     * If you want to search whole paths, use the #query_section() function instead.
     *
     * \param[in] name
     * The name of the sub-section which is to be returned.
     *
     * \returns
     * A pointer to the PropertyMap associated with \p name or \c nullptr if no section with that name exists.
     *
     * \note
     * This method does not allocate new memory, thus modifications to the returned PropertyMap affect the base
     * property map and the returned PropertyMap <b>must</b> not be deleted.
     */
    PropertyMap* get_sub_section(String name);

    /** \copydoc get_section() */
    const PropertyMap* get_sub_section(String name) const;

    /**
     * \brief Parses an object from the value of an entry
     *
     * If an entry for the given key exists, this function is equivalent to
     * <tt>this->get_entry(key).first.parse(x)</tt>.
     *
     * \param[in] key
     * The key of the entry whose value is to be parsed.
     *
     * \param[out] x
     * A reference to the object that is to be parsed into.
     *
     * \param[in] rtn_non_exist
     * The value that this function should return if there exisits no entry for the given key.
     *
     * \returns
     * - \c true, if the entry exists and was parsed successfully,
     * - \c false, if the entry exists but could not be parsed,
     * - \c rtn_non_exist, if the entry does not exist
     */
    template<typename T_>
    bool parse_entry(const String& key, T_& x, bool rtn_non_exist = true) const
    {
      std::pair<String, bool> it = this->get_entry(key);
      if(it.second)
        return it.first.parse(x);
      else
        return rtn_non_exist;
    }

    /**
     * \brief Returns a pointer to the parent section.
     *
     * \returns
     * A pointer to the parent section or \c nullptr if this section has no parent.
     */
    PropertyMap* get_parent()
    {
      return _parent;
    }

    /** \copydoc get_parent() */
    const PropertyMap* get_parent() const
    {
      return _parent;
    }

    /**
     * \brief Returns a pointer to the root section.
     */
    PropertyMap* get_root()
    {
      return (_parent != nullptr) ? _parent->get_root() : this;
    }

    /** \copydoc get_root() */
    const PropertyMap* get_root() const
    {
      return (_parent != nullptr) ? _parent->get_root() : this;
    }

    /**
     * \brief Returns a reference to the entry map.
     * \returns
     * A (const) reference to the entry map.
     */
    EntryMap& get_entry_map()
    {
      return _values;
    }

    /** \copydoc get_entry_map() */
    const EntryMap& get_entry_map() const
    {
      return _values;
    }

    /// Returns the first entry iterator.
    EntryIterator begin_entry()
    {
      return _values.begin();
    }

    /** \copydoc begin_entry() */
    ConstEntryIterator begin_entry() const
    {
      return _values.begin();
    }

    /// Returns the last entry iterator.
    EntryIterator end_entry()
    {
      return _values.end();
    }

    /** \copydoc end_entry() */
    ConstEntryIterator end_entry() const
    {
      return _values.end();
    }

    /// Returns the first section iterator.
    SectionIterator begin_section()
    {
      return _sections.begin();
    }

    /** \copydoc begin_section() */
    ConstSectionIterator begin_section() const
    {
      return _sections.begin();
    }

    /// Returns the last section iterator
    SectionIterator end_section()
    {
      return _sections.end();
    }

    /** \copydoc end_section() */
    ConstSectionIterator end_section() const
    {
      return _sections.end();
    }

    /**
     * \brief Merges another PropertyMap into \c this.
     *
     * \param[in] section
     * The section to be merged into \c this.
     *
     * \param[in] replace
     * Specifies the behaviour when conflicting key-value pairs are encountered. See add_entry() for more details.
     */
    void merge(const PropertyMap& section, bool replace = true);

    /**
     * \brief Parses a file in INI-format.
     *
     * \param[in] filename
     * The name of the file to be parsed.
     *
     * \param[in] replace
     * Specifies the behaviour when an entry with the same key already exists:
     *  - If \p replace = \c true, then the old value of the entry is replaced by the new \p value
     *  - If \p replace = \c false, then the old value of the entry is kept.
     *
     * \see PropertyMap::read(std::istream&)
     */
    void read(String filename, bool replace = true)
    {
      Dist::Comm comm(Dist::Comm::world());
      read(comm, filename, replace);
    }

    /**
     * \brief Parses a file in INI-format.
     *
     * \param[in] filename
     * The name of the file to be parsed.
     *
     * \param[in] replace
     * Specifies the behaviour when an entry with the same key already exists:
     *  - If \p replace = \c true, then the old value of the entry is replaced by the new \p value
     *  - If \p replace = \c false, then the old value of the entry is kept.
     *
     * \see PropertyMap::read(std::istream&)
     */
    void read(const Dist::Comm& comm, String filename, bool replace = true);

    /**
     * \brief Parses an input stream.
     *
     * \param[in] ifs
     * A reference to an input stream to be parsed.
     *
     * \param[in] replace
     * Specifies the behaviour when an entry with the same key already exists:
     *  - If \p replace = \c true, then the old value of the entry is replaced by the new \p value
     *  - If \p replace = \c false, then the old value of the entry is kept.
     */
    void read(std::istream& ifs, bool replace = true);

    /**
     * \brief Writes the property map into a file.
     *
     * This function writes the whole property map into a file, which can be parsed by the read() function.
     *
     * \param[in] filename
     * The name of the file into which to dump to.
     *
     * \see PropertyMap::write(std::ostream&, String::size_type) const
     */
    void write(String filename) const
    {
      Dist::Comm comm(Dist::Comm::world());
      write(comm, filename);
    }

    /**
     * \brief Writes the property map into a file.
     *
     * This function writes the whole property map into a file, which can be parsed by the read() function.
     *
     * \param[in] filename
     * The name of the file into which to dump to.
     *
     * \see PropertyMap::write(std::ostream&, String::size_type) const
     */
    void write(const Dist::Comm& comm, String filename) const;

    /**
     * \brief Writes the property map into an output stream.
     *
     * This function writes the whole property map into an output stream, which can be parsed.
     *
     * \param[in,out] os
     * A reference to an output stream to which to write to. It is silently assumed that the output stream is
     * open.
     *
     * \param[in] indent
     * Specifies the number of indenting whitespaces to be inserted at the beginning of each line.
     *
     * \note
     * The \p indent parameter is used internally for output formatting and does not need to be specified by the
     * caller.
     */
    void write(std::ostream& os, String::size_type indent = 0) const;
  }; // class PropertyMap
} // namespace FEAT

#endif // KERNEL_UTIL_PROPERTY_MAP_HPP
