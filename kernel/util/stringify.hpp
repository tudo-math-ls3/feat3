#pragma once
#ifndef UTIL_STRINGIFY_HHP
/// Header guard
#define UTIL_STRINGIFY_HHP 1

#include <sstream>
#include <string>

/**
 * \file stringify.hpp
 * \brief Methods to stringify different datatypes
 *
 * \author Dirk Ribbrock
 */

namespace Feast
{
    /**
     * Convert item to a string.
     * \param[in] item
     * The item to stringify
     */
    template <typename T_>
    std::string stringify(const T_ & item)
    {
        std::ostringstream s;
        s << item;
        return s.str();
    }

    /**
     * Convert item to a string (overload for std::string).
     * \param[in] item
     * The item to stringify
     */
    inline std::string stringify(const std::string & item)
    {
        return item;
    }

    /**
     * Convert item to a string (overload for char).
     * \param[in] item
     * The item to stringify
     */
    inline std::string stringify(const char & item)
    {
        return std::string(1, item);
    }

    /**
     * Convert item to a string (overload for unsigned char).
     * \param[in] item
     * The item to stringify
     */
    inline std::string stringify(const unsigned char & item)
    {
        return std::string(1, item);
    }

    /**
     * Convert item to a string (overload for bool).
     * \param[in] item
     * The item to stringify
     */
    inline std::string stringify(const bool & item)
    {
        return item ? "true" : "false";
    }

    /**
     * Convert item to a string (overload for char *, which isn't a
     * screwup like other pointers).
     * \param[in] item
     * The item to stringify
     */
    inline std::string stringify(const char * const item)
    {
        return std::string(item);
    }
}

#endif //UTIL_STRINGIFY_HHP
