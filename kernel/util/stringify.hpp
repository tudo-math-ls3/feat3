#ifndef UTIL_STRINGIFY_HHP
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
     */
    inline std::string stringify(const std::string & item)
    {
        return item;
    }

    /**
     * Convert item to a string (overload for char).
     */
    inline std::string stringify(const char & item)
    {
        return std::string(1, item);
    }

    /**
     * Convert item to a string (overload for unsigned char).
     */
    inline std::string stringify(const unsigned char & item)
    {
        return std::string(1, item);
    }

    /**
     * Convert item to a string (overload for bool).
     */
    inline std::string stringify(const bool & item)
    {
        return item ? "true" : "false";
    }

    /**
     * Convert item to a string (overload for char *, which isn't a
     * screwup like other pointers).
     */
    inline std::string stringify(const char * const item)
    {
        return std::string(item);
    }
}

#endif //UTIL_STRINGIFY_HHP
