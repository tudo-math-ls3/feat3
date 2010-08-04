#ifndef UTIL_STRINGIFY_HHP
#define UTIL_STRINGIFY_HHP 1

#include <sstream>
#include <string>
#include <tr1/memory>

namespace Feast
{
    /**
     * For use by stringify.
     */
    namespace stringify_internals
    {
        /**
         * Check that T_ is a sane type to be stringified.
         */
        template <typename T_>
        struct CheckType
        {
            /// Yes, we are a sane type.
            enum { value = 0 } Value;
        };

        /**
         * Check that T_ is a sane type to be stringified.
         */
        template <typename T_>
        struct CheckType<T_ *>
        {
            /// Yes, we are a sane type.
            enum { value = 0 } Value;
        };

        /**
         * Check that T_ is a sane type to be stringified, which it isn't
         * if it's a shared_ptr.
         */
        template <typename T_>
        struct CheckType<std::tr1::shared_ptr<T_> >
        {
        };
    }

    /**
     * Convert item to a string.
     */
    template <typename T_>
    std::string stringify(const T_ & item)
    {
        /* check that we're not trying to stringify a pointer or somesuch */

        /// \todo Evaluate check_for_stringify_silly_things.
        //int check_for_stringifying_silly_things = stringify_internals::CheckType<T_>::value;

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
