// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_DIST_HPP
#define KERNEL_UTIL_DIST_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/type_traits.hpp>

// includes, system
#include <vector>
#include <sstream>

// includes, MPI
#ifdef FEAT_HAVE_MPI
#include <mpi.h>
#endif // FEAT_HAVE_MPI

namespace FEAT
{
  /**
   * \brief Distributed Communication namespace
   *
   * This namespaces encapsulates the classes required for distributed communication
   * via the Message Passing Interface (MPI).
   */
  namespace Dist
  {
    /**
     * \brief Initializes the distributed communication system.
     *
     * This function is effectively a wrapper around the \c MPI_Init function.
     *
     * In addition, this function may perform further setup to initialize
     * additional datatypes, operations, etc.:
     *
     * - If FEAT_HAVE_QUADMATH is defined, this function will initialize the
     *   dt__float128 datatype.
     * - If FEAT_HAVE_HALFMATH is defined, this function will initialize the
     *   dt__half datatype.
     * - If FEAT_OVERRIDE_MPI_OPS is defined, this function will override
     *   the standard op_sum, op_max and op_min operations by custom
     *   implementations.
     *
     * \see \cite MPI31 Section 8.7, page 355
     */
    bool initialize(int& argc, char**& argv);

    /**
     * \brief Finalizes the distributed communication system.
     *
     * This function is effectively a wrapper around the \c MPI_Init function.
     *
     * In addition, this function may perform further cleanup to release
     * additionally defined datatype, operations, etc.
     *
     * \see \cite MPI31 Section 8.7, page 357
     */
    void finalize();

    /**
     * \brief Communication Datatype class
     *
     * This class effectively wraps around the \c MPI_Datatype handle.
     * Instances of this class are used to specify the datatype of messages.
     *
     * \author Peter Zajac
     */
    struct Datatype
    {
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// the MPI datatype handle
      MPI_Datatype dt;
      /// the size of the datatype in bytes
      std::size_t bytes;

      /**
       * \brief \c MPI_Datatype handle constructor
       *
       * \param[in] dt_
       * The \c MPI_Datatype handle that is to be wrapped.
       *
       * \param[in] bytes_
       * The size of the datatype in bytes.
       */
      explicit Datatype(MPI_Datatype dt_, std::size_t bytes_) : dt(dt_), bytes(bytes_) {}
#else
    private:
      int dt;
      std::size_t bytes;

    public:
      explicit Datatype(int dt_, std::size_t bytes_) : dt(dt_), bytes(bytes_) {}
#endif // FEAT_HAVE_MPI

    public:
      /// equality comparison operator
      bool operator==(const Datatype& other) const
      {
        return this->dt == other.dt;
      }

      /// inequality comparison operator
      bool operator!=(const Datatype& other) const
      {
        return this->dt != other.dt;
      }

      /**
       * \brief Returns the size of the datatype in bytes, as given in the constructor
       */
      std::size_t size() const
      {
        return bytes;
      }
    }; // class Datatype

    /// Datatype wrapper for \c MPI_BYTE
    extern const Datatype dt_byte;
    /// Datatype wrapper for \c MPI_CHAR
    extern const Datatype dt_char;
    /// Datatype wrapper for \c MPI_WCHAR
    extern const Datatype dt_wchar;
    /// Datatype wrapper for \c MPI_SIGNED_CHAR
    extern const Datatype dt_signed_char;
    /// Datatype wrapper for \c MPI_SHORT
    extern const Datatype dt_signed_short;
    /// Datatype wrapper for \c MPI_INT
    extern const Datatype dt_signed_int;
    /// Datatype wrapper for \c MPI_LONG
    extern const Datatype dt_signed_long;
    /// Datatype wrapper for \c MPI_LONG_LONG
    extern const Datatype dt_signed_long_long;
    /// Datatype wrapper for \c MPI_UNSIGNED_CHAR
    extern const Datatype dt_unsigned_char;
    /// Datatype wrapper for \c MPI_UNSIGNED_SHORT
    extern const Datatype dt_unsigned_short;
    /// Datatype wrapper for \c MPI_UNSIGNED
    extern const Datatype dt_unsigned_int;
    /// Datatype wrapper for \c MPI_UNSIGNED_LONG
    extern const Datatype dt_unsigned_long;
    /// Datatype wrapper for \c MPI_UNSIGNED_LONG_LONG
    extern const Datatype dt_unsigned_long_long;
    /// Datatype wrapper for \c MPI_FLOAT
    extern const Datatype dt_float;
    /// Datatype wrapper for \c MPI_DOUBLE
    extern const Datatype dt_double;
    /// Datatype wrapper for \c MPI_LONG_DOUBLE
    extern const Datatype dt_long_double;
    /// Datatype wrapper for \c MPI_INT8_T
    extern const Datatype dt_signed_int8;
    /// Datatype wrapper for \c MPI_INT16_T
    extern const Datatype dt_signed_int16;
    /// Datatype wrapper for \c MPI_INT32_T
    extern const Datatype dt_signed_int32;
    /// Datatype wrapper for \c MPI_INT64_T
    extern const Datatype dt_signed_int64;
    /// Datatype wrapper for \c MPI_UINT8_T
    extern const Datatype dt_unsigned_int8;
    /// Datatype wrapper for \c MPI_UINT16_T
    extern const Datatype dt_unsigned_int16;
    /// Datatype wrapper for \c MPI_UINT32_T
    extern const Datatype dt_unsigned_int32;
    /// Datatype wrapper for \c MPI_UINT64_T
    extern const Datatype dt_unsigned_int64;

#if defined(FEAT_HAVE_QUADMATH) || defined(DOXYGEN)
    /// custom Datatype for __float128
    extern const Datatype dt__float128;
#endif

#if defined(FEAT_HAVE_HALFMATH) || defined(DOXYGEN)
    /// custom Datatype for half_float::half
    extern const Datatype dt__half;
#endif

    /**
     * \brief Automatic Datatype deduction function template
     *
     * This function returns the Datatype object for any of the fundamental C/C++ datatypes.
     *
     * \tparam T_
     * The fundamental C/C++ datatype whose corresponding Datatype object is to be returned.
     *
     * \returns
     * A reference to the <c>dt_*</c> Datatype object representing the type \c T_.
     */
    template<typename T_> const Datatype& autotype();

    /// \cond internal
    template<> inline const Datatype& autotype<char>()              {return dt_char;}
    template<> inline const Datatype& autotype<wchar_t>()           {return dt_wchar;}
    template<> inline const Datatype& autotype<signed char>()       {return dt_signed_char;}
    template<> inline const Datatype& autotype<signed short>()      {return dt_signed_short;}
    template<> inline const Datatype& autotype<signed int>()        {return dt_signed_int;}
    template<> inline const Datatype& autotype<signed long>()       {return dt_signed_long;}
    template<> inline const Datatype& autotype<signed long long>()  {return dt_signed_long_long;}
    template<> inline const Datatype& autotype<unsigned char>()     {return dt_unsigned_char;}
    template<> inline const Datatype& autotype<unsigned short>()    {return dt_unsigned_short;}
    template<> inline const Datatype& autotype<unsigned int>()      {return dt_unsigned_int;}
    template<> inline const Datatype& autotype<unsigned long>()     {return dt_unsigned_long;}
    template<> inline const Datatype& autotype<unsigned long long>(){return dt_unsigned_long_long;}
    template<> inline const Datatype& autotype<float>()             {return dt_float;}
    template<> inline const Datatype& autotype<double>()            {return dt_double;}
    template<> inline const Datatype& autotype<long double>()       {return dt_long_double;}

#if defined(FEAT_HAVE_QUADMATH) || defined(DOXYGEN)
    template<> inline const Datatype& autotype<__float128>()        {return dt__float128;}
#endif

#if defined(FEAT_HAVE_HALFMATH) || defined(DOXYGEN)
    template<> inline const Datatype& autotype<half_float::half>()  {return dt__half;}
#endif

    /// \endcond

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    /**
     * \brief Communication Operation class
     *
     * This class effectively wraps around the \c MPI_Op handle.
     * Instances of this class are used to specify the operation for collective reduction messages.
     *
     * \author Peter Zajac
     */
    struct Operation
    {
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// the MPI operation handle
      MPI_Op op;

      /**
       * \brief \c MPI_Op handle constructor
       *
       * \param[in] op_
       * The \c MPI_Op handle that is to be wrapped.
       */
      explicit Operation(MPI_Op op_) : op(op_) {}
#else
    private:
      int op;
    public:
      explicit Operation(int op_) : op(op_) {}
#endif // FEAT_HAVE_MPI

    public:
      /// equality comparison operator
      bool operator==(const Operation& other) const
      {
        return this->op == other.op;
      }

      /// inequality comparison operator
      bool operator!=(const Operation& other) const
      {
        return this->op != other.op;
      }
    }; // class Operation

    /// Operation wrapper for \c MPI_SUM
    extern const Operation op_sum;
    /// Operation wrapper for \c MPI_MAX
    extern const Operation op_max;
    /// Operation wrapper for \c MPI_MIN
    extern const Operation op_min;

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    /**
     * \brief Communication Status class
     *
     * This class effectively wraps around the \c MPI_Status structure.
     *
     * \attention
     * This class is assumed to be a POD container for \c MPI_Status, i.e. you <b>must not</b>
     * add any other member variables or virtual functions to this class!
     *
     * \author Peter Zajac
     */
    class Status
    {
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
    public:
      /// the MPI status structure
      MPI_Status status;

      /// standard constructor
      Status() :
        status()
      {
        status.MPI_SOURCE = MPI_PROC_NULL;
      }

      /// status objects are copyable
      Status(const Status&) = default;
      /// status objects are copyable
      Status& operator=(const Status&) = default;

      /// \returns A pointer to the internal \c MPI_Status object.
      MPI_Status* mpi_status()
      {
        return &status;
      }

      /**
       * \brief Checks whether the source rank is \c MPI_PROC_NULL.
       *
       * \returns \c true, if the source rank is \c MPI_PROC_NULL, otherwise \c false.
       */
      bool is_null() const
      {
        return status.MPI_SOURCE == MPI_PROC_NULL;
      }

      /**
       * \returns The source rank, as stored in <c>status.MPI_SOURCE</c>
       */
      int source() const
      {
        return status.MPI_SOURCE;
      }

      /**
       * \returns The communication tag, as stored in <c>status.MPI_TAG</c>
       */
      int tag() const
      {
        return status.MPI_TAG;
      }

      /**
       * \returns The communication error status, as stored in <c>status.MPI_ERROR</c>
       */
      int error() const
      {
        return status.MPI_ERROR;
      }

      /**
       * \brief Returns the size of the message.
       *
       * This function effectively wraps around \c MPI_Get_count().
       *
       * \param[in] datatype
       * The datatype that the message is encoded in.
       *
       * \returns
       * The size of the message in instances of the corresponding datatype.
       */
      std::size_t get_count(const Datatype& datatype) const
      {
        int c(0);
        MPI_Get_count(&status, datatype.dt, &c);
        return std::size_t(c);
      }

      /// \returns The size of the message in bytes.
      std::size_t get_size() const
      {
        int c(0);
        MPI_Get_count(&status, MPI_BYTE, &c);
        return std::size_t(c);
      }
#else // non-MPI implementation
    public:
      /// status-dummy for non-MPI builds
      int status;

      Status() :
        status(0)
      {
      }
#endif // FEAT_HAVE_MPI
    }; // class Status

#ifdef FEAT_HAVE_MPI
    // Ensure that 'Status' is a POD container for 'MPI_Status',
    // as this is required by the RequestVector class.
    static_assert(sizeof(Status) == sizeof(MPI_Status), "invalid Status class size; class Status must be a POD");
#endif

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    /**
     * \brief Communication Request class
     *
     * This class effectively wraps around the \c MPI_Request handle.
     *
     * This class is move-constructible and move-assigneable, but objects of this class are non-copyable.
     *
     * \warning
     * This class \e intentionally does \b not free the internal MPI request upon destruction!\n
     * If you \e really want to free the request without waiting for it,
     * you have to do so manually by calling the #free() member function!
     *
     * \attention
     * This class is assumed to be a POD container for \c MPI_Request, i.e. you <b>must not</b>
     * add any other member variables or virtual functions to this class!
     *
     * \see \cite MPI31 Section 3.7
     *
     * \author Peter Zajac
     */
    class Request
    {
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
    public:
      /// our internal MPI request handle
      MPI_Request request;
#else
    private:
      /// request dummy for non-MPI builds
      int request;
#endif // FEAT_HAVE_MPI

    public:
      /**
       * \brief Standard constructor
       *
       * This constructor creates a (valid) null request,
       * i.e. a request which is already fulfilled.
       */
      Request();

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// MPI_Request constructor
      explicit Request(MPI_Request req_);

      /// \returns A pointer to the internal \c MPI_Request handle.
      MPI_Request* mpi_request()
      {
        return &request;
      }
#endif // FEAT_HAVE_MPI

      /// Request objects are non-copyable
      Request(const Request&) = delete;
      /// Request objects are non-copyable
      Request& operator=(const Request&) = delete;

      /**
       * \brief Move constructor
       *
       * This constructor moves the input request into the newly constructed request
       * and sets the input request \p other to a null request.
       *
       * \param[inout] other
       * The request that is to be moved.
       */
      Request(Request&& other);

      /**
       * \brief Move-assignment operator
       *
       * This operator moves the input request into this request object and sets the
       * input request \p other to a null request.
       *
       * \attention
       * The destination request represented by \p this must be a null request,
       * as this operator will fire an assertion failure otherwise!
       *
       * \param[inout] other
       * The request that is to be moved.
       *
       * \returns \p *this.
       */
      Request& operator=(Request&& other);

      /// non-virtual destructor
      ~Request()
      {
        XASSERTM(is_null(), "You are trying to destroy an active request!");
      }

      /**
       * \brief Checks whether this request is null.
       *
       * This function effectively checks whether the internal \c MPI_Request handle is equal to \c MPI_REQUEST_NULL.
       *
       * \returns \c true, if this request is null, otherwise \c false.
       *
       * \see \cite MPI31 Section 3.7.3, page 52
       */
      bool is_null() const;

      /**
       * \brief Frees the request.
       *
       * This function effectively calls \c MPI_Request_free() for the internal \c MPI_Request handle.
       *
       * \note
       * This function does nothing if the request is already null.
       *
       * \see \cite MPI31 Section 3.7.3, page 55
       */
      void free();

      /**
       * \brief Cancels the request.
       *
       * This function effectively calls \c MPI_Cancel() for the internal \c MPI_Request handle.
       *
       * \note
       * This function does nothing if the request is already null.
       *
       * \attention
       * This function does no free the request!
       *
       * \see \cite MPI31 Section 3.8.4, page 71
       */
      void cancel();

      /**
       * \brief Tests whether the request is fulfilled (or null) without blocking.
       *
       * This function effectively calls \c MPI_Test() for the internal \c MPI_Request handle.
       *
       * \note
       * For a null request, this function returns \c true, but leaves the \p status argument unmodified.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the request,
       * if it was fulfilled by this call.
       *
       * \returns \c true, if the request is fulfilled (or null), or \c false, if the request is still active.
       *
       * \see \cite MPI31 Section 3.7.3, page 54
       */
      bool test(Status& status);

      /**
       * \brief Tests whether the request is fulfilled (or null).
       *
       * This function calls #test(Status&) and discards the status object.
       *
       * \returns \c true, if the request is fulfilled (or null), or \c false, if the request is still active.
       */
      bool test()
      {
        Status status;
        return test(status);
      }

      /**
       * \brief Blocks until the request is fulfilled (or null).
       *
       * This function effectively calls \c MPI_Wait() for the internal \c MPI_Request handle.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the request,
       * if it was fulfilled by this call.
       *
       * \returns
       * \c true, if the request was fulfilled by this call, or
       * \c false, if the request was already null prior to this call.
       *
       * \see \cite MPI31 Section 3.7.3, page 53
       */
      bool wait(Status& status);

      /**
       * \brief Blocks until the request is fulfilled (or null).
       *
       * This function calls #wait(Status&) and discards the status object.
       *
       * \returns
       * \c true, if the request was fulfilled by this call, or
       * \c false, if the request was already null prior to this call.
       */
      bool wait()
      {
        Status status;
        return wait(status);
      }
    }; // class Request

#ifdef FEAT_HAVE_MPI
    // Ensure that 'Request' is a POD container for 'MPI_Request',
    // as this is required by the RequestVector class.
    static_assert(sizeof(Request) == sizeof(MPI_Request), "invalid Request class size; class Request must be a POD");
#endif

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    /**
     * \brief Communication Request vector class
     *
     * This class encapsulates vectors of Request and Status objects, which can be used to manage
     * a (related) set of requests. This class also offers the corresponding wait and test
     * wrapper functions to check for the completition of a single, multiple or all requests.
     *
     * This class is move-constructible and move-assigneable, but objects of this class are non-copyable.
     *
     * <b>Waiting for completition:</b>\n
     * There are two common scenarios regarding the completition of a request vector:
     * - Wait until all requests are fulfilled: this can be easily achieved by calling #wait_all()
     * - Wait until one request is fulfilled, then do some work that is related to the completition
     *   of this particular request, and repeat these two steps until all requests have been processed.
     *
     * The second of the above scenarios can be easily implemented by using a \c for or \c while loop
     * and the #wait_any() function:
       \code{.cpp}
       RequestVector reqs;

       ... // create some requests here

       for(std::size_t idx(0u); reqs.wait_any(idx); )
       {
         ... // request 'idx' has been fulfilled
       }

       // all requests have been processed
       \endcode
     *
     * \author Peter Zajac
     */
    class RequestVector
    {
    private:
      /// internal vector of Request objects
      std::vector<Request> _reqs;
      /// internal vector of Status objects
      std::vector<Status> _stats;

#if defined(FEAT_HAVE_MPI)
      int _isize() const
      {
        return int(_reqs.size());
      }

      MPI_Request* _reqs_array()
      {
        return reinterpret_cast<MPI_Request*>(_reqs.data());
      }

      MPI_Status* _stats_array()
      {
        return reinterpret_cast<MPI_Status*>(_stats.data());
      }
#endif // FEAT_HAVE_MPI

    public:
      /**
       * \brief Standard constructor
       *
       * This constructor creates an empty request vector.
       */
      RequestVector() :
        _reqs(),
        _stats()
      {
      }

      /**
       * \brief Sized constructor
       *
       * This constructor creates a request vector of the desired size consisting
       * of null requests and empty statuses.
       *
       * \param[in] size_
       * The number of null requests (and status objects) to allocate.
       */
      explicit RequestVector(std::size_t size_) :
        _reqs(size_),
        _stats(size_)
      {
      }

      /// request vectors are non-copyable
      RequestVector(const RequestVector&) = delete;
      /// request vectors are non-copyable
      RequestVector& operator=(const RequestVector&) = delete;

      /**
       * \brief Move constructor
       *
       * This constructor moves the input request vector into the newly constructed
       * request vector and clears the input request vector \p other.
       *
       * \param[inout] other
       * The source request vector that is to be moved.
       */
      RequestVector(RequestVector&& other) :
        _reqs(std::forward<std::vector<Request>>(other._reqs)),
        _stats(std::forward<std::vector<Status>>(other._stats))
      {
        other.clear();
      }

      /**
       * \brief Move-assignment operator
       *
       * This operator moves the input request vector into this request vector object
       * and clears the input request vector \p other.
       *
       * \attention
       * The destination request vector represented by \p this must be a null request,
       * i.e. it must not contain any active requests, as this operator will fire
       * an assertion failure otherwise!
       *
       * \param[inout] other
       * The source request vector that is to be moved.
       *
       * \returns \p *this.
       */
      RequestVector& operator=(RequestVector&& other)
      {
        if(this != &other)
        {
          XASSERT(is_null());
          _reqs = std::forward<std::vector<Request>>(other._reqs);
          _stats = std::forward<std::vector<Status>>(other._stats);
          other.clear();
        }
        return *this;
      }

      /// virtual destructor
      virtual ~RequestVector()
      {
        XASSERTM(is_null(), "You are trying to destroy an active request vector!");
      }

      /**
       * \brief Returns a (const) reference to a single request in the vector.
       *
       * \param[in] idx
       * The index of the request to be returned.
       *
       * \returns
       * A (const) reference to the desired Request object.
       */
      Request& get_request(std::size_t idx)
      {
        XASSERT(idx < _reqs.size());
        return _reqs[idx];
      }

      /** \copydoc get_request(std::size_t) */
      const Request& get_request(std::size_t idx) const
      {
        XASSERT(idx < _reqs.size());
        return _reqs[idx];
      }

      /**
       * \brief Returns a (const) reference to a single status in the vector.
       *
       * \param[in] idx
       * The index of the status to be returned.
       *
       * \returns
       * A (const) reference to the desired Status object.
       */
      Status& get_status(std::size_t idx)
      {
        XASSERT(idx < _stats.size());
        return _stats[idx];
      }

      /** \copydoc get_status(std::size_t) */
      const Status& get_status(std::size_t idx) const
      {
        XASSERT(idx < _stats.size());
        return _stats[idx];
      }

      /** \copydoc get_request(std::size_t) */
      Request& operator[](std::size_t idx)
      {
        XASSERT(idx < _reqs.size());
        return _reqs[idx];
      }

      /** \copydoc get_request(std::size_t) */
      const Request& operator[](std::size_t idx) const
      {
        XASSERT(idx < _reqs.size());
        return _reqs[idx];
      }

      /**
       * \brief Reserves sufficient space for a specified number of requests.
       *
       * This function effectively wraps around the <c>reserve()</c> function
       * of the internal <c>std::vector</c> objects.
       *
       * \param[in] size_
       * The minimum number of requests that space is to be reserved for.
       */
      void reserve(std::size_t size_)
      {
        if(size_ <= _reqs.size())
          return;
        _reqs.reserve(size_);
        _stats.reserve(size_);
      }

      /**
       * \brief Resizes the request vector.
       *
       * This function effectively wraps around the <c>resize()</c> function
       * of the internal <c>std::vector</c> objects.
       *
       * \attention
       * This function does not free requests that would be removed if the size of the request vector is
       * decreased! Instead, this function fires an assertion failure if a non-null request would be deleted.
       *
       * \param[in] size_
       * The new number of requests to be managed.
       */
      void resize(std::size_t size_)
      {
        // make sure that all to-be-removed requests are null
        for(std::size_t i(size_); i < _reqs.size(); ++i)
        {
          XASSERT(_reqs.at(i).is_null());
        }
        _reqs.resize(size_);
        _stats.resize(size_);
      }

      /**
       * \brief Compresses the request vector.
       *
       * This function removes all null requests (and their corresponding statuses) from the vector
       * and resizes the vector to the size of the remaining active requests.
       *
       * \note
       * You usually do not need this function, as the wait/test functions can perfectly work with
       * a vector of mixed null and active requests.
       *
       * \returns
       * The number of active requests remaining in the vector.
       */
      std::size_t compress()
      {
        std::size_t i(0), j(0), n = _reqs.size();
        while(i+j < n)
        {
          if(_reqs[i+j].is_null())
            ++j;
          else
          {
            _reqs[i] = std::move(_reqs[i+j]);
            _stats[i] = std::move(_stats[i+j]);
            ++i;
          }
        }
        _reqs.resize(n-j);
        _stats.resize(n-j);
        return n-j;
      }

      /**
       * \brief Inserts a new request at the end of the request vector.
       *
       * This function also reserves a new internal Status object for the request.
       *
       * \param[inout] request
       * The request that is to be inserted into the vector.
       * This object is set to a null request upon exit of this function.
       */
      void push_back(Request&& request)
      {
        _reqs.emplace_back(std::forward<Request>(request));
        _stats.emplace_back(Status());
      }

      /**
       * \brief Returns the total number of (both active and null) requests in the vector.
       *
       * \note
       * If you want to query the number of active (i.e. non-null) requests in the vector,
       * use the #active_count() function instead.
       *
       * \returns The total number of requests in the vector.
       */
      std::size_t size() const
      {
        return _reqs.size();
      }

      /**
       * \brief Returns the number of active requests in the vector.
       *
       * \note
       * If you want to query the total number of (active and null) requests in the vector,
       * use the #size() function instead.
       *
       * \returns The number of active requests in the vector.
       */
      std::size_t active_count() const
      {
        std::size_t count(0);
        for(const auto& r : _reqs)
        {
          if(!r.is_null())
            ++count;
        }
        return count;
      }

      /**
       * \brief Checks whether the request vector is empty.
       *
       * This function checks whether there are neither active nor null requests in the vector.
       *
       * \note
       * If you want to check whether there are no active requests in the vector,
       * use the #is_null() function instead.
       *
       * \returns \c true, if there are no requests in the vector, otherwise \c false.
       */
      bool empty() const
      {
        return _reqs.empty();
      }

      /**
       * \brief Checks whether all requests are null requests.
       *
       * This function checks whether are no active requests left in the vector.
       *
       * \note
       * If you want to check whether there are neither active nor null requests in the vector,
       * use the #empty() function instead.
       *
       * \returns \c true, if all requests are null, or
       * \c false, if there exists at least one active request in the vector.
       */
      bool is_null() const
      {
        for(const auto& r : _reqs)
        {
          if(!r.is_null())
            return false;
        }
        return true;
      }

      /**
       * \brief Frees all remaining active requests.
       *
       * This function frees all active requests, but does not resize or clear the vector.
       * Upon exit, the vector consists only of null requests (or no requests at all).
       *
       * \note
       * This function does nothing if there are no requests at all or if all requests are
       * already null requests.
       */
      void free()
      {
        for(auto& r : _reqs)
        {
          if(!r.is_null())
            r.free();
        }
      }

      /**
       * \brief Cancels all remaining active requests.
       *
       * This function cancels all active requests, but does not free them.
       *
       * \note
       * This function does nothing if there are no requests at all or if all requests are
       * already null requests.
       */
      void cancel()
      {
        for(auto& r : _reqs)
        {
          if(!r.is_null())
            r.cancel();
        }
      }

      /**
       * \brief Clears the request vector.
       *
       * This function sets the size of the request vector to zero.
       *
       * \attention
       * This function must not be called if there are still active requests in vector,
       * as this function will fire an assertion failure otherwise!
       *
       * \note
       * - If you want to remove all null requests from the vector, use the #compress() function instead.
       * - If you want to free all active requests in the vector, use the #free() function instead.
       * - If you want to cancel all active requests in the vector, use the #cancel() function instead.
       */
      void clear()
      {
        XASSERT(is_null());
        _reqs.clear();
        _stats.clear();
      }

      /**
       * \brief Tests whether a specific request is fulfilled (or null).
       *
       * \param[in] idx
       * The index of the request that is to be tested.
       *
       * \note
       * If this function returns \c true, the corresponding Status object can be queried by #get_status().
       *
       * \returns \c true, if the corresponding request is fulfilled, or \c false, if the request is still active.
       */
      bool test_for(std::size_t idx)
      {
        return get_request(idx).test(get_status(idx));
      }

      /**
       * \brief Tests whether a specific request is fulfilled (or null).
       *
       * \param[in] idx
       * The index of the request that is to be tested.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the request,
       * if it was fulfilled by this call.
       *
       * \attention
       * This function does \b not store the status of a successful test in the internal
       * status array of this RequestVector object for efficiency reasons!
       *
       * \returns \c true, if the corresponding request is fulfilled, or \c false, if the request is still active.
       */
      bool test_for(std::size_t idx, Status& status)
      {
        return get_request(idx).test(status);
      }

      /**
       * \brief Tests whether all active requests are fulfilled (or null).
       *
       * \note
       * If you want to check whether there are remaining active requests without changing the
       * state (i.e. fulfilling previously active requests), use the #is_null() function instead.
       *
       * \returns
       * \c true, if all requests are fulfilled, of \c false, if there at least one active request remains.
       *
       * \see \cite MPI31 Section 3.7.5, page 60
       */
      bool test_all();

      /**
       * \brief Tests whether one of the active requests has been fulfilled.
       *
       * \param[out] idx
       * The index of the previously active request that has been fulfilled by this call.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the active request
       * that has been fulfilled by this call.
       *
       * \attention
       * This function does \b not store the status of a successful test in the internal
       * status array of this RequestVector object for efficiency reasons!
       *
       * \returns
       * \c true, if a previously active request has been fulfilled by this call, or \c false,
       * if no active request has been fulfilled or if there were no active requests left
       * prior to the call of this function.
       *
       * \see \cite MPI31 Section 3.7.5, page 58
       */
      bool test_any(std::size_t& idx, Status& status);

      /**
       * \brief Tests whether one of the active requests has been fulfilled.
       *
       * \param[out] idx
       * The index of the previously active request that has been fulfilled by this call.
       *
       * \note
       * If this function returns \c true, the corresponding Status object can be queried by #get_status().
       *
       * \returns
       * \c true, if a previously active request has been fulfilled by this call, or \c false,
       * if no active request has been fulfilled or if there were no active requests left
       * prior to the call of this function.
       */
      bool test_any(std::size_t& idx)
      {
        Status status;
        if(!test_any(idx, status))
          return false;
        _stats[idx] = status;
        return true;
      }

      /**
       * \brief Tests whether one of the active requests has been fulfilled.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the active request
       * that has been fulfilled by this call.
       *
       * \attention
       * This function does \b not store the status of a successful test in the internal
       * status array of this RequestVector object for efficiency reasons!
       *
       * \returns
       * \c true, if a previously active request has been fulfilled by this call, or \c false,
       * if no active request has been fulfilled or if there were no active requests left
       * prior to the call of this function.
       */
      bool test_any(Status& status)
      {
        std::size_t idx;
        return test_any(idx, status);
      }

      /**
       * \brief Blocks until a specific request is fulfilled (or null).
       *
       * \param[in] idx
       * The index of the request that is to be waited for.
       *
       * \note
       * The corresponding Status object can be queried by #get_status().
       *
       * \returns
       * \c true, if the corresponding request was fulfilled by this call, or
       * \c false, if the corresponding request was already null prior to this call.
       */
      bool wait_for(std::size_t idx)
      {
        return get_request(idx).wait(get_status(idx));
      }

      /**
       * \brief Blocks until a specific request is fulfilled (or null).
       *
       * \param[in] idx
       * The index of the request that is to be waited for.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the request,
       * if it was fulfilled by this call.
       *
       * \attention
       * This function does \b not store the status of a successful wait in the internal
       * status array of this RequestVector object for efficiency reasons!
       *
       * \returns
       * \c true, if the corresponding request was fulfilled by this call, or
       * \c false, if the corresponding request was already null prior to this call.
       */
      bool wait_for(std::size_t idx, Status& status)
      {
        return get_request(idx).wait(status);
      }

      /**
       * \brief Blocks until all active requests are fulfilled.
       *
       * \note
       * This function does nothing if all requests are already null upon call.
       *
       * \see \cite MPI31 Section 3.7.5, page 59
       */
      void wait_all();

      /**
       * \brief Blocks until one of the active requests has been fulfilled.
       *
       * \param[out] idx
       * The index of the previously active request that has been fulfilled by this call.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the active request
       * that has been fulfilled by this call.
       *
       * \attention
       * This function does \b not store the status of a successful test in the internal
       * status array of this RequestVector object for efficiency reasons!
       *
       * \returns
       * \c true, if a previously active request has been fulfilled by this call, or \c false,
       * if there were no active requests left upon call of this function.
       *
       * \note
       * The return value of this function can be used as a break-condition in a \c while loop
       * as this function only returns \c false when all active requests have been fulfilled.
       *
       * \see \cite MPI31 Section 3.7.5, page 57
       */
      bool wait_any(std::size_t& idx, Status& status);

      /**
       * \brief Blocks until one of the active requests has been fulfilled.
       *
       * \param[out] idx
       * The index of the previously active request that has been fulfilled by this call.
       *
       * \note
       * If this function returns \c true, the corresponding Status object can be queried by #get_status().
       *
       * \returns
       * \c true, if a previously active request has been fulfilled by this call, or \c false,
       * if there were no active requests left upon call of this function.
       *
       * \note
       * The return value of this function can be used as a break-condition in a \c while loop
       * as this function only returns \c false when all active requests have been fulfilled.
       */
      bool wait_any(std::size_t& idx)
      {
        Status status;
        if(!wait_any(idx, status))
          return false;
        _stats[idx] = status;
        return true;
      }

      /**
       * \brief Blocks until one of the active requests has been fulfilled.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the active request
       * that has been fulfilled by this call.
       *
       * \attention
       * This function does \b not store the status of a successful test in the internal
       * status array of this RequestVector object for efficiency reasons!
       *
       * \returns
       * \c true, if a previously active request has been fulfilled by this call, or \c false,
       * if there were no active requests left upon call of this function.
       *
       * \note
       * The return value of this function can be used as a break-condition in a \c while loop
       * as this function only returns \c false when all active requests have been fulfilled.
       */
      bool wait_any(Status& status)
      {
        std::size_t idx;
        return wait_any(idx, status);
      }
    }; // class RequestVector

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    /**
     * \brief Communicator class
     *
     * This class effectively wraps around the \c MPI_Comm handle and contains wrapper functions
     * for many (but not all) MPI functions, which work on communicators, most notably all
     * message passing functions.
     *
     * \note
     * The documentation of this class is related to parallel MPI builds. For non-MPI builds, this
     * class offers corresponding serial implementations, which behave just as if they were called
     * for a MPI build with a world consisting of exactly 1 process. In many cases, the corresponding
     * serial version is an empty function, whereas in some other cases, a memcpy between the source
     * and receive buffers has to be performed.
     *
     * The wrapper functions provided by this class are written in accordance with the following guidelines:
     * - All function names correspond to their MPI function names in lowercase without the \c MPI_ prefix.
     * - All functions, which do not return a Request object, are of return type \c void.
     * - Function arguments have the same semantics and the same order as their MPI counterparts.
     *   In almost all cases, the function argument names also coincide.
     * - The \p comm argument is removed, as it is given by this Comm object itself.
     * - All size/count arguments are of type <c>std::size_t</c> rather than \c int.\n
     *   <b>Exception:</b> Arguments that are arrays of sizes/counts/displacements, are of type \c int
     *   for technical reasons, e.g. #allgatherv or #alltoallv.
     * - All \c MPI_Datatype handle arguments are replaced by references to Dist::Datatype objects.
     * - All \c MPI_Op handle arguments are replaced by references to Dist::Operation objects.
     * - All \c MPI_Status pointer arguments are replaced by references to Dist::Status objects.
     * - All \e output \c MPI_Request pointer arguments are removed, as nonblocking operations return
     *   a Dist::Request object representing the request.
     * - All \c tag arguments are defaulted to 0, unless there exists another non-defaulted argument
     *   following the \c tag argument, of course.
     *
     * Furthermore, this class offers additional overloads for many of the "pure" wrapper functions:
     * - For each function, which has at least one Datatype argument, there exists an overloaded
     *   template wrapper function, which automatically deducts the datatype(s) by using Dist::autotype.
     *   Of course, these overloads work only for arguments of fundamental C/C++ datatypes.
     *   For these overloads, the datatype arguments are removed.
     *
     * <b>Examples:</b>\n
     * 1. The blocking receive function has the following MPI prototype:
     *    \code{.cpp}
     *    int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
     *    \endcode
     *    The corresponding wrapper function for is the Dist::Comm::recv() member function:
     *    \code{.cpp}
     *    void recv(void* buf, std::size_t count, const Datatype& datatype, int source, int tag, Status& status) const
     *    \endcode
     *    The overloaded wrapper function template with automatic type deduction is:
     *    \code{.cpp}
     *    template<typename T_>
     *    void recv(T_* buffer, std::size_t count, int source, int tag, Status& status) const
     *    \endcode
     *
     * 2. The nonblocking reduce function has the following MPI prototype:
     *    \code{.cpp}
     *    int MPI_Ireduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
     *    \endcode
     *    The corresponding wrapper function for is the Dist::Comm::ireduce() member function:
     *    \code{.cpp}
     *    Request ireduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op, int root) const
     *    \endcode
     *    The overloaded wrapper function template with automatic type deduction is:
     *    \code{.cpp}
     *    template<typename T_>
     *    Request ireduce(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op, int root) const
     *    \endcode
     *
     * \author Peter Zajac
     */
    class Comm
    {
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
    public:
      /// our MPI communicator handle
      MPI_Comm comm;
#endif // FEAT_HAVE_MPI

    protected:
      /// our rank
      int _rank;
      /// the communicator size
      int _size;

#if !defined(FEAT_HAVE_MPI) && !defined(DOXYGEN)
      explicit Comm(int);
#endif

    public:
      /**
       * \brief Standard constructor
       *
       * This constructor creates a null communicator representing \c MPI_COMM_NULL.
       */
      Comm();

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /**
       * \brief \c MPI_Comm handle constructor
       *
       * This constructor creates a communicator for a given \c MPI_Comm handle.
       *
       * \param[in] comm_
       * The \c MPI_Comm handle that is to be encapsulated.
       *
       * \note
       * This class will automatically free the handle upon destruction by calling
       * \c MPI_Comm_free() unless \p comm_ is \c MPI_COMM_WORLD or \c MPI_COMM_NULL.
       */
      explicit Comm(MPI_Comm comm_);
#endif // FEAT_HAVE_MPI

      /// communicators are non-copyable
      Comm(const Comm&) = delete;
      /// communicators are non-copyable
      Comm& operator=(const Comm&) = delete;

      /**
       * \brief Move constructor
       *
       * This constructor moves the internal handle of the source communicator and
       * sets the source communicator to a null communicator.
       *
       * \note
       * If the source communicator is the world communicator, this constructor
       * will create a copy of the world communicator without modifying \p other.
       *
       * \param[inout] other
       * The source communicator that is to be moved.
       */
      Comm(Comm&& other);

      /**
       * \brief Move-assignment operator
       *
       * This operator moves the internal handle of the source communicator and
       * sets the source communicator to a null communicator.
       *
       * \attention
       * The destination communicator represented by \p this must be a null communicator,
       * as this operator will fire an assertion failure otherwise!
       *
       * \note
       * If the source communicator is the world communicator, this constructor
       * will create a copy of the world communicator without modifying \p other.
       *
       * \param[inout] other
       * The source communicator that is to be moved.
       */
      Comm& operator=(Comm&& other);

      /**
       * \brief virtual destructor
       *
       * This destructor frees the internal \c MPI_Comm handle by calling
       * \c MPI_Comm_free() unless it is \c MPI_COMM_WORLD or \c MPI_COMM_NULL.
       */
      virtual ~Comm();

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// \returns a const reference to the internal \c MPI_Comm handle
      const MPI_Comm& mpi_comm() const
      {
        return comm;
      }

      /// \returns a reference to the internal \c MPI_Comm handle
      MPI_Comm& mpi_comm()
      {
        return comm;
      }
#endif //FEAT_HAVE_MPI

      /**
       * \brief Returns a copy of the world communicator.
       */
      static Comm world();

      /**
       * \brief Returns a copy of the self communicator.
       */
      static Comm self();

      /**
       * \brief Returns a null communicator.
       */
      static Comm null();

      /**
       * \brief Checks whether this communicator is the world communicator.
       *
       * \returns \c true, if this communicator represents \c MPI_COMM_WORLD, otherwise \c false.
       */
      bool is_world() const;

      /**
       * \brief Checks whether this communicator is the self communicator.
       *
       * \returns \c true, if this communicator represents \c MPI_COMM_SELF, otherwise \c false.
       */
      bool is_self() const;

      /**
       * \brief Checks whether this communicator is a null communicator.
       *
       * \returns \c true, if this communicator represents \c MPI_COMM_NULL, otherwise \c false.
       */
      bool is_null() const;

      /**
       * \brief Returns the rank of this process in this communicator.
       *
       * \returns The rank of this process in this communicator.
       *
       * \see \cite MPI31 Section 6.4.1, page 236
       */
      int rank() const
      {
        return _rank;
      }

      /**
       * \brief Returns the size of this communicator.
       *
       * \returns The size of this communicator.
       *
       * \see \cite MPI31 Section 6.4.1, page 235
       */
      int size() const
      {
        return _size;
      }

      /**
       * \name Comm Creation
       */
      ///@{

      /**
       * \brief Creates a copy of this communicator.
       *
       * \see \cite MPI31, Section 6.4.2, page 238
       *
       * \returns A copy of this communicator.
       */
      Comm comm_dup() const;

      /**
       * \brief Creates a new sub-communicator from a strided range of ranks.
       *
       * This function is a short-cut, which creates a process group by using
       * \c MPI_Group_range_incl and then creates a corresponding new communicator
       * via \c MPI_Comm_create.\n
       * The ranks, which are contained in the new sub-communicator, are given
       * by <em>first, first+stride, first+2*stride, ..., first+(count-1)*stride</em>.
       *
       * \see \cite MPI31, Section 6.3.2, page 233
       * \see \cite MPI31, Section 6.4.2, page 240
       *
       * \param[in] count
       * The number desired ranks in the sub-communicator. Must be > 0.
       *
       * \param[in] first
       * The rank of the first process to be included in the sub-communicator.
       *
       * \param[in] stride
       * The stride for the rank range. Must be >= 1.
       *
       * \returns
       * A new sub-communicator for the range of processes. If this process is not
       * part of the new sub-communicator, the returned comm is a null communicator.
       */
      Comm comm_create_range_incl(int count, int first = 0, int stride = 1) const;

      /**
       * \brief Creates a new sub-communicator for a given set of ranks.
       *
       * This function is a short-cut, which creates a process group by using
       * \c MPI_Group_incl and then creates a corresponding new communicator
       * via \c MPI_Comm_create.\n
       *
       * \see \cite MPI31, Section 6.3.2, page 232
       * \see \cite MPI31, Section 6.4.2, page 240
       *
       * \param[in] n
       * The number of desired ranks in the sub-communicator. Must be > 0.
       *
       * \param[in] ranks
       * An array of unique ranks that should be part of the new sub-communicator.
       *
       * \returns
       * A new sub-communicator for the set of processes. If this process is not
       * part of the new sub-communicator, the returned comm is a null communicator.
       */
      Comm comm_create_incl(int n, const int* ranks) const;

      /**
       * \brief Creates a new sub-communicator by splitting this communicator.
       *
       * This functions splits this communicator into disjoint sub-communicators.
       *
       * \param[in] color
       * The color of this process.
       *
       * \param[in] key
       * The key for the ranking of the process.
       *
       * \see \cite MPI31, Section 6.4.2, page 244
       *
       * \returns
       * A new communicator for the set of processes of the same color.
       */
      Comm comm_split(int color, int key) const;

      ///@}

      /**
       * \name Barrier Synchronization
       */
      ///@{

      /**
       * \brief Blocking barrier
       *
       * \see \cite MPI31 Section 5.3, page 147
       */
      void barrier() const;

      /**
       * \brief Nonblocking barrier
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.1, page 198
       */
      Request ibarrier() const;

      // end of barrier synchronization group
      ///@}

      /**
       * \name Point-To-Point Communication (Blocking)
       */
      ///@{

      /**
       * \brief Blocking Send
       *
       * \param[in] buffer
       * The send buffer for the operation.
       *
       * \param[in] count
       * The size of the send buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[in] dest
       * The rank of the destination process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \see \cite MPI31 Section 3.2.1, page 24
       */
      void send(const void* buffer, std::size_t count, const Datatype& datatype, int dest, int tag = 0) const;

      /**
       * \brief Blocking Send
       *
       * This function automatically deducts the datatype of the send buffer (if possible).
       *
       * \param[in] buffer
       * The send buffer for the operation.
       *
       * \param[in] count
       * The size of the send buffer in datatype objects.
       *
       * \param[in] dest
       * The rank of the destination process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \see \cite MPI31 Section 3.2.1, page 24
       */
      template<typename T_>
      void send(const T_* buffer, std::size_t count, int dest, int tag = 0) const
      {
        send(buffer, count, autotype<T_>(), dest, tag);
      }

      /**
       * \brief Blocking Receive
       *
       * \param[in] buffer
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \param[in] source
       * The rank of the source process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the receive operation.
       *
       * \see \cite MPI31 Section 3.2.4, page 28
       */
      void recv(void* buffer, std::size_t count, const Datatype& datatype, int source, int tag, Status& status) const;

      /**
       * \brief Blocking Receive
       *
       * This function discards the Status object of the receive operation.
       *
       * \param[in] buffer
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \param[in] source
       * The rank of the source process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \see \cite MPI31 Section 3.2.4, page 28
       */
      void recv(void* buffer, std::size_t count, const Datatype& datatype, int source, int tag = 0) const
      {
        Status status;
        recv(buffer, count, datatype, source, tag, status);
      }

      /**
       * \brief Blocking Receive
       *
       * This function automatically deducts the datatype of the receive buffer (if possible).
       *
       * \param[in] buffer
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] source
       * The rank of the source process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \param[out] status
       * A reference to the Status object that receives information about the receive operation.
       *
       * \see \cite MPI31 Section 3.2.4, page 28
       */
      template<typename T_>
      void recv(T_* buffer, std::size_t count, int source, int tag, Status& status) const
      {
        recv(buffer, count, autotype<T_>(), source, tag, status);
      }

      /**
       * \brief Blocking Receive
       *
       * This function automatically deducts the datatype of the receive buffer (if possible).
       * This function discards the Status object of the receive operation.
       *
       * \param[in] buffer
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] source
       * The rank of the source process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \see \cite MPI31 Section 3.2.4, page 28
       */
      template<typename T_>
      void recv(T_* buffer, std::size_t count, int source, int tag = 0) const
      {
        Status status;
        recv(buffer, count, autotype<T_>(), source, tag, status);
      }

      // end of blocking point-to-point group
      ///@}

      /**
       * \name Point-To-Point Communication (Nonblocking)
       */
      ///@{

      /**
       * \brief Nonblocking Send
       *
       * \param[in] buffer
       * The send buffer for the operation.
       *
       * \param[in] count
       * The size of the send buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[in] dest
       * The rank of the destination process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 3.7.2, page 49
       */
      Request isend(const void* buffer, std::size_t count, const Datatype& datatype, int dest, int tag = 0) const;

      /**
       * \brief Nonblocking Send
       *
       * This function automatically deducts the datatype of the send buffer (if possible).
       *
       * \param[in] buffer
       * The send buffer for the operation.
       *
       * \param[in] count
       * The size of the send buffer in datatype objects.
       *
       * \param[in] dest
       * The rank of the destination process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 3.7.2, page 49
       */
      template<typename T_>
      Request isend(const T_* buffer, std::size_t count, int dest, int tag = 0) const
      {
        return isend(buffer, count, autotype<T_>(), dest, tag);
      }

      /**
       * \brief Nonblocking Receive
       *
       * \param[in] buffer
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \param[in] source
       * The rank of the source process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 3.7.2, page 51
       */
      Request irecv(void* buffer, std::size_t count, const Datatype& datatype, int source, int tag = 0) const;

      /**
       * \brief Nonblocking Receive
       *
       * This function automatically deducts the datatype of the receive buffer (if possible).
       *
       * \param[in] buffer
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] source
       * The rank of the source process.
       *
       * \param[in] tag
       * The tag for the message.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 3.7.2, page 51
       */
      template<typename T_>
      Request irecv(T_* buffer, std::size_t count, int source, int tag = 0) const
      {
        return irecv(buffer, count, autotype<T_>(), source, tag);
      }

      // end of nonblocking point-to-point group
      ///@}

      /**
       * \name Broadcasts
       */
      ///@{

      /**
       * \brief Blocking broadcast
       *
       * \param[inout] buffer
       * The send/receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send/receive buffer contents.
       *
       * \param[in] root
       * The rank of the root process that is sending the broadcast.
       *
       * \see \cite MPI31 Section 5.4, page 148
       */
      void bcast(void* buffer, std::size_t count, const Datatype& datatype, int root) const;

      /**
       * \brief Blocking broadcast
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[inout] buffer
       * The send/receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] root
       * The rank of the root process that is sending the broadcast.
       *
       * \see \cite MPI31 Section 5.4, page 148
       */
      template<typename T_>
      void bcast(T_* buffer, std::size_t count, int root) const
      {
        bcast(buffer, count, autotype<T_>(), root);
      }

      /**
       * \brief Nonblocking broadcast
       *
       * \param[inout] buffer
       * The send/receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send/receive buffer contents.
       *
       * \param[in] root
       * The rank of the root process that is sending the broadcast.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.2, page 199
       */
      Request ibcast(void* buffer, std::size_t count, const Datatype& datatype, int root) const;

      /**
       * \brief Nonblocking broadcast
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[inout] buffer
       * The send/receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] root
       * The rank of the root process that is sending the broadcast.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.2, page 199
       */
      template<typename T_>
      Request ibcast(T_* buffer, std::size_t count, int root) const
      {
        return ibcast(buffer, count, autotype<T_>(), root);
      }

      // end of broadcast group
      ///@}

      /**
       * \name Gather and Scatter
       */
      ///@{

      /**
       * \brief Blocking gather
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \param[in] root
       * The rank of the root process that is gathering the data.
       *
       * \see \cite MPI31 Section 5.5, page 149
       */
      void gather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const;

      /**
       * \brief Blocking gather
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \param[in] root
       * The rank of the root process that is gathering the data.
       *
       * \see \cite MPI31 Section 5.5, page 149
       */
      template<typename ST_, typename RT_>
      void gather(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, std::size_t recvcount, int root) const
      {
        gather(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcount, autotype<RT_>(), root);
      }

      /**
       * \brief Nonblocking gather
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \param[in] root
       * The rank of the root process that is gathering the data.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.3, page 200
       */
      Request igather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const;

      /**
       * \brief Nonblocking gather
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \param[in] root
       * The rank of the root process that is gathering the data.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.3, page 200
       */
      template<typename ST_, typename RT_>
      Request igather(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, std::size_t recvcount, int root) const
      {
        return igather(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcount, autotype<RT_>(), root);
      }

      /**
       * \brief Blocking scatter
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The number of datatype objects send to \e each process.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \param[in] root
       * The rank of the root process that is scattering the data.
       *
       * \see \cite MPI31 Section 5.6, page 159
       */
      void scatter(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const;

      /**
       * \brief Blocking scatter
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The number of datatype objects send to \e each process.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] root
       * The rank of the root process that is scattering the data.
       *
       * \see \cite MPI31 Section 5.6, page 159
       */
      template<typename ST_, typename RT_>
      void scatter(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, std::size_t recvcount, int root) const
      {
        scatter(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcount, autotype<RT_>(), root);
      }

      /**
       * \brief Nonblocking scatter
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The number of datatype objects send to \e each process.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \param[in] root
       * The rank of the root process that is scattering the data.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.4, page 202
       */
      Request iscatter(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const;

      /**
       * \brief Nonblocking scatter
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The number of datatype objects send to \e each process.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] root
       * The rank of the root process that is scattering the data.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.4, page 202
       */
      template<typename ST_, typename RT_>
      Request iscatter(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, std::size_t recvcount, int root) const
      {
        return iscatter(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcount, autotype<RT_>(), root);
      }

      /**
       * \brief Blocking gather-to-all
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \see \cite MPI31 Section 5.7, page 165
       */
      void allgather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const;

      /**
       * \brief Blocking gather-to-all
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \see \cite MPI31 Section 5.7, page 165
       */
      template<typename ST_, typename RT_>
      void allgather(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, std::size_t recvcount) const
      {
        allgather(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcount, autotype<RT_>());
      }

      /**
       * \brief Nonblocking gather-to-all
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.5, page 204
       */
      Request iallgather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const;

      /**
       * \brief Nonblocking gather-to-all
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.5, page 204
       */
      template<typename ST_, typename RT_>
      Request iallgather(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, std::size_t recvcount) const
      {
        return iallgather(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcount, autotype<RT_>());
      }

      /**
       * \brief Blocking gather-to-all
       *
       * \attention
       * In contrast to most other functions, this function uses \c int rather than <c>std::size_t</c>
       * as the type for the receive counts and displacements for technical reasons!
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcounts
       * The sizes of the receive buffers in datatype objects.
       *
       * \param[in] displs
       * The displacements of the receive buffers in datatype objects.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \see \cite MPI31 Section 5.7, page 166
       */
      void allgatherv(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, const int* recvcounts, const int* displs, const Datatype& recvtype) const;

      /**
       * \brief Blocking gather-to-all
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \attention
       * In contrast to most other functions, this function uses \c int rather than <c>std::size_t</c>
       * as the type for the receive counts and displacements for technical reasons!
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The size of the send buffer in datatype objects.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcounts
       * The sizes of the receive buffers in datatype objects.
       *
       * \param[in] displs
       * The displacements of the receive buffers in datatype objects.
       *
       * \see \cite MPI31 Section 5.7, page 166
       */
      template<typename ST_, typename RT_>
      void allgatherv(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, const int* recvcounts, const int* displs) const
      {
        allgatherv(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcounts, displs, autotype<RT_>());
      }

      /**
       * \brief Blocking All-to-All Scatter/Gather
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The number of datatype objects send to \e each process.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The size of the receive buffer in datatype objects.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \see \cite MPI31 Section 5.8, page 168
       */
      void alltoall(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const;

      /**
       * \brief Blocking All-to-All Scatter/Gather
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The number of datatype objects send to \e each process.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \see \cite MPI31 Section 5.8, page 168
       */
      template<typename ST_, typename RT_>
      void alltoall(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, std::size_t recvcount) const
      {
        alltoall(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcount, autotype<RT_>());
      }

      /**
       * \brief Nonblocking All-to-All Scatter/Gather
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The number of datatype objects send to \e each process.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.6, page 206
       */
      Request ialltoall(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const;

      /**
       * \brief Nonblocking All-to-All Scatter/Gather
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcount
       * The number of datatype objects send to \e each process.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcount
       * The number of datatype objects received from \e each process.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.6, page 206
       */
      template<typename ST_, typename RT_>
      Request ialltoall(const ST_* sendbuf, std::size_t sendcount, RT_* recvbuf, std::size_t recvcount) const
      {
        return ialltoall(sendbuf, sendcount, autotype<ST_>(), recvbuf, recvcount, autotype<RT_>());
      }

      /**
       * \brief Blocking All-to-All Scatter/Gather
       *
       * \attention
       * In contrast to most other functions, this function uses \c int rather than <c>std::size_t</c>
       * as the type for the send/receive counts and displacements for technical reasons!
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcounts
       * The number of datatype objects send to \e each process.
       *
       * \param[in] sdispls
       * The displacements of the send buffers in datatype objects.
       *
       * \param[in] sendtype
       * A reference to the Datatype object representing the send buffer contents.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcounts
       * The number of datatype objects received from \e each process.
       *
       * \param[in] rdispls
       * The displacements of the receive buffers in datatype objects.
       *
       * \param[in] recvtype
       * A reference to the Datatype object representing the receive buffer contents.
       *
       * \see \cite MPI31 Section 5.8, page 170
       */
      void alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, const Datatype& sendtype, void* recvbuf, const int* recvcounts, const int* rdispls, const Datatype& recvtype) const;

      /**
       * \brief Blocking All-to-All Scatter/Gather
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \attention
       * In contrast to most other functions, this function uses \c int rather than <c>std::size_t</c>
       * as the type for the send/receive counts and displacements for technical reasons!
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[in] sendcounts
       * The number of datatype objects send to \e each process.
       *
       * \param[in] sdispls
       * The displacements of the send buffers in datatype objects.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] recvcounts
       * The number of datatype objects received from \e each process.
       *
       * \param[in] rdispls
       * The displacements of the receive buffers in datatype objects.
       *
       * \see \cite MPI31 Section 5.8, page 170
       */
      template<typename ST_, typename RT_>
      void alltoallv(const ST_* sendbuf, const int* sendcounts, const int* sdispls, RT_* recvbuf, const int* recvcounts, const int* rdispls) const
      {
        alltoallv(sendbuf, sendcounts, sdispls, autotype<ST_>(), recvbuf, recvcounts, rdispls, autotype<RT_>());
      }

      // end of gather/scatter group
      ///@}

      /**
       * \name Reductions and Scans
       */
      ///@{

      /**
       * \brief Blocking Reduce
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send/receive buffer contents.
       *
       * \param[in] op
       * A reference to the Operation object representing the reduction operation.
       *
       * \param[in] root
       * The rank of the root process that receives the reduction result.
       *
       * \see \cite MPI31 Section 5.9.1, page 174
       */
      void reduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op, int root) const;

      /**
       * \brief Blocking Reduce
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] op
       * A reference to the Operation object representing the reduction operation.
       *
       * \param[in] root
       * The rank of the root process that receives the reduction result.
       *
       * \see \cite MPI31 Section 5.9.1, page 174
       */
      template<typename T_>
      void reduce(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op, int root) const
      {
        reduce(sendbuf, recvbuf, count, autotype<T_>(), op, root);
      }

      /**
       * \brief Nonblocking Reduce
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send/receive buffer contents.
       *
       * \param[in] op
       * A reference to the Operation object representing the reduction operation.
       *
       * \param[in] root
       * The rank of the root process that receives the reduction result.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.7, page 209
       */
      Request ireduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op, int root) const;

      /**
       * \brief Nonblocking Reduce
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] op
       * A reference to the Operation object representing the reduction operation.
       *
       * \param[in] root
       * The rank of the root process that receives the reduction result.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.7, page 209
       */
      template<typename T_>
      Request ireduce(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op, int root) const
      {
        return ireduce(sendbuf, recvbuf, count, autotype<T_>(), op, root);
      }

      /**
       * \brief Blocking All-Reduce
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send/receive buffer contents.
       *
       * \param[in] op
       * A reference to the Operation object representing the reduction operation.
       *
       * \see \cite MPI31 Section 5.9.6, page 187
       */
      void allreduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const;

      /**
       * \brief Blocking All-Reduce
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] op
       * A reference to the Operation object representing the reduction operation.
       *
       * \see \cite MPI31 Section 5.9.6, page 187
       */
      template<typename T_>
      void allreduce(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op) const
      {
        allreduce(sendbuf, recvbuf, count, autotype<T_>(), op);
      }

      /**
       * \brief Nonblocking All-Reduce
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send/receive buffer contents.
       *
       * \param[in] op
       * A reference to the Operation object representing the reduction operation.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.8, page 210
       */
      Request iallreduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const;

      /**
       * \brief Nonblocking All-Reduce
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] op
       * A reference to the Operation object representing the reduction operation.
       *
       * \returns A request object for the operation.
       *
       * \see \cite MPI31 Section 5.12.8, page 210
       */
      template<typename T_>
      Request iallreduce(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op) const
      {
        return iallreduce(sendbuf, recvbuf, count, autotype<T_>(), op);
      }

      /**
       * \brief Blocking Inclusive Scan
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send/receive buffer contents.
       *
       * \param[in] op
       * A reference to the Operation object representing the scan operation.
       *
       * \see \cite MPI31 Section 5.11.1, page 193
       */
      void scan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const;

      /**
       * \brief Blocking Inclusive Scan
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] op
       * A reference to the Operation object representing the scan operation.
       *
       * \see \cite MPI31 Section 5.11.1, page 193
       */
      template<typename T_>
      void scan(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op) const
      {
        scan(sendbuf, recvbuf, count, autotype<T_>(), op);
      }

      /*Request iscan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const;

      template<typename T_>
      Request iscan(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op) const
      {
        return iscan(sendbuf, recvbuf, count, autotype<T_>(), op);
      }*/

      /**
       * \brief Blocking Exclusive Scan
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] datatype
       * A reference to the Datatype object representing the send/receive buffer contents.
       *
       * \param[in] op
       * A reference to the Operation object representing the scan operation.
       *
       * \see \cite MPI31 Section 5.11.2, page 194
       */
      void exscan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const;

      /**
       * \brief Blocking Exclusive Scan
       *
       * This function automatically deducts the datatype of the send/receive buffer(s) (if possible).
       *
       * \note
       * \p sendbuf and \p recvbuf are allowed to be identical; in this case, the
       * \p sendbuf argument of the MPI function call is set to \c MPI_IN_PLACE.
       *
       * \param[in] sendbuf
       * The send buffer for the operation.
       *
       * \param[out] recvbuf
       * The receive buffer for the operation.
       *
       * \param[in] count
       * The size of the send/receive buffer in datatype objects.
       *
       * \param[in] op
       * A reference to the Operation object representing the scan operation.
       *
       * \see \cite MPI31 Section 5.11.2, page 194
       */
      template<typename T_>
      void exscan(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op) const
      {
        exscan(sendbuf, recvbuf, count, autotype<T_>(), op);
      }

      /*Request iexscan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const;

      template<typename T_>
      Request iexscan(const T_* sendbuf, T_* recvbuf, std::size_t count, const Operation& op) const
      {
        return iexscan(sendbuf, recvbuf, count, autotype<T_>(), op);
      }*/

      // end of reductions group
      ///@}

      /**
       * \name Extended and derived communications
       *
       * The following function are not mere wrappers for MPI functions,
       * but implement more complex utility tasks.
       */
      ///@{

      /**
       * \brief Blocking broadcast of a std::stringstream
       *
       * \param[inout] stream
       * The stream that is to be broadcasted.
       *
       * \param[in] root
       * The rank of the root process that is sending the broadcast.
       */
      void bcast_stringstream(std::stringstream& stream, int root = 0) const;

      /**
       * \brief Blocking broadcast of a BinaryStream
       *
       * \param[inout] stream
       * The stream that is to be broadcasted.
       *
       * \param[in] root
       * The rank of the root process that is sending the broadcast.
       */
      void bcast_binarystream(BinaryStream& stream, int root = 0) const;

      /**
       * \brief Prints a message line to an output stream
       *
       * \param[in] os
       * The output stream to write to. Must be a valid stream on the root process.
       *
       * \param[in] msg
       * The message that is to be written to \p os.
       * This function automatically appends a new-line after the message.
       *
       * \param[in] root
       * The root process that should perform the print.
       */
      void print(std::ostream& os, const String& msg, int root = 0) const;

      /**
       * \brief Prints a message line to the standard output stream \c cout
       *
       * \param[in] msg
       * The message that is to be written to \p cout.
       * This function automatically appends a new-line after the message.
       *
       * \param[in] root
       * The root process that should perform the print.
       */
      void print(const String& msg, int root = 0) const
      {
        print(std::cout, msg, root);
      }

      /**
       * \brief Prints the ordered messages of all processes to an output stream
       *
       * \param[in] os
       * The output stream to write to. Must be a valid stream on the root process.
       *
       * \param[in] msg
       * The message that is to be written to \p os.
       * This function automatically appends a new-line after the message (of each rank)
       * and prefixes each output line with the rank of the author process.
       *
       * \param[in] root
       * The root process that should collect the messages and perform the print.
       */
      void allprint(std::ostream& os, const String& msg, int root = 0) const;

      /**
       * \brief Prints the ordered messages of all processes to the standard output stream \c cout
       *
       * \param[in] msg
       * The message that is to be written to \p cout.
       * This function automatically appends a new-line after the message (of each rank)
       * and prefixes each output line with the rank of the author process.
       *
       * \param[in] root
       * The root process that should collect the messages and perform the print.
       */
      void allprint(const String& msg, int root = 0) const
      {
        allprint(std::cout, msg, root);
      }

      // end of extended comm group
      ///@}
    }; // class Comm
  } // namespace Dist
} // namespace FEAT

#endif // KERNEL_UTIL_DIST_HPP
