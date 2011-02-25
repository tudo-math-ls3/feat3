#pragma once
#ifndef UTIL_INSTANTIATION_POLICY_HPP
/// Header guard
#define UTIL_INSTANTIATION_POLICY_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  /**
  * \{
  * \name instantiation policy tags
  * \author Dirk Ribbrock
  *
  * Possible tags for InstantiationPolicy.
  */

  /// Tag for an instantiation policy that does not allow copying.
  struct NonCopyable;

  /// Tag for an instantiation policy that does not allow instantiation at all.
  struct NonInstantiable;

  /// Tag for an instantiation policy that does not allow more than one instance of the class.
  struct Singleton;

  /// \}


  /**
  * \brief utility base class to restrict instantiation behaviour of its descendants
  *
  * \author Dirk Ribbrock
  */
  template<
    typename T_,
    typename Method_>
  class InstantiationPolicy;


  /**
  * \brief base class to prevent copying
  *
  * \author Dirk Ribbrock
  */
  template<typename T_>
  class InstantiationPolicy<T_, NonCopyable>
  {

  private:

    /// Unwanted copy constructor: Do not implement!
    InstantiationPolicy(const InstantiationPolicy &);


    /// Unwanted copy assignment operator: Do not implement!
    InstantiationPolicy & operator= (const InstantiationPolicy &);


  public:

    /// Default constructor.
    InstantiationPolicy()
    {
    }
  };


  /**
  * \brief base class to prevent instantiation
  *
  * \author Dirk Ribbrock
  */
  template<typename T_>
  class InstantiationPolicy<T_, NonInstantiable>
  {

  private:

    /// Unwanted copy constructor: Do not implement!
    InstantiationPolicy(const InstantiationPolicy &);

    /// Unwanted copy assignment operator: Do not implement!
    InstantiationPolicy & operator= (const InstantiationPolicy &);

    /// Unwanted default constructor: Do not implement!
    InstantiationPolicy();
  };


  /**
  * \brief base class to prevent multiple instantiation
  *
  * \author Dirk Ribbrock
  */
  template <typename T_>
  class InstantiationPolicy<T_, Singleton>
  {

  private:

    class DeleteOnDestruction;

    friend class DeleteOnDestruction;

    /// Returns a pointer to our instance pointer.
    static T_ * * _instance_ptr()
    {
      static T_ * instance(nullptr);
      static DeleteOnDestruction delete_instance(&instance);

      return &instance;
    }


    /// Deletes the object of T_ that is pointed at by ptr.
    static void _delete(T_ * const ptr)
    {
      delete ptr;
    }


    /// Unwanted copy constructor: Do not implement!
    InstantiationPolicy(const InstantiationPolicy &);


    /// Unwanted copy assignment operator: Do not implement!
    InstantiationPolicy & operator= (const InstantiationPolicy &);


  protected:
    /// Default constructor.
    InstantiationPolicy()
    {
    }


  public:
    /// Returns the instance.
    static T_ * instance()
    {
      T_ * * instance_ptr(_instance_ptr());

      if (nullptr == *instance_ptr)
      {
        /// \todo Make thread safe
        //static Mutex m;
        //Lock l(m);

        instance_ptr = _instance_ptr();

        if (nullptr == *instance_ptr)
        {
          *instance_ptr = new T_;
        }
      }
      return *instance_ptr;
    }
  };



  /**
  * \brief allows automatic deletion of the destructed object
  *
  * \author Dirk Ribbrock
  */
  template <typename T_>
  class InstantiationPolicy<T_, Singleton>::DeleteOnDestruction
  {
  private:

    /// Pointer to our managed object
    T_ * * const _ptr;


  public:

    /// Constructor
    DeleteOnDestruction(T_ * * const ptr)
      : _ptr(ptr)
    {
    }


    /// Destructor
    ~DeleteOnDestruction()
    {
      InstantiationPolicy<T_, Singleton>::_delete(*_ptr);
      *_ptr = nullptr;
    }
  };
} // namespace FEAST

#endif // UTIL_INSTANTIATION_POLICY_HPP
