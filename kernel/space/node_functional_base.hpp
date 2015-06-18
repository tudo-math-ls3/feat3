#pragma once
#ifndef KERNEL_SPACE_NODE_FUNCTIONAL_BASE_HPP
#define KERNEL_SPACE_NODE_FUNCTIONAL_BASE_HPP 1

// includes, FEAST
#include <kernel/space/base.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Node-functional base class template
     *
     * This class acts as a base-class and interface documentation for Finite-Element Node-Functional implementations.
     *
     * \tparam Space_
     * The finite-element space that this node-functional is used by.
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      typename DataType_>
    class NodeFunctionalBase
    {
    public:
      /// space typedef
      typedef Space_ SpaceType;
      /// data type
      typedef DataType_ DataType;

#ifdef DOXYGEN
      /// specifies the maximum number of assigned DOFs
      static constexpr Index max_assigned_dofs = ...
#endif // DOXYGEN

    protected:
      /// space reference
      const SpaceType& _space;
      /// currently active cell index
      Index _cell_index;

      /// protected constructor
      explicit NodeFunctionalBase(const SpaceType& space) :
        _space(space),
        _cell_index(~Index(0))
      {
      }

    public:
      /**
       * \brief Prepares the node-functional for a given cell.
       *
       * \param[in] cell_index
       * The index of the cell that is to be used by the node-functional.
       */
      void prepare(Index cell_index)
      {
        _cell_index = cell_index;
      }

      /**
       * \brief Releases the node-functional from the current cell.
       */
      void finish()
      {
        _cell_index = ~Index(0);
      }

#ifdef DOXYGEN
      /**
       * \brief Returns the number of assigned dofs on the current cell.
       */
      int get_num_assigned_dofs() const;

      /**
       * \brief Evaluation operator.
       *
       * This operator evaluates an assigned dof on the current cell.
       *
       * \param[in] assign_idx
       * The index of the assigned dofs that is to be evaluated.
       *
       * \returns
       * The value of the node functional applied onto the functor.
       */
      template<
        typename NodeData_,
        typename Function_>
      void operator()(NodeData_& node_data, const Function_& function) const;
#endif // DOXYGEN
    }; // class NodeFunctionalBase<...>

    /**
     * \brief Null-Node-Functional class template
     *
     * This class implements the node-functional interface for an empty node functional set.
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      typename DataType_>
    class NodeFunctionalNull :
      public NodeFunctionalBase<Space_, DataType_>
    {
    public:
      static constexpr Index max_assigned_dofs = Index(0);

    private:
      /// base-class typedef
      typedef NodeFunctionalBase<Space_, DataType_> BaseClass;

    public:
      /// constructor
      explicit NodeFunctionalNull(const Space_& space) :
        BaseClass(space)
      {
      }

      /** \copydoc NodeFunctionalBase::get_num_assigned_dofs() */
      int get_num_assigned_dofs() const
      {
        return max_assigned_dofs;
      }

      /** \copydoc NodeFunctionalBase::operator()() */
      template<
        typename NodeData_,
        typename Function_>
      void operator()(NodeData_&, const Function_&) const
      {
        throw InternalError("invalid call of NodeFunctionalNull::operator()()");
      }
    }; // class NodeFunctionalNull<...>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_NODE_FUNCTIONAL_BASE_HPP
