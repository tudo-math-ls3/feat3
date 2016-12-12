#pragma once
#ifndef KERNEL_LAFEM_SLIP_FILTER_HPP
#define KERNEL_LAFEM_SLIP_FILTER_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/lafem/arch/slip_filter.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Slip Filter class template.
     *
     * \tparam Mem_
     * Memory architecture
     *
     * \tparam DT_
     * Data type, i.e. double
     *
     * \tparam IT_
     * Index type, i.e. unsigned int
     *
     * \tparam BlockSize_
     * Size of the blocks, i.e. 2 for a filter acting on a velocity field in 2d
     *
     * An application of this filter sets a given vector component, i.e. it sets the normal component of a given
     * velocity vector to zero. This is for enforcing boundary conditions like
     * \f[ u \cdot \nu = g ~\mathrm{on}~ \Gamma \f]
     * in the strong sense. An application of the slip filter \f$ f_s \f$ is can be written as
     * \f[ f_s(u,g) = u - (u \cdot \nu) \nu + g \nu. \f]
     *
     * In the Finite Element case, this boundary condition can only be enforced in the strong sense if the FE space
     * \f$ V_h \f$ is at least \f$ \mathcal{C}^0 \f$. One possible approximation is the condition that
     * \f[ \forall j \in \Gamma_h: \mathcal{N}_j(f_{s,h} \cdot \nu_h) = \mathcal{N}_j(g), \f]
     * where \f$ \mathcal{N}_j \f$ is the node functional to the Dof \p j associated with \f$ \Gamma_h, f_{s,h} \f$
     * is a discrete slip filter mapping and \f$ \nu_h \f$ is some discrete approximation of \f$ \nu \f$. Assuming
     * \f$ g = 0 \f$ for simplicity, if \f$v_j\f$ is a basis function of \f$ V_h \f$, possible choice is
     * \f[ f_{s,h}(v_j) = v_j - \frac{\mathcal{N}_j(v_j \cdot \nu_h)}{\mathcal{N}_j (\nu_h \cdot \nu_h)} \nu_h. \f]
     * In the case of conforming Lagrange elements, this reduces to
     * \f[ f_{s,h}(v_h)(x_j) = v_h(x_j) - \frac{v_h \cdot \nu_h}{\nu_h \cdot \nu_h}(x_j) \nu(x_j) \f]
     * for all Lagrange nodes \f$ x_j \in \Gamma \f$.
     *
     * At the moment, the SlipFilterAssembler is only implemented for Lagrange 1 functions and P1/Q1 transformations.
     * For conforming Lagrange \p p spaces, interpolation of the outer unit normal field to the Lagrange \p p space
     * needs to be implemented. For other function spaces (especially nonconforming ones) it not clear if the above
     * approximation is sensible in any way.
     *
     * Even for spaces with enough regularity (at least \f$ \mathcal{C}^0 \f$), the application of Node functionals
     * to FE functions of a different space is missing.
     *
     * Also missing: Filtering of matrices.
     *
     * \see Assembly::SlipFilterAssembler
     *
     * \author Jordi Paul
     */
    template<
      typename Mem_,
      typename DT_,
      typename IT_,
      int BlockSize_>
    class SlipFilter
    {
      public:
        /// mem-type typedef
        typedef Mem_ MemType;
        /// data-type typedef
        typedef DT_ DataType;
        /// index-type typedef
        typedef IT_ IndexType;
        /// The block size
        static constexpr int BlockSize = BlockSize_;
        /// Value type
        typedef Tiny::Vector<DataType, BlockSize> ValueType;
        /// Our supported vector type
        typedef DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_> VectorType;

        /// Our 'base' class type
        template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
        using FilterType = class SlipFilter<Mem2_, DT2_, IT2_, BlockSize_>;

        /// For creating filter with different Mem, DT, IT
        template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
        using FilterTypeByMDI = FilterType<Mem2_, DT2_, IT2_>;

        static_assert(BlockSize > 1, "BlockSize has to be >= 2 in SlipFilter!");

      private:
        /// This will contain the pointwise outer unit normal in all vertices
        SparseVectorBlocked<Mem_, DT_, IT_, BlockSize_> _nu;
        /// This will contain the data for filtering
        SparseVectorBlocked<Mem_, DT_, IT_, BlockSize_> _sv;

      public:
        /// default constructor
        SlipFilter() :
          _nu(),
          _sv()
          {
          }

        /**
         * \brief Constructor.
         *
         * \param[in] num_vertices
         * Number of vertices in the mesh to save the outer unit normal field at.
         *
         * \param[in] num_dofs
         * Number of DoFs that will be filtered.
         *
         */
        explicit SlipFilter(Index num_vertices, Index num_dofs) :
          _nu(num_vertices),
          _sv(num_dofs)
          {
          }

        /// move-ctor
        SlipFilter(SlipFilter && other) :
          _nu(std::move(other._nu)),
          _sv(std::move(other._sv))
          {
          }

        /// move-assignment operator
        SlipFilter & operator=(SlipFilter && other)
        {
          if(this != &other)
          {
            _sv = std::forward<decltype(other._sv)>(other._sv);
            _nu = std::forward<decltype(other._nu)>(other._nu);
          }
          return *this;
        }

        /// virtual destructor
        virtual ~SlipFilter()
        {
        }

        /// \brief Creates a clone of itself
        SlipFilter clone(CloneMode clone_mode = CloneMode::Deep) const
        {
          SlipFilter other;
          other.clone(*this, clone_mode);
          return std::move(other);
        }

        /// \brief Clones data from another SlipFilter
        void clone(const SlipFilter & other, CloneMode clone_mode = CloneMode::Deep)
        {
          _nu.clone(other.get_nu(), clone_mode);
          _sv.clone(other.get_filter_vector(), clone_mode);
        }

        /// \brief Converts data from another UnitFilter
        template<typename Mem2_, typename DT2_, typename IT2_, int BS_>
        void convert(const SlipFilter<Mem2_, DT2_, IT2_, BS_>& other)
        {
          _nu.convert(other.get_nu());
          _sv.convert(other.get_filter_vector());
        }

        /// \brief Clears the underlying data (namely the SparseVector)
        void clear()
        {
          _nu.clear();
          _sv.clear();
        }

        /// \cond internal
        SparseVectorBlocked<Mem_, DT_, IT_, BlockSize>& get_filter_vector()
        {
          return _sv;
        }
        const SparseVectorBlocked<Mem_, DT_, IT_, BlockSize>& get_filter_vector() const
        {
          return _sv;
        }

        SparseVectorBlocked<Mem_, DT_, IT_, BlockSize>& get_nu()
        {
          return _nu;
        }
        const SparseVectorBlocked<Mem_, DT_, IT_, BlockSize>& get_nu() const
        {
          return _nu;
        }
        /// \endcond

        /**
         * \brief Adds one element to the filter
         *
         * \param[in] idx Index where to add
         * \param[in] val Value to add
         *
         **/
        void add(IndexType idx, ValueType val)
        {
          _sv(idx, val);
        }

        /// \returns The total size of the filter.
        Index size() const
        {
          return _sv.size();
        }

        /// \returns The number of entries in the filter.
        Index used_elements() const
        {
          return _sv.used_elements();
        }

        /// \returns The index array.
        IT_* get_indices()
        {
          return _sv.indices();
        }

        /// \returns The index array.
        const IT_* get_indices() const
        {
          return _sv.indices();
        }

        /// \returns The value array.
        ValueType* get_values()
        {
          return _sv.elements();
        }

        /// \returns The value array.
        const ValueType* get_values() const
        {
          return _sv.elements();
        }

        /**
         * \brief Applies the filter onto the right-hand-side vector.
         *
         * \param[in,out] vector
         * A reference to the right-hand-side vector to be filtered.
         */
        void filter_rhs(VectorType& vector) const
        {
          if(_sv.empty())
            return;
          if(_sv.size() != vector.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");
          if(_sv.used_elements() > Index(0))
            Arch::SlipFilter<Mem_>::template filter_rhs<DT_, IT_, BlockSize_>
              (vector.template elements<Perspective::pod>(), _sv.template elements<Perspective::pod>(), _sv.indices(), _sv.used_elements());
        }

        /**
         * \brief Applies the filter onto the solution vector.
         *
         * \param[in,out] vector
         * A reference to the solution vector to be filtered.
         */
        void filter_sol(VectorType& vector) const
        {
          if(_sv.empty())
            return;
          // same as rhs
          filter_rhs(vector);
        }

        /**
         * \brief Applies the filter onto a defect vector.
         *
         * \param[in,out] vector
         * A reference to the defect vector to be filtered.
         */
        void filter_def(VectorType& vector) const
        {
          if(_sv.empty())
            return;
          if(_sv.size() != vector.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");
          if(_sv.used_elements() > Index(0))
            Arch::SlipFilter<Mem_>::template filter_rhs<DT_, IT_, BlockSize_>
              (vector.template elements<Perspective::pod>(), _sv.template elements<Perspective::pod>(), _sv.indices(), _sv.used_elements() );
        }

        /**
         * \brief Applies the filter onto a correction vector.
         *
         * \param[in,out] vector
         * A reference to the correction vector to be filtered.
         */
        void filter_cor(VectorType& vector) const
        {
          // same as def
          filter_def(vector);
        }
    }; // class SlipFilter<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_SLIP_FILTER_HPP
