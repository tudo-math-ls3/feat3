#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_LVLSET_2D_P1_HPP
#define KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_LVLSET_2D_P1_HPP 1

#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1.hpp>
// For HeavisideReg
#include <kernel/assembly/common_functions.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal

    /**
     * \brief Class template for Rumpf functionals with levelset term, 2d P1
     **/
    template<typename MemoryType_, typename DataType_>
    class RumpfFunctionalLevelset < MemoryType_, DataType_, Shape::Simplex<2> > :
    public RumpfFunctional < MemoryType_, DataType_, Shape::Simplex<2> >
    {
      public:
        /// Shape type of the underlying transformation
        typedef Shape::Simplex<2> ShapeType;
        /// Our base class
        typedef RumpfFunctional<MemoryType_, DataType_, ShapeType> BaseClass;

        /// Factor for the levelset penalty term
        DataType_ fac_lvlset;
        /// Factor for making the regularised Heaviside function steeper
        static DataType_ heaviside_reg_fac(){ return DataType_(10); };
        //static constexpr DataType_ heaviside_reg_fac = DataType_(1);

      public:
        /**
         * \brief Constructor
         **/
        RumpfFunctionalLevelset(
          const DataType_ fac_norm_,
          const DataType_ fac_det_,
          const DataType_ fac_cof_,
          const DataType_ fac_reg_) :
          BaseClass( fac_norm_,
          fac_det_,
          fac_cof_,
          fac_reg_),
          fac_lvlset(DataType_(10))
          {
          }

        /**
         * \brief Computes the additional levelset penalty term for the Rumpf functional
         **/
        template<typename Tl_>
        DataType_ compute_lvlset_penalty(const Tl_& lvlset_vals)
        {
          DataType_ penalty(0);
          for(Index i(0); i < Index(3); ++i)
          {
            for(Index j(0); j < i; ++j)
              penalty += FEAST::Assembly::Common::template HeavisideRegStatic<DataType_>::eval(-heaviside_reg_fac()*lvlset_vals(i)*lvlset_vals(j));
          }

          return penalty;
        } // compute_lvlset_penalty

        /**
         * \brief Adds the gradient of the additional levelset penalty term
         **/
        template<typename Tl_, typename Tlg_, typename Tgrad_>
        void add_lvlset_penalty_grad(const Tl_& lvlset_vals, const Tlg_& lvlset_grad_vals, Tgrad_& grad, DataType_& lvlset_constraint_last)
        {
          // Compute local gradient
          for(Index i(0); i < Index(3); ++i)
          {
            for(Index j(0); j < i; ++j)
            {
              auto lvlset_prod = -heaviside_reg_fac()*lvlset_vals(i)*lvlset_vals(j);
              // Derivative of the heaviside function
              auto heaviside_der = FEAST::Assembly::Common::template HeavisideRegStatic<DataType_>::der_x(lvlset_prod);
              for(Index d(0); d < Index(2); ++d)
                grad(d,i) -= heaviside_reg_fac() * fac_lvlset * lvlset_constraint_last * heaviside_der *
                  (lvlset_grad_vals(d,i) * lvlset_vals(j));
            }
          }

        } // add_levelset_penalty_grad

    }; // class RumpfFunctionalLevelset
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_LVLSET_2D_P1_HPP
