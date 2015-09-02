#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONALS_LVLSET
#define KERNEL_MESHOPT_RUMPF_FUNCTIONALS_LVLSET 1

#include <kernel/base_header.hpp>
#include <kernel/assembly/common_functions.hpp> // For HeavisideReg
#include <kernel/geometry/intern/face_index_mapping.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /// \cond internal

    /**
     * \brief Class template for Rumpf functionals with levelset term
     **/
    template<typename DataType_, typename ShapeType_>
    class RumpfFunctionalLevelset
    {
      public:
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type
        typedef ShapeType_ ShapeType;

        /// Type that maps face-local to cell-local vertex indices
        typedef Geometry::Intern::FaceIndexMapping<ShapeType, 1, 0> FimType;

        /// Factor for the levelset penalty term
        DataType fac_lvlset;
        /// Factor for making the regularised Heaviside function steeper
        static DataType heaviside_reg_fac(){ return DataType(1); };
        //static constexpr DataType heaviside_reg_fac = DataType(1);

      public:
        /**
         * \brief Constructor
         */
        RumpfFunctionalLevelset() :
          fac_lvlset(DataType(1))
          {
          }

        ///**
        // * \brief Copy constructor
        // */
        //RumpfFunctionalLevelset(const RumpfFunctionalLevelset& other) :
        //  fac_lvlset(other.fac_lvlset)
        //  {
        //  }

        /**
         * \brief Destructor
         */
        ~RumpfFunctionalLevelset()
        {
        }

        /**
         * \brief Computes the additional levelset penalty term for the Rumpf functional
         **/
        template<typename Tl_>
        DataType compute_lvlset_penalty(const Tl_& lvlset_vals) const
        {
          DataType penalty(0);
          // Loop over all edges of the current cell
          for(int edge(0); edge < Shape::FaceTraits<ShapeType,1>::count; ++edge)
          {
            // These are the indices of the vertices on edge
            int i(FimType::map(edge,0));
            int j(FimType::map(edge,1));

            penalty += FEAST::Assembly::Common::template HeavisideRegStatic<DataType>::eval(-heaviside_reg_fac()
                * lvlset_vals(i) * lvlset_vals(j));
          }

          // This is the version for penalising the diagonal cuts as well
          /*
             DataType penalty(0);
             for(Index i(0); i < Index(4); ++i)
             {
             for(Index j(0); j < i; ++j)
             penalty += FEAST::Assembly::Common::template HeavisideRegStatic<DataType>::eval(-heaviside_reg_fac*lvlset_vals(i)*lvlset_vals(j));
             }
             */

          return penalty;
        }

        /**
         * \brief Adds the gradient of the additional levelset penalty term
         **/
        template<typename Tl_, typename Tlg_, typename Tgrad_>
        void add_lvlset_penalty_grad(
          const Tl_& lvlset_vals,
          const Tlg_& lvlset_grad_vals,
          Tgrad_& grad,
          DataType lvlset_constraint_last)
          {
            // Compute local gradient
            for(int edge(0); edge < Shape::FaceTraits<ShapeType,1>::count; ++edge)
            {
              // These are the indices of the vertices on edge
              int i(FimType::map(edge,0));
              int j(FimType::map(edge,1));

              auto lvlset_prod = -heaviside_reg_fac() * lvlset_vals(i) * lvlset_vals(j);
              // Derivative of the heaviside function
              auto heaviside_der = FEAST::Assembly::Common::template HeavisideRegStatic<DataType>::der_x(lvlset_prod);
              for(int d(0); d < ShapeType::dimension; ++d)
              {
                grad(i,d) -= heaviside_reg_fac() * fac_lvlset * lvlset_constraint_last *
                  (heaviside_der * lvlset_grad_vals(i,d) * lvlset_vals(j));
                grad(j,d) -= heaviside_reg_fac() * fac_lvlset * lvlset_constraint_last *
                  (heaviside_der * lvlset_grad_vals(j,d) * lvlset_vals(i));
              }
            }

            // This is the version for penalising the diagonal cuts as well
            /*
            // Compute local gradient
            for(Index i(0); i < Index(4); ++i)
            {
            for(Index j(0); j < i; ++j)
            {
            auto lvlset_prod = -heaviside_reg_fac*lvlset_vals(i)*lvlset_vals(j);
            // Derivative of the heaviside function
            auto heaviside_der = FEAST::Assembly::Common::template HeavisideRegStatic<DataType>::der_x(lvlset_prod);
            for(Index d(0); d < 2; ++d)
            grad(d,i) -= heaviside_reg_fac*fac_lvlset*Math::sqr(lvlset_constraint_last) * (heaviside_der * lvlset_grad_vals(d,i) * lvlset_vals(j));
            }
            }
            */
          }

    }; // class RumpfFunctionalLevelset
    /// \endcond
  } // namespace Meshopt
} // namespace FEAST

#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONALS_LVLSET
