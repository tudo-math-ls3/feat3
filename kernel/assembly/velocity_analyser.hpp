#pragma once
#ifndef KERNEL_ASSEMBLY_VELOCITY_ANALYSER_HPP
#define KERNEL_ASSEMBLY_VELOCITY_ANALYSER_HPP 1

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Vector_>
      struct MultiGather;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Velocity field information structure
     *
     * This structure encapsulates various information which is computed and returned
     * by the VelocityAnalser class.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, int dim_>
    struct VelocityInfo
    {
      /**
       * \brief Velocity field L2-norm
       *
       * This entry contains the L2-norm of the velocity field:
       * \f[ \|u\|_{\mathcal{L}^2} = \Big(\int_\Omega \sum_{i=1}^d (u_i)^2\Big)^{\frac{1}{2}}\f]
       */
      DataType_ norm_l2;

      /**
       * \brief Velocity field H1-semi-norm
       *
       * This entry contains the H1-semi-norm of the velocity field:
       * \f[ |u|_{\mathcal{H}^^2} = \Big(\int_\Omega \sum_{i=1}^d (\nabla u_i)^2\Big)^{\frac{1}{2}}\f]
       */
      DataType_ norm_h1;

      /**
       * \brief Velocity field divergence L2-norm
       *
       * This entry contains the L2-norm of the divergence of the velocity field:
       * \f[ \| \nabla \dot u\|_{\mathcal{L}^2} = \Big(\int_\Omega \big(\sum_{i=1}^d \partial_i u_i\big)^2\Big)^{\frac{1}{2}}\f]
       */
      DataType_ divergence;

      /**
       * \brief Velocity field vorticity L2-norm
       *
       * This entry contains the L2-norm of the vorticity of the velocity field:
       * \f[ \| \nabla \times u\|_{\mathcal{L}^2} \f]
       */
      DataType_ vorticity;

      /**
       * \brief Velocity field components L2-norm
       *
       * This entry contains the L2-norms of all components of the velocity field.
       */
      Tiny::Vector<DataType_, dim_> norm_l2_comp;

      /**
       * \brief Velocity field components H1-semi-norm
       *
       * This entry contains the H1-semi-norms of all components of the velocity field.
       */
      Tiny::Vector<DataType_, dim_> norm_h1_comp;

      /// default CTOR
      VelocityInfo() :
        norm_l2(DataType_(0)),
        norm_h1(DataType_(0)),
        divergence(DataType_(0)),
        vorticity(DataType_(0)),
        norm_l2_comp(DataType_(0)),
        norm_h1_comp(DataType_(0))
      {
      }

      /// prints the info to an output stream
      friend std::ostream& operator<<(std::ostream& os, const VelocityInfo& vi)
      {
        // print L2-norms
        os << "L2-Norm...: " << scientify(vi.norm_l2) << " [";
        for(int i(0); i < dim_; ++i)
          os << " " << scientify(vi.norm_l2_comp[Index(i)]);
        os << " ]" << std::endl;

        // print L2-norms
        os << "H1-Norm...: " << scientify(vi.norm_h1) << " [";
        for(int i(0); i < dim_; ++i)
          os << " " << scientify(vi.norm_h1_comp[Index(i)]);
        os << " ]" << std::endl;

        // print divergence and vorticity norms
        os << "Divergence: " << scientify(vi.divergence) << std::endl;
        os << "Vorticity.: " << scientify(vi.vorticity) << std::endl;

        return os;
      }
    }; // struct VelocityInfo

    /**
     * \brief Velocity field analyser class
     *
     * This class can be used for post-processing of a 2D or 3D velocity vector field.
     * For more information about the computed quantities, see the VelocityInfo struct.
     * \author Peter Zajac
     */
    class VelocityAnalyser
    {
    private:
      /// \cond internal
      struct TrafoConfig : public Trafo::ConfigBase
      {
        static constexpr bool need_jac_det = true;
      };
      struct SpaceConfig : public Space::ConfigBase
      {
        static constexpr bool need_value = true;
        static constexpr bool need_grad = true;
      };

      template<typename DT_, int sm_, int sn_>
      static DT_ calc_vorticity(const Tiny::Matrix<DT_, 2, 2, sm_, sn_>& ders)
      {
        // 2D: (dx u2 - dy u1)^2
        return Math::sqr(ders[0][1] - ders[1][0]);
      }

      template<typename DT_, int sm_, int sn_>
      static DT_ calc_vorticity(const Tiny::Matrix<DT_, 3, 3, sm_, sn_>& ders)
      {
        // 3D: (dy u3 - dz u2)^2 + (dz u1 - dx u3)^2 + (dx u2 - dy u1)
        return
          Math::sqr(ders[2][1] - ders[1][2]) +
          Math::sqr(ders[0][2] - ders[2][0]) +
          Math::sqr(ders[0][1] - ders[1][0]);
      }
      /// \endcond

    public:
      /**
       * \brief Performs the analysis of a velocity field.
       *
       * \param[in] vector
       * The vector representing the velocity field to be analysed.
       *
       * \param[in] space
       * The finite element space for the velocity.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \returns
       * A VelocityInfo structure containing the computed information.
       */
      template<typename DataType_, typename IndexType_, int dim_, typename Space_, typename CubatureFactory_>
      static VelocityInfo<DataType_, dim_> compute(
        const LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>, dim_>& vector,
        const Space_& space, const CubatureFactory_& cubature_factory)
      {
        // first of all, verify the dimensions
        static_assert(Space_::shape_dim == dim_, "invalid velocity field dimension");

        /// space type
        typedef Space_ SpaceType;

        typedef LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>, dim_> VectorType;

        /// assembly traits
        typedef AsmTraits1<DataType_, SpaceType, TrafoConfig, SpaceConfig> AsmTraits;

        /// data type
        typedef typename AsmTraits::DataType DataType;

        // create the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);
        typedef typename AsmTraits::SpaceEvaluator  SpaceEvaluator;

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create matrix scatter-axpy
        Intern::MultiGather<VectorType> multi_gather(vector);

        // initialise result
        VelocityInfo<DataType_, dim_> info;

        // local velocity values
        Tiny::Matrix<DataType_, dim_, SpaceEvaluator::max_local_dofs> basis_val;

        // vector field values
        Tiny::Vector<DataType_, dim_> vals;

        // vector field derivatives
        Tiny::Matrix<DataType_, dim_, dim_> ders;

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // fetch local basis values
          basis_val.format();
          multi_gather.gather(basis_val, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // fetch number of local dofs
          Index num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(Index k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // compute basis function values and derivatives
            vals.format();
            ders.format();
            for(Index i(0); i < num_loc_dofs; ++i)
            {
              for(Index bi(0); bi < Index(dim_); ++bi)
              {
                vals[bi] += basis_val[bi][i] * space_data.phi[i].value;
                for(Index bj(0); bj < Index(dim_); ++bj)
                {
                  ders[bi][bj] += basis_val[bi][i] * space_data.phi[i].grad[bj];
                }
              }
            }

            // compute integration weight
            DataType omega = trafo_data.jac_det * cubature_rule.get_weight(k);

            // update info
            for(Index i(0); i < Index(dim_); ++i)
            {
              // update L2 component norm
              info.norm_l2_comp[i] += omega * Math::sqr(vals[i]);

              // update H1 component norm
              info.norm_h1_comp[i] += omega * Tiny::dot(ders[i], ders[i]);

            }

            // update divergence
            info.divergence += omega * Math::sqr(ders.trace());

            // update vorticity
            info.vorticity += omega * calc_vorticity(ders);

            // continue with next cubature point
          }

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // finally, compute the rest
        for(Index i(0); i < Index(dim_); ++i)
        {
          info.norm_l2 += info.norm_l2_comp[i];
          info.norm_h1 += info.norm_h1_comp[i];
          info.norm_l2_comp[i] = Math::sqrt(info.norm_l2_comp[i]);
          info.norm_h1_comp[i] = Math::sqrt(info.norm_h1_comp[i]);
        }

        // take the roots
        info.norm_l2 = Math::sqrt(info.norm_l2);
        info.norm_h1 = Math::sqrt(info.norm_h1);
        info.divergence = Math::sqrt(info.divergence);
        info.vorticity  = Math::sqrt(info.vorticity);

        // okay
        return info;
      }
    }; // class VelocityAnalyser

    /// \cond internal
    namespace Intern
    {
      template<typename DT_, typename IT_, int dim_>
      struct MultiGather<LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>, dim_> >
      {
        typedef MultiGather<LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>, dim_-1> > RestClass;

        LAFEM::GatherAxpy<LAFEM::DenseVector<Mem::Main, DT_, IT_> > _first_gather;
        RestClass _rest_gather;

        MultiGather(const LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>, dim_>& vector) :
          _first_gather(vector.first()), _rest_gather(vector.rest())
        {
        }

        template<int m_, int n_, int sm_, int sn_, typename DofMapping_>
        void gather(
          Tiny::Matrix<DT_, m_, n_, sm_, sn_>& vals,
          const DofMapping_& dof_mapping,
          Index offset = Index(0))
        {
          _first_gather(vals[offset], dof_mapping);
          _rest_gather.gather(vals, dof_mapping, offset+1);
        }
      };

      template<typename DT_, typename IT_>
      struct MultiGather<LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>, 1> >
      {
        LAFEM::GatherAxpy<LAFEM::DenseVector<Mem::Main, DT_, IT_>> _first_gather;

        MultiGather(const LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>, 1>& vector) :
          _first_gather(vector.first())
        {
        }

        template<int m_, int n_, int sm_, int sn_, typename DofMapping_>
        void gather(
          Tiny::Matrix<DT_, m_, n_, sm_, sn_>& vals,
          const DofMapping_& dof_mapping,
          Index offset = Index(0))
        {
          _first_gather(vals[offset], dof_mapping);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_VELOCITY_ANALYSER_HPP
