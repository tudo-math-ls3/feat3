// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_VELOCITY_ANALYSER_HPP
#define KERNEL_ASSEMBLY_VELOCITY_ANALYSER_HPP 1

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/dist.hpp>

namespace FEAT
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
      /// The dimension of the analysed velocity field
      static constexpr int dim = dim_;

      /**
       * \brief Velocity field H0-norm
       *
       * This entry contains the H0-norm of the velocity field:
       * \f[ \|u\|_{\mathcal{H}^0} = \Big(\int_\Omega \sum_{i=1}^d (u_i)^2\Big)^{\frac{1}{2}}\f]
       */
      DataType_ norm_h0;

      /**
       * \brief Velocity field H1-semi-norm
       *
       * This entry contains the H1-semi-norm of the velocity field:
       * \f[ |u|_{\mathcal{H}^1} = \Big(\int_\Omega \sum_{i=1}^d (\nabla u_i)^2\Big)^{\frac{1}{2}}\f]
       */
      DataType_ norm_h1;

      /**
       * \brief Velocity field divergence H0-norm
       *
       * This entry contains the H0-norm of the divergence of the velocity field:
       * \f[ \| \nabla \cdot u\|_{\mathcal{H}^0} = \Big(\int_\Omega \big(\sum_{i=1}^d \partial_i u_i\big)^2\Big)^{\frac{1}{2}}\f]
       */
      DataType_ divergence;

      /**
       * \brief Velocity field vorticity H0-norm
       *
       * This entry contains the H0-norm of the curl of the velocity field:
       * \f[ \| \nabla \times u\|_{\mathcal{H}^0} \f]
       */
      DataType_ vorticity;

      /**
       * \brief Velocity field components H0-norm
       *
       * This entry contains the H0-norms of all components of the velocity field.
       */
      Tiny::Vector<DataType_, dim_> norm_h0_comp;

      /**
       * \brief Velocity field components H1-semi-norm
       *
       * This entry contains the H1-semi-norms of all components of the velocity field.
       */
      Tiny::Vector<DataType_, dim_> norm_h1_comp;

      /// default CTOR
      VelocityInfo() :
        norm_h0(DataType_(0)),
        norm_h1(DataType_(0)),
        divergence(DataType_(0)),
        vorticity(DataType_(0)),
        norm_h0_comp(DataType_(0)),
        norm_h1_comp(DataType_(0))
      {
      }

      /**
       * \brief Synchronizes the velocity information over a communicator
       *
       * This function sums up the velocity information of all patches in a
       * parallel simulation to obtain the information for the global mesh.
       *
       * \param[in] comm
       * The communication over which to synchronize.
       */
      void synchronize(const Dist::Comm& comm)
      {
        DataType_ verr[2*dim_+4] =
        {
          Math::sqr(norm_h0),
          Math::sqr(norm_h1),
          Math::sqr(divergence),
          Math::sqr(vorticity),
        };
        for(int i(0); i < dim_; ++i)
        {
          verr[4+     i] = Math::sqr(norm_h0_comp[i]);
          verr[4+dim_+i] = Math::sqr(norm_h1_comp[i]);
        }

        comm.allreduce(verr, verr, std::size_t(2*dim_+4), Dist::op_sum);

        norm_h0    = Math::sqrt(verr[0]);
        norm_h1    = Math::sqrt(verr[1]);
        divergence = Math::sqrt(verr[2]);
        vorticity  = Math::sqrt(verr[3]);
        for(int i(0); i < dim_; ++i)
        {
          norm_h0_comp[i] = Math::sqrt(verr[4+     i]);
          norm_h1_comp[i] = Math::sqrt(verr[4+dim_+i]);
        }
      }


      /**
       * \brief Formats the velocity information as a string.
       *
       * \param[in] precision
       * The precision for floating point values.
       *
       * \param[in] pad_size
       * The leading string padding size. Should be >= 10 to achieve vertical alignment.
       *
       * \param[in] pad_char
       * The leading string padding character.
       *
       * \returns
       * A string containing the formatted velocity field info.
       */
      String format_string(int precision = 0, std::size_t pad_size = 10u, char pad_char = '.') const
      {
        String s;
        s += String("H0-Norm").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(norm_h0, precision) + " [";
        for(int i(0); i < dim_; ++i)
          (s += " ") += stringify_fp_sci(norm_h0_comp[i], precision);
        s += " ]\n";
        s += String("H1-Norm").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(norm_h1, precision) + " [";
        for(int i(0); i < dim_; ++i)
          (s += " ") += stringify_fp_sci(norm_h1_comp[i], precision);
        s += " ]\n";
        s += String("Divergence").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(divergence, precision) + "\n";
        s += String("Vorticity").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(vorticity, precision) + "\n";
        return s;
      }

      /// prints the info to an output stream
      friend std::ostream& operator<<(std::ostream& os, const VelocityInfo& vi)
      {
        // print H0-norms
        os << "H0-Norm...: " << stringify_fp_sci(vi.norm_h0) << " [";
        for(int i(0); i < dim_; ++i)
          os << " " << stringify_fp_sci(vi.norm_h0_comp[i]);
        os << " ]" << std::endl;

        // print H1-norms
        os << "H1-Norm...: " << stringify_fp_sci(vi.norm_h1) << " [";
        for(int i(0); i < dim_; ++i)
          os << " " << stringify_fp_sci(vi.norm_h1_comp[i]);
        os << " ]" << std::endl;

        // print divergence and vorticity norms
        os << "Divergence: " << stringify_fp_sci(vi.divergence) << std::endl;
        os << "Vorticity.: " << stringify_fp_sci(vi.vorticity) << std::endl;

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
      template<typename DT_, int sm_, int sn_>
      static DT_ calc_vorticity(const Tiny::Matrix<DT_, 2, 2, sm_, sn_>& jac)
      {
        // 2D: (dx u2 - dy u1)^2
        // J_ij = d/d_j f_i ==> (J_01 - J_10)^2
        return Math::sqr(jac[1][0] - jac[0][1]);
      }

      template<typename DT_, int sm_, int sn_>
      static DT_ calc_vorticity(const Tiny::Matrix<DT_, 3, 3, sm_, sn_>& jac)
      {
        // 3D: (dy u3 - dz u2)^2 + (dz u1 - dx u3)^2 + (dx u2 - dy u1)
        return
          Math::sqr(jac[2][1] - jac[1][2]) + // dy u3 - dz u2 = J_21 - J_12
          Math::sqr(jac[0][2] - jac[2][0]) + // dz u1 - dx u3 = J_02 - J_20
          Math::sqr(jac[1][0] - jac[0][1]);  // dx u2 - dy u1 = J_10 - J_01
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
        typedef AsmTraits1<DataType_, SpaceType, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

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

        // initialize result
        VelocityInfo<DataType_, dim_> info;

        // local velocity values
        Tiny::Matrix<DataType_, dim_, SpaceEvaluator::max_local_dofs> basis_val;

        // vector field values
        Tiny::Vector<DataType_, dim_> val;

        // vector field derivatives
        Tiny::Matrix<DataType_, dim_, dim_> jac;

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // initialize dof-mapping
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
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // compute basis function values and derivatives
            val.format();
            jac.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              for(int bi(0); bi < dim_; ++bi)
              {
                val[bi] += basis_val[bi][i] * space_data.phi[i].value;
                // J_ij = d/d_j f_i  <==> J = (f_1,...,f_n)^T * (d/d_1, ..., d/d_n)
                for(int bj(0); bj < dim_; ++bj)
                {
                  // J_ij = d/d_j f_i
                  jac[bi][bj] += basis_val[bi][i] * space_data.phi[i].grad[bj];
                }
              }
            }

            // compute integration weight
            DataType omega = trafo_data.jac_det * cubature_rule.get_weight(k);

            // update info
            for(int i(0); i < dim_; ++i)
            {
              // update H0 component norm
              info.norm_h0_comp[i] += omega * Math::sqr(val[i]);

              // update H1 component norm
              info.norm_h1_comp[i] += omega * Tiny::dot(jac[i], jac[i]);
            }

            // update divergence
            info.divergence += omega * Math::sqr(jac.trace());

            // update vorticity
            info.vorticity += omega * calc_vorticity(jac);

            // continue with next cubature point
          }

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // finally, compute the rest
        for(int i(0); i < dim_; ++i)
        {
          info.norm_h0 += info.norm_h0_comp[i];
          info.norm_h1 += info.norm_h1_comp[i];
          info.norm_h0_comp[i] = Math::sqrt(info.norm_h0_comp[i]);
          info.norm_h1_comp[i] = Math::sqrt(info.norm_h1_comp[i]);
        }

        // take the roots
        info.norm_h0 = Math::sqrt(info.norm_h0);
        info.norm_h1 = Math::sqrt(info.norm_h1);
        info.divergence = Math::sqrt(info.divergence);
        info.vorticity  = Math::sqrt(info.vorticity);

        // okay
        return info;
      }

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
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const Space_& space, const CubatureFactory_& cubature_factory)
      {
        // first of all, verify the dimensions
        static_assert(Space_::shape_dim == dim_, "invalid velocity field dimension");

        // validate vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");

        /// space type
        typedef Space_ SpaceType;

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        /// assembly traits
        typedef AsmTraits1<DataType_, SpaceType, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

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
        typename VectorType::GatherAxpy gather(vector);

        // initialize result
        VelocityInfo<DataType_, dim_> info;

        // local velocity values
        Tiny::Vector<Tiny::Vector<DataType, dim_>, SpaceEvaluator::max_local_dofs> basis_val;

        // vector field value
        Tiny::Vector<DataType_, dim_> val;

        // vector field first derivatives as jacobian matrix
        Tiny::Matrix<DataType_, dim_, dim_> jac;

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch local basis values
          basis_val.format();
          gather(basis_val, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // compute basis function values and derivatives
            val.format();
            jac.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              val.axpy(space_data.phi[i].value, basis_val[i]);
              // J_ij = d/d_j f_i  <==> J = (f_1,...,f_n)^T * (d/d_1, ..., d/d_n)
              jac.add_outer_product(basis_val[i], space_data.phi[i].grad);
            }

            // compute integration weight
            const DataType omega = trafo_data.jac_det * cubature_rule.get_weight(k);

            // update norms
            for(int i(0); i < dim_; ++i)
            {
              // update H0 component norm
              info.norm_h0_comp[i] += omega * Math::sqr(val[i]);

              // update H1 component norm
              info.norm_h1_comp[i] += omega * Tiny::dot(jac[i], jac[i]);
            }

            // update divergence
            info.divergence += omega * Math::sqr(jac.trace());

            // update vorticity
            info.vorticity += omega * calc_vorticity(jac);

            // continue with next cubature point
          }

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // finally, compute the rest
        for(int i(0); i < dim_; ++i)
        {
          info.norm_h0 += info.norm_h0_comp[i];
          info.norm_h1 += info.norm_h1_comp[i];
          info.norm_h0_comp[i] = Math::sqrt(info.norm_h0_comp[i]);
          info.norm_h1_comp[i] = Math::sqrt(info.norm_h1_comp[i]);
        }

        // take the roots
        info.norm_h0 = Math::sqrt(info.norm_h0);
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

        typename LAFEM::DenseVector<Mem::Main, DT_, IT_>::GatherAxpy _first_gather;
        RestClass _rest_gather;

        MultiGather(const LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>, dim_>& vector) :
          _first_gather(vector.first()), _rest_gather(vector.rest())
        {
        }

        template<int m_, int n_, int sm_, int sn_, typename DofMapping_>
        void gather(
          Tiny::Matrix<DT_, m_, n_, sm_, sn_>& vals,
          const DofMapping_& dof_mapping,
          int offset = 0)
        {
          _first_gather(vals[offset], dof_mapping);
          _rest_gather.gather(vals, dof_mapping, offset+1);
        }
      };

      template<typename DT_, typename IT_>
      struct MultiGather<LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>, 1> >
      {
        typename LAFEM::DenseVector<Mem::Main, DT_, IT_>::GatherAxpy _first_gather;

        MultiGather(const LAFEM::PowerVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>, 1>& vector) :
          _first_gather(vector.first())
        {
        }

        template<int m_, int n_, int sm_, int sn_, typename DofMapping_>
        void gather(
          Tiny::Matrix<DT_, m_, n_, sm_, sn_>& vals,
          const DofMapping_& dof_mapping,
          int offset = 0)
        {
          _first_gather(vals[offset], dof_mapping);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_VELOCITY_ANALYSER_HPP
