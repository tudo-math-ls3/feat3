#pragma once
#ifndef KERNEL_LAFEM_BICGSTAB_HPP
#define KERNEL_LAFEM_BICGSTAB_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/preconditioner.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief BiCGStab solver implementation
     *
     * This is just a temporary affair until there are some real(TM) solvers
     *
     **/
    struct BiCGStab
    {
      /**
       * \brief Applies the solver to a given system of linear equations.
       **/
      template <typename MT_, typename VT_>
      static void value(VT_ & x, const MT_ & A, const VT_ & b, Preconditioner<MT_, VT_> & precon, Index max_iters, typename VT_::DataType eps_relative, bool verbose = true)
      {
        // Each *_tilde vector is a correction, all others are defects.
        typedef typename VT_::DataType DT_;

        DT_ defnorm, defnorm_0, defnorm_00(DT_(1e14));
        Index iter = 0;
        DT_ rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
        //bool early_exit = 0;
        bool restarted = false;
        bool converged = 0;

        VT_ r(A.create_vector_l());
        VT_ r_tilde(A.create_vector_l());
        VT_ r_tilde_0(A.create_vector_l());
        VT_ p_tilde(A.create_vector_l());
        VT_ v(A.create_vector_l());
        VT_ v_tilde(A.create_vector_l());
        VT_ s(A.create_vector_l());
        VT_ s_tilde(A.create_vector_l());
        VT_ t(A.create_vector_l());
        VT_ t_tilde(A.create_vector_l());

        do
        {
          //r.defect(b, A, x);
          A.apply(r, x, b, -DT_(1));
          defnorm_0 = r.norm2();
          defnorm = defnorm_0;
          precon.apply(r_tilde_0, r);

          if (restarted == false)
          {
            defnorm_00 = defnorm_0;
          }
          r_tilde.copy(r_tilde_0);
          p_tilde.copy(r_tilde_0);

          rho_tilde = r_tilde_0.dot(r_tilde_0);

          // main BiCGStab loop
          do
          {
            iter = iter + 1;

            //v.product_matvec(A, p_tilde);
            A.apply(v, p_tilde);
            precon.apply(v_tilde, v);

            gamma_tilde = v_tilde.dot(r_tilde_0);

            if (Math::abs(gamma_tilde) < Math::abs(rho_tilde)*1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 1" << std::endl;
              break;
            }

            alpha_tilde = rho_tilde / gamma_tilde;

            if ((Math::abs(alpha_tilde) * v_tilde.norm2()) / defnorm < 1e-5)
            {
              restarted = true;;
              //std::cout << "Breakpoint 2" << std::endl;
              // \TODO warum ist das break hier nicht aktiv?
              //break;
            }

            DT_ malpha_tilde(-alpha_tilde);
            s.axpy(v, r, malpha_tilde);

            defnorm = s.norm2();
            if (defnorm < eps_relative * defnorm_00)
            {
              x.axpy(p_tilde, x, alpha_tilde);

              //early_exit = 1;
              converged = 1;
              //std::cout << "Breakpoint 3 (converged)" << std::endl;
              break;
            }
            s_tilde.axpy(v_tilde, r_tilde, malpha_tilde);

            //t.product_matvec(A, s_tilde);
            A.apply(t, s_tilde);

            precon.apply(t_tilde, t);

            gamma_tilde = t_tilde.dot(t_tilde);
            omega_tilde = t_tilde.dot(s_tilde);

            if (Math::abs(gamma_tilde) < Math::abs(omega_tilde) * 1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 4" << std::endl;
              break;
            }
            omega_tilde = omega_tilde / gamma_tilde;

            x.axpy(s_tilde, x, omega_tilde);
            x.axpy(p_tilde, x, alpha_tilde);

            DT_ momega_tilde(-omega_tilde);
            r.axpy(t, s, momega_tilde);

            defnorm = r.norm2();
            if (defnorm < eps_relative * defnorm_00)
            {
              converged = 1;
              //std::cout << "Breakpoint 5 (converged)" << std::endl;
              break;
            }

            r_tilde.axpy(t_tilde, s_tilde, momega_tilde);

            rho_tilde_old = rho_tilde;
            rho_tilde = r_tilde.dot(r_tilde_0);

            beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

            p_tilde.axpy(v_tilde, p_tilde, momega_tilde);
            p_tilde.scale(p_tilde, beta_tilde);
            p_tilde.axpy(p_tilde, r_tilde);

          } while (iter <= max_iters);

        } while (converged == 0 && iter < max_iters);


        Index used_iters = 0;
        used_iters = iter + 1;
        if (verbose)
          std::cout<<"Initial defect norm: " << defnorm_00 << ", Final defect norm: " << defnorm << ", used iters: " << used_iters << std::endl;
      }

      /**
       * \brief Variant with filter
       *
       * Copy-pasted from the original above
       *
       * \author Jordi Paul
       **/
      template <typename MT_, typename VT_, typename FT_>
      static void value(VT_ & x, const MT_ & A, const VT_ & b, const FT_ & filter, Preconditioner<MT_, VT_> & precon, Index max_iters, typename VT_::DataType eps_relative, bool verbose = true)
      {
        // Each *_tilde vector is a correction, all others are defects.
        typedef typename VT_::DataType DT_;

        DT_ defnorm, defnorm_0, defnorm_00(DT_(1e14));
        Index iter = 0;
        DT_ rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
        //bool early_exit = 0;
        bool restarted = false;
        bool converged = 0;

        VT_ r(A.create_vector_l());
        VT_ r_tilde(A.create_vector_l());
        VT_ r_tilde_0(A.create_vector_l());
        VT_ p_tilde(A.create_vector_l());
        VT_ v(A.create_vector_l());
        VT_ v_tilde(A.create_vector_l());
        VT_ s(A.create_vector_l());
        VT_ s_tilde(A.create_vector_l());
        VT_ t(A.create_vector_l());
        VT_ t_tilde(A.create_vector_l());

        do
        {
          //rdefect(b, A, x);
          Aapply(r, x, b, -DT_(1));
          filterfilter_def(r);

          defnorm_0 = r.norm2();
          defnorm = defnorm_0;
          precon.apply(r_tilde_0, r);
          filterfilter_cor(r_tilde_0);

          if (restarted == false)
          {
            defnorm_00 = defnorm_0;
          }
          r_tilde.copy(r_tilde_0);
          p_tilde.copy(r_tilde_0);

          rho_tilde = r_tilde_0dot(r_tilde_0);

          // main BiCGStab loop
          do
          {
            iter = iter + 1;

            //vproduct_matvec(A, p_tilde);
            Aapply(v, p_tilde);
            filterfilter_def(v);
            precon.apply(v_tilde, v);
            filterfilter_cor(v_tilde);

            gamma_tilde = v_tildedot(r_tilde_0);

            if (Math::abs(gamma_tilde) < Math::abs(rho_tilde)*1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 1" << std::endl;
              break;
            }

            alpha_tilde = rho_tilde / gamma_tilde;

            if ((Math::abs(alpha_tilde) * v_tilde.norm2()) / defnorm < 1e-5)
            {
              restarted = true;;
              //std::cout << "Breakpoint 2" << std::endl;
              // \TODO warum ist das break hier nicht aktiv?
              //break;
            }

            DT_ malpha_tilde(-alpha_tilde);
            saxpy(v, r, malpha_tilde);

            defnorm = s.norm2();
            if (defnorm < eps_relative * defnorm_00)
            {
              xaxpy(p_tilde, x, alpha_tilde);

              //early_exit = 1;
              converged = 1;
              //std::cout << "Breakpoint 3 (converged)" << std::endl;
              break;
            }
            s_tildeaxpy(v_tilde, r_tilde, malpha_tilde);

            //tproduct_matvec(A, s_tilde);
            Aapply(t, s_tilde);
            filterfilter_def(t);

            precon.apply(t_tilde, t);
            filterfilter_cor(t_tilde);

            gamma_tilde = t_tildedot(t_tilde);
            omega_tilde = t_tildedot(s_tilde);

            if (Math::abs(gamma_tilde) < Math::abs(omega_tilde) * 1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 4" << std::endl;
              break;
            }
            omega_tilde = omega_tilde / gamma_tilde;

            xaxpy(s_tilde, x, omega_tilde);
            xaxpy(p_tilde, x, alpha_tilde);

            DT_ momega_tilde(-omega_tilde);
            raxpy(t, s, momega_tilde);

            defnorm = r.norm2();
            if (defnorm < eps_relative * defnorm_00)
            {
              converged = 1;
              //std::cout << "Breakpoint 5 (converged)" << std::endl;
              break;
            }

            r_tildeaxpy(t_tilde, s_tilde, momega_tilde);

            rho_tilde_old = rho_tilde;
            rho_tilde = r_tildedot(r_tilde_0);

            beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

            p_tildeaxpy(v_tilde, p_tilde, momega_tilde);
            p_tildescale(p_tilde, beta_tilde);
            p_tildeaxpy(p_tilde, r_tilde);

          } while (iter <= max_iters);

        } while (converged == 0 && iter < max_iters);


        Index used_iters = 0;
        used_iters = iter + 1;
        if (verbose)
          std::cout<<"Initial defect norm: " << defnorm_00 << ", Final defect norm: " << defnorm << ", used iters: " << used_iters << std::endl;
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_BICGSTAB_HPP
