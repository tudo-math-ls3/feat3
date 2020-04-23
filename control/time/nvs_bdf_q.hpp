// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_TIME_NVS_BDF_Q_HPP
#define CONTROL_TIME_NVS_BDF_Q_HPP 1
#include <kernel/base_header.hpp>
#include <control/time/base.hpp>

#include "vector"

namespace FEAT
{
  namespace Control
  {
    namespace Time
    {
      /**
       * \brief Set of coefficients for operator splitting time discretisations based on BDF(q)
       *
       * \tparam DT_
       * The floating point precision
       *
       * This sounds more sophisticated than it is, as BDF(1) is the Implicit Euler scheme and BDF(2) is the only
       * other scheme supported until this class is extended to support higher order.
       *
       * For the incompressible Navier-Stokes equations
       * \f{align*}{
       *   \frac{\partial u}{\partial t} + (u \cdot \nabla ) u
       *   - \nabla \cdot \left(\frac{1}{Re} \nabla u - p I_d \right) & = f & \mathrm{in~} \Omega,
       *   u_{\partial \Omega = u_D} \\
       *   \nabla \cdot u & = 0 & \mathrm{in~} \Omega
       * \f}
       * with appropriate boundary conditions, a family of operator splitting techniques called
       * <em> projection methods </em> exists, where the equations are split into several, easier to handle parts
       * which are then solved successively. This usually introduces a <em> splitting </em> error.
       *
       * As a first step, replace the time derivative \f$\frac{\partial u}{\partial t}\f$ by a discrete approximation
       * \f$ \frac{D }{\delta_t}u \f$. For this class, the <em> backward differentiation formula </em> of order
       * \f$ q \f$ (BDF(q)) is used, where only \f$ q=1,2\f$ are implemented here.
       *
       * \f{align*}{
       *   \frac{D }{\delta_t} u^k & = \frac{1}{\delta_t} (\beta_q u^{k} - \sum_{k=0}^{q-1} \beta_k u^{k-1-q} ) \\
       *   & =
       *   \begin{cases}
       *     \frac{1}{\delta_t} (u^k - u^{k-1}), & q = 1 \\
       *     \frac{1}{\delta_t} \left(\frac{3}{2} u^k - (2 u^{k-1} - \frac{1}{2}u^{k-2}) \right), & q = 2
       *   \end{cases}
       * \f}
       *
       * The scheme splits the solution process into computing a <em> tentative velocity </em> \f$ \tilde{u} \f$
       * which is not divergence-free but satisfies the boundary condition \f$ \tilde{u}^* = u_D \f$ on
       * \f$ \partial \Omega\f$. It is then  projected to the space of (weakly) divergence free functions to obtain
       * the solution \f$u\f$, which then only satisfies the boundary condition \f$ u \cdot \nu = u_D \cdot \nu \f$
       * on \f$ \partial \Omega \f$.
       *
       * It is possible to skip computing \f$ u \f$ explicitly. For the convection term, extrapolation in time can be
       * used for linearisation and using \f$ u \f$ over \f$ \tilde{u} \f$ does not have any advantages. The only
       * other term using \f$ u \f$ is then the time discretisation, where the projection update can be carried out
       * implicitly by using the auxillary variable \f$ \phi \f$. Assuming that
       * \f$ \nabla \cdot (u \cdot \nabla) u = 0 \f$ and \f$ \nabla \cdot \Delta u =0 \f$, applying the divergence
       * operator to the momentum equation yields
       * \f{align*}{
       *   \nabla \cdot \frac{D }{\delta_t} u^k - \Delta p & = 0
       * \f}
       *
       * To linearise the convection term, the extrapolated velocity
       * \f$ u^{*,k} = 2\tilde{u}^{k-1} - \tilde{u}^{k-2}\f$
       * is used.
       *
       * The pressure in the momentum equation is also extrapolated with \f$ r \f$ from the previous time steps:
       * \f[
       *   p^{*,k} =
       *   \begin{cases}
       *     0, & r = 0 \\
       *     p^{k-1}, & r = 1 \\
       *     2 p^{k-1} - p^{k-2}, & r = 2
       *   \end{cases}
       * \f]
       *
       * The scheme is then as follow:
       * \f{align*}{
       *   \frac{D }{\delta_t} \tilde{u}^k + ( u^* \nabla ) u^k - \nabla \cdot
       *   \left( \nabla \tilde{u}^k - p^{*,k} I_d - \sum_{j=0}^{q-1} \frac{\beta_j}{\beta_q} \phi^{k-1-j} I_d \right) & = f^k,
       *   \tilde{u}^k_{| \partial \Omega} = u^k_D \\
       *   - \Delta \phi^k & = - \frac{\beta_q}{\delta_t} \nabla \cdot \tilde{u}^k,
       *   \partial_\nu \phi^k_{| \partial \Omega} = 0, \int_{\Omega} \phi^k dx = 0 \\
       *  p^k & = p^{k-1} + \phi^k - \frac{\chi}{Re} \nabla \cdot \tilde{u}^k
       *   \f}
       * The parameter \f$ \chi = 1 \f$ means the (better) rotational version of the scheme is used, otherwise
       * \f$ \chi = 0 \f$.
       *
       * Some special cases are the Chorin-Temam (or non-incremental first order pressure correction) scheme
       * ( \f$ (q,r) = (1,0) \f$ ), the Van Kan (or incremental first order pressure correction) scheme
       * ( \f$ (q,r) = (1,1) \f$ ).
       *
       * If the deformation tensor based bilinear form is used for the viscous term (meaning \f$ D(u) \f$ instead
       * of \f$ \nabla u \f$), the rotational coefficient gains an additional factor 2.
       *
       * If do-nothing boundaries for the momentum equation exist, the auxillary variable \f$ \phi \f$ needs Robin
       * boundary conditions.
       *
       * If the domain moves with the ALE velocity \f$ w \f$, the convection term should read
       * \f[
       *   ((u^{*,k} - w) \cdot \nabla) \tilde{u}^k
       * \f]
       *
       * For details (for the Stokes problem) see \cite GMS06.
       *
       * \author Jordi Paul
       *
       */
       template<typename DT_>
       class NvsBdfQ
       {
         public:
           /// The floating point precision
           typedef DT_ DataType;

         private:
           /// The k in BDF(k), has to be 1 or 2
           Index _num_steps;
           /// If k > 1, the methods performs this many steps as BDF(1) before switching to BDF(k)
           Index _num_steps_startup;
           /// How many backward steps should be taken for extrapolating the pressure
           Index _p_extrapolation_steps;
           /// Same in the startup phase
           Index _p_extrapolation_steps_startup;
           /// The Reynolds number
           const DataType _reynolds;
           /// The time step size
           DataType _delta_t;
           /// Shall we use the deformation tensor based bilinear form for the momentum equation?
           bool _use_deformation;
           /// Shall we use the rotational form of the operator splitting?
           bool _use_rotational_form;
           /// Are we still in the startup phase?
           bool _startup;

         public:
           /// 3 coefficients for the reaction, convection and diffusion terms for the left hand side of the momentum
           /// equation
           std::vector<DataType> coeff_lhs_v;
           /// 3 coefficients for the reaction, convection and diffusion terms for the right hand side of the momentum
           /// equation depending on the last computed velocity
           std::vector<DataType> coeff_rhs_v_1;
           /// k coeffcients for the auxillary variable phi enforcing the weak-divergence free condition
           std::vector<DataType> coeff_rhs_v_phi;
           /// \f$ n_p \f$ coefficients for extrapolating the pressure from previous time steps
           std::vector<DataType> coeff_extrapolation_p;

           /// Coefficient for the rhs for the momentum equation
           DataType coeff_rhs_f;
           /// Coefficient for the div(u) term for the projection step
           DataType coeff_rhs_proj_D_v;
           /// Coefficient the auxillary Poisson problem
           DataType coeff_rhs_phi;
           /// Coefficient for the reaction term depending on u[k-2]
           DataType coeff_rhs_v_2;
           /// Coefficient for the pressure gradient for the rhs of the momentum equation
           DataType coeff_rhs_v_p;

           /// How to handle the ALE term
           TermHandling ale_handling;
           /// How to handle the convection term
           TermHandling conv_handling;
           /// How to handle the viscous term
           TermHandling visc_handling;

           /**
            * \brief Constructor using a PropertyMap
            *
            * \param[in] config_section
            * The appropriate configuration section
            *
            * \param[in] reynolds
            * The Reynolds number
            *
            * \param[in] use_deformation
            * Use the deformation tensor based bilinear form for the viscous term?
            *
            * \param[in] startup
            * Set up as BDF(1) until finish_startup() is called?
            */
           NvsBdfQ(PropertyMap* config_section, const DataType reynolds, const bool use_deformation,
           const bool startup):
             _num_steps(~Index(0)),
             _num_steps_startup(1),
             _p_extrapolation_steps(~Index(0)),
             _p_extrapolation_steps_startup(~Index(0)),
             _reynolds(reynolds),
             _delta_t(0),
             _use_deformation(use_deformation),
             _use_rotational_form(false),
             _startup(startup),
             coeff_lhs_v(3),
             coeff_rhs_v_1(3),
             coeff_rhs_v_phi(0),
             coeff_extrapolation_p(0),
             coeff_rhs_f(0),
             coeff_rhs_proj_D_v(0),
             coeff_rhs_phi(0),
             coeff_rhs_v_2(0),
             coeff_rhs_v_p(0),
             ale_handling(TermHandling::off),
             conv_handling(TermHandling::off),
             visc_handling(TermHandling::off)
             {
               XASSERT(config_section != nullptr);
               // Get time discretisation order
               auto num_steps_p = config_section->query("num_steps");
               XASSERTM(num_steps_p.second, "TimeDiscretisation section is missing the mandatory num_steps entry!");
               _num_steps = Index(std::stoul(num_steps_p.first));
               XASSERTM(_num_steps == Index(1) || _num_steps == Index(2),"Only BDF1 and BDF2 are implemented!");

               // Get the order for extrapolation of the pressure in time
               auto p_extrapolation_p = config_section->query("p_extrapolation_steps");
               XASSERTM(p_extrapolation_p.second,
               "TimeDiscretisation section is missing the mandatory p_extrapolation_steps entry!");
               _p_extrapolation_steps = Index(std::stoul(p_extrapolation_p.first));
               XASSERTM(_p_extrapolation_steps < Index(3),"p extrapolation is only implemented for order 0, 1, 2!");
               _p_extrapolation_steps_startup = Math::min(_p_extrapolation_steps, Index(1));

               // Use the rotational form?
               auto use_rotational_form_p = config_section->query("use_rotational_form");
               XASSERTM(use_rotational_form_p.second,
               "TimeDiscretisation section is missing the mandatory use_rotational_form entry!");
               XASSERTM(std::stoi(use_rotational_form_p.first) == 0 || std::stoi(use_rotational_form_p.first) == 1,
               "use_rotational has to be set to 0 or 1");
               _use_rotational_form = std::stoi(use_rotational_form_p.first) == 1;

               // Get timestep size
               auto delta_t_p = config_section->query("delta_t");
               XASSERTM(delta_t_p.second, "TimeDiscretisation section is missing the mandatory delta_t entry!");
               _delta_t = DataType(std::stod(delta_t_p.first));
               XASSERT(_delta_t > DataType(0));

               // Use ale (meaning solve Navier-Stokes instead of plain Stokes)?
               auto ale_p = config_section->query("ALE");
               XASSERTM(ale_p.second,
               "TimeDiscretisation section is missing the mandatory ALE entry!");
               ale_p.first.parse(ale_handling);

               // Use convection (meaning solve Navier-Stokes instead of plain Stokes)?
               auto convection_p = config_section->query("convection");
               XASSERTM(convection_p.second,
               "TimeDiscretisation section is missing the mandatory convection entry!");
               convection_p.first.parse(conv_handling);

               // Use viscous (meaning solve Navier-Stokes instead of plain Stokes)?
               auto viscous_p = config_section->query("viscous");
               XASSERTM(viscous_p.second,
               "TimeDiscretisation section is missing the mandatory viscous entry!");
               viscous_p.first.parse(visc_handling);

               coeff_rhs_v_phi.clear();
               coeff_rhs_v_phi = std::vector<DataType>(_num_steps, DataType(0));

               coeff_extrapolation_p.clear();
               coeff_extrapolation_p = std::vector<DataType>(_p_extrapolation_steps, DataType(0));

               if(startup)
               {
                 set_coeffs(_num_steps_startup, _p_extrapolation_steps_startup);
               }
               else
               {
                 set_coeffs(_num_steps, _p_extrapolation_steps);
               }
             }

           /**
            * \brief Virtual destructor
            */
           virtual ~NvsBdfQ()
           {
           }

           /**
            * \brief Returns the time step size
            *
            * \returns The time step size
            */
           DataType delta_t()
           {
             return _delta_t;
           }

           /**
            * \brief Are we using the rotational form?
            *
            * \returns \c true if the rotational form is used.
            */
           bool is_rotational()
           {
             return _use_rotational_form;
           }

           /**
            * \brief Returns the current number of steps
            *
            * This is different from the maximum number of steps in the startup phase!
            *
            * \returns The current number of steps
            */
           Index get_num_steps()
           {
             return _startup ? _num_steps_startup : _num_steps;
           }

           /**
            * \brief Returns the maximum number of steps
            *
            * This is different from the current number of steps in the startup phase!
            *
            * \returns The maximum number of  steps
            */
           Index get_max_num_steps()
           {
             return Math::max(_num_steps, _num_steps_startup);
           }

           /**
            * \brief Returns the maximum number of pressure extrapolation steps
            *
            * This is different from the current number of pressure extrapolation steps in the startup phase!
            *
            * \returns The maximum number of pressure extrapolation steps
            */
           Index get_max_p_extrapolation_steps()
           {
             return Math::max(_p_extrapolation_steps, _p_extrapolation_steps_startup);
           }

           /**
            * \brief Returns the current number of pressure extrapolation steps
            *
            * This is different from the maximum number of pressure extrapolation steps in the startup phase!
            *
            * \returns The current number of pressure extrapolation steps
            */
           Index get_p_extrapolation_steps()
           {
             return _startup ? _p_extrapolation_steps_startup : _p_extrapolation_steps;
           }

           /**
            * \brief Finishes the startup phase
            *
            * This will set the coefficients again.
            */
           void finish_startup()
           {
             _startup = false;

             set_coeffs(_num_steps, _p_extrapolation_steps);
           }

           /**
            * \brief Sets the coefficients
            *
            * \param[in] num_steps
            * The \f$ q \f$
            *
            * \param[in] p_extrapolation_steps
            * The \f$ r \f$
            */
           void set_coeffs(const Index num_steps, const Index p_extrapolation_steps)
           {

             if(num_steps == Index(1))
             {
               // Mass matrix
               coeff_lhs_v.at(0) = DataType(1);
               // Convection matrix
               coeff_lhs_v.at(1) = (conv_handling == TermHandling::impl || ale_handling == TermHandling::impl)
                 ? _delta_t : DataType(0);
               // Stiffness matrix
               coeff_lhs_v.at(2) = (visc_handling == TermHandling::impl) ? _delta_t/_reynolds : DataType(0);

               // Mass matrix
               coeff_rhs_v_1.at(0) = DataType(1);
               // Convection matrix
               coeff_rhs_v_1.at(1) = (conv_handling == TermHandling::expl || ale_handling == TermHandling::expl)
                 ? -_delta_t : DataType(0);
               // Stiffness matrix
               coeff_rhs_v_1.at(2) = (visc_handling == TermHandling::expl) ? -_delta_t/_reynolds : DataType(0);

               // Mass matrix coefficient for previous time step vector
               coeff_rhs_v_2 = DataType(0);

               coeff_rhs_f = _delta_t;

               // coefficient for -(p*[k], div phi) on the rhs for the velocity eq.
               coeff_rhs_v_p = -_delta_t;
               // coefficient for -(phi[k], div phi) on the rhs for the velocity eq.
               coeff_rhs_v_phi.at(0) = -_delta_t;

               coeff_rhs_phi = DataType(1)/_delta_t;

             }
             else if(num_steps == Index(2))
             {
               // Mass matrix
               coeff_lhs_v.at(0) = DataType(1);
               // Convection matrix
               coeff_lhs_v.at(1) = (conv_handling == TermHandling::impl) ? DataType(2)/DataType(3)*_delta_t : DataType(0);
               // Stiffness matrix
               coeff_lhs_v.at(2) = (visc_handling == TermHandling::impl || ale_handling == TermHandling::impl)
                 ? DataType(2)/DataType(3)*_delta_t/_reynolds : DataType(0);

               // Mass matrix
               coeff_rhs_v_1.at(0) = DataType(4)/DataType(3);
               // Convection matrix
               coeff_rhs_v_1.at(1) = (conv_handling == TermHandling::expl || ale_handling == TermHandling::expl)
                 ? -DataType(2)/DataType(3)*_delta_t : DataType(0);
               // Stiffness matrix
               coeff_rhs_v_1.at(2) =
                 (visc_handling == TermHandling::expl) ? -DataType(2)/DataType(3)*_delta_t/_reynolds : DataType(0);

               coeff_rhs_f = DataType(2)/DataType(3)*_delta_t;

               // Mass matrix coefficient for previous time step vector
               coeff_rhs_v_2 = -DataType(1)/DataType(3);

               // coefficient for -(p*[k], div phi) on the rhs for the velocity eq.
               coeff_rhs_v_p = -DataType(2)/DataType(3)*_delta_t;

               // coefficient for -(phi[k], div phi) on the rhs for the velocity eq.
               coeff_rhs_v_phi.at(0) =  -DataType(8)/DataType(9)*_delta_t;
               // coefficient for -(phi[k-1], div phi) on the rhs for the velocity eq.
               coeff_rhs_v_phi.at(1) = DataType(2)/DataType(9)*_delta_t;

               coeff_rhs_phi = DataType(3)/(DataType(2)*_delta_t);

             }
             else
             {
               XABORTM("Invalid number of BDF steps: "+stringify(num_steps));
             }

             // Coefficient for the div v term for the projection step
             coeff_rhs_proj_D_v = _use_rotational_form ? DataType(1)/_reynolds : DataType(0);

             if(_use_deformation)
             {
               coeff_rhs_proj_D_v *= DataType(2);
             }

             // Coefficients for the extrapolation of the pressure
             if(p_extrapolation_steps == Index(1))
             {
               coeff_extrapolation_p.at(0) = DataType(1);
             }
             else if(p_extrapolation_steps == Index(2))
             {
               coeff_extrapolation_p.at(0) = DataType(2);
               coeff_extrapolation_p.at(1) = -DataType(1);
             }
             else
             {
               XABORTM("Invalid number of pressure extrapolation steps: "+stringify(p_extrapolation_steps));
             }

           }
       };
    } // namespace Time
  } // namespace Control
} // name FEAT

#endif // CONTROL_TIME_NVS_BDF_Q_HPP
