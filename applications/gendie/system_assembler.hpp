// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once
#include <kernel/util/property_map.hpp>
#include <applications/gendie/materials.hpp>
#include <kernel/assembly/burgers_velo_material_assembly_job.hpp>
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/util/stop_watch.hpp>

#include "format_helper.hpp"

#include <memory>

namespace Gendie
{
  // CRTP interface for Flow System Assembler
  template<typename Derived_, typename Domain_>
  class SystemFlowAssemblerBaseCRTP
  {
  public:
    typedef Domain_ DomainType;
    typedef typename Domain_::WeightType SystemDataType;
    static constexpr int padlen = 30;
    static constexpr char pc = '.';
    mutable FEAT::StopWatch watch_total, watch_restrict, watch_restrict_alloc;
    DomainType* domain;
    SystemDataType stiffness_factor = SystemDataType(1);
    SystemDataType theta = SystemDataType(0);
    FEAT::PreferredBackend asm_backend = FEAT::PreferredBackend::generic;


  protected:
    /**
      * \brief Casts \c this to its true type
      *
      * \returns A Derived_ reference to \c this
      */
    Derived_& cast() {return static_cast<Derived_&>(*this);}

    /// \copydoc cast()
    const Derived_& cast() const {return static_cast<const Derived_&>(*this);}

    bool _parse([[maybe_unused]] const FEAT::PropertyMap* prop)
    {
      return true;
    }

    void _set_domain(DomainType& domain_)
    {
      domain = &domain_;
    }

    void _init()
    {

    }

    void _done()
    {

    }

    void _set_jacobian(double DOXY(weight))
    {
      XABORTM("Has to be implemented by derived class");
    }

    template<typename LevelType_, typename VecType_>
    void _assemble_matrices(std::deque<std::shared_ptr<LevelType_>>& DOXY(system_levels), const VecType_& DOXY(vec_sol)) const
    {
      XABORTM("Specialization has to be assembled by base class");
    }

    template<typename LevelType_>
    void _assemble_matrices_isotrop(std::deque<std::shared_ptr<LevelType_>>& DOXY(system_levels)) const
    {
      XABORTM("Specialization has to be assembled by base class");
    }

    void _set_theta(SystemDataType _theta)
    {
      this->theta = _theta;
    }

    template<typename LevelType_, typename VecType_>
    void _fill_convections(const std::deque<std::shared_ptr<LevelType_>>& system_levels, std::vector<typename LevelType_::GlobalVeloVector>& conv_vectors, const VecType_& vec_sol) const
    {
      this->watch_restrict_alloc.start();
      typedef typename LevelType_::GlobalVeloVector ConvVecType;
      constexpr bool same_vecs = std::is_same_v<VecType_, ConvVecType>;
      if constexpr (same_vecs)
      {
        conv_vectors.push_back(ConvVecType(system_levels.front()->gate_velo, vec_sol.local().template at<0>(), FEAT::LAFEM::CloneMode::Shallow));
      }
      else
      {
        conv_vectors.push_back(system_levels.front()->create_global_vector_velo());
        conv_vectors.back().local().convert(vec_sol.local().template at<0>());
      }
      for(std::size_t i = 1; i < system_levels.size(); ++i)
      {
        conv_vectors.push_back(system_levels.at(i)->create_global_vector_velo());
        conv_vectors.back().format();
      }
      this->watch_restrict_alloc.stop();
    }

    template<typename LevelType_>
    void _restrict_convections(const std::deque<std::shared_ptr<LevelType_>>& system_levels, std::vector<typename LevelType_::GlobalVeloVector>& conv_vectors) const
    {
      this->watch_restrict.start();
      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        // do we need to participate in the communication?
        if((i+1) >= this->domain->size_virtual())
          break;
        // does this process have another system level?
        if((i+1) < system_levels.size())
        {
          const auto& vec_fine = conv_vectors.at(i);
          // create a coarse mesh velocity vector
          auto& vec_crs = conv_vectors.at(i+1);

          // truncate fine mesh velocity vector
          system_levels.at(i)->transfer_velo.trunc(vec_fine, vec_crs);
        }
        else
        {
          const auto& vec_fine = conv_vectors.at(i);
          // this process is a child, so send truncation to parent
          system_levels.at(i)->transfer_velo.trunc_send(vec_fine);
        }
      }
      this->watch_restrict.stop();
    }

    void _set_backend() const
    {
      FEAT::Backend::set_preferred_backend(this->asm_backend);
    }

    FEAT::String _format_timings() const
    {
      return FEAT::String{};
    }

    public:

    SystemFlowAssemblerBaseCRTP() = default;

    SystemFlowAssemblerBaseCRTP(DomainType& domain_, const FEAT::PropertyMap* prop)
     : domain(&domain_)
    {
      _parse(prop);
    }

    virtual ~SystemFlowAssemblerBaseCRTP() = default;

    SystemFlowAssemblerBaseCRTP(const SystemFlowAssemblerBaseCRTP&) = default;
    SystemFlowAssemblerBaseCRTP(SystemFlowAssemblerBaseCRTP&&) noexcept = default;
    SystemFlowAssemblerBaseCRTP& operator=(const SystemFlowAssemblerBaseCRTP&) = default;
    SystemFlowAssemblerBaseCRTP& operator=(SystemFlowAssemblerBaseCRTP&&) noexcept = default;

    void init()
    {
      watch_total.start();
      this->cast()._init();
      watch_total.stop();
    }

    void done()
    {
      watch_total.start();
      this->cast()._done();
      watch_total.stop();
    }

    /**
     * \brief Assembles the matrices defined by the system level
     *
     * \tparam VecType_ A global vectortype convertable to the inner system vectortype
     *
     * \param[in] vec_sol The local solution vector used as convection for assembly
     */
    template<typename LevelType_, typename VecType_>
    void assemble_matrices(std::deque<std::shared_ptr<LevelType_>>& system_levels, const VecType_& vec_sol) const
    {
      set_backend();
      watch_total.start();
      this->cast()._assemble_matrices(system_levels, vec_sol);
      watch_total.stop();
    }

    template<typename LevelType_>
    void assemble_matrices_isotrop(std::deque<std::shared_ptr<LevelType_>>& system_levels) const
    {
      set_backend();
      watch_total.start();
      this->cast()._assemble_matrices_isotrop(system_levels);
      watch_total.stop();
    }

    auto& get_domain()
    {
      return *domain;
    }

    /**
     * \brief Sets assembler such that the jacobian (if available) of the underlying system is assembled, weighted by a given parameter
     *
     * \param[in] weight Assembles the jacobian with a specific weight. How this parameter is interpreted exactly is specific to the actual system being assembled
     *                   but in general, weight == 0 means the standard system is assembled, weight > 0 some form ob (Pseudo)-Newton is assembled
     */
    void set_jacobian(double weight)
    {
      this->cast()._set_jacobian(weight);
    }

    SystemDataType calc_stiffness_factor() const
    {
      this->cast()._calc_stiffness_factor();
    }

    void set_stiffness_factor(SystemDataType factor)
    {
      this->stiffness_factor = factor;
    }

    SystemDataType get_stiffness_factor() const
    {
      return this->stiffness_factor;
    }

    void set_theta(SystemDataType _theta)
    {
      return this->cast()._set_theta(_theta);
    }

    void set_backend() const
    {
      this->_set_backend();
    }

    FEAT::String format_timings() const
    {
      return this->cast()._format_timings();
    }

    double get_total_time_elapsed() const
    {
      return this->watch_total.elapsed();
    }

  }; //class SystemFlowAssemblerBaseCRTP

  /**
   * \brief Standard Burgers Carreau assembler
   *
   */
  template<typename Domain_>
  class CarreauSingleMatSystemFlowAssembler : public SystemFlowAssemblerBaseCRTP<CarreauSingleMatSystemFlowAssembler<Domain_>, Domain_>
  {
  public:
    typedef SystemFlowAssemblerBaseCRTP<CarreauSingleMatSystemFlowAssembler, Domain_> BaseClass;
    typedef typename BaseClass::DomainType DomainType;
    typedef typename BaseClass::SystemDataType SystemDataType;
    friend BaseClass;

    SystemDataType reg_eps = 1E-100;
    bool navier = false;

  protected:

    bool _parse([[maybe_unused]] const FEAT::PropertyMap* prop)
    {

      return true;
    }

    void _init()
    {
    }

    void _done()
    {
    }

    void _set_jacobian(double weight)
    {
      _weight_jacobian = weight;
    }

    SystemDataType _calc_stiffness_factor() const
    {
      return this->_material->get_viscosity_model()->get_data()[0] * SystemDataType(1E3) / _mesh_scale;
    }

    /// always assume that matrix structure is preassambled
    template<typename LevelType_, typename VecType_>
    void _assemble_matrices(std::deque<std::shared_ptr<LevelType_>>& system_levels, const VecType_& vec_sol) const
    {
      typedef typename LevelType_::GlobalVeloVector ConvVecType;
      typedef typename LevelType_::DataType AssemblyDataType;

      // constexpr bool same_vecs = std::is_same_v<VecType_, ConvVecType>;
      // assemble convections
      std::vector<ConvVecType> conv_vecs;
      conv_vecs.reserve(system_levels.size());
      this->_fill_convections(system_levels, conv_vecs, vec_sol);
      this->_restrict_convections(system_levels, conv_vecs);

      this->watch_assembly.start();
      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        if(i == 0) this->watch_assembly_fine.start();
        // setup
        auto& cur_sys = *system_levels[i];
        auto& cur_dom = *this->domain->at(i);
        ASSERTM(!cur_sys.matrix_a.local().empty(), "A Matrix not allocated");
        ASSERTM(!cur_sys.matrix_b.local().empty(), "B Matrix not allocated");
        ASSERTM(!cur_sys.matrix_d.local().empty(), "D Matrix not allocated");
        FEAT::Assembly::BurgersVeloMaterialBlockedMatrixAssemblyJob burgers_mat_job(
          cur_sys.matrix_a.local(), conv_vecs[i].local(), cur_dom.space_velo,
          this->_cubature, this->_material->get_visc_func(AssemblyDataType(_mesh_scale)), this->_material->get_visc_der_func(AssemblyDataType(_mesh_scale))
        );
        {
          const auto visc_model_data = this->_material->get_viscosity_model()->get_data();
          burgers_mat_job.deformation = true;
          burgers_mat_job.nu = AssemblyDataType(visc_model_data[0] * SystemDataType(1E3) / this->_mesh_scale);
          burgers_mat_job.aT = AssemblyDataType(this->_material->get_abiat_aT(this->_temperature));
          burgers_mat_job.sd_delta = AssemblyDataType(0);
          burgers_mat_job.sd_nu = AssemblyDataType(visc_model_data[0] * SystemDataType(1E3) / this->_mesh_scale);
          burgers_mat_job.beta = navier ? AssemblyDataType(_material->get_density_gram_per_unit(this->_mesh_scale)) : AssemblyDataType(0);
          burgers_mat_job.reg_eps = AssemblyDataType(reg_eps);
          burgers_mat_job.frechet_beta = AssemblyDataType(SystemDataType(burgers_mat_job.beta) * this->_weight_jacobian);
          burgers_mat_job.frechet_material = AssemblyDataType(SystemDataType(1) * this->_weight_jacobian);
          burgers_mat_job.theta = AssemblyDataType(this->theta);
          burgers_mat_job.alpha = AssemblyDataType(SystemDataType(1) / this->stiffness_factor);
        }

        cur_sys.matrix_a.local().format();

        // actual assembly
        cur_dom.domain_asm.assemble(burgers_mat_job);
        if(i == 0) this->watch_assembly_fine.stop();
      }
      this->watch_assembly.stop();
    }

    template<typename LevelType_>
    void _assemble_matrices_isotrop(std::deque<std::shared_ptr<LevelType_>>& system_levels) const
    {
      typedef typename LevelType_::DataType AssemblyDataType;
      this->watch_assembly.start();

      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        if(i == 0) this->watch_assembly_fine.start();
        // setup
        auto& cur_sys = *system_levels[i];
        auto& cur_dom = *this->domain->at(i);
        ASSERTM(!cur_sys.matrix_a.local().empty(), "A Matrix not allocated");
        ASSERTM(!cur_sys.matrix_b.local().empty(), "B Matrix not allocated");
        ASSERTM(!cur_sys.matrix_d.local().empty(), "D Matrix not allocated");
        Assembly::Common::LaplaceOperatorBlocked<LevelType_::dim> laplace_op;

        cur_sys.matrix_a.local().format();

        // actual assembly
        Assembly::assemble_bilinear_operator_matrix_1(cur_dom.domain_asm, cur_sys.matrix_a.local(), laplace_op, cur_dom.space_velo, this->_cubature,
          AssemblyDataType(this->_material->get_viscosity_model()->get_data()[0]) * AssemblyDataType(this->_material->get_viscosity_model()->get_visco_scaling_factor(AssemblyDataType(_mesh_scale))) *
                 AssemblyDataType(SystemDataType(1) / this->stiffness_factor));
        if(i == 0) this->watch_assembly_fine.stop();
      }
      this->watch_assembly.stop();
    }

    FEAT::String _format_timings() const
    {
      enum _TimeID
      {
        total = 0,
        assembly = 1,
        assembly_fine = 2,
        restrict = 3,
        restrict_alloc = 4,
        num_entries = 5
      };
      FEAT::String s;
      double timings_max[_TimeID::num_entries], timings_min[_TimeID::num_entries], timings[_TimeID::num_entries];
      timings_min[_TimeID::total] = timings_max[_TimeID::total] = timings[_TimeID::total] = this->watch_total.elapsed();
      timings_min[_TimeID::assembly] = timings_max[_TimeID::assembly] = timings[_TimeID::assembly] = this->watch_assembly.elapsed();
      timings_min[_TimeID::assembly_fine] = timings_max[_TimeID::assembly_fine] = timings[_TimeID::assembly_fine] = this->watch_assembly_fine.elapsed();
      timings_min[_TimeID::restrict] = timings_max[_TimeID::restrict] = timings[_TimeID::restrict] = this->watch_restrict.elapsed();
      timings_min[_TimeID::restrict_alloc] = timings_max[_TimeID::restrict_alloc] = timings[_TimeID::restrict_alloc] = this->watch_restrict_alloc.elapsed();
      // sync our timings
      this->domain->comm().allreduce(timings_max, timings_max, _TimeID::num_entries, FEAT::Dist::op_max);
      this->domain->comm().allreduce(timings_min, timings_min, _TimeID::num_entries, FEAT::Dist::op_min);

      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Timings Single Material System Assembler----------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += format_subtime_mm("Allocate Restrictions", timings[_TimeID::restrict_alloc], timings[_TimeID::total], timings_min[_TimeID::restrict_alloc], timings_max[_TimeID::restrict_alloc], this->padlen);
      s += format_subtime_mm("Restrict Convection", timings[_TimeID::restrict], timings[_TimeID::total], timings_min[_TimeID::restrict], timings_max[_TimeID::restrict], this->padlen);
      s += format_subtime_mm("Assemble Systems", timings[_TimeID::assembly], timings[_TimeID::total], timings_min[_TimeID::assembly], timings_max[_TimeID::assembly], this->padlen);
      s += format_subtime_mm("Assemble Fine System", timings[_TimeID::assembly_fine], timings[_TimeID::total], timings_min[_TimeID::assembly_fine], timings_max[_TimeID::assembly_fine], this->padlen);
      s += format_subtime_mm("System Assembly Total Time", timings[_TimeID::total], timings[_TimeID::total], timings_min[_TimeID::total], timings_max[_TimeID::total], this->padlen);

      return s;
    }

    const Material<SystemDataType>* _material;
    mutable FEAT::StopWatch watch_assembly, watch_assembly_fine;
    String _cubature = "gauss-legendre:5";
    SystemDataType _temperature;
    SystemDataType _mesh_scale;
    SystemDataType _weight_jacobian;


  public:
    CarreauSingleMatSystemFlowAssembler() = default;

    CarreauSingleMatSystemFlowAssembler(DomainType& domain_, const FEAT::PropertyMap* prop, const Material<SystemDataType>& material, SystemDataType temperature, SystemDataType mesh_scale)
      : BaseClass(domain_, prop),
      _material(&material),
      _temperature(temperature),
      _mesh_scale(mesh_scale),
      _weight_jacobian(SystemDataType(0))
    {
      _parse(prop);
    }

    virtual ~CarreauSingleMatSystemFlowAssembler() = default;

    CarreauSingleMatSystemFlowAssembler(const CarreauSingleMatSystemFlowAssembler&) = default;
    CarreauSingleMatSystemFlowAssembler(CarreauSingleMatSystemFlowAssembler&&) noexcept = default;
    CarreauSingleMatSystemFlowAssembler& operator=(const CarreauSingleMatSystemFlowAssembler&) = default;
    CarreauSingleMatSystemFlowAssembler& operator=(CarreauSingleMatSystemFlowAssembler&&) noexcept = default;

    void set_material(const Material<SystemDataType>& material, SystemDataType mesh_scale)
    {
      _material = &material;
      _mesh_scale = mesh_scale;
    }

    void set_cubature(const String& cubature)
    {
      _cubature = cubature;
    }

  }; //class CarrauSystemFlowAssembler<...>
}
