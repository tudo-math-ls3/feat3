// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once
#include <kernel/util/property_map.hpp>
#include <applications/gendie/materials.hpp>
#include <kernel/assembly/burgers_velo_material_assembly_job.hpp>
#include <kernel/util/stop_watch.hpp>
#include "format_helper.hpp"

#include <memory>

namespace Gendie
{
  // CRTP interface for Flow System Assembler
  template<typename Derived_, typename Domain_, typename SystemLevel_>
  class DefectFlowAssemblerBaseCRTP
  {
  public:
    typedef Domain_ DomainType;
    typedef SystemLevel_ LevelType;
    typedef typename LevelType::GlobalSystemVector VectorType;
    typedef typename LevelType::GlobalVeloVector ConvVectorType;
    typedef typename VectorType::DataType DataType;
    static constexpr int padlen = 30;
    static constexpr char pc = '.';
    mutable FEAT::StopWatch watch_total, watch_init, watch_create_vector;
    DomainType* domain;
    std::shared_ptr<LevelType> system_level;
    DataType stiffness_factor = DataType(1);
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

    void _set_system(const std::shared_ptr<LevelType> system_level_)
    {
      system_level = system_level_;
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

    void _assemble_defect(VectorType& DOXY(vec_def), const VectorType& DOXY(vec_primal), const VectorType& DOXY(vec_conv), const VectorType& DOXY(vec_rhs)) const
    {
      XABORTM("Specialization has to be implemented by derived class");
    }

    VectorType _create_vector() const
    {
      XABORTM("To be implemented by derived class");
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

    DefectFlowAssemblerBaseCRTP() = default;

    DefectFlowAssemblerBaseCRTP(DomainType& domain_, const std::shared_ptr<LevelType>& system_level_, const FEAT::PropertyMap* prop)
     : domain(&domain_),
       system_level(system_level_)
    {
      _parse(prop);
    }

    virtual ~DefectFlowAssemblerBaseCRTP() = default;

    DefectFlowAssemblerBaseCRTP(const DefectFlowAssemblerBaseCRTP&) = default;
    DefectFlowAssemblerBaseCRTP(DefectFlowAssemblerBaseCRTP&&) noexcept = default;
    DefectFlowAssemblerBaseCRTP& operator=(const DefectFlowAssemblerBaseCRTP&) = default;
    DefectFlowAssemblerBaseCRTP& operator=(DefectFlowAssemblerBaseCRTP&&) noexcept = default;

    void init()
    {
      watch_total.start();
      watch_init.start();
      this->cast()._init();
      watch_init.stop();
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
     * \param[out] vec_def The defect vector to be assembled. Values are ignored.
     * \param[in] vec_primal The synced 1 primal vector that is integrated.
     * \param[in] vec_conv The synced 1 convection vector used.
     * \param[in] vec_rhs The synced 1 rhs vector.
     */
    void assemble_defect(VectorType& vec_def, const VectorType& vec_primal, const VectorType& vec_conv, const VectorType& vec_rhs) const
    {
      set_backend();
      watch_total.start();
      this->cast()._assemble_defect(vec_def, vec_primal, vec_conv, vec_rhs);
      watch_total.stop();
    }

    void assemble_unsynced_defect(VectorType& vec_def, const VectorType& vec_primal, const VectorType& vec_conv) const
    {
      set_backend();
      watch_total.start();
      this->cast()._assemble_unsynced_defect(vec_def, vec_primal, vec_conv);
      watch_total.stop();
    }

    auto& get_system_level()
    {
      return system_level;
    }

    auto& get_domain()
    {
      return *domain;
    }

    VectorType create_vector() const
    {
      watch_total.start();
      watch_create_vector.start();
      auto tmp = this->cast()._create_vector();
      watch_create_vector.stop();
      watch_total.stop();
      return tmp;
    }

    void set_system(const std::shared_ptr<LevelType>& system_level_)
    {
      this->_set_system(system_level_);
    }

    void set_backend() const
    {
      this->_set_backend();
    }

    DataType calc_stiffness_factor() const
    {
      return this->cast()._calc_stiffness_factor();
    }

    void set_stiffness_factor(DataType factor)
    {
      this->stiffness_factor = factor;
    }

    DataType get_stiffness_factor() const
    {
      return this->stiffness_factor;
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
  template<typename Domain_, typename SystemLevel_, bool asm_fbm_defect_ = true>
  class CarreauSingleMatDefectFlowAssembler : public DefectFlowAssemblerBaseCRTP<CarreauSingleMatDefectFlowAssembler<Domain_, SystemLevel_>, Domain_, SystemLevel_>
  {
  public:
    typedef DefectFlowAssemblerBaseCRTP<CarreauSingleMatDefectFlowAssembler, Domain_, SystemLevel_> BaseClass;
    typedef typename BaseClass::DomainType DomainType;
    typedef typename BaseClass::LevelType LevelType;
    typedef typename BaseClass::VectorType VectorType;
    typedef typename BaseClass::DataType DataType;
    friend BaseClass;

    static constexpr bool assemble_fbm_defect = asm_fbm_defect_;
    mutable FEAT::StopWatch watch_assembly, watch_fbm_filter, watch_grad_apply, watch_sync_and_rhs;

    DataType reg_eps = 1E-100;
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

    VectorType _create_vector() const
    {
      return this->system_level->create_global_vector_sys();
    }

    DataType _calc_stiffness_factor() const
    {
      return this->_material->get_viscosity_model()->get_data()[0] * DataType(1E3) / _mesh_scale;
    }

    //the front system has to have B and D matrix correctly assembled
    void _assemble_defect(VectorType& vec_def, const VectorType& vec_primal, const VectorType& vec_conv, const VectorType& vec_rhs) const
    {
      watch_assembly.start();
      // setup
      const auto& cur_sys = *this->system_level;
      auto& cur_dom = *this->domain->front();
      FEAT::Assembly::BurgersVeloMaterialBlockedVectorAssemblyJob burgers_vec_job(
        vec_def.local().template at<0>(), vec_primal.local().template at<0>(), vec_conv.local().template at<0>(), cur_dom.space_velo,
        _cubature, this->_material->get_visc_func(_mesh_scale), this->_material->get_visc_der_func(_mesh_scale)
      );
      {
        const auto visc_model_data = this->_material->get_viscosity_model()->get_data();
        burgers_vec_job.deformation = true;
        burgers_vec_job.nu = visc_model_data[0] * DataType(1E3) / _mesh_scale;
        burgers_vec_job.aT = this->_material->get_abiat_aT(this->_temperature);
        burgers_vec_job.beta = navier ? this->_material->get_density_gram_per_unit(_mesh_scale) : DataType(0);
        burgers_vec_job.reg_eps = reg_eps;
        burgers_vec_job.alpha = DataType(-1) / this->stiffness_factor;
      }
      vec_def.format();
      // actual assembly
      cur_dom.domain_asm.assemble(burgers_vec_job);
      watch_assembly.stop();
      watch_fbm_filter.start();
      if constexpr(assemble_fbm_defect)
        cur_sys.apply_fbm_filter_to_def(vec_def, vec_primal, -DataType(1)/*/this->stiffness_factor*/);
      watch_fbm_filter.stop();

      watch_grad_apply.start();
      // remaining defect assembly
      cur_sys.matrix_sys.local().block_b().apply(
        vec_def.local().template at<0>(), vec_primal.local().template at<1>(), vec_def.local().template at<0>(), DataType(-1));
      cur_sys.matrix_sys.local().block_d().apply(
        vec_def.local().template at<1>(), vec_primal.local().template at<0>(), vec_def.local().template at<1>(), DataType(-1));
      watch_grad_apply.stop();
      // in theory, one could save the unsynced/unfiltered defect vector now, which we will ignore here
      // synchronize
      watch_sync_and_rhs.start();
      vec_def.sync_0();
      // and now add our rhs, if not empty
      if(!vec_rhs.local().size())
      {
        vec_def.local().template at<0>().axpy(vec_rhs.local().template at<0>(), DataType(1)/this->stiffness_factor);
        vec_def.local().template at<1>().axpy(vec_rhs.local().template at<1>());
        // vec_def.axpy(vec_rhs);
      }
      watch_sync_and_rhs.stop();

    }

    //the front system has to have B and D matrix correctly assembled
    void _assemble_unsynced_defect(VectorType& vec_def, const VectorType& vec_primal, const VectorType& vec_conv) const
    {
      watch_assembly.start();
      // setup
      const auto& cur_sys = *this->system_level;
      auto& cur_dom = *this->domain->front();
      FEAT::Assembly::BurgersVeloMaterialBlockedVectorAssemblyJob burgers_vec_job(
        vec_def.local().template at<0>(), vec_primal.local().template at<0>(), vec_conv.local().template at<0>(), cur_dom.space_velo,
        _cubature, this->_material->get_visc_func(_mesh_scale), this->_material->get_visc_der_func(_mesh_scale)
      );
      {
        const auto visc_model_data = this->_material->get_viscosity_model()->get_data();
        burgers_vec_job.deformation = true;
        burgers_vec_job.nu = visc_model_data[0] * DataType(1E3) / _mesh_scale;
        burgers_vec_job.aT = this->_material->get_abiat_aT(this->_temperature);
        burgers_vec_job.beta = navier ? this->_material->get_density_gram_per_unit(_mesh_scale) : DataType(0);
        burgers_vec_job.reg_eps = reg_eps;
        burgers_vec_job.alpha = DataType(-1) / this->stiffness_factor;
      }
      vec_def.format();
      // actual assembly
      cur_dom.domain_asm.assemble(burgers_vec_job);
      watch_assembly.stop();
      watch_fbm_filter.start();
      if constexpr(assemble_fbm_defect)
        cur_sys.apply_fbm_filter_to_def(vec_def, vec_primal, -DataType(1));
      watch_fbm_filter.stop();

      watch_grad_apply.start();
      // remaining defect assembly
      cur_sys.matrix_sys.local().block_b().apply(
        vec_def.local().template at<0>(), vec_primal.local().template at<1>(), vec_def.local().template at<0>(), DataType(-1));
      cur_sys.matrix_sys.local().block_d().apply(
        vec_def.local().template at<1>(), vec_primal.local().template at<0>(), vec_def.local().template at<1>(), DataType(-1));
      watch_grad_apply.stop();
    }

    FEAT::String _format_timings() const
    {
      enum _TimeID
      {
        total = 0,
        assembly = 1,
        fbm_filter = 2,
        grad_apply = 3,
        sync_and_rhs = 4,
        init = 5,
        create_vector = 6,
        num_entries = 7
      };
      FEAT::String s;
      double timings_max[_TimeID::num_entries], timings_min[_TimeID::num_entries], timings[_TimeID::num_entries];
      timings_min[_TimeID::total] = timings_max[_TimeID::total] = timings[_TimeID::total] = this->watch_total.elapsed();
      timings_min[_TimeID::assembly] = timings_max[_TimeID::assembly] = timings[_TimeID::assembly] = this->watch_assembly.elapsed();
      timings_min[_TimeID::fbm_filter] = timings_max[_TimeID::fbm_filter] = timings[_TimeID::fbm_filter] = this->watch_fbm_filter.elapsed();
      timings_min[_TimeID::grad_apply] = timings_max[_TimeID::grad_apply] = timings[_TimeID::grad_apply] = this->watch_grad_apply.elapsed();
      timings_min[_TimeID::sync_and_rhs] = timings_max[_TimeID::sync_and_rhs] = timings[_TimeID::sync_and_rhs] = this->watch_sync_and_rhs.elapsed();
      timings_min[_TimeID::init] = timings_max[_TimeID::init] = timings[_TimeID::init] = this->watch_init.elapsed();
      timings_min[_TimeID::create_vector] = timings_max[_TimeID::create_vector] = timings[_TimeID::create_vector] = this->watch_create_vector.elapsed();
      // sync our timings
      this->domain->comm().allreduce(timings_max, timings_max, _TimeID::num_entries, FEAT::Dist::op_max);
      this->domain->comm().allreduce(timings_min, timings_min, _TimeID::num_entries, FEAT::Dist::op_min);

      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Timings Single Material Defect Assembler----------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += format_subtime_mm("Sync and add RHS", timings[_TimeID::sync_and_rhs], timings[_TimeID::total], timings_min[_TimeID::sync_and_rhs], timings_max[_TimeID::sync_and_rhs], this->padlen);
      s += format_subtime_mm("Apply Gradient Matrix", timings[_TimeID::grad_apply], timings[_TimeID::total], timings_min[_TimeID::grad_apply], timings_max[_TimeID::grad_apply], this->padlen);
      s += format_subtime_mm("Apply FBM Filter", timings[_TimeID::fbm_filter], timings[_TimeID::total], timings_min[_TimeID::fbm_filter], timings_max[_TimeID::fbm_filter], this->padlen);
      s += format_subtime_mm("Assemble Defect", timings[_TimeID::assembly], timings[_TimeID::total], timings_min[_TimeID::assembly], timings_max[_TimeID::assembly], this->padlen);
      s += format_subtime_mm("Init Assembler", timings[_TimeID::init], timings[_TimeID::total], timings_min[_TimeID::init], timings_max[_TimeID::init], this->padlen);
      s += format_subtime_mm("Create Defect Vector", timings[_TimeID::create_vector], timings[_TimeID::total], timings_min[_TimeID::create_vector], timings_max[_TimeID::create_vector], this->padlen);
      s += format_subtime_mm("Defect Assembly Total Time", timings[_TimeID::total], timings[_TimeID::total], timings_min[_TimeID::total], timings_max[_TimeID::total], this->padlen);

      return s;
    }

    const Material<DataType>* _material;
    String _cubature = "gauss-legendre:5";
    DataType _temperature;
    DataType _mesh_scale;


  public:
    CarreauSingleMatDefectFlowAssembler() = default;

    CarreauSingleMatDefectFlowAssembler(DomainType& domain_, const std::shared_ptr<SystemLevel_>& system_level_, const FEAT::PropertyMap* prop, const Material<DataType>& material, DataType temperature, DataType mesh_scale)
      : BaseClass(domain_, system_level_, prop),
      _material(&material),
      _temperature(temperature),
      _mesh_scale(mesh_scale)
    {
      _parse(prop);
    }

    virtual ~CarreauSingleMatDefectFlowAssembler() = default;

    CarreauSingleMatDefectFlowAssembler(const CarreauSingleMatDefectFlowAssembler&) = default;
    CarreauSingleMatDefectFlowAssembler(CarreauSingleMatDefectFlowAssembler&&) noexcept = default;
    CarreauSingleMatDefectFlowAssembler& operator=(const CarreauSingleMatDefectFlowAssembler&) = default;
    CarreauSingleMatDefectFlowAssembler& operator=(CarreauSingleMatDefectFlowAssembler&&) noexcept = default;

    void set_material(const Material<DataType>& material, DataType mesh_scale)
    {
      _material = &material;
      _mesh_scale = mesh_scale;
    }

    void set_cubature(const String& cubature)
    {
      _cubature = cubature;
    }

  }; //class CarrauDefectFlowAssembler<...>
}
