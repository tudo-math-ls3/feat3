// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include <kernel/runtime.hpp>
#include <kernel/util/arg_parser.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/voxel_assembly/voxel_assembly_common.hpp>


#include <array>
#include <functional>

namespace Gendie
{
  using namespace FEAT;

  enum class MaterialType : std::uint8_t
  {
    abc_model = 0,
    carreau = 1,
    carreau_yasuda = 2,
    power_law = 3,
    bingham = 4,
    pseudo_power_law = 5
  };

  enum class TemperatureModelType : std::uint8_t
  {
    isotherm = 0,
    c1c2 = 1,
    tbts = 2,
    melt_tbts = 3,
    etb = 4
  };

  std::ostream& operator<<(std::ostream& os, MaterialType ty)
  {
    switch(ty)
    {
      case Gendie::MaterialType::abc_model:
        return os << "ABC Model";
      case Gendie::MaterialType::carreau:
        return os << "Carreau Model";
      case Gendie::MaterialType::carreau_yasuda:
        return os << "Carreau Yasuda Model";
      case Gendie::MaterialType::power_law:
        return os << "Power Law Model";
      case Gendie::MaterialType::pseudo_power_law:
        return  os << "Pseudo Power Law Model";
      case Gendie::MaterialType::bingham:
        return os << "Bingham Model";
    }
    return os;
  }

  std::ostream& operator<<(std::ostream& os, TemperatureModelType ty)
  {
    switch(ty)
    {
      case Gendie::TemperatureModelType::isotherm:
        return os << "Isotherm";
      case Gendie::TemperatureModelType::c1c2:
        return os << "C1C2";
      case Gendie::TemperatureModelType::tbts:
        return os << "TBTS";
      case Gendie::TemperatureModelType::melt_tbts:
        return os << "Melt TBTS";
      case Gendie::TemperatureModelType::etb:
        return os << "ETB";
    }
    return os << "Unkown Model";
  }

  namespace Intern
  {
    inline bool try_parse_temperature_option(const String& option, TemperatureModelType& in_type)
    {
      if(option.compare_no_case("isotherm") == 0)
      {
        in_type = TemperatureModelType::isotherm;
        return true;
      }
      else if(option.compare_no_case("tbts") == 0)
      {
        in_type = TemperatureModelType::tbts;
        return true;
      }
      else if(option.compare_no_case("etb") == 0)
      {
        in_type = TemperatureModelType::etb;
        return true;
      }
      else if(option.compare_no_case("c1c2") == 0)
      {
        in_type = TemperatureModelType::c1c2;
        return true;
      }
      return false;
    }

    static bool is_temperature_model(const String& option)
    {
      TemperatureModelType type = TemperatureModelType::isotherm;
      return try_parse_temperature_option(option, type);
    }

    inline bool try_parse_material(const String& material_string, MaterialType& mat_type)
    {
      if(material_string.compare_no_case("carreau") == 0)
      {
        mat_type = MaterialType::carreau;
        return true;
      }
      else if(material_string.compare_no_case("carreau_yasuda") == 0)
      {
        mat_type = MaterialType::carreau_yasuda;
        return true;
      }
      else if(material_string.compare_no_case("abc_model") == 0)
      {
        mat_type = MaterialType::abc_model;
        return true;
      }
      else if(material_string.compare_no_case("powerlaw") == 0)
      {
        mat_type = MaterialType::power_law;
        return true;
      }
      else if(material_string.compare_no_case("bingham") == 0)
      {
        mat_type = MaterialType::bingham;
        return true;
      }
      return false;
    }

    static inline bool is_viscosity_model(const String& option)
    {
      MaterialType type = MaterialType::carreau;
      return try_parse_material(option, type);
    }

    static inline bool is_tbts(const String& option)
    {
      return option.compare_no_case("tbts") == 0;
    }

    static inline bool is_tbts_or_c1c2(const String& option)
    {
      return option.compare_no_case("tbts") == 0 || option.compare_no_case("c1c2") == 0;
    }

    static inline bool needs_viscosity_value_c(const String& option)
    {
      return option.compare_no_case("carreau") == 0 || option.compare_no_case("carreau_yasuda") == 0 || option.compare_no_case("abc_model") == 0 || option.compare_no_case("bingham") == 0;
    }

    static inline bool needs_viscosity_value_d(const String& option)
    {
      return option.compare_no_case("carreau_yasuda") == 0;
    }

    template<typename DT_>
    struct MaterialSectionBaseParser : public ArgParser
    {
      using ArgParser::parameter;
      Parameter<DT_> density = parameter(DT_(0))
                                .property("density")
                                .required()
                                .name("Density");

      Parameter<DT_> temperature_reference_temperature = parameter(DT_(0))
                                                          .property("temperature_model/reference_temperature")
                                                          .name("Reference Temperature");
      Parameter<DT_> temperature_value_a = parameter(DT_(0))
                                            .property("temperature_model/value_a")
                                            .name("Temperature Value A");
      Parameter<String> temperature_model = parameter(String(""))
                                              .property("temperature_model/model")
                                              .required()
                                              .name("Temperature Model")
                                              .help_string("Valid options are: isotherm, tbts, c1c2 (later maybe: melt_tbts, etb)")
                                              .validator(is_temperature_model)
                                              .needs_if(temperature_value_a, is_tbts)
                                              .needs_if(temperature_reference_temperature, is_tbts_or_c1c2);

      Parameter<String> viscosity_model = parameter(String(""))
                                            .property("viscosity_model/model")
                                            .required()
                                            .name("Viscosity Model")
                                            .help_string("Valid options are: carreau, carreau_yasuda, abc_model, powerlaw, bingham")
                                            .validator(is_viscosity_model);

      TemperatureModelType temperature_type() const
      {
        TemperatureModelType tmp = TemperatureModelType::isotherm;
        const bool success = try_parse_temperature_option(*temperature_model, tmp);
        XASSERTM(success, "Could not parse temperature model");
        return tmp;
      }
    };

    template<typename DT_>
    struct MaterialSectionParser : public ArgParser
    {
      using ArgParser::parameter;
      Parameter<DT_> viscosity_value_a = parameter(DT_(0))
                                           .property("viscosity_model/value_a")
                                           .name("Viscosity Model Value A")
                                           .required();
      Parameter<DT_> viscosity_value_b = parameter(DT_(0))
                                           .property("viscosity_model/value_b")
                                           .name("Viscosity Model Value B")
                                           .required();
      Parameter<DT_> viscosity_value_c = parameter(DT_(0))
                                           .property("viscosity_model/value_c")
                                           .name("Viscosity Model Value C");
      Parameter<DT_> viscosity_value_d = parameter(DT_(0))
                                           .property("viscosity_model/value_d")
                                           .name("Viscosity Model Value D");
      Parameter<String> viscosity_model = parameter(String(""))
                                            .property("viscosity_model/model")
                                            .required()
                                            .name("Viscosity Model")
                                            .help_string("Valid options are: carreau (A B C), carreau_yasuda (A B C D), abc_model (A B C), powerlaw (A B), bingham (A B C)")
                                            .validator(is_viscosity_model)
                                            .needs_if(viscosity_value_c, needs_viscosity_value_c)
                                            .needs_if(viscosity_value_d, needs_viscosity_value_d);
    };

  }

  template<typename DT_>
  struct ViscoModelBase
  {
    typedef DT_ DataType;
    template<typename DT2_> using ViscoFuncT = std::function<DT2_(DT2_, DT2_)>;
    template<typename DT2_> using ViscoDerFuncT = std::function<DT2_(DT2_, DT2_)>;
    static constexpr int data_len = 4;
    static constexpr DT_ min_shear = 1E-2;
    static constexpr DT_ max_shear = 1E+5;
    std::array<DT_, data_len> data;
    MaterialType type;
    //formatting variables
    static constexpr int prec = 2;
    static constexpr int width = 5;
    static constexpr bool clamp_shearrate = true;

    template<typename... T>
    ViscoModelBase(T... in) : data({in...}) {}

    ViscoModelBase() = default;
    ViscoModelBase(const ViscoModelBase&) = delete;
    ViscoModelBase& operator=(const ViscoModelBase&) = delete;
    ViscoModelBase(ViscoModelBase&&) = default;
    ViscoModelBase& operator=(ViscoModelBase&&) = default;

    virtual ~ViscoModelBase(){}

    virtual ViscoFuncT<float> get_visc_func(float mesh_scale) const = 0;
    virtual ViscoDerFuncT<float> get_visc_der_func(float mesh_scale) const = 0;
    virtual ViscoFuncT<double> get_visc_func(double mesh_scale) const = 0;
    virtual ViscoDerFuncT<double> get_visc_der_func(double mesh_scale) const = 0;
  #ifdef FEAT_HAVE_QUADMATH
    virtual ViscoFuncT<__float128> get_visc_func(__float128 mesh_scale) const = 0;
    virtual ViscoDerFuncT<__float128> get_visc_der_func(__float128 mesh_scale) const = 0;
  #endif
  #ifdef FEAT_HAVE_HALFMATH
    virtual ViscoFuncT<half> get_visc_func(half mesh_scale) const = 0;
    virtual ViscoDerFuncT<half> get_visc_der_func(half mesh_scale) const = 0;
  #endif

    template<typename DT2_>
    static DT2_ get_visco_scaling_factor(DT2_ mesh_scale)
    {
      // mesh_scale -> 1000 for mm, 100 for cm, ...
      return DT2_(1E3)/mesh_scale;
    }

    MaterialType get_type() const {return type;}
    std::array<DT_, data_len> get_data() const {return data;}

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const
    {
      String s;
      s += String("Material Type").pad_back(padlen, pc) + ": " + stringify(type);

      return s;
    }
  };

  template<typename DT_>
  struct CarreauModelBase : public ViscoModelBase<DT_>
  {
    typedef DT_ DataType;
    typedef ViscoModelBase<DT_> BaseClass;
    template<typename T_> using ViscoFuncT = typename BaseClass::template ViscoFuncT<T_>;
    template<typename T_> using ViscoDerFuncT = typename BaseClass::template ViscoDerFuncT<T_>;

    // static constexpr bool clamp_shearrate = BaseClass::clamp_shearrate;
    using BaseClass::clamp_shearrate;

    template<typename DT2_>
    ViscoFuncT<DT2_> _get_visc_func(DT2_ mesh_scale) const
    {
      const DT2_ visco_scaling_factor = DT2_(this->get_visco_scaling_factor(mesh_scale));
      if constexpr(clamp_shearrate)
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor), min_shear=DT2_(this->min_shear), max_shear=DT2_(this->max_shear)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * Math::pow(DT2_(1) + Math::pow(b * aT * Math::clamp(gamma_dot, min_shear, max_shear), d), (c - DT2_(1))/d);
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * Math::pow(DT2_(1) + Math::pow(b * aT * gamma_dot, d), (c - DT2_(1))/d);
        };
      }
    }

    template<typename DT2_>
    ViscoDerFuncT<DT2_> _get_visc_der_func(DT2_ mesh_scale) const
    {
      const DT2_ visco_scaling_factor = DT2_(this->get_visco_scaling_factor(mesh_scale));
      if constexpr(clamp_shearrate)
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor), min_shear=DT2_(this->min_shear), max_shear=DT2_(this->max_shear)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          gamma_dot = Math::clamp(gamma_dot, min_shear, max_shear);
          return aT * a * visco_scale * ((c-DT2_(1.0)) / d) *
                  Math::pow( DT2_(1.0) + Math::pow(b * gamma_dot * aT, d), ((c - DT2_(1.0) - d) / d)) *
                  d * b * aT * Math::pow(b * gamma_dot * aT, d-DT2_(1.0));
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * ((c-DT2_(1.0)) / d) *
                  Math::pow( DT2_(1.0) + Math::pow(b * gamma_dot * aT, d), ((c - DT2_(1.0) - d) / d)) *
                  d * b * aT * Math::pow(b * gamma_dot * aT, d-DT2_(1.0));
        };
      }
    }
    // context
    // data[0] -> mu_0 in Pa.s i.e. kg/m/s
    // data[1] -> reciprocal velocity
    // data[2] -> exponent
    // data[3] -> plateau index
    template<typename... T>
    CarreauModelBase(T... data) : BaseClass(data...){
      this->type = MaterialType::carreau;
    }

    CarreauModelBase() = default;
    CarreauModelBase(const CarreauModelBase&) = delete;
    CarreauModelBase& operator=(const CarreauModelBase&) = delete;
    CarreauModelBase(CarreauModelBase&&) = default;
    CarreauModelBase& operator=(CarreauModelBase&&) = default;
    virtual ~CarreauModelBase(){}


    ViscoFuncT<float> get_visc_func(float mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<float> get_visc_der_func(float mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }

    ViscoFuncT<double> get_visc_func(double mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<double> get_visc_der_func(double mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }

  #ifdef FEAT_HAVE_QUADMATH
    ViscoFuncT<__float128> get_visc_func(__float128 mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<__float128> get_visc_der_func(__float128 mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }
  #endif
  #ifdef FEAT_HAVE_HALFMATH
    ViscoFuncT<half> get_visc_func(half mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<half> get_visc_der_func(half mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }
  #endif

  }; // class CarreauModelBase

  template<typename DT_>
  struct ViscosityDataCarreauModel : public CarreauModelBase<DT_>
  {
    typedef DT_ DataType;
    typedef CarreauModelBase<DT_> BaseClass;
    // context
    // data[0] -> mu_0   //in Pa.s i.e kg/m/s
    // data[1] -> reciprocal velocity
    // data[2] -> exponent

    template<typename... T>
    ViscosityDataCarreauModel(T... data) : BaseClass(data...){
      this->type = MaterialType::carreau;
    }

    ViscosityDataCarreauModel(const FEAT::PropertyMap& config)
     : BaseClass()
    {
      this->type = MaterialType::carreau;
      Intern::MaterialSectionParser<DataType> parser;
      if(!parser.parse(config))
      {
        for(const auto& err : parser.errors())
        {
          std::cerr << err << "\n";
        }
        XABORTM("Could not parse material");
      }
      this->data[0] = *parser.viscosity_value_a;
      this->data[1] = *parser.viscosity_value_b;
      this->data[2] = DataType(1) - *parser.viscosity_value_c;
      this->data[3] = DataType(1);
    }

    ViscosityDataCarreauModel() = default;
    ViscosityDataCarreauModel(const ViscosityDataCarreauModel&) = delete;
    ViscosityDataCarreauModel& operator=(const ViscosityDataCarreauModel&) = delete;
    ViscosityDataCarreauModel(ViscosityDataCarreauModel&&) = default;
    ViscosityDataCarreauModel& operator=(ViscosityDataCarreauModel&&) = default;
    virtual ~ViscosityDataCarreauModel(){}

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Mu_0").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[0], this->prec, this->width) + "\n";
      s += String("Recip Velocity lambda").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[1], this->prec, this->width) + "\n";
      s += String("Exponent n").pad_back(padlen, pc) + ": " + stringify_fp_fix(-this->data[2]+DataType(1), this->prec, this->width);
      return s;
    }


  };

  template<typename DT_>
  struct ViscosityDataCarreauYasudaModel : public CarreauModelBase<DT_>
  {
    typedef DT_ DataType;
    typedef CarreauModelBase<DT_> BaseClass;
    // context
    // data[0] -> mu_0
    // data[1] -> reciprocal velocity
    // data[2] -> exponent
    // data[3] -> plateau index

    template<typename... T>
    ViscosityDataCarreauYasudaModel(T... data) : BaseClass(data...){
      this->type = MaterialType::carreau_yasuda;
    }

    ViscosityDataCarreauYasudaModel(const FEAT::PropertyMap& config)
     : BaseClass()
    {
      this->type = MaterialType::carreau_yasuda;
      Intern::MaterialSectionParser<DataType> parser;
      if(!parser.parse(config))
      {
        for(const auto& err : parser.errors())
        {
          std::cerr << err << "\n";
        }
        XABORTM("Could not parse material");
      }
      this->data[0] = *parser.viscosity_value_a;
      this->data[1] = *parser.viscosity_value_b;
      this->data[2] = *parser.viscosity_value_c;
      this->data[3] = *parser.viscosity_value_d;
    }

    ViscosityDataCarreauYasudaModel() = default;
    ViscosityDataCarreauYasudaModel(const ViscosityDataCarreauYasudaModel&) = delete;
    ViscosityDataCarreauYasudaModel& operator=(const ViscosityDataCarreauYasudaModel&) = delete;
    ViscosityDataCarreauYasudaModel(ViscosityDataCarreauYasudaModel&&) = default;
    ViscosityDataCarreauYasudaModel& operator=(ViscosityDataCarreauYasudaModel&&) = default;
    virtual ~ViscosityDataCarreauYasudaModel(){}

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Mu_0").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[0], this->prec, this->width) + "\n";
      s += String("Recip Velocity lambda").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[1], this->prec, this->width) + "\n";
      s += String("Exponent n").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[2], this->prec, this->width) + "\n";
      s += String("Plateau Index").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[3], this->prec, this->width);
      return s;
    }
  };

  template<typename DT_>
  struct ViscosityDataPseudoPowerLawModel : public CarreauModelBase<DT_>
  {
    typedef DT_ DataType;
    typedef CarreauModelBase<DT_> BaseClass;

    template<typename T_>
    using ViscoFuncT = typename BaseClass::template ViscoFuncT<T_>;
    template<typename T_>
    using ViscoDerFuncT = typename BaseClass::template ViscoFuncT<T_>;
    /// clipping value used for shearrate
    DataType clipping_value = DataType(1E-3);
    // context
    // data[0] -> mu_0   //in Pa.s i.e kg/m/s
    // data[1] -> reciprocal velocity
    // data[2] -> exponent

    template<typename... T>
    ViscosityDataPseudoPowerLawModel(T... data) : BaseClass(data...){
      this->type = MaterialType::carreau;
    }

    ViscosityDataPseudoPowerLawModel(const PropertyMap& config)
    {
      this->type = MaterialType::pseudo_power_law;
      {
        Intern::MaterialSectionParser<DataType> parser;
        if(!parser.parse(config))
        {
          for(const auto& err : parser.errors())
          {
            std::cerr << err << "\n";
          }
          XABORTM("Could not parse material");
        }
        const DataType clipp_val_exp = FEAT::Math::pow(clipping_value, *parser.viscosity_value_b-DataType(1));
        this->data[0] = *parser.viscosity_value_a * clipp_val_exp;
        this->data[1] = DataType(1)/clipping_value;
        this->data[2] = *parser.viscosity_value_b;
        this->data[3] = DataType(1);
      }
    }

    ViscosityDataPseudoPowerLawModel() = default;
    ViscosityDataPseudoPowerLawModel(const ViscosityDataPseudoPowerLawModel&) = delete;
    ViscosityDataPseudoPowerLawModel& operator=(const ViscosityDataPseudoPowerLawModel&) = delete;
    ViscosityDataPseudoPowerLawModel(ViscosityDataPseudoPowerLawModel&&) = default;
    ViscosityDataPseudoPowerLawModel& operator=(ViscosityDataPseudoPowerLawModel&&) = default;
    virtual ~ViscosityDataPseudoPowerLawModel(){}

    using BaseClass::clamp_shearrate;

    template<typename DT2_>
    ViscoFuncT<DT2_> _get_visc_func(DT2_ mesh_scale) const
    {
      const DT2_ visco_scaling_factor = DT2_(this->get_visco_scaling_factor(mesh_scale));
      if constexpr(clamp_shearrate)
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor), min_shear=DT2_(this->min_shear), max_shear=DT2_(this->max_shear)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * Math::pow(DT2_(1) + Math::pow(b * Math::clamp(gamma_dot, min_shear, max_shear), d), (c - DT2_(1))/d);
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * Math::pow(DT2_(1) + Math::pow(b * gamma_dot, d), (c - DT2_(1))/d);
        };
      }
    }

    template<typename DT2_>
    ViscoDerFuncT<DT2_> _get_visc_der_func(DT2_ mesh_scale) const
    {
      const DT2_ visco_scaling_factor = DT2_(this->get_visco_scaling_factor(mesh_scale));
      if constexpr(clamp_shearrate)
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor), min_shear=DT2_(this->min_shear), max_shear=DT2_(this->max_shear)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          gamma_dot = Math::clamp(gamma_dot, min_shear, max_shear);
          return aT * a * visco_scale * ((c-DT2_(1.0)) / d) *
                  Math::pow( DT2_(1.0) + Math::pow(b * gamma_dot, d), ((c - DT2_(1.0) - d) / d)) *
                  d * b * Math::pow(b * gamma_dot, d-DT2_(1.0));
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * ((c-DT2_(1.0)) / d) *
                  Math::pow( DT2_(1.0) + Math::pow(b * gamma_dot, d), ((c - DT2_(1.0) - d) / d)) *
                  d * b * Math::pow(b * gamma_dot, d-DT2_(1.0));
        };
      }
    }
    ViscoFuncT<float> get_visc_func(float mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<float> get_visc_der_func(float mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }

    ViscoFuncT<double> get_visc_func(double mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<double> get_visc_der_func(double mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }

  #ifdef FEAT_HAVE_QUADMATH
    ViscoFuncT<__float128> get_visc_func(__float128 mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<__float128> get_visc_der_func(__float128 mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }
  #endif
  #ifdef FEAT_HAVE_HALFMATH
    ViscoFuncT<half> get_visc_func(half mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<half> get_visc_der_func(half mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }
  #endif

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Mu_0").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[0], this->prec, this->width) + "\n";
      s += String("Recip Velocity lambda").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[1], this->prec, this->width) + "\n";
      s += String("Exponent n").pad_back(padlen, pc) + ": " + stringify_fp_fix(-this->data[2]+DataType(1), this->prec, this->width);
      return s;
    }
  };

  template<typename DT_>
  struct ViscosityDataPowerLawModel: public ViscoModelBase<DT_>
  {
    typedef DT_ DataType;
    typedef ViscoModelBase<DT_> BaseClass;
    template<typename T_> using ViscoFuncT = typename BaseClass::template ViscoFuncT<T_>;
    template<typename T_> using ViscoDerFuncT = typename BaseClass::template ViscoDerFuncT<T_>;
    // context
    // data[0] -> mu_0
    // data[1] -> powerlaw exponent

    static constexpr bool p_clamp_shearrate = true;

    template<typename... T>
    ViscosityDataPowerLawModel(T... data) : BaseClass(data...){
      this->type = MaterialType::power_law;
    }

    ViscosityDataPowerLawModel(const FEAT::PropertyMap& config)
     : BaseClass()
    {
      this->type = MaterialType::power_law;
      Intern::MaterialSectionParser<DataType> parser;
      if(!parser.parse(config))
      {
        for(const auto& err : parser.errors())
        {
          std::cerr << err << "\n";
        }
        XABORTM("Could not parse material");
      }
      this->data[0] = *parser.viscosity_value_a;
      this->data[1] = *parser.viscosity_value_b;
      this->data[2] = DataType(0);
      this->data[3] = DataType(0);
    }

    ViscosityDataPowerLawModel() = default;
    ViscosityDataPowerLawModel(const ViscosityDataPowerLawModel&) = delete;
    ViscosityDataPowerLawModel& operator=(const ViscosityDataPowerLawModel&) = delete;
    ViscosityDataPowerLawModel(ViscosityDataPowerLawModel&&) = default;
    ViscosityDataPowerLawModel& operator=(ViscosityDataPowerLawModel&&) = default;
    virtual ~ViscosityDataPowerLawModel(){}

    template<typename DT2_>
    ViscoFuncT<DT2_> _get_visc_func(DT2_ mesh_scale) const
    {
      const DT2_ visco_scaling_factor = DT2_(this->get_visco_scaling_factor(mesh_scale));
      if constexpr(p_clamp_shearrate)
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), visco_scale=DT2_(visco_scaling_factor), min_shear=DT2_(this->min_shear), max_shear=DT2_(this->max_shear)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * Math::pow(Math::clamp(gamma_dot, min_shear, max_shear), -(DT2_(1)-b));
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * Math::pow(gamma_dot, -(DT2_(1)-b));
        };
      }
    }

    template<typename DT2_>
    ViscoDerFuncT<DT2_> _get_visc_der_func(DT2_ mesh_scale) const
    {
      const DT2_ visco_scaling_factor = DT2_(this->get_visco_scaling_factor(mesh_scale));
      if constexpr(p_clamp_shearrate)
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), visco_scale=DT2_(visco_scaling_factor), min_shear=DT2_(this->min_shear), max_shear=DT2_(this->max_shear)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * (b-DT2_(1)) * Math::pow(Math::clamp(gamma_dot, min_shear, max_shear), -(DT2_(2)-b));
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * (b-DT2_(1)) * Math::pow(gamma_dot, -(DT2_(2)-b));
        };
      }
    }

    ViscoFuncT<float> get_visc_func(float mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<float> get_visc_der_func(float mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }

    ViscoFuncT<double> get_visc_func(double mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<double> get_visc_der_func(double mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }

  #ifdef FEAT_HAVE_QUADMATH
    ViscoFuncT<__float128> get_visc_func(__float128 mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<__float128> get_visc_der_func(__float128 mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }
  #endif
  #ifdef FEAT_HAVE_HALFMATH
    ViscoFuncT<half> get_visc_func(half mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<half> get_visc_der_func(half mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }
  #endif

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Mu_0").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[0], this->prec, this->width) + "\n";
      s += String("Powerlaw Exponent").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[1], this->prec, this->width);
      return s;
    }
  };

  template<typename DT_>
  struct ViscosityDataABCModel : CarreauModelBase<DT_>
  {
    typedef DT_ DataType;
    typedef CarreauModelBase<DT_> BaseClass;
    // context
    // data[0] -> A
    // data[1] -> B
    // data[2] -> C

    template<typename... T>
    ViscosityDataABCModel(T... data) : BaseClass(data...){
      this->type = MaterialType::abc_model;
    }

    ViscosityDataABCModel(const FEAT::PropertyMap& config)
     : BaseClass()
    {
      this->type = MaterialType::abc_model;
      Intern::MaterialSectionParser<DataType> parser;
      if(!parser.parse(config))
      {
        for(const auto& err : parser.errors())
        {
          std::cerr << err << "\n";
        }
        XABORTM("Could not parse material");
      }
      this->data[0] = *parser.viscosity_value_a;
      this->data[1] = *parser.viscosity_value_b;
      this->data[2] = *parser.viscosity_value_c + DataType(1);
      this->data[3] = DataType(1);
    }

    ViscosityDataABCModel() = default;
    ViscosityDataABCModel(const ViscosityDataABCModel&) = delete;
    ViscosityDataABCModel& operator=(const ViscosityDataABCModel&) = delete;
    ViscosityDataABCModel(ViscosityDataABCModel&&) = default;
    ViscosityDataABCModel& operator=(ViscosityDataABCModel&&) = default;
    virtual ~ViscosityDataABCModel(){}

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("A").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[0], this->prec, this->width) + "\n";
      s += String("B").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[1], this->prec, this->width) + "\n";
      s += String("C").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[2]-DataType(1), this->prec, this->width);
      return s;
    }
  };

  template<typename DT_>
  struct ViscosityDataBinghamModel : public ViscoModelBase<DT_>
  {
    typedef DT_ DataType;
    typedef ViscoModelBase<DT_> BaseClass;
    template<typename T_> using ViscoFuncT = typename BaseClass::template ViscoFuncT<T_>;
    template<typename T_> using ViscoDerFuncT = typename BaseClass::template ViscoDerFuncT<T_>;
    using BaseClass::clamp_shearrate;
    // context
    // data[0] -> mu_0
    // data[1] -> b
    // data[2] -> dynamic viscosity

    template<typename... T>
    ViscosityDataBinghamModel(T... data) : BaseClass(data...){
      this->type = MaterialType::power_law;
    }

    ViscosityDataBinghamModel(const FEAT::PropertyMap& config)
     : BaseClass()
    {
      this->type = MaterialType::bingham;
      Intern::MaterialSectionParser<DataType> parser;
      if(!parser.parse(config))
      {
        for(const auto& err : parser.errors())
        {
          std::cerr << err << "\n";
        }
        XABORTM("Could not parse material");
      }
      this->data[0] = *parser.viscosity_value_a;
      this->data[1] = *parser.viscosity_value_b;
      this->data[2] = *parser.viscosity_value_c;
      this->data[3] = DataType(0);
    }

    ViscosityDataBinghamModel() = default;
    ViscosityDataBinghamModel(const ViscosityDataBinghamModel&) = delete;
    ViscosityDataBinghamModel& operator=(const ViscosityDataBinghamModel&) = delete;
    ViscosityDataBinghamModel(ViscosityDataBinghamModel&&) = default;
    ViscosityDataBinghamModel& operator=(ViscosityDataBinghamModel&&) = default;
    virtual ~ViscosityDataBinghamModel(){}

    template<typename DT2_>
    ViscoFuncT<DT2_> _get_visc_func(DT2_ mesh_scale) const
    {
      const DT2_ visco_scaling_factor = DT2_(this->get_visco_scaling_factor(mesh_scale));
      if constexpr(clamp_shearrate)
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), visco_scale=DT2_(visco_scaling_factor), min_shear=DT2_(this->min_shear), max_shear=DT2_(this->max_shear)](DT2_ gamma_dot, DT2_) -> DT2_
        {
          return a * visco_scale + c * visco_scale / (b + Math::clamp(gamma_dot, min_shear, max_shear));
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_) -> DT2_
        {
          return a * visco_scale + c * visco_scale / (b + gamma_dot);
        };
      }
    }

    template<typename DT2_>
    ViscoDerFuncT<DT2_> _get_visc_der_func(DT2_ mesh_scale) const
    {
      const DT2_ visco_scaling_factor = DT2_(this->get_visco_scaling_factor(mesh_scale));
      if constexpr(clamp_shearrate)
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), visco_scale=DT2_(visco_scaling_factor), min_shear=DT2_(this->min_shear), max_shear=DT2_(this->max_shear)](DT2_ gamma_dot, DT2_) -> DT2_
        {
          return -c * visco_scale / Math::sqr(b + Math::clamp(gamma_dot, min_shear, max_shear));
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_) -> DT2_
        {
          return -c * visco_scale / Math::sqr(b + gamma_dot);
        };
      }
    }

    ViscoFuncT<float> get_visc_func(float mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<float> get_visc_der_func(float mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }

    ViscoFuncT<double> get_visc_func(double mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<double> get_visc_der_func(double mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }

  #ifdef FEAT_HAVE_QUADMATH
    ViscoFuncT<__float128> get_visc_func(__float128 mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<__float128> get_visc_der_func(__float128 mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }
  #endif
  #ifdef FEAT_HAVE_HALFMATH
    ViscoFuncT<half> get_visc_func(half mesh_scale) const override
    {
      return _get_visc_func(mesh_scale);
    }

    ViscoDerFuncT<half> get_visc_der_func(half mesh_scale) const override
    {
      return _get_visc_der_func(mesh_scale);
    }
  #endif

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Mu_0").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[0], this->prec, this->width) + "\n";
      s += String("b").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[1], this->prec, this->width);
      s += String("Dynamic viscosity").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[2], this->prec, this->width);
      return s;
    }

  };

  namespace Intern
  {
    template<typename DT_>
    std::unique_ptr<ViscoModelBase<DT_>> create_viscosity_model(const FEAT::PropertyMap& config)
    {
      typedef DT_ DataType;
      MaterialSectionBaseParser<DT_> base_parser;
      if(!base_parser.parse(config))
      {
        for(const auto& err : base_parser.errors())
        {
          std::cerr << err << "\n";
        }
        XABORTM("Parsing of Visco Mat went wrong");
      }
      Gendie::MaterialType model_type = Gendie::MaterialType::carreau;
      if(!try_parse_material(*base_parser.viscosity_model, model_type))
      {
        XABORTM(String("Unkown MaterialType ") + *base_parser.viscosity_model);
      }
      switch(model_type)
      {
        case Gendie::MaterialType::abc_model:
          return std::make_unique<ViscosityDataABCModel<DataType>>(config);
        case Gendie::MaterialType::carreau:
          return std::make_unique<ViscosityDataCarreauModel<DataType>>(config);
        case Gendie::MaterialType::carreau_yasuda:
          return std::make_unique<ViscosityDataCarreauYasudaModel<DataType>>(config);
        case Gendie::MaterialType::power_law:
          // return std::make_unique<ViscosityDataPowerLawModel<DataType>>(config); // reroute to pseudo powerlaw
        case Gendie::MaterialType::pseudo_power_law:
          return std::make_unique<ViscosityDataPseudoPowerLawModel<DataType>>(config);
        case Gendie::MaterialType::bingham:
          return std::make_unique<ViscosityDataBinghamModel<DataType>>(config);
      }
      XABORTM("Unsupported viscosity model");
      return nullptr;
    }
  }


  template<typename DT_>
  class Material
  {
  public:
    typedef DT_ DataType;
    static constexpr int data_size_temp = 4;
  private:
    std::unique_ptr<ViscoModelBase<DT_>> _visc_mat_p;
    std::array<DataType, data_size_temp> _temp_mat_data;
    /// density in kg/m^3
    DataType _density;
    TemperatureModelType _temp_type;

  protected:
    static constexpr std::size_t padlen = std::size_t(30);

    DataType _get_pseudo_at(DataType temp) const
    {
      switch (_temp_type)
      {
        case Gendie::TemperatureModelType::isotherm:
        {
          return DataType(1);
        }
        case Gendie::TemperatureModelType::c1c2:
        {
          auto lambda = [TB=_temp_mat_data[0], C1=_temp_mat_data[2], C2=_temp_mat_data[3]](DataType target_temp)
                        {
                          return Math::exp(- C1*(target_temp - TB)/(C2+target_temp-TB));
                        };
          return lambda(temp);

        }
        case Gendie::TemperatureModelType::tbts:
        {
          auto lambda = [TB=_temp_mat_data[0], TS=_temp_mat_data[1], C1=_temp_mat_data[2], C2=_temp_mat_data[3]](DataType target_temp)
                        {
                          return Math::pow(DataType(10),
                                C1*(TB-TS)/(C2+TB-TS) - C1*(target_temp - TS)/(C2+target_temp-TS));
                        };
          return lambda(temp);
        }
        case Gendie::TemperatureModelType::etb:
        {
          // ETB model: aT = exp(aux), aux = (E/R) * (1/(T+273.15) - 1/(T_B+273.15))
          // arr[0] = TB, arr[1] = E
          auto lambda = [TB=_temp_mat_data[0], E=_temp_mat_data[1]](DataType target_temp)
                        {
                          const DataType R = DataType(8.314462618); // Universal gas constant in J/(mol·K)
                          const DataType T_kelvin = target_temp + DataType(273.15);
                          const DataType TB_kelvin = TB + DataType(273.15);
                          const DataType aux = (E / R) * (DataType(1.0) / T_kelvin - DataType(1.0) / TB_kelvin);
                          return Math::exp(aux);
                        };
          return lambda(temp);
        }
        default:
        {
          XABORTM("Not implemented temperature model for _get_pseudo_at");
          return DataType(0);
        }
      }

    }

  public:
    Material() = default;
    Material(const Material&) = delete;
    Material& operator=(const Material&) = delete;
    Material(Material&&) = default;
    Material& operator=(Material&&) = default;
    ~Material() = default;

    Material(const PropertyMap& config)
    {
      Intern::MaterialSectionBaseParser<DataType> parser;
      if(!parser.parse(config))
      {
        for(const auto& err : parser.errors())
        {
          std::cerr << err << "\n";
        }
        XABORTM("Could not parse base material");
      }
      _density = *parser.density;
      _temp_type = parser.temperature_type();
      switch(_temp_type)
      {
        case Gendie::TemperatureModelType::isotherm:
          _temp_mat_data = {DataType(0), DataType(0), DataType(0), DataType(0)};
          break;
        case Gendie::TemperatureModelType::tbts:
          _temp_mat_data = {*parser.temperature_reference_temperature, *parser.temperature_value_a, DataType(8.86), DataType(101.6)};
          break;
        case Gendie::TemperatureModelType::c1c2:
          _temp_mat_data = {*parser.temperature_reference_temperature, DataType(0), DataType(8.86), DataType(101.6)};
          break;
        case Gendie::TemperatureModelType::etb:
          _temp_mat_data = {*parser.temperature_reference_temperature, *parser.temperature_value_a, DataType(0), DataType(0)};
          break;
        default:
          XABORTM("Unsupported temperature model");
      }
      _visc_mat_p = Intern::create_viscosity_model<DataType>(config);
    }

    /**
     * \brief Ctor writing visco mat ptr explicitly
     *
     */
    Material(DataType density, std::unique_ptr<ViscoModelBase<DataType>>&& visco_model) : _visc_mat_p(std::move(visco_model)), _density(density),  _temp_type(TemperatureModelType::isotherm)
    {
      _temp_mat_data[0] = DataType(200.);
    }

    template<typename DT2_>
    auto get_visc_func(DT2_ mesh_scaling_factor) const
    {
      return _visc_mat_p->get_visc_func(mesh_scaling_factor);
    }

    template<typename DT2_>
    auto get_visc_der_func(DT2_ mesh_scaling_factor) const
    {
      return _visc_mat_p->get_visc_der_func(mesh_scaling_factor);
    }

    void set_viscosity_model(DataType density, std::unique_ptr<ViscoModelBase<DataType>>&& visco_model)
    {
      _density = density;
      _visc_mat_p = std::move(visco_model);
    }


    const ViscoModelBase<DataType>* get_viscosity_model() const
    {
      return _visc_mat_p.get();
    }

    MaterialType get_current_visc_model_type() const
    {
      return _visc_mat_p->get_type();
    }

    DataType get_density() const
    {
      return _density;
    }

    DataType get_density_gram_per_unit(const DataType mesh_scaling_factor = 1E3) const
    {
      const DataType trans_val = DataType(1E3)/Math::cub(mesh_scaling_factor);
      return _density*trans_val;
    }

    DataType get_abiat_aT(DataType temperature) const
    {
      return this->_get_pseudo_at(temperature);
    }

    String format_string(int material_id = 0) const
    {
      String s;
      const char pc = '.';

      s = "Material " + stringify(material_id)  + " parameters:\n";
      s += String("Density ").pad_back(padlen, pc) + ": " + stringify_fp_fix(_density, 3, 10) + " kg/m^3\n";
      s += String("Temperature Model").pad_back(padlen, pc) + ": " + stringify(_temp_type) + "\n";
      s += String("TB").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temp_mat_data[0], 3, 10) + "\n";
      s += String("TS").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temp_mat_data[1], 3, 10) + "\n";
      s += String("C1").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temp_mat_data[2], 3, 10) + "\n";
      s += String("C2").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temp_mat_data[3], 3, 10) + "\n";
      s += _visc_mat_p->format_string(padlen, pc);

      return s;
    }

  }; //class Material<DT_>

  enum class MaterialParserError : std::uint8_t
  {
    success = 0,
    no_material_config = 1,
    no_materials = 2,
    invalid_material = 3
  };

  static inline const char* parse_error(MaterialParserError t)
  {
    switch(t)
    {
      case MaterialParserError::success:
        return "Success";
      case MaterialParserError::no_material_config:
        return "No Material Config Found";
      case MaterialParserError::no_materials:
        return "No Material Entry found in config";
      case MaterialParserError::invalid_material:
        return "Failed to parse material entry";
    }
    return "Unknown";
  }

  template<typename DT_>
  MaterialParserError parse_materials(std::vector<Material<DT_>>& materials, const FEAT::PropertyMap* mat_config, String& material_info)
  {
    XASSERTM(mat_config, "No material config list found");
    if(!mat_config)
      return MaterialParserError::no_material_config;
    std::size_t i = 0u;
    while(true)
    {
      const auto* cur_mat = mat_config->query_section(stringify(i));
      if(cur_mat == nullptr)
      {
        XASSERTM(i > 0u, "Did not querry any materials!");
        return i > 0u ? MaterialParserError::success : MaterialParserError::no_materials;
      }
      materials.emplace_back(*cur_mat);
      material_info += materials.back().format_string(int(i));
      ++i;
    }

    return MaterialParserError::success;
  }
}
