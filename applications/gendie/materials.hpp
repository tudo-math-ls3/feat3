// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include <kernel/runtime.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/voxel_assembly/voxel_assembly_common.hpp>
#include "parsing_helper.hpp"


#include <array>
#include <functional>

namespace Gendie
{
  using namespace FEAT;

  enum MaterialType
  {
    abc_model = 0,
    carreau = 1,
    carreau_yasuda = 2,
    power_law = 3,
    bingham = 4
  };

  enum TemperatureModelType
  {
    isotherm = 0,
    c1c2 = 1,
    tbts = 2,
    melt_tbts = 3,
    etb = 4
  };

  String stringify_enum(MaterialType ty)
  {
    switch(ty)
    {
      case Gendie::MaterialType::abc_model:
        return String("ABC Model");
      case Gendie::MaterialType::carreau:
        return String("Carreau Model");
      case Gendie::MaterialType::carreau_yasuda:
        return String("Carreau Yasuda Model");
      case Gendie::MaterialType::power_law:
        return String("Power Law Model");
      case Gendie::MaterialType::bingham:
        return String("Bingham Model");
    }
    return String("Unkown Model");
  }

  String stringify_enum(TemperatureModelType ty)
  {
    switch(ty)
    {
      case Gendie::TemperatureModelType::isotherm:
        return String("Isotherm");
      case Gendie::TemperatureModelType::c1c2:
        return String("C1C2");
      case Gendie::TemperatureModelType::tbts:
        return String("TBTS");
      case Gendie::TemperatureModelType::melt_tbts:
        return String("Melt TBTS");
      case Gendie::TemperatureModelType::etb:
        return String("ETB");
    }
    return String("Unkown Model");
  }

  std::ostream& operator<<(std::ostream& os, MaterialType mat)
  {
    os << stringify_enum(mat);
    return os;
  }

  std::ostream& operator<<(std::ostream& os, TemperatureModelType tmp_model)
  {
    os << stringify_enum(tmp_model);
    return os;
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
    static constexpr bool clamp_shearrate = false;

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

    static DataType get_visco_scaling_factor(DataType mesh_scale)
    {
      // mesh_scale -> 1000 for mm, 100 for cm, ...
      return DataType(1E3)/mesh_scale;
    }

    MaterialType get_type() const {return type;}
    std::array<DT_, data_len> get_data() const {return data;}
    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const
    {
      String s;
      s += String("Material Type").pad_back(padlen, pc) + ": " + stringify_enum(type);

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
                  Math::pow( DT2_(1.0) + Math::pow(b * gamma_dot, d), ((c - DT2_(1.0) - d) / d)) *
                  d * b * aT * Math::pow(b * gamma_dot * aT, d-DT2_(1.0));
        };
      }
      else
      {
        return [a=DT2_(this->data[0]), b=DT2_(this->data[1]), c=DT2_(this->data[2]), d=DT2_(this->data[3]), visco_scale=DT2_(visco_scaling_factor)](DT2_ gamma_dot, DT2_ aT) -> DT2_
        {
          return aT * a * visco_scale * ((c-DT2_(1.0)) / d) *
                  Math::pow( DT2_(1.0) + Math::pow(b * gamma_dot, d), ((c - DT2_(1.0) - d) / d)) *
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

    ViscosityDataCarreauModel(const PropertyMap* prop)
    {
      this->type = MaterialType::carreau;
      {
        DataType a,b,c;
        PARSE_PROP_OPTION(prop, "value_a", a, false);
        PARSE_PROP_OPTION(prop, "value_b", b, false);
        PARSE_PROP_OPTION(prop, "value_c", c, false);
        this->data[0] = a;
        this->data[1] = b;
        this->data[2] = -c+DataType(1);
        this->data[3] = DataType(1);
      }
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

    ViscosityDataCarreauYasudaModel() = default;
    ViscosityDataCarreauYasudaModel(const ViscosityDataCarreauYasudaModel&) = delete;
    ViscosityDataCarreauYasudaModel& operator=(const ViscosityDataCarreauYasudaModel&) = delete;
    ViscosityDataCarreauYasudaModel(ViscosityDataCarreauYasudaModel&&) = default;
    ViscosityDataCarreauYasudaModel& operator=(ViscosityDataCarreauYasudaModel&&) = default;
    virtual ~ViscosityDataCarreauYasudaModel(){}

    ViscosityDataCarreauYasudaModel(const PropertyMap* prop)
    {
      this->type = MaterialType::carreau_yasuda;
      {
        PARSE_PROP_OPTION(prop, "value_a", this->data[0], false);
        PARSE_PROP_OPTION(prop, "value_b", this->data[1], false);
        PARSE_PROP_OPTION(prop, "value_c", this->data[2], false);
        PARSE_PROP_OPTION(prop, "value_d", this->data[3], false);
      }
    }

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

    ViscosityDataPowerLawModel(const PropertyMap* prop)
    {
      this->type = MaterialType::power_law;
      {
        PARSE_PROP_OPTION(prop, "value_a", this->data[0], false);
        PARSE_PROP_OPTION(prop, "value_b", this->data[1], false);
      }
    }

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Mu_0").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[0], this->prec, this->width) + "\n";
      s += String("Powerlaw Exponent").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[1], this->prec, this->width);
      return s;
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

  };

  template<typename DT_>
  struct ViscosityDataABCModel : CarreauModelBase<DT_>
  {
    typedef DT_ DataType;
    typedef CarreauModelBase<DT_> BaseClass;
    // context
    // data[0] -> mu_0
    // data[1] -> reci
    // data[2] -> exponent

    template<typename... T>
    ViscosityDataABCModel(T... data) : BaseClass(data...){
      this->type = MaterialType::abc_model;
    }

    ViscosityDataABCModel() = default;
    ViscosityDataABCModel(const ViscosityDataABCModel&) = delete;
    ViscosityDataABCModel& operator=(const ViscosityDataABCModel&) = delete;
    ViscosityDataABCModel(ViscosityDataABCModel&&) = default;
    ViscosityDataABCModel& operator=(ViscosityDataABCModel&&) = default;
    virtual ~ViscosityDataABCModel(){}

    ViscosityDataABCModel(const PropertyMap* prop)
    {
      this->type = MaterialType::abc_model;
      {
        DataType a,b,c;
        PARSE_PROP_OPTION(prop, "value_a", a, false);
        PARSE_PROP_OPTION(prop, "value_b", b, false);
        PARSE_PROP_OPTION(prop, "value_c", c, false);
        this->data[0] = a;
        this->data[1] = b;
        this->data[2] = c+DataType(1);
        this->data[3] = DataType(1);
      }
    }

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
    // data[1] -> powerlaw exponent

    template<typename... T>
    ViscosityDataBinghamModel(T... data) : BaseClass(data...){
      this->type = MaterialType::power_law;
    }

    ViscosityDataBinghamModel(const PropertyMap* prop)
    {
      this->type = MaterialType::bingham;
      {
        PARSE_PROP_OPTION(prop, "value_a", this->data[0], false);
        PARSE_PROP_OPTION(prop, "value_b", this->data[1], false);
        PARSE_PROP_OPTION(prop, "value_c", this->data[2], false);
      }
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

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Mu_0").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[0], this->prec, this->width) + "\n";
      s += String("b").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[1], this->prec, this->width);
      s += String("Dynamic viscosity").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->data[2], this->prec, this->width);
      return s;
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

  };

  namespace Intern
  {
    bool parse_temperature_option(const std::pair<String, bool>& option, TemperatureModelType& in_type)
    {
      if(!option.second)
        return false;
      if(option.first.compare_no_case("isotherm") == 0)
      {
        in_type = TemperatureModelType::isotherm;
        return true;
      }
      else if(option.first.compare_no_case("tbts") == 0)
      {
        in_type = TemperatureModelType::tbts;
        return true;
      }
      XABORTM("No known temperature type " + option.first);
      return false;
    }

    bool parse_material_option(const std::pair<String, bool>& option, MaterialType& in_type)
    {
      if(!option.second)
        return false;
      if(option.first.compare_no_case("carreau") == 0)
      {
        in_type = MaterialType::carreau;
        return true;
      }
      else if (option.first.compare_no_case("carreau_yasuda") == 0)
      {
        in_type = MaterialType::carreau_yasuda;
        return true;
      }
      else if(option.first.compare_no_case("abc_model") == 0)
      {
        in_type = MaterialType::abc_model;
        return true;
      }
      else if(option.first.compare_no_case("powerlaw") == 0)
      {
        in_type = MaterialType::power_law;
        return true;
      }
      else if(option.first.compare_no_case("bingham") == 0)
      {
        in_type = MaterialType::bingham;
        return true;
      }
      XABORTM("No known model type: " + option.first);
      return false;
    }

    template<typename DT_>
    std::unique_ptr<ViscoModelBase<DT_>> create_viscosity_model(const FEAT::PropertyMap* prop)
    {
      typedef DT_ DataType;
      MaterialType model_type(MaterialType::carreau);
      std::unique_ptr<ViscoModelBase<DT_>> ptr;
      parse_material_option(prop->get_entry("model"), model_type);
      switch(model_type)
      {
        case Gendie::MaterialType::abc_model:
        {
          return std::make_unique<ViscosityDataABCModel<DataType>>(prop);
        }
        case Gendie::MaterialType::carreau:
        {
          return std::make_unique<ViscosityDataCarreauModel<DataType>>(prop);
        }
        case Gendie::MaterialType::carreau_yasuda:
        {
          return std::make_unique<ViscosityDataCarreauYasudaModel<DataType>>(prop);
        }
        case Gendie::MaterialType::power_law:
        {
          return std::make_unique<ViscosityDataPowerLawModel<DataType>>(prop);
        }
        case Gendie::MaterialType::bingham:
        {
          return std::make_unique<ViscosityDataBinghamModel<DataType>>(prop);
        }
        default:
        {
          XABORTM("Not implemented " + prop->get_entry("model").first + " yet");
        }
      }
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

    void _fill_temperature_model(const FEAT::PropertyMap* prop, TemperatureModelType& tmp, std::array<DataType, data_size_temp>& arr)
    {
      if(!Intern::parse_temperature_option(prop->get_entry("model"), tmp))
      {
        XABORTM("Could not parse the temperature model");
        return;
      }
      switch(tmp)
      {
        case Gendie::TemperatureModelType::isotherm:
        {
          {
            // TODO: Necessary to parse anything?
            // DataType ref_tmp(DataType(0));
            // PARSE_PROP_OPTION(prop, "reference_temperature", ref_tmp, false);
            // XASSERTM(ref_tmp > DataType(-273.15), "Temperature lower than absolute zero");
            arr[0]=DataType(0);
            arr[1]=DataType(0);
            arr[2]=DataType(0);
            arr[3]=DataType(0);
          }
          return;
        }
        case Gendie::TemperatureModelType::tbts:
        {
          // tmp[0] = TB,  tmp[1] = TS, tmp[2] = C1,  tmp[3] = C2
          {
            DataType TB(DataType(0));
            PARSE_PROP_OPTION(prop, "reference_temperature", TB, false);
            XASSERTM(TB > DataType(-273.15), "Constant parameter value_a gre");
            arr[0]=TB;
          }
          {
            DataType TS(DataType(0));
            PARSE_PROP_OPTION(prop, "value_a", TS, false);
            XASSERTM(TS > DataType(-273.15), "Constant parameter value_a gre");
            arr[1]=TS;
          }
          {
            DataType C1(DataType(8.86));
            arr[2]=C1;
          }
          {
            DataType C2(DataType(101.6));
            arr[3]=C2;
          }
          return;
        }
        case Gendie::TemperatureModelType::c1c2:
        {
          // tmp[0] = TB,  tmp[1] = TS, tmp[2] = C1,  tmp[3] = C2
          {
            DataType TB(DataType(0));
            PARSE_PROP_OPTION(prop, "reference_temperature", TB, false);
            XASSERTM(TB > DataType(-273.15), "Constant parameter value_a gre");
            arr[0]=TB;
          }
          {
            DataType TS(arr[1]);
            arr[1]=TS;
          }
          {
            DataType C1(DataType(8.86));
            arr[2]=C1;
          }
          {
            DataType C2(DataType(101.6));
            arr[3]=C2;
          }
          return;
        }
        default:
        {
          XABORTM("Not implemented " + prop->get_entry("model").first + " yet");
        }
      }
    }

    void _fill_viscosity_model(const FEAT::PropertyMap* prop)
    {
      this->_visc_mat_p = Intern::create_viscosity_model<DataType>(prop);
    }

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
        default:
        {
          auto lambda = [TB=_temp_mat_data[0], TS=_temp_mat_data[1], C1=_temp_mat_data[2], C2=_temp_mat_data[3]](DataType target_temp)
                        {
                          return Math::pow(DataType(10),
                                C1*(TB-TS)/(C2+TB-TS) - C1*(target_temp - TS)/(C2+target_temp-TS));
                        };
          return lambda(temp);
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


    /**
     * \brief ctor using a property map to parse the agruments
     *
     * \param[in] prop Property map with the following build
     *                 density: float
     *                 /temperature_model
     *                    model: String
     *                    ... temperature model parameters
     *                 /viscosity_model
     *                    model: String
     *                    ... viscosity model parameters
     */
    Material(const FEAT::PropertyMap* prop)
    {
      PARSE_PROP_OPTION(prop, "density", _density, false);
      _fill_temperature_model(prop->get_sub_section("temperature_model"), _temp_type, _temp_mat_data);
      _fill_viscosity_model(prop->get_sub_section("viscosity_model"));
    }

    /**
     * \brief Ctor writing visco mat ptr explicitly
     *
     */
    Material(DataType density, std::unique_ptr<ViscoModelBase<DataType>>&& visco_model) : _visc_mat_p(std::move(visco_model)), _density(density),  _temp_type(isotherm)
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
      s += String("Temperature Model").pad_back(padlen, pc) + ": " + stringify_enum(_temp_type) + "\n";
      s += String("TB").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temp_mat_data[0], 3, 10) + "\n";
      s += String("TS").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temp_mat_data[1], 3, 10) + "\n";
      s += String("C1").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temp_mat_data[2], 3, 10) + "\n";
      s += String("C2").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temp_mat_data[3], 3, 10) + "\n";
      s += _visc_mat_p->format_string(padlen, pc);

      return s;

    }
  }; //class Material<DT_>
}
