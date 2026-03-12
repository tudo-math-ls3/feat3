// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// This small application prints a human readable table with the shearrate
// to viscosity mappings of materials used in the gendie application
// Simply call ./material-viscosity-table --problem-file <problem_params.ini>

#include <kernel/runtime.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/assertion.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "materials.hpp"
#include "parsing_helper.hpp"
#include <kernel/util/simple_arg_parser.hpp>

using namespace FEAT;
using namespace Gendie;

// Function to print viscosity table for given material
template<typename DT_>
void print_viscosity_table(Material<DT_>& material, const std::vector<DT_>& shear_rates, const std::vector<DT_>& temperatures)
{
  String out;
  out += String("\nViscosity Table\n") + "=========================================\n";

  const int lw = 20;

  // Print header
  out += String("Shear Rate").pad_back(lw);
  for (const auto& temp : temperatures)
  {
    out += (String("Temp ") + stringify_fp_fix(temp, 1)).pad_back(lw);
  }
  out += "\n";

  // Print separator
  out += String("----------").pad_back(lw);
  for (size_t i = 0; i < temperatures.size(); ++i)
  {
    out += String("----------").pad_back(lw);
  }
  out += "\n";

  // Print viscosity values
  auto viscosity_func = material.get_visc_func(DT_(1E+3));

  for (const auto& gamma_dot : shear_rates)
  {
    out += stringify_fp_sci(gamma_dot, 2).pad_back(lw);
    for (const auto& temp : temperatures)
    {
      DT_ aT = material.get_abiat_aT(temp);
      DT_ viscosity = viscosity_func(gamma_dot, aT);
      out += stringify_fp_sci(viscosity, 4).pad_back(lw);
    }
    out += "\n";
  }
  out += "\n";

  std::cout << out;
}

// Function to print material info
template<typename DT_>
void print_material_info(Material<DT_>& material)
{
  std::cout << "Material Information\n";
  std::cout << "====================\n";
  std::cout << material.format_string() << "\n";
}

int main(int argc, char** argv)
{
  // Initialize FEAT runtime
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  SimpleArgParser args(argc, argv);

  args.support("problem-file", "The ini file containing materials under the section [problem_params/list_materials/k]");
  args.support("shearrate-range", "<lower> <higher> The range for the shearrates, from 1e-<lower> to 1e-<higher>, defaults to -2 6");
  args.support("shearrate-inc", "<increment> The multiplicative factor to increase the shearrate, Defaults to 10");
  args.support("temp-range", "<lower> <higher> The range for the temperatures in K, from <lower> <higher>, defaults from 200 400");
  args.support("temp-inc", "<increment> The additive factor to increase the temperature, Defaults to 40");

  args.query_unsupported();

  if(args.check("problem-file") <= 0)
  {
    XABORTM("You have to provide at least a problem file to be parsed");
  }

  try
  {
    // Create a property map and load from config file
    PropertyMap tmp_map;
    String filename = args.query("problem-file")->second.front();
    std::printf("Using property file %s \n", filename.c_str());
    tmp_map.read(filename);

    auto prop_map = tmp_map.treeify_structures();

    // Define shear rates and temperatures to test
    std::vector<double> shear_rates;
    {
      double start = 1E-2;
      double stop = 1E6;
      double increment = args.parse_default("shearrate-inc", double(10));

      if(args.check("shearrate-range") > 1)
      {
        double a = 1.;
        double b = 2.;
        auto opt = args.query("shearrate-range")->second;
        opt.at(0).parse(a);
        opt.at(1).parse(b);
        start = Math::pow(10, a);
        stop = Math::pow(10, b);
      }
      while(start <= stop)
      {
        shear_rates.emplace_back(start);
        start *= increment;
      }
    }
    std::vector<double> temperatures;
    {
      double start = 200;
      double stop = 400;
      double increment = args.parse_default("temp-inc", double(40));

      if(args.check("temp-range") > 1)
      {
        double a = 1.;
        double b = 2.;
        auto opt = args.query("temp-range")->second;
        opt.at(0).parse(a);
        opt.at(1).parse(b);
        start = a;
        stop = b;
      }
      while(start <= stop)
      {
        temperatures.emplace_back(start);
        start += increment;
      }
    }

    // Create the material using the property map
    std::vector<Material<double>> materials;
    parse_materials(materials, prop_map->get_sub_section("problem_params")->get_sub_section("list_materials"));

    for(auto& mat : materials)
    {
      // Print material information
      print_material_info(mat);

      // Print viscosity table
      print_viscosity_table(mat, shear_rates, temperatures);
    }

  }
  catch (const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}