// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#define FEAT_CUBATURE_TENSOR_PREFIX 1
#define FEAT_CUBATURE_SCALAR_PREFIX 1

#include <kernel/runtime.hpp>
#include <kernel/cubature/scalar/dynamic_factory.hpp>
#include <kernel/cubature/dynamic_factory.hpp>

#include <iostream>
#include <cstdio>
#include <set>
#include <map>

using namespace FEAT;
using namespace FEAT::Cubature;

typedef std::set<String> StringSet;
typedef std::map<String,String> StringMap;

template<
  typename Factory_,
  typename Functor_,
  bool variadic_ = (Factory_::variadic != 0)>
class AvailFunctorHelper;

template<
  typename Factory_,
  typename Functor_>
class AvailFunctorHelper<Factory_, Functor_, false>
{
private:
  Functor_& _functor;

public:
  explicit AvailFunctorHelper(Functor_& functor) :
    _functor(functor)
  {
    _functor.add_name(Factory_::name());
  }

  void alias(const String& name)
  {
    _functor.add_alias(name, Factory_::name());
  }
};

template<
  typename Factory_,
  typename Functor_>
class AvailFunctorHelper<Factory_, Functor_, true>
{
private:
  Functor_& _functor;

public:
  explicit AvailFunctorHelper(Functor_& functor) :
    _functor(functor)
  {
    _functor.add_name(Factory_::name() + ":<"
      + stringify(int(Factory_::min_points)) + "-"
      + stringify(int(Factory_::max_points)) + ">");
  }

  void alias(const String& name, int num_points)
  {
    _functor.add_alias(name, Factory_::name() + ":" + stringify(num_points));
  }
};

class AvailFunctor
{
private:
  StringSet& _names;
  StringMap& _aliases;

public:
  explicit AvailFunctor(StringSet& names, StringMap& aliases) :
    _names(names),
    _aliases(aliases)
  {
  }

  template<typename Factory_>
  void factory()
  {
    AvailFunctorHelper<Factory_, AvailFunctor> functor(*this);
    Factory_::alias(functor);
  }

  void add_name(const String& name)
  {
    _names.insert(name);
  }

  void add_alias(const String& alias, const String& name)
  {
    _aliases.insert(std::make_pair(alias, name));
  }
};

// Prints all available scalar cubature names
void print_avail_scalar()
{
  StringSet names;
  StringMap aliases;
  AvailFunctor functor(names, aliases);
  Scalar::FactoryWrapper::factory(functor);

  std::cout << "Available scalar cubature names:" << "\n";
  for(StringSet::iterator it(names.begin()); it != names.end(); ++it)
  {
    std::cout << *it << "\n";
  }

  std::cout << "\n" << "Available scalar cubature aliases:" << "\n";
  for(StringMap::iterator it(aliases.begin()); it != aliases.end(); ++it)
  {
    std::cout << it->first << " [ " << it->second << " ]" << "\n";
  }
  std::cout << "\n";
}

// Prints the cubature points for a specific scalar cubature rule
void test_dynamic_scalar(const String& name)
{
  // create rule
  Scalar::Rule<> rule;
  Scalar::DynamicFactory factory(name);
  if(!factory.create(rule))
  {
    std::cout << "ERROR: No scalar cubature rule named '" << name << "' found!" << "\n";
    return;
  }

  // print output
  std::cout << "\nPrinting scalar cubature rule '" << rule.get_name() << "':" << "\n";
  printf("      Weight          Coord\n");
  for(int i(0); i < rule.get_num_points(); ++i)
  {
    printf("%3i: %15.12f %15.12f\n", int(i), rule.get_weight(i), rule.get_coord(i));
  }
}

// Prints all available cubature names for a specific shape
template<typename Shape_>
void print_avail()
{
  StringSet names;
  StringMap aliases;
  AvailFunctor functor(names, aliases);
  FactoryWrapper<Shape_>::factory_no_refine(functor);

  std::cout << "Available cubature names for '" << Shape_::name() << "' :" << "\n";
  for(StringSet::iterator it(names.begin()); it != names.end(); ++it)
  {
    std::cout << *it << "\n";
  }

  std::cout << "\n" << "Available cubature aliases for '" << Shape_::name() << "' :" << "\n";
  for(StringMap::iterator it(aliases.begin()); it != aliases.end(); ++it)
  {
    std::cout << it->first << " [ " << it->second << " ]" << "\n";
  }
  std::cout << "\n";
}

// Prints the auto-degree aliases for a specific shape
template<typename Shape_>
void print_auto_degree()
{
  int max_degree = AutoAlias<Shape_>::max_auto_degree;

  std::cout << "\nPrinting auto-degree mappings for '" << Shape_::name() << "':" << "\n";
  for(int i(1); i <= max_degree; ++i)
  {
    std::cout << "auto-degree:" << i << " -> " << AutoAlias<Shape_>::map("auto-degree:" + stringify(i)) << "\n";
  }
}

// Prints the cubature points for a specific cubature rule and a specific shape
template<typename Shape_>
void test_dynamic(const String& name)
{
  if(name.compare_no_case("auto-degree") == 0)
  {
    print_auto_degree<Shape_>();
    return;
  }

  // create rule
  Rule<Shape_> rule;
  if(!DynamicFactory::create(rule, name))
  {
    std::cout << "ERROR: No " + Shape_::name() + " cubature rule named '" << name << "' found!" << "\n";
    return;
  }

  // print output
  std::cout << "\nPrinting " + Shape_::name() + " cubature rule '" << rule.get_name() << "':" << "\n";
  printf("      Weight          X-Coord");
  if(Shape_::dimension >= 2)
    printf("         Y-Coord");
  if(Shape_::dimension >= 3)
    printf("         Z-Coord");
  printf("\n");

  for(int i(0); i < rule.get_num_points(); ++i)
  {
    printf("%3i: %15.12f", int(i), rule.get_weight(i));
    for(int j(0); j < Shape_::dimension; ++j)
      printf(" %15.12f", rule.get_coord(i,j));
    printf("\n");
  }
}

// Prints all available cubature rules for all shapes
void print_avail_all()
{
  print_avail_scalar();
  print_avail< Shape::Simplex<1> >();
  print_avail< Shape::Simplex<2> >();
  print_avail< Shape::Simplex<3> >();
  print_avail< Shape::Hypercube<1> >();
  print_avail< Shape::Hypercube<2> >();
  print_avail< Shape::Hypercube<3> >();
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  if(argc < 2)
  {
    // print help message
    std::cout << "FEAT Dynamic Cubature List Tool" << "\n" << "\n";
    //            ---------1---------2---------3---------4---------5---------6---------7---------8
    //            123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    std::cout << "This tool can be used to display and debug all cubature rules currently" << "\n";
    std::cout << "supported by the FEAT::Cubature::DynamicFactory class." << "\n";
    std::cout << "\n";
    std::cout << "    USAGE:    cub-list <shape> [<name>]" << "\n";
    std::cout << "\n";
    std::cout << "where <shape> is either 's1', 's2', 's3', 'h1', 'h2' or 'h3' identifying the" << "\n";
    std::cout << "Simplex<n> or Hypercube<n> shape. Moreover, 'scalar' will refer to scalar rules." << "\n";
    std::cout << "\n";
    std::cout << "If both arguments are given, this tool will print out the weights and point" << "\n";
    std::cout << "coordinates of the cubature rule identified by <name>." << "\n";
    std::cout << "\n";
    std::cout << "If only the first argument is given, this tool will print a list of all" << "\n";
    std::cout << "available cubature rule names for the specified shape, including its aliases." << "\n";
    std::cout << "If <shape> is 'all', this tools prints out the cubature rules for all shapes." << "\n";
    std::cout << "\n" << "\n";
    std::cout << "Further information:" << "\n";
    std::cout << "--------------------" << "\n";
    std::cout << "There are two types of cubature rules: variadic and non-variadic ones." << "\n";
    std::cout << "A variaric cubature rule is parameterized in the number of desired cubature" << "\n";
    std::cout << "points. In this case the point-count parameter is appended by a colon to the" << "\n";
    std::cout << "cubature rule name; e.g. 'gauss-legendre:3', identifies the 3-point" << "\n";
    std::cout << "Gauss-Legendre rule." << "\n";
    std::cout << "For variadic cubature rules, this tool appends the minimum and maximum" << "\n";
    std::cout << "point count to the rule name, e.g. 'gauss-legendre:<1-5>' means that the" << "\n";
    std::cout << "DynamicFactory can generate Gauss-Legendre rules for at least 1 and at most" << "\n";
    std::cout << "5 points." << "\n";
    std::cout << "\n";
    std::cout << "For the Simplex<1> shape, cubature rules, which are generated from scalar rules," << "\n";
    std::cout << "are prefixed by 'scalar:'. For Hypercube<n> shapes, the prefix 'tensor:' states" << "\n";
    std::cout << "that the cubature rule is generated by a tensor-product of a scalar rule." << "\n";
    std::cout << "The 'scalar:' and 'tensor:' prefixes must not be used for the DynamicFactory" << "\n";
    std::cout << "by applications; these prefixes are only used by this tool for convenience." << "\n";
    std::cout << "\n";
    std::cout << "Several cubature rules have aliases, i.e. one and the same cubature rule may" << "\n";
    std::cout << "be identified by several names. This tool lists all available aliases and" << "\n";
    std::cout << "appends the 'real' name in square brackets. Example: The output line" << "\n";
    std::cout << "   simpson [ newton-cotes-closed:3 ]" << "\n";
    std::cout << "indicates that 'simpson' is an alias for the rule 'newton-cotes-closed:3'." << "\n";
    std::cout << "\n" << "\n";
    std::cout << "Automatic Rules:" << "\n";
    std::cout << "----------------" << "\n";
    std::cout << "In addition to the actual cubature rules and their aliases, there exist aliases" << "\n";
    std::cout << "called 'auto-aliased rules'. These auto-aliases are implemented for each shape" << "\n";
    std::cout << "and will try to select the most appropriate cubature rule fulfilling the desired" << "\n";
    std::cout << "requirements, if possible." << "\n";
    std::cout << "\n";
    std::cout << "Currently, the only supported auto-aliased rule is the 'auto-degree:<d>' alias," << "\n";
    std::cout << "which refers to the least-point inner cubature rule capable of integrating" << "\n";
    std::cout << "polynomials up to degree <d> exactly or to the highest-order cubature rule" << "\n";
    std::cout << "available if no exact cubature rule is available. Example: the 'auto-degree:2'" << "\n";
    std::cout << "alias refers to 'tensor:gauss-legendre:2' for Hypercube<n> shapes." << "\n";
    std::cout << "\n";
    std::cout << "A list of all auto-degree mappings supported by a shape can be printed by" << "\n";
    std::cout << "specifying 'auto-degree' as the cubature name, e.g. the command" << "\n";
    std::cout << "   cub-list h2 auto-degree" << "\n";
    std::cout << "will print all auto-degree aliases for quadrilaterals." << "\n";
    std::cout << "\n" << "\n";
    std::cout << "Refined Rules:" << "\n";
    std::cout << "--------------" << "\n";
    std::cout << "Each Simplex and Hypercube cubature rule may also be refined by prefixing the" << "\n";
    std::cout << "cubature name by 'refine*n:', where 'n' is a positive integer specifying the" << "\n";
    std::cout << "desired number of refinement steps. This allows to generate cubature rules of" << "\n";
    std::cout << "higher precision (but not higher order) from other cubature rules." << "\n";
    std::cout << "Note: 'refine:' is identical to 'refine*1:'." << "\n";
    return 0;
  }

  if(argc > 2)
  {
    if(String(argv[1]).compare_no_case("scalar") == 0)
      test_dynamic_scalar(argv[1]);
    else if(String(argv[1]).compare_no_case("s1") == 0)
      test_dynamic< Shape::Simplex<1> >(argv[2]);
    else if(String(argv[1]).compare_no_case("s2") == 0)
      test_dynamic< Shape::Simplex<2> >(argv[2]);
    else if(String(argv[1]).compare_no_case("s3") == 0)
      test_dynamic< Shape::Simplex<3> >(argv[2]);
    else if(String(argv[1]).compare_no_case("h1") == 0)
      test_dynamic< Shape::Hypercube<1> >(argv[2]);
    else if(String(argv[1]).compare_no_case("h2") == 0)
      test_dynamic< Shape::Hypercube<2> >(argv[2]);
    else if(String(argv[1]).compare_no_case("h3") == 0)
      test_dynamic< Shape::Hypercube<3> >(argv[2]);
    else
      std::cout << "ERROR: Unknown Shape '" << argv[1] << "'" << "\n";
  }
  else
  {
    if(String(argv[1]).compare_no_case("all") == 0)
      print_avail_all();
    else if(String(argv[1]).compare_no_case("scalar") == 0)
      print_avail_scalar();
    else if(String(argv[1]).compare_no_case("s1") == 0)
      print_avail< Shape::Simplex<1> >();
    else if(String(argv[1]).compare_no_case("s2") == 0)
      print_avail< Shape::Simplex<2> >();
    else if(String(argv[1]).compare_no_case("s3") == 0)
      print_avail< Shape::Simplex<3> >();
    else if(String(argv[1]).compare_no_case("h1") == 0)
      print_avail< Shape::Hypercube<1> >();
    else if(String(argv[1]).compare_no_case("h2") == 0)
      print_avail< Shape::Hypercube<2> >();
    else if(String(argv[1]).compare_no_case("h3") == 0)
      print_avail< Shape::Hypercube<3> >();
    else
      std::cout << "ERROR: Unknown Shape '" << argv[1] << "'" << "\n";
  }

  return 0;
}
