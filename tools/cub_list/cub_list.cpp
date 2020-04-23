// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#define FEAT_CUBATURE_TENSOR_PREFIX 1
#define FEAT_CUBATURE_SCALAR_PREFIX 1
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

  std::cout << "Available scalar cubature names:" << std::endl;
  for(StringSet::iterator it(names.begin()); it != names.end(); ++it)
  {
    std::cout << *it << std::endl;
  }

  std::cout << std::endl << "Available scalar cubature aliases:" << std::endl;
  for(StringMap::iterator it(aliases.begin()); it != aliases.end(); ++it)
  {
    std::cout << it->first << " [ " << it->second << " ]" << std::endl;
  }
  std::cout << std::endl;
}

// Prints the cubature points for a specific scalar cubature rule
void test_dynamic_scalar(const String& name)
{
  // create rule
  Scalar::Rule<> rule;
  Scalar::DynamicFactory factory(name);
  if(!factory.create(rule))
  {
    std::cout << "ERROR: No scalar cubature rule named '" << name << "' found!" << std::endl;
    return;
  }

  // print output
  std::cout << "\nPrinting scalar cubature rule '" << rule.get_name() << "':" << std::endl;
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

  std::cout << "Available cubature names for '" << Shape_::name() << "' :" << std::endl;
  for(StringSet::iterator it(names.begin()); it != names.end(); ++it)
  {
    std::cout << *it << std::endl;
  }

  std::cout << std::endl << "Available cubature aliases for '" << Shape_::name() << "' :" << std::endl;
  for(StringMap::iterator it(aliases.begin()); it != aliases.end(); ++it)
  {
    std::cout << it->first << " [ " << it->second << " ]" << std::endl;
  }
  std::cout << std::endl;
}

// Prints the auto-degree aliases for a specific shape
template<typename Shape_>
void print_auto_degree()
{
  int max_degree = AutoAlias<Shape_>::max_auto_degree;

  std::cout << "\nPrinting auto-degree mappings for '" << Shape_::name() << "':" << std::endl;
  for(int i(1); i <= max_degree; ++i)
  {
    std::cout << "auto-degree:" << i << " -> " << AutoAlias<Shape_>::map("auto-degree:" + stringify(i)) << std::endl;
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
    std::cout << "ERROR: No " + Shape_::name() + " cubature rule named '" << name << "' found!" << std::endl;
    return;
  }

  // print output
  std::cout << "\nPrinting " + Shape_::name() + " cubature rule '" << rule.get_name() << "':" << std::endl;
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
  if(argc < 2)
  {
    // print help message
    std::cout << "FEAT Dynamic Cubature List Tool" << std::endl << std::endl;
    //            ---------1---------2---------3---------4---------5---------6---------7---------8
    //            123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    std::cout << "This tool can be used to display and debug all cubature rules currently" << std::endl;
    std::cout << "supported by the FEAT::Cubature::DynamicFactory class." << std::endl;
    std::cout << std::endl;
    std::cout << "    USAGE:    cub-list <shape> [<name>]" << std::endl;
    std::cout << std::endl;
    std::cout << "where <shape> is either 's1', 's2', 's3', 'h1', 'h2' or 'h3' identifying the" << std::endl;
    std::cout << "Simplex<n> or Hypercube<n> shape. Moreover, 'scalar' will refer to scalar rules." << std::endl;
    std::cout << std::endl;
    std::cout << "If both arguments are given, this tool will print out the weights and point" << std::endl;
    std::cout << "coordinates of the cubature rule identified by <name>." << std::endl;
    std::cout << std::endl;
    std::cout << "If only the first argument is given, this tool will print a list of all" << std::endl;
    std::cout << "available cubature rule names for the specified shape, including its aliases." << std::endl;
    std::cout << "If <shape> is 'all', this tools prints out the cubature rules for all shapes." << std::endl;
    std::cout << std::endl << std::endl;
    std::cout << "Further information:" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "There are two types of cubature rules: variadic and non-variadic ones." << std::endl;
    std::cout << "A variaric cubature rule is parameterised in the number of desired cubature" << std::endl;
    std::cout << "points. In this case the point-count parameter is appended by a colon to the" << std::endl;
    std::cout << "cubature rule name; e.g. 'gauss-legendre:3', identifies the 3-point" << std::endl;
    std::cout << "Gauss-Legendre rule." << std::endl;
    std::cout << "For variadic cubature rules, this tool appends the minimum and maximum" << std::endl;
    std::cout << "point count to the rule name, e.g. 'gauss-legendre:<1-5>' means that the" << std::endl;
    std::cout << "DynamicFactory can generate Gauss-Legendre rules for at least 1 and at most" << std::endl;
    std::cout << "5 points." << std::endl;
    std::cout << std::endl;
    std::cout << "For the Simplex<1> shape, cubature rules, which are generated from scalar rules," << std::endl;
    std::cout << "are prefixed by 'scalar:'. For Hypercube<n> shapes, the prefix 'tensor:' states" << std::endl;
    std::cout << "that the cubature rule is generated by a tensor-product of a scalar rule." << std::endl;
    std::cout << "The 'scalar:' and 'tensor:' prefixes must not be used for the DynamicFactory" << std::endl;
    std::cout << "by applications; these prefixes are only used by this tool for convenience." << std::endl;
    std::cout << std::endl;
    std::cout << "Several cubature rules have aliases, i.e. one and the same cubature rule may" << std::endl;
    std::cout << "be identified by several names. This tool lists all available aliases and" << std::endl;
    std::cout << "appends the 'real' name in square brackets. Example: The output line" << std::endl;
    std::cout << "   simpson [ newton-cotes-closed:3 ]" << std::endl;
    std::cout << "indicates that 'simpson' is an alias for the rule 'newton-cotes-closed:3'." << std::endl;
    std::cout << std::endl << std::endl;
    std::cout << "Automatic Rules:" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout << "In addition to the actual cubature rules and their aliases, there exist aliases" << std::endl;
    std::cout << "called 'auto-aliased rules'. These auto-aliases are implemented for each shape" << std::endl;
    std::cout << "and will try to select the most appropriate cubature rule fulfilling the desired" << std::endl;
    std::cout << "requirements, if possible." << std::endl;
    std::cout << std::endl;
    std::cout << "Currently, the only supported auto-aliased rule is the 'auto-degree:<d>' alias," << std::endl;
    std::cout << "which refers to the least-point inner cubature rule capable of integrating" << std::endl;
    std::cout << "polynomials up to degree <d> exactly or to the highest-order cubature rule" << std::endl;
    std::cout << "available if no exact cubature rule is available. Example: the 'auto-degree:2'" << std::endl;
    std::cout << "alias refers to 'tensor:gauss-legendre:2' for Hypercube<n> shapes." << std::endl;
    std::cout << std::endl;
    std::cout << "A list of all auto-degree mappings supported by a shape can be printed by" << std::endl;
    std::cout << "specifying 'auto-degree' as the cubature name, e.g. the command" << std::endl;
    std::cout << "   cub-list h2 auto-degree" << std::endl;
    std::cout << "will print all auto-degree aliases for quadrilaterals." << std::endl;
    std::cout << std::endl << std::endl;
    std::cout << "Refined Rules:" << std::endl;
    std::cout << "--------------" << std::endl;
    std::cout << "Each Simplex and Hypercube cubature rule may also be refined by prefixing the" << std::endl;
    std::cout << "cubature name by 'refine*n:', where 'n' is a positive integer specifying the" << std::endl;
    std::cout << "desired number of refinement steps. This allows to generate cubature rules of" << std::endl;
    std::cout << "higher precision (but not higher order) from other cubature rules." << std::endl;
    std::cout << "Note: 'refine:' is identical to 'refine*1:'." << std::endl;
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
      std::cout << "ERROR: Unknown Shape '" << argv[1] << "'" << std::endl;
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
      std::cout << "ERROR: Unknown Shape '" << argv[1] << "'" << std::endl;
  }

  return 0;
}
