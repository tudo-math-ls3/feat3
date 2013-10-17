#define FEAST_CUBATURE_TENSOR_PREFIX 1
#define FEAST_CUBATURE_SCALAR_PREFIX 1
#include <kernel/cubature/scalar/dynamic_factory.hpp>
#include <kernel/cubature/dynamic_factory.hpp>

#include <iostream>
#include <cstdio>
#include <set>
#include <map>

using namespace FEAST;
using namespace FEAST::Cubature;

typedef std::set<String> StringSet;
typedef std::map<String,String> StringMap;

#include "avail_functor.hpp"

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

void test_dynamic_scalar(const String& name)
{
  // create rule
  Scalar::Rule<> rule;
  Scalar::DynamicFactory factory(name);
  //if(!Scalar::DynamicFactory::create(rule, name))
  if(!factory.create(rule))
  {
    std::cout << "\nNo scalar cubature rule named '" << name << "' found!" << std::endl;
    return;
  }

  // print output
  std::cout << "\nPrinting scalar cubature rule '" << rule.get_name() << "':" << std::endl;
  printf("      Weight          Coord\n");
  for(Index i(0); i < rule.get_num_points(); ++i)
  {
    printf("%3i: %15.12f %15.12f\n", int(i), rule.get_weight(i), rule.get_coord(i));
  }
}

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
  std::cout << std::endl;}

template<typename Shape_>
void test_dynamic(const String& name)
{
  // create rule
  Rule<Shape_> rule;
  if(!DynamicFactory::create(rule, name))
  {
    std::cout << "\nNo " + Shape_::name() + " cubature rule named '" << name << "' found!" << std::endl;
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

  for(Index i(0); i < rule.get_num_points(); ++i)
  {
    printf("%3i: %15.12f", int(i), rule.get_weight(i));
    for(Index j(0); j < Index(Shape_::dimension); ++j)
      printf(" %15.12f", rule.get_coord(i,j));
    printf("\n");
  }
}

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
  if(argc > 2)
  {
    if(String(argv[1]).compare_no_case("-s1") == 0)
      test_dynamic< Shape::Simplex<1> >(argv[2]);
    else if(String(argv[1]).compare_no_case("-s2") == 0)
      test_dynamic< Shape::Simplex<2> >(argv[2]);
    else if(String(argv[1]).compare_no_case("-s3") == 0)
      test_dynamic< Shape::Simplex<3> >(argv[2]);
    else if(String(argv[1]).compare_no_case("-h1") == 0)
      test_dynamic< Shape::Hypercube<1> >(argv[2]);
    else if(String(argv[1]).compare_no_case("-h2") == 0)
      test_dynamic< Shape::Hypercube<2> >(argv[2]);
    else if(String(argv[1]).compare_no_case("-h3") == 0)
      test_dynamic< Shape::Hypercube<3> >(argv[2]);
  }
  else if(argc > 1)
  {
    if(String(argv[1]).compare_no_case("-s1") == 0)
      print_avail< Shape::Simplex<1> >();
    else if(String(argv[1]).compare_no_case("-s2") == 0)
      print_avail< Shape::Simplex<2> >();
    else if(String(argv[1]).compare_no_case("-s3") == 0)
      print_avail< Shape::Simplex<3> >();
    else if(String(argv[1]).compare_no_case("-h1") == 0)
      print_avail< Shape::Hypercube<1> >();
    else if(String(argv[1]).compare_no_case("-h2") == 0)
      print_avail< Shape::Hypercube<2> >();
    else if(String(argv[1]).compare_no_case("-h3") == 0)
      print_avail< Shape::Hypercube<3> >();
    else
      test_dynamic_scalar(argv[1]);
  }
  else
    print_avail_all();

  return 0;
}
