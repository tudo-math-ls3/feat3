#define FEAST_CUBATURE_TENSOR_PREFIX 1
#include <kernel/cubature/scalar/dynamic_factory.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <iostream>
#include <cstdio>

#define SEP std::cout << "************************************************************" \
  "************************************************************" << std::endl

#define SEM std::cout << "------------------------------------------------------------" \
  "------------------------------------------------------------" << std::endl

#define SEC(x) SEP; std::cout << (x) << std::endl; SEM

using namespace FEAST;
using namespace FEAST::Cubature;

void print_avail_scalar(bool aliases = true)
{
  std::cout << "Available scalar cubature names:" << std::endl;
  Scalar::DynamicFactory<>::print_avail(true);
  std::cout << std::endl;
}

template<typename Shape_>
void print_avail()
{
  typedef DynamicFactory<Shape_> Factory;

  std::cout << "Available cubature names for '" << Shape_::name() << "' :" << std::endl;
  DynamicFactory<Shape_>::print_avail(true);
  std::cout << std::endl;
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

void test_dynamic_scalar(const String& name)
{
  // create rule
  Scalar::Rule<> rule;
  if(!Scalar::DynamicFactory<>::create(rule, name))
  {
    std::cout << "\nNo scalar cubature rule named '" << name << "' found!" << std::endl;
    return;
  }

  // print output
  std::cout << "\nPrinting scalar cubature rule '" << rule.get_name() << "':" << std::endl;
  for(Index i(0); i < rule.get_num_points(); ++i)
  {
    printf("%3i: %15.12f %15.12f\n", i, rule.get_weight(i), rule.get_coord(i));
  }
}

template<typename Shape_>
void test_dynamic(const String& name)
{
  // create rule
  Rule<Shape_> rule;
  if(!DynamicFactory<Shape_>::create(rule, name))
  {
    std::cout << "\nNo " + Shape_::name() + " cubature rule named '" << name << "' found!" << std::endl;
    return;
  }

  // print output
  std::cout << "\nPrinting " + Shape_::name() + " cubature rule '" << rule.get_name() << "':" << std::endl;
  for(Index i(0); i < rule.get_num_points(); ++i)
  {
    printf("%3i: %15.12f", i, rule.get_weight(i));
    for(Index j(0); j < Shape_::dimension; ++j)
      printf(" %15.12f", rule.get_coord(i,j));
    printf("\n");
  }
}

void test_dynamic_all(const String& name)
{
  test_dynamic_scalar(name);
  test_dynamic< Shape::Simplex<1> >(name);
  test_dynamic< Shape::Simplex<2> >(name);
  test_dynamic< Shape::Simplex<3> >(name);
  test_dynamic< Shape::Hypercube<1> >(name);
  test_dynamic< Shape::Hypercube<2> >(name);
  test_dynamic< Shape::Hypercube<3> >(name);
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
    test_dynamic_scalar(argv[1]);
  else
    print_avail_all();

  //test_dynamic< Shape::Hypercube<2> >("barycentre");
  //test_dynamic< Shape::Hypercube<2> >("refine*2:barycentre");

  return 0;
}