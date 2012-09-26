#define FEAST_CUBATURE_TENSOR_PREFIX
#include <kernel/cubature/scalar/dynamic_factory.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <iostream>

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

void test_scalar()
{
  typedef Scalar::Rule<> RuleType;
  typedef Scalar::DynamicFactory<> Factory;
  typedef Scalar::DriverFactory<Scalar::PulcherrimaDriver> Pulcherrima;
  typedef Scalar::DriverFactory<Scalar::GaussLegendreDriver> GaussLegendre;

  RuleType pulcherrima_1(Pulcherrima::create());
  RuleType gauss2_1(GaussLegendre::create(2));

  RuleType pulcherrima_2(Factory::create("pulcherrima"));
  RuleType gauss2_2(Factory::create("gauss-legendre:2"));
  RuleType midpoint(Factory::create("midpoint"));
}

void test_quad()
{
  typedef Rule<Shape::Quadrilateral> RuleType;
  typedef DynamicFactory<Shape::Quadrilateral> Factory;
  typedef DriverFactory<BarycentreDriver,Shape::Quadrilateral> Barycentre;
  typedef DriverFactory<TrapezoidalDriver,Shape::Quadrilateral> Trapezoidal;

  RuleType barycentre_1(Barycentre::create());
  RuleType trapezoidal_1(Trapezoidal::create());

  RuleType barycentre_2(Factory::create("barycentre"));
  RuleType trapezoidal_2(Factory::create("trapezoidal"));
}

void test_tria()
{
  typedef Rule<Shape::Triangle> RuleType;
  typedef DynamicFactory<Shape::Triangle> Factory;
  typedef DriverFactory<BarycentreDriver,Shape::Triangle> Barycentre;
  typedef DriverFactory<TrapezoidalDriver,Shape::Triangle> Trapezoidal;

  RuleType barycentre_1(Barycentre::create());
  RuleType trapezoidal_1(Trapezoidal::create());

  RuleType barycentre_2(Factory::create("barycentre"));
  RuleType trapezoidal_2(Factory::create("trapezoidal"));
}

int main(int argc, char* argv[])
{
  print_avail_all();
  return 0;
}
