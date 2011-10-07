// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/pretty_printer.hpp>
#include <kernel/util/assertion.hpp>
#include <test_system/test_system.hpp>

// includes, system
#include <iostream>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
* \brief testing the pretty printer
*
* \test
* This tests the functionality of the pretty printer.
*
* \tparam Tag_
* description missing
*
* \tparam DT_
* description missing
*
* \author Hilmar Wobker
*/
template<
  typename Tag_,
  typename DT_>
class PrettyPrinterTest
  : public TaggedTest<Tag_, DT_>
{

public:

  /// CTOR
  PrettyPrinterTest()
    : TaggedTest<Tag_, DT_>("pretty_printer_test")
  {
  }


  /// main routine
  void run() const
  {
    CONTEXT("PrettyPrinterTest::run()");

    // test the PrettyPrinter
    String prefix(String("Proc 42"));
    PrettyPrinter pp(40, '#', prefix + " ");
    pp.add_line_sep();
    pp.add_line_centered("Testing the pretty printer");
    pp.add_line_sep();
    pp.add_line("left bla blub");
    pp.add_line("too long too long too long too long too long too long");
    pp.add_line_no_right_delim("also too long too long too long too long too long too long");
    pp.add_line("left bla blub");
    pp.add_line_sep();

    // print it like this
    std::cout << pp.block();

    // modify the pretty printer a little bit:
    pp.reset_block();
    pp.set_width(60);
    pp.set_prefix("WUPP ");
    pp.set_delim('*');
    pp.add_line_sep();
    pp.add_line_centered("Testing the pretty printer once again");
    pp.add_line_sep();
    pp.add_line_centered("bla");
    pp.add_line_centered("bla bla");
    pp.add_line_centered("bla bla bla");
    pp.add_line_centered("bla bla");
    pp.add_line_centered("bla");
    pp.add_line_sep();

    // or print it like this
    pp.print(std::cout);
  } // run()
}; // PrettyPrinterTest

// create test instance
PrettyPrinterTest<Nil, Nil> pretty_printer_test;
