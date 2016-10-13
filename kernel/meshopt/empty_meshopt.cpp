#ifndef DOXYGEN

// Dummy class instance to silence ipo linker optimization warnings about empty libmeshopt
class ipo_foobar_meshopt
{
public:
  int i;
  ipo_foobar_meshopt() :
    i(0)
  {
    (void)i;
  }
} ipo_barfoo_meshopt;

#endif // DOXYGEN
