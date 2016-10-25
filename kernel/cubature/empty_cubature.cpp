#ifndef DOXYGEN

// Dummy class instance to silence ipo linker optimization warnings about empty libcubature
class ipo_foobar_cubature
{
public:
  int i;
  ipo_foobar_cubature() :
    i(0)
  {
    (void)i;
  }
} ipo_barfoo_cubature;

#endif // DOXYGEN
