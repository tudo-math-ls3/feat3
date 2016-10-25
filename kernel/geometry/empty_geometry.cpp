#ifndef DOXYGEN

// Dummy class instance to silence ipo linker optimization warnings about empty libgeometry
class ipo_foobar_geometry
{
public:
  int i;
  ipo_foobar_geometry() :
    i(0)
  {
    (void)i;
  }
} ipo_barfoo_geometry;

#endif // DOXYGEN
