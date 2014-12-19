#ifndef DOXYGEN

// dummy class instance to silence ipo linker optimization warnings about empty liblafem
class ipo_foobar
{
public:
  int i;
  ipo_foobar() :
    i(0)
  {
    (void)i;
  }
} ipo_barfoo;

#endif // DOXYGEN
