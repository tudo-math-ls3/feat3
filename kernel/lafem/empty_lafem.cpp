#ifndef DOXYGEN

// dummy class instance to silence ipo linker optimization warnings about empty liblafem
class ipo_foobar_lafem
{
public:
  int i;
  ipo_foobar_lafem() :
    i(0)
  {
    (void)i;
  }
} ipo_barfoo_lafem;

#endif // DOXYGEN
