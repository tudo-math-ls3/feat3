#ifndef DOXYGEN

// dummy class instance to silence ipo linker optimization warnings about empty liblafem_arch
class ipo_foobar_lafem_arch
{
public:
  int i;
  ipo_foobar_lafem_arch() :
    i(0)
  {
    (void)i;
  }
} ipo_barfoo_lafem_arch;

#endif // DOXYGEN
