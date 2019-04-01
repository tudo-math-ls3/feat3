// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
