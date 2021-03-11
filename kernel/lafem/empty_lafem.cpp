// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
