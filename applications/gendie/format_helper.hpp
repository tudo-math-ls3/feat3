#pragma once
#include <kernel/util/string.hpp>

namespace Gendie
{
  inline FEAT::String format_submemory_mm_t(const FEAT::String& s, unsigned long long total, unsigned long long mmin, unsigned long long mmax, unsigned long long padlen)
  {
    double mb = double(total) / (1024.0*1024.0*1024.0);
    double mmi = double(mmin) / (1024.0*1024.0*1024.0);
    double mma = double(mmax) / (1024.0*1024.0*1024.0);
    return s.pad_back(padlen, '.') + ": " + FEAT::stringify_fp_fix(mb, 3, 12) + " GB { " + FEAT::stringify_fp_fix(mmi, 3, 12) + " GB / " + FEAT::stringify_fp_fix(mma, 3, 12) + " GB }\n";
  }

  inline FEAT::String format_submemory_mm_t_kb(const FEAT::String& s, unsigned long long total, unsigned long long mmin, unsigned long long mmax, unsigned long long padlen)
  {
    double mb = double(total) / (1024.0);
    double mmi = double(mmin) / (1024.0);
    double mma = double(mmax) / (1024.0);
    return s.pad_back(padlen, '.') + ": " + FEAT::stringify_fp_fix(mb, 3, 12) + " KB { " + FEAT::stringify_fp_fix(mmi, 3, 12) + " KB / " + FEAT::stringify_fp_fix(mma, 3, 12) + " KB }\n";
  }

  inline FEAT::String format_subtime_mm(const FEAT::String& s, double t, double total, double tmin, double tmax, unsigned long long padlen)
  {
    return s.pad_back(padlen, '.') + ": " + FEAT::stringify_fp_fix(t, 3, 10) + " sec [" + FEAT::stringify_fp_fix(100.0*t/total, 3, 7)
      + "% ] { " + FEAT::stringify_fp_fix(tmin, 3, 10) + " / " + FEAT::stringify_fp_fix(tmax, 3, 10) + " }\n";
  }

  template <typename IT_>
  inline FEAT::String format_index_mm(const FEAT::String& s, IT_ total, IT_ tmin, IT_ tmax, unsigned long long padlen)
  {
    return s.pad_back(padlen, '.') + ": " + FEAT::stringify(total) + " { " + FEAT::stringify(tmin) + " / " + FEAT::stringify(tmax) + " }\n";
  }
}