#include <kernel/util/assertion.hpp>
#include <kernel/util/runtime.hpp>

#include <cstdio>

namespace FEAT
{
  void assertion(
    bool expr,
    const char * const expr_str,
    const char * const func,
    const char * const file,
    const int line,
    const char * const msg)
  {
    // alright?
    if(expr)
      return;

    // write error message if available
    if(msg != nullptr)
      fprintf(stderr, "\nERROR: ASSERTION FAILED: %s\n", msg);
    else
      fprintf(stderr, "\nERROR: ASSERTION FAILED\n");

    // write basic information
    fprintf(stderr, "Expression: %s\n", expr_str);
    fprintf(stderr, "Function..: %s\n", func);
    fprintf(stderr, "File......: %s\n", file);
    fprintf(stderr, "Line......: %i\n", line);

    // flush stderr
    fflush(stderr);

    // abort execution;
    // this may also write a call-stack dump if possible
    Runtime::abort();
  }
} // namespace FEAT
