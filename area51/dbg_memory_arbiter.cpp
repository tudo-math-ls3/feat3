#include <kernel/runtime.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/memory_aux.hpp>

using namespace FEAT;

void dump_arbiter(Memory::Arbiter& arbiter)
{
  Memory::TypedView<unsigned char> v(arbiter.view(Memory::Location::main, Memory::Access::read));
  //const unsigned char* c = (const unsigned char*)v.raw_r();
  const unsigned char* c = v.get_r();
  for(std::size_t i(0); i < v.bytes(); ++i, ++c)
  {
    printf("%02X", int(*c));
    if(((i+1u) % 32) == 0)
      printf("\n");
    else if(((i+1u) % 8) == 0)
      printf(" ");
  }
  printf("\n");
}

void dump_mapped()
{
  printf("Rsrv: %10zu Bytes\nMain: %10zu Bytes\nCUDA: %10zu Bytes\n\n",
    Memory::Manager::bytes_reserved(),
    Memory::Manager::bytes_allocated_main(),
    Memory::Manager::bytes_allocated_cuda());
}

int main(int argc, char** argv)
{
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  Memory::Arbiter ma(128);



  {
    Memory::View v = ma.view(Memory::Location::main, Memory::Access::read_write);

    Memory::memset_main(v.raw_w(), (int)0xDEADBEEF, v.bytes());
  }

  dump_mapped();

  ma.format(Memory::Location::main);

  dump_arbiter(ma);
  dump_mapped();

  ma.format_random(Memory::Location::main);

  dump_arbiter(ma);
  dump_mapped();

#ifdef FEAT_HAVE_CUDA
  {
    Memory::Arbiter va = ma.attach(36, 48);
    Memory::View v = va.view(Memory::Location::cuda, Memory::Access::read_write);
    v.format_to_zero();
  }
#endif

  dump_arbiter(ma);
  dump_mapped();

  {
    Memory::Arbiter va = ma.attach(24, 16);
    Memory::TypedView<std::uint32_t> v = va.typed_view<std::uint32_t>(Memory::Location::main, Memory::Access::read_write);
    for(std::size_t i = 0; i < v.size(); ++i)
      v[i] = 0xBEBADEC0;
  }

  {
    Memory::Arbiter va = ma.attach(40, 16);
    Memory::TypedView<std::uint32_t> v = va.typed_view<std::uint32_t>(Memory::Location::main, Memory::Access::read_write);
    for(std::size_t i = 0; i < v.size(); ++i)
      v[i] = 0xEFBEADDE;
  }

  Memory::Manager::free_all_unmapped(Memory::Location::cuda);

  dump_arbiter(ma);
  dump_mapped();

  {
    Memory::Arbiter va1 = ma.attach(32, 16);
    Memory::Arbiter va2 = ma.attach(8, 72);
    Memory::TypedView<std::uint32_t> v1 = va1.typed_view<std::uint32_t>(Memory::Location::main, Memory::Access::read);
    Memory::TypedView<std::uint64_t> v2 = va2.typed_view<std::uint64_t>(Memory::Location::main, Memory::Access::read);
  }

  //Memory::Arbiter::release(h);
  ma.release();

  dump_mapped();

  return 0;
}
