#pragma once
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/string.hpp>
namespace Gendie
{
  #define PARSE_PROP_OPTION(prop_map, key, arg, not_found_rtn) XASSERTM(prop_map->parse_entry(key, arg, not_found_rtn), FEAT::String("Could not parse ") + key)

  inline bool check_for_config_option(const std::pair<FEAT::String, bool>& config_query, bool default_rtn = false)
  {
    if(!config_query.second)
      return default_rtn;
    return !(config_query.first.empty() || (config_query.first.compare_no_case(FEAT::String("false")) == 0) || (config_query.first.compare_no_case(FEAT::String("off")) == 0));
  }

  inline FEAT::PreferredBackend parse_backend(const FEAT::String& backend_string)
  {
    if(backend_string.compare_no_case("cuda") == 0)
    {
      return FEAT::PreferredBackend::cuda;
    }
    else if(backend_string.compare_no_case("mkl") == 0)
    {
      return FEAT::PreferredBackend::mkl;
    }
    else if(backend_string.compare_no_case("generic") == 0)
    {
      return FEAT::PreferredBackend::generic;
    }
    XABORTM("Unknown backend " + backend_string);
    return FEAT::PreferredBackend::generic;

  }
}