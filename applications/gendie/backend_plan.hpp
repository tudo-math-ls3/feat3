// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/backend.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/cuda_util.hpp>

#include <deque>
#include <memory>
#include <vector>

namespace Gendie
{
  using namespace FEAT;

  namespace Intern
  {
    /// validates a backend option string
    inline bool is_backend(const FEAT::String& option)
    {
      return (option.compare_no_case("cuda") == 0) || (option.compare_no_case("mkl") == 0) || (option.compare_no_case("generic") == 0);
    }

    /// converts a backend option string to the corresponding PreferredBackend enum
    inline FEAT::PreferredBackend backend_from_string(const FEAT::String& option)
    {
      if(option.compare_no_case("cuda") == 0)
        return FEAT::PreferredBackend::cuda;
      if(option.compare_no_case("mkl") == 0)
        return FEAT::PreferredBackend::mkl;
      if(option.compare_no_case("generic") == 0)
        return FEAT::PreferredBackend::generic;
      XABORTM("Unknown backend " + option);
      return FEAT::PreferredBackend::generic;
    }

    /**
     * \brief Per-level compute-backend selector.
     *
     * Resolves, for a given multigrid level, which compute backend should be used. The choice can be
     * made directly (an explicit per-level list, index 0 = finest level) or heuristically (a level
     * size threshold: levels with a size metric >= \c threshold use \c heuristic_high, else
     * \c heuristic_low). If neither applies for a level, the \c fallback backend is used.
     *
     * Resolution precedence for level \c i with size metric \c level_size:
     *  1. \c per_level[i] if the list has an entry for level \c i,
     *  2. threshold heuristic if \c threshold > 0,
     *  3. \c fallback.
     */
    struct BackendSelector
    {
      /// explicit per-level backends (index 0 = finest); may be shorter than the number of levels
      std::vector<FEAT::PreferredBackend> per_level;
      /// level-size threshold enabling the heuristic (0 = disabled)
      FEAT::Index threshold = FEAT::Index(0);
      /// backend used by the heuristic for levels with size >= threshold
      FEAT::PreferredBackend heuristic_high = FEAT::PreferredBackend::cuda;
      /// backend used by the heuristic for levels with size < threshold
      FEAT::PreferredBackend heuristic_low = FEAT::PreferredBackend::generic;
      /// fallback backend when neither list nor heuristic applies
      FEAT::PreferredBackend fallback = FEAT::PreferredBackend::generic;

      /// resolves the backend for a level given a (cheap, local) size metric
      FEAT::PreferredBackend backend_for(std::size_t level, FEAT::Index level_size) const
      {
        if(level < per_level.size())
          return per_level.at(level);
        if(threshold > FEAT::Index(0))
          return (level_size >= threshold) ? heuristic_high : heuristic_low;
        return fallback;
      }

      /// returns a compact, human-readable description of this selector (for info output)
      FEAT::String display() const
      {
        FEAT::String s;
        if(!per_level.empty())
        {
          s += "list[";
          for(std::size_t i(0); i < per_level.size(); ++i)
            s += (i == 0 ? "" : " ") + FEAT::stringify(per_level.at(i));
          s += "] ";
        }
        if(threshold > FEAT::Index(0))
          s += "threshold(>=" + FEAT::stringify(threshold) + " -> " + FEAT::stringify(heuristic_high)
             + ", else " + FEAT::stringify(heuristic_low) + ") ";
        s += "fallback(" + FEAT::stringify(fallback) + ")";
        return s;
      }
    };

    /**
     * \brief Builds a BackendSelector from an explicit list of backend tokens and a threshold.
     *
     * \param[in] tokens Whitespace-parsed per-level backend tokens (index 0 = finest); validated.
     * \param[in] threshold Level-size threshold (0 disables the heuristic).
     * \param[in] fallback The fallback backend when neither list nor heuristic applies.
     */
    inline BackendSelector make_backend_selector(const std::deque<FEAT::String>& tokens,
      FEAT::Index threshold, FEAT::PreferredBackend fallback)
    {
      BackendSelector sel;
      sel.fallback = fallback;
      sel.threshold = threshold;
      for(const auto& token : tokens)
      {
        if(token.empty())
          continue;
        XASSERTM(is_backend(token), "Invalid per-level backend token '" + token + "'");
        sel.per_level.push_back(backend_from_string(token));
      }
      return sel;
    }

    /**
     * \brief Builds a BackendSelector directly from a PropertyMap.
     *
     * \param[in] prop The (root) property map to read from; may be \c nullptr.
     * \param[in] list_key The key path of the whitespace-separated per-level backend list
     *   (e.g. "backend-parameters/assembly-backends"). Empty tokens are ignored.
     * \param[in] threshold_key The key path of the level-size threshold (Index; 0 disables heuristic).
     * \param[in] fallback The fallback backend when neither list nor heuristic applies.
     */
    inline BackendSelector parse_backend_selector(const FEAT::PropertyMap* prop,
      const FEAT::String& list_key, const FEAT::String& threshold_key, FEAT::PreferredBackend fallback)
    {
      if(prop == nullptr)
      {
        BackendSelector sel;
        sel.fallback = fallback;
        return sel;
      }
      FEAT::Index threshold(0);
      prop->parse_entry(threshold_key, threshold, true);
      return make_backend_selector(prop->query(list_key, "").split_by_whitespaces(), threshold, fallback);
    }

    /**
     * \brief Solver decorator that pins the compute backend around every operation.
     *
     * Wraps an inner solver and sets the preferred backend before delegating each
     * SolverBase operation. This lets a MultiGridHierarchy (which otherwise runs the whole cycle
     * under one global backend) use a distinct backend for each level's smoother / coarse solver.
     *
     * \note The backend is set (not restored) before delegating; inter-level transfers therefore run
     * under the most-recently-pinned level backend.
     */
    template<typename Vector_>
    class BackendPinnedSolver : public FEAT::Solver::SolverBase<Vector_>
    {
    public:
      typedef Vector_ VectorType;
      typedef FEAT::Solver::SolverBase<Vector_> InnerType;

    private:
      std::shared_ptr<InnerType> _inner;
      FEAT::PreferredBackend _backend;
      /// precomputed NVTX range-name prefix, e.g. "Pin[cuda] " (built once to avoid per-call formatting)
      FEAT::String _nvtx_prefix;

    public:
      explicit BackendPinnedSolver(std::shared_ptr<InnerType> inner, FEAT::PreferredBackend backend) :
        _inner(std::move(inner)),
        _backend(backend),
        _nvtx_prefix("Pin[" + FEAT::stringify(backend) + "] ")
      {
      }

      FEAT::PreferredBackend get_backend() const
      {
        return _backend;
      }

      const std::shared_ptr<InnerType>& get_inner() const
      {
        return _inner;
      }

      virtual String name() const override
      {
        return "BackendPinned[" + FEAT::stringify(_backend) + "](" + _inner->name() + ")";
      }

      virtual void init_symbolic() override
      {
        FEAT::Util::ScopedNvtxRange _nvtx((_nvtx_prefix + "initS").c_str());
        FEAT::Backend::set_preferred_backend(_backend);
        _inner->init_symbolic();
      }

      virtual void init_numeric() override
      {
        FEAT::Util::ScopedNvtxRange _nvtx((_nvtx_prefix + "initN").c_str());
        FEAT::Backend::set_preferred_backend(_backend);
        _inner->init_numeric();
      }

      virtual void done_numeric() override
      {
        FEAT::Util::ScopedNvtxRange _nvtx((_nvtx_prefix + "doneN").c_str());
        FEAT::Backend::set_preferred_backend(_backend);
        _inner->done_numeric();
      }

      virtual void done_symbolic() override
      {
        FEAT::Util::ScopedNvtxRange _nvtx((_nvtx_prefix + "doneS").c_str());
        FEAT::Backend::set_preferred_backend(_backend);
        _inner->done_symbolic();
      }

      virtual FEAT::Solver::Status apply(Vector_& vec_cor, const Vector_& vec_def) override
      {
        FEAT::Util::ScopedNvtxRange _nvtx((_nvtx_prefix + "apply").c_str());
        FEAT::Backend::set_preferred_backend(_backend);
        return _inner->apply(vec_cor, vec_def);
      }
    }; // class BackendPinnedSolver

    /// creates a BackendPinnedSolver wrapping the given inner solver
    template<typename Vector_>
    inline std::shared_ptr<BackendPinnedSolver<Vector_>> new_backend_pinned_solver(
      std::shared_ptr<FEAT::Solver::SolverBase<Vector_>> inner, FEAT::PreferredBackend backend)
    {
      return std::make_shared<BackendPinnedSolver<Vector_>>(std::move(inner), backend);
    }
  } // namespace Intern
} // namespace Gendie
