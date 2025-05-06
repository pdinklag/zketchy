#pragma once

#include <pm.hpp>

namespace zk::internal {

#ifdef BENCHMARK
    using Result = pm::Result;
    using MemoryTimePhase = pm::MemoryTimePhase;
    using TimePhase = pm::TimePhase;

    constexpr bool do_benchmark = true;
#else
    using Result = pm::NoopResult;
    using MemoryTimePhase = pm::NoopPhase;
    using TimePhase = pm::NoopPhase;

    constexpr bool do_benchmark = false;
#endif

}
