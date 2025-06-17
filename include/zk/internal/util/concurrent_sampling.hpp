#pragma once

#include <algorithm>
#include <concepts>

#include <allocator/alignedallocator.hpp>
#include <data-structures/hash_table_mods.hpp>
#include <utils/hash/murmur2_hash.hpp>
#include <data-structures/table_config.hpp>

namespace zk::internal {

template<std::unsigned_integral Fingerprint, std::unsigned_integral Index>
class ConcurrentSampling {
public:
    using Map = typename growt::table_config<Fingerprint, Index, utils_tm::hash_tm::murmur2_hash, growt::AlignedAllocator<>, hmod::growable>::table_type;

    struct UpdateLeftmost {
        using mapped_type = Index;

        mapped_type operator()(mapped_type& lhs, const mapped_type& rhs) const { return lhs = std::min(lhs, rhs); }
    };
};

}
