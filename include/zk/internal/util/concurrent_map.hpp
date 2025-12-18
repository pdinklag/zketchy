#pragma once

#include <algorithm>
#include <concepts>

#include <allocator/alignedallocator.hpp>
#include <data-structures/hash_table_mods.hpp>
#include <utils/hash/murmur2_hash.hpp>
#include <data-structures/table_config.hpp>

namespace zk::internal {

template<typename Key, typename Value>
using ConcurrentMap = typename growt::table_config<Key, Value, utils_tm::hash_tm::murmur2_hash, growt::AlignedAllocator<>, hmod::growable>::table_type;

}
