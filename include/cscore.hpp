#include <iopp/file_input_stream.hpp>
#include "internal/util/si_iec_literals.hpp"

namespace zk {

double cscore(iopp::FileInputStream& in, size_t const pattern_len, size_t const sampling_exp, size_t const block_size = 64_Ki, size_t const buffer_size = 32_Mi);

}
