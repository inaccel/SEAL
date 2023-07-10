#pragma once

#include <seal/util/xcl2.hpp>

namespace seal
{
    typedef struct
    {
        cl::Buffer key1;
        cl::Buffer key2;
        cl::Buffer key3;
    } OCLPublicKey;
} // namespace seal
