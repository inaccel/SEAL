#pragma once

#include "seal/modulus.h"
#include "seal/util/uintarithsmallmod.h"
#include <cstdint>

typedef struct
{
    uint64_t key1 : 52;
    uint64_t key2 : 52;
    uint64_t key3 : 52;
    uint64_t key4 : 52;
    uint64_t key5 : 48;
} __attribute__((packed)) DyadmultKeys1_t;

/// @brief
/// Struct DyadmultKeys2_t
/// @param[in] key1-6 stores the bits of compressed switch key data
typedef struct
{
    uint64_t key1 : 4;
    uint64_t key2 : 52;
    uint64_t key3 : 52;
    uint64_t key4 : 52;
    uint64_t key5 : 52;
    uint64_t key6 : 44;
} __attribute__((packed)) DyadmultKeys2_t;

/// @brief
/// Struct DyadmultKeys3_t
/// @param[in] key1-5 stores the bits of compressed switch key data
typedef struct
{
    uint64_t key1 : 8;
    uint64_t key2 : 52;
    uint64_t key3 : 52;
    uint64_t key4 : 52;
    uint64_t key5 : 52;
    uint64_t NOT_USED : 40;
} __attribute__((packed)) DyadmultKeys3_t;

typedef struct
{
    uint64_t data[8][4];
} moduli_t;

namespace seal
{
    namespace util
    {
        namespace keyswitch
        {
            void loadTwiddleFactors(
                uint64_t n, uint64_t key_modulus_size, const Modulus *moduli, uint64_t *root_of_unity_powers_ptr);

            void buildModulusMeta(
                uint64_t key_modulus_size, const Modulus *moduli, const MultiplyUIntModOperand *modswitch_factors,
                moduli_t *modulus_meta);

            void buildInvnMeta(
                uint64_t n, uint64_t key_modulus_size, const Modulus *moduli, const uint64_t *root_of_unity_powers_ptr,
                moduli_t *invn);

            void loadKeys(
                uint64_t n, uint64_t decomp_modulus_size, uint64_t key_modulus_size, const uint64_t **k_switch_keys,
                DyadmultKeys1_t *key_vector1, DyadmultKeys2_t *key_vector2, DyadmultKeys3_t *key_vector3);

            void readOutput(
                uint64_t n, uint64_t decomp_modulus_size, const Modulus *moduli, const uint64_t *fpga_result,
                uint64_t *result, uint64_t batch_size);
        } // namespace keyswitch
    } // namespace util
} // namespace seal
