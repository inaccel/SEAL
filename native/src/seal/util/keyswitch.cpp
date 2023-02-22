#include "seal/util/keyswitch.h"
#include <stdexcept>
#include "hexl/hexl.hpp"

#define BIT_MASK(BITS) ((1UL << BITS) - 1)
#define MAX_DEGREE 16384

using namespace std;
using namespace intel::hexl;

namespace seal
{
    namespace util
    {
        namespace keyswitch
        {
            void ComputeRootOfUnityPowers(
                uint64_t m_q, uint64_t m_degree, uint64_t m_degree_bits, uint64_t m_w,
                uint64_t *inv_root_of_unity_powers, uint64_t *precon64_inv_root_of_unity_powers,
                uint64_t *root_of_unity_powers, uint64_t *precon64_root_of_unity_powers)
            {
                uint64_t inv_root_of_unity_powers_pre[MAX_DEGREE];

                // 64-bit preconditioning
                root_of_unity_powers[0] = 1;
                inv_root_of_unity_powers_pre[0] = 1;
                uint64_t idx = 0;
                uint64_t prev_idx = idx;

                for (size_t i = 1; i < m_degree; i++)
                {
                    idx = ReverseBits(i, m_degree_bits);
                    root_of_unity_powers[idx] = MultiplyMod(root_of_unity_powers[prev_idx], m_w, m_q);
                    inv_root_of_unity_powers_pre[idx] = InverseMod(root_of_unity_powers[idx], m_q);

                    prev_idx = idx;
                }

                precon64_root_of_unity_powers[0] = 0;
                for (size_t i = 1; i < m_degree; i++)
                {
                    precon64_root_of_unity_powers[i] = MultiplyFactor(root_of_unity_powers[i], 64, m_q).BarrettFactor();
                }

                idx = 0;

                for (size_t m = (m_degree >> 1); m > 0; m >>= 1)
                {
                    for (size_t i = 0; i < m; i++)
                    {
                        inv_root_of_unity_powers[idx] = inv_root_of_unity_powers_pre[m + i];
                        idx++;
                    }
                }

                inv_root_of_unity_powers[m_degree - 1] = 0;

                for (uint64_t i = 0; i < m_degree; i++)
                {
                    precon64_inv_root_of_unity_powers[i] =
                        MultiplyFactor(inv_root_of_unity_powers[i], 64, m_q).BarrettFactor();
                }
            }

            void loadTwiddleFactors(
                uint64_t n, uint64_t key_modulus_size, const Modulus *moduli, uint64_t *root_of_unity_powers_ptr)
            {
                for (uint64_t i = 0; i < key_modulus_size; i++)
                {
                    ComputeRootOfUnityPowers(
                        moduli[i].value(), n, Log2(n), MinimalPrimitiveRoot(2 * n, moduli[i].value()),
                        root_of_unity_powers_ptr + i * n * 4, root_of_unity_powers_ptr + i * n * 4 + n,
                        root_of_unity_powers_ptr + i * n * 4 + n * 2, root_of_unity_powers_ptr + i * n * 4 + n * 3);
                }
            }

            uint64_t precomputeModulusK(uint64_t modulus)
            {
                uint64_t k = 0;
                for (uint64_t i = 64; i > 0; i--)
                {
                    if ((1UL << i) >= modulus)
                    {
                        k = i;
                    }
                }
                return k;
            }

            void buildModulusMeta(
                uint64_t key_modulus_size, const Modulus *moduli, const MultiplyUIntModOperand *modswitch_factors,
                moduli_t *modulus_meta)
            {
                for (uint64_t i = 0; i < key_modulus_size; i++)
                {
                    (*modulus_meta).data[i][0] = moduli[i].value();
                    (*modulus_meta).data[i][1] = MultiplyFactor(1, 64, moduli[i].value()).BarrettFactor();
                    uint64_t modulus = moduli[i].value();
                    uint64_t twice_modulus = 2 * modulus;
                    uint64_t four_times_modulus = 4 * modulus;
                    uint64_t arg2 = modswitch_factors[i].operand;
                    const int InputModFactor = 8;
                    arg2 = ReduceMod<InputModFactor>(arg2, modulus, &twice_modulus, &four_times_modulus);
                    (*modulus_meta).data[i][2] = arg2;
                    uint64_t k = precomputeModulusK(moduli[i].value());
                    __extension__ __int128 a = 1;
                    uint64_t r = (uint64_t)((a << (2 * k)) / moduli[i].value());
                    (*modulus_meta).data[i][3] = (r << 8) | k;
                }
            }

            void buildInvnMeta(
                uint64_t n, uint64_t key_modulus_size, const Modulus *moduli, const uint64_t *root_of_unity_powers_ptr,
                moduli_t *invn)
            {
                for (uint64_t i = 0; i < key_modulus_size; i++)
                {
                    uint64_t inv_n = InverseMod(n, moduli[i].value());
                    uint64_t W_op = root_of_unity_powers_ptr[i * n * 4 + n - 1];
                    uint64_t inv_nw = MultiplyMod(inv_n, W_op, moduli[i].value());
                    uint64_t y_barrett_n = DivideUInt128UInt64Lo(inv_n, 0, moduli[i].value());
                    uint64_t y_barrett_nw = DivideUInt128UInt64Lo(inv_nw, 0, moduli[i].value());
                    (*invn).data[i][0] = inv_n;
                    unsigned long k = precomputeModulusK(moduli[i].value());
                    __extension__ __int128 a = 1;
                    uint64_t r = (uint64_t)((a << (2 * k)) / moduli[i].value());
                    (*invn).data[i][1] = (r << 8) | k;
                    (*invn).data[i][2] = y_barrett_n;
                    (*invn).data[i][3] = y_barrett_nw;
                }
            }

            void loadKeys(
                uint64_t n, uint64_t decomp_modulus_size, uint64_t key_modulus_size, const uint64_t **k_switch_keys,
                DyadmultKeys1_t *key_vector1, DyadmultKeys2_t *key_vector2, DyadmultKeys3_t *key_vector3)
            {
                size_t kv_idx = 0;
                for (uint64_t k = 0; k < decomp_modulus_size; k++)
                {
                    for (uint64_t j = 0; j < n; j++, kv_idx++)
                    {
                        for (uint64_t i = 0; i < key_modulus_size; i++)
                        {
                            uint64_t key1 = k_switch_keys[k][i * n + j];
                            uint64_t key2 = k_switch_keys[k][(i + key_modulus_size) * n + j];
                            if (i == 0)
                            {
                                key_vector1[kv_idx].key1 = key1 & BIT_MASK(52);
                                key_vector1[kv_idx].key2 = key2 & BIT_MASK(52);
                            }
                            else if (i == 1)
                            {
                                key_vector1[kv_idx].key3 = key1 & BIT_MASK(52);
                                key_vector1[kv_idx].key4 = key2 & BIT_MASK(52);
                            }
                            else if (i == 2)
                            {
                                key_vector1[kv_idx].key5 = key1 & BIT_MASK(48);
                                key_vector2[kv_idx].key1 = (key1 >> 48) & BIT_MASK(4);
                                key_vector2[kv_idx].key2 = key2 & BIT_MASK(52);
                            }
                            else if (i == 3)
                            {
                                key_vector2[kv_idx].key3 = key1 & BIT_MASK(52);
                                key_vector2[kv_idx].key4 = key2 & BIT_MASK(52);
                            }
                            else if (i == 4)
                            {
                                key_vector2[kv_idx].key5 = key1 & BIT_MASK(52);
                                key_vector2[kv_idx].key6 = key2 & BIT_MASK(44);
                                key_vector3[kv_idx].key1 = (key2 >> 44) & BIT_MASK(8);
                            }
                            else if (i == 5)
                            {
                                key_vector3[kv_idx].key2 = key1 & BIT_MASK(52);
                                key_vector3[kv_idx].key3 = key2 & BIT_MASK(52);
                            }
                            else if (i == 6)
                            {
                                key_vector3[kv_idx].key4 = key1 & BIT_MASK(52);
                                key_vector3[kv_idx].key5 = key2 & BIT_MASK(52);
                                key_vector3[kv_idx].NOT_USED = 0;
                            }
                            else
                            {
                                throw invalid_argument("keys not supported");
                            }
                        }
                    }
                }
            }

            void readOutput(
                uint64_t n, uint64_t decomp_modulus_size, const Modulus *moduli, const uint64_t *fpga_result,
                uint64_t *result, uint64_t batch_size)
            {
                size_t output_size = n * decomp_modulus_size * 2;
                for (size_t off = 0; off < batch_size * output_size; off += output_size)
                {
                    for (size_t i = 0; i < decomp_modulus_size; i++)
                    {
                        uint64_t modulus = moduli[i].value();
                        for (size_t j = 0; j < n; j++)
                        {
                            size_t k = i * n + j;
                            result[off + k] += fpga_result[off + 2 * k];
                            if (result[off + k] >= modulus)
                                result[off + k] -= modulus;

                            size_t idx = k + n * decomp_modulus_size;
                            result[off + idx] += fpga_result[off + 2 * k + 1];
                            if (result[off + idx] >= modulus)
                                result[off + idx] -= modulus;
                        }
                    }
                }
            }
        } // namespace keyswitch
    } // namespace util
} // namespace seal
