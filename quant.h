#pragma once

#include <vector>
#include <cstdint>
#include <optional>
#include <array>

#include "sound_unit.h"

namespace atrac3 {

    constexpr std::array<uint8_t, 7> ATRAC3_HUFF_TAB_SIZES = { 9, 5, 7, 9, 15, 31, 63 };
    constexpr std::array<uint8_t, 8> ATRAC3_CLC_LENGTH_TAB = { 0, 4, 3, 3, 4, 4, 5, 6 };
    constexpr std::array<int8_t, 4> ATRAC3_MANTISSA_CLC_TAB = { 0, 1, -2, -1 };
    constexpr std::array<int8_t, 18> ATRAC3_MANTISSA_VLC_TAB = { 0, 0, 0, 1, 0, -1, 1, 0, -1, 0, 1, 1, 1, -1, -1, 1, -1, -1 };
    constexpr std::array<std::pair<uint8_t, uint8_t>, 139> ATRAC3_HUFF_TABS = { {
        {31, 1}, {32, 3}, {33, 3}, {34, 4}, {35, 4}, {36, 5}, {37, 5}, {38, 5}, {39, 5},
        {31, 1}, {32, 3}, {30, 3}, {33, 3}, {29, 3},
        {31, 1}, {32, 3}, {30, 3}, {33, 4}, {29, 4}, {34, 4}, {28, 4},
        {31, 1}, {32, 3}, {30, 3}, {33, 4}, {29, 4}, {34, 5}, {28, 5}, {35, 5}, {27, 5},
        {31, 2}, {32, 3}, {30, 3}, {33, 4}, {29, 4}, {34, 4}, {28, 4}, {38, 4}, {24, 4},
        {35, 5}, {27, 5}, {36, 6}, {26, 6}, {37, 6}, {25, 6},
        {31, 3}, {32, 4}, {30, 4}, {33, 4}, {29, 4}, {34, 4}, {28, 4}, {46, 4}, {16, 4},
        {35, 5}, {27, 5}, {36, 5}, {26, 5}, {37, 5}, {25, 5}, {38, 6}, {24, 6}, {39, 6},
        {23, 6}, {40, 6}, {22, 6}, {41, 6}, {21, 6}, {42, 7}, {20, 7}, {43, 7}, {19, 7},
        {44, 7}, {18, 7}, {45, 7}, {17, 7},
        {31, 3}, {62, 4}, {0, 4}, {32, 5}, {30, 5}, {33, 5}, {29, 5}, {34, 5}, {28, 5},
        {35, 5}, {27, 5}, {36, 5}, {26, 5}, {37, 6}, {25, 6}, {38, 6}, {24, 6}, {39, 6},
        {23, 6}, {40, 6}, {22, 6}, {41, 6}, {21, 6}, {42, 6}, {20, 6}, {43, 6}, {19, 6},
        {44, 6}, {18, 6}, {45, 7}, {17, 7}, {46, 7}, {16, 7}, {47, 7}, {15, 7}, {48, 7},
        {14, 7}, {49, 7}, {13, 7}, {50, 7}, {12, 7}, {51, 7}, {11, 7}, {52, 8}, {10, 8},
        {53, 8}, {9, 8}, {54, 8}, {8, 8}, {55, 8}, {7, 8}, {56, 8}, {6, 8}, {57, 8},
        {5, 8}, {58, 8}, {4, 8}, {59, 8}, {3, 8}, {60, 8}, {2, 8}, {61, 8}, {1, 8}
    } };

    constexpr std::array<float, 8> ATRAC3_INV_MAX_QUANT = {
        0.0f,
        1.0f / 1.5f,
        1.0f / 2.5f,
        1.0f / 3.5f,
        1.0f / 4.5f,
        1.0f / 7.5f,
        1.0f / 15.5f,
        1.0f / 31.5f,
    };

    constexpr std::array<size_t, 33> ATRAC3_SUBBAND_TAB = {
        0, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160, 176, 192, 224, 256, 288, 320,
        352, 384, 416, 448, 480, 512, 576, 640, 704, 768, 896, 1024,
    };

    constexpr std::array<int32_t, 32> ENERGY_THRESHOLD = {
        7, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 4, 5, 5, 5, 6,
        6, 7, 7, 8, 10, 13, 17, 22, 28, 35, 49, 74, 109, 155, 250, 441,
    };

    constexpr std::array<int32_t, 8> CODEBOOK_BOUNDARY = { 1, 2, 2, 2, 4, 6, 6, 40 };
    constexpr std::array<int32_t, 8> COST_REDUCTION = { 6, 40, 40, 60, 76, 60, 60, 100 };
    constexpr std::array<int32_t, 8> BASE_COST_SCALAR = { 100, 15, 20, 25, 29, 35, 45, 55 };
    constexpr std::array<int32_t, 8> SF_THRESHOLD_OFFSET = { 55, 3, 5, 7, 9, 12, 15, 18 };

    constexpr std::array<int32_t, 18> LF_OVERSHOOT_WEIGHTS_Q8 = {
        256, 188, 97, 97, 97, 97, 97, 97,
        74, 74, 74, 74, 74, 74, 74, 74, 74, 74
    };

    constexpr std::array<int32_t, 5> TONAL_SPREAD_Q0 = { -225, -266, -307, -317, -1024 };

    struct QuantizedSubband {
        uint8_t table_index = 0;
        bool has_scale_factor = false;
        uint8_t scale_factor_index = 0;
        std::vector<int8_t> mantissas;
        size_t payload_bits = 0;
        float mse = 0.0f;
        float max_abs_err = 0.0f;

        std::vector<float> dequantized(size_t len) const;
        static QuantizedSubband uncoded(const float* coeffs, size_t len);
    };

    struct SearchOptions {
        float lambda = 0.0001f;
        size_t max_candidates_per_band = 64;
        bool use_fixed_schedule = false;
        std::optional<size_t> target_bits;
        std::array<bool, 32> tonal_marked_subbands = {};
    };

    struct TonalExtractionResult {
        std::vector<TonalComponent> tonal_components;
        TonalCodingModeSelector tonal_mode_selector = TonalCodingModeSelector::AllVlc;
        size_t tonal_bits = 0;
        uint8_t coded_qmf_bands = 1;
        std::array<bool, 32> tonal_subbands = {};
    };

    struct SpectrumEncoding {
        SpectralUnit spectral_unit;
        std::vector<QuantizedSubband> quantized_subbands;
        std::vector<float> reconstructed;
        size_t payload_bits = 0;
        float mse = 0.0f;
    };

    float scale_factor(uint8_t index);
    uint8_t optimal_sf_index_for_peak(float peak, uint8_t selector);
    QuantizedSubband choose_subband_encoding(const float* coefficients, size_t len, CodingMode coding_mode, const SearchOptions& options);
    SpectrumEncoding build_spectral_unit(const float* coefficients, size_t len, CodingMode coding_mode, const SearchOptions& options);
    ChannelSoundUnit build_basic_sound_unit_from_encoding(const SpectrumEncoding& encoding);
    uint8_t fast_peak_to_sf_index(const float* coefficients, size_t len);
    std::vector<uint8_t> compute_group_sf_indices(const float* coefficients, size_t len);
    int32_t estimate_band_bit_cost(const uint8_t* group_sf_indices, size_t group_count, uint8_t tbl_index, int32_t sf_index, size_t subband_width);
    void encode_mantissas_to_writer(class BitWriter& writer, uint8_t selector, CodingMode coding_mode, const int8_t* mantissas, size_t count);
    size_t coded_qmf_bands_for_subband_count(size_t subband_count);
    size_t fixed_sound_unit_bits(size_t subband_count);
    size_t candidate_total_bits(const QuantizedSubband* candidate);
    TonalExtractionResult extract_tonal_components(float* residual, size_t residual_len, size_t budget_bits, uint8_t coded_qmf_bands, CodingMode coding_mode, size_t max_entries);

}