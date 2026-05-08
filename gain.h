#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "common.h"
#include "sound_unit.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace atrac3 {

    constexpr size_t GAIN_HISTORY_SLOTS = 32;
    constexpr size_t GAIN_CURVE_SLOTS = 64;
    constexpr size_t GAIN_CURVE_SAMPLES = 256;
    constexpr size_t WINDOW_EDGE_ZERO_SAMPLES = 32;

    constexpr uint8_t UNITY_GAIN_LEVEL_CODE = 4;
    constexpr size_t GAIN_LEVEL_CODE_COUNT = 16;

    constexpr std::array<int32_t, GAIN_LEVEL_CODE_COUNT> GAIN_LEVEL_EXPONENTS = {
        -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    };

    constexpr std::array<float, 12> GAIN_INTERPOLATION_STEPS = {
        0.5946045f, 0.7070923f, 0.84088135f, 0.35354614f, 0.5f, 0.7070923f,
        0.2102356f, 0.35354614f, 0.5946045f, 0.125f, 0.25f, 0.5f
    };

    enum class DecoderWindowKind {
        Full,
        ZeroTail,
        ZeroHead,
        ZeroEdges
    };

    struct GainCurve {
        std::array<float, GAIN_CURVE_SAMPLES> samples;
        size_t first_change_sample;
    };

    inline std::array<float, GAIN_CURVE_SLOTS> combined_gain_history(
        const std::array<float, GAIN_HISTORY_SLOTS>& previous,
        const std::array<float, GAIN_HISTORY_SLOTS>& current)
    {
        std::array<float, GAIN_CURVE_SLOTS> out{};
        std::copy(previous.begin(), previous.end(), out.begin());
        std::copy(current.begin(), current.end(), out.begin() + GAIN_HISTORY_SLOTS);
        return out;
    }

    inline std::array<float, GAIN_CURVE_SLOTS> combined_gain_profile(
        const std::array<float, GAIN_HISTORY_SLOTS>& current,
        const std::array<float, GAIN_HISTORY_SLOTS>& previous)
    {
        std::array<float, GAIN_CURVE_SLOTS> out{};
        std::copy(current.begin(), current.end(), out.begin());
        std::copy(previous.begin(), previous.end(), out.begin() + GAIN_HISTORY_SLOTS);
        return out;
    }

    inline std::array<float, GAIN_HISTORY_SLOTS> estimate_envelope_slots(const std::vector<float>& samples)
    {
        if (samples.size() != SAMPLES_PER_BAND) {
            throw std::invalid_argument("gain envelope expects correct number of band samples");
        }

        std::array<float, GAIN_HISTORY_SLOTS> out{};
        size_t chunk_size = SAMPLES_PER_BAND / GAIN_HISTORY_SLOTS;
        for (size_t slot_index = 0; slot_index < GAIN_HISTORY_SLOTS; ++slot_index) {
            float max_val = 0.0f;
            for (size_t i = 0; i < chunk_size; ++i) {
                float val = std::abs(samples[slot_index * chunk_size + i]);
                if (val > max_val) {
                    max_val = val;
                }
            }
            out[slot_index] = max_val;
        }
        return out;
    }

    inline DecoderWindowKind decoder_window_kind(bool flag_a_active, bool flag_b_active)
    {
        if (!flag_a_active && !flag_b_active) {
            return DecoderWindowKind::Full;
        }
        else if (!flag_a_active && flag_b_active) {
            return DecoderWindowKind::ZeroTail;
        }
        else if (flag_a_active && !flag_b_active) {
            return DecoderWindowKind::ZeroHead;
        }
        else {
            return DecoderWindowKind::ZeroEdges;
        }
    }

    inline std::array<float, GAIN_CURVE_SAMPLES> decoder_window_table(DecoderWindowKind kind)
    {
        std::array<float, GAIN_CURVE_SAMPLES> out{};
        for (size_t index = 0; index < GAIN_CURVE_SAMPLES; ++index) {
            out[index] = static_cast<float>(
                std::sin(((static_cast<double>(index) + 0.5) * M_PI) / static_cast<double>(GAIN_CURVE_SAMPLES))
                );
        }

        switch (kind) {
        case DecoderWindowKind::Full:
            break;
        case DecoderWindowKind::ZeroTail:
            std::fill(out.begin() + (GAIN_CURVE_SAMPLES - WINDOW_EDGE_ZERO_SAMPLES), out.end(), 0.0f);
            break;
        case DecoderWindowKind::ZeroHead:
            std::fill(out.begin(), out.begin() + WINDOW_EDGE_ZERO_SAMPLES, 0.0f);
            break;
        case DecoderWindowKind::ZeroEdges:
            std::fill(out.begin(), out.begin() + WINDOW_EDGE_ZERO_SAMPLES, 0.0f);
            std::fill(out.begin() + (GAIN_CURVE_SAMPLES - WINDOW_EDGE_ZERO_SAMPLES), out.end(), 0.0f);
            break;
        }
        return out;
    }

    GainCurve build_gain_curve(const GainBand& current, const GainBand& previous);

    GainBand estimate_gain_band(
        const std::array<float, GAIN_HISTORY_SLOTS>& current,
        const std::array<float, GAIN_HISTORY_SLOTS>& previous,
        size_t band_index,
        float history_peak_state);

}