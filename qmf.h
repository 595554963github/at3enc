#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include "common.h"
#include "gain.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace atrac3 {

    constexpr std::array<float, 24> EXE_QMF_48TAP_HALF = {
        0.000015258789f, 0.00009608304f, 0.00005861498f, -0.00031435178f,
        -0.00025285266f, 0.00089026295f, 0.00054333656f, -0.002123024f,
        -0.00081761716f, 0.004399848f, 0.00078923843f, -0.008183379f,
        0.00006384667f, 0.014029815f, -0.002570447f, -0.022687245f,
        0.008143067f, 0.03558198f, -0.019632578f, -0.056703273f,
        0.045504123f, 0.10373335f, -0.13785878f, -0.48455644f
    };

    constexpr std::array<float, 24> LEGACY_QMF_48TAP_HALF = {
        -0.00001461907f, -0.00009205479f, -0.00005615757f, 0.00030117269f,
        0.0002422519f, -0.00085293897f, -0.0005205574f, 0.0020340169f,
        0.0007833389f, -0.0042153862f, -0.0007561499f, 0.007840294f,
        -0.00006116992f, -0.01344162f, 0.002462682f, 0.021736089f,
        -0.007801671f, -0.03409022f, 0.01880949f, 0.05432601f,
        -0.04359638f, -0.09938437f, 0.1320791f, 0.4642416f
    };

    inline std::array<float, 48> mirrored_qmf_window() {
        std::array<float, 48> out{};
        for (size_t i = 0; i < 24; ++i) {
            float v = EXE_QMF_48TAP_HALF[i] * 2.0f;
            out[i] = v;
            out[47 - i] = v;
        }
        return out;
    }

    constexpr size_t DIRECT_QMF_STATE_SAMPLES = 138;
    constexpr size_t DIRECT_QMF_STAGE2_HISTORY_SAMPLES = 132;
    constexpr size_t DIRECT_QMF_SCRATCH_SAMPLES = 1024 + DIRECT_QMF_STATE_SAMPLES;
    constexpr float DEFAULT_DIRECT_QMF_OUTPUT_GAIN = 4.0f;

    constexpr std::array<size_t, 24> STAGE1_EVEN_COEFF_ORDER = {
        46, 2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 44, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40
    };

    constexpr std::array<size_t, 24> STAGE1_ODD_COEFF_ORDER = {
        47, 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 45, 1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41
    };

    constexpr std::array<int, 24> STAGE1_PHASE0_OFFSETS = {
        6, -38, -34, -30, -26, -22, -18, -14, -10, -6, -2, 2, 4, -40, -36, -32, -28, -24, -20, -16, -12, -8, -4, 0
    };

    constexpr std::array<int, 24> STAGE1_PHASE1_OFFSETS = {
        7, -37, -33, -29, -25, -21, -17, -13, -9, -5, -1, 3, 5, -39, -35, -31, -27, -23, -19, -15, -11, -7, -3, 1
    };

    constexpr std::array<int, 24> STAGE1_PHASE2_OFFSETS = {
        8, -36, -32, -28, -24, -20, -16, -12, -8, -4, 0, 4, 6, -38, -34, -30, -26, -22, -18, -14, -10, -6, -2, 2
    };

    constexpr std::array<int, 24> STAGE1_PHASE3_OFFSETS = {
        9, -35, -31, -27, -23, -19, -15, -11, -7, -3, 1, 5, 7, -37, -33, -29, -25, -21, -17, -13, -9, -5, -1, 3
    };

    constexpr std::array<float, 4> STAGE2_SIGN_MASKS = { -0.0f, 0.0f, -0.0f, 0.0f };

    inline std::array<float, 48> exe_qmf_window() {
        std::array<float, 48> out{};
        for (size_t i = 0; i < 24; ++i) {
            out[i] = EXE_QMF_48TAP_HALF[i];
            out[47 - i] = EXE_QMF_48TAP_HALF[i];
        }
        return out;
    }

    inline std::array<std::array<float, 24>, 4> build_direct_stage2_coeffs() {
        std::array<std::array<float, 24>, 4> out{};
        auto full = exe_qmf_window();
        for (int lane = 0; lane < 4; ++lane) {
            int parity = lane & 1;
            out[lane][0] = full[46 + parity];
            for (int tap = 0; tap < 23; ++tap) {
                out[lane][tap + 1] = full[parity + tap * 2];
            }
        }
        return out;
    }

    struct FourBandFrame {
        std::array<std::vector<float>, 4> bands;
        std::array<float, 1024> interleaved;
    };

    inline std::array<std::array<float, 32>, 4> estimate_envelopes_from_interleaved(
        const std::array<float, 1024>& interleaved) {
        std::array<std::array<float, 32>, 4> out{};
        size_t slots_per_band = 256 / 8;

        for (size_t slot_index = 0; slot_index < slots_per_band; ++slot_index) {
            size_t slot_base = slot_index * 8 * 4;
            for (size_t sample_offset = 0; sample_offset < 8; ++sample_offset) {
                size_t sample_base = slot_base + sample_offset * 4;
                for (size_t band_index = 0; band_index < 4; ++band_index) {
                    float abs_value = std::abs(interleaved[sample_base + band_index]);
                    if (abs_value > out[band_index][slot_index]) {
                        out[band_index][slot_index] = abs_value;
                    }
                }
            }
        }

        return out;
    }

    class TwoBandQmf {
    public:
        TwoBandQmf() : delay_{}, window_{} {
            window_ = exe_qmf_window();
        }

        std::pair<std::vector<float>, std::vector<float>> split_block(const std::vector<float>& input) {
            if ((input.size() & 1) != 0) {
                throw std::invalid_argument("QMF input must contain even number of samples");
            }

            std::vector<float> history;
            history.reserve(delay_.size() + input.size());
            history.insert(history.end(), delay_.begin(), delay_.end());
            history.insert(history.end(), input.begin(), input.end());

            size_t out_len = input.size() / 2;
            std::vector<float> low(out_len);
            std::vector<float> high(out_len);

            for (size_t k = 0; k < out_len; ++k) {
                size_t base = k * 2;
                float even_acc = 0.0f;
                float odd_acc = 0.0f;

                for (size_t tap = 0; tap < 24; ++tap) {
                    even_acc += history[base + tap * 2] * window_[tap * 2];
                    odd_acc += history[base + tap * 2 + 1] * window_[tap * 2 + 1];
                }

                low[k] = even_acc + odd_acc;
                high[k] = odd_acc - even_acc;
            }

            size_t delay_len = delay_.size();
            for (size_t i = 0; i < delay_len; ++i) {
                delay_[i] = history[history.size() - delay_len + i];
            }

            return { std::move(low), std::move(high) };
        }

    private:
        std::array<float, 46> delay_;
        std::array<float, 48> window_;
    };

    class FourBandQmf {
    public:
        FourBandQmf() : root_(), low_split_(), high_split_(), direct_state_{} {
            window_ = exe_qmf_window();
            direct_stage2_coeffs_ = build_direct_stage2_coeffs();
        }

        std::array<std::vector<float>, 4> split_frame(const std::vector<float>& frame) {
            return split_frame_with_layout(frame).bands;
        }

        FourBandFrame split_frame_with_layout(const std::vector<float>& frame) {
            if (frame.size() != 1024) {
                throw std::invalid_argument("ATRAC3 4-band analysis expects 1024 samples");
            }

            return split_frame_direct(frame);
        }

    private:
        FourBandFrame split_frame_direct(const std::vector<float>& frame) {
            std::array<float, DIRECT_QMF_SCRATCH_SAMPLES> scratch{};
            for (size_t i = 0; i < DIRECT_QMF_STATE_SAMPLES; ++i) {
                scratch[i] = direct_state_[i];
            }

            std::array<float, 1024> interleaved{};
            size_t cursor = DIRECT_QMF_STAGE2_HISTORY_SAMPLES;

            for (size_t sample_idx = 0; sample_idx < 256; ++sample_idx) {
                size_t input_base = sample_idx * 4;
                scratch[cursor + 6] = frame[input_base];
                scratch[cursor + 7] = frame[input_base + 1];
                scratch[cursor + 8] = frame[input_base + 2];
                scratch[cursor + 9] = frame[input_base + 3];

                float phase0 = direct_dot(scratch, cursor, STAGE1_EVEN_COEFF_ORDER, STAGE1_PHASE0_OFFSETS);
                float phase1 = direct_dot(scratch, cursor, STAGE1_ODD_COEFF_ORDER, STAGE1_PHASE1_OFFSETS);
                float phase2 = direct_dot(scratch, cursor, STAGE1_EVEN_COEFF_ORDER, STAGE1_PHASE2_OFFSETS);
                float phase3 = direct_dot(scratch, cursor, STAGE1_ODD_COEFF_ORDER, STAGE1_PHASE3_OFFSETS);

                float stage1_a = phase0 + phase1;
                float stage1_b = phase2 + phase3;
                float stage1_c = phase1 + (-phase0);
                float stage1_d = phase3 + (-phase2);

                scratch[cursor - 40] = stage1_a;
                scratch[cursor - 39] = stage1_b;
                scratch[cursor - 38] = stage1_c;
                scratch[cursor - 37] = stage1_d;

                float band_mix0 = direct_stage2_dot(scratch, cursor, stage1_a, 0);
                float band_mix1 = direct_stage2_dot(scratch, cursor, stage1_b, 1);
                float band_mix2 = direct_stage2_dot(scratch, cursor, stage1_c, 2);
                float band_mix3 = direct_stage2_dot(scratch, cursor, stage1_d, 3);

                interleaved[input_base] = (band_mix0 + apply_sign_mask(band_mix1, STAGE2_SIGN_MASKS[1])) * DEFAULT_DIRECT_QMF_OUTPUT_GAIN;
                interleaved[input_base + 1] = (band_mix1 + apply_sign_mask(band_mix0, STAGE2_SIGN_MASKS[0])) * DEFAULT_DIRECT_QMF_OUTPUT_GAIN;
                interleaved[input_base + 2] = (band_mix3 + apply_sign_mask(band_mix2, STAGE2_SIGN_MASKS[2])) * DEFAULT_DIRECT_QMF_OUTPUT_GAIN;
                interleaved[input_base + 3] = (band_mix2 + apply_sign_mask(band_mix3, STAGE2_SIGN_MASKS[3])) * DEFAULT_DIRECT_QMF_OUTPUT_GAIN;

                cursor += 4;
            }

            for (size_t i = 0; i < DIRECT_QMF_STATE_SAMPLES; ++i) {
                direct_state_[i] = scratch[1024 + i];
            }

            std::array<std::vector<float>, 4> bands;
            for (size_t i = 0; i < 4; ++i) {
                bands[i].resize(256);
            }
            for (size_t sample_idx = 0; sample_idx < 256; ++sample_idx) {
                size_t base = sample_idx * 4;
                bands[0][sample_idx] = interleaved[base];
                bands[1][sample_idx] = interleaved[base + 1];
                bands[2][sample_idx] = interleaved[base + 2];
                bands[3][sample_idx] = interleaved[base + 3];
            }

            return { bands, interleaved };
        }

        float direct_dot(const std::array<float, DIRECT_QMF_SCRATCH_SAMPLES>& scratch, size_t cursor,
            const std::array<size_t, 24>& coeff_order,
            const std::array<int, 24>& offsets) const {
            float sum = 0.0f;
            for (size_t i = 0; i < 24; ++i) {
                int offset = static_cast<int>(cursor) + offsets[i];
                if (offset >= 0 && static_cast<size_t>(offset) < scratch.size()) {
                    sum += window_[coeff_order[i]] * scratch[static_cast<size_t>(offset)];
                }
            }
            return sum;
        }

        float direct_stage2_dot(const std::array<float, DIRECT_QMF_SCRATCH_SAMPLES>& scratch, size_t cursor,
            float current, int lane) const {
            float acc = direct_stage2_coeffs_[lane][0] * current;
            for (int tap = 0; tap < 23; ++tap) {
                int offset = static_cast<int>(cursor) - static_cast<int>(DIRECT_QMF_STAGE2_HISTORY_SAMPLES) + lane + tap * 4;
                if (offset >= 0 && static_cast<size_t>(offset) < scratch.size()) {
                    acc += direct_stage2_coeffs_[lane][tap + 1] * scratch[static_cast<size_t>(offset)];
                }
            }
            return acc;
        }

        float apply_sign_mask(float value, float sign) const {
            union FloatBits {
                float f;
                uint32_t u;
            };
            FloatBits vb;
            FloatBits sb;
            vb.f = value;
            sb.f = sign;
            FloatBits result;
            result.u = vb.u ^ sb.u;
            return result.f;
        }

        TwoBandQmf root_;
        TwoBandQmf low_split_;
        TwoBandQmf high_split_;
        std::array<float, DIRECT_QMF_STATE_SAMPLES> direct_state_;
        std::array<float, 48> window_;
        std::array<std::array<float, 24>, 4> direct_stage2_coeffs_;
    };

}