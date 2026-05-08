#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace atrac3 {

    constexpr size_t MDCT_COEFFS_PER_BAND = 256;
    constexpr size_t MDCT_INPUT_SAMPLES = MDCT_COEFFS_PER_BAND * 2;

    inline std::array<float, MDCT_COEFFS_PER_BAND> atrac3_analysis_window_half() {
        std::array<float, MDCT_COEFFS_PER_BAND> out{};
        for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
            double phase = ((static_cast<double>(i) + 0.5) / static_cast<double>(MDCT_COEFFS_PER_BAND) - 0.5) * M_PI;
            out[i] = static_cast<float>((std::sin(phase) + 1.0) * 0.5);
        }
        return out;
    }

    inline std::array<float, MDCT_INPUT_SAMPLES> symmetric_window_from_half(const std::array<float, MDCT_COEFFS_PER_BAND>& half) {
        std::array<float, MDCT_INPUT_SAMPLES> out{};
        for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
            out[i] = half[i];
            out[MDCT_COEFFS_PER_BAND + i] = half[MDCT_COEFFS_PER_BAND - 1 - i];
        }
        return out;
    }

    class Mdct256 {
    public:
        Mdct256() : Mdct256(symmetric_window_from_half(atrac3_analysis_window_half())) {}

        explicit Mdct256(const std::array<float, MDCT_INPUT_SAMPLES>& window)
            : analysis_window_(window)
        {
            build_pre_rotation();
            build_fft();
            build_post_rotation();
            build_bit_reverse();
        }

        std::array<float, MDCT_COEFFS_PER_BAND> forward(const std::array<float, MDCT_INPUT_SAMPLES>& input) const {
            std::array<float, 513> scratch{};
            std::array<bool, 128> visited{};
            visited.fill(false);

            for (int index = 0; index < 64; ++index) {
                scratch[257 + index] = -(analysis_window_[384 + index * 2] * input[384 + index * 2])
                    - (analysis_window_[383 - index * 2] * input[383 - index * 2]);
            }
            for (int index = 0; index < 128; ++index) {
                scratch[321 + index] = analysis_window_[index * 2] * input[index * 2]
                    - analysis_window_[255 - index * 2] * input[255 - index * 2];
            }
            for (int index = 0; index < 64; ++index) {
                scratch[449 + index] = analysis_window_[511 - index * 2] * input[511 - index * 2]
                    + analysis_window_[256 + index * 2] * input[256 + index * 2];
            }

            for (int index = 0; index < 128; ++index) {
                int source = 257 + index * 2;
                float left = scratch[source];
                float right = scratch[source + 1];
                float cos_val = pre_rotation_cos_[index];
                float neg_sin = pre_rotation_neg_sin_[index];

                scratch[1 + index] = cos_val * left - neg_sin * right;
                scratch[129 + index] = cos_val * right + neg_sin * left;
            }

            for (int index = 0; index < 128; ++index) {
                if (visited[index]) continue;
                size_t target = bit_reverse_[index];
                float real_val = scratch[1 + target];
                float imag_val = scratch[129 + target];
                scratch[1 + target] = scratch[1 + index];
                scratch[129 + target] = scratch[129 + index];
                scratch[1 + index] = real_val;
                scratch[129 + index] = imag_val;
                visited[target] = true;
            }

            size_t span = 64;
            size_t stage_width = 2;
            for (int stage = 0; stage < 7; ++stage) {
                size_t half_width = stage_width / 2;
                size_t block_start = 0;
                size_t block_end = half_width;

                for (size_t block = 0; block < span; ++block) {
                    size_t twiddle = 0;
                    size_t upper = block_start;
                    size_t lower = block_end;

                    for (size_t i = 0; i < half_width; ++i) {
                        float rotated_real = fft_cos_[twiddle] * scratch[1 + lower]
                            - fft_neg_sin_[twiddle] * scratch[129 + lower];
                        float rotated_imag = fft_cos_[twiddle] * scratch[129 + lower]
                            + fft_neg_sin_[twiddle] * scratch[1 + lower];

                        float upper_real = scratch[1 + upper];
                        float upper_imag = scratch[129 + upper];
                        scratch[1 + lower] = upper_real - rotated_real;
                        scratch[129 + lower] = upper_imag - rotated_imag;
                        scratch[1 + upper] = upper_real + rotated_real;
                        scratch[129 + upper] = upper_imag + rotated_imag;

                        twiddle += span;
                        upper += 1;
                        lower += 1;
                    }

                    block_start += stage_width;
                    block_end += stage_width;
                }

                span /= 2;
                stage_width *= 2;
            }

            for (int index = 0; index < 128; ++index) {
                float left_real = scratch[256 - index];
                float right_real = scratch[1 + index];
                float right_imag = scratch[129 + index];
                float left_imag = scratch[128 - index];

                scratch[257 + index] = post_alpha_minus_[index] * left_real
                    + post_beta_minus_[index] * right_real
                    + post_alpha_plus_[index] * right_imag
                    + post_beta_plus_[index] * left_imag;
                scratch[512 - index] = left_real * post_beta_plus_[index]
                    + (post_alpha_plus_[index] * right_real - left_imag * post_alpha_minus_[index])
                    - post_beta_minus_[index] * right_imag;
            }

            std::array<float, MDCT_COEFFS_PER_BAND> output{};
            for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
                output[i] = scratch[257 + i] * (1.0f / 128.0f);
            }
            return output;
        }

        std::array<float, MDCT_COEFFS_PER_BAND> forward_reference(const std::array<float, MDCT_INPUT_SAMPLES>& input) const {
            double n = static_cast<double>(MDCT_COEFFS_PER_BAND);
            std::array<float, MDCT_COEFFS_PER_BAND> output{};

            for (int k = 0; k < MDCT_COEFFS_PER_BAND; ++k) {
                double acc = 0.0;
                for (size_t idx = 0; idx < MDCT_INPUT_SAMPLES; ++idx) {
                    double windowed = static_cast<double>(input[idx]) * static_cast<double>(analysis_window_[idx]);
                    double angle = (M_PI / n) * ((static_cast<double>(idx) + 0.5 + n / 2.0) * (static_cast<double>(k) + 0.5));
                    acc += windowed * std::cos(angle);
                }
                output[k] = static_cast<float>(acc);
            }
            return output;
        }

    private:
        void build_pre_rotation() {
            constexpr size_t PRE_ROTATION_BINS = 128;
            pre_rotation_cos_.resize(PRE_ROTATION_BINS);
            pre_rotation_neg_sin_.resize(PRE_ROTATION_BINS);
            for (size_t i = 0; i < PRE_ROTATION_BINS; ++i) {
                double phase = static_cast<double>(i) * M_PI / PRE_ROTATION_BINS;
                pre_rotation_cos_[i] = static_cast<float>(std::cos(phase));
                pre_rotation_neg_sin_[i] = static_cast<float>(-std::sin(phase));
            }
        }

        void build_fft() {
            constexpr size_t FFT_TWIDDLES = 64;
            fft_cos_.resize(FFT_TWIDDLES);
            fft_neg_sin_.resize(FFT_TWIDDLES);
            for (size_t i = 0; i < FFT_TWIDDLES; ++i) {
                double phase = static_cast<double>(i) * M_PI / FFT_TWIDDLES;
                fft_cos_[i] = static_cast<float>(std::cos(phase));
                fft_neg_sin_[i] = static_cast<float>(-std::sin(phase));
            }
        }

        void build_post_rotation() {
            constexpr size_t PRE_ROTATION_BINS = 128;
            post_alpha_minus_.resize(PRE_ROTATION_BINS);
            post_beta_minus_.resize(PRE_ROTATION_BINS);
            post_beta_plus_.resize(PRE_ROTATION_BINS);
            post_alpha_plus_.resize(PRE_ROTATION_BINS);

            for (size_t i = 0; i < PRE_ROTATION_BINS; ++i) {
                double alpha = static_cast<double>(2 * i + 1) * M_PI / 1024.0;
                double beta = 5.0 * alpha;

                post_alpha_minus_[i] = static_cast<float>(0.5 * (std::cos(beta) - std::sin(alpha)));
                post_alpha_plus_[i] = static_cast<float>(0.5 * (std::cos(beta) + std::sin(alpha)));
                post_beta_minus_[i] = static_cast<float>(0.5 * (std::cos(alpha) - std::sin(beta)));
                post_beta_plus_[i] = static_cast<float>(0.5 * (std::cos(alpha) + std::sin(beta)));
            }
        }

        void build_bit_reverse() {
            constexpr size_t PRE_ROTATION_BINS = 128;
            bit_reverse_.resize(PRE_ROTATION_BINS);
            for (size_t i = 0; i < PRE_ROTATION_BINS; ++i) {
                size_t reversed = 0;
                for (int bit = 0; bit < 7; ++bit) {
                    if (i & (1ULL << bit)) {
                        reversed |= (1ULL << (6 - bit));
                    }
                }
                bit_reverse_[i] = reversed;
            }
        }

        std::array<float, MDCT_INPUT_SAMPLES> analysis_window_;
        std::vector<float> pre_rotation_cos_;
        std::vector<float> pre_rotation_neg_sin_;
        std::vector<float> fft_cos_;
        std::vector<float> fft_neg_sin_;
        std::vector<float> post_alpha_minus_;
        std::vector<float> post_beta_minus_;
        std::vector<float> post_beta_plus_;
        std::vector<float> post_alpha_plus_;
        std::vector<size_t> bit_reverse_;
    };

}
