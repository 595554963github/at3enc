#include "synthesis.h"
#include "qmf.h"

namespace atrac3 {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

    static std::array<float, IMDCT_OUTPUT_SAMPLES> atrac3_decoder_window() {
        std::array<float, IMDCT_OUTPUT_SAMPLES> out = {};

        for (int i = 0; i < 128; ++i) {
            int j = 255 - i;
            double wi = std::sin(((static_cast<double>(i) + 0.5) / 256.0 - 0.5) * M_PI) + 1.0;
            double wj = std::sin(((static_cast<double>(j) + 0.5) / 256.0 - 0.5) * M_PI) + 1.0;
            double w = 0.5 * (wi * wi + wj * wj);

            out[i] = static_cast<float>(wi / w);
            out[511 - i] = out[i];
            out[j] = static_cast<float>(wj / w);
            out[511 - j] = out[j];
        }

        return out;
    }

    Imdct256::Imdct256() : window(atrac3_decoder_window()), scale(IMDCT_SCALE) {}

    std::array<float, IMDCT_OUTPUT_SAMPLES> Imdct256::inverse(const std::array<float, MDCT_COEFFS_PER_BAND>& input) const {
        double n = static_cast<double>(MDCT_COEFFS_PER_BAND);
        std::array<float, IMDCT_OUTPUT_SAMPLES> out = {};

        for (size_t sample_index = 0; sample_index < IMDCT_OUTPUT_SAMPLES; ++sample_index) {
            double acc = 0.0;
            double sample_idx_d = static_cast<double>(sample_index);
            for (size_t coeff_index = 0; coeff_index < MDCT_COEFFS_PER_BAND; ++coeff_index) {
                double angle = (M_PI / n) * (sample_idx_d + 0.5 + n / 2.0) * (static_cast<double>(coeff_index) + 0.5);
                acc += static_cast<double>(input[coeff_index]) * std::cos(angle);
            }
            out[sample_index] = static_cast<float>(acc * static_cast<double>(scale) * static_cast<double>(window[sample_index]));
        }

        return out;
    }

    IqmfState::IqmfState() : delay{}, window(mirrored_qmf_window()) {}

    std::vector<float> IqmfState::synthesize(const float* inlo, const float* inhi, size_t len) {
        if (len % 2 != 0) {
            throw std::runtime_error("iQMF input length must be even");
        }

        size_t n_in = len;
        std::vector<float> temp(46 + n_in * 2, 0.0f);
        for (size_t i = 0; i < 46; ++i) {
            temp[i] = delay[i];
        }

        float* p3 = temp.data() + 46;
        for (size_t i = 0; i < n_in; i += 2) {
            p3[2 * i] = inlo[i] + inhi[i];
            p3[2 * i + 1] = inlo[i] - inhi[i];
            if (i + 1 < n_in) {
                p3[2 * i + 2] = inlo[i + 1] + inhi[i + 1];
                p3[2 * i + 3] = inlo[i + 1] - inhi[i + 1];
            }
        }

        std::vector<float> out(n_in * 2, 0.0f);
        size_t start = 0;
        size_t out_index = 0;

        for (size_t i = 0; i < n_in; ++i) {
            float s1 = 0.0f;
            float s2 = 0.0f;

            for (size_t tap = 0; tap < 48; tap += 2) {
                s1 += temp[start + tap] * window[tap];
                s2 += temp[start + tap + 1] * window[tap + 1];
            }

            out[out_index] = s2;
            out[out_index + 1] = s1;
            start += 2;
            out_index += 2;
        }

        for (size_t i = 0; i < 46; ++i) {
            delay[i] = temp[n_in * 2 + i];
        }

        return out;
    }

    SynthesisChannel::SynthesisChannel() : imdct_overlap{}, qmf_low{}, qmf_high{}, qmf_root() {}

    Atrac3Synthesis::Atrac3Synthesis(size_t channel_count)
        : imdct(), channels(channel_count) {
    }

    std::vector<std::vector<float>> Atrac3Synthesis::synthesize_frame(const std::vector<const float*>& coefficients) {
        if (coefficients.size() != channels.size()) {
            throw std::runtime_error("channel count mismatch");
        }

        std::vector<std::vector<float>> result;
        result.reserve(channels.size());

        for (size_t channel_index = 0; channel_index < channels.size(); ++channel_index) {
            result.push_back(synthesize_channel(channel_index, coefficients[channel_index]));
        }

        return result;
    }

    std::vector<float> Atrac3Synthesis::synthesize_channel(size_t channel_index, const float* coefficients) {
        if (channel_index >= channels.size()) {
            throw std::runtime_error("channel index out of range");
        }

        SynthesisChannel& channel = channels[channel_index];
        std::array<std::array<float, SAMPLES_PER_BAND>, QMF_BANDS> bands{};

        for (size_t band_index = 0; band_index < QMF_BANDS; ++band_index) {
            size_t start = band_index * SAMPLES_PER_BAND;
            std::array<float, MDCT_COEFFS_PER_BAND> band_coefficients{};
            for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
                band_coefficients[i] = coefficients[start + i];
            }

            if (band_index & 1) {
                std::reverse(band_coefficients.begin(), band_coefficients.end());
            }

            std::array<float, IMDCT_OUTPUT_SAMPLES> imdct_output = imdct.inverse(band_coefficients);

            for (size_t sample_index = 0; sample_index < SAMPLES_PER_BAND; ++sample_index) {
                bands[band_index][sample_index] = imdct_output[sample_index] + channel.imdct_overlap[band_index][sample_index];
            }
            for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
                channel.imdct_overlap[band_index][i] = imdct_output[SAMPLES_PER_BAND + i];
            }
        }

        std::vector<float> low_half = channel.qmf_low.synthesize(bands[0].data(), bands[1].data(), SAMPLES_PER_BAND);
        std::vector<float> high_half = channel.qmf_high.synthesize(bands[3].data(), bands[2].data(), SAMPLES_PER_BAND);

        return channel.qmf_root.synthesize(low_half.data(), high_half.data(), low_half.size() / 2);
    }

} // namespace atrac3