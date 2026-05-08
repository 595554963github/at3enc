#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "common.h"
#include "mdct.h"
#include "qmf.h"

namespace atrac3 {

    constexpr size_t IMDCT_OUTPUT_SAMPLES = 512;
    constexpr float IMDCT_SCALE = 2.0f / static_cast<float>(MDCT_COEFFS_PER_BAND);

    struct Imdct256 {
        std::array<float, IMDCT_OUTPUT_SAMPLES> window;
        float scale;

        Imdct256();

        std::array<float, IMDCT_OUTPUT_SAMPLES> inverse(const std::array<float, MDCT_COEFFS_PER_BAND>& input) const;
    };

    struct IqmfState {
        std::array<float, 46> delay;
        std::array<float, 48> window;

        IqmfState();

        std::vector<float> synthesize(const float* inlo, const float* inhi, size_t len);
    };

    struct SynthesisChannel {
        std::array<std::array<float, MDCT_COEFFS_PER_BAND>, QMF_BANDS> imdct_overlap;
        IqmfState qmf_low;
        IqmfState qmf_high;
        IqmfState qmf_root;

        SynthesisChannel();
    };

    struct Atrac3Synthesis {
        Imdct256 imdct;
        std::vector<SynthesisChannel> channels;

        Atrac3Synthesis(size_t channel_count);

        std::vector<std::vector<float>> synthesize_frame(const std::vector<const float*>& coefficients);

    private:
        std::vector<float> synthesize_channel(size_t channel_index, const float* coefficients);
    };

}