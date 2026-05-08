#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <optional>
#include <algorithm>
#include <functional>
#include "common.h"
#include "bitstream.h"
#include "sound_unit.h"
#include "mdct.h"
#include "qmf.h"
#include "quant.h"
#include "gain.h"

namespace atrac3 {

    struct PrototypeOptions {
        CodingMode coding_mode;
        float lambda;
        std::optional<size_t> frame_limit;
        size_t start_frame;
        size_t flush_frames;
        std::optional<size_t> target_bits_per_channel;
        bool joint_stereo = false;

        PrototypeOptions() : coding_mode(CodingMode::Clc), lambda(0.0001f), start_frame(0), flush_frames(0), joint_stereo(false) {}
    };

    struct PrototypeFrameChannel {
        ChannelSoundUnit sound_unit;
        SpectrumEncoding spectrum;
        size_t bit_len;
        std::vector<uint8_t> bytes;
    };

    struct PrototypeFrame {
        std::vector<PrototypeFrameChannel> channels;
        size_t bit_len;
        std::vector<uint8_t> bytes;
    };

    struct PrototypeEncodeResult {
        uint32_t sample_rate;
        size_t channel_count;
        size_t frame_count;
        std::vector<PrototypeFrame> frames;
        std::vector<uint8_t> bytes;
        std::vector<float> original_samples;
        PrototypeOptions options;
    };

    struct GainInspectionBand {
        std::array<float, GAIN_HISTORY_SLOTS> current_envelope;
        std::array<float, GAIN_HISTORY_SLOTS> previous_envelope;
        GainBand gain_band;
    };

    struct GainInspectionChannel {
        std::vector<GainInspectionBand> bands;
    };

    struct AnalyzedChannel {
        std::vector<float> coefficients;
        std::vector<GainBand> gain_bands;
    };

    class PrototypeEncoder {
    public:
        PrototypeEncoder(size_t channel_count);

        PrototypeEncodeResult encode_wav(const float* samples, size_t sample_count, uint32_t sample_rate, size_t channels, const PrototypeOptions& options);

        PrototypeFrame encode_frame(const std::vector<const float*>& channels, CodingMode coding_mode, const SearchOptions& search);

        std::vector<GainInspectionChannel> inspect_gain_frame(const std::vector<std::vector<float>>& channels);

        std::vector<std::vector<float>> analyze_frame_coefficients(const std::vector<std::vector<float>>& channels);

        void reset();

        static void set_quantizer_compat_gain(float value);
        static void set_apply_odd_band_reverse(bool value);
        static void set_apply_gain_estimation(bool value);
        static void set_analysis_sample_offset(int value);
        static void set_swap_gain_curve_order(bool value);
        static void set_use_reference_mdct(bool value);

    private:
        AnalyzedChannel analyze_channel_for_encoding(size_t channel_index, const std::vector<float>& samples);

        std::vector<float> analyze_channel_raw(size_t channel_index, const std::vector<float>& samples);

        AnalyzedChannel analyze_channel(size_t channel_index, const std::vector<float>& samples, bool encode_gain);

        PrototypeFrame encode_analyzed_frame(const std::vector<AnalyzedChannel>& analysis_channels, CodingMode coding_mode, const SearchOptions& search);

        std::vector<FourBandQmf> qmf;
        Mdct256 mdct;
        std::vector<std::array<std::array<float, MDCT_COEFFS_PER_BAND>, 4>> overlap;
        std::vector<std::vector<GainBand>> previous_gain_bands;
        std::vector<std::array<std::array<float, GAIN_HISTORY_SLOTS>, 4>> previous_envelopes;
        std::vector<std::array<float, 4>> previous_peak_state;
        std::vector<AnalyzedChannel> pending_analysis;
    };

}