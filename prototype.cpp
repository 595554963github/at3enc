#include "prototype.h"
#include "wav.h"
#include "bitstream.h"
#include "sound_unit.h"
#include "quant.h"
#include "gain.h"
#include "container.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cstring>

namespace atrac3 {

static float g_quantizer_compat_gain = 7500.0f;
static bool g_apply_odd_band_reverse = true;
static bool g_apply_gain_estimation = false;
static int g_analysis_sample_offset = 69;
static bool g_swap_gain_curve_order = false;
static bool g_use_reference_mdct = false;

static std::vector<float> slice_frame(const std::vector<float>& channel, size_t start) {
    std::vector<float> frame(1024, 0.0f);
    int shifted_start = static_cast<int>(start) + g_analysis_sample_offset;
    for (size_t dst_index = 0; dst_index < frame.size(); ++dst_index) {
        int src_index = shifted_start + static_cast<int>(dst_index);
        if (src_index >= 0 && src_index < static_cast<int>(channel.size())) {
            frame[dst_index] = channel[src_index];
        }
    }
    return frame;
}

static std::array<float, MDCT_COEFFS_PER_BAND> compensate_band_samples(
    const std::array<float, MDCT_COEFFS_PER_BAND>& band_samples,
    const std::array<float, GAIN_CURVE_SAMPLES>& curve) {
    std::array<float, MDCT_COEFFS_PER_BAND> out;
    for (size_t index = 0; index < MDCT_COEFFS_PER_BAND; ++index) {
        float gain = curve[index];
        if (std::abs(gain) > 1e-9f) {
            out[index] = band_samples[index] / gain;
        } else {
            out[index] = band_samples[index];
        }
    }
    return out;
}

static size_t minimum_subband_count_for_target_bits(size_t target_bits) {
    if (target_bits >= 1536) {
        for (size_t i = 0; i < 33; ++i) {
            if (ATRAC3_SUBBAND_TAB[i] >= 768) {
                return i;
            }
        }
    }
    return 0;
}

static void pad_spectral_unit(SpectrumEncoding& encoding, size_t min_subband_count) {
    while (encoding.spectral_unit.subbands.size() < min_subband_count) {
        SpectralSubband subband;
        subband.table_index = 0;
        subband.scale_factor_index = nullptr;
        encoding.spectral_unit.subbands.push_back(subband);
    }
}

static PrototypeFrameChannel encode_channel(ChannelSoundUnit sound_unit, SpectrumEncoding spectrum) {
    PrototypeFrameChannel result;
    result.sound_unit = sound_unit;
    result.spectrum = spectrum;
    BitWriter writer;
    sound_unit.write_to(writer);
    result.bit_len = writer.bits_written();
    writer.byte_align_zero();
    result.bytes = writer.into_bytes();
    return result;
}

PrototypeEncoder::PrototypeEncoder(size_t channel_count)
    : qmf(channel_count), mdct(),
    overlap(channel_count), previous_gain_bands(channel_count),
    previous_envelopes(channel_count), previous_peak_state(channel_count), pending_analysis(channel_count) {
    for (size_t i = 0; i < channel_count; ++i) {
        overlap[i] = std::array<std::array<float, MDCT_COEFFS_PER_BAND>, 4>();
        for (int j = 0; j < 4; ++j) {
            overlap[i][j].fill(0.0f);
        }
        previous_gain_bands[i].resize(4);
        previous_envelopes[i] = std::array<std::array<float, GAIN_HISTORY_SLOTS>, 4>();
        for (int j = 0; j < 4; ++j) {
            previous_envelopes[i][j].fill(0.0f);
        }
        previous_peak_state[i].fill(0.0f);
        pending_analysis[i].coefficients.resize(1024, 0.0f);
        pending_analysis[i].gain_bands.resize(4);
    }
}

PrototypeEncodeResult PrototypeEncoder::encode_wav(const float* samples, size_t sample_count, uint32_t sample_rate, size_t channels, const PrototypeOptions& options) {
    size_t channel_count = channels;

    std::vector<std::vector<float>> channel_samples(channels);
    for (size_t i = 0; i < channel_count; ++i) {
        channel_samples[i].resize(sample_count / channel_count);
        for (size_t j = 0; j < channel_samples[i].size(); ++j) {
            channel_samples[i][j] = samples[j * channel_count + i];
        }
    }

    size_t pcm_frames = sample_count / channels;
    size_t start_sample = options.start_frame * 1024;
    size_t remaining_pcm_frames = (pcm_frames > start_sample) ? (pcm_frames - start_sample) : 0;
    size_t input_frames = (remaining_pcm_frames + 1023) / 1024;
    size_t available_frames = input_frames + options.flush_frames;
    size_t frame_count = options.frame_limit.value_or(available_frames);
    if (frame_count > available_frames) {
        frame_count = available_frames;
    }

    std::vector<std::vector<AnalyzedChannel>> all_analyses;
    all_analyses.reserve(frame_count);

    for (size_t frame_index = 0; frame_index < frame_count; ++frame_index) {
        size_t start = start_sample + frame_index * 1024;
        std::vector<AnalyzedChannel> analyses;
        analyses.reserve(channel_count);
        for (size_t ch = 0; ch < channel_count; ++ch) {
            std::vector<float> frame_samples = slice_frame(channel_samples[ch], start);
            analyses.push_back(analyze_channel_for_encoding(ch, frame_samples));
        }
        pending_analysis = analyses;
        all_analyses.push_back(analyses);
    }

    SearchOptions search_opts;
    search_opts.lambda = options.lambda;
    search_opts.target_bits = options.target_bits_per_channel;
    search_opts.max_candidates_per_band = 64;
    for (int i = 0; i < 32; ++i) {
        search_opts.tonal_marked_subbands[i] = false;
    }

    std::vector<PrototypeFrame> frames;
    frames.reserve(frame_count);
    for (const auto& analysis : all_analyses) {
        frames.push_back(encode_analyzed_frame(analysis, options.coding_mode, search_opts));
    }

    std::vector<uint8_t> all_bytes;
    for (const auto& f : frames) {
        all_bytes.insert(all_bytes.end(), f.bytes.begin(), f.bytes.end());
    }

    std::vector<float> original_samples(samples, samples + sample_count);

    PrototypeEncodeResult result;
    result.sample_rate = sample_rate;
    result.channel_count = channel_count;
    result.frame_count = frame_count;
    result.frames = frames;
    result.bytes = all_bytes;
    result.original_samples = original_samples;
    result.options = options;
    return result;
}

PrototypeFrame PrototypeEncoder::encode_frame(const std::vector<const float*>& channels, CodingMode coding_mode, const SearchOptions& search) {
    if (channels.size() != qmf.size()) {
        throw std::runtime_error("channel count mismatch");
    }

    std::vector<AnalyzedChannel> current_analysis;
    current_analysis.reserve(channels.size());
    for (size_t channel_index = 0; channel_index < channels.size(); ++channel_index) {
        std::vector<float> samples_vec(1024);
        for (size_t i = 0; i < 1024; ++i) {
            samples_vec[i] = channels[channel_index][i];
        }
        current_analysis.push_back(analyze_channel_for_encoding(channel_index, samples_vec));
    }
    pending_analysis = current_analysis;
    return encode_analyzed_frame(current_analysis, coding_mode, search);
}

PrototypeFrame PrototypeEncoder::encode_analyzed_frame(const std::vector<AnalyzedChannel>& analysis_channels, CodingMode coding_mode, const SearchOptions& search) {
    if (analysis_channels.size() != qmf.size()) {
        throw std::runtime_error("analysis channel count mismatch");
    }

    std::vector<PrototypeFrameChannel> frame_channels;
    frame_channels.reserve(analysis_channels.size());
    BitWriter frame_writer;

    bool joint_stereo = (analysis_channels.size() == 2) && 
                        search.target_bits.has_value() && 
                        search.target_bits.value() < 1000;

    for (size_t ch = 0; ch < analysis_channels.size(); ++ch) {
        const auto& analysis = analysis_channels[ch];
        size_t gain_payload_bits = 0;
        for (const auto& band : analysis.gain_bands) {
            gain_payload_bits += band.points.size() * 9;
        }
        size_t min_subband_count = minimum_subband_count_for_target_bits(search.target_bits.value_or(0));
        size_t padded_skip_bits = (min_subband_count > 0) ? (min_subband_count - 1) * 3 : 0;
        size_t base_target = search.target_bits.value_or(1536);
        size_t spectral_budget = (base_target > gain_payload_bits + padded_skip_bits) ?
            (base_target - gain_payload_bits - padded_skip_bits) : 0;

        uint8_t coded_qmf_target = 3;
        std::vector<float> residual = analysis.coefficients;

        if (joint_stereo && ch == 1) {
            const auto& left_analysis = analysis_channels[0];
            for (size_t i = 0; i < residual.size() && i < left_analysis.coefficients.size(); ++i) {
                residual[i] = analysis.coefficients[i] - left_analysis.coefficients[i];
            }
        }

        TonalExtractionResult tonal_result = extract_tonal_components(
            residual.data(), residual.size(), spectral_budget, coded_qmf_target, coding_mode, 4
        );

        SearchOptions adjusted_search = search;
        adjusted_search.target_bits = (spectral_budget > tonal_result.tonal_bits) ?
            std::optional<size_t>(spectral_budget - tonal_result.tonal_bits) : std::optional<size_t>(0);
        for (int i = 0; i < 32; ++i) {
            adjusted_search.tonal_marked_subbands[i] = tonal_result.tonal_subbands[i];
        }

        SpectrumEncoding spectrum = build_spectral_unit(residual.data(), residual.size(), coding_mode, adjusted_search);
        pad_spectral_unit(spectrum, min_subband_count);

        ChannelSoundUnit sound_unit = build_basic_sound_unit_from_encoding(spectrum);

        size_t coded_qmf_bands = std::min(static_cast<size_t>(sound_unit.coded_qmf_bands), static_cast<size_t>(coded_qmf_target));
        sound_unit.coded_qmf_bands = static_cast<uint8_t>(coded_qmf_bands);
        sound_unit.gain_bands.resize(coded_qmf_bands);
        for (size_t i = 0; i < coded_qmf_bands; ++i) {
            sound_unit.gain_bands[i] = analysis.gain_bands[i];
        }

        sound_unit.tonal_mode_selector = tonal_result.tonal_mode_selector;
        sound_unit.tonal_components.reserve(tonal_result.tonal_components.size());
        for (const auto& tc : tonal_result.tonal_components) {
            TonalComponent c = tc;
            c.band_flags.resize(coded_qmf_bands, false);
            c.cells.resize(coded_qmf_bands * 4, TonalCell());
            sound_unit.tonal_components.push_back(c);
        }

        try {
            size_t total = sound_unit.bit_len();
            if (total > base_target) {
                sound_unit.tonal_mode_selector = TonalCodingModeSelector::AllVlc;
                sound_unit.tonal_components.clear();
                adjusted_search.target_bits = spectral_budget;
                SpectrumEncoding s2 = build_spectral_unit(
                    analysis.coefficients.data(),
                    analysis.coefficients.size(),
                    coding_mode,
                    adjusted_search
                );
                pad_spectral_unit(s2, min_subband_count);
                sound_unit.spectrum = s2.spectral_unit;
                spectrum = s2;
            }
        } catch (...) {
        }

        frame_channels.push_back(encode_channel(sound_unit, spectrum));
        frame_channels.back().sound_unit.write_to(frame_writer);
    }

    PrototypeFrame frame;
    frame.channels = std::move(frame_channels);
    frame.bit_len = frame_writer.bits_written();
    frame_writer.byte_align_zero();
    frame.bytes = frame_writer.into_bytes();
    return frame;
}

AnalyzedChannel PrototypeEncoder::analyze_channel_for_encoding(size_t channel_index, const std::vector<float>& samples) {
    return analyze_channel(channel_index, samples, true);
}

std::vector<float> PrototypeEncoder::analyze_channel_raw(size_t channel_index, const std::vector<float>& samples) {
    return analyze_channel(channel_index, samples, false).coefficients;
}

AnalyzedChannel PrototypeEncoder::analyze_channel(size_t channel_index, const std::vector<float>& samples, bool encode_gain) {
    if (samples.size() != 1024) {
        throw std::runtime_error("expected 1024 samples");
    }

    AnalyzedChannel result;
    result.coefficients.resize(1024, 0.0f);
    result.gain_bands.resize(4);

    FourBandFrame frame = qmf[channel_index].split_frame_with_layout(samples);
    auto envelopes = estimate_envelopes_from_interleaved(frame.interleaved);

    bool gain_enabled = encode_gain && g_apply_gain_estimation;

    for (size_t band_index = 0; band_index < 4; ++band_index) {
        std::array<float, GAIN_HISTORY_SLOTS> envelope;
        for (size_t i = 0; i < GAIN_HISTORY_SLOTS; ++i) {
            envelope[i] = envelopes[band_index][i];
        }

        std::array<float, GAIN_HISTORY_SLOTS> previous_envelope;
        for (size_t i = 0; i < GAIN_HISTORY_SLOTS; ++i) {
            previous_envelope[i] = previous_envelopes[channel_index][band_index][i];
        }

        previous_envelopes[channel_index][band_index] = envelope;

        float history_peak_state = previous_peak_state[channel_index][band_index];

        GainBand gain_band;
        if (gain_enabled) {
            gain_band = estimate_gain_band(envelope, previous_envelope, band_index, history_peak_state);
        }
        result.gain_bands[band_index] = gain_band;

        std::array<float, MDCT_COEFFS_PER_BAND> band_samples;
        for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
            band_samples[i] = frame.bands[band_index][i];
        }

        std::array<float, MDCT_COEFFS_PER_BAND> analysis_samples = band_samples;

        if (gain_enabled) {
            GainCurve curve;
            if (g_swap_gain_curve_order) {
                curve = build_gain_curve(previous_gain_bands[channel_index][band_index], gain_band);
            } else {
                curve = build_gain_curve(gain_band, previous_gain_bands[channel_index][band_index]);
            }
            analysis_samples = compensate_band_samples(band_samples, curve.samples);
        }

        std::array<float, MDCT_INPUT_SAMPLES> mdct_input;
        for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
            mdct_input[i] = overlap[channel_index][band_index][i];
        }
        for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
            mdct_input[MDCT_COEFFS_PER_BAND + i] = analysis_samples[i];
        }

        for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
            overlap[channel_index][band_index][i] = analysis_samples[i];
        }

        std::array<float, MDCT_COEFFS_PER_BAND> band_coefficients;
        if (g_use_reference_mdct) {
            band_coefficients = mdct.forward_reference(mdct_input);
        } else {
            band_coefficients = mdct.forward(mdct_input);
            for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
                band_coefficients[i] *= g_quantizer_compat_gain;
            }
        }

        if (g_apply_odd_band_reverse && (band_index & 1)) {
            std::reverse(band_coefficients.begin(), band_coefficients.end());
        }

        size_t start = band_index * MDCT_COEFFS_PER_BAND;
        for (size_t i = 0; i < MDCT_COEFFS_PER_BAND; ++i) {
            result.coefficients[start + i] = band_coefficients[i];
        }
    }

    if (gain_enabled) {
        for (size_t band_index = 0; band_index < 4; ++band_index) {
            float max_val = 0.0f;
            for (size_t i = 0; i < GAIN_HISTORY_SLOTS; ++i) {
                max_val = std::max(max_val, previous_envelopes[channel_index][band_index][i]);
            }
            previous_peak_state[channel_index][band_index] = max_val;
        }
        previous_gain_bands[channel_index] = result.gain_bands;
    }

    return result;
}

void PrototypeEncoder::reset() {
    for (size_t i = 0; i < qmf.size(); ++i) {
        qmf[i] = FourBandQmf();
    }
    for (size_t i = 0; i < overlap.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            overlap[i][j].fill(0.0f);
            previous_envelopes[i][j].fill(0.0f);
            previous_peak_state[i][j] = 0.0f;
            previous_gain_bands[i][j] = GainBand();
            pending_analysis[i].coefficients.assign(1024, 0.0f);
            pending_analysis[i].gain_bands.assign(4, GainBand());
        }
    }
}

void PrototypeEncoder::set_quantizer_compat_gain(float value) {
    g_quantizer_compat_gain = value;
}

void PrototypeEncoder::set_apply_odd_band_reverse(bool value) {
    g_apply_odd_band_reverse = value;
}

void PrototypeEncoder::set_apply_gain_estimation(bool value) {
    g_apply_gain_estimation = value;
}

void PrototypeEncoder::set_analysis_sample_offset(int value) {
    g_analysis_sample_offset = value;
}

void PrototypeEncoder::set_swap_gain_curve_order(bool value) {
    g_swap_gain_curve_order = value;
}

void PrototypeEncoder::set_use_reference_mdct(bool value) {
    g_use_reference_mdct = value;
}

std::vector<GainInspectionChannel> PrototypeEncoder::inspect_gain_frame(const std::vector<std::vector<float>>& channels) {
    if (channels.size() != qmf.size()) {
        throw std::runtime_error("channel count mismatch");
    }

    std::vector<GainInspectionChannel> debug_channels;
    debug_channels.reserve(channels.size());

    for (size_t channel_index = 0; channel_index < channels.size(); ++channel_index) {
        const auto& samples = channels[channel_index];
        if (samples.size() != 1024) {
            throw std::runtime_error("expected 1024 samples");
        }

        FourBandFrame frame = qmf[channel_index].split_frame_with_layout(samples);
        auto envelopes = estimate_envelopes_from_interleaved(frame.interleaved);
        std::array<std::array<float, GAIN_HISTORY_SLOTS>, 4> current_envelopes;
        std::vector<GainBand> current_gain_bands(4);
        std::vector<GainInspectionBand> debug_bands;
        debug_bands.reserve(4);

        for (size_t band_index = 0; band_index < 4; ++band_index) {
            auto previous_envelope = previous_envelopes[channel_index][band_index];
            auto current_envelope = envelopes[band_index];
            float history_peak_state = previous_peak_state[channel_index][band_index];
            GainBand gain_band;
            if (g_apply_gain_estimation) {
                gain_band = estimate_gain_band(current_envelope, previous_envelope, band_index, history_peak_state);
            }

            current_envelopes[band_index] = current_envelope;
            current_gain_bands[band_index] = gain_band;

            float max_val = 0.0f;
            for (size_t i = 0; i < GAIN_HISTORY_SLOTS; ++i) {
                max_val = std::max(max_val, previous_envelope[i]);
            }
            previous_peak_state[channel_index][band_index] = max_val;

            GainInspectionBand debug_band;
            debug_band.current_envelope = current_envelope;
            debug_band.previous_envelope = previous_envelope;
            debug_band.gain_band = gain_band;
            debug_bands.push_back(debug_band);
        }

        previous_envelopes[channel_index] = current_envelopes;
        previous_gain_bands[channel_index] = current_gain_bands;

        GainInspectionChannel debug_channel;
        debug_channel.bands = std::move(debug_bands);
        debug_channels.push_back(debug_channel);
    }

    return debug_channels;
}

std::vector<std::vector<float>> PrototypeEncoder::analyze_frame_coefficients(const std::vector<std::vector<float>>& channels) {
    if (channels.size() != qmf.size()) {
        throw std::runtime_error("channel count mismatch");
    }

    std::vector<std::vector<float>> result;
    result.reserve(channels.size());

    for (size_t channel_index = 0; channel_index < channels.size(); ++channel_index) {
        const auto& samples = channels[channel_index];
        result.push_back(analyze_channel_raw(channel_index, samples));
    }

    return result;
}

}
