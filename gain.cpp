#include "gain.h"
#include <algorithm>
#include <cmath>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

namespace atrac3 {

    constexpr float GAIN_TRIGGER_RATIO = 1.85f;
    constexpr float GAIN_MIN_ABS_LEVEL = 1e-4f;

    static float gain_exponent_to_scale(int32_t exponent)
    {
        if (exponent < 0) {
            return 1.0f / static_cast<float>(1u << static_cast<uint32_t>(-exponent));
        }
        else {
            return static_cast<float>(1u << static_cast<uint32_t>(exponent));
        }
    }

    static size_t first_non_unity_sample(const std::array<float, GAIN_CURVE_SAMPLES>& samples)
    {
        for (size_t i = 0; i < GAIN_CURVE_SAMPLES; ++i) {
            if (std::abs(samples[i] - 1.0f) > 1e-6f) {
                return i;
            }
        }
        return GAIN_CURVE_SAMPLES - 1;
    }

    static size_t band_history_slots(size_t band_index)
    {
        return (8 - std::min(band_index, static_cast<size_t>(7))) * 8;
    }

    static std::array<float, 8> coarse_history_maxima(const std::array<float, GAIN_HISTORY_SLOTS>& previous)
    {
        std::array<float, 8> coarse{};
        for (size_t group = 0; group < 8; ++group) {
            float max_val = 0.0f;
            for (size_t i = 0; i < 4; ++i) {
                float val = previous[group * 4 + i];
                if (val > max_val) {
                    max_val = val;
                }
            }
            coarse[group] = max_val;
        }
        return coarse;
    }

    static float gain_threshold_multiplier(int32_t mode)
    {
        return (mode == -1) ? 2.0f : GAIN_TRIGGER_RATIO;
    }

    static float forward_threshold_multiplier(int32_t mode)
    {
        return (mode == -1) ? 2.0f : 1.6f;
    }

    static size_t coarse_location(
        const std::array<float, GAIN_CURVE_SLOTS>& history,
        size_t coarse_index,
        float threshold)
    {
        size_t base = coarse_index * 4;
        size_t location = base;
        if (coarse_index != 0 &&
            (base + 4 < history.size() ? history[base + 4] : 0.0f) < threshold &&
            history[base + 3] < threshold)
        {
            location = base + 3;
            if (history[base + 2] < threshold) {
                location = base + 2;
            }
        }
        return location;
    }

    static int32_t gain_step_from_ratio(float peak, float baseline)
    {
        float ratio = (peak / std::max(baseline, GAIN_MIN_ABS_LEVEL)) * static_cast<float>(M_SQRT2);
        union { float f; uint32_t u; } conv;
        conv.f = std::max(ratio, 1.0f);
        uint32_t bits = conv.u;
        return std::max(static_cast<int32_t>((bits >> 23) - 127), 0);
    }

    static size_t merge_backward_points(
        std::array<int32_t, 8>& positions,
        std::array<int32_t, 8>& level_deltas,
        size_t forward_count,
        size_t backward_start)
    {
        size_t write_index = forward_count;
        size_t read_index = backward_start;
        while (read_index < 7) {
            positions[write_index] = positions[read_index];
            level_deltas[write_index] = level_deltas[read_index];
            write_index++;
            read_index++;
        }
        return write_index;
    }

    GainCurve build_gain_curve(const GainBand& current, const GainBand& previous)
    {
        std::array<int32_t, GAIN_CURVE_SLOTS> slot_exponents{};

        size_t previous_slot = 0;
        for (const auto& point : previous.points) {
            int32_t exponent = GAIN_LEVEL_EXPONENTS[point.level];
            size_t end_slot = std::min(static_cast<size_t>(point.location) + GAIN_HISTORY_SLOTS, GAIN_CURVE_SLOTS - 1);
            while (previous_slot <= end_slot) {
                slot_exponents[previous_slot] = exponent;
                previous_slot++;
            }
            if (previous_slot == GAIN_CURVE_SLOTS) {
                break;
            }
        }

        size_t current_slot = 0;
        for (const auto& point : current.points) {
            int32_t exponent = GAIN_LEVEL_EXPONENTS[point.level];
            size_t end_slot = std::min(static_cast<size_t>(point.location), GAIN_CURVE_SLOTS - 1);
            while (current_slot <= end_slot) {
                slot_exponents[current_slot] += exponent;
                current_slot++;
            }
            if (current_slot == GAIN_CURVE_SLOTS) {
                break;
            }
        }

        std::array<float, GAIN_CURVE_SAMPLES> samples{};
        int32_t previous_exponent = slot_exponents[GAIN_CURVE_SLOTS - 1];
        float current_gain = gain_exponent_to_scale(previous_exponent);

        for (size_t slot_index = GAIN_CURVE_SLOTS - 1; slot_index < GAIN_CURVE_SLOTS; slot_index--) {
            int32_t exponent = slot_exponents[slot_index];
            size_t base = slot_index * 4;
            if (exponent == previous_exponent) {
                samples[base] = current_gain;
                samples[base + 1] = current_gain;
                samples[base + 2] = current_gain;
                samples[base + 3] = current_gain;
            }
            else if (exponent < previous_exponent) {
                float old_gain = gain_exponent_to_scale(previous_exponent);
                size_t interp_index = std::min(
                    static_cast<size_t>(previous_exponent - exponent - 1) * 3,
                    GAIN_INTERPOLATION_STEPS.size() - 3
                );
                samples[base + 1] = old_gain * GAIN_INTERPOLATION_STEPS[interp_index];
                samples[base + 2] = old_gain * GAIN_INTERPOLATION_STEPS[interp_index + 1];
                samples[base + 3] = old_gain * GAIN_INTERPOLATION_STEPS[interp_index + 2];
                current_gain = gain_exponent_to_scale(exponent);
                samples[base] = current_gain;
            }
            else {
                current_gain = gain_exponent_to_scale(exponent);
                size_t interp_index = std::min(
                    static_cast<size_t>(exponent - previous_exponent - 1) * 3,
                    GAIN_INTERPOLATION_STEPS.size() - 3
                );
                samples[base + 1] = current_gain * GAIN_INTERPOLATION_STEPS[interp_index + 2];
                samples[base + 2] = current_gain * GAIN_INTERPOLATION_STEPS[interp_index + 1];
                samples[base + 3] = current_gain * GAIN_INTERPOLATION_STEPS[interp_index];
                samples[base] = current_gain;
            }
            previous_exponent = exponent;
        }

        size_t first_change_sample = first_non_unity_sample(samples);

        GainCurve result;
        result.samples = samples;
        result.first_change_sample = std::min(first_change_sample, GAIN_CURVE_SAMPLES - 1);
        return result;
    }

    GainBand estimate_gain_band(
        const std::array<float, GAIN_HISTORY_SLOTS>& current,
        const std::array<float, GAIN_HISTORY_SLOTS>& previous,
        size_t band_index,
        float history_peak_state)
    {
        auto history = combined_gain_profile(current, previous);
        int32_t mode = static_cast<int32_t>(std::min(band_index, static_cast<size_t>(7)));
        auto coarse = coarse_history_maxima(current);

        float scan_max = coarse[7];
        size_t limit = band_history_slots(band_index);
        for (size_t i = GAIN_HISTORY_SLOTS; i < limit; ++i) {
            if (i >= history.size()) break;
            float val = history[i];
            if (val > scan_max) {
                scan_max = val;
            }
        }
        scan_max = std::max(scan_max, GAIN_MIN_ABS_LEVEL);

        float threshold_mul = gain_threshold_multiplier(mode);
        float threshold = scan_max * threshold_mul;
        int32_t remaining_steps = 4;
        std::array<int32_t, 8> positions;
        std::array<int32_t, 8> level_deltas;
        positions.fill(32);
        level_deltas.fill(0);
        size_t coarse_insert = 7;

        for (size_t coarse_index = 7; coarse_index < 8; coarse_index--) {
            float peak = coarse[coarse_index];
            if (scan_max <= peak) {
                if (peak > GAIN_MIN_ABS_LEVEL && peak > threshold) {
                    coarse_insert--;
                    positions[coarse_insert] = static_cast<int32_t>(coarse_location(history, coarse_index, threshold));
                    int32_t step = std::min(gain_step_from_ratio(peak, scan_max), remaining_steps);
                    remaining_steps -= step;
                    level_deltas[coarse_insert] = -step;

                    if (remaining_steps < 1 || coarse_insert == 5) {
                        break;
                    }
                }

                threshold = peak * threshold_mul;
                scan_max = peak;
            }
        }

        size_t backward_start = coarse_insert;
        size_t forward_count = 0;
        int32_t additional_steps = 15 - remaining_steps;
        if (additional_steps > 0) {
            float running_peak = std::max(history_peak_state, std::max(history[0], GAIN_MIN_ABS_LEVEL));
            float running_threshold = running_peak * forward_threshold_multiplier(mode);
            size_t scan_limit = std::min(static_cast<size_t>(positions[backward_start]), static_cast<size_t>(32));

            for (size_t slot = 0; slot < scan_limit; ++slot) {
                if (slot + 1 >= history.size()) break;
                float value = history[slot + 1];
                if (value < running_peak) {
                    continue;
                }
                if (value <= GAIN_MIN_ABS_LEVEL || value <= running_threshold) {
                    running_threshold = value * forward_threshold_multiplier(mode);
                    running_peak = value;
                    continue;
                }

                positions[forward_count] = static_cast<int32_t>(slot);
                int32_t step = gain_step_from_ratio(value, running_peak);
                if (forward_count > 0 &&
                    positions[forward_count - 1] == static_cast<int32_t>(slot) - 1 &&
                    level_deltas[forward_count - 1] <= step)
                {
                    forward_count--;
                    additional_steps += level_deltas[forward_count];
                    step += level_deltas[forward_count];
                }
                if (step > additional_steps) {
                    step = additional_steps;
                }
                additional_steps -= step;
                level_deltas[forward_count] = step;
                forward_count++;

                if (forward_count == backward_start || additional_steps < 1) {
                    break;
                }

                running_threshold = value * forward_threshold_multiplier(mode);
                running_peak = value;
            }
        }

        size_t total_points = merge_backward_points(positions, level_deltas, forward_count, backward_start);
        if (total_points == 0) {
            return GainBand();
        }

        int32_t running_level = static_cast<int32_t>(UNITY_GAIN_LEVEL_CODE);
        for (size_t index = total_points - 1; index < total_points; index--) {
            running_level += level_deltas[index];
            level_deltas[index] = running_level;
        }

        GainBand band;
        for (size_t index = 0; index < total_points; ++index) {
            int32_t level_clamped = level_deltas[index];
            if (level_clamped < 0) level_clamped = 0;
            if (level_clamped > static_cast<int32_t>(GAIN_LEVEL_CODE_COUNT - 1)) level_clamped = GAIN_LEVEL_CODE_COUNT - 1;
            uint8_t level = static_cast<uint8_t>(level_clamped);

            int32_t pos_clamped = positions[index];
            if (pos_clamped < 0) pos_clamped = 0;
            if (pos_clamped > 31) pos_clamped = 31;
            uint8_t location = static_cast<uint8_t>(pos_clamped);

            band.points.push_back(GainPoint{ level, location });
        }

        return band;
    }

} // namespace atrac3