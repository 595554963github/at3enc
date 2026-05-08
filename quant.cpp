#include "quant.h"
#include "bitstream.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>
#include <map>
#include <stdexcept>

namespace atrac3 {

struct HuffmanEntry {
    int8_t symbol;
    uint32_t code;
    uint8_t bits;
};

static std::vector<std::vector<HuffmanEntry>> get_huffman_codebooks() {
    static std::vector<std::vector<HuffmanEntry>> codebooks;
    if (codebooks.empty()) {
        size_t offset = 0;
        for (size_t i = 0; i < ATRAC3_HUFF_TAB_SIZES.size(); ++i) {
            size_t size = ATRAC3_HUFF_TAB_SIZES[i];
            size_t next = offset + size;
            
            auto build_canonical_codebook = [](const std::vector<std::pair<uint8_t, uint8_t>>& raw) {
                std::vector<HuffmanEntry> result;
                if (raw.empty()) return result;
                
                size_t max_bits = 0;
                for (const auto& entry : raw) {
                    max_bits = std::max(max_bits, (size_t)entry.second);
                }
                
                std::vector<uint32_t> counts(max_bits + 1, 0);
                for (const auto& entry : raw) {
                    counts[entry.second]++;
                }
                
                std::vector<uint32_t> next_codes(max_bits + 1, 0);
                uint32_t code = 0;
                for (size_t bits = 1; bits <= max_bits; ++bits) {
                    code = (code + counts[bits - 1]) << 1;
                    next_codes[bits] = code;
                }
                
                for (const auto& entry : raw) {
                    uint8_t bits = entry.second;
                    uint32_t code = next_codes[bits];
                    next_codes[bits]++;
                    result.push_back({static_cast<int8_t>(entry.first - 31), code, bits});
                }
                
                return result;
            };
            
            std::vector<std::pair<uint8_t, uint8_t>> raw;
            for (size_t j = offset; j < next; ++j) {
                raw.push_back(ATRAC3_HUFF_TABS[j]);
            }
            codebooks.push_back(build_canonical_codebook(raw));
            offset = next;
        }
    }
    return codebooks;
}

std::vector<float> QuantizedSubband::dequantized(size_t len) const {
    std::vector<float> res(len, 0.0f);
    if (table_index == 0 || !has_scale_factor) {
        return res;
    }
    float sf = scale_factor(scale_factor_index);
    float inv_max = ATRAC3_INV_MAX_QUANT[table_index];
    float scale = sf * inv_max;
    size_t n = std::min(mantissas.size(), len);
    for (size_t i = 0; i < n; ++i) {
        res[i] = mantissas[i] * scale;
    }
    return res;
}

QuantizedSubband QuantizedSubband::uncoded(const float* coeffs, size_t len) {
    QuantizedSubband q;
    q.table_index = 0;
    q.has_scale_factor = false;
    q.scale_factor_index = 0;
    q.payload_bits = 0;
    q.mse = 0.0f;
    for (size_t i = 0; i < len; ++i) {
        float f = coeffs ? coeffs[i] : 0.0f;
        q.mse += f * f;
    }
    if (len > 0) {
        q.mse /= (float)len;
    }
    return q;
}

static const std::array<std::array<int8_t, 2>, 9> ATRAC3_MANTISSA_VLC_PAIRS = {{
    {0, 0}, {0, 1}, {0, -1}, {1, 0}, {-1, 0}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}
}};

float scale_factor(uint8_t index) {
    return std::pow(2.0f, (index - 15.0f) / 3.0f);
}

static uint32_t twos_complement_bits(int32_t value, uint8_t bits) {
    if (bits == 32) {
        return static_cast<uint32_t>(value);
    } else {
        uint32_t mask = (1u << bits) - 1;
        return static_cast<uint32_t>(value) & mask;
    }
}

static uint8_t clc_symbol_index(int8_t value) {
    for (size_t i = 0; i < ATRAC3_MANTISSA_CLC_TAB.size(); ++i) {
        if (ATRAC3_MANTISSA_CLC_TAB[i] == value) {
            return static_cast<uint8_t>(i);
        }
    }
    throw std::runtime_error("mantissa not found in CLC table");
}

static int8_t vlc_pair_symbol(int8_t left, int8_t right) {
    for (int i = 0; i < 9; ++i) {
        if (ATRAC3_MANTISSA_VLC_PAIRS[i][0] == left && ATRAC3_MANTISSA_VLC_PAIRS[i][1] == right) {
            return static_cast<int8_t>(i);
        }
    }
    throw std::runtime_error("pair not found in VLC pairs");
}

static const HuffmanEntry* find_huffman_entry(uint8_t selector, int8_t symbol) {
    auto codebooks = get_huffman_codebooks();
    if (selector < 1 || selector > 7) return nullptr;
    const auto& codebook = codebooks[selector - 1];
    for (const auto& entry : codebook) {
        if (entry.symbol == symbol) {
            return &entry;
        }
    }
    return nullptr;
}

static void encode_mantissas_to_payload(uint8_t selector, CodingMode coding_mode, const std::vector<int8_t>& mantissas, RawBitPayload& payload) {
    if (selector < 1 || selector > 7) {
        return;
    }
    
    if (coding_mode == CodingMode::Clc) {
        if (selector == 1) {
            for (size_t i = 0; i + 1 < mantissas.size(); i += 2) {
                uint8_t hi = clc_symbol_index(mantissas[i]);
                uint8_t lo = clc_symbol_index(mantissas[i + 1]);
                payload.push_bits((hi << 2) | lo, ATRAC3_CLC_LENGTH_TAB[1]);
            }
        } else {
            uint8_t width = ATRAC3_CLC_LENGTH_TAB[selector];
            for (int8_t m : mantissas) {
                uint32_t val = twos_complement_bits(m, width);
                payload.push_bits(val, width);
            }
        }
    } else {
        if (selector == 1) {
            for (size_t i = 0; i + 1 < mantissas.size(); i += 2) {
                int8_t symbol = vlc_pair_symbol(mantissas[i], mantissas[i + 1]);
                const HuffmanEntry* entry = find_huffman_entry(1, symbol);
                if (entry) {
                    payload.push_bits(entry->code, entry->bits);
                }
            }
        } else {
            for (int8_t m : mantissas) {
                const HuffmanEntry* entry = find_huffman_entry(selector, m);
                if (entry) {
                    payload.push_bits(entry->code, entry->bits);
                }
            }
        }
    }
}

static int32_t round_ties_even(float x) {
    float floor_x = std::floor(x);
    float frac = x - floor_x;
    if (frac < 0.5f) return (int32_t)floor_x;
    if (frac > 0.5f) return (int32_t)floor_x + 1;
    int32_t lower = (int32_t)floor_x;
    return (lower & 1) == 0 ? lower : lower + 1;
}

static std::vector<int8_t> quantize_selector1_clc(const float* coefficients, size_t len, float scale) {
    std::vector<int8_t> mantissas;
    mantissas.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        float val = coefficients[i] / scale;
        int8_t best = 0;
        float best_err = 1e30f;
        for (int8_t candidate : ATRAC3_MANTISSA_CLC_TAB) {
            float err = std::abs(val - candidate);
            if (err < best_err) {
                best_err = err;
                best = candidate;
            }
        }
        mantissas.push_back(best);
    }
    return mantissas;
}

static std::vector<int8_t> quantize_selector1_vlc(const float* coefficients, size_t len, float scale) {
    std::vector<int8_t> mantissas;
    mantissas.reserve(len);
    for (size_t i = 0; i + 1 < len; i += 2) {
        float left = coefficients[i] / scale;
        float right = coefficients[i + 1] / scale;
        
        std::array<int8_t, 2> best_pair = {0, 0};
        float best_err = 1e30f;
        
        for (const auto& pair : ATRAC3_MANTISSA_VLC_PAIRS) {
            float err_left = std::abs(left - pair[0]);
            float err_right = std::abs(right - pair[1]);
            float err = err_left * err_left + err_right * err_right;
            if (err < best_err) {
                best_err = err;
                best_pair = pair;
            }
        }
        
        mantissas.push_back(best_pair[0]);
        mantissas.push_back(best_pair[1]);
    }
    return mantissas;
}

static std::vector<int8_t> quantize_signed_clc(const float* coefficients, size_t len, uint8_t selector, float scale) {
    uint8_t width = ATRAC3_CLC_LENGTH_TAB[selector];
    int32_t min_value = -(1i32 << (width - 1));
    int32_t max_value = (1i32 << (width - 1)) - 1;
    
    std::vector<int8_t> mantissas;
    mantissas.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        int32_t val = round_ties_even(coefficients[i] / scale);
        val = std::max(min_value, std::min(max_value, val));
        mantissas.push_back(static_cast<int8_t>(val));
    }
    return mantissas;
}

static std::vector<int8_t> quantize_vlc(const float* coefficients, size_t len, uint8_t selector, float scale) {
    auto codebooks = get_huffman_codebooks();
    if (selector < 1 || selector > 7) {
        return {};
    }
    
    const auto& codebook = codebooks[selector - 1];
    std::vector<int8_t> mantissas;
    mantissas.reserve(len);
    
    for (size_t i = 0; i < len; ++i) {
        float val = coefficients[i] / scale;
        int8_t best = 0;
        float best_err = 1e30f;
        
        for (const auto& entry : codebook) {
            float err = std::abs(val - entry.symbol);
            if (err < best_err) {
                best_err = err;
                best = entry.symbol;
            }
        }
        
        mantissas.push_back(best);
    }
    
    return mantissas;
}

static QuantizedSubband quantize_subband_internal(const float* coefficients, size_t len, uint8_t selector, uint8_t sf_index, CodingMode coding_mode) {
    float sf = scale_factor(sf_index);
    float inv_max_q = ATRAC3_INV_MAX_QUANT[selector];
    float scale = sf * inv_max_q;
    
    std::vector<int8_t> mantissas;
    
    if (coding_mode == CodingMode::Clc) {
        if (selector == 1) {
            mantissas = quantize_selector1_clc(coefficients, len, scale);
        } else {
            mantissas = quantize_signed_clc(coefficients, len, selector, scale);
        }
    } else {
        if (selector == 1) {
            mantissas = quantize_selector1_vlc(coefficients, len, scale);
        } else {
            mantissas = quantize_vlc(coefficients, len, selector, scale);
        }
    }
    
    int final_sf = sf_index;
    std::vector<int8_t> final_mantissas = mantissas;
    float final_scale = scale;
    
    if (selector >= 2 && !final_mantissas.empty()) {
        while (final_sf <= 60) {
            bool any_odd = false;
            bool any_nonzero = false;
            for (int8_t m : final_mantissas) {
                if (m & 1) any_odd = true;
                if (m != 0) any_nonzero = true;
            }
            if (any_odd || !any_nonzero) break;
            for (int8_t& m : final_mantissas) m /= 2;
            final_sf += 3;
            final_scale = scale_factor(final_sf) * inv_max_q;
        }
    }
    
    size_t payload_bits = 0;
    {
        BitWriter w;
        if (coding_mode == CodingMode::Clc) {
            if (selector == 1) {
                for (size_t i = 0; i + 1 < final_mantissas.size(); i += 2) {
                    uint8_t hi = clc_symbol_index(final_mantissas[i]);
                    uint8_t lo = clc_symbol_index(final_mantissas[i + 1]);
                    w.write_bits((hi << 2) | lo, ATRAC3_CLC_LENGTH_TAB[1]);
                }
            } else {
                uint8_t width = ATRAC3_CLC_LENGTH_TAB[selector];
                for (int8_t m : final_mantissas) {
                    uint32_t val = twos_complement_bits(m, width);
                    w.write_bits(val, width);
                }
            }
        } else {
            if (selector == 1) {
                for (size_t i = 0; i + 1 < final_mantissas.size(); i += 2) {
                    int8_t symbol = vlc_pair_symbol(final_mantissas[i], final_mantissas[i + 1]);
                    const HuffmanEntry* entry = find_huffman_entry(1, symbol);
                    if (entry) {
                        w.write_bits(entry->code, entry->bits);
                    }
                }
            } else {
                for (int8_t m : final_mantissas) {
                    const HuffmanEntry* entry = find_huffman_entry(selector, m);
                    if (entry) {
                        w.write_bits(entry->code, entry->bits);
                    }
                }
            }
        }
        payload_bits = w.bits_written();
    }
    
    float mse = 0.0f;
    float max_abs_err = 0.0f;
    for (size_t i = 0; i < len; ++i) {
        float recon = final_mantissas[i] * final_scale;
        float err = coefficients[i] - recon;
        mse += err * err;
        max_abs_err = std::max(max_abs_err, std::abs(err));
    }
    mse /= len;
    
    QuantizedSubband q;
    q.table_index = selector;
    q.has_scale_factor = true;
    q.scale_factor_index = static_cast<uint8_t>(final_sf);
    q.mantissas = final_mantissas;
    q.payload_bits = payload_bits;
    q.mse = mse;
    q.max_abs_err = max_abs_err;
    return q;
}

uint8_t optimal_sf_index_for_peak(float peak, uint8_t selector) {
    if (peak <= 0.0f) return 0;
    float max_mantissa = 0.0f;
    switch (selector) {
        case 1: max_mantissa = 1.0f; break;
        case 2: max_mantissa = 2.0f; break;
        case 3: max_mantissa = 3.0f; break;
        case 4: max_mantissa = 4.0f; break;
        case 5: max_mantissa = 7.0f; break;
        case 6: max_mantissa = 15.0f; break;
        case 7: max_mantissa = 31.0f; break;
        default: return 0;
    }
    float imq = ATRAC3_INV_MAX_QUANT[selector];
    float needed_sf = peak / (max_mantissa * imq);
    float log_sf = std::log2(std::max(needed_sf, 1e-10f));
    int index = static_cast<int>(std::ceil(log_sf * 3.0f + 15.0f));
    if (index < 0) index = 0;
    if (index > 63) index = 63;
    return static_cast<uint8_t>(index);
}

QuantizedSubband choose_subband_encoding(const float* coefficients, size_t len, CodingMode coding_mode, const SearchOptions& options) {
    if (len == 0) return QuantizedSubband::uncoded(coefficients, len);

    QuantizedSubband best = QuantizedSubband::uncoded(coefficients, len);
    float best_score = best.mse;

    float peak = 0.0f;
    for (size_t i = 0; i < len; ++i) peak = std::max(peak, std::abs(coefficients[i]));
    
    for (uint8_t selector = 1; selector <= 7; ++selector) {
        uint8_t center = optimal_sf_index_for_peak(peak, selector);
        int lo = std::max(0, (int)center - 6);
        int hi = std::min(63, (int)center + 6);

        for (int sf = lo; sf <= hi; ++sf) {
            QuantizedSubband candidate = quantize_subband_internal(coefficients, len, selector, (uint8_t)sf, coding_mode);
            float score = candidate.mse + options.lambda * candidate.payload_bits;
            if (score < best_score - 1e-12f) {
                best = candidate;
                best_score = score;
            }
        }
    }
    return best;
}

struct BudgetSolution {
    std::vector<QuantizedSubband> selected;
    size_t used_bits;
    float mse;
};

static BudgetSolution solve_band_budget(const std::vector<std::vector<QuantizedSubband>>& candidates, size_t target_bits) {
    if (candidates.empty()) return { {}, 0, 1e30f };

    size_t coded_qmf_bands = coded_qmf_bands_for_subband_count(candidates.size());
    size_t fixed_bits = 6 + 2 + coded_qmf_bands * 3 + 5 + 5 + 1;

    if (fixed_bits > target_bits) return { {}, 0, 1e30f };

    size_t band_budget = target_bits - fixed_bits;
    size_t band_count = candidates.size();
    size_t state_count = band_budget + 1;

    std::vector<float> costs((band_count + 1) * state_count, 1e30f);
    std::vector<int> parents(band_count * state_count, -1);
    std::vector<int> choices(band_count * state_count, -1);
    costs[0] = 0.0f;

    for (size_t band = 0; band < band_count; ++band) {
        size_t current = band * state_count;
        size_t next = (band + 1) * state_count;
        for (size_t used = 0; used <= band_budget; ++used) {
            float cur_cost = costs[current + used];
            if (cur_cost >= 1e29f) continue;

            for (size_t ci = 0; ci < candidates[band].size(); ++ci) {
                const auto& cand = candidates[band][ci];
                size_t cand_bits;
                if (cand.table_index == 0) {
                    cand_bits = 3;
                } else {
                    cand_bits = 3 + 6 + cand.payload_bits;
                }
                size_t next_bits = used + cand_bits;
                if (next_bits > band_budget) continue;

                float next_cost = cur_cost + cand.mse;
                size_t slot = next + next_bits;
                if (next_cost < costs[slot] - 1e-12f) {
                    costs[slot] = next_cost;
                    parents[band * state_count + next_bits] = static_cast<int>(used);
                    choices[band * state_count + next_bits] = static_cast<int>(ci);
                }
            }
        }
    }

    size_t final_offset = band_count * state_count;
    size_t best_bits = band_budget;
    float best_cost = costs[final_offset + band_budget];

    for (size_t used = 0; used <= band_budget; ++used) {
        float cost = costs[final_offset + used];
        if (cost < best_cost - 1e-12f) {
            best_cost = cost;
            best_bits = used;
        }
    }

    if (best_cost >= 1e29f) return { {}, 0, 1e30f };

    std::vector<QuantizedSubband> selected;
    size_t used = best_bits;
    for (int band = static_cast<int>(band_count) - 1; band >= 0; --band) {
        int ci = choices[band * state_count + used];
        if (ci < 0) {
            for (size_t i = 0; i < band_count; ++i) {
                selected.push_back(QuantizedSubband::uncoded(nullptr, 0));
            }
            return { selected, fixed_bits, best_cost };
        }
        selected.insert(selected.begin(), candidates[band][ci]);
        used = parents[band * state_count + used];
    }

    size_t total_bits = fixed_bits;
    for (const auto& s : selected) {
        if (s.table_index == 0) {
            total_bits += 3;
        } else {
            total_bits += 3 + 6 + s.payload_bits;
        }
    }

    return { selected, total_bits, best_cost };
}

static SpectrumEncoding build_spectral_unit_budgeted(
    const float* coefficients,
    size_t len,
    CodingMode coding_mode,
    const SearchOptions& search,
    size_t target_bits) {
    if (len != 1024) throw std::runtime_error("expected 1024 coefficients");
    if (target_bits == 0) throw std::runtime_error("target_bits must be positive");

    std::vector<std::vector<uint8_t>> group_sf;
    group_sf.reserve(32);
    std::array<int32_t, 32> band_peak_sf = {};
    std::array<int32_t, 32> band_energy_sum = {};
    size_t num_active_bands = 0;

    for (int band = 0; band < 32; ++band) {
        size_t start = ATRAC3_SUBBAND_TAB[band];
        size_t end = ATRAC3_SUBBAND_TAB[band + 1];
        auto groups = compute_group_sf_indices(coefficients + start, end - start);

        int32_t peak = 0;
        int32_t energy = 0;
        for (uint8_t sf : groups) {
            int32_t s = static_cast<int32_t>(sf);
            peak = std::max(peak, s);
            energy += s;
        }

        if (energy < ENERGY_THRESHOLD[band] && peak < 3) {
            peak = 0;
            energy = 0;
        }

        band_peak_sf[band] = peak;
        band_energy_sum[band] = energy;
        if (peak > 0) {
            num_active_bands = band + 1;
        }
        group_sf.push_back(groups);
    }

    if (num_active_bands == 0) {
        num_active_bands = 1;
    }

    size_t fixed_overhead = fixed_sound_unit_bits(num_active_bands);
    int32_t available_bits = static_cast<int32_t>(target_bits) - static_cast<int32_t>(fixed_overhead);
    if (available_bits < 0) available_bits = 0;

    std::array<uint8_t, 32> tbl_indices = {};
    std::array<int32_t, 32> sf_indices = {};
    std::array<std::array<int32_t, 8>, 32> cost_at = {};

    for (size_t band = 0; band < num_active_bands; ++band) {
        if (band_peak_sf[band] == 0) continue;
        sf_indices[band] = band_peak_sf[band];
        size_t width = ATRAC3_SUBBAND_TAB[band + 1] - ATRAC3_SUBBAND_TAB[band];
        for (uint8_t tbl = 1; tbl <= 7; ++tbl) {
            cost_at[band][tbl] = estimate_band_bit_cost(
                group_sf[band].data(),
                group_sf[band].size(),
                tbl,
                sf_indices[band],
                width);
        }
    }

    int32_t budget_10x = available_bits * 10;

    std::array<int32_t, 32> effective_peak = band_peak_sf;
    for (size_t band = 0; band < num_active_bands; ++band) {
        if (!search.tonal_marked_subbands[band]) {
            continue;
        }
        int32_t tbl_guess = ((band_peak_sf[band] + 4) / 8);
        if (tbl_guess < 1) tbl_guess = 1;
        if (tbl_guess > 7) tbl_guess = 7;
        size_t spread_idx = ((tbl_guess - 1) >> 1);
        spread_idx = std::min(spread_idx, (size_t)4);
        int32_t factor = std::abs(TONAL_SPREAD_Q0[spread_idx]);
        int32_t lower_peak = (band > 0) ? band_peak_sf[band - 1] : 0;
        int32_t upper_peak = (band + 1 < num_active_bands) ? band_peak_sf[band + 1] : 0;

        int32_t self_delta = (factor * (band_peak_sf[band] + upper_peak)) >> 8;
        effective_peak[band] = std::min(effective_peak[band] + self_delta, (int32_t)63);
        if (band + 1 < num_active_bands) {
            int32_t up_delta = (factor * (band_peak_sf[band] + upper_peak)) >> 10;
            effective_peak[band + 1] = std::min(effective_peak[band + 1] + up_delta, (int32_t)63);
        }
        if (band > 0) {
            int32_t dn_delta = (factor * (band_peak_sf[band] + lower_peak)) >> 11;
            effective_peak[band - 1] = std::min(effective_peak[band - 1] + dn_delta, (int32_t)63);
        }
    }

    int32_t total_cost = 0;
    for (size_t band = 0; band < num_active_bands; ++band) {
        if (band_peak_sf[band] == 0) continue;
        int32_t initial = ((effective_peak[band] + 4) / 8);
        if (initial < 1) initial = 1;
        if (initial > 7) initial = 7;
        tbl_indices[band] = static_cast<uint8_t>(initial);
        total_cost += cost_at[band][initial];
    }

    while (total_cost > budget_10x) {
        int best_band = -1;
        int32_t best_savings = 0;
        float best_score = -1e30f;

        for (size_t band = 0; band < num_active_bands; ++band) {
            uint8_t current_tbl = tbl_indices[band];
            if (current_tbl <= 1) continue;
            uint8_t next_tbl = current_tbl - 1;
            int32_t savings = cost_at[band][current_tbl] - cost_at[band][next_tbl];
            if (savings <= 0) continue;

            float weight = (band < 18) ? (float)LF_OVERSHOOT_WEIGHTS_Q8[band] : (float)LF_OVERSHOOT_WEIGHTS_Q8[17];
            float score = weight * (float)savings / ((float)band_peak_sf[band] + 1.0f);
            if (score > best_score) {
                best_score = score;
                best_band = static_cast<int>(band);
                best_savings = savings;
            }
        }

        if (best_band >= 0) {
            total_cost -= best_savings;
            tbl_indices[best_band] -= 1;
        } else {
            break;
        }
    }

    while (true) {
        int best_band = -1;
        float best_efficiency = -1e30f;
        int32_t best_delta = 0;

        for (size_t band = 0; band < num_active_bands; ++band) {
            uint8_t current_tbl = tbl_indices[band];
            if (current_tbl == 0 || current_tbl >= 7) continue;
            uint8_t next_tbl = current_tbl + 1;
            int32_t delta = cost_at[band][next_tbl] - cost_at[band][current_tbl];
            if (delta <= 0 || total_cost + delta > budget_10x) continue;
            float efficiency = (float)effective_peak[band] / (float)delta;
            if (efficiency > best_efficiency) {
                best_efficiency = efficiency;
                best_band = static_cast<int>(band);
                best_delta = delta;
            }
        }

        if (best_band >= 0) {
            tbl_indices[best_band] += 1;
            total_cost += best_delta;
        } else {
            break;
        }
    }

    std::vector<QuantizedSubband> quantized_subbands;
    quantized_subbands.reserve(num_active_bands);
    std::vector<float> reconstructed(len, 0.0f);
    std::vector<SpectralSubband> spectral_subbands;
    spectral_subbands.reserve(num_active_bands);
    size_t payload_bits = 0;
    size_t used_bits = fixed_overhead;

    for (size_t band = 0; band < num_active_bands; ++band) {
        size_t start = ATRAC3_SUBBAND_TAB[band];
        size_t end = ATRAC3_SUBBAND_TAB[band + 1];
        const float* slice = coefficients + start;
        size_t band_len = end - start;

        if (tbl_indices[band] == 0 || band_peak_sf[band] == 0) {
            QuantizedSubband uncoded = QuantizedSubband::uncoded(slice, band_len);
            SpectralSubband subband;
            subband.table_index = 0;
            subband.scale_factor_index = nullptr;
            spectral_subbands.push_back(subband);
            quantized_subbands.push_back(uncoded);
            used_bits += 3;
            continue;
        }

        uint8_t selector = tbl_indices[band];
        float peak = 0.0f;
        for (size_t i = 0; i < band_len; ++i) {
            peak = std::max(peak, std::abs(slice[i]));
        }
        uint8_t sf_center = optimal_sf_index_for_peak(peak, selector);
        QuantizedSubband best;
        float best_score = 1e30f;
        bool has_best = false;

        for (int8_t delta = 0; delta <= 3; ++delta) {
            int32_t sf_val = static_cast<int32_t>(sf_center) + delta;
            if (sf_val < 0) sf_val = 0;
            if (sf_val > 63) sf_val = 63;
            uint8_t sf_try = static_cast<uint8_t>(sf_val);
            try {
                QuantizedSubband candidate = quantize_subband_internal(slice, band_len, selector, sf_try, coding_mode);
                size_t bits = candidate_total_bits(&candidate);
                if (used_bits + bits <= target_bits) {
                    float score = candidate.max_abs_err;
                    if (!has_best || score < best_score - 1e-12f) {
                        best = candidate;
                        best_score = score;
                        has_best = true;
                    }
                }
            } catch (...) {}
        }

        if (!has_best) {
            for (uint8_t fallback_sel = selector - 1; fallback_sel >= 1; --fallback_sel) {
                uint8_t sf_fb = optimal_sf_index_for_peak(peak, fallback_sel);
                try {
                    QuantizedSubband candidate = quantize_subband_internal(slice, band_len, fallback_sel, sf_fb, coding_mode);
                    size_t bits = candidate_total_bits(&candidate);
                    if (used_bits + bits <= target_bits) {
                        best = candidate;
                        has_best = true;
                        break;
                    }
                } catch (...) {}
            }
        }

        if (!has_best) {
            best = QuantizedSubband::uncoded(slice, band_len);
        }

        used_bits += candidate_total_bits(&best);
        auto recon_band = best.dequantized(band_len);
        for (size_t i = 0; i < band_len; ++i) {
            reconstructed[start + i] = recon_band[i];
        }
        payload_bits += best.payload_bits;
        SpectralSubband subband;
        subband.table_index = best.table_index;
        if (best.has_scale_factor && best.table_index != 0) {
            subband.scale_factor_index = new uint8_t(best.scale_factor_index);
            encode_mantissas_to_payload(best.table_index, coding_mode, best.mantissas, subband.payload);
        } else {
            subband.scale_factor_index = nullptr;
        }
        spectral_subbands.push_back(subband);
        quantized_subbands.push_back(best);
    }

    size_t surplus = target_bits - used_bits;
    for (int pass = 0; pass < 4; ++pass) {
        if (surplus < 20) break;
        std::vector<size_t> order;
        for (size_t b = 0; b < quantized_subbands.size(); ++b) {
            uint8_t t = quantized_subbands[b].table_index;
            if (t > 0 && t < 7) {
                order.push_back(b);
            }
        }
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            uint8_t sa = quantized_subbands[a].has_scale_factor ? quantized_subbands[a].scale_factor_index : 0;
            uint8_t sb = quantized_subbands[b].has_scale_factor ? quantized_subbands[b].scale_factor_index : 0;
            return sb < sa;
        });

        size_t upgrades_this_pass = 0;
        for (size_t band : order) {
            if (surplus < 20) break;
            size_t start = ATRAC3_SUBBAND_TAB[band];
            size_t end = ATRAC3_SUBBAND_TAB[band + 1];
            const float* slice = coefficients + start;
            size_t band_len = end - start;
            uint8_t cur_tbl = quantized_subbands[band].table_index;
            if (cur_tbl >= 7) continue;
            uint8_t new_tbl = cur_tbl + 1;
            float peak = 0.0f;
            for (size_t i = 0; i < band_len; ++i) {
                peak = std::max(peak, std::abs(slice[i]));
            }
            uint8_t sf_center = optimal_sf_index_for_peak(peak, new_tbl);
            QuantizedSubband best_upgrade;
            float best_upgrade_score = quantized_subbands[band].mse;
            bool has_upgrade = false;

            for (int8_t delta = -2; delta <= 2; ++delta) {
                int32_t sf_val = static_cast<int32_t>(sf_center) + delta;
                if (sf_val < 0) sf_val = 0;
                if (sf_val > 63) sf_val = 63;
                uint8_t sf_try = static_cast<uint8_t>(sf_val);
                try {
                    QuantizedSubband cand = quantize_subband_internal(slice, band_len, new_tbl, sf_try, coding_mode);
                    size_t extra = candidate_total_bits(&cand) - candidate_total_bits(&quantized_subbands[band]);
                    if (extra <= surplus && cand.mse < best_upgrade_score - 1e-12f) {
                        best_upgrade = cand;
                        best_upgrade_score = cand.mse;
                        has_upgrade = true;
                    }
                } catch (...) {}
            }

            if (has_upgrade) {
                size_t extra = candidate_total_bits(&best_upgrade) - candidate_total_bits(&quantized_subbands[band]);
                surplus = surplus > extra ? surplus - extra : 0;
                used_bits += extra;
                auto recon = best_upgrade.dequantized(band_len);
                for (size_t i = 0; i < band_len; ++i) {
                    reconstructed[start + i] = recon[i];
                }
                payload_bits = payload_bits + best_upgrade.payload_bits - quantized_subbands[band].payload_bits;
                SpectralSubband subband;
                subband.table_index = best_upgrade.table_index;
                if (best_upgrade.has_scale_factor && best_upgrade.table_index != 0) {
                    subband.scale_factor_index = new uint8_t(best_upgrade.scale_factor_index);
                    encode_mantissas_to_payload(best_upgrade.table_index, coding_mode, best_upgrade.mantissas, subband.payload);
                } else {
                    subband.scale_factor_index = nullptr;
                }
                spectral_subbands[band] = subband;
                quantized_subbands[band] = best_upgrade;
                upgrades_this_pass += 1;
            }
        }
        if (upgrades_this_pass == 0) break;
    }

    while (spectral_subbands.size() > 1) {
        auto& last = spectral_subbands.back();
        if (last.table_index == 0) {
            spectral_subbands.pop_back();
        } else {
            break;
        }
    }

    for (size_t band = num_active_bands; band < 32; ++band) {
        size_t start = ATRAC3_SUBBAND_TAB[band];
        size_t end = ATRAC3_SUBBAND_TAB[band + 1];
        reconstructed[start] = 0.0f;
        for (size_t i = start + 1; i < end; ++i) reconstructed[i] = 0.0f;
        quantized_subbands.push_back(QuantizedSubband::uncoded(coefficients + start, end - start));
    }

    if (spectral_subbands.empty()) {
        spectral_subbands.push_back(SpectralSubband());
        spectral_subbands[0].table_index = 0;
        spectral_subbands[0].scale_factor_index = nullptr;
    }

    float mse = 0.0f;
    for (size_t i = 0; i < len; ++i) {
        float err = coefficients[i] - reconstructed[i];
        mse += err * err;
    }
    mse /= len;

    SpectrumEncoding result;
    result.spectral_unit.coding_mode = coding_mode;
    result.spectral_unit.subbands = spectral_subbands;
    result.quantized_subbands = quantized_subbands;
    result.reconstructed = reconstructed;
    result.payload_bits = payload_bits;
    result.mse = mse;
    return result;
}

SpectrumEncoding build_spectral_unit(const float* coefficients, size_t len, CodingMode coding_mode, const SearchOptions& options) {
    if (options.target_bits.has_value() && options.target_bits.value() > 0) {
        return build_spectral_unit_budgeted(coefficients, len, coding_mode, options, options.target_bits.value());
    }

    if (len != 1024) throw std::runtime_error("expected 1024 coefficients");

    std::vector<std::vector<QuantizedSubband>> candidates(32);

    float peak = 0.0f;

    for (int band = 0; band < 32; ++band) {
        size_t start = ATRAC3_SUBBAND_TAB[band];
        size_t end = ATRAC3_SUBBAND_TAB[band + 1];
        size_t band_len = end - start;

        if (band_len == 0) continue;

        candidates[band].push_back(QuantizedSubband::uncoded(coefficients + start, band_len));

        for (size_t i = 0; i < band_len; ++i) {
            peak = std::max(peak, std::abs(coefficients[start + i]));
        }

        for (uint8_t selector = 1; selector <= 7; ++selector) {
            uint8_t center = optimal_sf_index_for_peak(peak, selector);
            int lo = std::max(0, (int)center - 6);
            int hi = std::min(63, (int)center + 6);
            for (int sf = lo; sf <= hi; ++sf) {
                auto cand = quantize_subband_internal(coefficients + start, band_len, selector, (uint8_t)sf, coding_mode);
                bool duplicate = false;
                for (const auto& existing : candidates[band]) {
                    size_t existing_bits = (existing.table_index == 0) ? 3 : (3 + 6 + existing.payload_bits);
                    size_t cand_bits = (cand.table_index == 0) ? 3 : (3 + 6 + cand.payload_bits);
                    if (existing_bits == cand_bits && existing.mse <= cand.mse + 1e-12f) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    candidates[band].push_back(cand);
                }
            }
        }
        if (candidates[band].size() > 64) {
            std::sort(candidates[band].begin(), candidates[band].end(),
                [](const QuantizedSubband& a, const QuantizedSubband& b) {
                    size_t bits_a = (a.table_index == 0) ? 3 : (3 + 6 + a.payload_bits);
                    size_t bits_b = (b.table_index == 0) ? 3 : (3 + 6 + b.payload_bits);
                    if (bits_a != bits_b) return bits_a < bits_b;
                    return a.mse < b.mse;
                });
            candidates[band].resize(64);
        }
    }

    auto solution = solve_band_budget(candidates, 1536);

    if (!solution.selected.empty() && solution.used_bits > 0) {
        std::vector<QuantizedSubband> quantized_subbands(32);
        std::vector<float> reconstructed(len, 0.0f);
        std::vector<SpectralSubband> spectral_subbands;
        size_t payload_bits = 0;
        size_t last_coded = 0;

        for (size_t band = 0; band < solution.selected.size(); ++band) {
            size_t start = ATRAC3_SUBBAND_TAB[band];
            size_t end = ATRAC3_SUBBAND_TAB[band + 1];
            size_t band_len = end - start;
            if (band_len > 0) {
                auto recon = solution.selected[band].dequantized(band_len);
                for (size_t i = 0; i < band_len; ++i) {
                    reconstructed[start + i] = recon[i];
                }
            }
            quantized_subbands[band] = solution.selected[band];
            if (solution.selected[band].table_index != 0) {
                last_coded = band;
            }
        }

        for (size_t band = 0; band <= last_coded; ++band) {
            SpectralSubband subband;
            subband.table_index = quantized_subbands[band].table_index;
            if (quantized_subbands[band].has_scale_factor && quantized_subbands[band].table_index != 0) {
                subband.scale_factor_index = new uint8_t(quantized_subbands[band].scale_factor_index);
                encode_mantissas_to_payload(quantized_subbands[band].table_index, coding_mode, quantized_subbands[band].mantissas, subband.payload);
            } else {
                subband.scale_factor_index = nullptr;
            }
            spectral_subbands.push_back(subband);
            payload_bits += quantized_subbands[band].payload_bits;
        }

        if (spectral_subbands.empty()) {
            spectral_subbands.push_back(SpectralSubband());
            spectral_subbands[0].table_index = 0;
            spectral_subbands[0].scale_factor_index = nullptr;
        }

        float mse = 0.0f;
        for (size_t i = 0; i < len; ++i) {
            float err = coefficients[i] - reconstructed[i];
            mse += err * err;
        }
        mse /= len;

        SpectrumEncoding result;
        result.spectral_unit.coding_mode = coding_mode;
        result.spectral_unit.subbands = spectral_subbands;
        result.quantized_subbands = quantized_subbands;
        result.reconstructed = reconstructed;
        result.payload_bits = payload_bits;
        result.mse = mse;
        return result;
    }

    std::vector<QuantizedSubband> quantized_subbands(32);
    std::vector<float> reconstructed(len, 0.0f);
    size_t last_coded = 0;

    for (int band = 0; band < 32; ++band) {
        size_t start = ATRAC3_SUBBAND_TAB[band];
        size_t end = ATRAC3_SUBBAND_TAB[band + 1];
        size_t band_len = end - start;
        if (band_len == 0) continue;

        auto quantized = choose_subband_encoding(coefficients + start, band_len, coding_mode, options);
        auto recon = quantized.dequantized(band_len);
        for (size_t i = 0; i < band_len; ++i) {
            reconstructed[start + i] = recon[i];
        }
        quantized_subbands[band] = quantized;
        if (quantized.table_index != 0) {
            last_coded = band;
        }
    }

    std::vector<SpectralSubband> spectral_subbands;
    size_t payload_bits = 0;
    for (size_t band = 0; band <= last_coded; ++band) {
        SpectralSubband subband;
        subband.table_index = quantized_subbands[band].table_index;
        if (quantized_subbands[band].has_scale_factor && quantized_subbands[band].table_index != 0) {
            subband.scale_factor_index = new uint8_t(quantized_subbands[band].scale_factor_index);
            encode_mantissas_to_payload(quantized_subbands[band].table_index, coding_mode, quantized_subbands[band].mantissas, subband.payload);
        } else {
            subband.scale_factor_index = nullptr;
        }
        spectral_subbands.push_back(subband);
        payload_bits += quantized_subbands[band].payload_bits;
    }

    if (spectral_subbands.empty()) {
        spectral_subbands.push_back(SpectralSubband());
        spectral_subbands[0].table_index = 0;
        spectral_subbands[0].scale_factor_index = nullptr;
    }

    float mse = 0.0f;
    for (size_t i = 0; i < len; ++i) {
        float err = coefficients[i] - reconstructed[i];
        mse += err * err;
    }
    mse /= len;

    SpectrumEncoding result;
    result.spectral_unit.coding_mode = coding_mode;
    result.spectral_unit.subbands = spectral_subbands;
    result.quantized_subbands = quantized_subbands;
    result.reconstructed = reconstructed;
    result.payload_bits = payload_bits;
    result.mse = mse;
    return result;
}

ChannelSoundUnit build_basic_sound_unit_from_encoding(const SpectrumEncoding& encoding) {
    size_t subband_count = encoding.spectral_unit.subbands.size();
    size_t coded_qmf_bands = 1;
    if (subband_count > 8) coded_qmf_bands = 2;
    if (subband_count > 16) coded_qmf_bands = 3;
    if (subband_count > 24) coded_qmf_bands = 4;

    ChannelSoundUnit unit;
    unit.coded_qmf_bands = static_cast<uint8_t>(coded_qmf_bands);
    unit.gain_bands.resize(coded_qmf_bands);
    unit.tonal_mode_selector = TonalCodingModeSelector::AllVlc;
    unit.tonal_components.clear();
    unit.spectrum = encoding.spectral_unit;
    return unit;
}

uint8_t fast_peak_to_sf_index(const float* coefficients, size_t len) {
    uint32_t max_abs_bits = 0;
    for (size_t i = 0; i < len; ++i) {
        float c = coefficients[i];
        uint32_t bits = reinterpret_cast<const uint32_t&>(c);
        uint32_t doubled = bits << 1;  // shift away sign, |c|*2
        if (doubled > max_abs_bits) {
            max_abs_bits = doubled;
        }
    }
    if (max_abs_bits == 0) {
        return 0;
    }
    int32_t exponent = (max_abs_bits >> 24); // 8-bit exponent of doubled value
    uint32_t mantissa = max_abs_bits & 0x00FFFFFF; // 24-bit mantissa
    
    int32_t sf_index = 3 * exponent - 364;
    
    if (mantissa > 0x00965FE9) {
        sf_index += 1;
    } else if (mantissa < 0x00428A30) {
        sf_index -= 1;
    }
    
    if (sf_index < 0 || sf_index > 63) {
        return 0;
    }
    return static_cast<uint8_t>(sf_index);
}

std::vector<uint8_t> compute_group_sf_indices(const float* coefficients, size_t len) {
    std::vector<uint8_t> result;
    for (size_t i = 0; i < len; i += 4) {
        size_t glen = std::min((size_t)4, len - i);
        result.push_back(fast_peak_to_sf_index(coefficients + i, glen));
    }
    return result;
}

int32_t estimate_band_bit_cost(const uint8_t* group_sf_indices, size_t group_count, uint8_t tbl_index, int32_t sf_index, size_t subband_width) {
    if (tbl_index == 0) return 0;
    int32_t threshold = sf_index - SF_THRESHOLD_OFFSET[tbl_index];
    int32_t cost = BASE_COST_SCALAR[tbl_index] * static_cast<int32_t>(subband_width) + 60;
    for (size_t i = 0; i < group_count; ++i) {
        if (static_cast<int32_t>(group_sf_indices[i]) < threshold) {
            cost -= COST_REDUCTION[tbl_index];
        }
    }
    if (cost < 0) cost = 0;
    return cost;
}

void encode_mantissas_to_writer(class BitWriter& writer, uint8_t selector, CodingMode coding_mode, const int8_t* mantissas, size_t count) {
    if (selector < 1 || selector > 7) return;
    
    if (coding_mode == CodingMode::Clc) {
        if (selector == 1) {
            for (size_t i = 0; i + 1 < count; i += 2) {
                uint8_t hi = clc_symbol_index(mantissas[i]);
                uint8_t lo = clc_symbol_index(mantissas[i + 1]);
                writer.write_bits((hi << 2) | lo, ATRAC3_CLC_LENGTH_TAB[1]);
            }
        } else {
            uint8_t width = ATRAC3_CLC_LENGTH_TAB[selector];
            for (size_t i = 0; i < count; ++i) {
                uint32_t val = twos_complement_bits(mantissas[i], width);
                writer.write_bits(val, width);
            }
        }
    } else {
        if (selector == 1) {
            for (size_t i = 0; i + 1 < count; i += 2) {
                int8_t symbol = vlc_pair_symbol(mantissas[i], mantissas[i + 1]);
                const HuffmanEntry* entry = find_huffman_entry(1, symbol);
                if (entry) {
                    writer.write_bits(entry->code, entry->bits);
                }
            }
        } else {
            for (size_t i = 0; i < count; ++i) {
                const HuffmanEntry* entry = find_huffman_entry(selector, mantissas[i]);
                if (entry) {
                    writer.write_bits(entry->code, entry->bits);
                }
            }
        }
    }
}

size_t coded_qmf_bands_for_subband_count(size_t subband_count) {
    size_t last_end = ATRAC3_SUBBAND_TAB[std::min(subband_count, (size_t)32)];
    size_t bands = ((last_end - 1) >> 8) + 1;
    return std::clamp(bands, (size_t)1, (size_t)4);
}

size_t fixed_sound_unit_bits(size_t subband_count) {
    size_t cqb = coded_qmf_bands_for_subband_count(subband_count);
    return 6 + 2 + cqb * 3 + 5 + 5 + 1;
}

size_t candidate_total_bits(const QuantizedSubband* candidate) {
    if (candidate->table_index == 0) return 3;
    return 3 + 6 + candidate->payload_bits;
}

static int subband_index_for_position(size_t pos) {
    for (int band = 0; band < 32; ++band) {
        if (pos >= ATRAC3_SUBBAND_TAB[band] && pos < ATRAC3_SUBBAND_TAB[band + 1]) {
            return band;
        }
    }
    return -1;
}

TonalExtractionResult extract_tonal_components(float* residual, size_t residual_len, size_t budget_bits, uint8_t coded_qmf_bands, CodingMode coding_mode, size_t max_entries) {
    TonalExtractionResult result;
    result.tonal_mode_selector = TonalCodingModeSelector::AllVlc;
    result.tonal_bits = 0;
    result.coded_qmf_bands = coded_qmf_bands;
    for (int i = 0; i < 32; ++i) result.tonal_subbands[i] = false;
    
    size_t qmf_bands = coded_qmf_bands;
    size_t total_cells = qmf_bands * 4;
    size_t spectral_end = std::min(qmf_bands * 256, residual_len);
    float abs_threshold = 2.0f;
    uint8_t qstep = 7;
    uint8_t clc_width = ATRAC3_CLC_LENGTH_TAB[qstep];
    float imq = ATRAC3_INV_MAX_QUANT[qstep];
    int32_t max_mantissa = (1i32 << (clc_width - 1)) - 1;
    int32_t min_mantissa = -(1i32 << (clc_width - 1));
    
    size_t entry_bits = 12 + static_cast<size_t>(clc_width) * 4;
    size_t new_band_bits = 12;
    size_t base_bits = 5 + 2 + qmf_bands + 3 + 3;
    size_t tonal_budget = budget_bits / 4;
    
    if (tonal_budget < base_bits + entry_bits + new_band_bits) {
        return result;
    }
    
    std::vector<std::vector<TonalEntry>> cells(total_cells);
    std::vector<bool> band_active(qmf_bands, false);
    size_t total_bits_used = base_bits;
    size_t total_entries_count = 0;
    
    std::vector<std::pair<size_t, float>> candidates;
    size_t pos = 0;
    while (pos + 3 < spectral_end) {
        size_t end = std::min(pos + 4, spectral_end);
        float peak_val = 0.0f;
        for (size_t i = pos; i < end; ++i) {
            peak_val = std::max(peak_val, std::abs(residual[i]));
        }
        if (peak_val >= abs_threshold) {
            candidates.push_back({pos, peak_val});
        }
        pos += 4;
    }
    
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });
    
    for (const auto& [pos, _] : candidates) {
        if (total_entries_count >= max_entries) break;
        
        size_t end = std::min(pos + 4, spectral_end);
        float peak_val = 0.0f;
        for (size_t i = pos; i < end; ++i) {
            peak_val = std::max(peak_val, std::abs(residual[i]));
        }
        if (peak_val < abs_threshold) continue;
        
        uint8_t sf_idx = optimal_sf_index_for_peak(peak_val, qstep);
        float sf_val = scale_factor(sf_idx);
        float scale = sf_val * imq;
        if (scale < 1e-12f) continue;
        
        std::array<int32_t, 4> mantissas = {0, 0, 0, 0};
        bool all_zero = true;
        for (size_t i = 0; i < 4; ++i) {
            if (pos + i < spectral_end) {
                mantissas[i] = round_ties_even(residual[pos + i] / scale);
                mantissas[i] = std::clamp(mantissas[i], min_mantissa, max_mantissa);
                if (mantissas[i] != 0) all_zero = false;
            }
        }
        if (all_zero) continue;
        
        size_t qmf_band_idx = pos / 256;
        size_t cell_idx = pos / 64;
        if (qmf_band_idx >= qmf_bands || cell_idx >= total_cells) continue;
        if (cells[cell_idx].size() >= 7) continue;
        
        size_t band_cost = band_active[qmf_band_idx] ? 0 : new_band_bits;
        if (total_bits_used + band_cost + entry_bits > tonal_budget) continue;
        
        for (size_t i = 0; i < 4; ++i) {
            if (pos + i < spectral_end) {
                residual[pos + i] -= static_cast<float>(mantissas[i]) * scale;
            }
        }
        
        TonalEntry entry;
        entry.scale_factor_index = sf_idx;
        entry.position = static_cast<uint8_t>(pos % 64);
        for (int32_t m : mantissas) {
            entry.payload.push_bits(twos_complement_bits(m, clc_width), clc_width);
        }
        cells[cell_idx].push_back(entry);
        
        if (!band_active[qmf_band_idx]) {
            band_active[qmf_band_idx] = true;
            total_bits_used += new_band_bits;
        }
        total_bits_used += entry_bits;
        total_entries_count += 1;
        
        int subband = subband_index_for_position(pos);
        if (subband >= 0 && subband < 32) {
            result.tonal_subbands[subband] = true;
        }
    }
    
    if (total_entries_count == 0) {
        return result;
    }
    
    TonalComponent component;
    component.band_flags = band_active;
    component.coded_values_minus_one = 3;
    component.quant_step_index = qstep;
    component.coding_mode = nullptr;
    for (const auto& cell_entries : cells) {
        TonalCell cell;
        cell.entries = cell_entries;
        component.cells.push_back(cell);
    }
    
    result.tonal_mode_selector = TonalCodingModeSelector::AllClc;
    result.tonal_components.push_back(component);
    result.tonal_bits = total_bits_used;
    
    return result;
}

}
