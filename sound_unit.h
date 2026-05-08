#pragma once

#include <vector>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include "bitstream.h"
#include "common.h"

namespace atrac3 {

    constexpr uint8_t SOUND_UNIT_ID = 0x28;
    constexpr size_t MAX_CODED_QMF_BANDS = 4;
    constexpr size_t MAX_CODED_SUBBANDS = 32;
    constexpr size_t MAX_TONAL_COMPONENTS = 64;
    constexpr size_t TONAL_CELLS_PER_QMF_BAND = 4;
    constexpr size_t MAX_GAIN_POINTS = 7;
    constexpr size_t MAX_TONAL_ENTRIES_PER_CELL = 7;

    enum class CodingMode {
        Vlc = 0,
        Clc = 1,
    };

    inline uint32_t CodingMode_bit(CodingMode self) {
        return static_cast<uint32_t>(self);
    }

    enum class SpectralTableKind {
        Skip,
        Pairwise,
        Single,
    };

    enum class TonalCodingModeSelector {
        AllVlc = 0,
        AllClc = 1,
        PerComponent = 3,
    };

    inline uint32_t TonalCodingModeSelector_bits(TonalCodingModeSelector self) {
        return static_cast<uint32_t>(self);
    }

    inline CodingMode* TonalCodingModeSelector_shared_mode(TonalCodingModeSelector self) {
        static CodingMode vlc = CodingMode::Vlc;
        static CodingMode clc = CodingMode::Clc;
        switch (self) {
        case TonalCodingModeSelector::AllVlc: return &vlc;
        case TonalCodingModeSelector::AllClc: return &clc;
        default: return nullptr;
        }
    }

    struct GainPoint {
        uint8_t level;
        uint8_t location;

        void write_to(BitWriter& writer) const {
            writer.write_bits(level, 4);
            writer.write_bits(location, 5);
        }
    };

    struct BitChunk {
        uint32_t value;
        uint8_t bits;
    };

    struct RawBitPayload {
        std::vector<BitChunk> chunks;

        void push_bits(uint32_t value, uint8_t bits) {
            chunks.push_back({ value, bits });
        }

        size_t bit_len() const {
            size_t len = 0;
            for (const auto& chunk : chunks) {
                len += chunk.bits;
            }
            return len;
        }

        void write_to(BitWriter& writer) const {
            for (const auto& chunk : chunks) {
                writer.write_bits(chunk.value, chunk.bits);
            }
        }
    };

    struct GainBand {
        std::vector<GainPoint> points;

        GainBand() = default;

        void write_to(BitWriter& writer) const {
            writer.write_bits(static_cast<uint32_t>(points.size()), 3);
            for (const auto& point : points) {
                point.write_to(writer);
            }
        }
    };

    struct TonalEntry {
        uint8_t scale_factor_index = 0;
        uint8_t position = 0;
        RawBitPayload payload;

        void write_to(BitWriter& writer) const {
            writer.write_bits(scale_factor_index, 6);
            writer.write_bits(position, 6);
            payload.write_to(writer);
        }
    };

    struct TonalCell {
        std::vector<TonalEntry> entries;

        TonalCell() = default;

        void write_to(BitWriter& writer) const {
            writer.write_bits(static_cast<uint32_t>(entries.size()), 3);
            for (const auto& entry : entries) {
                entry.write_to(writer);
            }
        }
    };

    struct TonalComponent {
        std::vector<bool> band_flags;
        uint8_t coded_values_minus_one;
        uint8_t quant_step_index;
        CodingMode* coding_mode;
        std::vector<TonalCell> cells;

        TonalComponent() : coding_mode(nullptr), coded_values_minus_one(0), quant_step_index(2) {}

        ~TonalComponent() {}

        void write_to(BitWriter& writer, size_t coded_qmf_bands, TonalCodingModeSelector selector) const {
            for (bool flag : band_flags) {
                writer.write_bit(flag);
            }
            writer.write_bits(coded_values_minus_one, 3);
            writer.write_bits(quant_step_index, 3);
            if (selector == TonalCodingModeSelector::PerComponent) {
                writer.write_bits(CodingMode_bit(*coding_mode), 1);
            }
            for (size_t cell_index = 0; cell_index < cells.size(); ++cell_index) {
                if (!band_flags[cell_index / TONAL_CELLS_PER_QMF_BAND]) {
                    continue;
                }
                cells[cell_index].write_to(writer);
            }
        }
    };

    struct SpectralSubband {
        uint8_t table_index;
        uint8_t* scale_factor_index;
        RawBitPayload payload;

        SpectralSubband() : table_index(0), scale_factor_index(nullptr) {}

        ~SpectralSubband() {
            if (scale_factor_index != nullptr) {
                delete scale_factor_index;
            }
        }

        SpectralSubband(const SpectralSubband& other)
            : table_index(other.table_index), scale_factor_index(nullptr), payload(other.payload) {
            if (other.scale_factor_index != nullptr) {
                scale_factor_index = new uint8_t(*other.scale_factor_index);
            }
        }

        SpectralSubband& operator=(const SpectralSubband& other) {
            if (this != &other) {
                if (scale_factor_index != nullptr) {
                    delete scale_factor_index;
                }
                table_index = other.table_index;
                payload = other.payload;
                if (other.scale_factor_index != nullptr) {
                    scale_factor_index = new uint8_t(*other.scale_factor_index);
                }
                else {
                    scale_factor_index = nullptr;
                }
            }
            return *this;
        }
    };

    struct SpectralUnit {
        CodingMode coding_mode;
        std::vector<SpectralSubband> subbands;

        SpectralUnit() : coding_mode(CodingMode::Vlc) {}

        void write_to(BitWriter& writer) const;
    };

    struct ChannelSoundUnit {
        uint8_t coded_qmf_bands;
        std::vector<GainBand> gain_bands;
        TonalCodingModeSelector tonal_mode_selector;
        std::vector<TonalComponent> tonal_components;
        SpectralUnit spectrum;

        ChannelSoundUnit() : coded_qmf_bands(1), tonal_mode_selector(TonalCodingModeSelector::AllVlc) {
            gain_bands.push_back(GainBand());
        }

        void write_to(BitWriter& writer) const;

        size_t bit_len() const {
            BitWriter writer;
            write_to(writer);
            return writer.bits_written();
        }
    };

}