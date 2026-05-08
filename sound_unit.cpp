#include "sound_unit.h"

namespace atrac3 {

    void SpectralUnit::write_to(BitWriter& writer) const {
        if (subbands.empty() || subbands.size() > 32) {
            throw std::runtime_error("spectral unit has invalid number of coded subbands");
        }

        writer.write_bits(static_cast<uint32_t>(subbands.size() - 1), 5);
        writer.write_bits(CodingMode_bit(coding_mode), 1);

        for (const auto& subband : subbands) {
            writer.write_bits(subband.table_index, 3);
        }

        for (const auto& subband : subbands) {
            if (subband.scale_factor_index != nullptr) {
                writer.write_bits(*subband.scale_factor_index, 6);
            }
        }

        for (const auto& subband : subbands) {
            if (subband.table_index != 0) {
                subband.payload.write_to(writer);
            }
        }
    }

    void ChannelSoundUnit::write_to(BitWriter& writer) const {
        writer.write_bits(SOUND_UNIT_ID, 6);
        writer.write_bits(coded_qmf_bands - 1, 2);

        for (const auto& band : gain_bands) {
            band.write_to(writer);
        }

        writer.write_bits(static_cast<uint32_t>(tonal_components.size()), 5);
        if (!tonal_components.empty()) {
            writer.write_bits(TonalCodingModeSelector_bits(tonal_mode_selector), 2);
            for (const auto& component : tonal_components) {
                component.write_to(writer, coded_qmf_bands, tonal_mode_selector);
            }
        }

        spectrum.write_to(writer);
    }

}