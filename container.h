#pragma once

#include <vector>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <optional>
#include "common.h"
#include "prototype.h"
#include "sound_unit.h"

namespace atrac3 {

    constexpr uint32_t ATRAC3_TARGET_SAMPLE_RATE = 44100;

    enum class Atrac3Bitrate {
        Kbps66 = 0x00,
        Kbps105 = 0x01,
        Kbps132 = 0x02,
        Kbps165 = 0x03,
    };

    struct Atrac3ContainerOptions {
        std::optional<Atrac3Bitrate> bitrate;
        Atrac3ContainerOptions() : bitrate(std::nullopt) {}
    };

    struct Atrac3Container {
        uint32_t sample_rate;
        uint16_t channel_count;
        uint32_t frame_count;
        Atrac3Bitrate bitrate;
        std::vector<uint8_t> frames;

        std::vector<uint8_t> encode() const;
        static Atrac3Container decode(const std::vector<uint8_t>& data);
    };

    inline uint16_t Atrac3Bitrate_block_align(Atrac3Bitrate self, uint16_t channels) {
        uint16_t per_channel;
        switch (self) {
        case Atrac3Bitrate::Kbps66: per_channel = 96; break;
        case Atrac3Bitrate::Kbps105: per_channel = 152; break;
        case Atrac3Bitrate::Kbps132: per_channel = 192; break;
        case Atrac3Bitrate::Kbps165: per_channel = 248; break;
        default: per_channel = 192;
        }
        return per_channel * channels;
    }

    inline size_t Atrac3Bitrate_bits_per_channel_frame(Atrac3Bitrate self, uint16_t channels) {
        uint16_t block_align = Atrac3Bitrate_block_align(self, channels);
        return (block_align * 8) / channels;
    }

    inline Atrac3Bitrate Atrac3Bitrate_from_bits_per_channel_frame(size_t bits) {
        switch (bits) {
        case 768: return Atrac3Bitrate::Kbps66;
        case 1216: return Atrac3Bitrate::Kbps105;
        case 1536: return Atrac3Bitrate::Kbps132;
        case 1984: return Atrac3Bitrate::Kbps165;
        default: throw std::runtime_error("no ATRAC3 bitrate with this bit budget per channel frame");
        }
    }

    inline uint32_t Atrac3Bitrate_kbps(Atrac3Bitrate self, uint16_t channels) {
        if (channels == 1) {
            switch (self) {
            case Atrac3Bitrate::Kbps66: return 33;
            case Atrac3Bitrate::Kbps105: return 52;
            case Atrac3Bitrate::Kbps132: return 66;
            case Atrac3Bitrate::Kbps165: return 83;
            default: throw std::runtime_error("invalid ATRAC3 bitrate");
            }
        }
        else {
            switch (self) {
            case Atrac3Bitrate::Kbps66: return 66;
            case Atrac3Bitrate::Kbps105: return 105;
            case Atrac3Bitrate::Kbps132: return 132;
            case Atrac3Bitrate::Kbps165: return 165;
            default: throw std::runtime_error("invalid ATRAC3 bitrate");
            }
        }
    }

    constexpr uint32_t FOURCC_WAVE = 0x45564157;
    constexpr uint32_t FOURCC_fmt_ = 0x20746d66;
    constexpr uint32_t FOURCC_data = 0x61746164;
    constexpr uint32_t FOURCC_fact = 0x74636166;
    constexpr uint32_t FOURCC_RIFF = 0x46464952;
    constexpr uint16_t WAVE_FORMAT_ATRAC3 = 0x0270;
    constexpr size_t ATRAC3_WAVE_EXTRA_SIZE = 14;
    constexpr size_t ATRAC3_FACT_SIZE = 4;

    struct WaveFormatEx {
        uint16_t w_format_tag;
        uint16_t n_channels;
        uint32_t n_samples_per_sec;
        uint32_t n_avg_bytes_per_sec;
        uint16_t n_block_align;
        uint16_t w_bits_per_sample;
        uint16_t cb_size;
    };

    struct RiffChunkHeader {
        uint32_t id;
        uint32_t size;
    };

    struct WaveFormatExExtraAtrac3 {
        uint16_t unknown1;
        uint16_t unknown2;
        uint32_t channel_config;
        uint16_t unknown3;
        uint16_t joint_stereo;
        uint16_t unknown4;
    };

    inline void push_u16_le(std::vector<uint8_t>& bytes, uint16_t value) {
        bytes.push_back(static_cast<uint8_t>(value & 0xFF));
        bytes.push_back(static_cast<uint8_t>((value >> 8) & 0xFF));
    }

    inline void push_u32_le(std::vector<uint8_t>& bytes, uint32_t value) {
        bytes.push_back(static_cast<uint8_t>(value & 0xFF));
        bytes.push_back(static_cast<uint8_t>((value >> 8) & 0xFF));
        bytes.push_back(static_cast<uint8_t>((value >> 16) & 0xFF));
        bytes.push_back(static_cast<uint8_t>((value >> 24) & 0xFF));
    }

    inline Atrac3Bitrate choose_smallest_fitting_bitrate(size_t max_channel_bytes, uint16_t channels) {
        std::vector<Atrac3Bitrate> bitrates = {
            Atrac3Bitrate::Kbps66,
            Atrac3Bitrate::Kbps105,
            Atrac3Bitrate::Kbps132,
            Atrac3Bitrate::Kbps165
        };
        for (Atrac3Bitrate bitrate : bitrates) {
            size_t slot_size = Atrac3Bitrate_block_align(bitrate, channels) / channels;
            if (max_channel_bytes <= slot_size) {
                return bitrate;
            }
        }
        throw std::runtime_error("encoded channel sound unit exceeds supported ATRAC3 budgets");
    }

    inline uint16_t coding_mode_extradata(const PrototypeEncodeResult& encoded) {
        return 0;
    }

    inline std::vector<uint8_t> wrap_prototype_in_riff_at3(
        const PrototypeEncodeResult& encoded,
        const Atrac3ContainerOptions& options) {

        if (encoded.channel_count < 1 || encoded.channel_count > 2) {
            throw std::runtime_error("prototype ATRAC3 container currently supports 1 or 2 channels");
        }
        if (encoded.frame_count == 0) {
            throw std::runtime_error("cannot wrap an empty encode result");
        }

        size_t max_channel_bytes = 0;
        for (const auto& frame : encoded.frames) {
            for (const auto& channel : frame.channels) {
                if (channel.bytes.size() > max_channel_bytes) {
                    max_channel_bytes = channel.bytes.size();
                }
            }
        }

        Atrac3Bitrate bitrate;
        if (options.bitrate.has_value()) {
            bitrate = options.bitrate.value();
        }
        else {
            bitrate = choose_smallest_fitting_bitrate(max_channel_bytes, static_cast<uint16_t>(encoded.channel_count));
        }

        uint16_t block_align = Atrac3Bitrate_block_align(bitrate, static_cast<uint16_t>(encoded.channel_count));
        size_t channel_slot_size = block_align / encoded.channel_count;
        uint32_t avg_bytes_per_sec = ((static_cast<uint32_t>(block_align) * ATRAC3_TARGET_SAMPLE_RATE) + (SAMPLES_PER_FRAME / 2)) / SAMPLES_PER_FRAME;
        uint32_t sample_count = static_cast<uint32_t>(encoded.frame_count) * SAMPLES_PER_FRAME;
        uint16_t coding_mode = coding_mode_extradata(encoded);
        uint32_t data_size = static_cast<uint32_t>(encoded.frame_count) * static_cast<uint32_t>(block_align);

        uint16_t joint_stereo_flag = 0;
        if (encoded.channel_count == 2 && bitrate == Atrac3Bitrate::Kbps66) {
            joint_stereo_flag = 1;
        }

        std::vector<uint8_t> bytes;
        bytes.reserve(76 + data_size);

        uint32_t riff_size = 4 + (8 + 32) + (8 + 8) + (8 + data_size);

        bytes.insert(bytes.end(), { 'R', 'I', 'F', 'F' });
        push_u32_le(bytes, riff_size);
        bytes.insert(bytes.end(), { 'W', 'A', 'V', 'E' });

        bytes.insert(bytes.end(), { 'f', 'm', 't', ' ' });
        push_u32_le(bytes, 32);
        push_u16_le(bytes, WAVE_FORMAT_ATRAC3);
        push_u16_le(bytes, static_cast<uint16_t>(encoded.channel_count));
        push_u32_le(bytes, ATRAC3_TARGET_SAMPLE_RATE);
        push_u32_le(bytes, avg_bytes_per_sec);
        push_u16_le(bytes, block_align);
        push_u16_le(bytes, 0);
        push_u16_le(bytes, ATRAC3_WAVE_EXTRA_SIZE);
        push_u16_le(bytes, 1);
        push_u32_le(bytes, 0x1000);
        push_u16_le(bytes, coding_mode);
        push_u16_le(bytes, joint_stereo_flag);
        push_u16_le(bytes, 1);
        push_u16_le(bytes, 0);

        bytes.insert(bytes.end(), { 'f', 'a', 'c', 't' });
        push_u32_le(bytes, 8);
        push_u32_le(bytes, sample_count);
        push_u32_le(bytes, SAMPLES_PER_FRAME);

        bytes.insert(bytes.end(), { 'd', 'a', 't', 'a' });
        push_u32_le(bytes, data_size);

        for (const auto& frame : encoded.frames) {
            if (frame.channels.size() != encoded.channel_count) {
                throw std::runtime_error("frame channel count mismatch");
            }
            for (const auto& channel : frame.channels) {
                size_t write_size = std::min(channel.bytes.size(), channel_slot_size);
                bytes.insert(bytes.end(), channel.bytes.begin(), channel.bytes.begin() + write_size);
                size_t padding = channel_slot_size - write_size;
                for (size_t i = 0; i < padding; ++i) {
                    bytes.push_back(0);
                }
            }
        }

        return bytes;
    }

    inline std::vector<uint8_t> Atrac3Container::encode() const {
        std::vector<uint8_t> out;
        out.reserve(12 + 8 + sizeof(WaveFormatEx) + ATRAC3_WAVE_EXTRA_SIZE + 8 + ATRAC3_FACT_SIZE + 8 + frames.size());

        size_t bits_per_channel_frame = Atrac3Bitrate_bits_per_channel_frame(bitrate, channel_count);
        size_t bytes_per_channel_frame = (bits_per_channel_frame + 7) / 8;
        size_t bytes_per_frame = bytes_per_channel_frame * channel_count;

        WaveFormatEx fmt_header = {};
        fmt_header.w_format_tag = WAVE_FORMAT_ATRAC3;
        fmt_header.n_channels = channel_count;
        fmt_header.n_samples_per_sec = sample_rate;
        fmt_header.n_avg_bytes_per_sec = static_cast<uint32_t>(bytes_per_frame * (static_cast<size_t>(sample_rate) / SAMPLES_PER_FRAME));
        fmt_header.n_block_align = static_cast<uint16_t>(bytes_per_frame);
        fmt_header.w_bits_per_sample = 0;
        fmt_header.cb_size = ATRAC3_WAVE_EXTRA_SIZE;

        WaveFormatExExtraAtrac3 fmt_extra = {};
        fmt_extra.unknown1 = 0;
        fmt_extra.unknown2 = static_cast<uint16_t>(bits_per_channel_frame);
        fmt_extra.channel_config = 0;
        fmt_extra.unknown3 = 1;
        fmt_extra.joint_stereo = 0;
        fmt_extra.unknown4 = 0;

        size_t fact_payload_size = 4;
        uint32_t fact_payload = static_cast<uint32_t>(frame_count);

        size_t data_payload_size = frames.size();

        size_t riff_size = 4 + (8 + sizeof(WaveFormatEx) + ATRAC3_WAVE_EXTRA_SIZE) + (8 + fact_payload_size) + (8 + data_payload_size);

        out.resize(12);
        std::memcpy(&out[0], &FOURCC_RIFF, 4);
        std::memcpy(&out[4], &riff_size, 4);
        std::memcpy(&out[8], &FOURCC_WAVE, 4);

        size_t fmt_size = sizeof(WaveFormatEx) + ATRAC3_WAVE_EXTRA_SIZE;
        out.resize(out.size() + 8 + fmt_size);
        size_t offset = 12;
        std::memcpy(&out[offset], &FOURCC_fmt_, 4);
        std::memcpy(&out[offset + 4], &fmt_size, 4);
        offset += 8;
        std::memcpy(&out[offset], &fmt_header, sizeof(fmt_header));
        offset += sizeof(fmt_header);
        std::memcpy(&out[offset], &fmt_extra, sizeof(fmt_extra));
        offset += sizeof(fmt_extra);

        out.resize(out.size() + 8 + fact_payload_size);
        std::memcpy(&out[offset], &FOURCC_fact, 4);
        std::memcpy(&out[offset + 4], &fact_payload_size, 4);
        offset += 8;
        std::memcpy(&out[offset], &fact_payload, 4);
        offset += 4;

        out.resize(out.size() + 8 + data_payload_size);
        std::memcpy(&out[offset], &FOURCC_data, 4);
        std::memcpy(&out[offset + 4], &data_payload_size, 4);
        offset += 8;
        std::memcpy(&out[offset], frames.data(), frames.size());

        return out;
    }

    inline Atrac3Container Atrac3Container::decode(const std::vector<uint8_t>& data) {
        size_t offset = 0;

        if (data.size() < 12) {
            throw std::runtime_error("not a WAVE file");
        }

        uint32_t riff_id, wave_id;
        std::memcpy(&riff_id, &data[offset], 4);
        offset += 4;

        uint32_t riff_size;
        std::memcpy(&riff_size, &data[offset], 4);
        offset += 4;

        std::memcpy(&wave_id, &data[offset], 4);
        offset += 4;

        if (riff_id != FOURCC_RIFF || wave_id != FOURCC_WAVE) {
            throw std::runtime_error("not a WAVE file");
        }

        std::optional<WaveFormatEx> fmt_header;
        std::optional<WaveFormatExExtraAtrac3> fmt_extra;
        std::optional<uint32_t> fact_payload;
        std::vector<uint8_t> data_payload;

        while (offset + 8 <= data.size()) {
            RiffChunkHeader chunk;
            std::memcpy(&chunk.id, &data[offset], 4);
            offset += 4;
            std::memcpy(&chunk.size, &data[offset], 4);
            offset += 4;

            size_t chunk_end = offset + chunk.size;
            if (chunk_end > data.size()) {
                chunk_end = data.size();
            }

            if (chunk.id == FOURCC_fmt_) {
                if (chunk.size < sizeof(WaveFormatEx)) {
                    throw std::runtime_error("not a WAVE file");
                }
                WaveFormatEx fmt;
                std::memcpy(&fmt, &data[offset], sizeof(fmt));
                offset += sizeof(fmt);

                if (fmt.w_format_tag != WAVE_FORMAT_ATRAC3) {
                    throw std::runtime_error("not an ATRAC3 file");
                }

                if (fmt.cb_size != ATRAC3_WAVE_EXTRA_SIZE) {
                    throw std::runtime_error("not an ATRAC3 file");
                }

                if (chunk.size < sizeof(WaveFormatEx) + ATRAC3_WAVE_EXTRA_SIZE) {
                    throw std::runtime_error("not an ATRAC3 file");
                }

                WaveFormatExExtraAtrac3 extra;
                std::memcpy(&extra, &data[offset], sizeof(extra));
                offset += sizeof(extra);

                fmt_header = fmt;
                fmt_extra = extra;
            }
            else if (chunk.id == FOURCC_fact) {
                if (chunk.size < 4) {
                    throw std::runtime_error("not an ATRAC3 file");
                }
                uint32_t fact;
                std::memcpy(&fact, &data[offset], 4);
                offset += 4;
                fact_payload = fact;
            }
            else if (chunk.id == FOURCC_data) {
                data_payload.resize(chunk.size);
                std::memcpy(data_payload.data(), &data[offset], chunk.size);
                offset += chunk.size;
            }
            else {
                offset = chunk_end;
            }

            offset = chunk_end;
        }

        if (!fmt_header || !fmt_extra || !fact_payload || data_payload.empty()) {
            throw std::runtime_error("not an ATRAC3 file");
        }

        return Atrac3Container{
            fmt_header->n_samples_per_sec,
            fmt_header->n_channels,
            *fact_payload,
            Atrac3Bitrate_from_bits_per_channel_frame(fmt_extra->unknown2),
            std::move(data_payload),
        };
    }

}