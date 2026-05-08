#pragma once

#include <vector>
#include <cstdint>
#include <stdexcept>
#include <fstream>
#include <cstring>

namespace atrac3 {

    struct WavData {
        uint32_t sample_rate;
        uint16_t channels;
        std::vector<float> samples;

        size_t frames() const {
            return samples.size() / channels;
        }

        std::vector<float> channel_samples(size_t channel) const {
            if (channel >= channels) {
                throw std::out_of_range("channel out of range");
            }
            std::vector<float> result;
            result.reserve(frames());
            for (size_t i = 0; i < samples.size(); i += channels) {
                result.push_back(samples[i + channel]);
            }
            return result;
        }
    };

    inline WavData read_wav(const std::string& path) {
        std::ifstream file(path, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open WAV file");
        }

        char riff[4];
        file.read(riff, 4);
        if (std::memcmp(riff, "RIFF", 4) != 0) {
            throw std::runtime_error("Not a WAV file");
        }

        uint32_t file_size;
        file.read(reinterpret_cast<char*>(&file_size), 4);

        char wave[4];
        file.read(wave, 4);
        if (std::memcmp(wave, "WAVE", 4) != 0) {
            throw std::runtime_error("Not a WAV file");
        }

        uint16_t channels = 0;
        uint32_t sample_rate = 0;
        uint16_t bits_per_sample = 0;
        uint16_t audio_format = 0;
        std::vector<uint8_t> data;

        while (file.good()) {
            char chunk_id[4];
            file.read(chunk_id, 4);
            if (!file.good()) break;

            uint32_t chunk_size;
            file.read(reinterpret_cast<char*>(&chunk_size), 4);

            if (std::memcmp(chunk_id, "fmt ", 4) == 0) {
                file.read(reinterpret_cast<char*>(&audio_format), 2);
                file.read(reinterpret_cast<char*>(&channels), 2);
                file.read(reinterpret_cast<char*>(&sample_rate), 4);
                uint32_t byte_rate;
                file.read(reinterpret_cast<char*>(&byte_rate), 4);
                uint16_t block_align;
                file.read(reinterpret_cast<char*>(&block_align), 2);
                file.read(reinterpret_cast<char*>(&bits_per_sample), 2);

                if (chunk_size > 16) {
                    file.ignore(chunk_size - 16);
                }
            }
            else if (std::memcmp(chunk_id, "data", 4) == 0) {
                data.resize(chunk_size);
                file.read(reinterpret_cast<char*>(data.data()), chunk_size);
            }
            else {
                file.ignore(chunk_size);
            }
        }

        if (channels != 1 && channels != 2) {
            throw std::runtime_error("Only mono and stereo WAV supported");
        }

        if (audio_format != 1 && audio_format != 3) {
            throw std::runtime_error("Unsupported WAV format, only PCM and float are supported");
        }

        std::vector<float> samples;
        size_t sample_count = data.size() / (bits_per_sample / 8);

        if (sample_count % channels != 0) {
            throw std::runtime_error("Sample count not multiple of channels");
        }

        samples.reserve(sample_count);

        if (bits_per_sample == 16 && audio_format == 1) {
            const int16_t* int_samples = reinterpret_cast<const int16_t*>(data.data());
            for (size_t i = 0; i < sample_count; ++i) {
                samples.push_back(static_cast<float>(int_samples[i]) / 32768.0f);
            }
        }
        else if (bits_per_sample == 32 && audio_format == 1) {
            const int32_t* int_samples = reinterpret_cast<const int32_t*>(data.data());
            for (size_t i = 0; i < sample_count; ++i) {
                samples.push_back(static_cast<float>(int_samples[i]) / 2147483648.0f);
            }
        }
        else if (bits_per_sample == 32 && audio_format == 3) {
            const float* float_samples = reinterpret_cast<const float*>(data.data());
            for (size_t i = 0; i < sample_count; ++i) {
                samples.push_back(float_samples[i]);
            }
        }
        else {
            throw std::runtime_error("Unsupported bit depth");
        }

        return WavData{ sample_rate, channels, std::move(samples) };
    }

} // namespace atrac3