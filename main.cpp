#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <optional>
#include <vector>
#include <cmath>
#include <algorithm>
#include "wav.h"
#include "prototype.h"
#include "container.h"

using namespace atrac3;

std::vector<float> resample_audio(const std::vector<float>& samples, int src_sample_rate, int dst_sample_rate, int channels) {
    if (src_sample_rate == dst_sample_rate) {
        return samples;
    }

    std::cout << "Resampling from " << src_sample_rate << " Hz to " << dst_sample_rate << " Hz..." << std::endl;

    double factor = static_cast<double>(dst_sample_rate) / src_sample_rate;
    size_t src_samples_per_channel = samples.size() / channels;
    size_t dst_samples_per_channel = static_cast<size_t>(src_samples_per_channel * factor);
    
    std::vector<float> result(dst_samples_per_channel * channels);

    for (int ch = 0; ch < channels; ch++) {
        for (size_t i = 0; i < dst_samples_per_channel; i++) {
            double position = static_cast<double>(i) / factor;
            size_t index = static_cast<size_t>(position);
            double frac = position - index;

            float y0, y1, y2, y3;
            
            if (index == 0) {
                y0 = samples[ch];
                y1 = samples[ch];
                y2 = samples[ch + channels];
                y3 = samples[ch + 2 * channels];
            } else if (index >= src_samples_per_channel - 1) {
                y0 = samples[ch + (src_samples_per_channel - 2) * channels];
                y1 = samples[ch + (src_samples_per_channel - 1) * channels];
                y2 = samples[ch + (src_samples_per_channel - 1) * channels];
                y3 = samples[ch + (src_samples_per_channel - 1) * channels];
            } else {
                y0 = samples[ch + (index - 1) * channels];
                y1 = samples[ch + index * channels];
                y2 = samples[ch + (index + 1) * channels];
                y3 = samples[ch + (index + 2) * channels];
            }

            double a0 = 2.0 * y1;
            double a1 = -y0 + y2;
            double a2 = 2.0 * y0 - 5.0 * y1 + 4.0 * y2 - y3;
            double a3 = -y0 + 3.0 * y1 - 3.0 * y2 + y3;

            double v = (a0 + frac * (a1 + frac * (a2 + frac * a3))) * 0.5;
            result[i * channels + ch] = static_cast<float>(v);
        }
    }

    return result;
}

struct CliOptions {
    std::string input_path;
    std::string output_path;
    std::optional<Atrac3Bitrate> bitrate;
    bool encode_mode;
};

CliOptions parse_args(int argc, char* argv[]) {
    CliOptions opts{};
    opts.encode_mode = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-e" || arg == "--encode") {
            opts.encode_mode = true;
        }
        else if (arg == "-b" || arg == "--bitrate") {
            if (i + 1 < argc) {
                int br = std::stoi(argv[++i]);
                switch (br) {
                case 66: opts.bitrate = Atrac3Bitrate::Kbps66; break;
                case 105: opts.bitrate = Atrac3Bitrate::Kbps105; break;
                case 132: opts.bitrate = Atrac3Bitrate::Kbps132; break;
                default:
                    std::cerr << "Unsupported bitrate " << br << ", using 132 kbps" << std::endl;
                    opts.bitrate = Atrac3Bitrate::Kbps132;
                }
            }
        }
        else if (arg[0] != '-') {
            if (opts.input_path.empty()) {
                opts.input_path = arg;
            }
            else if (opts.output_path.empty()) {
                opts.output_path = arg;
            }
        }
    }

    if (opts.input_path.empty()) {
        throw std::runtime_error("No input file specified");
    }

    if (opts.output_path.empty()) {
        size_t dot_pos = opts.input_path.rfind('.');
        if (dot_pos != std::string::npos) {
            opts.output_path = opts.input_path.substr(0, dot_pos) + ".at3";
        }
        else {
            opts.output_path = opts.input_path + ".at3";
        }
    }

    return opts;
}

void encode_file(const CliOptions& opts) {
    std::cout << "Reading WAV file: " << opts.input_path << std::endl;
    auto wav = read_wav(opts.input_path);

    std::cout << "WAV file info:" << std::endl;
    std::cout << "  Sample rate: " << wav.sample_rate << " Hz" << std::endl;
    std::cout << "  Channels: " << wav.channels << std::endl;

    std::vector<float> samples = wav.samples;
    int sample_rate = wav.sample_rate;

    if (wav.sample_rate != ATRAC3_TARGET_SAMPLE_RATE) {
        samples = resample_audio(wav.samples, wav.sample_rate, ATRAC3_TARGET_SAMPLE_RATE, wav.channels);
        sample_rate = ATRAC3_TARGET_SAMPLE_RATE;
    }

    size_t pcm_frames = samples.size() / wav.channels;
    size_t atrac3_frames = (pcm_frames + 1023) / 1024;
    std::cout << "  ATRAC3 Frames: " << atrac3_frames << std::endl;

    Atrac3Bitrate actual_bitrate = Atrac3Bitrate::Kbps132;
    if (opts.bitrate.has_value()) {
        actual_bitrate = opts.bitrate.value();
    }

    size_t total_bits_per_frame = Atrac3Bitrate_block_align(actual_bitrate, wav.channels) * 8;
    size_t target_bits_per_channel = total_bits_per_frame / wav.channels;

    std::cout << "  Target bits per channel frame: " << target_bits_per_channel << std::endl;
    std::cout << "  Block align: " << Atrac3Bitrate_block_align(actual_bitrate, wav.channels) << std::endl;

    PrototypeEncoder encoder(wav.channels);
    PrototypeOptions proto_opts;
    proto_opts.coding_mode = CodingMode::Clc;
    proto_opts.lambda = 0.0001f;
    proto_opts.frame_limit = atrac3_frames;
    proto_opts.start_frame = 0;
    proto_opts.flush_frames = 0;
    proto_opts.target_bits_per_channel = target_bits_per_channel;
    proto_opts.joint_stereo = (wav.channels == 2 && actual_bitrate == Atrac3Bitrate::Kbps66);

    auto result = encoder.encode_wav(
        samples.data(),
        samples.size(),
        sample_rate,
        wav.channels,
        proto_opts
    );

    std::cout << "Encoded " << result.frame_count << " frames" << std::endl;

    Atrac3ContainerOptions container_opts;
    container_opts.bitrate = actual_bitrate;

    auto at3_data = wrap_prototype_in_riff_at3(result, container_opts);

    std::cout << "Writing output file: " << opts.output_path << std::endl;
    std::ofstream out_file(opts.output_path, std::ios::binary);
    if (!out_file.is_open()) {
        throw std::runtime_error("Failed to open output file");
    }
    out_file.write(reinterpret_cast<const char*>(at3_data.data()), at3_data.size());
    out_file.close();

    std::cout << "Encoding complete!" << std::endl;
    std::cout << "Output: " << opts.output_path << std::endl;
    std::cout << "File size: " << at3_data.size() << " bytes" << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "ATRAC3 Encoder" << std::endl;

    std::string input_wav;
    int bitrate = 132;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-e" && i + 1 < argc) {
            input_wav = argv[++i];
        }
        else if (arg == "-b" && i + 1 < argc) {
            bitrate = std::stoi(argv[++i]);
        }
        else if (arg == "-h") {
            std::cout << "Usage: at3codec -e input.wav [-b 66|105|132]" << std::endl;
            return 0;
        }
    }

    if (input_wav.empty()) {
        std::cerr << "Error: No input file" << std::endl;
        return 1;
    }

    std::string output_at3 = input_wav.substr(0, input_wav.rfind('.')) + ".at3";

    auto wav = read_wav(input_wav);
    std::vector<float> samples = wav.samples;
    int sample_rate = wav.sample_rate;

    if (wav.sample_rate != 44100) {
        samples = resample_audio(wav.samples, wav.sample_rate, 44100, wav.channels);
        sample_rate = 44100;
    }

    Atrac3Bitrate at3_bitrate;
    switch (bitrate) {
    case 66: at3_bitrate = Atrac3Bitrate::Kbps66; break;
    case 105: at3_bitrate = Atrac3Bitrate::Kbps105; break;
    default: at3_bitrate = Atrac3Bitrate::Kbps132;
    }

    size_t pcm_frames = samples.size() / wav.channels;
    size_t atrac3_frames = (pcm_frames + 1023) / 1024;
    size_t total_bits_per_frame = Atrac3Bitrate_block_align(at3_bitrate, wav.channels) * 8;
    size_t target_bits_per_channel = total_bits_per_frame / wav.channels;

    PrototypeEncoder encoder(wav.channels);
    PrototypeOptions proto_opts;
    proto_opts.coding_mode = CodingMode::Clc;
    proto_opts.lambda = 0.0001f;
    proto_opts.frame_limit = atrac3_frames;
    proto_opts.target_bits_per_channel = target_bits_per_channel;
    proto_opts.joint_stereo = (wav.channels == 2 && at3_bitrate == Atrac3Bitrate::Kbps66);

    auto result = encoder.encode_wav(samples.data(), samples.size(), sample_rate, wav.channels, proto_opts);

    Atrac3ContainerOptions container_opts;
    container_opts.bitrate = at3_bitrate;
    auto at3_data = wrap_prototype_in_riff_at3(result, container_opts);

    std::ofstream out_file(output_at3, std::ios::binary);
    out_file.write(reinterpret_cast<const char*>(at3_data.data()), at3_data.size());

    std::cout << "Done: " << output_at3 << std::endl;
    return 0;
}