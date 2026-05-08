#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace atrac3 {

constexpr size_t SAMPLES_PER_FRAME = 1024;
constexpr size_t QMF_BANDS = 4;
constexpr size_t SAMPLES_PER_BAND = SAMPLES_PER_FRAME / QMF_BANDS;

static std::array<float, 256> encode_window() {
    std::array<float, 256> w;
    for (int i = 0; i < 256; ++i) {
        float t = ((static_cast<float>(i) + 0.5f) / 256.0f) - 0.5f;
        w[i] = std::sin(t * static_cast<float>(M_PI)) + 1.0f;
    }
    return w;
}

static const std::array<float, 256> ENCODE_WINDOW = encode_window();

} // namespace atrac3
