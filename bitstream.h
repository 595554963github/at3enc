#pragma once

#include <vector>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace atrac3 {

    class BitWriter {
    public:
        BitWriter() : bit_ptr(0) {}

        void write_bits(uint32_t val, int n) {
            for (int i = 0; i < n; i++) {
                size_t bit_idx = 7 - (bit_ptr % 8);
                size_t byte_idx = bit_ptr / 8;

                if (byte_idx >= (int)data.size()) {
                    data.push_back(0);
                }

                int bit = (val >> (n - 1 - i)) & 1;
                if (bit) {
                    data[byte_idx] |= 1 << bit_idx;
                }
                bit_ptr++;
            }
        }

        void write_bit(bool bit) {
            write_bits(bit ? 1 : 0, 1);
        }

        void byte_align_zero() {
            while ((bit_ptr & 7) != 0) {
                write_bit(false);
            }
        }

        const std::vector<uint8_t>& flush() const {
            return data;
        }

        std::vector<uint8_t> into_bytes() {
            return std::move(data);
        }

        size_t bits_written() const {
            return bit_ptr;
        }

    private:
        std::vector<uint8_t> data;
        size_t bit_ptr = 0;
    };

    class BitReader {
    public:
        explicit BitReader(const std::vector<uint8_t>& bytes)
            : data(bytes.data()), len(bytes.size()), bit_ptr(0) {
        }

        uint32_t read_bits(int n) {
            uint32_t val = 0;
            for (int i = 0; i < n; i++) {
                size_t byte_idx = bit_ptr / 8;
                size_t bit_idx = 7 - (bit_ptr % 8);
                int bit = (data[byte_idx] >> bit_idx) & 1;
                val = (val << 1) | bit;
                bit_ptr++;
            }
            return val;
        }

        size_t bit_pos() const { return bit_ptr; }

    private:
        const uint8_t* data;
        size_t len;
        size_t bit_ptr;
    };

}