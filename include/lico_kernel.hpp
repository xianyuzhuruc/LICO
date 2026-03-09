// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <climits>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <fitting-method/optimal_pla.hpp>
#include <fitting-method/swing_pla.hpp>
#include <config.hpp>
#if RESIDUAL_COMPRESS
#include <lico_fastpfor.hpp>
#endif


namespace lico
{
    #define BIT_CEIL(x) ((x) < 2 ? 1u : 1u << (64u - __builtin_clzll((x) - 1)))
    #define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

    template <typename K>
    class LICO {

    public:
        typedef int64_t Simd_Value;
        typedef int32_t Correction_Value;
        typedef int32_t Intercept_Value;
        typedef uint32_t Covered_Value;

        typedef int32_t Segment_Value;

        uint64_t n;                           ///< The number of elements this index was built on.
        uint32_t Epsilon_Data;                ///< The epsilon value used to build the index.actual it's uint32_t
        uint8_t bpc_epsilon;

        std::vector<Segment_Value> seg_s_i;

        std::vector<Segment_Value> seg_delta_x;
        std::vector<Segment_Value> seg_delta_y;
        std::vector<Segment_Value> seg_y_b;
        std::vector<Segment_Value> seg_x_b;
        std::vector<Covered_Value> seg_covered;


        sdsl::int_vector<64> seg_covered_compress;
        sdsl::int_vector<64> seg_delta_y_compress;
        sdsl::int_vector<64> seg_delta_x_compress;
        sdsl::int_vector<64> seg_y_b_compress;
        sdsl::int_vector<64> seg_x_b_compress;


        sdsl::int_vector<64> corrections_none; // corrections for saving, each value <= (epsilon / 4)
        sdsl::bit_vector signs_none; // signs compress for saving
        std::vector<uint32_t> corrections_compress;

        std::vector<Correction_Value> corrections_vector; // corrections for decode

        Segment_Value first_offset = 0;

        uint32_t segments_size;
        uint8_t bpc_covered, bpc_delta_y, bpc_delta_x, bpc_y_b, bpc_x_b;

        /// Sentinel value to avoid bounds checking.
        static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity() : std::numeric_limits<K>::max();

        using position_type = typename std::conditional_t<sizeof(K) <= 4, uint32_t, uint64_t>;

#ifdef OptimalPLA
        using canonical_segment = typename internal_optimal::OptimalPiecewiseLinearModel<position_type, K>::CanonicalSegment;
#elifdef SwingFilterPLA
        using canonical_segment = typename internal_swing::SwingPiecewiseLinearModel<position_type, K>::CanonicalSegment;
#endif

        uint64_t errorPointCount;

        LICO(): n(0), errorPointCount(0) {}
        LICO(uint32_t epsilon) {
            bpc_epsilon = BIT_WIDTH(epsilon);
            Epsilon_Data = epsilon;
        }

        explicit LICO(const std::vector<K>& data, uint32_t epsilon) : LICO(data.begin(), data.end(), epsilon) {}

        template <typename RandomIt>
        LICO(RandomIt begin, RandomIt end, uint32_t epsilon, Segment_Value first_offset = 0):
            first_offset(first_offset),
            n(std::distance(begin, end)),
            errorPointCount(0){
            Epsilon_Data = epsilon;
            bpc_epsilon = BIT_WIDTH(epsilon);
            corrections_vector.resize(n);

#if RESIDUAL_COMPRESS
#else
            corrections_none.resize(((n * bpc_epsilon + 63) >> 6)); // set 3 to avoid insEufficient memory
            signs_none.resize(n);
#endif

            build(begin, end, epsilon);

#if RESIDUAL_COMPRESS
            if (residual_compress_type == "fastpfor")
                corrections_compress = compress_residuals_fastpfor(corrections_vector);
#endif

            segments_compress();
        }

        template <typename RandomIt>
        void build(RandomIt begin, RandomIt end, uint64_t epsilon) {

            auto n = (uint64_t) std::distance(begin, end);
            if (n == 0) return;

            if (*std::prev(--end) == sentinel)
                throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

            std::vector<canonical_segment> canonical_segments;
            canonical_segments.reserve(epsilon > 0 ? n / (epsilon * epsilon) : n / 8);

            auto in_fun = [begin](auto i) { return std::pair<position_type, K>(i, begin[i]); };
            auto out_fun = [&canonical_segments](auto cs) { canonical_segments.push_back(cs); };

#ifdef OptimalPLA
            auto n_segments = internal_optimal::make_segmentation_par(n, epsilon, in_fun, out_fun);
#elifdef SwingFilterPLA
            auto n_segments = internal_swing::make_segmentation_par(n, epsilon, in_fun, out_fun);
#endif

            for (auto it = canonical_segments.begin(); it < canonical_segments.end(); ++it) {
                auto i = it -> get_first_x();
                auto j = std::next(it) != canonical_segments.end() ? std::next(it) -> get_first_x() : n;
                build_segments(*it, n, i, j, begin); // build the segment
            }
        }

        template <typename RandomIt>
        void build_segments(const canonical_segment& cs, uint64_t n, uint64_t i, uint64_t j, RandomIt data) {

            uint32_t first = cs.get_first_x();
            if (first == n) return;

            auto [delta_y, delta_x, y_b, x_b] = cs.get_lico_segment();
            x_b += first_offset; // important to add offset
            if (first == n - 1) {
                seg_covered.push_back(1);
                seg_delta_x.push_back(1);
                seg_delta_y.push_back(0);
                seg_x_b.push_back(first + first_offset);
                seg_y_b.push_back(data[first]);


#if RESIDUAL_COMPRESS
                corrections_vector[first] = 0;
#else
                set_correction(corrections_none.data(), first, 0, 0, signs_none);
#endif
                return;
            }


            seg_covered.push_back(j - i);
            seg_s_i.push_back(first);
            seg_delta_x.push_back(delta_x);
            seg_delta_y.push_back(delta_y);
            seg_x_b.push_back(x_b);
            seg_y_b.push_back(y_b);

//
#if RESIDUAL_COMPRESS
            int32_t last_correction = 0;
            for (Covered_Value p = first; p < j; p++) {
                // int64_t error = static_cast<int64_t> (data[p]) - seg_approximate(p, first, cs_exponent, cs_significand, cs_intercept);
                int64_t error = static_cast<int64_t> (data[p]) -  seg_approximate(delta_y, delta_x,  y_b, x_b, p + first_offset);
                int32_t error_diff = error - last_correction;
                corrections_vector[p] = error_diff;
                last_correction = error;
            }
#else
            uint64_t* corrections_ptr = corrections_none.data();
            sdsl::bit_vector& signs_ptr = signs_none;

            for (Covered_Value p = first; p < j; p++) {
                int64_t error = static_cast<int64_t> (data[p]) -  seg_approximate(delta_y, delta_x,  y_b, x_b, p + first_offset);
                uint8_t sign_value = error >= 0 ? 0 : 1;
                error = std::abs(error);
                if (error <= Epsilon_Data)
                    set_correction(corrections_ptr, p, error, sign_value, signs_ptr);
                else
                    std::cerr << "Error: error epsilon: " << error << " " << Epsilon_Data << " " << (sign_value == 0 ? 1 : -1) << " " << p << " " << data[p]  << std::endl;
            }
#endif
        }

        int64_t seg_approximate(Segment_Value delta_y, Segment_Value delta_x, Segment_Value y_b, Segment_Value x_b, Segment_Value j) {
            int64_t result = static_cast<int64_t> (delta_y) * (j - x_b) + (delta_x >> 1) * ((j > x_b) - (j < x_b));
            result /= delta_x;
            result += y_b;
            return result;
        }

        void segments_compress() {
            segments_size = seg_covered.size();
            auto max_iter_1 = std::max_element(seg_covered.begin(), seg_covered.end());
            bpc_covered = BIT_WIDTH(*max_iter_1);
            auto max_iter_2 = std::max_element(seg_delta_y.begin(), seg_delta_y.end());
            bpc_delta_y = BIT_WIDTH(*max_iter_2);
            auto max_iter_3 = std::max_element(seg_delta_x.begin(), seg_delta_x.end());
            bpc_delta_x = BIT_WIDTH(*max_iter_3);
            auto max_iter_4 = std::max_element(seg_y_b.begin(), seg_y_b.end());
            bpc_y_b = BIT_WIDTH(*max_iter_4);
            auto max_iter_5 = std::max_element(seg_x_b.begin(), seg_x_b.end());
            bpc_x_b = BIT_WIDTH(*max_iter_5);


            seg_covered_compress.resize((segments_size * bpc_covered + 63) >> 6);
            seg_delta_y_compress.resize((segments_size * bpc_delta_y + 63) >> 6);
            seg_delta_x_compress.resize((segments_size * bpc_delta_x + 63) >> 6);
            seg_y_b_compress.resize((segments_size * bpc_y_b + 63) >> 6);
            seg_x_b_compress.resize((segments_size * bpc_x_b + 63) >> 6);

            // for (int i = 1; i < seg_covered.size(); i++) {
            //     if (seg_intercept[i] <= seg_intercept[i - 1]) {
            //         std::cerr << "Error: intercept not increasing: " << seg_intercept[i] << " " << seg_intercept[i - 1] << std::endl;
            //     }
            // }

            for (int i = 0; i < seg_covered.size(); i++) {
                auto covered = seg_covered[i];
                set_segment_compress(seg_covered_compress.data(), i, covered, bpc_covered);
                auto delta_y = seg_delta_y[i];
                set_segment_compress(seg_delta_y_compress.data(), i, delta_y, bpc_delta_y);
                auto delta_x = seg_delta_x[i];
                set_segment_compress(seg_delta_x_compress.data(), i, delta_x, bpc_delta_x);
                auto y_b = seg_y_b[i];
                set_segment_compress(seg_y_b_compress.data(), i, y_b, bpc_y_b);
                auto x_b = seg_x_b[i];
                set_segment_compress(seg_x_b_compress.data(), i, x_b, bpc_x_b);
            }
            // clear the original segments
            seg_covered.clear();
            seg_delta_y.clear();
            seg_delta_x.clear();
            seg_y_b.clear();
            seg_x_b.clear();

            corrections_vector.clear();
            std::vector<Covered_Value>().swap(seg_covered);
            std::vector<Segment_Value>().swap(seg_delta_y);
            std::vector<Segment_Value>().swap(seg_delta_x);
            std::vector<Segment_Value>().swap(seg_y_b);
            std::vector<Segment_Value>().swap(seg_x_b);
            std::vector<Correction_Value>().swap(corrections_vector);
        }


        void free_memory() {
            std::vector<Covered_Value>().swap(seg_covered);
            std::vector<Segment_Value>().swap(seg_delta_y);
            std::vector<Segment_Value>().swap(seg_delta_x);
            std::vector<Segment_Value>().swap(seg_y_b);
            std::vector<Segment_Value>().swap(seg_x_b);

            std::vector<Correction_Value>().swap(corrections_vector);
        }



        // int64_t seg_approximate(uint32_t i, uint32_t first, uint8_t exponent, int64_t significand, int32_t intercept) const {
        //     return (int64_t(significand * (i - first)) >> exponent) + intercept;
        // }

        uint64_t segment_slope_significand_max() const {
            auto max_iter = std::max_element(seg_y_b.begin(), seg_y_b.end());
            uint64_t max_slope_significand = *max_iter;
            return max_slope_significand;
        }

        uint32_t segment_slope_exponent_max() const {
            auto max_iter = std::max_element(seg_delta_y.begin(), seg_delta_y.end());
            uint32_t max_slope_exponent = *max_iter;
            return max_slope_exponent;
        }

        uint64_t size() const {
            return n;
        }
        

        uint64_t total_size_in_bytes() const {
            return segment_size_in_bytes() + corrections_size_in_bytes() + signs_size_in_bytes();
        }

        uint64_t ground_truth_build_size_in_bytes() const {
            return (double(segments_size) * 20.0)  + sizeof(uint8_t) + corrections_size_in_bytes() + signs_size_in_bytes();
        }

        uint64_t segment_size_in_bytes() const {
            return ((seg_covered_compress.bit_size() + seg_delta_y_compress.bit_size() + seg_delta_x_compress.bit_size() + seg_y_b_compress.bit_size() + seg_x_b_compress.bit_size()) / CHAR_BIT) + sizeof(uint8_t); //here uint8_t represents the byte need for Epsilon_Data, actually it's only need 8 bit, but we store it as u32 for convenient
        }

        uint64_t corrections_size_in_bytes() const {
#if RESIDUAL_COMPRESS
            return corrections_compress.size() * 4;
#else
            return (corrections_none.bit_size()) / CHAR_BIT;
#endif
        }

        uint64_t signs_size_in_bytes() const {
#if RESIDUAL_COMPRESS
            return 0;
#else
            return (signs_none.bit_size()) / CHAR_BIT;
#endif
        }


        void report_residual_list(std::ofstream &file) {
#if RESIDUAL_COMPRESS
            file << n << std::endl;
            uint32_t first = 0;
            for (uint32_t i = 0; i < segments_size; i++) {
                auto covered = seg_covered[i];
                int32_t last_correction = 0;
                // file << covered << std::endl;
                for (int j = first; j < covered + first; j++) {
                    auto delta_residual = corrections_vector[j];
                    auto true_residual = last_correction + delta_residual;
                    last_correction = true_residual;
                    file << true_residual << " " << delta_residual << " " << zigzag(delta_residual) << std::endl;
                }
                first += covered;
            }
#else
            throw std::runtime_error("don't support to report uncompress residuals");
#endif
        }

        void segment_init() {
            seg_covered.resize(segments_size);
            seg_delta_y.resize(segments_size);
            seg_delta_x.resize(segments_size);
            seg_y_b.resize(segments_size);
            seg_x_b.resize(segments_size);


            for (int i = 0; i < segments_size; i++) {
                seg_covered[i] = get_segment_compress(seg_covered_compress.data(), i, bpc_covered);
                seg_delta_y[i] = get_segment_compress(seg_delta_y_compress.data(), i, bpc_delta_y);
                seg_delta_x[i] = get_segment_compress(seg_delta_x_compress.data(), i, bpc_delta_x);
                seg_y_b[i] = get_segment_compress(seg_y_b_compress.data(), i, bpc_y_b);
                seg_x_b[i] = get_segment_compress(seg_x_b_compress.data(), i, bpc_x_b);
            }
        }

        void residual_init() {
#if RESIDUAL_COMPRESS
            if (residual_compress_type == "fastpfor") {
                std::vector<uint32_t> uncompressed_output(n);
                corrections_vector.resize(n);
                decompress_residuals_fastpfor(corrections_compress, n, uncompressed_output, corrections_vector);
            } else {
                throw std::runtime_error("Not support uncompress residuals");
            }
#else
            corrections_vector.resize(n);
            int32_t * corrections_pointer = corrections_vector.data();
            uint32_t first = 0;
            for(uint32_t i = 0; i < segments_size; i++) {
                auto covered = seg_covered[i];
                int32_t last_correction = 0;
                for (Covered_Value j = first; j < first + covered; ++j) {
                    int32_t correction_varify = get_correction(corrections_none.data(), j, signs_none);
#ifndef __CUDACC__
                    *corrections_pointer++ = correction_varify - last_correction;
                    last_correction = correction_varify;
#else
                    std::cerr << "None" << std::endl;
                    *corrections_pointer++ = correction_varify;
#endif
                }
                first += covered;
            }
#endif
        }

        void normal_init(){
            segment_init();
            residual_init();
        }

        void normal_clean() {
            seg_covered.clear();
            seg_delta_x.clear();
            seg_y_b.clear();
            seg_delta_y.clear();
            corrections_vector.clear();
            std::vector<Covered_Value>().swap(seg_covered);
            std::vector<Segment_Value>().swap(seg_delta_y);
            std::vector<Segment_Value>().swap(seg_delta_x);
            std::vector<Segment_Value>().swap(seg_y_b);
            std::vector<Segment_Value>().swap(seg_x_b);
            std::vector<Correction_Value>().swap(corrections_vector);
        }

        std::vector<K> scalar_decode() {
            std::vector<K> output;
            output.resize(n);
            uint32_t pointer = 0;

            residual_init();

            int32_t last_first = 0;
            for (uint32_t i= 0; i < segments_size; i++) {
                const auto covered = seg_covered[i];
                const auto delta_y = seg_delta_y[i];
                const auto delta_x = seg_delta_x[i];
                const auto delta_x_divide = delta_x >> 1;
                const auto y_b = seg_y_b[i];
                const auto x_b = seg_x_b[i];

                int32_t last_correction = y_b +  corrections_vector[pointer];

                int32_t j = last_first + first_offset;
                int64_t result = (static_cast<int64_t> (delta_y) * (j - x_b) + (delta_x_divide) * (j > x_b ? 1 : j == x_b ? 0 : -1));
                output[pointer++] = result / delta_x + last_correction;
                j++;
                while (j < last_first + covered + first_offset) {
                    last_correction = last_correction + corrections_vector[pointer];

                    result += delta_y + (delta_x_divide) * (j == x_b ? 1 : j == x_b + 1 ? 1 : 0) ;

                    output[pointer++] = result / delta_x + last_correction;
                    j++;
                }
                last_first += covered;
            }

            return output;
        }

        uint64_t get_correction_bit_offset(uint64_t i, uint32_t bpc) {
            return i * bpc;
        }

        void set_segment_compress(uint64_t* compress_val, uint64_t i, uint64_t value, uint8_t bpc) {
                uint64_t idx = get_correction_bit_offset(i, bpc);
                sdsl::bits::write_int(compress_val + (idx >> 6u), value, idx & 0x3f, bpc);
            }

        int64_t get_segment_compress(const uint64_t* compress_val, uint64_t i, uint8_t bpc) {
                uint64_t idx = get_correction_bit_offset(i, bpc);
                uint64_t val = sdsl::bits::read_int(compress_val + (idx >> 6u), idx & 0x3f, bpc);
                return val;
            }


        void set_correction(uint64_t* corrections_compress_val, uint64_t i, uint64_t value, uint8_t sign_value, sdsl::bit_vector& signs_ptr) {
            uint64_t idx = get_correction_bit_offset(i, bpc_epsilon);
            sdsl::bits::write_int(corrections_compress_val + (idx >> 6u), value, idx & 0x3f, bpc_epsilon);
            signs_ptr[i] = sign_value & 0x01;
        }


        int64_t get_correction(const uint64_t* corrections_compress_val, int64_t i, const sdsl::bit_vector& signs_ptr)  {
            uint64_t idx = get_correction_bit_offset(i, bpc_epsilon);
            uint64_t correction = sdsl::bits::read_int(corrections_compress_val + (idx >> 6u), idx & 0x3f, bpc_epsilon);
            if (correction == 0 && signs_ptr[i] == 1)
                return Epsilon_Data + 1;
            return signs_ptr[i] == 0 ? correction : -correction;
        }

     };

}
