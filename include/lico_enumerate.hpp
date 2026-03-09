#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <config.hpp>
#include <lico_kernel.hpp>
#include <tools.hpp>
#include <libdivide.h>

#if RESIDUAL_COMPRESS
#include <lico_fastpfor.hpp>
#endif




namespace lico_sequence {
    using LICOIndex = lico::LICO<uint32_t>;
    typedef int32_t Segment_Value;
    typedef int64_t Simd_Value;
    typedef uint64_t Slope_Value;
    typedef int32_t Correction_Value;


    struct segment{

        Segment_Value covered; // covered
        Segment_Value delta_y;
        Segment_Value delta_x;
        Segment_Value y_b; // 32 bits
        Segment_Value x_b;
        Segment_Value first;

        inline segment(Segment_Value first, Segment_Value delta_y, Segment_Value delta_x, Segment_Value y_b, Segment_Value x_b, Segment_Value covered) :
                first(first), delta_y(delta_y), delta_x(delta_x), y_b(y_b), x_b(x_b), covered(covered) {}
    };

#if USE_HUGEPAGE
    template <typename T>
    class HugePageAllocator {
        constexpr static size_t PAGE_SIZE = 2 * 1024 * 1024; // 2MB
    public:
        using value_type = T;

        HugePageAllocator() = default;
        ~HugePageAllocator() = default;

        static T* allocate(std::size_t n) {
            const size_t size = (n * sizeof(T) + PAGE_SIZE - 1) / PAGE_SIZE * PAGE_SIZE;
            // const size_t size = n * sizeof(T);
            void* ptr = mmap(nullptr, size,
                             PROT_WRITE | PROT_READ,
                             MAP_SHARED | MAP_ANONYMOUS | MAP_HUGETLB,
                             -1, 0);
            if (ptr == MAP_FAILED) {
                std::cerr << "Failed to allocate huge pages: " << std::strerror(errno) << std::endl;
                throw std::bad_alloc();
            }
            return static_cast<T*>(ptr);
        }

        static void deallocate(T* p, std::size_t n) noexcept {
            const size_t size = (n * sizeof(T) + PAGE_SIZE - 1) / PAGE_SIZE * PAGE_SIZE;
            munmap(p, size);
        }

        bool operator == (const HugePageAllocator&) const { return true; }
        bool operator != (const HugePageAllocator&) const { return false; }
    };
#endif

    template <typename K>
    class lico_enumerator{

        public:

        uint64_t n;

        std::vector<segment> segments;

        std::vector<K> segment_first_values; // size == segments.size()

        std::vector<K> segment_max;

        std::vector<uint64_t> parted_sizes;

        std::vector<Correction_Value> corrections_vector; // corrections for decode

        std::vector<std::vector<uint32_t>> corrections_compress_fastpfor;

        std::vector<uint32_t> block_sizes;

        void load_block_size(std::vector<uint32_t> all_block_sizes) {
            this -> block_sizes = std::move(all_block_sizes);
        }

        void load_residuals(uint64_t data_size, std::vector<Correction_Value> corrections_vector) {
            this -> n = data_size;
            this -> corrections_vector = std::move(corrections_vector);
        }

        void load_residuals_fastpfor(uint64_t data_size, std::vector<std::vector<uint32_t>> compress_fastpfor, std::vector<uint64_t> parted_sizes) {
            this -> n = data_size;
            this -> corrections_compress_fastpfor = std::move(compress_fastpfor);
            this -> parted_sizes = std::move(parted_sizes);
        }


        // used for Query Test
        K current_value = INT_MAX;
        K next_first_value = INT_MAX;
        Correction_Value current_correction = 0;
        Segment_Value current_pos = 0;
        Segment_Value current_segment = 0;
        int32_t total_segment_size = 0;

#if USE_HUGEPAGE
        std::vector<K, HugePageAllocator<K>> current_value_vector;
#else
        std::vector<K> current_value_vector;
#endif


        void query_init(const std::string decode_type, const std::string query_type) {
            if (decode_type == "normal") {
                if (query_type == "intersection") {
                    current_pos = 0;
                    current_segment = 0;
                    total_skip = 0;
                } else if (query_type == "union") {
                    current_pos = 0;
                    current_segment = 0;
                    current_value = INT_MAX - 1;
                    total_skip = 0;
                }
            } else if (decode_type == "simd") {
                simd_init();
                std::vector<Correction_Value> ().swap(corrections_vector);
                std::vector<segment> ().swap(segments);
                current_pos = 0;
                current_segment = 0;
                current_value = INT_MAX - 1;
                total_skip = 0;
            }
        }

        void decode_query(K* output, const std::string decode_type) {
            if (decode_type == "normal") {
                normal_decode(output);
            } else if (decode_type == "simd") {
                simd_decode_512i_2d1s(output);
            }
        }

        long double total_skip = 0;


        K* next_pointer;

        void next_init() {
            simd_init();
            std::vector<Correction_Value> ().swap(corrections_vector);
            std::vector<segment> ().swap(segments);

            current_pos = 0;
            current_value_vector.resize(n + 1);

            simd_decode_512i_2d1s(current_value_vector.data());

            current_value_vector[n] = INT_MAX;
            next_pointer = current_value_vector.data();
            current_value = *next_pointer;
        }

        K docid() {
            return current_value;
        }

        void next() {
            current_value = *++next_pointer;
        }


        int64_t current_numerator;
        using branchfree_t = libdivide::branchfree_divider<int64_t>; // use libdivide to speed the division

        uint64_t residuals_resume_duration;

        void build_segment_first_values() { // pre-decode the first value of each segment
            segment_first_values.clear();
            segment_first_values.resize(segments.size() + 1);

            K* out = segment_first_values.data();
            size_t posi = 0;

            auto start = std::chrono::high_resolution_clock::now();
            for (const auto& s : segments) {
                const int32_t j_x_diff = s.first - s.x_b;
                const int64_t delta_x_add = (s.delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
                const int64_t numerator = static_cast<int64_t>(s.delta_y) * j_x_diff + delta_x_add;
        #if RESIDUAL_COMPRESS
                out[posi++] = static_cast<K>(numerator / s.delta_x + s.y_b + unzigzag_int32(corrections_vector[s.first]));
        #else
                out[posi++] = static_cast<K>(numerator / s.delta_x + s.y_b + corrections_vector[s.first]);
        #endif
            }
            out[posi] = INT_MAX; // sentinel

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
        }

        void resume_residuals() {
#if RESIDUAL_COMPRESS // Notably, LICO do not need resume residuals from gaps indeed, it's only necessary for LICO++.
            int32_t i = 0;
            __m512i v_correction;
            auto start = std::chrono::high_resolution_clock::now();
            for (; i + 16 < corrections_vector.size(); i += 16) {
                v_correction = _mm512_loadu_epi32(&corrections_vector[i]);
                v_correction = _mm512_xor_si512(_mm512_srli_epi32(v_correction, 1), _mm512_srai_epi32(_mm512_slli_epi32(v_correction, 31), 31)); // (u >> 1) ^ -(u & 1)
                _mm512_storeu_epi32(&corrections_vector[i], v_correction);
            }
            for (; i < corrections_vector.size(); i ++) {
                corrections_vector[i] = unzigzag_int32(corrections_vector[i]);
            }

            Correction_Value* __restrict correct_pointer = corrections_vector.data();
            for (const auto& s : segments) {
                for (int32_t j = 1; j < s.covered; j++) {
                    *(correct_pointer + 1) += *correct_pointer;
                    correct_pointer++;
                }
                correct_pointer++;
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
            residuals_resume_duration = total_duration;
#else // To simplify coding, we unify the sequential decoding pipline using the gaps of residuals, but this transformation is not necessary for LICO. So we exclude this cost in LICO.
            Correction_Value* __restrict correct_pointer = corrections_vector.data();
            for (const auto& s : segments) {
                for (int32_t j = 1; j < s.covered; j++) {
                    *(correct_pointer + 1) += *correct_pointer;
                    correct_pointer++;
                }
                correct_pointer++;
            }
#endif
        }

        void next_geq_init() {
            residuals_decode();
            build_segment_first_values();
            resume_residuals();

            total_segment_size = static_cast<int32_t>(segments.size());

            current_pos = 0;
            current_segment = 0;
            const segment* seg = &segments[0];
            current_correction = seg -> y_b;

            const int32_t j_x_diff = seg -> first - seg -> x_b;
            const int64_t delta_x_add = (seg->delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
            current_numerator = static_cast<int64_t>(seg->delta_y) * j_x_diff + delta_x_add;

            current_value = segment_first_values[0];
            next_first_value = segment_first_values[1];
        }


        void random_geq_init() {
            residuals_decode();
            build_segment_first_values();
            resume_residuals();

            total_segment_size = static_cast<int32_t>(segments.size());

            current_segment = 0;
        }

        std::vector<K> segment_first_positions;

        void build_segment_first_positions() { // pre-decode the first value of each segment
            segment_first_positions.clear();
            segment_first_positions.resize(segments.size() + 1);

            K* out = segment_first_positions.data();

            auto start = std::chrono::high_resolution_clock::now();
            for (const auto& s : segments) {
                *out++ = s.first;
            }

            *out = INT_MAX;
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
        }

        void random_access_init() {
            residuals_decode();
            resume_residuals();
            build_segment_first_positions();
        }

        K random_access(K position) {
            // const auto seg = std::lower_bound(segments.data(), segments.data() + segments.size(), position,
                // [](const segment& seg, const K& pos) {return seg.first < pos;});

            const auto it = std::upper_bound(segment_first_positions.begin(), segment_first_positions.end(), position) - segment_first_positions.begin() - 1;

            const segment* seg = &segments[it];

            // assert(it < segment_first_positions.size());
            const int32_t j_x_diff = position - seg -> x_b;
            const int64_t numerator = seg -> delta_y * j_x_diff + (seg -> delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
            return static_cast<K>(numerator / seg -> delta_x + seg -> y_b + corrections_vector[position]);

        }


        K random_geq(K posting_value) {
            const K* sfv = segment_first_values.data();
            const K* it = std::upper_bound(sfv, sfv + total_segment_size, posting_value);
            current_segment = std::max(static_cast<int32_t>(it - sfv) - 1, 0);

            // 获取当前段
            const segment* seg = &segments[current_segment];
            const int32_t j_start = seg -> first;
            const int32_t j_end = j_start + seg -> covered;

            const int64_t dy = seg -> delta_y;
            const int64_t dx_div = seg -> delta_x;
            const int64_t dx_half64 = seg -> delta_x >> 1;
            const int32_t xb = seg -> x_b;


            auto decode_value = [&](const int32_t j) -> K {
                const int32_t j_x_diff = j - xb;
                const int64_t numerator = dy * j_x_diff + dx_half64 * ((j_x_diff > 0) - (j_x_diff < 0));
                const Correction_Value last_correction = seg -> y_b + corrections_vector[j];
                return static_cast<K>(numerator / dx_div + last_correction);
            };

            K result;

            // binary_search
            int32_t left = j_start, right = j_end;

            while (left < right) {
                int32_t mid = (left + right) >> 1;
                const K mid_value = decode_value(mid);
                if (mid_value < posting_value) {
                    left = mid + 1;
                } else {
                    result = mid_value;
                    right = mid;
                }
            }

            if (left >= j_end) {
                return segment_first_values[current_segment + 1];
            }

            return result;
        }

        inline K nextgeq(const K posting_value) {
            if (current_value >= posting_value) {
                return current_value;
            }

            const K* sfv = segment_first_values.data();

            if (posting_value >= next_first_value) {
                const int32_t seg_idx = current_segment;
                int32_t i = seg_idx + 1;
                const int32_t limit = std::min<int32_t>(seg_idx + 1 + 8, total_segment_size);

                while (i < limit && posting_value >= sfv[i]) {
                    ++i;
                }

                if (i < total_segment_size && posting_value < sfv[i]) {
                    current_segment = i - 1;
                } else {
                    const K* begin = sfv + i;
                    const K* end = sfv + total_segment_size;
                    const K* it = std::upper_bound(begin, end, posting_value);
                    current_segment = static_cast<int32_t>(it - sfv) - 1;
                }

                current_pos = segments[current_segment].first;
                next_first_value = sfv[current_segment + 1];
            }

            // binary search
            const segment* seg = &segments[current_segment];

            const int32_t j_end = seg -> first + seg -> covered;
            const int32_t xb = seg -> x_b;
            const int32_t yb = seg -> y_b;

            const int64_t dy = seg -> delta_y;
            const int64_t dx_div = seg -> delta_x;
            const int64_t dx_half64 = seg -> delta_x >> 1;

            int32_t left = current_pos, right = j_end, mid;

            K mid_value;
            int32_t j_x_diff;
            int64_t numerator;
            Correction_Value correction;

            // if (left > xb) {
            //     while (left < right) {
            //         mid = (left + right) >> 1;
            //
            //         numerator = dy * (mid - xb) + dx_half64;
            //         correction  = yb + corrections_vector[mid];
            //         mid_value = static_cast<K>(static_cast<Correction_Value> (numerator / dx_div) + correction);
            //
            //         if (mid_value < posting_value) {
            //             left = mid + 1;
            //         } else {
            //             current_value = mid_value;
            //             right = mid;
            //         }
            //     }
            // } else {
            while (left < right) {
                mid = (left + right) >> 1;

                j_x_diff = mid - xb;
                numerator = dy * j_x_diff + dx_half64 * ((j_x_diff > 0) - (j_x_diff < 0));
                correction = yb + corrections_vector[mid];
                mid_value = static_cast<K>(static_cast<Correction_Value> (numerator / dx_div) + correction);

                if (mid_value < posting_value) {
                    left = mid + 1;
                } else {
                    current_value = mid_value;
                    right = mid;
                }
            }
            // }

            if (left < j_end) {
                current_pos = left + 1; // next_point is left + 1
            } else {
                current_pos = j_end;
                current_value = next_first_value;
                ++current_segment;
                next_first_value = sfv[current_segment + 1];
            }

            return current_value;
        }

        K random_geq_naive(K posting_value) {
            const K* sfv = segment_first_values.data();
            const K* it = std::upper_bound(sfv, sfv + total_segment_size, posting_value);
            current_segment = std::max(static_cast<int32_t>(it - sfv) - 1, 0);

            if (current_segment >= total_segment_size) {
                return INT_MAX;
            }

            // internal segment scan
            const segment* seg = &segments[current_segment];
            int32_t j = seg -> first;
            const int32_t j_end   = j + seg -> covered;

            const int32_t xb = seg -> x_b;
            const int32_t xb_minus_1 = xb - 1;

            const int64_t dy = seg -> delta_y;
            const branchfree_t dx_div(seg -> delta_x);
            const int64_t dx_half64 = seg -> delta_x >> 1;

            const int32_t j_x_diff = j - seg -> x_b;
            int64_t numerator = dy * j_x_diff + dx_half64 * ((j_x_diff > 0) - (j_x_diff < 0));

            Correction_Value* corr_ptr = &corrections_vector[j];
            Correction_Value last_correction = seg -> y_b;

            {
                const int32_t p1_end = std::min<int32_t>(j_end, xb_minus_1);
                while (j < p1_end) {
#if RESIDUAL_COMPRESS
                    last_correction += unzigzag_int32(*corr_ptr++);
#else
                    last_correction += *corr_ptr++;
#endif
                    const K v = static_cast<K>(numerator / dx_div + last_correction);
                    numerator += dy;
                    if (v >= posting_value) {
                        return v;
                    }
                    ++j;
                }
            }

            // special process xb and xb - 1. Notably, here j is j + 1 for numerator, so we need check xb and xb - 1
            if ((j == xb_minus_1 || j == xb)) {
                for (int k = 0; k < 2  && (j == xb_minus_1 || j == xb); ++k) {
#if RESIDUAL_COMPRESS
                    last_correction += unzigzag_int32(*corr_ptr++);
#else
                    last_correction += *corr_ptr++;
#endif
                    const K v = static_cast<K>(numerator / dx_div + last_correction);
                    numerator += dy + dx_half64; // 两个特殊点都要额外加 dx_half
                    if (v >= posting_value) {
                        return v;
                    }
                    ++j;
                }
            }

            // last part
            while (j < j_end) {
#if RESIDUAL_COMPRESS
                last_correction += unzigzag_int32(*corr_ptr++);
#else
                last_correction += *corr_ptr++;
#endif
                const K v = static_cast<K>(numerator / dx_div + last_correction);
                numerator += dy;
                if (v >= posting_value) {
                    return v;
                }
                ++j;
            }

            return segment_first_values[current_segment + 1];
        }

        inline K nextgeq_naive(K posting_value) {
            if (current_value >= posting_value) return current_value;

            const K* sfv = segment_first_values.data();

            // small sequential step
            if (posting_value >= next_first_value) {
                const int32_t seg_idx = current_segment;
                int32_t i = seg_idx + 1;
                const int32_t limit = std::min<int32_t>(seg_idx + 1 + 8, total_segment_size); // 8 times, if we still don't found the suitable segment, then use upper_bound

                while (i < limit && posting_value >= sfv[i]) {
                    ++i;
                }

                if (i < total_segment_size && posting_value < sfv[i]) {
                    current_segment = i - 1;
                } else {
                    const K* begin = sfv + i;
                    const K* end   = sfv + total_segment_size;
                    const K* it    = std::upper_bound(begin, end, posting_value);
                    current_segment = static_cast<int32_t>(it - sfv) - 1;

                    if (current_segment >= total_segment_size) {
                        current_pos = n;
                        current_value = INT_MAX;
                        return current_value;
                    }
                }

                const segment* seg = &segments[current_segment];
                current_pos = seg -> first;
                current_correction = seg -> y_b;
                const int32_t j_x_diff = seg -> first - seg -> x_b;
                const int64_t delta_x_add = (seg->delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
                current_numerator = static_cast<int64_t>(seg->delta_y) * j_x_diff + delta_x_add;

                next_first_value = sfv[current_segment + 1];
            }

            // internal segment scan
            while (true) {
                const segment* seg = &segments[current_segment];
                const int32_t j_end   = seg -> first + seg -> covered;

                if (current_pos >= j_end) {
                    ++current_segment;
                    if (current_segment >= total_segment_size) {
                        current_pos = n;
                        // current_value = INT_MAX;
                        return INT_MAX;
                    }

                    const segment* nseg = &segments[current_segment];
                    current_pos = nseg -> first;
                    current_correction = nseg -> y_b;
                    const int32_t j_x_diff = nseg -> first - nseg -> x_b;
                    const int64_t delta_x_add = (nseg -> delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
                    current_numerator = static_cast<int64_t>(nseg -> delta_y) * j_x_diff + delta_x_add;
                    next_first_value = sfv[current_segment + 1];
                    continue;
                }

                const int32_t xb = seg -> x_b;
                const int32_t xb_minus_1 = xb - 1;
                const int64_t dy = seg -> delta_y;
                const int64_t dx_half64 = static_cast<int64_t>(seg -> delta_x >> 1);
                // const int64_t dx_div = seg -> delta_x;
                const branchfree_t dx_div(seg -> delta_x);

                const Correction_Value* corr_ptr = &corrections_vector[current_pos];
                const Correction_Value* corr_end = &corrections_vector[j_end];

                int32_t j = current_pos;

                {
                    const int32_t p1_end = std::min<int32_t>(j_end, xb_minus_1);
                    while (j < p1_end) {
        #if RESIDUAL_COMPRESS
                        current_correction += unzigzag_int32(*corr_ptr++);
        #else
                        current_correction += *corr_ptr++;
        #endif
                        const K v = static_cast<K>(current_numerator / dx_div + current_correction);
                        current_numerator += dy;
                        if (v >= posting_value) {
                            current_value = v;
                            current_pos = ++j;
                            return current_value;
                        }
                        ++j;
                    }
                }

                // special process xb and xb - 1. Notably, here j is j + 1 for numerator, so we need check xb and xb - 1
                if ((j == xb_minus_1 || j == xb)) {
                    // if ((j == xb_minus_1 || j == xb) && j < j_end && corr_ptr < corr_end) {
                    for (int k = 0; k < 2  && (j == xb_minus_1 || j == xb); ++k) {
        #if RESIDUAL_COMPRESS
                        current_correction += unzigzag_int32(*corr_ptr++);
        #else
                        current_correction += *corr_ptr++;
        #endif
                        const K v = static_cast<K>(current_numerator / dx_div + current_correction);
                        current_numerator += dy + dx_half64; // 两个特殊点都要额外加 dx_half
                        if (v >= posting_value) {
                            current_value = v;
                            current_pos = ++j;
                            return current_value;
                        }
                        ++j;
                    }
                }

                // last part
                while (j < j_end) {
        #if RESIDUAL_COMPRESS
                    current_correction += unzigzag_int32(*corr_ptr++);
        #else
                    current_correction += *corr_ptr++;
        #endif
                    const K v = static_cast<K>(current_numerator / dx_div + current_correction);
                    current_numerator += dy;
                    if (v >= posting_value) {
                        current_value = v;
                        current_pos = ++j;
                        return current_value;
                    }
                    ++j;
                }

                // go to next segment
                ++current_segment;
                if (current_segment >= total_segment_size) {
                    current_pos = n;
                    return INT_MAX;
                }

                const segment* nseg = &segments[current_segment];
                current_pos = nseg->first;
                current_correction = nseg->y_b;
                {
                    const int32_t j_x_diff2 = nseg->first - nseg->x_b;
                    const int64_t delta_x_add2 = (nseg->delta_x >> 1) * ((j_x_diff2 > 0) - (j_x_diff2 < 0));
                    current_numerator = static_cast<int64_t>(nseg->delta_y) * j_x_diff2 + delta_x_add2;
                }
                next_first_value = sfv[current_segment + 1];
            }
        }

        void residuals_decode() {
#if RESIDUAL_COMPRESS
            total_duration = 0;
            if (residual_compress_type == "fastpfor") {
                std::vector<uint32_t> uncompressed_output;
                uncompressed_output.resize(n);

                for (int32_t i = 0; i < 3; i++)
                    for (int32_t j = 0; j < n; j++)
                        uncompressed_output[j] = j + 1 + i; // just for warm up

                auto start = std::chrono::high_resolution_clock::now();

                decompress_residuals_fastpfor_parted(corrections_compress_fastpfor, uncompressed_output, parted_sizes);

                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                total_duration = duration.count();

                // actually, we can get int32 vector directly by the fastpfor decoding, beacuse 0 <= zigzag(residuals) <= 2*epsilon, and our epsilons are all smaller than (INT_MAX / 2), but dut to the FastPFor library only support ouput uint32_t type, so we cast them into int32_t
                corrections_vector.resize(n);
                for (int32_t i = 0; i < n; i++) {
                    corrections_vector[i] = static_cast<Correction_Value> (uncompressed_output[i]);
                }
            }
#endif
            // uncompress residuals is not need to be decoded
        }

        // the naive scale decode algorithm
        void normal_decode_naive(K* __restrict output) {
            Correction_Value* __restrict correct_pointer = corrections_vector.data();
            const auto end_iter = segments.end();

            auto start = std::chrono::high_resolution_clock::now();
            for (auto it = segments.begin(); it != end_iter; ++it) {
                const auto& seg = *it;
                const auto covered = seg.covered;
                const int64_t delta_y = seg.delta_y;
                const int64_t delta_x = seg.delta_x;
                const auto delta_x_divide = delta_x >> 1;
                const auto y_b = seg.y_b;
                const auto x_b = seg.x_b;

                int32_t last_correction = y_b +  (*correct_pointer++);

                int32_t j = seg.first;
                const int32_t j_end = j + covered;
                int64_t result = (static_cast<int64_t> (delta_y) * (j - x_b) + (delta_x_divide) * ((j > x_b) - (j < x_b)));
                *output++ = result / delta_x + last_correction;
                j++;
                while (j < j_end) {
                    last_correction = last_correction + (*correct_pointer++);

                    result += delta_y + (delta_x_divide) * ((j == x_b) + ( j == x_b + 1)) ;

                    *output++ = result / delta_x + last_correction;
                    j++;
                }
            }

            auto end = std::chrono::high_resolution_clock::now();
            total_duration += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }

        // the scale decode algorithm
        void normal_decode(K* __restrict output) {

            Correction_Value* __restrict correct_pointer = corrections_vector.data();
            const auto end_iter = segments.end();

            auto start = std::chrono::high_resolution_clock::now();
            for (auto it = segments.begin(); it != end_iter; ++it) {
                const auto& seg = *it;

                const int32_t covered = seg.covered;
                if (__builtin_expect(covered <= 0, 0)) { continue; }

                const int64_t dy = seg.delta_y;
                const int64_t dx = seg.delta_x; // 假定 dx > 0
                const int64_t half = (dx >> 1); // floor(dx/2)
                const int32_t xb = seg.x_b;

                int32_t last_correction = seg.y_b;
#if RESIDUAL_COMPRESS
                last_correction += (unzigzag_int32(*correct_pointer++));
#else
                last_correction += (*correct_pointer++);
#endif

                int32_t j = seg.first;
                // K* __restrict out = output + j;
                const int32_t j_end = j + covered;

                // numerator = dy * (j - xb) + half * sign(j - xb)
                int64_t num = dy * (int64_t(j) - int64_t(xb)) + half * ((int64_t)(j > xb) - (int64_t)(j < xb));

                int32_t q = num / dx;
                int32_t r = num - q * dx;

                *output++ = static_cast<K>(q + last_correction);
                ++j;

                if (__builtin_expect(j >= j_end, 0)) { continue; }

                const int32_t step_q = dy / dx;
                const int32_t step_r = dy - step_q * dx;

                // first part: j < xb, without half, the numerator is negative
                const int32_t e1 = (xb > j) ? ((j_end < xb) ? j_end : xb) : j;
                while (j < e1) {
#if RESIDUAL_COMPRESS
                    last_correction += (unzigzag_int32(*correct_pointer++));
#else
                    last_correction += (*correct_pointer++);
#endif
                    q += step_q;
                    r += step_r;

                    if (r > 0) {
                        r -= dx;
                        ++q;
                    }

                    *output++ = static_cast<K>(q + last_correction);
                    ++j;
                }

                if (j < j_end && j == xb) {
#if RESIDUAL_COMPRESS
                    last_correction += (unzigzag_int32(*correct_pointer++));
#else
                    last_correction += (*correct_pointer++);
#endif
                    *output++ = static_cast<K>(last_correction);
                    ++j;
                }

                // one more half step, if j = xb + 1. Due to j is added sequentially, so we do not need to double check where j is equal to xb + 1
                if (j < j_end) {
#if RESIDUAL_COMPRESS
                    last_correction += (unzigzag_int32(*correct_pointer++));
#else
                    last_correction += (*correct_pointer++);
#endif
                    q = step_q;
                    r = step_r + half;

                    if (__builtin_expect(r >= dx, 0)) {
                        r -= dx;
                        ++q;
                    }

                    *output++ = static_cast<K>(q + last_correction);
                    ++j;
                }

                // last part, without half, the numerator is positive
                while (j < j_end) {
#if RESIDUAL_COMPRESS
                    last_correction += (unzigzag_int32(*correct_pointer++));
#else
                    last_correction += (*correct_pointer++);
#endif
                    q += step_q;
                    r += step_r;

                    if (__builtin_expect(r >= dx, 0)) {
                        r -= dx;
                        ++q;
                    }

                    *output++ = static_cast<K>(q + last_correction);
                    ++j;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            total_duration += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }

        // used for SIMD
        std::vector<segment> segments_sort; // resorted segments
        std::vector<Segment_Value*> delta_y_simd;
        std::vector<Segment_Value*> delta_x_simd;
        std::vector<Segment_Value*> y_b_simd;
        std::vector<Segment_Value*> x_b_simd;
#if USE_HUGEPAGE
        std::vector<Correction_Value, HugePageAllocator<Correction_Value>> corrections_simd;
#else
        std::vector<Correction_Value> corrections_simd;
#endif
        std::vector<Correction_Value> corrections_tail;
        std::vector<Segment_Value*> first_simd;
        std::vector<Segment_Value*> list_covered_simd;
        std::vector<Segment_Value> group_min_covered_simd;

        uint64_t total_calculated = 0;
        uint64_t conversion_time = 0;
        uint64_t total_duration = 0;

        uint64_t idx = 0;

        template <typename T>
        T* aligned_new(uint64_t num_elements) {
            void* ptr = std::aligned_alloc(align_val, num_elements * sizeof(T));
            if (!ptr) throw std::bad_alloc();
            return static_cast<T*>(ptr);
        }

        template <typename T>
        static void aligned_delete(T* ptr) {
            std::free(ptr);
        }

#if USE_HUGEPAGE
        constexpr static size_t HUGE_PAGE_SIZE = 2 * 1024 * 1024;
        template <typename T>
        T* aligned_new_huge(uint64_t num_elements) {
            const size_t size = (num_elements * sizeof(T) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE;
            void* ptr = mmap(nullptr, size,
                             PROT_WRITE | PROT_READ,
                             MAP_SHARED | MAP_ANONYMOUS | MAP_HUGETLB,
                             -1, 0);
            if (ptr == MAP_FAILED) {
                std::cerr << "Failed to allocate huge pages: " << std::strerror(errno) << " " << size << " " << num_elements << std::endl;
                throw std::bad_alloc();
            }
            return static_cast<T*>(ptr);
        }

        template <typename T>
        static void aligned_delete_huge(T* p, std::size_t num_elements) noexcept {
            const size_t size = (num_elements * sizeof(T) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE;
            munmap(p, size);
        }
#endif

        void free_memory(std::string decode_type = "simd") {
            if (decode_type == "simd") {
                for (auto i = 0; i < delta_y_simd.size(); i++) {
                    aligned_delete(delta_y_simd[i]);
                    aligned_delete(delta_x_simd[i]);
                    aligned_delete(y_b_simd[i]);
                    aligned_delete(x_b_simd[i]);
                    aligned_delete(first_simd[i]);
                    aligned_delete(list_covered_simd[i]);
                }
                std::vector<Segment_Value> ().swap(group_min_covered_simd);
                std::vector<Correction_Value> ().swap(corrections_tail);
                std::vector<Segment_Value*> ().swap(delta_y_simd);
                std::vector<Segment_Value*> ().swap(delta_x_simd);
                std::vector<Segment_Value*> ().swap(y_b_simd);
                std::vector<Segment_Value*> ().swap(first_simd);
                std::vector<Segment_Value*> ().swap(list_covered_simd);
                std::vector<segment> ().swap(segments_sort);
#if USE_HUGEPAGE
                std::vector<K, HugePageAllocator<K>> ().swap(current_value_vector);
                std::vector<Correction_Value, HugePageAllocator<Correction_Value>> ().swap(corrections_simd);
#else
                std::vector<Correction_Value> ().swap(corrections_simd);
                std::vector<K> ().swap(current_value_vector);
#endif
            }

            std::vector<Correction_Value> ().swap(corrections_vector);
            std::vector<segment> ().swap(segments);
        }

        constexpr static uint32_t simd_limit = 4; // the minimum length in one simd group
        constexpr static double simd_group_limit = 0.0;
#if SIMD_512_1D1S
        constexpr static int32_t simd_group_length = 16;
        constexpr static K align_val = 64; // for avx512: 16 * 4 bytes
#elif SIMD_256_2D1S
        constexpr static int32_t simd_group_length = 4;
        constexpr static K align_val = 16;
#else
        constexpr static int32_t simd_group_length = 8;
        constexpr static K align_val = 32; // for avx512: 16 * 4 bytes
#endif

        void memory_layout() {
            std::vector<segment> tmp = segments;
            std::sort(tmp.begin(), tmp.end(), [](const segment& a, const segment& b) { return a.covered > b.covered;});

            std::vector<segment> simd;
            std::vector<segment> tail;
            simd.reserve(tmp.size());
            tail.reserve(tmp.size());

            const size_t n = tmp.size();
            size_t i = 0;

            while (i < n) {
                if (n - i < static_cast<size_t>(simd_group_length)) {
                    tail.insert(tail.end(), tmp.begin() + static_cast<ptrdiff_t>(i), tmp.end());
                    break;
                }

                //  [i, i + simd_group_length)
                const size_t head = i;
                const size_t end  = i + static_cast<size_t>(simd_group_length) - 1;

                const double head_cov = static_cast<double>(tmp[head].covered);
                const double tail_cov = static_cast<double>(tmp[end].covered);

                if (tail_cov >= head_cov * simd_group_limit && tail_cov >= simd_limit) {
                    simd.insert(simd.end(),
                                tmp.begin() + static_cast<ptrdiff_t>(i),
                                tmp.begin() + static_cast<ptrdiff_t>(i + simd_group_length));
                    i += static_cast<size_t>(simd_group_length);
                } else {
                    tail.push_back(tmp[i]);
                    ++i;
                }
            }

            std::sort(tail.begin(), tail.end(), [](const segment& a, const segment& b) { return a.first < b.first; });

            segments_sort.clear();
            segments_sort.insert(segments_sort.end(), simd.begin(), simd.end());
            idx = simd.size();
            segments_sort.insert(segments_sort.end(), tail.begin(), tail.end());
        }

        // our Memory Re-Layout
        void simd_init() {
            idx = 0;
            memory_layout();

            Segment_Value min_cover = INT_MAX;
            // Segment_Value max_min_covered = 0;

            uint64_t true_idx = 0;
            // for (auto it = segments_sort.begin(); it + simd_group_length <= segments_sort.end(); it = it + simd_group_length) {
            assert(idx % simd_group_length == 0);
            for (auto it = segments_sort.begin(); it +simd_group_length <= segments_sort.begin() + idx; it = it + simd_group_length) {
                alignas(align_val) Segment_Value *delta_y_simd_tmp = aligned_new<Segment_Value>(simd_group_length);
                alignas(align_val) Segment_Value *delta_x_simd_tmp = aligned_new<Segment_Value>(simd_group_length);
                alignas(align_val) Segment_Value *y_b_simd_tmp = aligned_new<Segment_Value>(simd_group_length);
                alignas(align_val) Segment_Value *x_b_simd_tmp = aligned_new<Segment_Value>(simd_group_length);
                alignas(align_val) Segment_Value *covered_tmp = aligned_new<Segment_Value>(simd_group_length);
                alignas(align_val) Segment_Value *first_tmp = aligned_new<Segment_Value>(simd_group_length);

                // std::vector<segment> simd_segments(it, it + simd_group_length);
                true_idx += simd_group_length;

                int i = 0;
                for (auto its = it; its != it + simd_group_length; ++its, ++i) {
                    auto covered = its -> covered;
                    delta_y_simd_tmp[i] = static_cast<Segment_Value>(its -> delta_y);
                    delta_x_simd_tmp[i] = static_cast<Segment_Value>(its -> delta_x);
                    y_b_simd_tmp[i] = static_cast<Segment_Value>(its -> y_b);
                    x_b_simd_tmp[i] = static_cast<Segment_Value>(its -> x_b);
                    covered_tmp[i] = static_cast<Segment_Value> (covered);
                    first_tmp[i] = static_cast<Segment_Value> (its -> first);
                    min_cover = min_cover < covered ? min_cover : covered;
                }

                delta_y_simd.emplace_back(delta_y_simd_tmp);
                delta_x_simd.emplace_back(delta_x_simd_tmp);
                y_b_simd.emplace_back(y_b_simd_tmp);
                x_b_simd.emplace_back(x_b_simd_tmp);
                first_simd.emplace_back(first_tmp);
                list_covered_simd.emplace_back(covered_tmp);
                // max_min_covered = min_cover;
                // max_min_covered = min_cover - min_cover % 2;
                group_min_covered_simd.emplace_back(min_cover);
            }

            assert(idx == true_idx);

            residuals_decode();
            create_corrections();
            create_corrections_tail();
        }

        void create_corrections() {
            total_calculated = 0;
            for (const auto length: group_min_covered_simd)
                total_calculated += length * simd_group_length;

            corrections_simd.resize(total_calculated, 0);
            uint64_t corrections_pointer = 0;
            for (int i = 0;i < group_min_covered_simd.size(); i++) {
                Segment_Value group_min_covered = group_min_covered_simd[i];
                Segment_Value *first_tmp = first_simd[i];
                for (Segment_Value j = 0; j < group_min_covered; j++) {
                    for (Segment_Value k = 0; k < simd_group_length; k++) {
                        corrections_simd[corrections_pointer++] = corrections_vector[first_tmp[k] + j];
                    }
                }
            }
        }

        void create_corrections_tail() {
            corrections_tail.resize(n - total_calculated);
             Segment_Value correct_pointers = 0;
             for(int i = 0; i < group_min_covered_simd.size(); i++) {
                 Segment_Value *first_tmp = first_simd[i];
                 Segment_Value *list_covered = list_covered_simd[i];
                 Segment_Value group_min_covered = group_min_covered_simd[i];
                 for (K k = 0; k < simd_group_length; k++) {
                     if (group_min_covered < list_covered[k]) {
                         for (Segment_Value pos = group_min_covered; pos < list_covered[k]; pos++){
                             corrections_tail[correct_pointers++] = corrections_vector[first_tmp[k] + pos];
                         }
                     }
                 }
             }

             for(auto it = segments_sort.begin() + idx; it < segments_sort.end(); it++) {
                 for (Segment_Value pos = it -> first; pos < (it -> first + it -> covered); pos++){
                     corrections_tail[correct_pointers++] = corrections_vector[pos];
                 }
             }
        }

#if SIMD_512_1D1S
        void simd_decode_512i_1d1s(K* __restrict output) {
            Correction_Value* __restrict correct_pointer = corrections_tail.data();
            alignas(align_val) Correction_Value *vec_16x32int = aligned_new<Correction_Value>(16);
            alignas(align_val) Correction_Value *vec_16x32int_q = aligned_new<Correction_Value>(16);
            alignas(align_val) Correction_Value *vec_16x32int_r = aligned_new<Correction_Value>(16);
            alignas(align_val) Correction_Value *vec_16x32int_step_q = aligned_new<Correction_Value>(16);
            alignas(align_val) Correction_Value *vec_16x32int_step_r = aligned_new<Correction_Value>(16);
            const Correction_Value *corrections_p = corrections_simd.data();
            const int cover_length_size = group_min_covered_simd.size();
            K* p_out[16];

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < cover_length_size; i++) {

                Segment_Value  *first_tmp = first_simd[i];

                for (int j = 0; j < 16; j++) {
                    p_out[j] = output + first_tmp[j];
                }

                // #Block 1 predictions, we should calculate the first numerator
                // # 1 Load
                const __m512i v_zeros = _mm512_setzero_si512();
                const __m512i v_ones = _mm512_set1_epi32(1);
                __m512i v_xb  = _mm512_load_epi32(x_b_simd[i]);
                __m512i v_xb_plus1 = _mm512_add_epi32(v_xb, v_ones);
                __m512i v_j   = _mm512_load_epi32(first_simd[i]);
                __m512i v_dx   = _mm512_load_epi32(delta_x_simd[i]);
                __m512i v_dxd =  _mm512_srli_epi32(v_dx, 1); // delta_x / 2

                // #1 Calculate numerator
                // numerator = dy * (j - xb) + half * sign(j - xb)
                for (int k = 0; k < 16; k++) {
                    const int64_t dy = delta_y_simd[i][k];
                    const int64_t dx = delta_x_simd[i][k];
                    const int64_t first_sub_xb = first_simd[i][k] - x_b_simd[i][k];
                    const int64_t num = dy * (first_sub_xb) + (dx >> 1) * ((first_sub_xb > 0) - (first_sub_xb < 0));
                    vec_16x32int_q[k] = num / dx;
                    vec_16x32int_r[k] = num - vec_16x32int_q[k] * dx;
                    vec_16x32int_step_q[k] = dy / dx;
                    vec_16x32int_step_r[k] = dy - vec_16x32int_step_q[k] * dx;
                }

                // #1 Calculate slope item
                __m512i v_q = _mm512_load_epi32(vec_16x32int_q);
                __m512i v_r = _mm512_load_epi32(vec_16x32int_r);
                __m512i v_step_q = _mm512_load_epi32(vec_16x32int_step_q);
                __m512i v_step_r = _mm512_load_epi32(vec_16x32int_step_r);

                // #1 residual block
                __m512i v_yb = _mm512_load_epi32(y_b_simd[i]); // y_b
                __m512i v_correction = _mm512_loadu_epi32(corrections_p); // residuals
                corrections_p += 16;
#if RESIDUAL_COMPRESS
                __m512i mask = _mm512_srai_epi32(_mm512_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm512_xor_si512(_mm512_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm512_add_epi32(v_yb, v_correction); // y_b + residuals
                __m512i v_int32 = _mm512_add_epi32(v_q, v_yb); // q + residuals + y_b

                // #1 store back
                _mm512_i32scatter_epi32(output, v_j, v_int32, 4);

                // #2 predications
                // j = j +1
                v_j = _mm512_add_epi32(v_j, v_ones);
                __mmask16 m_positive = _mm512_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                //  q += step_q, r += step_r
                v_q = _mm512_add_epi32(v_q, v_step_q);
                v_r = _mm512_add_epi32(v_r, v_step_r);

                // (j == x_b or j == x_b+1), r += half
                __mmask16 m_ind =  _mm512_cmpeq_epi32_mask(v_j, v_xb) |  _mm512_cmpeq_epi32_mask(v_j, v_xb_plus1);
                v_r = _mm512_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                __mmask16 m_r_exceed = _mm512_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                __mmask16 m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                m_r_exceed = _mm512_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                m_process = m_process | (m_positive & m_r_exceed);

                v_r = _mm512_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                v_q = _mm512_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1


                // #2 residual block
                v_correction = _mm512_loadu_epi32(corrections_p);
                corrections_p += 16;
#if RESIDUAL_COMPRESS
                mask = _mm512_srai_epi32(_mm512_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm512_xor_si512(_mm512_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm512_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                // q + y_b + residuals
                v_int32 = _mm512_add_epi32(v_q, v_yb);

                // store back
                _mm512_i32scatter_epi32(output, v_j, v_int32, 4);

                const Segment_Value cover_length_simd = group_min_covered_simd[i];

                for (int32_t j_step = 2; j_step < cover_length_simd; j_step++) {
                    // #Decode Unit 1
                    // j = j +1
                    v_j = _mm512_add_epi32(v_j, v_ones);
                    m_positive = _mm512_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                    //  q += step_q, r += step_r
                    v_q = _mm512_add_epi32(v_q, v_step_q);
                    v_r = _mm512_add_epi32(v_r, v_step_r);

                    // if (j == x_b or j == x_b+1), r += half
                    m_ind =  _mm512_cmpeq_epi32_mask(v_j, v_xb) |  _mm512_cmpeq_epi32_mask(v_j, v_xb_plus1);
                    if (__builtin_expect(m_ind != 0, 0))
                        v_r = _mm512_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                    // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                    m_r_exceed = _mm512_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                    m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                    m_r_exceed = _mm512_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                    m_process = m_process | (m_positive & m_r_exceed); // (j <= xb && r > 0) || (j > xb && r >= dx)

                    v_r = _mm512_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                    v_q = _mm512_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1

                    // residual block
                    v_correction = _mm512_loadu_epi32(corrections_p);
                    corrections_p += 16;
#if RESIDUAL_COMPRESS
                    mask = _mm512_srai_epi32(_mm512_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm512_xor_si512(_mm512_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm512_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm512_add_epi32(v_q, v_yb);

                    // #Store Unit 1
                    _mm512_i32scatter_epi32(output, v_j, v_int32, 4);
                }

                // last residuals
                _mm512_store_epi32(vec_16x32int, v_yb);
                _mm512_store_epi32(vec_16x32int_q, v_q);
                _mm512_store_epi32(vec_16x32int_r, v_r);
                _mm512_store_epi32(vec_16x32int_step_q, v_step_q);
                _mm512_store_epi32(vec_16x32int_step_r, v_step_r);


                for (Segment_Value k = 0; k < 16; k++) {
                    const Segment_Value covered = list_covered_simd[i][k];
                    if (cover_length_simd < covered) {
                        int32_t j = first_tmp[k] + cover_length_simd;
                        K* __restrict out = output + j;
                        const int32_t dx = delta_x_simd[i][k];
                        const int32_t half = dx >> 1; // delta_x / 2
                        const int32_t xb = x_b_simd[i][k];
                        const int32_t j_end = first_tmp[k] + covered;

                        // step delta
                        const int32_t step_r = vec_16x32int_step_r[k];
                        const int32_t step_q = vec_16x32int_step_q[k];

                        int32_t r = vec_16x32int_r[k];
                        int32_t q = vec_16x32int_q[k];
                        int32_t last_correction = vec_16x32int[k];

                        // first part: j < xb, without half, the numerator is negative
                        const int32_t checkpoint = (xb > j) ? ((j_end < xb) ? j_end : xb) : j;

                        // std::cerr << j << " " << xb << " " << checkpoint << "\n";
                        while (__builtin_expect(j < checkpoint, 0)) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q += step_q;
                            r += step_r;

                            if (__builtin_expect(r > 0, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }

                        // q = 0 when j == xb
                        if (j == xb && j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            *out++ = static_cast<K>(last_correction);
                            ++j;
                        }

                        // one more half step, if j = xb + 1. Notably, we need to double check if j == xb + 1 here, because we do not know where the current j locates, it's different to normal decode.
                        if (j == xb + 1 && j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q = step_q;
                            r = step_r + half;

                            if (__builtin_expect(r >= dx, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }

                        // last part, without half, the numerator is positive
                        while (j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q += step_q;
                            r += step_r;

                            if (__builtin_expect(r >= dx, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }
                    } else {
                        break;
                    }
                }
            }

            simd_tail_process(output, correct_pointer);

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
            aligned_delete(vec_16x32int);
            aligned_delete(vec_16x32int_q);
            aligned_delete(vec_16x32int_r);
            aligned_delete(vec_16x32int_step_q);
            aligned_delete(vec_16x32int_step_r);
        }
#endif

#if SIMD_256_1D1S
        void simd_decode_256i_1d1s(K* __restrict output) {
            Correction_Value* __restrict correct_pointer = corrections_tail.data();
            alignas(align_val) Correction_Value *vec_8x32int = aligned_new<Correction_Value>(8);
            alignas(align_val) Correction_Value *vec_8x32int_q = aligned_new<Correction_Value>(8);
            alignas(align_val) Correction_Value *vec_8x32int_r = aligned_new<Correction_Value>(8);
            alignas(align_val) Correction_Value *vec_8x32int_step_q = aligned_new<Correction_Value>(8);
            alignas(align_val) Correction_Value *vec_8x32int_step_r = aligned_new<Correction_Value>(8);
            const Correction_Value *corrections_p = corrections_simd.data();
            const int cover_length_size = group_min_covered_simd.size();
            K* p_out[8];

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < cover_length_size; i++) {

                Segment_Value  *first_tmp = first_simd[i];

                for (int j = 0; j < 8; j++) {
                    p_out[j] = output + first_tmp[j];
                }

                // #Block 1 predictions, we should calculate the first numerator
                // # 1 Load
                const __m256i v_zeros = _mm256_setzero_si256();
                const __m256i v_ones = _mm256_set1_epi32(1);
                __m256i v_xb  = _mm256_load_epi32(x_b_simd[i]);
                __m256i v_xb_plus1 = _mm256_add_epi32(v_xb, v_ones);
                __m256i v_j   = _mm256_load_epi32(first_simd[i]);
                __m256i v_dx   = _mm256_load_epi32(delta_x_simd[i]);
                __m256i v_dxd =  _mm256_srli_epi32(v_dx, 1); // delta_x / 2

                // #1 Calculate numerator
                // numerator = dy * (j - xb) + half * sign(j - xb)
                for (int k = 0; k < 8; k++) {
                    const int64_t dy = delta_y_simd[i][k];
                    const int64_t dx = delta_x_simd[i][k];
                    const int64_t first_sub_xb = first_simd[i][k] - x_b_simd[i][k];
                    const int64_t num = dy * (first_sub_xb) + (dx >> 1) * ((first_sub_xb > 0) - (first_sub_xb < 0));
                    vec_8x32int_q[k] = num / dx;
                    vec_8x32int_r[k] = num - vec_8x32int_q[k] * dx;
                    vec_8x32int_step_q[k] = dy / dx;
                    vec_8x32int_step_r[k] = dy - vec_8x32int_step_q[k] * dx;
                }

                // #1 Calculate slope item
                __m256i v_q = _mm256_load_epi32(vec_8x32int_q);
                __m256i v_r = _mm256_load_epi32(vec_8x32int_r);
                __m256i v_step_q = _mm256_load_epi32(vec_8x32int_step_q);
                __m256i v_step_r = _mm256_load_epi32(vec_8x32int_step_r);

                // #1 residual block
                __m256i v_yb = _mm256_load_epi32(y_b_simd[i]); // y_b
                __m256i v_correction = _mm256_loadu_epi32(corrections_p); // residuals
                corrections_p += 8;
#if RESIDUAL_COMPRESS
                __m256i mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals
                __m256i v_int32 = _mm256_add_epi32(v_q, v_yb); // q + residuals + y_b

                // #1 store back
                for (int k = 0; k < 8; k++) {
                    *p_out[k] = _mm256_extract_epi32(v_int32, k);
                    *p_out[k]++;
                }


                // #2 predications
                // j = j +1
                v_j = _mm256_add_epi32(v_j, v_ones);
                __mmask8 m_positive = _mm256_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                //  q += step_q, r += step_r
                v_q = _mm256_add_epi32(v_q, v_step_q);
                v_r = _mm256_add_epi32(v_r, v_step_r);

                // (j == x_b or j == x_b+1), r += half
                __mmask8 m_ind =  _mm256_cmpeq_epi32_mask(v_j, v_xb) |  _mm256_cmpeq_epi32_mask(v_j, v_xb_plus1);
                v_r = _mm256_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                __mmask8 m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                __mmask8 m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                m_process = m_process | (m_positive & m_r_exceed);

                v_r = _mm256_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                v_q = _mm256_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1


                // #2 residual block
                v_correction = _mm256_loadu_epi32(corrections_p);
                corrections_p += 8;
#if RESIDUAL_COMPRESS
                mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                // q + y_b + residuals
                v_int32 = _mm256_add_epi32(v_q, v_yb);

                // store back
                for (int k = 0; k < 8; k++) {
                    *p_out[k] = _mm256_extract_epi32(v_int32, k);
                    *p_out[k]++;
                }

                const Segment_Value cover_length_simd = group_min_covered_simd[i];

                Segment_Value j_step = 2;
                // const Segment_Value two_unit = cover_length_simd - cover_length_simd % 2;
                for (; j_step < cover_length_simd; j_step++) {
                    // #Decode Unit 1
                    // j = j +1
                    v_j = _mm256_add_epi32(v_j, v_ones);
                    m_positive = _mm256_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                    //  q += step_q, r += step_r
                    v_q = _mm256_add_epi32(v_q, v_step_q);
                    v_r = _mm256_add_epi32(v_r, v_step_r);

                    // if (j == x_b or j == x_b+1), r += half
                    m_ind =  _mm256_cmpeq_epi32_mask(v_j, v_xb) |  _mm256_cmpeq_epi32_mask(v_j, v_xb_plus1);
                    if (__builtin_expect(m_ind != 0, 0))
                        v_r = _mm256_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                    // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                    m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                    m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                    m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                    m_process = m_process | (m_positive & m_r_exceed); // (j <= xb && r > 0) || (j > xb && r >= dx)

                    v_r = _mm256_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                    v_q = _mm256_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1

                    // residual block
                    v_correction = _mm256_loadu_epi32(corrections_p);
                    corrections_p += 8;
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm256_add_epi32(v_q, v_yb);

                    // #Store Unit 1
                    for (int k = 0; k < 8; k++) {
                        *p_out[k] = _mm256_extract_epi32(v_int32, k);
                        *p_out[k]++;
                    }
                }

                // last residuals
                _mm256_store_epi32(vec_8x32int, v_yb);
                _mm256_store_epi32(vec_8x32int_q, v_q);
                _mm256_store_epi32(vec_8x32int_r, v_r);
                _mm256_store_epi32(vec_8x32int_step_q, v_step_q);
                _mm256_store_epi32(vec_8x32int_step_r, v_step_r);


                for (Segment_Value k = 0; k < 8; k++) {
                    const Segment_Value covered = list_covered_simd[i][k];
                    if (cover_length_simd < covered) {
                        int32_t j = first_tmp[k] + cover_length_simd;
                        K* __restrict out = output + j;
                        const int32_t dx = delta_x_simd[i][k];
                        const int32_t half = dx >> 1; // delta_x / 2
                        const int32_t xb = x_b_simd[i][k];
                        const int32_t j_end = first_tmp[k] + covered;
                        // step delta

                        const int32_t step_r = vec_8x32int_step_r[k];
                        const int32_t step_q = vec_8x32int_step_q[k];

                        int32_t r = vec_8x32int_r[k];
                        int32_t q = vec_8x32int_q[k];
                        int32_t last_correction = vec_8x32int[k];

                        // first part: j < xb, without half, the numerator is negative
                        const int32_t checkpoint = (xb > j) ? ((j_end < xb) ? j_end : xb) : j;

                        // std::cerr << j << " " << xb << " " << checkpoint << "\n";
                        while (__builtin_expect(j < checkpoint, 0)) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q += step_q;
                            r += step_r;

                            if (__builtin_expect(r > 0, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }

                        // q = 0 when j == xb
                        if (j == xb && j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            *out++ = static_cast<K>(last_correction);
                            ++j;
                        }

                        // one more half step, if j = xb + 1. Notably, we need to double check if j == xb + 1 here, because we do not know where the current j locates, it's different to normal decode.
                        if (j == xb + 1 && j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q = step_q;
                            r = step_r + half;

                            if (__builtin_expect(r >= dx, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }

                        // last part, without half, the numerator is positive
                        while (j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q += step_q;
                            r += step_r;

                            if (__builtin_expect(r >= dx, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }
                    } else {
                        break;
                    }
                }
            }

            simd_tail_process(output, correct_pointer);

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
            aligned_delete(vec_8x32int);
            aligned_delete(vec_8x32int_q);
            aligned_delete(vec_8x32int_r);
            aligned_delete(vec_8x32int_step_q);
            aligned_delete(vec_8x32int_step_r);
        }
#endif

#if SIMD_256_2D1S
        void simd_decode_256i_2d1s(K* __restrict output) {
            Correction_Value* __restrict correct_pointer = corrections_tail.data();
            // __m512i rearrange_idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
            __m256i rearrange_idx = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
            alignas(align_val) Correction_Value *vec_4x32int = aligned_new<Correction_Value>(4);
            alignas(align_val) Correction_Value *vec_4x32int_q = aligned_new<Correction_Value>(4);
            alignas(align_val) Correction_Value *vec_4x32int_r = aligned_new<Correction_Value>(4);
            alignas(align_val) Correction_Value *vec_4x32int_step_q = aligned_new<Correction_Value>(4);
            alignas(align_val) Correction_Value *vec_4x32int_step_r = aligned_new<Correction_Value>(4);
            const Correction_Value *corrections_p = corrections_simd.data();
            const int cover_length_size = group_min_covered_simd.size();

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < cover_length_size; i++) {
                Segment_Value  *first_tmp = first_simd[i];
                K *p0 = output + first_tmp[0];
                K *p1 = output + first_tmp[1];
                K *p2 = output + first_tmp[2];
                K *p3 = output + first_tmp[3];

                // #Block 1 predictions, we should calculate the first numerator
                // # 1 Load
                const __m128i v_zeros = _mm_setzero_si128();
                const __m128i v_ones = _mm_set1_epi32(1);
                __m128i v_xb  = _mm_load_epi32(x_b_simd[i]);
                __m128i v_xb_plus1 = _mm_add_epi32(v_xb, v_ones);
                __m128i v_j   = _mm_load_epi32(first_tmp);
                __m128i v_dx   = _mm_load_epi32(delta_x_simd[i]);
                __m128i v_dxd =  _mm_srli_epi32(v_dx, 1); // delta_x / 2

                // #1 Calculate numerator
                // numerator = dy * (j - xb) + half * sign(j - xb)
                for (int k = 0; k < 4; k++) {
                    const int64_t dy = delta_y_simd[i][k];
                    const int64_t dx = delta_x_simd[i][k];
                    const int64_t first_sub_xb = first_simd[i][k] - x_b_simd[i][k];
                    const int64_t num = dy * (first_sub_xb) + (dx >> 1) * ((first_sub_xb > 0) - (first_sub_xb < 0));
                    vec_4x32int_q[k] = num / dx;
                    vec_4x32int_r[k] = num - vec_4x32int_q[k] * dx;
                    vec_4x32int_step_q[k] = dy / dx;
                    vec_4x32int_step_r[k] = dy - vec_4x32int_step_q[k] * dx;
                }

                __m128i v_q = _mm_load_epi32(vec_4x32int_q);
                __m128i v_r = _mm_load_epi32(vec_4x32int_r);
                __m128i v_step_q = _mm_load_epi32(vec_4x32int_step_q);
                __m128i v_step_r = _mm_load_epi32(vec_4x32int_step_r);

                // #1 residual block
                __m128i v_yb = _mm_load_epi32(y_b_simd[i]); // y_b
                __m128i v_correction = _mm_loadu_epi32(corrections_p); // residuals
                corrections_p += 4;
#if RESIDUAL_COMPRESS
                __m128i mask = _mm_srai_epi32(_mm_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm_xor_si128(_mm_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm_add_epi32(v_yb, v_correction); // y_b + residuals
                __m128i v_int32 = _mm_add_epi32(v_q, v_yb); // q + residuals + y_b

                // #1 store to lower 128 bits
                __m256i v_result = _mm256_castsi128_si256(v_int32); // lower 128 bits


                // #2 predications
                // j = j +1
                v_j = _mm_add_epi32(v_j, v_ones);
                __mmask8 m_positive = _mm_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                //  q += step_q, r += step_r
                v_q = _mm_add_epi32(v_q, v_step_q);
                v_r = _mm_add_epi32(v_r, v_step_r);

                // (j == x_b or j == x_b+1), r += half
                __mmask8 m_ind =  _mm_cmpeq_epi32_mask(v_j, v_xb) |  _mm_cmpeq_epi32_mask(v_j, v_xb_plus1);
                v_r = _mm_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                __mmask8 m_r_exceed = _mm_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                __mmask8 m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                m_r_exceed = _mm_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                m_process = m_process | (m_positive & m_r_exceed);

                v_r = _mm_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                v_q = _mm_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1


                // #2 residual block
                v_correction = _mm_loadu_epi32(corrections_p);
                corrections_p += 4;
#if RESIDUAL_COMPRESS
                mask = _mm_srai_epi32(_mm_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm_xor_si128(_mm_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                // q + y_b + residuals
                v_int32 = _mm_add_epi32(v_q, v_yb);

                // store back
                v_result = _mm256_inserti64x2(v_result, v_int32, 1); // upper 256 bits
                v_result = _mm256_permutexvar_epi32(rearrange_idx, v_result); // rearrange

                __m128i t0 = _mm256_extracti32x4_epi32(v_result, 0);
                __m128i t1 = _mm256_extracti32x4_epi32(v_result, 1);

                *((int64_t*)(p0)) = _mm_extract_epi64(t0, 0);  // a0, a1
                *((int64_t*)(p1)) = _mm_extract_epi64(t0, 1);  // b0, b1
                *((int64_t*)(p2)) = _mm_extract_epi64(t1, 0);  // c0, c1
                *((int64_t*)(p3)) = _mm_extract_epi64(t1, 1);  // d0, d1

                p0 = p0 + 2;
                p1 = p1 + 2;
                p2 = p2 + 2;
                p3 = p3 + 2;

                const Segment_Value cover_length_simd = group_min_covered_simd[i];

                Segment_Value j_step = 2;
                const Segment_Value two_unit = cover_length_simd - cover_length_simd % 2;
                for (; j_step < two_unit; j_step += 2) {
                    // #Decode Unit 1
                    // j = j +1
                    v_j = _mm_add_epi32(v_j, v_ones);
                    m_positive = _mm_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                    //  q += step_q, r += step_r
                    v_q = _mm_add_epi32(v_q, v_step_q);
                    v_r = _mm_add_epi32(v_r, v_step_r);

                    // if (j == x_b or j == x_b+1), r += half
                    m_ind =  _mm_cmpeq_epi32_mask(v_j, v_xb) |  _mm_cmpeq_epi32_mask(v_j, v_xb_plus1);
                    if (__builtin_expect(m_ind != 0, 0))
                        v_r = _mm_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                    // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                    m_r_exceed = _mm_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                    m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                    m_r_exceed = _mm_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                    m_process = m_process | (m_positive & m_r_exceed); // (j <= xb && r > 0) || (j > xb && r >= dx)

                    v_r = _mm_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                    v_q = _mm_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1

                    // residual block
                    v_correction = _mm_loadu_epi32(corrections_p);
                    corrections_p += 4;
#if RESIDUAL_COMPRESS
                    mask = _mm_srai_epi32(_mm_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm_xor_si128(_mm_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm_add_epi32(v_q, v_yb);

                    // tmp store to lower 256 bits
                    v_result = _mm256_castsi128_si256(v_int32);

                    // ——————————————————————————————————————————————
                    // #Decode Unit 2
                    // j = j +1
                    v_j = _mm_add_epi32(v_j, v_ones);
                    m_positive = _mm_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                    //  q += step_q, r += step_r
                    v_q = _mm_add_epi32(v_q, v_step_q);
                    v_r = _mm_add_epi32(v_r, v_step_r);

                    // if (j == x_b or j == x_b+1), r += half
                    m_ind =  _mm_cmpeq_epi32_mask(v_j, v_xb) |  _mm_cmpeq_epi32_mask(v_j, v_xb_plus1);
                    if (__builtin_expect(m_ind != 0, 0))
                        v_r = _mm_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                    // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                    m_r_exceed = _mm_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                    m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                    m_r_exceed = _mm_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                    m_process = m_process | (m_positive & m_r_exceed); // (j <= xb && r > 0) || (j > xb && r >= dx)

                    v_r = _mm_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                    v_q = _mm_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1

                    // residual block
                    v_correction = _mm_loadu_epi32(corrections_p);
                    corrections_p += 4;
#if RESIDUAL_COMPRESS
                    mask = _mm_srai_epi32(_mm_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm_xor_si128(_mm_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm_add_epi32(v_q, v_yb);

                    // tmp store to upper 128bits
                    v_result = _mm256_inserti64x2(v_result, v_int32, 1); // upper 128 bits

                    // ——————————————————————————————————————————————
                    // #Store Unit 1
                    v_result = _mm256_permutexvar_epi32(rearrange_idx, v_result); // rearrange

                    t0 = _mm256_extracti32x4_epi32(v_result, 0);
                    t1 = _mm256_extracti32x4_epi32(v_result, 1);


                    *((int64_t*)(p0)) = _mm_extract_epi64(t0, 0);  // a0, a1
                    *((int64_t*)(p1)) = _mm_extract_epi64(t0, 1);  // b0, b1
                    *((int64_t*)(p2)) = _mm_extract_epi64(t1, 0);  // c0, c1
                    *((int64_t*)(p3)) = _mm_extract_epi64(t1, 1);  // d0, d1

                    p0 = p0 + 2;
                    p1 = p1 + 2;
                    p2 = p2 + 2;
                    p3 = p3 + 2;
                }

                if (two_unit < cover_length_simd) {
                    // #Last Decode Unit
                    // j = j +1
                    v_j = _mm_add_epi32(v_j, v_ones);
                    m_positive = _mm_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                    //  q += step_q, r += step_r
                    v_q = _mm_add_epi32(v_q, v_step_q);
                    v_r = _mm_add_epi32(v_r, v_step_r);

                    // if (j == x_b or j == x_b+1), r += half
                    m_ind =  _mm_cmpeq_epi32_mask(v_j, v_xb) |  _mm_cmpeq_epi32_mask(v_j, v_xb_plus1);
                    if (__builtin_expect(m_ind != 0, 0))
                        v_r = _mm_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                    // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                    m_r_exceed = _mm_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                    m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                    m_r_exceed = _mm_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                    m_process = m_process | (m_positive & m_r_exceed); // (j <= xb && r > 0) || (j > xb && r >= dx)

                    v_r = _mm_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                    v_q = _mm_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1

                    // residual block
                    v_correction = _mm_loadu_epi32(corrections_p);
                    corrections_p += 4;
#if RESIDUAL_COMPRESS
                    mask = _mm_srai_epi32(_mm_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm_xor_si128(_mm_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // q + y_b + residuals
                    v_int32 = _mm_add_epi32(v_q, v_yb);

                    _mm_store_epi32(vec_4x32int, v_int32);

                    *p0 = vec_4x32int[0];
                    *p1 = vec_4x32int[1];
                    *p2 = vec_4x32int[2];
                    *p3 = vec_4x32int[3];
                }

                // last residuals
                _mm_store_epi32(vec_4x32int, v_yb);
                _mm_store_epi32(vec_4x32int_q, v_q);
                _mm_store_epi32(vec_4x32int_r, v_r);
                _mm_store_epi32(vec_4x32int_step_q, v_step_q);
                _mm_store_epi32(vec_4x32int_step_r, v_step_r);


                for (Segment_Value k = 0; k < 4; k++) {
                    const Segment_Value covered = list_covered_simd[i][k];
                    if (cover_length_simd < covered) {
                        int32_t j = first_tmp[k] + cover_length_simd;
                        K* __restrict out = output + j;

                        const int64_t dx = delta_x_simd[i][k];
                        const int32_t half = dx >> 1; // delta_x / 2
                        const int32_t xb = x_b_simd[i][k];
                        const int32_t j_end = first_tmp[k] + covered;
                        // step delta
                        const int32_t step_r = vec_4x32int_step_r[k];
                        const int32_t step_q = vec_4x32int_step_q[k];

                        int32_t r = vec_4x32int_r[k];
                        int32_t q = vec_4x32int_q[k];
                        int32_t last_correction = vec_4x32int[k];

                        // first part: j < xb, without half, the numerator is negative
                        const int32_t checkpoint = (xb > j) ? ((j_end < xb) ? j_end : xb) : j;

                        while (__builtin_expect(j < checkpoint, 0)) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q += step_q;
                            r += step_r;

                            if (__builtin_expect(r > 0, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }

                        // q = 0 when j == xb
                        if (j == xb && j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            *out++ = static_cast<K>(last_correction);
                            ++j;
                        }

                        // one more half step, if j = xb + 1. Notably, we need to double check if j == xb + 1 here, because we do not know where the current j locates, it's different to normal decode.
                        if (j == xb + 1 && j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q = step_q;
                            r = step_r + half;

                            if (__builtin_expect(r >= dx, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }

                        // last part, without half, the numerator is positive
                        while (j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q += step_q;
                            r += step_r;

                            if (__builtin_expect(r >= dx, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }
                    } else {
                        break;
                    }
                }
            }

            simd_tail_process(output, correct_pointer);

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
            aligned_delete(vec_4x32int);
            aligned_delete(vec_4x32int_q);
            aligned_delete(vec_4x32int_r);
            aligned_delete(vec_4x32int_step_q);
            aligned_delete(vec_4x32int_step_r);
        }
#endif

#if SIMD_512_2D1S
        void simd_decode_512i_2d1s(K* __restrict output) { // with i32scatter_epi64 and rearrange. we store two consecutive int32 as one int64
            Correction_Value* __restrict correct_pointer = corrections_tail.data();
            __m512i rearrange_idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
            alignas(align_val) Correction_Value *vec_8x32int = aligned_new<Correction_Value>(8);
            alignas(align_val) Correction_Value *vec_8x32int_q = aligned_new<Correction_Value>(8);
            alignas(align_val) Correction_Value *vec_8x32int_r = aligned_new<Correction_Value>(8);
            alignas(align_val) Correction_Value *vec_8x32int_step_q = aligned_new<Correction_Value>(8);
            alignas(align_val) Correction_Value *vec_8x32int_step_r = aligned_new<Correction_Value>(8);
            const Correction_Value *corrections_p = corrections_simd.data();
            const int cover_length_size = group_min_covered_simd.size();

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < cover_length_size; i++) {

                Segment_Value  *first_tmp = first_simd[i];

                // #Block 1 predictions, we should calculate the first numerator
                // # 1 Load
                const __m256i v_zeros = _mm256_setzero_si256();
                const __m256i v_ones = _mm256_set1_epi32(1);
                const __m256i v_xb  = _mm256_load_epi32(x_b_simd[i]);
                const __m256i v_xb_plus1 = _mm256_add_epi32(v_xb, v_ones);
                const __m256i v_dx   = _mm256_load_epi32(delta_x_simd[i]);
                const __m256i v_dxd =  _mm256_srli_epi32(v_dx, 1); // delta_x / 2

                __m256i v_j   = _mm256_load_epi32(first_tmp);

                const __m256i address_offset = v_j;

                K* base_address = output;

                // #1 Calculate numerator
                // numerator = dy * (j - xb) + half * sign(j - xb)
                for (int k = 0; k < 8; k++) {
                    const int64_t dy = delta_y_simd[i][k];
                    const int64_t dx = delta_x_simd[i][k];
                    const int64_t first_sub_xb = first_simd[i][k] - x_b_simd[i][k];
                    const int64_t num = dy * (first_sub_xb) + (dx >> 1) * ((first_sub_xb > 0) - (first_sub_xb < 0));
                    vec_8x32int_q[k] = num / dx;
                    vec_8x32int_r[k] = num - vec_8x32int_q[k] * dx;
                    vec_8x32int_step_q[k] = dy / dx;
                    vec_8x32int_step_r[k] = dy - vec_8x32int_step_q[k] * dx;
                }

                __m256i v_q = _mm256_load_epi32(vec_8x32int_q);
                __m256i v_r = _mm256_load_epi32(vec_8x32int_r);
                __m256i v_step_q = _mm256_load_epi32(vec_8x32int_step_q);
                __m256i v_step_r = _mm256_load_epi32(vec_8x32int_step_r);

                // #1 residual block
                __m256i v_yb = _mm256_load_epi32(y_b_simd[i]); // y_b
                __m256i v_correction = _mm256_loadu_epi32(corrections_p); // residuals
                corrections_p += 8;
#if RESIDUAL_COMPRESS
                __m256i mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals
                __m256i v_int32 = _mm256_add_epi32(v_q, v_yb); // q + residuals + y_b

                // #1 store to lower 256 bits
                __m512i v_result = _mm512_castsi256_si512(v_int32); // lower 256 bits


                // #2 predications
                // j = j +1
                v_j = _mm256_add_epi32(v_j, v_ones);
                __mmask8 m_positive = _mm256_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                //  q += step_q, r += step_r
                v_q = _mm256_add_epi32(v_q, v_step_q);
                v_r = _mm256_add_epi32(v_r, v_step_r);

                // (j == x_b or j == x_b+1), r += half
                __mmask8 m_ind =  _mm256_cmpeq_epi32_mask(v_j, v_xb) |  _mm256_cmpeq_epi32_mask(v_j, v_xb_plus1);
                v_r = _mm256_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                __mmask8 m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                __mmask8 m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                m_process = m_process | (m_positive & m_r_exceed);

                v_r = _mm256_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                v_q = _mm256_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1


                // #2 residual block
                v_correction = _mm256_loadu_epi32(corrections_p);
                corrections_p += 8;
#if RESIDUAL_COMPRESS
                mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                // q + y_b + residuals
                v_int32 = _mm256_add_epi32(v_q, v_yb);

                // store back
                v_result = _mm512_inserti64x4(v_result, v_int32, 1); // upper 256 bits
                v_result = _mm512_permutexvar_epi32(rearrange_idx, v_result); // rearrange

                _mm512_i32scatter_epi64(base_address, address_offset, v_result, 4); // the scale is 4, because output arrary is int32
                base_address += 2;

                const Segment_Value cover_length_simd = group_min_covered_simd[i];

                Segment_Value j_step = 2;
                const Segment_Value two_unit = cover_length_simd - cover_length_simd % 2;
                for (; j_step < two_unit; j_step += 2) {
                    // #Decode Unit 1
                    // j = j +1
                    v_j = _mm256_add_epi32(v_j, v_ones);
                    m_positive = _mm256_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                    //  q += step_q, r += step_r
                    v_q = _mm256_add_epi32(v_q, v_step_q);
                    v_r = _mm256_add_epi32(v_r, v_step_r);

                    // if (j == x_b or j == x_b+1), r += half
                    m_ind =  _mm256_cmpeq_epi32_mask(v_j, v_xb) |  _mm256_cmpeq_epi32_mask(v_j, v_xb_plus1);
                    if (__builtin_expect(m_ind != 0, 0))
                        v_r = _mm256_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                    // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                    m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                    m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                    m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                    m_process = m_process | (m_positive & m_r_exceed); // (j <= xb && r > 0) || (j > xb && r >= dx)

                    v_r = _mm256_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                    v_q = _mm256_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1

                    // residual block
                    v_correction = _mm256_loadu_epi32(corrections_p);
                    corrections_p += 8;
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm256_add_epi32(v_q, v_yb);

                    // tmp store to lower 256 bits
                    v_result = _mm512_castsi256_si512(v_int32);

                    // ——————————————————————————————————————————————
                    // #Decode Unit 2
                    // j = j +1
                    v_j = _mm256_add_epi32(v_j, v_ones);
                    m_positive = _mm256_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                    //  q += step_q, r += step_r
                    v_q = _mm256_add_epi32(v_q, v_step_q);
                    v_r = _mm256_add_epi32(v_r, v_step_r);

                    // if (j == x_b or j == x_b+1), r += half
                    m_ind =  _mm256_cmpeq_epi32_mask(v_j, v_xb) |  _mm256_cmpeq_epi32_mask(v_j, v_xb_plus1);
                    if (__builtin_expect(m_ind != 0, 0))
                        v_r = _mm256_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                    // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                    m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                    m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                    m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                    m_process = m_process | (m_positive & m_r_exceed); // (j <= xb && r > 0) || (j > xb && r >= dx)

                    v_r = _mm256_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                    v_q = _mm256_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1

                    // residual block
                    v_correction = _mm256_loadu_epi32(corrections_p);
                    corrections_p += 8;
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm256_add_epi32(v_q, v_yb);

                    // tmp store to upper 256bits
                    v_result = _mm512_inserti64x4(v_result, v_int32, 1); // upper 256 bits

                    // ——————————————————————————————————————————————
                    // #Store Unit 1
                    v_result = _mm512_permutexvar_epi32(rearrange_idx, v_result); // rearrange

                    _mm512_i32scatter_epi64(base_address, address_offset, v_result, 4);
                    base_address += 2;
                }

                if (two_unit < cover_length_simd) {
                    // #Last Decode Unit
                    // j = j +1
                    v_j = _mm256_add_epi32(v_j, v_ones);
                    m_positive = _mm256_cmp_epi32_mask(v_j, v_xb, _MM_CMPINT_GT); // j > xb mask

                    //  q += step_q, r += step_r
                    v_q = _mm256_add_epi32(v_q, v_step_q);
                    v_r = _mm256_add_epi32(v_r, v_step_r);

                    // if (j == x_b or j == x_b+1), r += half
                    m_ind =  _mm256_cmpeq_epi32_mask(v_j, v_xb) |  _mm256_cmpeq_epi32_mask(v_j, v_xb_plus1);
                    if (__builtin_expect(m_ind != 0, 0))
                        v_r = _mm256_mask_add_epi32(v_r, m_ind, v_r, v_dxd);

                    // j <= xb，if r > 0 then q++ r-=dx ; j > xb, if v_r >= dx then q++ r-=dx
                    m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_zeros, _MM_CMPINT_GT); // r > 0 mask
                    m_process = ~m_positive & m_r_exceed; // j <= xb && r > 0
                    m_r_exceed = _mm256_cmp_epi32_mask(v_r, v_dx, _MM_CMPINT_GE); // r >= dx mask
                    m_process = m_process | (m_positive & m_r_exceed); // (j <= xb && r > 0) || (j > xb && r >= dx)

                    v_r = _mm256_mask_sub_epi32(v_r, m_process, v_r, v_dx); // r -= dx
                    v_q = _mm256_mask_add_epi32(v_q, m_process, v_q, v_ones); // q += 1

                    // residual block
                    v_correction = _mm256_loadu_epi32(corrections_p);
                    corrections_p += 8;
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // q + y_b + residuals
                    v_int32 = _mm256_add_epi32(v_q, v_yb);

                    _mm256_i32scatter_epi32(output, v_j, v_int32, 4);
                }

                // last residuals
                _mm256_store_epi32(vec_8x32int, v_yb);
                _mm256_store_epi32(vec_8x32int_q, v_q);
                _mm256_store_epi32(vec_8x32int_r, v_r);
                _mm256_store_epi32(vec_8x32int_step_q, v_step_q);
                _mm256_store_epi32(vec_8x32int_step_r, v_step_r);

                for (Segment_Value k = 0; k < 8; k++) {
                    const Segment_Value covered = list_covered_simd[i][k];
                    if (cover_length_simd < covered) {
                        int32_t j = first_tmp[k] + cover_length_simd;
                        K* __restrict out = output + j;

                        const int64_t dx = delta_x_simd[i][k];
                        const int32_t half = dx >> 1; // delta_x / 2
                        const int32_t xb = x_b_simd[i][k];
                        const int32_t j_end = first_tmp[k] + covered;
                        // step delta
                        const int32_t step_r = vec_8x32int_step_r[k];
                        const int32_t step_q = vec_8x32int_step_q[k];

                        int32_t r = vec_8x32int_r[k];
                        int32_t q = vec_8x32int_q[k];
                        int32_t last_correction = vec_8x32int[k];

                        // first part: j < xb, without half, the numerator is negative
                        const int32_t checkpoint = (xb > j) ? ((j_end < xb) ? j_end : xb) : j;

                        while (__builtin_expect(j < checkpoint, 0)) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q += step_q;
                            r += step_r;

                            if (__builtin_expect(r > 0, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }

                        // q = 0 when j == xb
                        if (j == xb && j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            *out++ = static_cast<K>(last_correction);
                            ++j;
                        }

                        // one more half step, if j = xb + 1. Notably, we need to double check if j == xb + 1 here, because we do not know where the current j locates, it's different to normal decode.
                        if (j == xb + 1 && j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q = step_q;
                            r = step_r + half;

                            if (__builtin_expect(r >= dx, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }

                        // last part, without half, the numerator is positive
                        while (j < j_end) {
#if RESIDUAL_COMPRESS
                            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                            last_correction += (*correct_pointer++);
#endif
                            q += step_q;
                            r += step_r;

                            if (__builtin_expect(r >= dx, 0)) {
                                r -= dx;
                                ++q;
                            }

                            *out++ = static_cast<K>(q + last_correction);
                            ++j;
                        }
                    } else {
                        break;
                    }
                }
            }

            simd_tail_process(output, correct_pointer);

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
            aligned_delete(vec_8x32int);
            aligned_delete(vec_8x32int_q);
            aligned_delete(vec_8x32int_r);
            aligned_delete(vec_8x32int_step_q);
            aligned_delete(vec_8x32int_step_r);
        }

#endif

        void simd_tail_process(K* __restrict output, Correction_Value* __restrict correct_pointer) {
            const auto end_iter = segments_sort.end();
            for (auto it = segments_sort.begin() + idx; it != end_iter; ++it) {
                const auto& seg = *it;

                const int32_t covered = seg.covered;
                if (__builtin_expect(covered <= 0, 0)) { continue; }

                const int64_t dy = seg.delta_y;
                const int64_t dx = seg.delta_x; // 假定 dx > 0
                const int32_t half = (dx >> 1); // floor(dx/2)
                const int32_t xb = seg.x_b;

                int32_t last_correction = seg.y_b;
#if RESIDUAL_COMPRESS
                last_correction += (unzigzag_int32(*correct_pointer++));
#else
                last_correction += (*correct_pointer++);
#endif

                int32_t j = seg.first;
                K* __restrict out = output + j;
                const int32_t j_end = j + covered;

                // numerator = dy * (j - xb) + half * sign(j - xb)
                int64_t num = dy * (int64_t(j) - int64_t(xb)) + (int64_t) half * ((int64_t)(j > xb) - (int64_t)(j < xb));

                // simulate division
                int32_t q = num / dx;
                int32_t r = num - q * dx;

                // the first value
                *out++ = static_cast<K>(q + last_correction);
                ++j;

                if (__builtin_expect(j >= j_end, 0)) { continue; }

                // step delta
                const int32_t step_q = dy / dx;
                const int32_t step_r = dy - step_q * dx;

                // first part: j < xb, without half, the numerator is negative
                const int32_t e1 = (xb > j) ? ((j_end < xb) ? j_end : xb) : j;
                while (j < e1) {
#if RESIDUAL_COMPRESS
                    last_correction += (unzigzag_int32(*correct_pointer++));
#else
                    last_correction += (*correct_pointer++);
#endif
                    q += step_q;
                    r += step_r;

                    if (r > 0) {
                        r -= dx;
                        ++q;
                    }

                    *out++ = static_cast<K>(q + last_correction);
                    ++j;
                }

                if (j < j_end && j == xb) {
#if RESIDUAL_COMPRESS
                    last_correction += (unzigzag_int32(*correct_pointer++));
#else
                    last_correction += (*correct_pointer++);
#endif
                    *out++ = static_cast<K>(last_correction);
                    ++j;
                }

                // one more half step, if j = xb + 1
                if (j < j_end) {
#if RESIDUAL_COMPRESS
            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                    last_correction += (*correct_pointer++);
#endif
                    q = step_q;
                    r = step_r + half;

                    if (__builtin_expect(r >= dx, 0)) {
                        r -= dx;
                        ++q;
                    }

                    *out++ = static_cast<K>(q + last_correction);
                    ++j;
                }

                // last part, without half, the numerator is positive
                while (j < j_end) {
#if RESIDUAL_COMPRESS
            last_correction += (unzigzag_int32(*correct_pointer++));
#else
                    last_correction += (*correct_pointer++);
#endif
                    q += step_q;
                    r += step_r;

                    if (__builtin_expect(r >= dx, 0)) {
                        r -= dx;
                        ++q;
                    }

                    *out++ = static_cast<K>(q + last_correction);
                    ++j;
                }
            }
        }
    };


    template <typename K>
    static lico_enumerator<K> create_enumerator_from_indexes(std::vector<LICOIndex> &indexes) {
        lico_enumerator<K> enumerator;
        std::vector<int32_t> all_corrections_vector;
        std::vector<uint64_t> all_parted_size;
        std::vector<uint32_t> block_sizes_vector;
        std::vector<std::vector<uint32_t>> all_corrections_compress_fastpfor;

        uint64_t total_data_size = 0;
        int32_t last_first = 0;

        for (auto& idx : indexes) {
            idx.segment_init();
            total_data_size += idx.n;
            all_parted_size.push_back(idx.n);
            block_sizes_vector.push_back(idx.segments_size);

            for (int i = 0; i < idx.segments_size; ++i) {
                auto first = last_first;
                auto delta_y = idx.seg_delta_y[i];
                auto delta_x = idx.seg_delta_x[i];
                auto y_b = idx.seg_y_b[i];
                auto x_b = idx.seg_x_b[i];
                auto covered = idx.seg_covered[i];
                last_first += covered;
                enumerator.segments.emplace_back(first, delta_y, delta_x, y_b, x_b, covered);
            }

#if RESIDUAL_COMPRESS
            if (residual_compress_type == "fastpfor") {
                all_corrections_compress_fastpfor.push_back(idx.corrections_compress);
            } else {
                throw std::runtime_error("residual_compress_type not recognised");
            }
#else
                idx.residual_init();
                all_corrections_vector.insert(all_corrections_vector.end(), idx.corrections_vector.begin(), idx.corrections_vector.end());
#endif
        }

#if RESIDUAL_COMPRESS
        if (residual_compress_type == "fastpfor") {
            enumerator.load_residuals_fastpfor(total_data_size, all_corrections_compress_fastpfor, all_parted_size);
        } else {
            throw std::runtime_error("residual_compress_type not recognised");
        }
#else
        enumerator.load_residuals(total_data_size, all_corrections_vector);
#endif

        enumerator.load_block_size(block_sizes_vector);

        return enumerator;
    }

    template <typename K>
    static lico_enumerator<K> create_enumerator_from_single_index(LICOIndex &idx) {
        lico_enumerator<K> enumerator;

        uint32_t last_first = 0;

        idx.segment_init();

        for (int i = 0; i < idx.segments_size; ++i) {
            auto first = last_first;
            auto delta_y = idx.seg_delta_y[i];
            auto delta_x = idx.seg_delta_x[i];
            auto y_b = idx.seg_y_b[i];
            auto x_b = idx.seg_x_b[i];
            auto covered = idx.seg_covered[i];
            last_first += covered;
            enumerator.segments.emplace_back(first, delta_y, delta_x, y_b, x_b, covered);
        }

#if RESIDUAL_COMPRESS
        if (residual_compress_type == "fastpfor") {
            std::vector<uint64_t> parted_size{idx.n};
            std::vector<std::vector<uint32_t>> compress_fastpfor{idx.corrections_compress};
            enumerator.load_residuals_fastpfor(idx.n, compress_fastpfor, parted_size);
        } else {
            throw std::runtime_error("residual_compress_type not recognised");
        }
#else
            idx.residual_init();
            enumerator.load_residuals(idx.n, idx.corrections_vector);
#endif
        enumerator.load_block_size(std::vector<uint32_t>{idx.segments_size});
        idx.free_memory();

        return enumerator;
    }

}

