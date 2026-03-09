#pragma once

#include <vp2intersect.h>
#include <vp2union.hpp>
#include <lico_kernel.hpp>
#include <lico_enumerate.hpp>
#include <mm_file.hpp>

namespace lico_sequence {
    template <typename K, uint64_t epsilon = 64> // K is uint32_t or uint64_t, Floating is unused
    class lico_querier{
        uint64_t data_size = 0;
        uint32_t query_num = 0;
        uint64_t query_time = 0;
        std::string input_basename = "";

        std::vector<std::vector<uint32_t>> read_query(const std::string& filename) {
            std::vector<std::vector<uint32_t>> idLists;
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return idLists; // Return an empty vector if the file could not be opened.
            }
            std::string line;
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                std::vector<uint32_t> ids;
                uint32_t id;
                while (iss >> id) {                 // Extract uint32_t from the line until no more can be found.
                    ids.push_back(id);
                }
                remove_duplicate_terms(ids);
                idLists.push_back(ids);
            }
            query_num = idLists[0].size();
            std::cout << "Total query sequences: " << idLists.size() << " Query num: " << (query_num > 5 ? 5 : query_num) << std::endl;
            file.close();
            return idLists;
        }

        std::vector<std::vector<uint32_t>> read_query_random(const std::string& filename) {
            std::vector<std::vector<uint32_t>> idLists;
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return idLists; // Return an empty vector if the file could not be opened.
            }
            std::string line;
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                std::vector<uint32_t> ids;
                uint32_t id;
                while (iss >> id) {                 // Extract uint32_t from the line until no more can be found.
                    ids.push_back(id);
                }
                // remove_duplicate_terms(ids);
                idLists.push_back(ids);
            }
            query_num = idLists[0].size();
            std::cout << "Total query sequences: " << idLists.size() << " Query num: " << (query_num > 5 ? 5 : query_num) << std::endl;
            file.close();
            return idLists;
        }

        std::vector<lico_enumerator <K>> load_model(std::vector<uint32_t> idx_list) {
            if (input_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                throw std::runtime_error("file format error");
            }

            std::vector<lico_enumerator <K>> index_sequences;

            std::ifstream in_header(input_basename + "idx.size", std::ios::binary);
            if (!in_header) {
                std::cerr << "Error: Cannot open idx.size for reading." << std::endl;
                throw std::runtime_error("File open error");
            }

            in_header.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            K index_num = 0;
            in_header.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));

            if (epsilon == 0) {
                for (K i = 0; i < index_num; ++i) {
                    size_t partition_size;
                    in_header.read(reinterpret_cast<char*>(&partition_size), sizeof(partition_size));
                    auto it = std::find(idx_list.begin(), idx_list.end(), i);
                    if (it != idx_list.end()) {
                        std::vector<LICOIndex> partition_index;
                        partition_index.reserve(partition_size);

                        for (size_t j = 0; j < partition_size; ++j) {
                            std::string filename = input_basename + std::to_string(i) + "_" + std::to_string(j) + ".idx";
                            std::ifstream in(filename, std::ios::binary);
                            if (!in) {
                                std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                                throw std::runtime_error("File open error");
                            }

                            uint32_t Epsilon_Data;
                            in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));

                            LICOIndex variant_index = lico::LICO<K>(Epsilon_Data);

                            read_index_data(in, variant_index);

                            partition_index.push_back(std::move(variant_index));
                            in.close();
                        }

                        lico_enumerator <K> enumerator_tmp = create_enumerator_from_indexes <K>(partition_index);

                        index_sequences.emplace_back(enumerator_tmp);
                    }
                }
                if (index_sequences.size() != idx_list.size()) {
                    throw std::runtime_error("index_sequences.size() != idx_list.size()");
                }
            }
            else {
                for (auto i : idx_list) {
                    if (i >= index_num) {
                        std::cerr << "Index count out of range: " << i << "Total index: " << index_num << std::endl;
                        return {};
                    }
                    std::string filename = input_basename + std::to_string(i) + ".idx";

                    std::ifstream in(filename, std::ios::binary);
                    if (!in) {
                        std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                        throw std::runtime_error("File open error");
                    }

                    uint32_t Epsilon_Data;
                    in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));

                    LICOIndex variant_index = lico::LICO<K>(Epsilon_Data);

                    read_index_data(in, variant_index);

                    in.close();

                    lico_enumerator <K> enumerator_tmp = create_enumerator_from_single_index<K>(variant_index);

                    index_sequences.emplace_back(enumerator_tmp);
                }
            }

            in_header.close();

            return index_sequences;
        }

        lico_enumerator <K> load_single_model(const uint32_t idx) {
            if (input_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                throw std::runtime_error("file format error");
            }

            lico_enumerator <K> index_sequence;

            std::ifstream in_header(input_basename + "idx.size", std::ios::binary);

            if (!in_header) {
                std::cerr << "Error: Cannot open idx.size for reading." << std::endl;
                throw std::runtime_error("File open error");
            }

            in_header.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            K index_num = 0;
            in_header.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));

            if (epsilon == 0) {
                if (idx > index_num)
                    throw std::runtime_error("LICO Index count out of range");
                for (K i = 0; i < index_num; ++i) {
                    size_t partition_size;
                    in_header.read(reinterpret_cast<char*>(&partition_size), sizeof(partition_size));
                    if (i == idx) {
                        std::vector<LICOIndex> partition_index;
                        partition_index.reserve(partition_size);

                        for (size_t j = 0; j < partition_size; ++j) {
                            std::string filename = input_basename + std::to_string(i) + "_" + std::to_string(j) + ".idx";
                            std::ifstream in(filename, std::ios::binary);
                            if (!in) {
                                std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                                throw std::runtime_error("File open error");
                            }

                            uint32_t Epsilon_Data;
                            in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));

                            LICOIndex variant_index = lico::LICO<K>(Epsilon_Data);

                            read_index_data(in, variant_index);

                            partition_index.push_back(std::move(variant_index));
                            in.close();
                        }

                        index_sequence = create_enumerator_from_indexes <K>(partition_index);
                    }
                }
            }
            else {
                if (idx >= index_num) {
                    std::cerr << "Index count out of range: " << idx << "Total index: " << index_num << std::endl;
                    return {};
                }
                std::string filename = input_basename + std::to_string(idx) + ".idx";

                std::ifstream in(filename, std::ios::binary);
                if (!in) {
                    std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                    throw std::runtime_error("File open error");
                }

                uint32_t Epsilon_Data;
                in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));

                LICOIndex variant_index = lico::LICO<K>(Epsilon_Data);

                read_index_data(in, variant_index);

                in.close();

                index_sequence = create_enumerator_from_single_index<K>(variant_index);
            }

            in_header.close();

            return index_sequence;
        }

        std::vector<K> load_real_data (const uint32_t idx) {
            mm::file_source<K> input(dataset_filename_debug.c_str(), mm::advice::sequential);
            K const* data = input.data();
            int list_idx = 0;
            for (size_t i = 2; i < input.size();){
                uint64_t n = data[i];
                if (list_idx == idx) {
                    std::vector<K> sequence(data + i + 1, data + i + n + 1);
                    // input.close();
                    return sequence;
                }
                list_idx++;
                i += n + 1;
            }
            input.close();
            return {};
        }

        static uint32_t intersect_u32_normal(uint32_t const *a, uint32_t const *b,uint32_t a_length, uint32_t b_length, uint32_t *out) {
            uint32_t* const intersect_start = out;

            const uint32_t* a_end = a + a_length;
            const uint32_t* b_end = b + b_length;

            while (a < a_end && b < b_end) {
                if (*a < *b) {
                    ++a;
                } else if (*a > *b) {
                    ++b;
                } else {
                    // Found a match
                    *out++ = *a;
                    ++a;
                    ++b;
                }
            }
            return out - intersect_start;
        }

        static uint32_t intersect_u32_simd(uint32_t const* a, uint32_t const* b, uint32_t a_length, uint32_t b_length, uint32_t *out) {
            uint32_t const* const a_end = a + a_length;
            uint32_t const* const b_end = b + b_length;
            uint32_t c = 0;
            union vec_t {
                __m512i zmm;
                uint32_t u32[16];
            } a_vec, b_vec;

            while (a + 16 < a_end && b + 16 < b_end) {
                a_vec.zmm = _mm512_loadu_si512((__m512i const*)a);
                b_vec.zmm = _mm512_loadu_si512((__m512i const*)b);

                // Intersecting registers with `_mm512_2intersect_epi16_mask` involves a lot of shuffling
                // and comparisons, so we want to avoid it if the slices don't overlap at all
                uint32_t a_min;
                uint32_t a_max = a_vec.u32[15];
                uint32_t b_min = b_vec.u32[0];
                uint32_t b_max = b_vec.u32[15];

                // If the slices don't overlap, advance the appropriate pointer
                while (a_max < b_min && a + 32 < a_end) {
                    a += 16;
                    a_vec.zmm = _mm512_loadu_si512((__m512i const*)a);
                    a_max = a_vec.u32[15];
                }
                a_min = a_vec.u32[0];
                while (b_max < a_min && b + 32 < b_end) {
                    b += 16;
                    b_vec.zmm = _mm512_loadu_si512((__m512i const*)b);
                    b_max = b_vec.u32[15];
                }
                b_min = b_vec.u32[0];

                // Now we are likely to have some overlap, so we can intersect the registers
                __mmask16 a_matches = _mm512_2intersect_epi32_mask(a_vec.zmm, b_vec.zmm);
                _mm512_mask_compressstoreu_epi32(out + c, a_matches, a_vec.zmm);


                c += _mm_popcnt_u32(a_matches); // The `_popcnt32` symbol isn't recognized by MSVC

                // Determine the number of entries to skip in each array, by comparing
                // every element in the vector with the last (largest) element in the other array
                __m512i a_last_broadcasted = _mm512_set1_epi32(*(int const*)&a_max);
                __m512i b_last_broadcasted = _mm512_set1_epi32(*(int const*)&b_max);
                __mmask16 a_step_mask = _mm512_cmple_epu32_mask(a_vec.zmm, b_last_broadcasted);
                __mmask16 b_step_mask = _mm512_cmple_epu32_mask(b_vec.zmm, a_last_broadcasted);
                a += 16 - __lzcnt16((uint16_t)a_step_mask);
                b += 16 - __lzcnt16((uint16_t)b_step_mask);
            }

            // Handle the tail:
            c += intersect_u32_normal(a, b, a_end - a, b_end - b, out + c);
            return c;
        }

        long double avg_skip = 0;
        long double avg_query_total_size = 0;
        long double avg_query_real_size = 0;

        void remove_duplicate_terms(std::vector<uint32_t>& terms) {
            std::sort(terms.begin(), terms.end());
            terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
        }

        void query_test_intersection(const std::vector<std::vector<uint32_t>> &query_list, const std::string &decode_type) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_round = 0;
            uint64_t avg_time_per = 0;
            K repeat_num = 1;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                avg_time_round = 0;
                for (auto &query : query_list) {
                    std::vector<lico_enumerator <K>> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const lico_enumerator <K> &a, const lico_enumerator <K> &b) {return a.n < b.n;});

                    // int branch_num = query.size() == 2 ? 2 : query.size() == 3 ? 3 : 4;
                    int branch_num = query.size() == 2 ? 2 : 3;

                    for (int i = 0; i < branch_num; i++) {
                        index_sequences[i].query_init(decode_type, "intersection");
                    }

                    uint32_t query_id_idx = 0;
                    uint32_t candidate_posting = 0;
                    uint32_t equal_result = 0;
                    uint32_t candidate_posting_tmp = 0;

                    K *intersection_result_p1, *intersection_result_p2;
#if USE_HUGEPAGE
                    std::vector<K, HugePageAllocator<K>> intersection_result_1(universe_size); // for simd
                    std::vector<K, HugePageAllocator<K>> intersection_result_2(universe_size);
#else
                    std::vector<K> intersection_result_1(universe_size); // for normal
                    std::vector<K> intersection_result_2(universe_size);
#endif
                    intersection_result_p1 = intersection_result_1.data();
                    intersection_result_p2 = intersection_result_2.data();


                    avg_time_per = 0;

                    if (query.size() == 1) {
                        throw std::runtime_error("query is empty");
                    }
                    else if (query.size() == 2) {
                        index_sequences[0].decode_query(intersection_result_p1, decode_type);
                        index_sequences[1].decode_query(intersection_result_p2, decode_type);
                        avg_time_per += index_sequences[0].total_duration;
                        avg_time_per += index_sequences[1].total_duration;

                        auto start = std::chrono::high_resolution_clock::now();
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, index_sequences[0].n,  index_sequences[1].n, intersection_result_p1);
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();

                    } else {
                        // first intersection
                        index_sequences[0].decode_query(intersection_result_p1, decode_type);
                        index_sequences[1].decode_query(intersection_result_p2, decode_type);
                        avg_time_per += index_sequences[0].total_duration;
                        avg_time_per += index_sequences[1].total_duration;

                        auto start = std::chrono::high_resolution_clock::now();
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, index_sequences[0].n, index_sequences[1].n, intersection_result_p1);
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();

                        // second intersection
                        index_sequences[2].decode_query(intersection_result_p2, decode_type);
                        avg_time_per += index_sequences[2].total_duration;

                        start = std::chrono::high_resolution_clock::now();
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, equal_result, index_sequences[2].n, intersection_result_p1);
                        end = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();

                        if (branch_num > 3) {
                            index_sequences[3].decode_query(intersection_result_p2, decode_type);
                            avg_time_per += index_sequences[3].total_duration;

                            start = std::chrono::high_resolution_clock::now();
                            equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, equal_result, index_sequences[3].n, intersection_result_p1);
                            end = std::chrono::high_resolution_clock::now();
                            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                            avg_time_per += duration.count();

                        }

                        // next_geq candidate
                        if (index_sequences.size() > branch_num) {
                            for (query_id_idx = 3; query_id_idx < index_sequences.size(); query_id_idx++) {
                                index_sequences[query_id_idx].next_geq_init();
                                avg_time_per += index_sequences[query_id_idx].total_duration;
                            }

                            uint32_t intersection_length = equal_result;
                            equal_result = 0;
                            uint32_t intersection_idx = 0;
                            bool found = true;

                            start = std::chrono::high_resolution_clock::now();

                            while (intersection_idx < intersection_length) {
                                candidate_posting = intersection_result_p1[intersection_idx];
                                for (query_id_idx = 3; query_id_idx < index_sequences.size(); query_id_idx++) {
                                    candidate_posting_tmp = index_sequences[query_id_idx].nextgeq(candidate_posting);
                                    if (candidate_posting_tmp != candidate_posting) {
                                        candidate_posting = candidate_posting_tmp;
                                        found = false;
                                        break;
                                    }
                                }
                                if (found) {
                                    intersection_result_p2[equal_result++] = candidate_posting;
                                    intersection_idx++;
                                } else {
                                    intersection_idx++;
                                    // intersection_idx = lower_bound(intersection_result_p1 + intersection_idx + 1, intersection_result_p1 + intersection_length, candidate_posting) - intersection_result_p1;
                                    while (intersection_idx < intersection_length && intersection_result_p1[intersection_idx] < candidate_posting)
                                        intersection_idx++;
                                    found = true;
                                }
                            }

                            end = std::chrono::high_resolution_clock::now();
                            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                            avg_time_per += duration.count();
                        }
                    }


                    if (repeat == 0) {
                        total += equal_result;
                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;
                    }

                    avg_time_round += avg_time_per;

                    if (decode_type == "simd") {
                        for (auto &enumerator : index_sequences) {
                            enumerator.free_memory();
                        }
                    }
                }
                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
                if (repeat > 0)
                    querys_per_times.push_back(avg_time_round);
            }

            std::sort(querys_per_times.begin(), querys_per_times.end());
            uint64_t avg_time = std::accumulate(querys_per_times.begin(), querys_per_times.end(), 0ULL);

            std::cerr << "Average query time: " <<  static_cast<long double> (avg_time) / query_list.size() / repeat_num / 1000 << ", " << "Median query time: " << static_cast<long double> (querys_per_times[querys_per_times.size() > 1 ? ((querys_per_times.size() + 1) / 2) : 0]) / query_list.size() / 1000 << std::endl;
            std::cerr << "Average skip rate: " << avg_skip / query_list.size() << ", Average query total size: " << avg_query_total_size / query_list.size() << ", Average query real size: " << avg_query_real_size / query_list.size() << std::endl;
        }

        void query_test_union(const std::vector<std::vector<uint32_t>> &query_list, const std::string &decode_type) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_per = 0;
            uint64_t avg_time_round = 0;
            K repeat_num = 1;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                avg_time_round = 0;
                for (auto &query : query_list) {
                    std::vector<lico_enumerator <K>> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const lico_enumerator <K> &a, const lico_enumerator <K> &b) {return a.n < b.n;});

                    for (int i = 0; i < query.size(); i++) {
                        index_sequences[i].query_init(decode_type, "union");
                    }

                    uint64_t equal_result = 0;
                    K *union_result_p1, *union_result_p2, *union_result_p3;
#if USE_HUGEPAGE
                    std::vector<K, HugePageAllocator<K>>  union_result_1(universe_size); // for simd
                    std::vector<K, HugePageAllocator<K>>  union_result_2(universe_size);
                    std::vector<K, HugePageAllocator<K>>  union_result_3(universe_size);
# else
                    std::vector<K> union_result_1(universe_size); // for normal
                    std::vector<K> union_result_2(universe_size);
                    std::vector<K> union_result_3(universe_size);
#endif
                    union_result_p1 = union_result_1.data();
                    union_result_p2 = union_result_2.data();
                    union_result_p3 = union_result_3.data();

                    avg_time_per = 0;

                    if (query.size() == 1) {
                        throw std::runtime_error("query is empty");
                    }
                    if (query.size() == 2) {
                        index_sequences[0].decode_query(union_result_p1, decode_type);
                        avg_time_per += index_sequences[0].total_duration;
                        index_sequences[1].decode_query(union_result_p2, decode_type);
                        avg_time_per += index_sequences[1].total_duration;

                        auto start = std::chrono::high_resolution_clock::now();
                        equal_result = union_u32_simd(union_result_p1, union_result_p2, index_sequences[0].n, index_sequences[1].n, union_result_p3);
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();
                    } else {
                        index_sequences[0].decode_query(union_result_p1, decode_type);
                        avg_time_per += index_sequences[0].total_duration;
                        index_sequences[1].decode_query(union_result_p2, decode_type);
                        avg_time_per += index_sequences[1].total_duration;

                        auto start = std::chrono::high_resolution_clock::now();
                        equal_result = union_u32_simd(union_result_p1, union_result_p2, index_sequences[0].n, index_sequences[1].n, union_result_p3);
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();

                        for (int i = 2; i < query.size(); i++) {
                            index_sequences[i].decode_query(union_result_p1, decode_type);
                            avg_time_per += index_sequences[i].total_duration;

                            start = std::chrono::high_resolution_clock::now();
                            equal_result = union_u32_simd(union_result_p3, union_result_p1, equal_result, index_sequences[i].n, union_result_p2);
                            end = std::chrono::high_resolution_clock::now();
                            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                            avg_time_per += duration.count();

                            K *union_result_tmp = union_result_p2;
                            union_result_p2 = union_result_p3;
                            union_result_p3 = union_result_tmp;
                        }
                    }

                    if (repeat == 0) {
                        total += equal_result;

                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;
                    }

                    avg_time_round += avg_time_per;
                    if (decode_type == "simd") {
                        for (auto &enumerator : index_sequences) {
                            enumerator.free_memory();
                        }
                    }
                }
                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
                if (repeat > 0)
                    querys_per_times.push_back(avg_time_round);
            }
            std::sort(querys_per_times.begin(), querys_per_times.end());
            uint64_t avg_time = std::accumulate(querys_per_times.begin(), querys_per_times.end(), 0LL);

            std::cerr << "Average query time: " <<  static_cast<long double> (avg_time) / query_list.size() / repeat_num / 1000 << ", " << "Median query time: " << static_cast<long double> (querys_per_times[querys_per_times.size() > 1 ? ((querys_per_times.size() + 1) / 2) : 0]) / query_list.size() / 1000 << std::endl;
            std::cerr << "Average skip rate: " << avg_skip / query_list.size() << ", Average query total size: " << avg_query_total_size / query_list.size() << ", Average query real size: " << avg_query_real_size / query_list.size() << std::endl;
        }

        void random_access_debug(const std::vector<std::vector<uint32_t>> &query_list) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_per = 0;
            uint64_t avg_time_round = 0;
            uint64_t total = 0;

            for (auto &query : query_list) {
                std::vector<K> real_sequence = load_real_data(query[0]);

#if USE_HUGEPAGE
                std::vector<K, HugePageAllocator<K>> result_1(universe_size); // for simd
# else
                std::vector<K> result_1(universe_size); // for normal
#endif
                K *result_p1 = result_1.data();
                const K* result_p1_start = result_p1;

                if (query.size() == 1) {
                    throw std::runtime_error("query is empty");
                }

                auto start = std::chrono::high_resolution_clock::now();
                for (int32_t i = 1; i < query.size(); i++) {
                    K posting_list = query[i];
                    K candidate = *(std::lower_bound(real_sequence.begin(), real_sequence.end(), posting_list));
                    if (candidate == posting_list) {
                        total++;
                        // std::cerr << candidate << " ";
                    }

                        // *result_p1++ = candidate;
                }
                // std::cerr << std::endl;
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                avg_time_per = duration.count();


                // total += result_p1 - result_p1_start;

                avg_time_round += avg_time_per;
            }
            std::cerr << "Total size: " << total << std::endl;
        }

        void query_test_random_access(const std::vector<std::vector<uint32_t>> &query_list, const std::string &decode_type) {
            std::vector<uint64_t> query_per_times;
            std::vector<uint64_t> residual_per_times;
            uint64_t avg_time_per = 0;
            uint64_t total_list_size = 0;
            uint64_t avg_time_round = 0, avg_residual_decode_round = 0;
            int32_t repeat_num = 3;

            for (int32_t repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                avg_time_round = 0;
                avg_residual_decode_round = 0;

                for (auto &query : query_list) {
                    lico_enumerator <K> index_sequence = load_single_model(query[0]);

                    index_sequence.random_access_init();

#if USE_HUGEPAGE
                    std::vector<K, HugePageAllocator<K>> result_1(universe_size); // for simd
# else
                    std::vector<K> result_1(universe_size); // for normal
#endif
                    K *result_p1 = result_1.data();
                    const K* result_p1_start = result_p1;

                    if (query.size() == 1) {
                        throw std::runtime_error("query is empty");
                    }

                    auto start = std::chrono::high_resolution_clock::now();
                    for (int32_t i = 1; i < query.size(); i++) {
                        K position = query[i];
                        *result_p1++ = index_sequence.random_access(position);
                    }
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                    avg_time_per = index_sequence.total_duration + duration.count();


                    if (repeat == 0) {
                        total += result_p1 - result_p1_start;
                        total_list_size += index_sequence.n;
                    }

                    avg_residual_decode_round += index_sequence.residuals_resume_duration;
                    avg_time_round += avg_time_per;
                }

                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
                if (repeat > 0) {
                    query_per_times.push_back(avg_time_round);
                    residual_per_times.push_back(avg_residual_decode_round);
                }

            }

            std::sort(query_per_times.begin(), query_per_times.end());
            uint64_t avg_time = std::accumulate(query_per_times.begin(), query_per_times.end(), 0LL);

            std::cerr << "Average residual decode per int time: " << static_cast<long double> (std::accumulate(residual_per_times.begin(), residual_per_times.end(), 0LL)) / query_list.size() / repeat_num * 1000 * 1000 / total_list_size << std::endl;
            std::cerr << "Average query time: " <<  static_cast<long double> (avg_time) / query_list.size() / repeat_num / 1000 << ", " << "Median query time: " << static_cast<long double> (query_per_times[query_per_times.size() > 1 ? ((query_per_times.size() + 1) / 2) : 0]) / query_list.size() / 1000 << std::endl;
        }

        static void intersection_candidate_test(std::vector<lico_enumerator <K>> &index_sequences, K* intersection_result_p1, uint32_t &equal_result, uint32_t &query_id_idx, uint32_t &candidate_posting_tmp, const uint32_t m) {
            uint32_t candidate_posting = index_sequences[0].docid();
            K* intersection_start = intersection_result_p1;

            while (candidate_posting < INT_MAX) {
                for (; query_id_idx < m; ++query_id_idx) { // find the same docID for all words
                    candidate_posting_tmp = index_sequences[query_id_idx].nextgeq(candidate_posting); // for word i  get next docID greater or equal to candidate
                    if (candidate_posting_tmp != candidate_posting) { // if exist a docID different from candidate
                        candidate_posting = candidate_posting_tmp; // candidate is the new docID
                        query_id_idx= 0; // i restarts
                        break; // break the loop and read the new candidate
                    }
                }
                if (query_id_idx == m) { // if all words have the same docID
                    *intersection_result_p1++ = candidate_posting; // add to the intersection
                    candidate_posting = index_sequences[0].nextgeq(candidate_posting + 1); // get a new candidate
                    query_id_idx = 1; // i restarts
                }

            }

            equal_result = intersection_result_p1 - intersection_start;
        }

    public:
        K universe_size;
        std::string dataset_filename_debug;
        void test_query(const std::string input_filename, const std::string &decode_type, const std::string &query_filename, const std::string query_type, const std::string &dataset_filename) {
            mm::file_source<K> input(dataset_filename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            assert(data[0] == 1);
            universe_size = data[1];
            input.close();
            dataset_filename_debug = dataset_filename;

            std::vector<std::vector<uint32_t>> query_list = read_query_random(query_filename);
            input_basename = input_filename;
            if (query_type == "AND")
                query_test_intersection(query_list, decode_type);
            else if (query_type == "OR")
                query_test_union(query_list, decode_type);
            else if (query_type == "RND")
                // random_access_debug(query_list);
                query_test_random_access(query_list, decode_type);
            else {
                std::cerr << "Error: query_type must be AND, OR, ANDB or ORB" << std::endl;
            }
        }
    };
}
