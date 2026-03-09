#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <perf_event.hpp>
#include <lico_kernel.hpp>
#include <lico_enumerate.hpp>

namespace lico_sequence {
    template <typename K, size_t epsilon = 64> // K is uint32_t or uint64_t
    class lico_decoder{

    public:
        uint64_t data_size = 0;
        uint64_t total_list = 0;

        uint64_t total_decode_time = 0;
        uint64_t max_decode_time = 0;
        uint64_t min_decode_time = UINT64_MAX - 1;

        uint64_t total_calculated = 0;
        uint64_t total_calculated_add = 0;
        uint64_t total_conversion_time = 0;
        uint64_t total_unequal = 0;

#ifdef PERF_OUTPUT
        std::vector<std::string> perf_header;
        std::vector<std::vector<double>> perf_data;
#endif

        void decode_test(lico_enumerator <K> &index, std::string &decode_type, bool warm_up = true) {
            total_list++;

            if (decode_type == "simd") {
                index.simd_init();

#if USE_HUGEPAGE
                std::vector<K, HugePageAllocator<K>> result1(index.n);
#else
                std::vector<K> result1(index.n);
#endif
                if (warm_up) {
                    for (int warm_time = 0; warm_time < 5; warm_time++) {
                        for (int i = 0; i < index.n; i ++) {
                            result1[i] = 0;
                        }
                    }
                }

#ifdef PERF_OUTPUT
                PerfEvent perf_event;
                perf_event.startCounters();
#endif

#if SIMD_512_1D1S
                index.simd_decode_512i_1d1s(result1.data());
#elif SIMD_512_2D1S
                index.simd_decode_512i_2d1s(result1.data());
#elif SIMD_256_2D1S
                index.simd_decode_256i_2d1s(result1.data());
#elif SIMD_256_1D1S
                index.simd_decode_256i_1d1s(result1.data());
#endif

#ifdef PERF_OUTPUT
                perf_event.stopCounters();
                std::stringstream header_out;
                std::stringstream data_out;
                PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
                perf_event.printReport(header_out, data_out, 1);
                std::vector<double> perf_data_tmp = parseDoubles(data_out.str(), ',');
                for (int i = 0; i < perf_data_tmp.size(); i++) {
                    perf_data[i].push_back(perf_data_tmp[i]);
                }
#endif

                total_calculated += index.total_calculated;
                if (index.total_duration > max_decode_time)
                    max_decode_time = index.total_duration;
                if (index.total_duration < min_decode_time)
                    min_decode_time = index.total_duration;

#if USE_HUGEPAGE
                std::vector<K, HugePageAllocator<K>> ().swap(result1);
#else
                std::vector<K> ().swap(result1);
#endif
            }
            else if (decode_type == "normal") {
                index.residuals_decode();

#if USE_HUGEPAGE
                std::vector<K, HugePageAllocator<K>> result2(index.n);
#else
                std::vector<K> result2(index.n);
#endif

                if (warm_up) {
                    for (int warm_time = 0; warm_time < 5; warm_time++) {
                        for (int i = 0; i < index.n; i ++) {
                            result2[i] = 0;
                        }
                    }
                }

#ifdef PERF_OUTPUT
                PerfEvent perf_event;
                perf_event.startCounters();
#endif

                index.normal_decode(result2.data());

#ifdef PERF_OUTPUT
                perf_event.stopCounters();
                std::stringstream header_out;
                std::stringstream data_out;
                PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
                perf_event.printReport(header_out, data_out, 1);
                std::vector<double> perf_data_tmp = parseDoubles(data_out.str(), ',');
                for (int i = 0; i < perf_data_tmp.size(); i++) {
                    perf_data[i].push_back(perf_data_tmp[i]);
                }
#endif

#if USE_HUGEPAGE
                std::vector<K, HugePageAllocator<K>> ().swap(result2);
#else
                std::vector<K> ().swap(result2);
#endif

                if (index.total_duration > max_decode_time)
                    max_decode_time = index.total_duration;
                if (index.total_duration < min_decode_time)
                    min_decode_time = index.total_duration;
            }

            index.free_memory(decode_type);
            data_size += index.n;
            total_decode_time += index.total_duration;
        }

        void result_statistic(std::string &decode_type) {
            std::cerr << "Total list: " << total_list << std::endl;
            if (decode_type == "simd" || decode_type == "simd_simple") {
                std::cerr << "Decode time 1 simd, average: " << total_decode_time / total_list  << ", max: " << max_decode_time << ", min: " << min_decode_time << " microseconds" <<std::endl;
                std::cerr << "Total calculated: " << total_calculated << " (" << double(total_calculated) / data_size << ") , Total calculated add: " << total_calculated_add << ", Total Conversion time: " << total_conversion_time << ", Total Unequal: " << total_unequal << std::endl;
            }
            if (decode_type == "normal") {
                std::cerr << "Decode time 2 normal, average: " << total_decode_time / total_list  << ", max: " << max_decode_time << ", min: " << min_decode_time << " microseconds" <<std::endl;
            }

            std::cerr << "Decode per integer: " << static_cast<long double> (total_decode_time) / data_size * 1000.0 << " nanoseconds" << std::endl;

#ifdef PERF_OUTPUT
            std::cerr << "Performance: " << std::endl;
            for (const auto &header : perf_header) {
                std::cerr << std::setw(15) << header;
            }
            std::cerr << std::endl;

            for (const auto &data : perf_data) {
                double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
                std::cerr << std::setw(15) << std::fixed << std::setprecision(2) << mean;
            }
            std::cerr << std::endl;
#endif

            std::cerr << std::endl;
        }

        void test_model(const std::string input_basename, std::string decode_type) {
#ifdef PERF_OUTPUT
            // init perf_header
            PerfEvent perf_event;
            perf_event.startCounters();
#endif

            if (input_basename.empty() || input_basename.back() != '/') {
                std::cerr << "Error: input_basename must end with '/'" << std::endl;
                return;
            }
            std::cerr << "Load index from: " << input_basename << std::endl;

#ifdef PERF_OUTPUT
            perf_event.stopCounters();
            std::stringstream header_out;
            std::stringstream data_out;
            PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
            perf_event.printReport(header_out, data_out, 1);
            perf_header = split_str(header_out.str(), ',');
            perf_data.resize(perf_header.size());
#endif

            std::ifstream in_header(input_basename + "idx.size", std::ios::binary);
            if (!in_header) {
                std::cerr << "Error: Cannot open idx.size for reading." << std::endl;
                return;
            }


            uint64_t data_size_tmp = 0;
            in_header.read(reinterpret_cast<char*>(&data_size_tmp), sizeof(data_size_tmp));
            K index_num = 0;
            in_header.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));

            if (epsilon == 0) {
                for (K i = 0; i < index_num; ++i) {
                    size_t partition_size;
                    in_header.read(reinterpret_cast<char*>(&partition_size), sizeof(size_t));

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

                        LICOIndex index = lico::LICO<K>(Epsilon_Data);

                        read_index_data<LICOIndex>(in, index);

                        partition_index.push_back(std::move(index));
                        in.close();
                    }

                    lico_enumerator <K> enumerator_tmp = create_enumerator_from_indexes <K>(partition_index);

                    decode_test(enumerator_tmp, decode_type);
                }
            }
            else {
                for (K i = 0; i < index_num; ++i) {
                    std::string filename = input_basename + std::to_string(i) + ".idx";
                    std::ifstream in(filename, std::ios::binary);

                    uint32_t Epsilon_Data;

                    while (in.peek() != EOF) {
                        in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));

                        LICOIndex index = lico::LICO<K>(Epsilon_Data);

                        read_index_data<LICOIndex>(in, index);

                        lico_enumerator <K> enumerator_tmp = create_enumerator_from_single_index <K> (index);

                        decode_test(enumerator_tmp, decode_type);
                    }

                    in.close();
                }

            }

            in_header.close();
            result_statistic(decode_type);
        }

    };
}
