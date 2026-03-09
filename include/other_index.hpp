#pragma once
#include <iostream>
#include <vector>
#include <chrono>
#include <climits>
#include <fstream>
#include <sstream>
#include <../external/mm_file/include/mm_file/mm_file.hpp>
#include <codecfactory.h>

typedef uint32_t K;

namespace other_sequence {

    static size_t block_size = 128;
    static uint32_t block_size_u32 = 128;

    static FastPForLib::CODECFactory factory_;
    static FastPForLib::IntegerCODEC& codec_ = *factory_.getFromName("simdfastpfor256");

    class fastpfor_list {
    public:
        fastpfor_list() : compressed_data_(), n(0) {}

        void encode(std::vector<K>& data) {
            if (data.empty()) {
                compressed_data_.clear();
                block_headers_.clear();
                compress_size_ = 0;
                n = 0;
                return;
            }

            n = data.size();
            compressed_data_.resize(n + 4096);

            K* compressed_data_begin = compressed_data_.data();
            K* compressed_data_pointer = compressed_data_.data();

            const uint64_t block_num = (data.size() + block_size - 1) / block_size;
            block_headers_.resize(block_num * 2 - 1); // n block_max, n-1 block_endpoint

            // generate gap vector
            std::vector<K> docs_buf(block_size);
            K* docs_it = data.data();
            auto start = std::chrono::high_resolution_clock::now();
            uint32_t last_doc(-1);
            for (uint64_t b = 0; b < block_num; b++) {
                size_t cur_block_size = std::min<size_t>(block_size, n - b * block_size);
                size_t true_compress_size = cur_block_size + 2048;
                for (size_t i = 0; i < cur_block_size; ++i) {
                    K doc(*docs_it++);
                    docs_buf[i] = doc - last_doc - 1;
                    last_doc = doc;
                }
                block_headers_[b] = last_doc; // max_doc of each block

                codec_.encodeArray(docs_buf.data(), cur_block_size, compressed_data_pointer, true_compress_size);
                compressed_data_pointer += true_compress_size;

                if (b != block_num - 1) {
                    block_headers_[block_num + b] = compressed_data_pointer - compressed_data_begin; // endpoint
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            total_duration_ = duration.count();

            compressed_data_.resize(compressed_data_pointer - compressed_data_begin);
            compressed_data_.shrink_to_fit();
            compress_size_ = compressed_data_.size() + block_headers_.size();
        }

        void decode(K* output) {
            const uint64_t block_num = (n + block_size - 1) / block_size;
            assert(block_num > 0);

            uint32_t* compressed_data_begin = compressed_data_.data();
            uint32_t* block_max_begin = block_headers_.data();
            uint32_t* block_endpoint_begin = block_max_begin + block_num;


            uint32_t endpoint = 0;
            uint32_t i = 0;

            for (; i != block_num - 1; ++i) {
                uint32_t cur_base = (i ? block_max_begin[i - 1] : uint32_t(-1)) + 1;
                uint32_t const* ptr = compressed_data_begin + endpoint;
                size_t compressed_length = block_endpoint_begin[i] - endpoint;

                codec_.decodeArray(ptr, compressed_length, output, block_size);

                output[0] += cur_base;
                for (uint32_t k = 1; k < block_size_u32; ++k) {
                    output[k] += output[k - 1] + 1;
                }

                endpoint = block_endpoint_begin[i];
                output += block_size;
            }

            // last block
            uint32_t cur_base = (i ? block_max_begin[i - 1] : uint32_t(-1)) + 1;
            uint32_t const* ptr = compressed_data_begin + endpoint;
            size_t compressed_length = compressed_data_.size() - endpoint;
            size_t last_block_size = n - (block_num - 1) * block_size;
            codec_.decodeArray(ptr, compressed_length, output, last_block_size);

            output[0] += cur_base;
            for (uint32_t k = 1; k < last_block_size; ++k) {
                output[k] += output[k - 1] + 1;
            }
        }

        const std::vector<uint32_t>& get_compressed_data() const {
            return compressed_data_;
        }

        const std::vector<uint32_t>& get_header_data() const {
            return block_headers_;
        }

        void set_compressed_data(std::vector<uint32_t> &data, size_t original_size) {
            compressed_data_ = std::move(data);
            n = original_size;
        }

        void set_header_data(std::vector<uint32_t> &data) {
            block_headers_ = std::move(data);
            compress_size_ = compressed_data_.size() + block_headers_.size();
        }

        size_t size() const {
            return (compressed_data_.size() + block_headers_.size()) * sizeof(uint32_t);
        }


        size_t original_size() const {
            return n;
        }

        bool empty() const {
            return compressed_data_.empty();
        }

        size_t time_cost_ns() const {
            return total_duration_;
        }

        void next_geq_init() {
            total_duration_ = 0;
            total_skip_ = 0;
            block_buf.resize(block_size_u32);
            total_block_num = (n + block_size_u32 - 1) / block_size_u32;
            total_endpoint_begin = block_headers_.data() + total_block_num;

            // decode the first block
            current_block = 0;

            auto start = std::chrono::high_resolution_clock::now();

            const uint32_t cur_block_size = current_block < (total_block_num - 1) ? block_size_u32 : n - (total_block_num - 1) * block_size_u32;
            uint32_t* ptr = compressed_data_.data();
            current_value_pointer = block_buf.data();
            current_remain = cur_block_size;
            current_block_max = block_headers_[current_block];

            size_t decode_block_size = cur_block_size;
            codec_.decodeArray(ptr, total_endpoint_begin[0], current_value_pointer, decode_block_size);

            uint32_t cur_base = uint32_t(-1) + 1;
            current_value_pointer[0] += cur_base;
            for (uint32_t k = 1; k < cur_block_size; ++k) {
                current_value_pointer[k] += current_value_pointer[k - 1] + 1;
            }

            current_value = *current_value_pointer;

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration_ = duration.count();
        }

        uint32_t nextgeq(uint32_t posting_value) {
            if (posting_value <= current_value) return current_value;

            if (posting_value <= current_block_max) { // in current block
                while (current_remain > 0) {
                    if (*current_value_pointer >= posting_value) {
                        current_value = *current_value_pointer++;
                        current_remain--;
                        return current_value;
                    }
                    current_remain--;
                    current_value_pointer++;
                }
            } else {
                auto end   = block_headers_.begin() + total_block_num;
                auto it = std::lower_bound(block_headers_.begin() + current_block + 1, end, posting_value);
                if (it == end) {
                    current_value = INT_MAX;
                    return current_value;
                }
                total_skip_ += current_remain + (it - block_headers_.begin() - (current_block + 1)) * block_size;

                current_block = it - block_headers_.begin();
                // decode one block
                const uint32_t cur_block_size = current_block < (total_block_num - 1) ? block_size_u32 : n - (total_block_num - 1) * block_size_u32;
                uint32_t* ptr = compressed_data_.data() + total_endpoint_begin[current_block - 1];
                current_value_pointer = block_buf.data();
                current_remain = cur_block_size;
                current_block_max = block_headers_[current_block];

                size_t compressed_length = current_block < total_block_num - 1 ? total_endpoint_begin[current_block] - total_endpoint_begin[current_block - 1] : compressed_data_.size() - total_endpoint_begin[current_block - 1];

                size_t decode_block_size = cur_block_size;
                codec_.decodeArray(ptr, compressed_length, current_value_pointer, decode_block_size);

                uint32_t cur_base = block_headers_[current_block - 1] + 1;
                current_value_pointer[0] += cur_base;
                for (uint32_t k = 1; k < cur_block_size; ++k) {
                    current_value_pointer[k] += current_value_pointer[k - 1] + 1;
                }

                while (current_remain > 0) {
                    if (*current_value_pointer >= posting_value) {
                        current_value = *current_value_pointer++;
                        current_remain--;
                        return current_value;
                    }
                    current_remain--;
                    current_value_pointer++;
                }
            }
            current_value = INT_MAX;
            return INT_MAX;
        }

        void next_init() {
            total_duration_ = 0;
            total_skip_ = 0;
            block_buf.resize(block_size_u32);
            total_block_num = (n + block_size - 1) / block_size;
            total_endpoint_begin = block_headers_.data() + total_block_num;

            // decode the first block
            current_block = 0;

            auto start = std::chrono::high_resolution_clock::now();

            const uint32_t cur_block_size = current_block < (total_block_num - 1) ? block_size : n - (total_block_num - 1) * block_size;
            uint32_t* ptr = compressed_data_.data();
            current_value_pointer = block_buf.data();
            current_remain = cur_block_size - 1;

            size_t decode_block_size = cur_block_size;
            codec_.decodeArray(ptr, total_endpoint_begin[0], current_value_pointer, decode_block_size);

            uint32_t cur_base = uint32_t(-1) + 1;
            current_value_pointer[0] += cur_base;
            for (uint32_t k = 1; k < cur_block_size; ++k) {
                current_value_pointer[k] += current_value_pointer[k - 1] + 1;
            }

            current_value = *current_value_pointer;

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration_ = duration.count();
        }

        void next() {
            if (current_remain > 0) {
                current_value = *(++current_value_pointer);
                current_remain--;
            } else {
                // decode next block
                current_block++;
                if (current_block < total_block_num) {
                    const uint32_t cur_block_size = current_block < (total_block_num - 1) ? block_size : n - (total_block_num - 1) * block_size;
                    uint32_t* ptr = compressed_data_.data() + total_endpoint_begin[current_block - 1];
                    current_value_pointer = block_buf.data();
                    size_t compressed_length = current_block < total_block_num - 1 ? total_endpoint_begin[current_block] - total_endpoint_begin[current_block - 1] : compressed_data_.size() - total_endpoint_begin[current_block - 1];

                    size_t decode_block_size = cur_block_size;
                    codec_.decodeArray(ptr, compressed_length, current_value_pointer, decode_block_size);

                    uint32_t cur_base = block_headers_[current_block - 1] + 1;
                    current_value_pointer[0] += cur_base;
                    for (uint32_t k = 1; k < cur_block_size; ++k) {
                        current_value_pointer[k] += current_value_pointer[k - 1] + 1;
                    }

                    current_remain = cur_block_size - 1;
                    current_value = *current_value_pointer;
                } else {
                    current_remain = 0;
                    current_value = INT_MAX;
                }
            }
        }

        uint32_t docid() {
            return current_value;
        }



    public:
        uint32_t current_value;
        uint32_t current_block;
        uint32_t current_block_max;
        uint32_t current_remain;
        uint32_t total_block_num;
        uint32_t* current_value_pointer;
        uint32_t* total_endpoint_begin;
        std::vector<uint32_t> block_buf;

        size_t total_skip_;

        std::vector<uint32_t> compressed_data_;
        std::vector<uint32_t> block_headers_;
        size_t n;
        size_t compress_size_;
        size_t total_duration_ = 0;


    };

    // builder
    class other_index_build {
    private:
        uint64_t data_size = 0;
        uint64_t data_compressed_size = 0;

        std::vector<fastpfor_list> encoded_lists;
        std::vector<uint64_t> list_lengths;

        bool use_u8 = false;

    public:
        void build_index(const std::string& input_basename, const std::string& index_type) {
            std::cerr << "Read File [Build]: " << input_basename << std::endl;
            mm::file_source<uint32_t> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            std::cerr << "Universe Size: " << data[1] << std::endl;

            K list_count = 0;
            for (size_t i = 2; i < input.size();) {
                uint32_t n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n + 1);
                data_size += n;
                list_lengths.push_back(n);
                list_count++;
                if (list_count % 100 == 0)
                    std::cerr << list_count << " ";
                if (index_type == "fastpfor") {
                    fastpfor_list list_encoder;
                    list_encoder.encode(sequence);
                    encoded_lists.push_back(std::move(list_encoder));
                }

                i += n + 1;
            }
        }

        void statistic_index(std::string output_basename = "") {
            if (output_basename.empty()) {
                std::cerr << "Warning: No output basename provided. Using stdout for statistics." << std::endl;
            }

            std::ofstream file;
            if (!output_basename.empty()) {
                file.open(output_basename + ".statistic_log.txt");
                if (!file.is_open()) {
                    std::cerr << "Error: Cannot open " << output_basename << ".statistic_log.txt for writing." << std::endl;
                    return;
                }
                file << "#list_id,original_elements,original_bytes,compressed_bytes,compression_ratio(%)\n";
            }

            uint64_t total_original_bytes = 0;
            uint64_t total_compressed_bytes = 0;
            uint64_t total_duration = 0;

            for (size_t i = 0; i < encoded_lists.size(); ++i) {
                const auto& list_encoder = encoded_lists[i];
                const size_t original_size = list_encoder.original_size();
                const size_t compressed_bytes = list_encoder.size();

                const size_t original_bytes = original_size * sizeof(K);

                total_original_bytes += original_bytes;
                total_compressed_bytes += compressed_bytes;
                total_duration += list_encoder.time_cost_ns();

                double compression_ratio = 0.0;
                if (original_bytes > 0) {
                    compression_ratio = static_cast<double>(original_bytes) / compressed_bytes * 100.0;
                }

                if (file.is_open()) {
                    file << i << " "
                         << original_size << " "
                         << original_bytes << " "
                         << compressed_bytes << " "
                         << compression_ratio << "%\n";
                }
            }

            double global_compression_ratio = 0.0;
            double bits_per_integer = 0.0;

            if (total_original_bytes > 0) {
                global_compression_ratio = static_cast<double>(total_original_bytes) /  total_compressed_bytes * 100.0;
            }

            if (data_size > 0) {
                bits_per_integer = static_cast<double>(total_compressed_bytes * 8) / data_size;
            }

            std::ostream& out = file.is_open() ? static_cast<std::ostream&>(file) : std::cout;

            out << "\n# Global Statistics\n";
            out << "Total original size (elements): " << data_size << "\n";
            out << "Total original bytes: " << total_original_bytes << " Bytes\n";
            out << "Total compressed bytes: " << total_compressed_bytes << " Bytes, " << static_cast<long double> (total_compressed_bytes) / 1024 / 1024 / 1024 <<  " GiB \n";
            out << "Global compression ratio: " << global_compression_ratio << "%\n";
            out << "Compression factor (original/compressed): " << (total_compressed_bytes > 0 ? static_cast<double>(total_original_bytes) / total_compressed_bytes : 0.0) << "%\n";
            out << "Average bits per integer (BPI): " << bits_per_integer << " bits\n";

            std::cerr << "\nBuild Time: " << static_cast<long double> (total_duration) / 1e9 << "\n";
            std::cerr << "\nTotal compressed bytes: " << total_compressed_bytes << " Bytes, " << static_cast<long double> (total_compressed_bytes) / 1024 / 1024 / 1024 << " GiB, Global compression ratio: " << global_compression_ratio << "%, " << "Average bits per integer (BPI): " << bits_per_integer << " bits\n";

            if (file.is_open()) {
                file.close();
                std::cerr << "Statistics saved to: " << output_basename << ".statistic_log.txt" << std::endl;
            }
        }

        void data_test(const std::string& input_basename, const std::string& index_type) {
            std::cerr << "Read File [Test]: " << input_basename << std::endl;
            mm::file_source<uint32_t> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            std::cerr << "Universe Size: " << data[1] << std::endl;

            size_t data_unequal_test = 0;
            K idx = 0;
            for (size_t i = 2; i < input.size();) {
                uint32_t n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n + 1);

                std::vector<K> result(sequence.size());

                encoded_lists[idx].decode(result.data());

                assert(result.size() == sequence.size());

                for (auto j = 0; j < result.size(); j++) {
                    if (sequence[j] != result[j]) {
                        data_unequal_test++;
                        // std::cerr << "Unequal Value: " << result_decode[j] << " " << sequence[j] << " " << posi << " " << j << std::endl;
                    }
                }

                i += n + 1;
                idx++;
            }
            std::cerr << "Unequal postings: " << data_unequal_test << std::endl << std::endl;
        }

        void decode_test() {
            uint64_t total_duration = 0;
            uint64_t total_integers = 0;
            for (uint32_t i = 0; i < encoded_lists.size(); i++) {
                auto &list = encoded_lists[i];
                total_integers += list.original_size();
                std::vector<K> result(list.original_size());
                for (int32_t j = 0; j < 3; j++)
                    for (int32_t k = 0; k < list.original_size(); k++)
                        result[k] = j + k;

                auto start = std::chrono::high_resolution_clock::now();
                list.decode(result.data());
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                total_duration += duration.count();
            }
            std::cerr << "Total duration: " << total_duration << std::endl;
            std::cerr << "Decode per int: " << static_cast<long double> (total_duration) / total_integers << " ns" << std::endl;
        }

        uint64_t index_size() {
            data_compressed_size = 0;
            for (const auto& lst : encoded_lists) {
                data_compressed_size += lst.size();
            }
            return data_compressed_size;
        }

        void save_model(const std::string& output_basename) {
            if (output_basename.empty() || output_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                return;
            }

            std::cerr << "Save index to: " << output_basename << std::endl;

            std::ofstream out_header(output_basename + "idx.size", std::ios::binary);
            if (!out_header) {
                std::cerr << "Error: Cannot open idx.size for writing." << std::endl;
                return;
            }

            out_header.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));
            use_u8 = false;
            out_header.write(reinterpret_cast<const char*>(&use_u8), sizeof(use_u8));

            size_t total_list_size = encoded_lists.size();
            out_header.write(reinterpret_cast<const char*>(&total_list_size), sizeof(size_t));

            for (size_t i = 0; i < total_list_size; ++i) {
                const auto& compressed_data = encoded_lists[i].get_compressed_data();
                const auto& block_header = encoded_lists[i].get_header_data();
                size_t comp_size = compressed_data.size();
                size_t header_size = block_header.size();
                size_t orig_size = encoded_lists[i].original_size();

                std::string filename = output_basename + std::to_string(i) + ".idx";
                std::ofstream out(filename, std::ios::binary);
                if (!out) {
                    std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
                    throw std::runtime_error("file open error");
                }

                out.write(reinterpret_cast<const char*>(&orig_size), sizeof(orig_size));
                out.write(reinterpret_cast<const char*>(&comp_size), sizeof(comp_size));
                out.write(reinterpret_cast<const char*>(&header_size), sizeof(header_size));
                out.write(reinterpret_cast<const char*>(compressed_data.data()), comp_size * sizeof(uint32_t));
                out.write(reinterpret_cast<const char*>(block_header.data()), header_size * sizeof(uint32_t));
                out.close();
            }
            out_header.close();
        }

        void load_model(const std::string& input_basename) {
            if (input_basename.empty() || input_basename.back() != '/') {
                std::cerr << "Error: input_basename must end with '/'" << std::endl;
                return;
            }

            std::cerr << "Load index from: " << input_basename << std::endl;

            std::ifstream in_header(input_basename + "idx.size", std::ios::binary);
            if (!in_header) {
                std::cerr << "Error: Cannot open idx.size for reading." << std::endl;
                return;
            }

            in_header.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            in_header.read(reinterpret_cast<char*>(&use_u8), sizeof(use_u8));

            size_t total_list_size;
            in_header.read(reinterpret_cast<char*>(&total_list_size), sizeof(size_t));
            in_header.close();

            encoded_lists.clear();
            encoded_lists.reserve(total_list_size);
            list_lengths.clear();
            list_lengths.reserve(total_list_size);

            for (size_t i = 0; i < total_list_size; ++i) {
                std::string filename = input_basename + std::to_string(i) + ".idx";
                std::ifstream in(filename, std::ios::binary);
                if (!in) {
                    std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                    continue;
                }

                size_t orig_size, comp_size, header_size;
                in.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
                in.read(reinterpret_cast<char*>(&comp_size), sizeof(comp_size));
                in.read(reinterpret_cast<char*>(&header_size), sizeof(header_size));

                std::vector<uint32_t> compressed_data(comp_size);
                in.read(reinterpret_cast<char*>(compressed_data.data()), comp_size * sizeof(uint32_t));
                std::vector<uint32_t> block_header(header_size);
                in.read(reinterpret_cast<char*>(block_header.data()), header_size * sizeof(uint32_t));
                in.close();

                fastpfor_list list;
                list.set_compressed_data(compressed_data, orig_size);
                list.set_header_data(block_header);
                encoded_lists.push_back(std::move(list));
                list_lengths.push_back(orig_size);
            }
        }

        void decode_list(size_t idx, K* output) {
            if (idx < encoded_lists.size()) {
                encoded_lists[idx].decode(output);
            }
        }
    };

    // decode
    class other_index_decode {
    private:
        uint64_t data_size = 0;
        uint64_t total_duration = 0;
        uint64_t total_integers = 0;

        std::vector<uint64_t> list_lengths;
        bool use_u8 = false;

    public:
        void decode_test(fastpfor_list &list) {
            total_integers += list.original_size();
            std::vector<K> result(list.original_size());
            for (int32_t j = 0; j < 3; j++)
                for (int32_t k = 0; k < list.original_size(); k++)
                    result[k] = j + k;

            auto start = std::chrono::high_resolution_clock::now();
            list.decode(result.data());
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            total_duration += duration.count();
        }

        void result_statistics() {
            std::cerr << "Total duration: " << total_duration << std::endl;
            std::cerr << "Decode per int: " << static_cast<long double> (total_duration) / total_integers << " ns" << std::endl;
        }

        void load_model(const std::string& input_basename) {
            if (input_basename.empty() || input_basename.back() != '/') {
                std::cerr << "Error: input_basename must end with '/'" << std::endl;
                return;
            }

            std::cerr << "Load index from: " << input_basename << std::endl;

            std::ifstream in_header(input_basename + "idx.size", std::ios::binary);
            if (!in_header) {
                std::cerr << "Error: Cannot open idx.size for reading." << std::endl;
                return;
            }

            in_header.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            in_header.read(reinterpret_cast<char*>(&use_u8), sizeof(use_u8));

            size_t total_list_size;
            in_header.read(reinterpret_cast<char*>(&total_list_size), sizeof(size_t));
            in_header.close();


            for (size_t i = 0; i < total_list_size; ++i) {
                std::string filename = input_basename + std::to_string(i) + ".idx";
                std::ifstream in(filename, std::ios::binary);
                if (!in) {
                    std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                    continue;
                }

                size_t orig_size, comp_size, header_size;
                in.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
                in.read(reinterpret_cast<char*>(&comp_size), sizeof(comp_size));
                in.read(reinterpret_cast<char*>(&header_size), sizeof(header_size));

                std::vector<uint32_t> compressed_data(comp_size);
                in.read(reinterpret_cast<char*>(compressed_data.data()), comp_size * sizeof(uint32_t));
                std::vector<uint32_t> block_header(header_size);
                in.read(reinterpret_cast<char*>(block_header.data()), header_size * sizeof(uint32_t));
                in.close();

                fastpfor_list list;
                list.set_compressed_data(compressed_data, orig_size);
                list.set_header_data(block_header);
                decode_test(list);
            }
        }
    };

    class other_index_query {
    private:
        uint64_t data_size = 0;
        uint64_t total_duration = 0;

        std::vector<uint64_t> list_lengths;

        bool use_u8 = false;

        long double avg_skip = 0;
        long double avg_query_total_size = 0;
        long double avg_query_real_size = 0;
        std::string input_basename;

    public:
        static void remove_duplicate_terms(std::vector<uint32_t>& terms) {
            std::sort(terms.begin(), terms.end());
            terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
        }

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
            uint32_t query_num = idLists[0].size();
            std::cout << "Total query sequences: " << idLists.size() << " Query num: " << (query_num > 5 ? 5 : query_num) << std::endl;
            file.close();
            return idLists;
        }

        static void intersection_candidate_test(std::vector<fastpfor_list> &index_sequences, K* intersection_result_p1, uint32_t &equal_result, uint32_t &query_id_idx, uint32_t &candidate_posting_tmp, const uint32_t m) {
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

        void query_test_intersection_benchmark(const std::vector<std::vector<uint32_t>> &query_list) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_per = 0;
            uint64_t avg_time_round = 0;
            K repeat_num = 5;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                int32_t query_id = 0;
                avg_time_round = 0;
                for (auto &query : query_list) {
                    query_id++;
                    if (query.size() < 2)
                        continue;

                    std::vector<fastpfor_list> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const fastpfor_list &a, const fastpfor_list &b) {return a.n < b.n;});

                    uint32_t query_id_idx = 0;
                    uint32_t equal_result = 0;
                    uint32_t candidate_posting_tmp = 0;
                    const uint32_t m = index_sequences.size();
                    int32_t intersection_size = query.size() == 2 ? index_sequences[1].n + 1 : index_sequences[2].n + 1;

                    K *intersection_result_p1;
                    // std::vector<K, HugePageAllocator<K>> intersection_result_1(intersection_size); // for simd
                    // std::vector<K, HugePageAllocator<K>> intersection_result_2(intersection_size);
                    std::vector<K> intersection_result_1(intersection_size); // for normal
                    // std::vector<K> intersection_result_2(intersection_size);
                    intersection_result_p1 = intersection_result_1.data();
                    // intersection_result_p2 = intersection_result_2.data();

                    // warm up
                    for (auto k = 0; k < 5; k++) {
                        for (auto i = 0;i < intersection_size; i++)
                            intersection_result_p1[i] = i + 1 + k;
                    }

                    avg_time_per = 0;

                    for (int32_t i = 0; i < index_sequences.size(); i++) {
                        index_sequences[i].next_geq_init();
                        avg_time_per += index_sequences[i].total_duration_;
                    }

                    auto start = std::chrono::high_resolution_clock::now();

                    intersection_candidate_test(index_sequences, intersection_result_p1, equal_result, query_id_idx, candidate_posting_tmp, m);
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                    avg_time_per += duration.count();

                    if (repeat == 0) {
                        total += equal_result;
                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip_;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;

                    }

                    avg_time_round += avg_time_per;

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

        static void union_candidate_test(std::vector<fastpfor_list> &index_sequences, K* result_p1, uint32_t &result, const uint32_t m) {
            uint32_t cur_doc = std::min_element(index_sequences.begin(), index_sequences.end(), [](fastpfor_list lhs, fastpfor_list rhs) {return lhs.docid() < rhs.docid();})->docid();
            K* result_start = result_p1;

            while (cur_doc < INT_MAX) {
                *result_p1++ = cur_doc;
                uint32_t next_doc = INT_MAX;
                for (int32_t i = 0; i < m; i++) {
                    if (index_sequences[i].docid() == cur_doc) {
                        index_sequences[i].next();
                    }
                    if (index_sequences[i].docid() < next_doc) {
                        next_doc = index_sequences[i].docid();
                    }
                }
                cur_doc = next_doc;
            }

            result = result_p1 - result_start;
        }

        void query_test_union_benchmark(const std::vector<std::vector<uint32_t>> &query_list) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_per = 0;
            uint64_t avg_time_round = 0;
            K repeat_num = 5;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                int32_t query_id = 0;
                avg_time_round = 0;
                for (auto &query : query_list) {
                    query_id++;
                    // std::cerr << query_id << " ";
                    if (query.size() < 2)
                        continue;

                    std::vector<fastpfor_list> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const fastpfor_list &a, const fastpfor_list &b) {return a.n < b.n;});


                    uint32_t equal_result = 0;
                    const uint32_t m = index_sequences.size();

                    K *union_result_p1;
                    std::vector<K> union_result_1(universe_size); // for normal
                    union_result_p1 = union_result_1.data();

                    // warm up
                    for (auto k = 0; k < 5; k++) {
                        for (auto i = 0;i < universe_size; i++)
                            union_result_p1[i] = i + 1 + k;
                    }

                    avg_time_per = 0;

                    for (int32_t i = 0; i < index_sequences.size(); i++) {
                        index_sequences[i].next_init();
                        avg_time_per += index_sequences[i].total_duration_;
                    }

                    auto start = std::chrono::high_resolution_clock::now();

                    union_candidate_test(index_sequences, union_result_p1, equal_result, m);

                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                    avg_time_per += duration.count();

                    if (repeat == 0) {
                        total += equal_result;
                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip_;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;

                    }

                    avg_time_round += avg_time_per;

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

        K universe_size;
        void test_query(const std::string &input_filename, const std::string &query_filename, const std::string query_type, const std::string &dataset_filename) {
            mm::file_source<K> input(dataset_filename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            assert(data[0] == 1);
            universe_size = data[1];
            input.close();

            std::vector<std::vector<uint32_t>> query_list = read_query(query_filename);
            input_basename = input_filename;
            if (query_type == "ANDB")
                query_test_intersection_benchmark(query_list);
            else if (query_type == "ORB")
                query_test_union_benchmark(query_list);
            else {
                std::cerr << "Error: query_type must be ANDB or ORB" << std::endl;
            }
        }

        std::vector<fastpfor_list> load_model(const std::vector<uint32_t> &idx_list) {
            if (input_basename.empty() || input_basename.back() != '/') {
                std::cerr << "Error: input_basename must end with '/'" << std::endl;
                return {};
            }

            std::ifstream in_header(input_basename + "idx.size", std::ios::binary);
            if (!in_header) {
                std::cerr << "Error: Cannot open idx.size for reading." << std::endl;
                return {};
            }

            in_header.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            in_header.read(reinterpret_cast<char*>(&use_u8), sizeof(use_u8));

            size_t total_list_size;
            in_header.read(reinterpret_cast<char*>(&total_list_size), sizeof(size_t));
            in_header.close();

            std::vector<fastpfor_list> result_list;
            result_list.reserve(idx_list.size());
            for (auto &i : idx_list) {
                std::string filename = input_basename + std::to_string(i) + ".idx";
                std::ifstream in(filename, std::ios::binary);
                if (!in) {
                    std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                    continue;
                }

                size_t orig_size, comp_size, header_size;
                in.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
                in.read(reinterpret_cast<char*>(&comp_size), sizeof(comp_size));
                in.read(reinterpret_cast<char*>(&header_size), sizeof(header_size));

                std::vector<uint32_t> compressed_data(comp_size);
                in.read(reinterpret_cast<char*>(compressed_data.data()), comp_size * sizeof(uint32_t));
                std::vector<uint32_t> block_header(header_size);
                in.read(reinterpret_cast<char*>(block_header.data()), header_size * sizeof(uint32_t));
                in.close();

                fastpfor_list list;
                list.set_compressed_data(compressed_data, orig_size);
                list.set_header_data(block_header);
                result_list.push_back(list);
            }

            return result_list;
        }
    };
} // namespace other_sequence