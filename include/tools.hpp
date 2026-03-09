# pragma once
#include <config.hpp>
#include <algorithm>
#include <fstream>
#include <sstream>

# ifndef  RESIDUAL_COMPRESS
# define RESIDUAL_COMPRESS 1
# endif

template <typename K> std::vector<K> load_data(const std::string &filename) {
    /* Open file. */
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open())
        exit(EXIT_FAILURE);

    /* Read number of keys. */
    K n_keys;
    in.read(reinterpret_cast<char*>(&n_keys), sizeof(K));

    /* Initialize vector. */
    std::vector<K> data;
    data.resize(n_keys);

    /* Read keys. */
    in.read(reinterpret_cast<char*>(data.data()), n_keys * sizeof(K));  
    in.close();
    std::sort(data.begin(), data.end());
    return data;
}

std::vector<std::string> split_str(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<double> parseDoubles(const std::string &line, char delimiter) {
    std::vector<double> values;
    std::vector<std::string> tokens = split_str(line, delimiter);

    for (const auto &token : tokens) {
        try {
            values.push_back(std::stod(token));
        } catch (const std::exception &e) {
            std::cerr << "Error parsing token: " << token << " (" << e.what() << ")\n";
        }
    }
    return values;
}

// ---------- sign <-> zigzag ----------
inline static uint32_t zigzag(int32_t x) {
    return (static_cast<uint32_t>(x) << 1) ^ static_cast<uint32_t>(x >> 31);
}

inline static int32_t unzigzag(uint32_t u) {
    return static_cast<int32_t>((u >> 1) ^ -(static_cast<int32_t>(u & 1)));
}

inline static int32_t unzigzag_int32(int32_t u) {
    return ((u >> 1) ^ -(u & 1));
}

template<typename IndexT>
static void write_index_data(std::ofstream& out, const IndexT& index) {
    out.write(reinterpret_cast<const char*>(&index.Epsilon_Data), sizeof(index.Epsilon_Data));
    out.write(reinterpret_cast<const char*>(&index.n), sizeof(index.n));

    // segments metadata
    out.write(reinterpret_cast<const char*>(&index.segments_size), sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(&index.bpc_covered), sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&index.bpc_delta_y), sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&index.bpc_delta_x), sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&index.bpc_y_b), sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(&index.bpc_x_b), sizeof(uint8_t));

    // seg arrays
    auto write_vector = [&out](const auto& vec) {
        size_t size = vec.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        if (size > 0) {
            out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(decltype(vec[0])));
        }
    };
    write_vector(index.seg_covered_compress);
    write_vector(index.seg_delta_y_compress);
    write_vector(index.seg_delta_x_compress);
    write_vector(index.seg_y_b_compress);
    write_vector(index.seg_x_b_compress);

#if RESIDUAL_COMPRESS
    write_vector(index.corrections_compress);
#else
    write_vector(index.corrections_none);
    auto write_bit_vector = [&out](const auto& vec) {
        size_t size = vec.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        if (size > 0) {
            out.write(reinterpret_cast<const char*>(vec.data()), (size + 7) / 8);
        }
    };
    write_bit_vector(index.signs_none);
#endif
}

template<typename IndexT>
static void read_index_data(std::ifstream& in, IndexT& index) {
    in.read(reinterpret_cast<char*>(&index.n), sizeof(index.n));

    in.read(reinterpret_cast<char*>(&index.segments_size), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&index.bpc_covered), sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&index.bpc_delta_y), sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&index.bpc_delta_x), sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&index.bpc_y_b), sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&index.bpc_x_b), sizeof(uint8_t));

    // seg arrays
    auto read_vector = [&in](auto& vec) {
        size_t size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        if (size > 0) {
            vec.resize(size);
            in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(decltype(vec[0])));
        }
    };

    read_vector(index.seg_covered_compress);
    read_vector(index.seg_delta_y_compress);
    read_vector(index.seg_delta_x_compress);
    read_vector(index.seg_y_b_compress);
    read_vector(index.seg_x_b_compress);

#if RESIDUAL_COMPRESS
    read_vector(index.corrections_compress);
#else
    read_vector(index.corrections_none);
    auto read_bit_vector = [&in](auto& vec) {
        size_t size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        if (size > 0) {
            vec.resize(size);
            in.read(reinterpret_cast<char*>(vec.data()), (size + 7) / 8);
        }
    };
    read_bit_vector(index.signs_none);
#endif
}


inline static int64_t sign_bit(int64_t x) {
    return static_cast<int64_t> (static_cast<uint64_t>(x) >> 63);
}