# pragma once
#include <config.hpp>
#if RESIDUAL_COMPRESS
#include <vector>
#include <tools.hpp>
#include <codecfactory.h>
#include <deltautil.h>

static FastPForLib::CODECFactory factory;
static FastPForLib::IntegerCODEC &codec = *factory.getFromName("simdfastpfor256");

static std::vector<uint32_t> compress_residuals_fastpfor(const std::vector<int32_t>& ints) {
    std::vector<uint32_t> zigzag_uints;
    zigzag_uints.reserve(ints.size());
    for (size_t i = 0; i < ints.size(); ++i)
        zigzag_uints.push_back(zigzag(ints[i]));

    std::vector<uint32_t> compressed_output(ints.size() + 1024);
    size_t compressed_size = compressed_output.size();

    codec.encodeArray(zigzag_uints.data(), zigzag_uints.size(), compressed_output.data(), compressed_size);
    compressed_output.resize(compressed_size);
    compressed_output.shrink_to_fit();
    return compressed_output;
}

static void decompress_residuals_fastpfor(std::vector<uint32_t>& uints, uint64_t n, std::vector<uint32_t> &uncompressed_output, std::vector<int32_t> &unzigzag_ints) {
    size_t uncompressed_size = uncompressed_output.size();
    codec.decodeArray(uints.data(), uints.size(), uncompressed_output.data(), uncompressed_size);
    assert(n == uncompressed_size);
    for (size_t i = 0; i < n; ++i)
        unzigzag_ints[i] = unzigzag(uncompressed_output[i]);

}

static void decompress_residuals_fastpfor_parted(std::vector<std::vector<uint32_t>>& uints, std::vector<uint32_t> &uncompressed_output, std::vector<uint64_t> &parted_sizes) {
    uint64_t decoded_size = 0;
    for (uint32_t i = 0; i < parted_sizes.size(); ++i) {
        size_t parted_uncompressed_size = parted_sizes[i];
        codec.decodeArray(uints[i].data(), uints[i].size(), uncompressed_output.data() + decoded_size, parted_uncompressed_size);
        decoded_size += parted_uncompressed_size;
    }
}
#endif
