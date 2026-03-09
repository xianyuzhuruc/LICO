#include <iostream>
#include <string>
#include <fstream>
#include <lico_decode.hpp>

std::string decode_type = "";

template <uint64_t epsilon=64>
void test_collection_lico(const std::string input_basename, const std::string output_basename) {

    typedef lico_sequence::lico_decoder<uint32_t, epsilon> PGM_INDEX_DECODER;
    PGM_INDEX_DECODER index;

    index.test_model(output_basename, decode_type);
}

int main(int argc, const char** argv)
{
   // cpu_set_t mask;
   // CPU_ZERO(&mask); // clear CPU mask
   // CPU_SET(16, &mask); // set CPU 16
   //
   // if (sched_setaffinity(0, sizeof(mask), &mask) == -1) { // 0 is the calling process
   //     perror("sched_setaffinity");
   //     return 1;
   // }

    int mandatory = 6;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <index_type> <collection_basename> <output_basename> <epsilon> <decode_type> <test_log> <residual_compress_type>" << std::endl;
        return 1;
    }

    const std::string index_type = argv[1];
    const std::string input_basename = argv[2];
    const std::string output_basename = argv[3];
    const std::string epsilonstr = argv[4];
    decode_type = argv[5];
    // residual_compress_type = argv[7];
    uint32_t epsilon = static_cast<uint32_t>(std::stoi(epsilonstr));

    if (index_type == "lico" || index_type == "lico_swing")
        switch (epsilon)
        {
            case 0: test_collection_lico<0> (input_basename, output_basename); break;
            case 1: test_collection_lico<1>(input_basename, output_basename); break;
            case 7: test_collection_lico<7>(input_basename, output_basename); break;
            case 15: test_collection_lico<15>(input_basename, output_basename); break;
            // case 16: test_collection_lico<16>(input_basename, output_basename); break;
            case 31: test_collection_lico<31>(input_basename, output_basename); break;
            // case 32: test_collection_lico<32>(input_basename, output_basename); break;
            case 63: test_collection_lico<63>(input_basename, output_basename); break;
            // case 64: test_collection_lico<64>(input_basename, output_basename); break;
            // case 126: test_collection_lico<126>(input_basename, output_basename); break;
            case 127: test_collection_lico<127>(input_basename, output_basename); break;
            // case 128: test_collection_lico<128>(input_basename, output_basename); break;
            case 255: test_collection_lico<255>(input_basename, output_basename); break;
            // case 256: test_collection_lico<256>(input_basename, output_basename); break;
            case 511: test_collection_lico<511>(input_basename, output_basename); break;
            // case 512: test_collection_lico<512>(input_basename, output_basename); break;
            case 1023: test_collection_lico<1023>(input_basename, output_basename); break;
            // case 1024: test_collection_lico<1024>(input_basename, output_basename); break;
            case 2047: test_collection_lico<2047>(input_basename, output_basename); break;
            // case 2048: test_collection_lico<2048>(input_basename, output_basename); break;
            case 4095: test_collection_lico<4095>(input_basename, output_basename); break;
            case 8191: test_collection_lico<8191>(input_basename, output_basename); break;
            case 16383: test_collection_lico<16383>(input_basename, output_basename); break;
            case 32767: test_collection_lico<32767>(input_basename, output_basename); break;
            case 65535: test_collection_lico<65535>(input_basename, output_basename); break;
            // case 131071: test_collection_lico<131071>(input_basename, output_basename); break;
            // case 262143: test_collection_lico<262143>(input_basename, output_basename); break;
            default: std::cerr << "Unsupported Epsilon Value: " << epsilon << std::endl; break;
        }
    else
        std::cerr << "ERROR: only support lico index " << std::endl;

    return 0;
}
