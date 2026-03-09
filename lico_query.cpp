#include <lico_query.hpp>

std::string log_filename = "";
std::string decode_type = "";
std::string query_filename = "";
std::string query_type = "";

template <uint64_t epsilon=64>
void test_collection_lico(const std::string input_basename, const std::string output_filename) {

    typedef lico_sequence::lico_querier<uint32_t, epsilon> PGM_INDEX_QUERIER;
    PGM_INDEX_QUERIER index;
    index.test_query(output_filename, decode_type, query_filename, query_type, input_basename + ".docs");

}

int main(int argc, const char** argv)
{
//    cpu_set_t mask;
//    CPU_ZERO(&mask); // clear CPU mask
//    CPU_SET(16, &mask); // set CPU 16
//
//    if (sched_setaffinity(0, sizeof(mask), &mask) == -1) { // 0 is the calling process
//        perror("sched_setaffinity");
//        return 1;
//    }

    int mandatory = 8;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <index_type> <collection_basename> <index_filename> <epsilon> <decode_type> <query_filename> <query_type> <log_filename> <residual_compress_type>" << std::endl;
        return 1;
    }

    const std::string index_type = argv[1];
    const std::string input_basename = argv[2];
    const std::string output_basename = argv[3];
    const std::string epsilon_str = argv[4];
    uint32_t epsilon = static_cast<uint32_t>(std::stoi(epsilon_str));
    decode_type = argv[5];
    query_filename = argv[6];
    query_type = argv[7];
    log_filename = argv[8];
    residual_compress_type = argv[9];

    if (index_type == "lico")
        switch (epsilon)
        {
            case 0: test_collection_lico<0>(input_basename, output_basename); break; // represents adopting data partition method
            case 1: test_collection_lico<1>(input_basename, output_basename); break;
            case 7: test_collection_lico<7>(input_basename, output_basename); break;
            case 15: test_collection_lico<15>(input_basename, output_basename); break;
            case 31: test_collection_lico<31>(input_basename, output_basename); break;
            case 63: test_collection_lico<63>(input_basename, output_basename); break;
            case 127: test_collection_lico<127>(input_basename, output_basename); break;
            case 255: test_collection_lico<255>(input_basename, output_basename); break;
            case 511: test_collection_lico<511>(input_basename, output_basename); break;
            case 1023: test_collection_lico<1023>(input_basename, output_basename); break;
            case 2047: test_collection_lico<2047>(input_basename, output_basename); break;
            case 4095: test_collection_lico<4095>(input_basename, output_basename); break;
            case 8191: test_collection_lico<8191>(input_basename, output_basename); break;
            case 16383: test_collection_lico<16383>(input_basename, output_basename); break;
            case 32767: test_collection_lico<32767>(input_basename, output_basename); break;
            case 65535: test_collection_lico<65535>(input_basename, output_basename); break;
            case 131071: test_collection_lico<131071>(input_basename, output_basename); break;
            case 262143: test_collection_lico<262143>(input_basename, output_basename); break;
            default: std::cerr << "Unsupported Epsilon Value: " << epsilon << std::endl; break;
        }
    else
        std::cerr << "ERROR: only support pgm index " << std::endl;
    return 0;
}
