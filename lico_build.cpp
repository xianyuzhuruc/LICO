#include <config.hpp>
#include <lico_build.hpp>

std::string log_filename = "";
std::string decode_type = "";

uint32_t lico_m = 0;
bool read_only = false;

template <uint64_t epsilon=64>
void build_collection_lico(const std::string input_basename, const std::string output_basename) {

    typedef lico_sequence::lico_builder<uint32_t, epsilon> PGM_INDEX_BUILDER;
    PGM_INDEX_BUILDER index;
    if (!read_only){
        index.build_model(input_basename + ".docs", lico_m);
        index.data_test(input_basename + ".docs");
        if(!log_filename.empty()) { // save covered and residual
            index.statistic_index(log_filename);
        }
        index.save_model(output_basename);

    } else {
        index.load_model(output_basename);
        // index.data_test(input_basename + ".docs");

        if(!log_filename.empty()) { // save covered and residual
            index.statistic_index(log_filename);
        }
    }
}

int main(int argc, const char** argv)
{
    std::ios::sync_with_stdio(0);
    int mandatory = 6;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <index_type> <collection_basename> <output_basename> <epsilon> <read_only> <decode_type> <test_log> <LICO_m> <residual_compress_type>" << std::endl;
        return 1;
    }

    const std::string index_type = argv[1];
    const std::string input_basename = argv[2];
    const std::string output_basename = argv[3];
    const std::string epsilonstr = argv[4];
    const std::string read_only_str = argv[5];
    decode_type = argv[6];
    log_filename = argv[7];
    const std::string lico_m_str = argv[8];
    const uint64_t epsilon = static_cast<uint64_t>(std::stoi(epsilonstr));
    lico_m = static_cast<uint32_t>(std::stoi(lico_m_str));
    residual_compress_type = argv[9];

    read_only = (read_only_str == "t");

    if (index_type == "lico" || index_type != "lico_swing" || index_type != "lico_greedy")  //compress docs
        switch (epsilon)
        {
            case 0: build_collection_lico<0>(input_basename, output_basename); break; // represents adopting data partition method
            case 1: build_collection_lico<1>(input_basename, output_basename); break;
            case 3: build_collection_lico<3>(input_basename, output_basename); break;
            case 7: build_collection_lico<7>(input_basename, output_basename); break;
            case 15: build_collection_lico<15>(input_basename, output_basename); break;
            case 31: build_collection_lico<31>(input_basename, output_basename); break;
            case 63: build_collection_lico<63>(input_basename, output_basename); break;
            case 127: build_collection_lico<127>(input_basename, output_basename); break;
            case 255: build_collection_lico<255>(input_basename, output_basename); break;
            case 511: build_collection_lico<511>(input_basename, output_basename); break;
            case 1023: build_collection_lico<1023>(input_basename, output_basename); break;
            case 2047: build_collection_lico<2047>(input_basename, output_basename); break;
            case 4095: build_collection_lico<4095>(input_basename, output_basename); break;
            case 8191: build_collection_lico<8191>(input_basename, output_basename); break;
            case 16383: build_collection_lico<16383>(input_basename, output_basename); break;
            case 32767: build_collection_lico<32767>(input_basename, output_basename); break;
            case 65535: build_collection_lico<65535>(input_basename, output_basename); break;
            case 131071: build_collection_lico<131071>(input_basename, output_basename); break;
            case 262143: build_collection_lico<262143>(input_basename, output_basename); break;

            default: std::cerr << "Unsupported Epsilon Value: " << epsilon << std::endl; break;
        }
    else
        std::cerr << "ERROR: only support pgm index " << std::endl;

    return 0;
}
