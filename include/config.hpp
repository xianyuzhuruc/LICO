#pragma once
#include <string>

#ifndef  RESIDUAL_COMPRESS
#define RESIDUAL_COMPRESS 0
#endif

#ifndef HUGEPAGE
#define USE_HUGEPAGE 0
#else
#define USE_HUGEPAGE 1
#endif

#ifndef SIMD_512_1D1S
#define SIMD_512_1D1S 0
#else
#define SIMD_512_1D1S 1
#endif

#ifndef SIMD_256_1D1S
#define SIMD_256_1D1S 0
#else
#define SIMD_256_1D1S 1
#endif

#ifndef SIMD_256_2D1S
#define SIMD_256_2D1S 0
#else
#define SIMD_256_2D1S 1
#endif

#ifndef SIMD_512_2D1S
#define SIMD_512_2D1S 0
#else
#define SIMD_512_2D1S 1
#endif



#ifndef SwingFilterPLA
    #define OptimalPLA
#endif

static constexpr size_t block_size_rice = 8192;

static std::string residual_compress_type = "fastpfor";

static bool use_huge_pages = true;

