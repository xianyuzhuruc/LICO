#pragma once
#include <vector>
#include <queue>
#include <tuple>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <climits>
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include <fstream>
#include <limits>
#include <cstring>

#if RESIDUAL_COMPRESS
static inline long double LICO_S_scale = 136*200L; // scaling for S in bits cost
#else
static inline long double LICO_S_scale = 136L; // scaling for S in bits cost
#endif
static size_t page_size;    // keys per page

#ifndef MIN_PAGES_PER_BLOCK
#define MIN_PAGES_PER_BLOCK 1
#endif

using gap_prefix_type = __uint128_t;

// Half-open index ranges everywhere: [start_idx, end_idx)
struct Page {
    size_t start_idx = 0; // inclusive key index
    size_t end_idx   = 0; // exclusive key index
};

struct Block {
    size_t start_idx = 0; // inclusive key index
    size_t end_idx   = 0; // exclusive key index
    size_t epsilon   = 1; // epsilon*
    long double cost = 0.0L; // bits at epsilon*
};

struct BlockStats {
    size_t n_keys = 0;  // t - s
    size_t n_gaps = 0;  // n_keys - 1
    long double mu  = 0.0L; // mean gap
    long double var = 0.0L; // unbiased gap variance
};

static inline long double safe_log2l(long double x) {
    if (!(x > 0.0L)) x = LDBL_MIN; // guard
    // log2l is preferred, but log2 works fine with long double on many libcs
    return std::log2(x);
}

static bool output_log = true;

class lico_partition {
public:
    using K = uint32_t;
    size_t global_m = 4096;       // target max blocks for greedy
#ifndef SYN_DATA
    long double LICO_C_coef = 0.109967L;
        // long double LICO_C_coef = 0.0000003520L;
//     Coverage (%) | Slope C | Intercept b |    R^2
// -------------|---------|-------------|----------
//           10 | 2.8227737825 |      114.62 |   0.5164
//           20 | 1.3718942841 |      207.17 |   0.6195
//           30 | 0.8261582635 |      335.99 |   0.7402
//           40 | 0.5118288497 |      620.18 |   0.8397
//           50 | 0.3474892392 |     1251.43 |   0.9138
//           60 | 0.2491987578 |     2718.44 |   0.9467
//           70 | 0.1845336882 |     5320.10 |   0.9717
//           80 | 0.1090760007 |    13849.29 |   0.9231
//           90 | 0.0448444793 |    38532.45 |   0.7585
//          100 | 0.0000003520 |   145026.96 |   0.0000
#else
    long double LICO_C_coef = 0.0000499143L;
    // long double LICO_C_coef = 0.000046L;
    // long double LICO_C_coef = 0.8958612L;
    // long double LICO_C_coef = 0.237791L;
#endif

    // long double LICO_C_coef = 0.0019967L;

    // Inputs / state
    size_t total_size = 0;  // number of keys
    std::vector<gap_prefix_type> gap_prefix_sum;     // size n, prefix on gaps
    std::vector<gap_prefix_type> gap_prefix_squares; // size n, prefix on gap^2

    // Outputs
    std::vector<Page> pages;   // paging over keys, half-open
    std::vector<Block> blocks; // result blocks in key indices [s,t)
    std::string partition_type = "optimal";
    long double optimal_cost = 0.0L;
    long double greedy_cost  = 0.0L;

    lico_partition() = default;

    explicit lico_partition(const std::vector<K>& data,  size_t global_m = 4096, double lico_c =0.109967, size_t t_page_size = 1024) : global_m(global_m) {
        page_size = t_page_size;
        LICO_C_coef = lico_c;
        if (output_log) {
                std::cerr << global_m << " " << page_size << " " << LICO_C_coef << " " << LICO_S_scale << "\n";
                output_log = false;
        }
        total_size = data.size();
        build_page(data);
        build_prefix(data);
    }

    ~lico_partition() {
        std::vector<Block>().swap(blocks);
    }

    void memory_clean() {
        std::vector<gap_prefix_type>().swap(gap_prefix_sum);
        std::vector<gap_prefix_type>().swap(gap_prefix_squares);
        std::vector<Page>().swap(pages);
    }

    void gap_prefix_clean() {
        std::vector<gap_prefix_type>().swap(gap_prefix_sum);
        std::vector<gap_prefix_type>().swap(gap_prefix_squares);
    }

    BlockStats stats(size_t s, size_t t) const {
        BlockStats r{};
        if (t <= s) return r;

        const size_t n_keys = t - s;
        r.n_keys = n_keys;

        if (n_keys < 2) return r;

        // gaps are between consecutive keys: indices [s+1 .. t-1]
        const size_t g0 = s;     // prefix index to subtract
        const size_t g1 = t - 1; // last gap index to include in prefix

        const long double sum =
            (long double)(gap_prefix_sum[g1] - gap_prefix_sum[g0]);
        const long double sum_sq =
            (long double)(gap_prefix_squares[g1] - gap_prefix_squares[g0]);

        const size_t n_gaps = n_keys - 1;
        r.n_gaps = n_gaps;

        const long double mu = sum / (long double)n_gaps;
        long double var = 0.0L;
        if (n_gaps >= 2) {
            var = (sum_sq - (long double)n_gaps * mu * mu) / (long double)(n_gaps - 1);
            if (var < 0.0L) var = 0.0L; // numeric guard
        }

        r.mu = mu;
        r.var = var;
        return r;
    }

    static inline int highest_bit_position(size_t x) {
        return (int)(sizeof(size_t) * CHAR_BIT - 1) - __builtin_clzl(x);
    }

    static inline size_t fill_all_bits_from_msb(size_t x) {
        if (x == 0) return 0;
        int k = highest_bit_position(x);
        return (1UL << (k + 1)) - 1;
    }

    size_t epsilon_star_from(size_t n_keys, long double sigma2) {
        if (n_keys < 2 || !(sigma2 > 0.0L)) return 1;
        const long double n = (long double)(n_keys - 1); // #gaps
        const long double factor =
            (2.0L * std::log(2.0L)) * LICO_C_coef * LICO_S_scale * sigma2 * (n / (n + 1.0L));
        if (!(factor > 0.0L)) return 1;

        long double eps = std::sqrt(factor);
        if (!std::isfinite(eps) || eps <= 1.0L)
            return 1;
        const long double eps_up = std::ceil(eps);
        if (eps_up > (long double)std::numeric_limits<size_t>::max())
            return std::numeric_limits<size_t>::max();
        size_t eps_size_t = (size_t) eps_up;

        // return eps_size_t; // no scale


        if (eps_size_t > 1)
            return fill_all_bits_from_msb(eps_size_t);
        else
            return 1; // epsilon* must be at least 1
    }

    size_t epsilon_star_range(size_t s, size_t t) {
        const auto st = stats(s, t);
        return epsilon_star_from(st.n_keys, st.var);
    }

    // S(P|ε*) ≈ n_gaps * ( 1.957 + 0.5 * log2(C*S*sigma^2) )
    long double bits_cost_star_from(size_t n_keys, long double sigma2) {
        if (n_keys < 2 || !(sigma2 > 0.0L)) return 0.0L;
        const size_t n_gaps = n_keys - 1;
        long double term = 1.957L + 0.5L * safe_log2l(LICO_C_coef * LICO_S_scale * sigma2);
        if (term < 0.0L) term = 0.0L; // never negative bits
        return (long double)n_gaps * term;
    }

    long double bits_cost_range(size_t s, size_t t) {
        const auto st = stats(s, t);
        return bits_cost_star_from(st.n_keys, st.var);
    }

    // Merge adjacent blocks that are contiguous and have identical epsilon.
    // If recompute_eps is true, we will also recompute epsilon* for the merged range;
    // otherwise we keep the shared epsilon as-is (since they matched).
    void merge_adjacent_same_epsilon(bool recompute_eps = false) {
        // return; // no scale

        if (blocks.size() < 2) return;

        // Ensure blocks are in ascending order by start index
        std::sort(blocks.begin(), blocks.end(),
                [](const Block& a, const Block& b){ return a.start_idx < b.start_idx; });

        std::vector<Block> merged;
        merged.reserve(blocks.size());

        Block cur = blocks[0];
        for (size_t i = 1; i < blocks.size(); ++i) {
            const Block& nxt = blocks[i];

            const bool contiguous = (cur.end_idx == nxt.start_idx);
            const bool same_eps   = (cur.epsilon == nxt.epsilon);

            if (contiguous && same_eps) {
                // Extend current block to include next
                cur.end_idx = nxt.end_idx;
                // keep epsilon as-is unless we choose to recompute after the loop
                // cost will be recomputed after the loop on the final merged ranges
            } else {
                // finalize current and start new
                merged.push_back(cur);
                cur = nxt;
            }
        }
        merged.push_back(cur);

        // Recompute cost (and optionally epsilon) for each merged block
        long double total_bits = 0.0L;
        for (auto &blk : merged) {
            if (recompute_eps) {
                blk.epsilon = epsilon_star_range(blk.start_idx, blk.end_idx);
            }
            blk.cost = bits_cost_range(blk.start_idx, blk.end_idx);
            total_bits += blk.cost;
        }

        blocks.swap(merged);

        // Update partition totals
        if (partition_type == "optimal") {
            optimal_cost = total_bits;
        } else if (partition_type == "greedy") {
            greedy_cost = total_bits;
        }
    }


    // ---------- Optimal (DP) over pages ----------
    void optimal_partition() {
        partition_type = "optimal";
        optimal_cost = 0.0L;
        blocks.clear();

        const size_t n = pages.size();
        if (n == 0) return;

        // Feasibility: each block must have at least MIN_PAGES_PER_BLOCK pages.
        const size_t max_blocks_by_pages =
            (MIN_PAGES_PER_BLOCK == 0) ? n : (n / MIN_PAGES_PER_BLOCK);

        size_t m = global_m;
        if (m == 0) m = 1;
        if (m > max_blocks_by_pages) {
            // std::cerr << "[optimal_partition] requested m=" << m
            //         << " > feasible " << max_blocks_by_pages
            //         << " under MIN_PAGES_PER_BLOCK=" << MIN_PAGES_PER_BLOCK
            //         << ". Using m=" << max_blocks_by_pages << " instead.\n";
            m = max_blocks_by_pages;
        }

        const long double INF = 1e300L;

        // dp[j][k] = min cost for first j pages with exactly k blocks
        std::vector<std::vector<long double>> dp(n + 1, std::vector<long double>(m + 1, INF));
        std::vector<std::vector<int>> prev_i(n + 1, std::vector<int>(m + 1, -1));
        std::vector<std::vector<size_t>> prev_eps(n + 1, std::vector<size_t>(m + 1, 0));

        dp[0][0] = 0.0L;

        // j: end page (exclusive), k: block count, i: previous cut
        for (size_t j = 1; j <= n; ++j) {
            // With exactly k blocks, we must be able to place k*MIN_PAGES_PER_BLOCK pages.
            const size_t k_max_here = std::min(m, j / std::max<size_t>(1, MIN_PAGES_PER_BLOCK));
            for (size_t k = 1; k <= k_max_here; ++k) {
                long double best = INF;
                int best_i = -1;
                size_t best_eps = 0;

                // i must allow k-1 blocks before i, and at least MIN_PAGES_PER_BLOCK for [i,j)
                const size_t i_min = (k - 1) * MIN_PAGES_PER_BLOCK;
                const size_t i_max = (j >= MIN_PAGES_PER_BLOCK) ? (j - MIN_PAGES_PER_BLOCK) : 0;
                if (i_min > i_max) continue;

                for (size_t i = i_min; i <= i_max; ++i) {
                    if (!(dp[i][k - 1] < INF)) continue;

                    const size_t s_key = pages[i].start_idx;
                    const size_t t_key = pages[j - 1].end_idx; // exclusive
                    const long double c = bits_cost_range(s_key, t_key);
                    const long double cand = dp[i][k - 1] + c;

                    if (cand < best) {
                        best = cand;
                        best_i = static_cast<int>(i);
                        best_eps = epsilon_star_range(s_key, t_key);
                    }
                }

                dp[j][k] = best;
                prev_i[j][k] = best_i;
                prev_eps[j][k] = best_eps;
            }
        }

        // EXACTLY m blocks: take dp[n][m]
        if (!(dp[n][m] < INF)) {
            // unreachable under constraints; emit one big block as a hard fallback
            std::cerr << "[optimal_partition] dp[n][m] unreachable; emitting single block fallback.\n";
            Block blk;
            blk.start_idx = pages[0].start_idx;
            blk.end_idx   = pages[n - 1].end_idx;
            blk.epsilon   = epsilon_star_range(blk.start_idx, blk.end_idx);
            blk.cost      = bits_cost_range(blk.start_idx, blk.end_idx);
            optimal_cost  = blk.cost;
            blocks.push_back(blk);
            return;
        }

        // Backtrack exactly m blocks
        std::vector<Block> rev;
        rev.reserve(m);

        int j = static_cast<int>(n);
        int k = static_cast<int>(m);
        while (k > 0 && j > 0) {
            int i = prev_i[j][k];
            if (i < 0) {
                // safety (shouldn't happen if dp[n][m] < INF)
                Block blk;
                blk.start_idx = pages[0].start_idx;
                blk.end_idx   = pages[j - 1].end_idx;
                blk.epsilon   = epsilon_star_range(blk.start_idx, blk.end_idx);
                blk.cost      = bits_cost_range(blk.start_idx, blk.end_idx);
                rev.push_back(blk);
                break;
            }

            Block blk;
            blk.start_idx = pages[(size_t)i].start_idx;
            blk.end_idx   = pages[(size_t)j - 1].end_idx; // exclusive
            blk.epsilon   = prev_eps[j][k];               // ε* chosen for [i,j)
            blk.cost      = bits_cost_range(blk.start_idx, blk.end_idx);
            rev.push_back(blk);

            j = i;
            --k;
        }

        std::reverse(rev.begin(), rev.end());
        blocks.swap(rev);
        optimal_cost = dp[n][m];

        // Optional: merge neighbors with identical epsilon
        merge_adjacent_same_epsilon(/*recompute_eps=*/false);
    }

    // ---------- Greedy (top-down) over pages ----------
    // Splits the most expensive interval until we have up to global_m blocks
    void greedy_partition() {
        partition_type = "greedy";
        greedy_cost = 0.0L;
        blocks.clear();

        using Interval = std::tuple<long double, size_t, size_t>; // (cost, s_page, t_page) with [s,t)
        auto cmp = [](const Interval& a, const Interval& b){ return std::get<0>(a) < std::get<0>(b); }; // max-heap
        std::priority_queue<Interval, std::vector<Interval>, decltype(cmp)> Q(cmp);

        const size_t n = pages.size();
        if (n == 0) return;

        // initial whole range
        {
            const size_t s_key = pages[0].start_idx;
            const size_t t_key = pages[n - 1].end_idx; // exclusive
            const long double full = bits_cost_range(s_key, t_key);
            Q.emplace(full, 0, n);
        }

        // helper to finalize a page interval [s,t)
        auto finalize_block = [&](size_t s, size_t t) {
            Block blk;
            blk.start_idx = pages[s].start_idx;
            blk.end_idx   = pages[t - 1].end_idx; // exclusive
            blk.epsilon   = epsilon_star_range(blk.start_idx, blk.end_idx);
            blk.cost      = bits_cost_range(blk.start_idx, blk.end_idx);
            greedy_cost  += blk.cost;
            blocks.push_back(blk);
        };

        if (LICO_C_coef != 0) {
            while (Q.size() < global_m && !Q.empty()) {
                auto [cur_cost, s, t] = Q.top();
                Q.pop();

                // number of pages in interval
                const size_t pn = t - s;
                // To split into [s,k) and [k,t), both must have at least MIN_PAGES_PER_BLOCK pages
                const size_t k_begin = s + MIN_PAGES_PER_BLOCK;
                const size_t k_end   = t - MIN_PAGES_PER_BLOCK;

                if (pn < 2 * MIN_PAGES_PER_BLOCK || k_begin > k_end) {
                    finalize_block(s, t);
                    continue;
                }

                // Find best split k
                size_t best_k = k_begin;
                long double best_sum = 1e300L;
                for (size_t k = k_begin; k <= k_end; ++k) {
                    const long double left  = bits_cost_range(pages[s].start_idx,     pages[k - 1].end_idx);
                    const long double right = bits_cost_range(pages[k].start_idx,     pages[t - 1].end_idx);
                    const long double sum_lr = left + right;
                    if (sum_lr < best_sum) { best_sum = sum_lr; best_k = k; }
                }

                // Push children intervals back to the heap
                {
                    const long double left_cost =
                        bits_cost_range(pages[s].start_idx, pages[best_k - 1].end_idx);
                    Q.emplace(left_cost, s, best_k);
                }
                {
                    const long double right_cost =
                        bits_cost_range(pages[best_k].start_idx, pages[t - 1].end_idx);
                    Q.emplace(right_cost, best_k, t);
                }
            }
        } else {
            LICO_C_coef = 0.001L;
        }


        // Finalize remaining intervals
        while (!Q.empty()) {
            auto [_, s, t] = Q.top();
            Q.pop();
            finalize_block(s, t);
        }

        // Sort blocks by key start
        std::sort(blocks.begin(), blocks.end(),
                  [](const Block& a, const Block& b){ return a.start_idx < b.start_idx; });
        // Merge neighbors sharing the same epsilon; keep epsilon as-is
        merge_adjacent_same_epsilon(/*recompute_eps=*/false);
    }

    // ---------- Reporting ----------
    void summarize(std::string outputbase_name){
        std::ofstream file(outputbase_name + ".internal_log.txt");
        file << "Partition type: " << partition_type << "\n";
        file << "Total Blocks: " << blocks.size() << "\n";
        long double total_bits = 0.0L;
        for (const auto &blk : blocks) {
            total_bits += blk.cost;
            file << "Block: [" << blk.start_idx << ", " << blk.end_idx
                      << "), eps*: " << blk.epsilon
                      << ", bits: " << (double)blk.cost << "\n";
        }
        file << "Total Bits: " << (double)total_bits << "\n";
    }

    void calculate_max_min_page_gap_variance(std::ofstream &file_gap_sta){
        long double max_gap_var = 0.0L;
        long double min_gap_var = 1e300L;
        for (const auto &page: pages) {
            auto page_stats = stats(page.start_idx, page.end_idx);
            if (page_stats.n_gaps == 0) continue;
            if (page_stats.var > max_gap_var) max_gap_var = page_stats.var;
            if (page_stats.var < min_gap_var) min_gap_var = page_stats.var;
        }
        file_gap_sta << "Max Page Gap Variance:\t" << (double)max_gap_var
                     << "\tMin Page Gap Variance:\t" << (double)min_gap_var
                     << "\tPage Number:\t" << pages.size() << "\t";
    }

    void summarize_gap_variance_per_page(std::ofstream &file_gap_sta) {
        file_gap_sta << pages.size() << "\n";
        for (const auto &page: pages) {
            auto page_stats = stats(page.start_idx, page.end_idx);
            // if (page_stats.n_gaps == 0) continue;
            file_gap_sta << page_stats.var << "\n";
        }
    }

private:
    // Build prefix sums of gaps (strictly increasing keys required)
    void build_prefix(const std::vector<K>& data) {
        const size_t n = data.size();
        gap_prefix_sum.assign(n, 0);
        gap_prefix_squares.assign(n, 0);
        if (n == 0) return;

        for (size_t i = 1; i < n; ++i) {
            if (data[i] <= data[i - 1])
                throw std::invalid_argument("Data must be strictly increasing (positive gaps).");
            const uint64_t g = data[i] - data[i - 1];
            gap_prefix_sum[i]      = gap_prefix_sum[i - 1] + (gap_prefix_type)g;
            gap_prefix_squares[i]  = gap_prefix_squares[i - 1]
                                   + (gap_prefix_type)g * (gap_prefix_type)g;
        }
    }

    // Build half-open pages over key indices
    void build_page(const std::vector<K>& data) {
        const size_t n = data.size();
        pages.clear();
        size_t s = 0;
        while (s < n) {
            size_t e = s + page_size;
            if (e > n) e = n; // exclusive
            pages.push_back({s, e});
            s = e;
        }
    }
};

