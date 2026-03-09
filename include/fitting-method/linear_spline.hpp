#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>

using K=uint32_t;

namespace ls {

    struct Segment {
        int32_t startIndex;
        double slope;           // stored parameter only
        Segment(int32_t startIndex, double slope): startIndex(startIndex), slope(slope) {}
    };

    struct FitResult {
        std::vector<Segment> segments;
        std::size_t segmentBytes;
        std::size_t errorBytes;
        int32_t epsilon;
        std::size_t residualViolations; // count of |residual| > epsilon after fitting
    };

    struct EncodedSpline {
        std::vector<Segment> segments;          // segment start indices and slopes
        K firstValue;                // seed value at index 0
        std::vector<std::int32_t> residuals;    // residuals for all non-start indices in order
        std::size_t n;                          // original array length
        int32_t epsilon;                   // max absolute residual bound
    };

    static inline std::size_t ceilDiv(std::size_t x, std::size_t y) {
        return x == 0 ? 0 : (x + y - 1) / y;
    }

    inline FitResult fitLinearSpline_swing(const std::vector<K>& data, int32_t epsilon, int32_t offset = 0) {
        FitResult result{};
        result.epsilon = epsilon;
        result.residualViolations = 0;

        const std::size_t n = data.size();
        if (n == 0) {
            result.segmentBytes = 0;
            result.errorBytes = 0;
            return result;
        }

        if (n == 1) {
            Segment s{offset, 0.0};
            result.segments.push_back(s);
            // Storage: slope (double) + startIndex (size_t). First startIndex could be implicit (0) if desired.
            result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::size_t));
            // Error storage: N * bits per value
            std::size_t bits = (epsilon <= 0) ? 0 : static_cast<std::size_t>(std::ceil(std::log2(2.0 * static_cast<double>(epsilon) + 1.0)));
            result.errorBytes = ceilDiv(n * bits, static_cast<std::size_t>(8));
            return result;
        }

        const double INF = std::numeric_limits<double>::infinity();

        int32_t s = 0; // current segment start index
        double minSlope = -INF;
        double maxSlope = INF;

        for (int32_t i = s + 1; i < n; ++i) {
            const int32_t dx = i - s;
            const double dy = static_cast<double>(data[i]) - static_cast<double>(data[s]);
            // Use a tightened bound so that integer rounding keeps residual within ±epsilon.
            // Since we anchor b at the segment start value exactly, the allowable deviation
            // at distance dx is ±epsLine.
            const double epsLine = std::max(0.0, static_cast<double>(epsilon));
            const double upper = (dy + epsLine) / static_cast<double>(dx);
            const double lower = (dy - epsLine) / static_cast<double>(dx);

            const double newMin = std::max(minSlope, lower);
            const double newMax = std::min(maxSlope, upper);

            if (newMin <= newMax) {
                minSlope = newMin;
                maxSlope = newMax;
                continue;
            }

            // Including point i would violate feasibility. Close segment at e = i-1
            const int32_t e = i - 1;
            const double a = (std::isfinite(minSlope) && std::isfinite(maxSlope)) ? (0.5 * (minSlope + maxSlope)) : 0.0;
            result.segments.push_back(Segment{s + offset, a});

            // Start new segment at e (overlap boundary)
            s = e;
            // Initialize bounds using current point i relative to new s
            const int32_t dx2 = i - s;
            const double dy2 = static_cast<double>(data[i]) - static_cast<double>(data[s]);
            const double upper2 = (dy2 + epsLine) / static_cast<double>(dx2);
            const double lower2 = (dy2 - epsLine) / static_cast<double>(dx2);
            minSlope = lower2;
            maxSlope = upper2;
        }

        // Finalize last segment from s to n-1
        {
            int32_t e = n - 1;
            const int32_t dx_e = e - s;
            double a = 0.0;
            if (dx_e > 0) {
                const double dy_e = static_cast<double>(data[e]) - static_cast<double>(data[s]);
                const double epsLine = std::max(0.0, static_cast<double>(epsilon) - 0.5);
                const double upper_e = (dy_e + epsLine) / static_cast<double>(dx_e);
                const double lower_e = (dy_e - epsLine) / static_cast<double>(dx_e);
                const double hi = std::min(maxSlope, upper_e);
                const double lo = std::max(minSlope, lower_e);
                a = (std::isfinite(lo) && std::isfinite(hi)) ? (0.5 * (lo + hi)) : 0.0;
            }
            result.segments.push_back(Segment{s + offset, a});
        }

        // Storage computations
        // Minimal storage per segment: slope (double) + startIndex (size_t). End index is inferred.
        result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(int32_t));

        // Error storage: bits needed to represent residual in [-epsilon, +epsilon]
        // If epsilon <= 0, assume no residual storage needed.
        // Integer rounding can add up to 0.5, so allow +/- (epsilon+1).
        // Residuals are guaranteed within ±epsilon after rounding by using tightened fit.
        std::size_t bitsPerValue = (epsilon <= 0) ? 0 : static_cast<std::size_t>(
            std::ceil(std::log2(2.0 * static_cast<double>(epsilon) + 1.0))
        );
        // Error storage: N * bits per value, as requested
        result.errorBytes = ceilDiv(n * bitsPerValue, 8);

        // Count residual violations using the overlap boundary convention used by encode/decode.
        // For each segment, take start s and end e = next.startIndex (or n-1 for last),
        // compute integer prediction with llround and residual r = data[i] - pred.
        // Count where |r| > epsilon.
        // std::size_t violations = 0;
        // const std::size_t segCount = result.segments.size();
        // for (std::size_t si = 0; si < segCount; ++si) {
        //     const auto& seg = result.segments[si];
        //     std::size_t s = seg.startIndex;
        //     std::size_t e = (si + 1 < segCount) ? result.segments[si + 1].startIndex : (n - 1);
        //     const std::int64_t y_s = data[s];
        //     const double b = static_cast<double>(y_s) - seg.slope * static_cast<double>(s);
        //     for (std::size_t i = s + 1; i <= e; ++i) {
        //         const double pred = seg.slope * static_cast<double>(i) + b;
        //         const std::int64_t p = static_cast<std::int64_t>(std::llround(pred));
        //         const std::int64_t r = data[i] - p;
        //         if (std::llabs(r) > epsilon) ++violations;
        //     }
        // }
        // result.residualViolations = violations;

        return result;
    }

    inline FitResult fitLinearSpline_rw(const std::vector<K>& data, int32_t epsilon, int32_t offset = 0) {
        FitResult result{};
        result.epsilon = epsilon;
        result.residualViolations = 0;

        const std::size_t n = data.size();
        if (n == 1) {
            result.segments.push_back(Segment{0 + offset, 0.0});
            // 每段最小存储：slope(double) + startIndex(int32_t)
            result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));
            // 残差位：覆盖 [-ε, +ε]，位数 ceil(log2(2ε+1))
            const int32_t eps = std::max<int32_t>(0, epsilon);
            const std::size_t bits = (eps <= 0) ? 0 : static_cast<std::size_t>(
                std::ceil(std::log2(2.0 * static_cast<double>(eps) + 1.0))
            );
            result.errorBytes = ceilDiv(n * bits, static_cast<std::size_t>(8));
            return result;
        }

        const int32_t eps = std::max<int32_t>(0, epsilon);

        std::size_t s = 0;
        while (s < n - 1) {
            // 用首两个点确定斜率（索引步长为 1）
            const std::int64_t ys  = static_cast<std::int64_t>(data[s]);
            const std::int64_t ys1 = static_cast<std::int64_t>(data[s + 1]);

            // 固定段斜率为首差：a = y_{s+1} - y_s
            const double a_dir = static_cast<double>(ys1 - ys);
            // 截距使直线通过 (s, y_s)：b = y_s - a*s
            const double b = static_cast<double>(ys) - a_dir * static_cast<double>(s);

            // 扫描后续点，逐点用该斜率预测并检查残差是否在 ±epsilon
            std::size_t lastGoodE = s + 1; // 至少包含 s+1（该点残差为 0）
            std::size_t j = s + 1;
            while (j < n) {
                const double pred = a_dir * static_cast<double>(j) + b;
                const std::int64_t p = static_cast<std::int64_t>(std::llround(pred));
                const std::int64_t r = static_cast<std::int64_t>(data[j]) - p;

                if (std::llabs(r) <= eps) {
                    lastGoodE = j;
                    ++j;
                } else {
                    break; // 首个超出 ±epsilon 的点
                }
            }

            // 输出一段，斜率固定为 a_dir
            result.segments.push_back(Segment{
                static_cast<std::int32_t>(s + offset),
                a_dir
            });

            // 边界规则：本段终点 = 下段起点（重叠边界）
            if (lastGoodE == n - 1) break;
            s = lastGoodE;
        }

        // 存储开销：每段 slope(double) + startIndex(int32_t)
        result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));

        // 残差位宽估算：覆盖 [-ε, +ε]
        const std::size_t bitsPerValue = (eps <= 0) ? 0 : static_cast<std::size_t>(std::ceil(std::log2(2.0 * static_cast<double>(eps) + 1.0)));
        result.errorBytes = ceilDiv(n * bitsPerValue, static_cast<std::size_t>(8));


        // std::size_t violations = 0;
        // const std::size_t segCount = result.segments.size();
        // for (std::size_t si = 0; si < segCount; ++si) {
        //     const auto& seg = result.segments[si];
        //     const std::size_t segStart = static_cast<std::size_t>(seg.startIndex);
        //     const std::size_t segEnd = (si + 1 < segCount)
        //         ? static_cast<std::size_t>(result.segments[si + 1].startIndex)
        //         : (n - 1);
        //
        //     const std::int64_t y_s = static_cast<std::int64_t>(data[segStart]);
        //     const double b = static_cast<double>(y_s) - seg.slope * static_cast<double>(segStart);
        //
        //     for (std::size_t i = segStart + 1; i <= segEnd; ++i) {
        //         const double pred = seg.slope * static_cast<double>(i) + b;
        //         const std::int64_t p = static_cast<std::int64_t>(std::llround(pred));
        //         const std::int64_t r = static_cast<std::int64_t>(data[i]) - p;
        //         if (std::llabs(r) > eps) ++violations;
        //     }
        // }
        // result.residualViolations = violations;

        return result;
    }


    inline FitResult fitLinearSpline_1(const std::vector<K>& data, int32_t epsilon, int32_t offset = 0) {
    FitResult result{};
    result.epsilon = epsilon;
    result.residualViolations = 0;

    const std::size_t n = data.size();
    if (n == 1) {
        result.segments.push_back(Segment{0 + offset, 0.0});
        result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));
        const int32_t eps = std::max<int32_t>(0, epsilon);
        const std::size_t bits = (eps <= 0) ? 0 : static_cast<std::size_t>(
            std::ceil(std::log2(2.0 * static_cast<double>(eps) + 1.0))
        );
        result.errorBytes = ceilDiv(n * bits, static_cast<std::size_t>(8));
        return result;
    }

    const double eps = static_cast<double>(std::max<int32_t>(0, epsilon));

    // ---- GreedySplineCorridor (absolute error corridor ±eps) ----
    // 输出的每个 Segment 存：段起点 index + 段 slope
    // 截距 b 在解码端由 (startIndex, y[startIndex], slope) 恢复：b = y_s - slope*s

    auto y = [&](std::size_t i) -> double {
        return static_cast<double>(static_cast<std::int64_t>(data[i]));
    };

    // 计算在基点 B 固定时，点 i 对 slope 的可行区间 [lo, hi]
    // 约束：| (yB + a*(i-B)) - y_i | <= eps
    // => y_i - eps <= yB + a*(i-B) <= y_i + eps
    // => (y_i - eps - yB)/(i-B) <= a <= (y_i + eps - yB)/(i-B)
    auto slopeIntervalFromPoint = [&](std::size_t B, std::size_t i, double& lo, double& hi) {
        const double dy = y(i) - y(B);
        const double dx = static_cast<double>(i - B);
        lo = (dy - eps) / dx;
        hi = (dy + eps) / dx;
        if (lo > hi) std::swap(lo, hi);
    };

    // 选取一个确定 slope（从当前可行区间里）
    auto chooseSlope = [](double lo, double hi) -> double {
        // 最简单稳妥：取中点（保证在可行区间内）
        return 0.5 * (lo + hi);
    };

    std::size_t B = 0;               // base point index
    std::size_t lastFeasible = 1;    // 最后一个仍可由当前 B 连接的点（至少 1）

    // 初始可行 slope 区间由点 1 给出
    double lo, hi;
    slopeIntervalFromPoint(B, 1, lo, hi);

    // 先占位：段起点 B=0 的 slope 之后才能确定（当段结束时才知道最终 corridor）
    // 为了不改变输出结构，我们在段结束时 push Segment{B, chosenSlope}
    while (true) {
        // 扫描 C，从 i=2 开始（论文从 3 到 n，因其 1-based；这里 0-based）
        std::size_t i = lastFeasible + 1;
        for (; i < n; ++i) {
            double cli, chi;
            slopeIntervalFromPoint(B, i, cli, chi);

            // 更新 corridor：与当前 [lo,hi] 取交集
            const double newLo = std::max(lo, cli);
            const double newHi = std::min(hi, chi);

            if (newLo <= newHi) {
                // 仍可达：收紧 corridor，继续
                lo = newLo;
                hi = newHi;
                lastFeasible = i;
            } else {
                // 不可达：按照论文，选 i-1 作为新基点（knot）
                break;
            }
        }

        // 当前段确定：从 corridor 中选一个 slope
        const double a = chooseSlope(lo, hi);
        result.segments.push_back(Segment{
            static_cast<std::int32_t>(B + offset),
            a
        });

        if (lastFeasible >= n - 1) {
            // 已经能到末尾：结束（末点隐含由最后一段覆盖/或由外部使用 segments 推断）
            break;
        }

        // 开启下一段：新基点 = lastFeasible （即 i-1）
        B = lastFeasible;

        // 至少要从 B 连到 B+1 才能形成下一段；若数据长度保证 n>=2，这里安全
        const std::size_t next = B + 1;
        if (next >= n) break;

        lastFeasible = next;
        slopeIntervalFromPoint(B, next, lo, hi);
    }

    // 存储开销：每段 slope(double) + startIndex(int32_t)
    result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));

    // 残差位宽估算：覆盖 [-ε, +ε]
    const int32_t eps_i = std::max<int32_t>(0, epsilon);
    const std::size_t bitsPerValue = (eps_i <= 0)
        ? 0
        : static_cast<std::size_t>(std::ceil(std::log2(2.0 * static_cast<double>(eps_i) + 1.0)));
    result.errorBytes = ceilDiv(n * bitsPerValue, static_cast<std::size_t>(8));

    return result;
}

    inline FitResult fitLinearSpline_no_can(const std::vector<K>& data, int32_t epsilon, int32_t offset = 0) {
    FitResult result{};
    result.epsilon = epsilon;
    result.residualViolations = 0;

    const std::size_t n = data.size();
    if (n == 1) {
        result.segments.push_back(Segment{0 + offset, 0.0});
        result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));
        const int32_t eps = std::max<int32_t>(0, epsilon);
        const std::size_t bits = (eps <= 0) ? 0 : static_cast<std::size_t>(
            std::ceil(std::log2(2.0 * static_cast<double>(eps) + 1.0))
        );
        result.errorBytes = ceilDiv(n * bits, static_cast<std::size_t>(8));
        return result;
    }

    const double eps = static_cast<double>(std::max<int32_t>(0, epsilon));

    auto y = [&](std::size_t i) -> double {
        return static_cast<double>(static_cast<std::int64_t>(data[i]));
    };

    // 从 B 到点 P 的“上界射线”斜率：连接 (B, y(B)) 到 (P, y(P)+eps)
    auto slopeUpper = [&](std::size_t B, std::size_t P) -> double {
        const double dx = static_cast<double>(P - B);
        return ( (y(P) + eps) - y(B) ) / dx;
    };

    // 从 B 到点 P 的“下界射线”斜率：连接 (B, y(B)) 到 (P, y(P)-eps)
    auto slopeLower = [&](std::size_t B, std::size_t P) -> double {
        const double dx = static_cast<double>(P - B);
        return ( (y(P) - eps) - y(B) ) / dx;
    };

    // 从 B 到 C 的中心线斜率（不加 eps），用于可达性判定（对应 BC）
    auto slopeCenter = [&](std::size_t B, std::size_t C) -> double {
        const double dx = static_cast<double>(C - B);
        return ( y(C) - y(B) ) / dx;
    };

    // 选定要写入 Segment 的 slope：从当前 corridor 的两条边界（BU、BL）取中点
    auto chooseSlopeFromCorridor = [&](std::size_t B, std::size_t U, std::size_t L) -> double {
        const double su = slopeUpper(B, U);
        const double sl = slopeLower(B, L);
        return 0.5 * (su + sl);
    };

    // ---- GreedySplineCorridor (Figure 1 原文写法：B/U/L) ----
    std::size_t B = 0;
    // R 从 S[1] 开始：你的结构里段用 startIndex 表示 knot；先不立刻 push，
    // 等确定 corridor 后再为该段选择 slope 并 push Segment{B, slope}
    // 初始化 U = S[2] + eps, L = S[2] - eps  (0-based 即点 1)
    std::size_t U = 1;
    std::size_t L = 1;

    // i = 3..n (1-based) => i = 2..n-1 (0-based)
    for (std::size_t i = 2; i < n; ++i) {
        const std::size_t C = i;

        // 判定：if BC is left of BU or right of BL
        // 用斜率比较实现：不可达当 slope(BC) > slope(BU) 或 slope(BC) < slope(BL)
        const double sBU = slopeUpper(B, U);
        const double sBL = slopeLower(B, L);
        const double sBC = slopeCenter(B, C);

        if (sBC > sBU || sBC < sBL) {
            // 违反走廊：B = S[i-1], R append <B>
            // 在你结构里：先把“上一段（旧 B）”输出
            result.segments.push_back(Segment{
                static_cast<std::int32_t>(B + offset),
                chooseSlopeFromCorridor(B, U, L)
            });

            // 更新 base point 为 i-1
            B = i - 1;

            // U = C + eps, L = C - eps（对应原文：U=C+, L=C-）
            U = C;
            L = C;
        } else {
            // 仍可达：计算 U' / L'（候选 bound 来自当前 C）
            const std::size_t Uprime = C;
            const std::size_t Lprime = C;

            // 原文：U = C + eps, L = C - eps（先设为候选）
            // 然后：
            // if BU is left of BU' then U = U
            // if BL is right of BL' then L = L
            //
            // 用斜率比较表达“更紧”：
            // - 上界更紧：slopeUpper(B, Uprime) < slopeUpper(B, U)  => U = Uprime
            // - 下界更紧：slopeLower(B, Lprime) > slopeLower(B, L)  => L = Lprime
            const double sBUprime = slopeUpper(B, Uprime);
            const double sBLprime = slopeLower(B, Lprime);

            if (sBUprime < sBU) U = Uprime;
            if (sBLprime > sBL) L = Lprime;
        }
    }

    // 原文末尾：R append <S[n]>
    // 你的结构不显式存最后点；但需要输出最后一段（最后一个 base point B 开始的段）
    result.segments.push_back(Segment{
        static_cast<std::int32_t>(B + offset),
        chooseSlopeFromCorridor(B, U, L)
    });

    // 存储开销：每段 slope(double) + startIndex(int32_t)
    result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));

    // 残差位宽估算：覆盖 [-ε, +ε]
    const int32_t eps_i = std::max<int32_t>(0, epsilon);
    const std::size_t bitsPerValue = (eps_i <= 0)
        ? 0
        : static_cast<std::size_t>(std::ceil(std::log2(2.0 * static_cast<double>(eps_i) + 1.0)));
    result.errorBytes = ceilDiv(n * bitsPerValue, static_cast<std::size_t>(8));

    return result;
}

    inline FitResult fitLinearSpline_no(const std::vector<K>& data, int32_t epsilon, int32_t offset = 0) {
    FitResult result{};
    result.epsilon = epsilon;
    result.residualViolations = 0;

    const std::size_t n = data.size();
    if (n == 1) {
        result.segments.push_back(Segment{0 + offset, 0.0});
        result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));
        const int32_t eps_i = std::max<int32_t>(0, epsilon);
        const std::size_t bits = (eps_i <= 0) ? 0 : static_cast<std::size_t>(
            std::ceil(std::log2(2.0 * static_cast<double>(eps_i) + 1.0))
        );
        result.errorBytes = ceilDiv(n * bits, static_cast<std::size_t>(8));
        return result;
    }

    const double eps = static_cast<double>(std::max<int32_t>(0, epsilon));

    // 预留一点空间，减少扩容；真实段数通常远小于 n
    result.segments.reserve(std::min<std::size_t>(n, 256));

    auto yd = [&](std::size_t i) -> double {
        return static_cast<double>(static_cast<std::int64_t>(data[i]));
    };

    // ---- GreedySplineCorridor (Figure 1 原文写法：B/U/L) ----
    std::size_t B = 0;
    double yB = yd(B);

    // 初始化 U = S[2] + eps, L = S[2] - eps (0-based: index 1)
    std::size_t U = 1;
    std::size_t L = 1;

    // 维护当前走廊边界斜率（相对于当前 B）
    // sBU = slope of line from (B, yB) to (U, y(U)+eps)
    // sBL = slope of line from (B, yB) to (L, y(L)-eps)
    {
        const double y1 = yd(1);
        const double invDx = 1.0; // (1 - 0) == 1
        const double dy = y1 - yB;
        // U=L=1 所以两条都可直接算
        // (dy ± eps) / 1
        // 上界：y1+eps - yB = dy + eps
        // 下界：y1-eps - yB = dy - eps
        // 注意：这里的斜率是从 B 出发的边界射线
        double sBU = (dy + eps) * invDx;
        double sBL = (dy - eps) * invDx;

        // i = 2..n-1
        for (std::size_t i = 2; i < n; ++i) {
            const std::size_t C = i;
            const double yC = yd(C);
            const double invDxC = 1.0 / static_cast<double>(C - B);

            // sBC = (yC - yB)/(C-B)
            const double sBC = (yC - yB) * invDxC;

            // 不可达：BC is left of BU OR right of BL
            // 实现约定：越过上界 => sBC > sBU，越过下界 => sBC < sBL
            if (sBC > sBU || sBC < sBL) {
                // 输出上一段：从当前 corridor 选 slope
                // chooseSlope = 0.5*(sBU + sBL)
                result.segments.push_back(Segment{
                    static_cast<std::int32_t>(B + offset),
                    0.5 * (sBU + sBL)
                });

                // 更新 base point 为 i-1
                B = i - 1;
                yB = yd(B);

                // 重置 U/L 为当前 C，并重置走廊为只由 C 约束
                U = C;
                L = C;

                const double invDxReset = 1.0 / static_cast<double>(C - B); // 这里 (C-(i-1)) == 1
                const double dyReset = yC - yB;
                sBU = (dyReset + eps) * invDxReset;
                sBL = (dyReset - eps) * invDxReset;
            } else {
                // 候选 bound 来自当前 C（U' = C, L' = C）
                const double dyC = yC - yB;
                const double sBUprime = (dyC + eps) * invDxC;
                const double sBLprime = (dyC - eps) * invDxC;

                // 上界更紧：sBUprime < sBU
                if (sBUprime < sBU) {
                    sBU = sBUprime;
                    U = C;
                }
                // 下界更紧：sBLprime > sBL
                if (sBLprime > sBL) {
                    sBL = sBLprime;
                    L = C;
                }
            }
        }

        // 输出最后一段
        result.segments.push_back(Segment{
            static_cast<std::int32_t>(B + offset),
            0.5 * (sBU + sBL)
        });
    }

    // 存储开销：每段 slope(double) + startIndex(int32_t)
    result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));

    // 残差位宽估算：覆盖 [-ε, +ε]
    const int32_t eps_i = std::max<int32_t>(0, epsilon);
    const std::size_t bitsPerValue = (eps_i <= 0)
        ? 0
        : static_cast<std::size_t>(std::ceil(std::log2(2.0 * static_cast<double>(eps_i) + 1.0)));
    result.errorBytes = ceilDiv(n * bitsPerValue, static_cast<std::size_t>(8));

        std::size_t violations = 0;
        const std::size_t segCount = result.segments.size();
        for (std::size_t si = 0; si < segCount; ++si) {
            const auto& seg = result.segments[si];
            const std::size_t segStart = static_cast<std::size_t>(seg.startIndex);
            const std::size_t segEnd = (si + 1 < segCount)
                ? static_cast<std::size_t>(result.segments[si + 1].startIndex)
                : (n - 1);

            const std::int64_t y_s = static_cast<std::int64_t>(data[segStart]);
            const double b = static_cast<double>(y_s) - seg.slope * static_cast<double>(segStart);

            for (std::size_t i = segStart + 1; i <= segEnd; ++i) {
                const double pred = seg.slope * static_cast<double>(i) + b;
                const std::int64_t p = static_cast<std::int64_t>(std::llround(pred));
                const std::int64_t r = static_cast<std::int64_t>(data[i]) - p;
                if (std::llabs(r) > eps) ++violations;
            }
        }
        result.residualViolations = violations;

    return result;
}

    inline FitResult fitLinearSpline(const std::vector<K>& data, int32_t epsilon, int32_t offset = 0) {
        FitResult result{};
        result.epsilon = epsilon;
        result.residualViolations = 0;

        const std::size_t n = data.size();
        if (n == 1) {
            result.segments.push_back(Segment{0 + offset, 0.0});
            result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));
            const int32_t eps_i = std::max<int32_t>(0, epsilon);
            const std::size_t bits = (eps_i <= 0) ? 0 : static_cast<std::size_t>(
                std::ceil(std::log2(2.0 * static_cast<double>(eps_i) + 1.0))
            );
            result.errorBytes = ceilDiv(n * bits, static_cast<std::size_t>(8));
            return result;
        }

        // 外部误差仍是 epsilon；但为了保证 llround 后仍满足 ±epsilon，
        // corridor 用 epsCorr = epsilon - 0.5（最小为 0）
        const double epsOut  = static_cast<double>(std::max<int32_t>(0, epsilon));
        const double epsCorr = std::max(0.0, epsOut);

        result.segments.reserve(std::min<std::size_t>(n, 256));

        auto yd = [&](std::size_t i) -> double {
            return static_cast<double>(static_cast<std::int64_t>(data[i]));
        };

        std::size_t B = 0;
        double yB = yd(B);

        // 初始化 U/L = 1
        std::size_t U = 1, L = 1;

        // 初始 corridor 边界斜率（用 epsCorr）
        const double y1 = yd(1);
        double sBU = ( (y1 + epsCorr) - yB ); // /1
        double sBL = ( (y1 - epsCorr) - yB ); // /1

        for (std::size_t i = 2; i < n; ++i) {
            const std::size_t C = i;
            const double yC = yd(C);
            const double invDxC = 1.0 / static_cast<double>(C - B);

            const double sBC = (yC - yB) * invDxC;

            if (sBC > sBU || sBC < sBL) {
                // 输出上一段（取 corridor 中点）
                result.segments.push_back(Segment{
                    static_cast<std::int32_t>(B + offset),
                    0.5 * (sBU + sBL)
                });

                // 新 base point = i-1
                B = i - 1;
                yB = yd(B);

                // 重置 U/L = C，并重置 corridor（注意这里 dx=1）
                U = C; L = C;

                const double invDxReset = 1.0 / static_cast<double>(C - B); // == 1
                const double dyReset = yC - yB;
                sBU = (dyReset + epsCorr) * invDxReset;
                sBL = (dyReset - epsCorr) * invDxReset;
            } else {
                // 候选 bound 来自 C
                const double dyC = yC - yB;
                const double sBUprime = (dyC + epsCorr) * invDxC;
                const double sBLprime = (dyC - epsCorr) * invDxC;

                if (sBUprime < sBU) { sBU = sBUprime; U = C; }
                if (sBLprime > sBL) { sBL = sBLprime; L = C; }
            }
        }

        // 最后一段
        result.segments.push_back(Segment{
            static_cast<std::int32_t>(B + offset),
            0.5 * (sBU + sBL)
        });

        result.segmentBytes = result.segments.size() * (sizeof(double) + sizeof(std::int32_t));

        const int32_t eps_i = std::max<int32_t>(0, epsilon);
        const std::size_t bitsPerValue = (eps_i <= 0)
            ? 0
            : static_cast<std::size_t>(std::ceil(std::log2(2.0 * static_cast<double>(eps_i) + 1.0)));
        result.errorBytes = ceilDiv(n * bitsPerValue, static_cast<std::size_t>(8));

        std::size_t violations = 0;
        const std::size_t segCount = result.segments.size();
        for (std::size_t si = 0; si < segCount; ++si) {
            const auto& seg = result.segments[si];
            const std::size_t segStart = static_cast<std::size_t>(seg.startIndex-offset);
            const std::size_t segEnd = (si + 1 < segCount) ? static_cast<std::size_t>(result.segments[si + 1].startIndex-offset) : (n - 1);

            const std::int64_t y_s = static_cast<std::int64_t>(data[segStart]);
            const double b = static_cast<double>(y_s) - seg.slope * static_cast<double>(segStart);

            for (std::size_t i = segStart + 1; i <= segEnd; ++i) {
                const double pred = seg.slope * static_cast<double>(i) + b;
                const std::int64_t p = static_cast<std::int64_t>(std::llround(pred));
                const std::int64_t r = static_cast<std::int64_t>(data[i]) - p;
                if (std::llabs(r) > epsilon) ++violations;
            }
        }
        result.residualViolations = violations;

        return result;
    }

    inline EncodedSpline encodeLinearSpline_parted(const std::vector<Segment>& segments, const std::vector<K> &data) {
        EncodedSpline enc{};
        enc.n = data.size();

        enc.segments = std::move(segments);
        enc.segments.emplace_back(Segment{static_cast<int32_t> (enc.n - 1), 0});
        enc.firstValue = data[0];

        // Generate residuals for all non-start indices
        enc.residuals.clear();
        enc.residuals.reserve(enc.n ? (enc.n - enc.segments.size()) : 0);

        std::size_t segCount = enc.segments.size();
        for (int32_t si = 0; si < segCount - 1; ++si) {
            const auto& seg = enc.segments[si];
            const int32_t s = seg.startIndex;
            const int32_t e = enc.segments[si + 1].startIndex;
            // Overlap boundary: next segment starts at current end index
            // int32_t e = (si + 1 < segCount) ? enc.segments[si + 1].startIndex : (enc.n - 1);

            // Intercept derived from actual value at start.
            const K y_s = data[s];
            const double b = static_cast<double>(y_s) - seg.slope * static_cast<double>(s);

            // Residuals for indices strictly after s (including boundary e)
            for (std::size_t i = s + 1; i <= e; ++i) {
                const double pred = seg.slope * static_cast<double>(i) + b;
                const int32_t p = static_cast<std::int32_t>(std::llround(pred));
                const int32_t r = data[i] - p;
                enc.residuals.push_back(r);
            }
        }

        return enc;
    }


    inline EncodedSpline encodeLinearSpline(const std::vector<K>& data, int32_t epsilon) {
        EncodedSpline enc{};
        enc.n = data.size();
        enc.epsilon = epsilon;
        if (data.empty()) {
            enc.firstValue = 0;
            return enc;
        }

        // Fit segments
        FitResult fit = fitLinearSpline(data, epsilon);
        enc.segments = std::move(fit.segments);
        enc.segments.emplace_back(Segment{static_cast<int32_t> (enc.n - 1), 0});
        enc.firstValue = data[0];

        // Generate residuals for all non-start indices
        enc.residuals.clear();
        enc.residuals.reserve(enc.n ? (enc.n - enc.segments.size()) : 0);

        std::size_t segCount = enc.segments.size();
        for (int32_t si = 0; si < segCount - 1; ++si) {
            const auto& seg = enc.segments[si];
            const int32_t s = seg.startIndex;
            const int32_t e = enc.segments[si + 1].startIndex;
            // Overlap boundary: next segment starts at current end index
            // int32_t e = (si + 1 < segCount) ? enc.segments[si + 1].startIndex : (enc.n - 1);

            // Intercept derived from actual value at start.
            const K y_s = data[s];
            const double b = static_cast<double>(y_s) - seg.slope * static_cast<double>(s);

            // Residuals for indices strictly after s (including boundary e)
            for (std::size_t i = s + 1; i <= e; ++i) {
                const double pred = seg.slope * static_cast<double>(i) + b;
                const std::int64_t p = static_cast<std::int64_t>(std::llround(pred));
                const std::int64_t r = data[i] - p;
                enc.residuals.push_back(r);
            }
        }

        return enc;
    }

    inline void decodeLinearSpline(const EncodedSpline& enc, K* out) {
        out[0] = enc.firstValue;

        std::size_t rpos = 0; // residual cursor
        const std::size_t segCount = enc.segments.size();
        for (int32_t si = 0; si < segCount - 1; ++si) {
            const auto& seg = enc.segments[si];

            const int32_t e = enc.segments[si + 1].startIndex;
            const int32_t s = seg.startIndex;

            // Use previously decoded boundary value as start anchor
            const K y_s = out[s];
            const double b = static_cast<double>(y_s) - seg.slope * static_cast<double>(s);

            for (int32_t i = s + 1; i <= e; ++i) {
                const double pred = seg.slope * static_cast<double>(i) + b;
                const int32_t p = static_cast<int32_t>(std::llround(pred));
                out[i] = p + enc.residuals[rpos++];
            }
        }
    }
} // namespace ls
