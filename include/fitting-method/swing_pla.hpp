#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#ifdef _OPENMP
// #include <omp.h>
#else
// #pragma message("Compilation with -fopenmp is optional but recommended")
#define omp_get_num_procs() 1
#define omp_get_max_threads() 1
#endif



namespace lico::internal_swing{
    template<typename T>
    using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>, long double,
                                                    std::conditional_t<(sizeof(T) < 8), int64_t, int64_t>>;

    template <typename X, typename Y>
    class SwingPiecewiseLinearModel {
    private:
        using SX = LargeSigned<X>;
        using SY = LargeSigned<Y>;

        struct Point {
            X x{};
            Y y{};

            // Get the slope between self and p
            long double operator-(const Point& p) const {
                if (p.x == x) {
                    throw std::invalid_argument("can not handle two same points");
                }
                long double dx = static_cast<long double>(x) - static_cast<long double>(p.x);
                long double dy = static_cast<long double>(y) - static_cast<long double>(p.y);

                long double slope = dy / dx;
                return slope;
            };
        };

        X first_x = 0;
        X last_x = 0;
        Point end[2];
        int points = 0;

        const Y epsilon;
        long double rho_max, rho_min;
        Point so;

    public:
        class CanonicalSegment;

        explicit SwingPiecewiseLinearModel(Y epsilon) : epsilon(epsilon) {
            if (epsilon < 0)
                throw std::invalid_argument("epsilon cannot be negative");
        }

        bool add_point(const X& x, const Y& y) {
            auto max_y = std::numeric_limits<Y>::max();
            auto min_y = std::numeric_limits<Y>::lowest();
            Point p1{x, y >= max_y - epsilon ? max_y : y + epsilon};
            Point p2{x, y <= min_y + epsilon ? min_y : y - epsilon};

            if (points == 0) {
                // It is the first point of the segment
                first_x = x;
                so = {x, y};
                points++;
                last_x = x;
                return true;
            }

            if (points == 1) {
                end[0] = p1;
                end[1] = p2;
                points++;
                rho_max = so - end[0];
                rho_min = so - end[1];
                last_x = x;
                return true;
            }

            long double slope1 = so - p2; // lower bound
            long double slope2 = so - p1; // upper bound
            bool outside_line1 = slope1 > rho_max;
            bool outside_line2 = slope2 < rho_min;

            if (outside_line1 || outside_line2) {
                points = 0;
                return false;
            }

            if (slope2 < rho_max) {
                end[0] = p1;
                rho_max = slope2;
            }

            if (slope1 > rho_min) {
                end[1] = p2;
                rho_min = slope1;
            }
            last_x = x;
            return true;
        }

        CanonicalSegment get_segment() {
            if (points == 1) {
                return CanonicalSegment(so, first_x, last_x);
            }
            return CanonicalSegment(so, end[0], end[1], first_x, last_x);
        }

        void reset() {
            points = 0;
        }
    };

    template <typename X, typename Y>
    class SwingPiecewiseLinearModel<X, Y>::CanonicalSegment {
    public:
        CanonicalSegment() = default;

        friend class SwingPiecewiseLinearModel;

        Point so;
        Point end[2];
        X first_x;
        X last_x;

        struct Segment {
            Segment(const long double& slope, const long double& intercept, const X& first_x, const X& last_x) :
                slope(slope), intercept(intercept), first_x(first_x), last_x(last_x) {};
            Segment() = default;

        public:
            long double slope;
            long double intercept;
            X first_x;
            X last_x;
        };

        CanonicalSegment(Point so, X first_x, X last_x) : so(so), first_x(first_x), last_x(last_x) {
            end[0] = so;
            end[1] = so;
        };

        CanonicalSegment(Point so, const Point& p1, const Point& p2, X first_x, X last_x)
            : so(so), end{p1, p2}, first_x(first_x), last_x(last_x) {};

        bool one_point() const {
            return so.x == end[0].x;
        }

        // float slope
        Segment get_canonical_segment(const X& origin) const {
            Point p1 = this->end[0];
            Point p2 = this->end[1];


            if (one_point())
                return {0, static_cast<long double>(end[0].y + end[1].y) / 2, first_x, last_x};

            auto [i_x, i_y] = so;


            long double rho_max = (static_cast<long double>(p1.y) - i_y) / (static_cast<long double>(p1.x) - i_x);

            long double rho_min = (static_cast<long double>(p2.y) - i_y) / (static_cast<long double>(p2.x) - i_x);

            long double slope = (rho_max + rho_min) / 2.;
            auto intercept = i_y - (i_x - origin) * slope;
            Segment result = Segment(slope, intercept, first_x, last_x);
            return result;
        }

        X get_first_x() const { return first_x; }


        // Here we choose the "max slope" bound through the pivot 'so' and the upper envelope end[0].
        std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> get_lico_segment() const {
            if (one_point()) {
                return {0, 0, 0, 0};
            }

            // Use upper bound (y+epsilon) at last_x as the max slope anchor
            const auto& p_up = end[0];

            // Compute deltas relative to pivot 'so'
            // Assumptions: X, Y are non-decreasing domain/codomain (e.g., x sorted, y is rank/index),
            // so delta_x >= 1 and delta_y >= 0 for multi-point segments.
            auto dx_signed = static_cast<LargeSigned<X>>(p_up.x) - static_cast<LargeSigned<X>>(so.x);
            auto dy_signed = static_cast<LargeSigned<Y>>(p_up.y) - static_cast<LargeSigned<Y>>(so.y);

            uint32_t delta_x = static_cast<uint32_t>(dx_signed);
            uint32_t delta_y = static_cast<uint32_t>(dy_signed);
            uint32_t x_b = static_cast<uint32_t>(so.x);
            uint32_t y_b = static_cast<uint32_t>(so.y);

            return {delta_y, delta_x, y_b, x_b};
        }
    };

    template <typename Fin, typename Fout>
    size_t make_segmentation(size_t n, double epsilon, Fin in, Fout out) {
        size_t c = 0;
        auto p = in(0);

        using X = typename std::invoke_result_t<Fin, size_t>::first_type;
        using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
        SwingPiecewiseLinearModel<X, Y> opt(epsilon);
        auto add_point = [&](X x, Y y) {
            if (!opt.add_point(x, y)) {
                out(opt.get_segment());
                opt.add_point(x, y);
                ++c;
            }
        };

        opt.add_point(p.first, p.second);

        for (size_t i = 1; i < n; ++i) {
            auto next_p = in(i);
            if (next_p.second == p.second) {
                if (next_p.second != in(i+1).second){
                    add_point(next_p.first, next_p.second);
                }
                continue;
            }
            p = next_p;
            add_point(p.first, p.second);
        }

        out(opt.get_segment());
        return ++c;
    }

    template <typename Fin, typename Fout>
    size_t make_segmentation_par(size_t n, size_t epsilon, Fin in, Fout out) {
        // auto parallelism = std::min<size_t>(omp_get_max_threads(), 20); // 20 is an arbitrary limit
        // auto chunk_size = n / parallelism;
        // auto c = 0ull;

        // if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation(n, epsilon, in, out);

//         using K = typename std::invoke_result_t<Fin, size_t>;
//         using canonical_segment = typename SwingPiecewiseLinearModel<K, int>::CanonicalSegment;
//         std::vector<std::vector<canonical_segment>> results(parallelism);
//
// #pragma omp parallel for reduction(+ : c) num_threads(parallelism)
//         for (auto i = 0; i < parallelism; ++i) {
//             auto first = i * chunk_size;
//             auto last = i == parallelism - 1 ? n : first + chunk_size;
//             if (first > 0) {
//                 for (; first < last; ++first)
//                     if (in(first) != in(first - 1))
//                         break;
//                 if (first == last)
//                     continue;
//             }
//
//             auto in_fun = [in](auto j) { return in(j); };
//             auto out_fun = [&results, i](const auto& cs) { results[i].emplace_back(cs); };
//             results[i].reserve(chunk_size / (epsilon > 0 ? epsilon * epsilon : 16));
//             c += make_segmentation(n, first, last, epsilon, in_fun, out_fun);
//         }
//
//         for (auto& v : results)
//             for (auto& cs : v)
//                 out(cs);
//
//         return c;
    }
}
