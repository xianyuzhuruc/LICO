#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <algorithm>


#ifdef _OPENMP
// #include <omp.h>
#else
// #pragma message("Compilation with -fopenmp is optional but recommended")
#define omp_get_num_procs() 1
#define omp_get_max_threads() 1
#endif


namespace lico::internal_greedy {
    template<typename T>
    using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>, long double,
                                                    std::conditional_t<(sizeof(T) < 8), int64_t, int64_t>>;

    template <typename X, typename Y>
    class GreedyPiecewiseLinearModel {
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

        X first_x = 0; // First x value of the segment
        X last_x = 0; // Last x value of the segment
        Point initial[2]; // Initial points of the segment
        Point end[2]; // End points of the segment
        int points = 0; // Number of points in the segment, used to check if this is the first point or not

        const Y epsilon; // Error bound for the PLA method
        long double rho_max, rho_min; // Maximum and minimum slopes of the segment
        std::pair<long double, long double> so; // Pivot point of the segment, used to calculate the intercept


    public:
        class CanonicalSegment;

        explicit GreedyPiecewiseLinearModel(Y epsilon) : epsilon(epsilon) {
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
                initial[0] = p1;
                initial[1] = p2;
                points++;
                last_x = x;
                return true;
            }

            if (points == 1) {
                // It is the second point of the segment
                end[0] = p1;
                end[1] = p2;
                points++;
                long double so_x = (static_cast<long double>(end[0].x) + initial[0].x) / 2.0;
                long double so_y = (static_cast<long double>(end[0].y) + initial[0].y + end[1].y + initial[1].y) / 4.0;
                rho_min = initial[0] - end[1];
                rho_max = initial[1] - end[0];
                so = {so_x, so_y};
                last_x = x;
                return true;
            }

            double slope1 = (static_cast<long double>(y - epsilon) - so.second) / (static_cast<long double>(x) - so.
                first);
            double slope2 = (static_cast<long double>(y + epsilon) - so.second) / (static_cast<long double>(x) - so.
                first);
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
                return CanonicalSegment(initial[0], initial[1], first_x, last_x);
            }
            return CanonicalSegment(so, end[0], end[1], first_x, last_x);
        }

        void reset() {
            points = 0;
        }
    };

    template <typename X, typename Y>
    class GreedyPiecewiseLinearModel<X, Y>::CanonicalSegment {
    public:
        CanonicalSegment() = default;

        friend class GreedyPiecewiseLinearModel;

        std::pair<long double, long double> so;
        Point end[2];
        X first_x;
        X last_x;

        struct Segment {
            Segment(const long double& slope, const long double& intercept, const X& first_x, const X& last_x) :
                slope(slope), intercept(intercept), first_x(first_x), last_x(last_x) {};
            Segment() = default;

        public:
            long double slope;
            int64_t intercept;
            X first_x;
            X last_x;
        };

        CanonicalSegment(const Point& p1, const Point& p2, X first_x, X last_x) : end{p1, p2}, first_x(first_x),
            last_x(last_x) {
            so = {(p1.x + p2.x) / 2, (p1.y + p2.y) / 2};
        };

        CanonicalSegment(std::pair<long double, long double> so, const Point& p1, const Point& p2, X first_x, X last_x)
            : so(so), end{p1, p2}, first_x(first_x), last_x(last_x) {};

        bool one_point() const {
            return so.first == end[0].x;
        }

        Segment get_canonical_segment(const X& origin) const {
            Point p1 = this->end[0];
            Point p2 = this->end[1];


            if (one_point())
                return {0, static_cast<long double>(end[0].y + end[1].y) / 2, first_x, last_x};

            auto [i_x, i_y] = so;

            long double rho_max = (static_cast<long double>(p1.y) - i_y) / (static_cast<long double>(p1.x) - i_x);
            long double rho_min = (static_cast<long double>(p2.y) - i_y) / (static_cast<long double>(p2.x) - i_x);

            long double slope = (rho_max + rho_min) / 2.;
            long double raw_intercept = i_y - slope * (i_x - static_cast<double>(origin));
            int64_t intercept = static_cast<int64_t>(std::llround(raw_intercept));
            Segment result = Segment(slope, intercept, first_x, last_x);
            return result;
        }

        X get_first_x() const { return first_x; }

        std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> get_lico_segment() const {
            if (one_point()) {
                return {0u, 0u, 0u, 0u};
            }

            // Upper envelope point at the segment end (last_x)
            const auto& p_up = end[0];

            uint32_t x_b = p_up.first;
            uint32_t y_b = p_up.second;

            uint32_t delta_x = static_cast<uint32_t>(x_b - so.first);
            uint32_t delta_y = static_cast<uint32_t>(y_b - so.second);

            return {delta_y, delta_x, y_b, x_b};
        }

    };


    template <typename Fin, typename Fout>
    size_t make_segmentation(size_t n, double epsilon, Fin in, Fout out) {
        size_t c = 0;
        auto p = in(0);

        using X = typename std::invoke_result_t<Fin, size_t>::first_type;
        using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
        GreedyPiecewiseLinearModel<X, Y> opt(epsilon);
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
        // printf("The number of threads is %d\n", parallelism);
        // auto chunk_size = n / parallelism;
        // auto c = 0ull;

        // if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation(n, epsilon, in, out);

//         using K = typename std::invoke_result_t<Fin, size_t>;
//         using canonical_segment = typename GreedyPiecewiseLinearModel<K, int>::CanonicalSegment;
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

} // namespace Greedy::internal
