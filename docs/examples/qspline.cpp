// clang++ qspline.cpp -o qspline --std=c++2a -stdlib=libc++

#include <cstdio>
#include <limits>
#include <optional>
#include <vector>
#include <string>
#include <fstream>

using number_t = float;
const number_t NaN = std::numeric_limits<number_t>::quiet_NaN();

// p: (x, y)
struct Point2 {
    number_t x;
    number_t y;
};

// S: a(x-xs)^2 + b(x-xs) + c
struct SectionCurve2 {
    number_t a;
    number_t b;
    number_t c;
    number_t xs;
    number_t xe;
};

static std::vector<Point2> read_points_from_csv_file(char const * path);
static std::vector<SectionCurve2> qspline_for_points(std::vector<Point2> const & points, std::optional<number_t> opt_a0);
static number_t y_on_qspline(std::vector<SectionCurve2> const & curves, number_t x);

int main(int argc, char ** argv) {
    if (argc < 2) {
        std::puts("usage: qspline source [a0]");
        return 0;
    }
    std::optional<number_t> opt_a0;
    if (argc >= 3) {
        number_t a0;
        if (std::sscanf(argv[2], "%f", &a0) != 1) {
            std::printf("Error %d\n", -__LINE__);
            return -1;
        }
        opt_a0.emplace(a0);
    }
    try {
        // input
        auto const points = read_points_from_csv_file(argv[1]);
        auto const curves = qspline_for_points(points, opt_a0);
        // output
        auto range = points.back().x - points.front().x;
        for (auto x = points.front().x; x < points.back().x; x += range / 100) {
            std::printf("%f,%f\n", x, y_on_qspline(curves, x));
        }
        std::printf("%f,%f\n", points.back().x, y_on_qspline(curves, points.back().x));
    } catch (int & error) {
        std::printf("Error %d\n", error);
        return -1;
    }

    return 0;
}

static std::vector<Point2> read_points_from_csv_file(char const * path) {
    std::vector<Point2> points;

    std::ifstream ifs(path);
    if (ifs.fail()) {
        throw -__LINE__;
    }

    std::string line;
    while (std::getline(ifs, line)) {
        decltype(points)::value_type point;
        if (std::sscanf(line.c_str(), "%f,%f", &point.x, &point.y) != 2) {
            throw -__LINE__;
        }
        points.push_back(point);
    }
    return points;
}

static std::vector<SectionCurve2> qspline_for_points(std::vector<Point2> const & points, std::optional<number_t> opt_a0) {
    std::vector<SectionCurve2> S;

    if (!(points.size() >= 2)) {
        throw -__LINE__;
    }
    for (auto i = 0; i < points.size() - 1; ++i) {
        if (!(points[i].x < points[i+1].x)) {
            throw -__LINE__;
        }
    }

    auto const N = points.size() - 1; // N+1 points => N curves

    // Prepare curves
    for (auto i = 0; i < N; ++i) {
        SectionCurve2 curve = {
            .a = NaN,
            .b = NaN,
            .c = NaN,
            .xs = points[i].x,
            .xe = points[i+1].x
        };
        S.push_back(curve);
    }
    // c_i
    for (auto i = 0; i < N; ++i) {
        S[i].c = points[i].y;
    }
    // a_i
    number_t a0;
    if (opt_a0.has_value()) {
        a0 = opt_a0.value();
    } else {
        number_t Uc = 1;
        number_t Vc = 0;
        number_t Wc = 0;
        for (auto i = 0; i < N; ++i) {
            number_t Uj = 1;
            number_t Vj = 0;
            number_t Wj = 0;
            for (auto j = i-1; j > 0; --j) {
                // m => -1
                number_t dxj_0_m = points[j].x - points[j-1].x;
                number_t dxj_1_0 = points[j+1].x - points[j].x;
                number_t dyj_0_m = points[j].y - points[j-1].y;
                number_t dyj_1_0 = points[j+1].y - points[j].y;
                number_t uj_m = - dxj_0_m / dxj_1_0;
                number_t vj_m = (dyj_1_0 / dxj_1_0 - dyj_0_m / dxj_0_m) / dxj_1_0;

                number_t Uj_m = Uj * (uj_m * uj_m);
                number_t Vj_m = Uj * (2 * uj_m * vj_m) + Vj * (uj_m);
                number_t Wj_m = Uj * (vj_m * vj_m) + Vj * (vj_m) + Wj;
                Uj = Uj_m;
                Vj = Vj_m;
                Wj = Wj_m;
            }
            Uc += Uj;
            Vc += Vj;
            Wc += Wj;
        }
        a0 = -Vc / (2 * Uc);
    }
    S[0].a = a0;
    for (auto i = 0; i < N-1; ++i) {
        number_t dxi_1_0 = points[i+1].x - points[i].x;
        number_t dxi_2_1 = points[i+2].x - points[i+1].x;
        number_t dyi_1_0 = points[i+1].y - points[i].y;
        number_t dyi_2_1 = points[i+2].y - points[i+1].y;
        number_t ui = - dxi_1_0 / dxi_2_1;
        number_t vi = (dyi_2_1 / dxi_2_1 - dyi_1_0 / dxi_1_0) / dxi_2_1;
        // a_{i+1} = u_i a_i + v_i
        S[i+1].a = ui * S[i].a + vi;
    }
    // b_i
    for (auto i = 0; i < N; ++i) {
        number_t dxi_1_0 = points[i+1].x - points[i].x;
        number_t dyi_1_0 = points[i+1].y - points[i].y;
        S[i].b = dyi_1_0 / dxi_1_0 - S[i].a * dxi_1_0;
    }

    return S;
}

static number_t y_on_qspline(std::vector<SectionCurve2> const & curves, number_t x) {
    for (auto curve: curves) {
        if (x >= curve.xs && x <= curve.xe) {
            auto nx = x - curve.xs;
            return curve.a * nx * nx + curve.b * nx + curve.c;
        }
    }
    return NaN;
}
