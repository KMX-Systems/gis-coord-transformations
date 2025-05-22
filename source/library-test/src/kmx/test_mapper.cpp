// Copyright (c) 2025 - present KMX Systems. All rights reserved.
#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <cmath>   // For std::floor, std::abs, std::isfinite
#include <cstdint> // For std::uint8_t, etc.
#include <format>  // For std::format, std::vformat, std::make_format_args, std::format_to C++20/23
#include <kmx/gis/uniform_grid_mapper.hpp>
#include <optional>
#include <stdexcept>
#include <string>        // For std::string
#include <type_traits>   // For std::is_floating_point_v, std::is_integral_v, std::is_unsigned_v
#include <unordered_map> // For testing
#include <utility>       // For std::pair

namespace kmx::gis
{
    // Test constants
    static constexpr double default_x_min_c = 280000.0;
    static constexpr double default_x_max_c = 840000.0;
    static constexpr double default_y_min_c = 230000.0;
    static constexpr double default_y_max_c = 850000.0;

    static constexpr coordinate::xy<double> default_min_corner_c(default_x_min_c, default_y_min_c);
    static constexpr coordinate::xy<double> default_max_corner_c(default_x_max_c, default_y_max_c);

    static constexpr std::uint8_t default_row_count_u8_c = 5u;
    static constexpr std::uint8_t default_col_count_u8_c = 10u;

    // For testing uniform_grid_mapper, which is in kmx::gis
    TEST_CASE("kmx::gis::uniform_grid_mapper functionality", "[uniform_grid_mapper]")
    {
        using test_point_type = coordinate::xy<double>;
        using test_mapper_defaults = uniform_grid_mapper<test_point_type, std::uint8_t, double>;
        using test_cell_type = typename test_mapper_defaults::cell_type;

        const test_point_type p0_0_d(0.0, 0.0);
        const test_point_type p10_10_d(10.0, 10.0);
        const test_point_type p0_10_d(0.0, 10.0);
        const test_point_type neg_width_min_d(10.0, 0.0);
        const test_point_type neg_width_max_d(0.0, 10.0);
        const test_point_type neg_height_min_d(0.0, 10.0);
        const test_point_type neg_height_max_d(0.0, 0.0);

        SECTION("Constructor validation (defaults)")
        {
            REQUIRE_NOTHROW(test_mapper_defaults(p0_0_d, p10_10_d, 1u, 1u));
            REQUIRE_THROWS_AS(test_mapper_defaults(p0_0_d, p10_10_d, 0u, 1u), std::invalid_argument);
            REQUIRE_THROWS_AS(test_mapper_defaults(p0_0_d, p10_10_d, 1u, 0u), std::invalid_argument);
            REQUIRE_THROWS_AS(test_mapper_defaults(neg_width_min_d, neg_width_max_d, 1u, 1u), std::invalid_argument);
            REQUIRE_THROWS_AS(test_mapper_defaults(neg_height_max_d, neg_height_min_d, 1u, 1u), std::invalid_argument);
            REQUIRE_THROWS_AS(test_mapper_defaults(p0_0_d, p0_10_d, 1u, 1u), std::invalid_argument);
        }

        const test_mapper_defaults mapper_defaults(default_min_corner_c, default_max_corner_c, default_row_count_u8_c,
                                                   default_col_count_u8_c);

        SECTION("Size method (defaults)")
        {
            const auto dims = mapper_defaults.size();
            REQUIRE(dims.first == default_row_count_u8_c);
            REQUIRE(dims.second == default_col_count_u8_c);
        }

        SECTION("Points within grid boundaries (defaults)")
        {
            REQUIRE(mapper_defaults[default_min_corner_c].value() == test_cell_type(0, 0));
            REQUIRE(mapper_defaults[default_max_corner_c].value() ==
                    test_cell_type(default_col_count_u8_c - 1, default_row_count_u8_c - 1));

            const test_point_type top_left(default_min_corner_c.x, default_max_corner_c.y);
            REQUIRE(mapper_defaults[top_left].value() == test_cell_type(0, default_row_count_u8_c - 1));

            const test_point_type bottom_right(default_max_corner_c.x, default_min_corner_c.y);
            REQUIRE(mapper_defaults[bottom_right].value() == test_cell_type(default_col_count_u8_c - 1, 0));

            const double cell_width = (default_max_corner_c.x - default_min_corner_c.x) / static_cast<double>(default_col_count_u8_c);
            const double cell_height = (default_max_corner_c.y - default_min_corner_c.y) / static_cast<double>(default_row_count_u8_c);
            const test_point_type inside_first_cell(default_min_corner_c.x + (0.1 * cell_width),
                                                    default_min_corner_c.y + (0.1 * cell_height));
            REQUIRE(mapper_defaults[inside_first_cell].value() == test_cell_type(0, 0));

            const test_point_type center_pt((default_min_corner_c.x + default_max_corner_c.x) / 2.0,
                                            (default_min_corner_c.y + default_max_corner_c.y) / 2.0);
            REQUIRE(mapper_defaults[center_pt].value() == test_cell_type(default_col_count_u8_c / 2, default_row_count_u8_c / 2));
        }

        SECTION("Points on internal grid lines (defaults)")
        {
            const double cell_width = (default_max_corner_c.x - default_min_corner_c.x) / static_cast<double>(default_col_count_u8_c);
            const double cell_height = (default_max_corner_c.y - default_min_corner_c.y) / static_cast<double>(default_row_count_u8_c);

            const test_point_type on_line1(default_min_corner_c.x + cell_width, default_min_corner_c.y + (0.5 * cell_height));
            REQUIRE(mapper_defaults[on_line1].value() == test_cell_type(1, 0));
            const test_point_type on_line2(default_min_corner_c.x + (0.5 * cell_width), default_min_corner_c.y + cell_height);
            REQUIRE(mapper_defaults[on_line2].value() == test_cell_type(0, 1));
            const test_point_type on_line3(default_min_corner_c.x + cell_width, default_min_corner_c.y + cell_height);
            REQUIRE(mapper_defaults[on_line3].value() == test_cell_type(1, 1));
        }

        SECTION("Points outside grid boundaries (defaults)")
        {
            REQUIRE_FALSE(mapper_defaults[coordinate::xy<double>(default_min_corner_c.x - 100.0, default_min_corner_c.y)].has_value());
            REQUIRE_FALSE(mapper_defaults[coordinate::xy<double>(default_max_corner_c.x + 100.0, default_max_corner_c.y)].has_value());
        }

        SECTION("Points very close to boundaries (within epsilon) (defaults)")
        {
            const double eps_test = 1e-10;
            REQUIRE(mapper_defaults[coordinate::xy<double>(default_min_corner_c.x - eps_test, default_min_corner_c.y)].value() ==
                    test_cell_type(0, 0));
            REQUIRE(mapper_defaults[coordinate::xy<double>(default_max_corner_c.x + eps_test, default_max_corner_c.y)].value() ==
                    test_cell_type(default_col_count_u8_c - 1, default_row_count_u8_c - 1));
        }

        SECTION("1x1 Grid (defaults)")
        {
            const test_mapper_defaults mapper_1x1(p0_0_d, {100.0, 50.0}, 1u, 1u);
            REQUIRE(mapper_1x1[p0_0_d].value() == test_cell_type(0, 0));
            REQUIRE(mapper_1x1[{50.0, 25.0}].value() == test_cell_type(0, 0));
            REQUIRE(mapper_1x1[{100.0, 50.0}].value() == test_cell_type(0, 0));
            REQUIRE_FALSE(mapper_1x1[{100.1, 50.0}].has_value());
        }

        using test_point_float_type = coordinate::xy<float>;
        using test_mapper_float_float = uniform_grid_mapper<test_point_float_type, std::uint8_t, float>;
        using test_cell_float_u8 = typename test_mapper_float_float::cell_type;

        const test_point_float_type min_corner_f(static_cast<float>(default_min_corner_c.x), static_cast<float>(default_min_corner_c.y));
        const test_point_float_type max_corner_f(static_cast<float>(default_max_corner_c.x), static_cast<float>(default_max_corner_c.y));

        SECTION("Functionality with float coords and float intermediate FP (float, uint8_t, float)")
        {
            const test_mapper_float_float mapper_float_processing(min_corner_f, max_corner_f, default_row_count_u8_c,
                                                                  default_col_count_u8_c);
            REQUIRE(mapper_float_processing[min_corner_f].value() == test_cell_float_u8(0, 0));
            const test_point_float_type on_boundary_f(
                min_corner_f.x + (max_corner_f.x - min_corner_f.x) / static_cast<float>(default_col_count_u8_c), min_corner_f.y);
            REQUIRE(mapper_float_processing[on_boundary_f].value() == test_cell_float_u8(1, 0));
        }

        SECTION("Usage in std::unordered_map")
        {
            std::unordered_map<test_cell_type, std::string> cell_data;

            test_cell_type cell_a(5, 1);
            test_cell_type cell_b(default_col_count_u8_c - 1, default_row_count_u8_c - 1);

            cell_data[cell_a] = "Data for cell A";
            cell_data[cell_b] = "Data for cell B";

            REQUIRE(cell_data.count(cell_a) == 1);
            REQUIRE(cell_data[cell_a] == "Data for cell A");
        }
    }

// Updated guard for constexpr tests:
// Primary needs are constexpr cmath (for std::floor) and is_constant_evaluated.
// Format is used by xy's printing, but not directly by the mapper's core constexpr logic being tested by static_asserts on operator[].
// However, if a static_assert tried to print an xy (e.g. in a custom message), then format's constexpr-ness would matter.
// For robustness, we can keep it, or make the guard even more specific to the features used in the static_assert block.
#if defined(__cpp_lib_is_constant_evaluated) && __cpp_lib_is_constant_evaluated >= 201811L && defined(__cpp_lib_constexpr_cmath) && \
    __cpp_lib_constexpr_cmath >= 202202L && defined(__cpp_lib_format) && __cpp_lib_format >= 201907L

    TEST_CASE("kmx::gis::uniform_grid_mapper constexpr evaluation", "[uniform_grid_mapper][constexpr]")
    {
        using test_point_type_cx = coordinate::xy<double>;
        using test_index_comp_type_cx = std::uint8_t;
        // Explicitly specify all template arguments for clarity in test
        using test_mapper_type_cx = uniform_grid_mapper<test_point_type_cx, test_index_comp_type_cx, double>;
        using test_cell_type_cx = typename test_mapper_type_cx::cell_type;

        static constexpr test_point_type_cx min_c(0.0, 0.0);
        static constexpr test_point_type_cx max_c(100.0, 50.0);
        static constexpr test_index_comp_type_cx row_count_cx = 5;
        static constexpr test_index_comp_type_cx col_count_cx = 10;

        static constexpr test_mapper_type_cx constexpr_mapper(min_c, max_c, row_count_cx, col_count_cx);

        static_assert(constexpr_mapper.size().first == row_count_cx);
        static_assert(constexpr_mapper.size().second == col_count_cx);

        static constexpr test_point_type_cx pt1_cx(0.0, 0.0);
        static constexpr std::optional<test_cell_type_cx> res1_cx = constexpr_mapper[pt1_cx];
        static_assert(res1_cx.has_value());
        static_assert(res1_cx.value().x == 0 && res1_cx.value().y == 0);
        static_assert(res1_cx.value() == test_cell_type_cx(0, 0));

        static constexpr test_point_type_cx pt2_cx(100.0, 50.0);
        static constexpr auto res2_cx = constexpr_mapper[pt2_cx];
        static_assert(res2_cx.has_value());
        static_assert(res2_cx.value().x == (col_count_cx - 1) && res2_cx.value().y == (row_count_cx - 1));

        static constexpr test_cell_type_cx cell_idx_a(1, 2);
        static constexpr test_cell_type_cx cell_idx_b(1, 3);
        static_assert(cell_idx_a < cell_idx_b);

        // Runtime check to ensure the test case itself executes if compiled
        const auto runtime_res1 = constexpr_mapper[pt1_cx];
        REQUIRE(runtime_res1.has_value());
        REQUIRE(runtime_res1.value() == test_cell_type_cx(0, 0));
    }
#else
    TEST_CASE("kmx::gis::uniform_grid_mapper constexpr evaluation (SKIPPED)", "[uniform_grid_mapper][constexpr][!shouldfail]")
    {
        std::string reason = "Skipping uniform_grid_mapper constexpr tests due to missing features: ";
        bool first = true;
        auto append_reason = [&](bool defined, long long val, long long required, const char* name)
        {
            if (!defined || (required > 0 && val < required))
            {
                if (!first)
                    reason += ", ";
                reason += name;
                if (defined)
                    reason += " (val=" + std::to_string(val) + " req>=" + std::to_string(required) + ")";
                else
                    reason += " (not defined)";
                first = false;
            }
        };

    #ifdef __cpp_lib_is_constant_evaluated
        append_reason(true, __cpp_lib_is_constant_evaluated, 201811L, "__cpp_lib_is_constant_evaluated");
    #else
        append_reason(false, 0, 201811L, "__cpp_lib_is_constant_evaluated");
    #endif

    #ifdef __cpp_lib_constexpr_cmath
        append_reason(true, __cpp_lib_constexpr_cmath, 202202L, "__cpp_lib_constexpr_cmath");
    #else
        append_reason(false, 0, 202202L, "__cpp_lib_constexpr_cmath");
    #endif

    #ifdef __cpp_lib_format
        append_reason(true, __cpp_lib_format, 201907L, "__cpp_lib_format"); // Basic format needed for xy printing if used in static_asserts
    #else
        append_reason(false, 0, 201907L, "__cpp_lib_format");
    #endif

        if (first)
            reason += "Unknown (all checked macros seem defined as expected)."; // Should not happen if skip occurs

        WARN(reason);
    }
#endif

} // namespace kmx::gis
