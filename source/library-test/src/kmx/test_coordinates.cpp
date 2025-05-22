// Copyright (c) 2025 - present KMX Systems. All rights reserved.
#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <cmath>      // For std::floor, std::abs, std::isfinite
#include <cstdint>    // For std::uint8_t, etc.
#include <functional> // For std::hash
#include <kmx/gis/coordinate/base.hpp>
#include <limits> // For std::numeric_limits
#include <stdexcept>

namespace kmx::gis::coordinate
{
    TEST_CASE("kmx::gis::coordinate::xy functionality (generic)", "[coordinate][xy]")
    {
        SECTION("Floating point xy (double)")
        {
            using coord_d = xy<double>;
            const coord_d p_zero;
            REQUIRE(p_zero.x == 0.0);
            REQUIRE(p_zero.y == 0.0);

            const coord_d p1(1.23456, 7.89123);
            REQUIRE(p1.x == Catch::Approx(1.23456));
            REQUIRE(p1.y == Catch::Approx(7.89123));

            REQUIRE_NOTHROW(p1.validate());
            const coord_d p_nan(std::numeric_limits<double>::quiet_NaN(), 1.0);
            REQUIRE_THROWS_AS(p_nan.validate(), std::invalid_argument);

            REQUIRE(p1.to_string_simple(2) == "1.23, 7.89");
            std::stringstream ss;
            ss << p1;
            REQUIRE(ss.str() == "xy(1.23456, 7.89123)");

            const coord_d p_add1(1.0, 2.0), p_add2(3.0, 0.5);
            REQUIRE(p_add1 + p_add2 == xy<double>(4.0, 2.5));
            REQUIRE(p_add1 * 2.0 == xy<double>(2.0, 4.0));

            const coord_d p_eq1(1.0, 2.0);
            const coord_d p_eq2(1.0, 2.0000000001);
            const coord_d p_eq3(1.0, 2.001);
            REQUIRE(p_eq1 == p_eq2);
            REQUIRE_FALSE(p_eq1 == p_eq3);

            const coord_d p_ord1(1.0, 2.0);
            const coord_d p_ord2(1.0, 3.0);
            const coord_d p_ord3(2.0, 1.0);
            REQUIRE(p_ord1 < p_ord2);
            REQUIRE(p_ord2 < p_ord3);
        }

        SECTION("Integral xy (uint8_t)")
        {
            using coord_u8 = xy<std::uint8_t>;
            const coord_u8 c_zero;
            REQUIRE(c_zero.x == 0);
            REQUIRE(c_zero.y == 0);

            const coord_u8 c1(10, 20);
            REQUIRE(c1.x == 10);
            REQUIRE(c1.y == 20);

            REQUIRE_NOTHROW(c1.validate());

            REQUIRE(c1.to_string_simple() == "10, 20");
            std::stringstream ss;
            ss << c1;
            REQUIRE(ss.str() == "xy(10, 20)");

            const coord_u8 c_add1(5, 3), c_add2(2, 7);
            REQUIRE(c_add1 + c_add2 == xy<std::uint8_t>(7, 10));
            REQUIRE(c_add1 * static_cast<std::uint8_t>(2) == xy<std::uint8_t>(10, 6));
            REQUIRE(c1 / static_cast<std::uint8_t>(3) == xy<std::uint8_t>(3, 6));
            REQUIRE_THROWS_AS(c1 / static_cast<std::uint8_t>(0), std::runtime_error);

            const coord_u8 c_eq1(10, 20);
            const coord_u8 c_eq2(10, 20);
            const coord_u8 c_eq3(11, 20);
            REQUIRE(c_eq1 == c_eq2);
            REQUIRE_FALSE(c_eq1 == c_eq3);
            REQUIRE(c_eq1 < c_eq3);
        }

        SECTION("Hashing xy")
        {
            using coord_d = xy<double>;
            using coord_u8 = xy<std::uint8_t>;

            std::hash<coord_d> hash_d;
            std::hash<coord_u8> hash_u8;

            const coord_d p1_d(1.0, 2.0);
            const coord_d p1_d_exact_copy = p1_d;
            REQUIRE(hash_d(p1_d) == hash_d(p1_d_exact_copy));

            const coord_u8 c1_u8(10, 20);
            const coord_u8 c2_u8(10, 20);
            REQUIRE(hash_u8(c1_u8) == hash_u8(c2_u8));
        }
    }

    TEST_CASE("kmx::gis::coordinate::xyz functionality", "[coordinate][xyz]")
    {
        using coord3d_d = xyz<double>; // Alias for xyz<double>
        using coord3d_f = xyz<float>;  // Alias for xyz<float>

        SECTION("Default constructor")
        {
            const coord3d_d p;
            REQUIRE(p.x == 0.0);
            REQUIRE(p.y == 0.0);
            REQUIRE(p.z == 0.0);
        }

        SECTION("Parameterized constructor")
        {
            const coord3d_d p(1.1, 2.2, 3.3);
            REQUIRE(p.x == Catch::Approx(1.1));
            REQUIRE(p.y == Catch::Approx(2.2));
            REQUIRE(p.z == Catch::Approx(3.3));
        }

        SECTION("Validation")
        {
            const coord3d_d p_valid(1.0, 2.0, 3.0);
            REQUIRE_NOTHROW(p_valid.validate());

            const coord3d_d p_nan_x(std::numeric_limits<double>::quiet_NaN(), 1.0, 2.0);
            REQUIRE_THROWS_AS(p_nan_x.validate(), std::invalid_argument);

            const coord3d_d p_nan_y(1.0, std::numeric_limits<double>::quiet_NaN(), 2.0);
            REQUIRE_THROWS_AS(p_nan_y.validate(), std::invalid_argument);

            const coord3d_d p_nan_z(1.0, 2.0, std::numeric_limits<double>::quiet_NaN());
            REQUIRE_THROWS_AS(p_nan_z.validate(), std::invalid_argument);

            const coord3d_d p_inf_x(std::numeric_limits<double>::infinity(), 1.0, 2.0);
            REQUIRE_THROWS_AS(p_inf_x.validate(), std::invalid_argument);
        }

        SECTION("String formatting (to_string_simple, print, println, operator<<)")
        {
            const coord3d_d p(1.234567, 8.9, -0.54321);

            REQUIRE(p.to_string_simple(2) == "1.23, 8.90, -0.54");          // Precision 2
            REQUIRE(p.to_string_simple(5) == "1.23457, 8.90000, -0.54321"); // Default precision 5, with rounding

            std::stringstream ss_print;
            p.print(ss_print, 3);
            REQUIRE(ss_print.str() == "1.235, 8.900, -0.543");

            std::stringstream ss_println;
            p.println(ss_println, 1);
            REQUIRE(ss_println.str() == "1.2, 8.9, -0.5\n"); // Includes newline

            std::stringstream ss_op_shift;
            ss_op_shift << p; // Uses friend operator<<, which formats as "xyz(x, y, z)" with precision 5
            REQUIRE(ss_op_shift.str() == "xyz(1.23457, 8.90000, -0.54321)");
        }

        SECTION("Spaceship operator and equality")
        {
            const coord3d_d p1(1.0, 2.0, 3.0);
            const coord3d_d p2(1.0, 2.0, 3.0);
            const coord3d_d p3_x(2.0, 2.0, 3.0);
            const coord3d_d p4_y(1.0, 3.0, 3.0);
            const coord3d_d p5_z(1.0, 2.0, 4.0);

            // Equality (defaulted by spaceship, so bitwise for floats)
            REQUIRE(p1 == p2);
            REQUIRE((p1 <=> p2) == std::strong_ordering::equal);

            REQUIRE(p1 != p3_x);
            REQUIRE(p1 != p4_y);
            REQUIRE(p1 != p5_z);

            // Ordering (lexicographical)
            REQUIRE(p1 < p3_x); // 1.0 < 2.0 (x differs)
            REQUIRE(p1 < p4_y); // x same, 2.0 < 3.0 (y differs)
            REQUIRE(p1 < p5_z); // x, y same, 3.0 < 4.0 (z differs)

            REQUIRE(p3_x > p1);
            REQUIRE(p4_y > p1);
            REQUIRE(p5_z > p1);

            REQUIRE(p1 <= p2);
            REQUIRE(p1 <= p5_z);
            REQUIRE(p5_z >= p1);
        }

        SECTION("Arithmetic operators")
        {
            const coord3d_d p1(1.0, 2.0, 3.0);
            const coord3d_d p2(0.5, 1.5, 2.5);
            coord3d_d p_mut = p1;

            // Addition
            REQUIRE(p1 + p2 == coord3d_d(1.5, 3.5, 5.5));
            p_mut = p1;
            p_mut += p2;
            REQUIRE(p_mut == coord3d_d(1.5, 3.5, 5.5));

            // Subtraction
            REQUIRE(p1 - p2 == coord3d_d(0.5, 0.5, 0.5));
            p_mut = p1;
            p_mut -= p2;
            REQUIRE(p_mut == coord3d_d(0.5, 0.5, 0.5));

            // Scalar Multiplication
            REQUIRE(p1 * 2.0 == coord3d_d(2.0, 4.0, 6.0));
            REQUIRE(2.5 * p1 == coord3d_d(2.5, 5.0, 7.5)); // Non-member
            p_mut = p1;
            p_mut *= 0.5;
            REQUIRE(p_mut == coord3d_d(0.5, 1.0, 1.5));

            // Scalar Division
            REQUIRE(p1 / 2.0 == coord3d_d(0.5, 1.0, 1.5));
            REQUIRE_THROWS_AS(p1 / 0.0, std::runtime_error);
            p_mut = p1;
            p_mut /= 4.0;
            REQUIRE(p_mut == coord3d_d(0.25, 0.5, 0.75));
            REQUIRE_THROWS_AS(p_mut /= 0.0, std::runtime_error);

            // Unary Minus
            REQUIRE(-p1 == coord3d_d(-1.0, -2.0, -3.0));
        }

        SECTION("Float type xyz")
        {
            const coord3d_f p1_f(1.0f, 2.5f, -0.5f);
            const coord3d_f p2_f(0.5f, 0.5f, 0.5f);

            REQUIRE(p1_f + p2_f == coord3d_f(1.5f, 3.0f, 0.0f));
            REQUIRE(p1_f * 2.0f == coord3d_f(2.0f, 5.0f, -1.0f));

            std::stringstream ss;
            ss << p1_f;
            // Check if formatting for float is reasonable
            // Default precision is 5 for operator<< in xyz
            // Exact string might depend on float to string conversion nuances
            // For "1.0f, 2.5f, -0.5f" with 5 places:
            REQUIRE(ss.str() == "xyz(1.00000, 2.50000, -0.50000)");
        }
    }
}
