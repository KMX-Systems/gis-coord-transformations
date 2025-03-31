// Copyright (c) 2025 - present KMX Systems. All rights reserved.
#define CATCH_CONFIG_MAIN
#include <array>
#include <catch2/catch_all.hpp>
// #include <iostream>
#include <kmx/gis.hpp>
#include <tuple>

namespace kmx::gis
{
    using T = double;
    using K = constants<T>;
    using conv = conversion<T>;
    using stereo70_params = stereo70::projection_params<T>;

    // Define tolerances for checks
    constexpr T stereo_tolerance_m = 0.001;
    constexpr T wgs84_lonlat_tolerance_deg = 0.0001;
    constexpr T wgs84_h_tolerance_m = 30;

    TEST_CASE("Projection Origin Mapping", "[projection][origin]")
    {
        constexpr T origin_tolerance_m = stereo70_params::origin_proj_tol;

        SECTION("Krasovsky Origin to Stereo70 False Origin")
        {
            geodetic_coord<T> kraso_origin = {.latitude = stereo70_params::lat0_deg, .longitude = stereo70_params::lon0_deg, .height = {}};

            // Use internal projection function directly for this test
            stereo70::coordinate<T> proj_origin = conv::krasovsky_geodetic_to_stereo70(kraso_origin);

            CHECK(std::abs(proj_origin.x - stereo70_params::fe) < origin_tolerance_m);
            CHECK(std::abs(proj_origin.y - stereo70_params::fn) < origin_tolerance_m);
        }

        SECTION("Stereo70 False Origin to Krasovsky Origin")
        {
            stereo70::coordinate<T> stereo_origin = {.x = stereo70_params::fe, .y = stereo70_params::fn, .z = {}};

            // Use internal inverse projection function directly
            geodetic_coord<T> kraso_origin_calc = conv::stereo70_to_krasovsky_geodetic(stereo_origin);

            CHECK(std::abs(kraso_origin_calc.longitude - stereo70_params::lat0_deg) < stereo70_params::origin_tol_deg);
            CHECK(std::abs(kraso_origin_calc.longitude - stereo70_params::lon0_deg) < stereo70_params::origin_tol_deg);
        }
    }

    TEST_CASE("WGS84 <-> Stereo70 Conversions", "[conversion 1]")
    {
#if 0
        // Print Helmert parameters being used for this test case
        const auto& params = conv::pulkovo58_to_wgs84;
        std::cout << "\nINFO: Testing with Helmert Parameters (Pulkovo 1942 -> WGS84):\n";
        std::cout << std::fixed << std::setprecision(4);
        std::cout << " dx=" << params.dx << " m, dy=" << params.dy << " m, dz=" << params.dz << " m\n";
        std::cout << " rx=" << params.rx_sec << "\", ry=" << params.ry_sec << "\", rz=" << params.rz_sec << "\"\n";
        std::cout << " ds=" << params.ds_ppm << " ppm\n";
        std::cout << std::defaultfloat; // Reset float format
#endif

        using input_expected_output_pair = std::tuple<wgs84::coordinate<T>, stereo70::coordinate<T>, const char*>;

        static constexpr std::array test_data {
            input_expected_output_pair {{.latitude = 46.76952896129325, .longitude = 23.589875634659435, .height = 346},
                                        {.x = 392437.167, .y = 586510.612, .z = 346},
                                        "Cluj"},
            input_expected_output_pair {{.latitude = 45.79759722839506, .longitude = 24.15208824953577, .height = 428},
                                        {.x = 434216.523, .y = 477889.546, .z = 428},
                                        "Sibiu"},
            input_expected_output_pair {{.latitude = 45.64191842559299, .longitude = 25.588845714078214, .height = 594},
                                        {.x = 546017.285, .y = 460409.378, .z = 594},
                                        "Brasov"}};

        for (const auto& [wgs_in, stereo_expected, location]: test_data)
        {
            INFO("Testing coordinates from " << location);

            SECTION("Forward conversion (WGS84 -> Stereo70)")
            {
                const auto stereo_calc = conv::wgs84_to_stereo70(wgs_in);
                CHECK(std::abs(stereo_calc.x - stereo_expected.x) < stereo_tolerance_m);
                CHECK(std::abs(stereo_calc.y - stereo_expected.y) < stereo_tolerance_m);
                CHECK(std::abs(stereo_calc.z - stereo_expected.z) < wgs84_h_tolerance_m);
            }

            SECTION("Round-trip conversion (WGS84 -> Stereo70 -> WGS84)")
            {
                const auto stereo_calc = conv::wgs84_to_stereo70(wgs_in);
                const auto wgs_rt = conv::stereo70_to_wgs84(stereo_calc);

                CHECK(std::abs(wgs_rt.latitude - wgs_in.latitude) < wgs84_lonlat_tolerance_deg);
                CHECK(std::abs(wgs_rt.latitude - wgs_in.latitude) < wgs84_lonlat_tolerance_deg);
                CHECK(std::abs(wgs_rt.height - wgs_in.height) < wgs84_h_tolerance_m);
            }
        }
    }

    TEST_CASE("Stereo70 -> WGS84 Conversions", "[conversion 2]")
    {
        using input_expected_output_pair = std::tuple<stereo70::coordinate<T>, wgs84::coordinate<T>, const char*>;

        static constexpr std::array test_data {
            input_expected_output_pair {
                {.x = 392434.50, .y = 586512.00, .z = 346}, {.latitude = 46.769533, .longitude = 23.589875, .height = 346}, "Cluj"},
            input_expected_output_pair {
                {.x = 434218, .y = 477894, .z = 428}, {.latitude = 45.797629, .longitude = 24.152137, .height = 428}, "Sibiu"},
            input_expected_output_pair {
                {.x = 546029.50, .y = 460415.00, .z = 594}, {.latitude = 45.641958, .longitude = 25.589031, .height = 594}, "Brasov"}};

        for (const auto& [stereo_in, wgs_expected, location]: test_data)
        {
            INFO("Testing coordinates from " << location);

            const auto wgs_calc = conv::stereo70_to_wgs84(stereo_in);

            INFO("Calculated: lat=" << wgs_calc.latitude << " lon=" << wgs_calc.longitude << " alt=" << wgs_calc.height);

            CHECK(std::abs(wgs_calc.latitude - wgs_expected.latitude) < wgs84_lonlat_tolerance_deg);
            CHECK(std::abs(wgs_calc.latitude - wgs_expected.latitude) < wgs84_lonlat_tolerance_deg);
            CHECK(std::abs(wgs_calc.height - wgs_expected.height) < wgs84_h_tolerance_m);
        }
    }

    TEST_CASE("Edge Cases", "[conversion][edge]")
    {
        SECTION("Invalid Input Handling")
        {
            // Test NaN inputs
            wgs84::coordinate<T> nan_coord {.latitude = std::numeric_limits<T>::quiet_NaN(), .longitude = 0, .height = 0};

            REQUIRE_THROWS_AS(conv::wgs84_to_stereo70(nan_coord), std::invalid_argument);

            // Test out of range latitude
            wgs84::coordinate<T> invalid_lat {.latitude = 100, .longitude = 0, .height = 0};

            REQUIRE_THROWS_AS(conv::wgs84_to_stereo70(invalid_lat), std::invalid_argument);
        }

        SECTION("Pole Handling")
        {
            // North Pole
            wgs84::coordinate<T> north_pole {.latitude = 90, .longitude = 0, .height = 0};

            REQUIRE_NOTHROW(conv::wgs84_to_stereo70(north_pole));

            // South Pole
            wgs84::coordinate<T> south_pole {.latitude = -90, .longitude = 0, .height = 0};

            REQUIRE_NOTHROW(conv::wgs84_to_stereo70(south_pole));
        }
    }
}
