// Copyright (c) 2025 - present KMX Systems. All rights reserved.
#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <kmx/gis.hpp>

namespace kmx::gis
{
    using Approx = Catch::Approx;

    using T = double;
    using K = constants<T>;
    using conv = conversion<T>;

    // Define tolerances for checks
    constexpr T stereo_tolerance_m = 1;
    constexpr T wgs84_lonlat_tolerance_deg = 0.0001;
    constexpr T wgs84_h_tolerance_m = 30;

#if 1
    TEST_CASE("Projection Origin Mapping", "[projection][origin]")
    {
        constexpr double origin_tolerance_m = conv::stereo70_origin_proj_tolerance_m;

        SECTION("Krasovsky Origin to Stereo70 False Origin")
        {
            geodetic_coord<T> kraso_origin = {.latitude = conv::stereo70_lat0_deg, .longitude = conv::stereo70_lon0_deg, .height = K::zero};

            // Use internal projection function directly for this test
            stereo70::coordinate<T> proj_origin = conv::krasovsky_geodetic_to_stereo70(kraso_origin);

            CHECK(proj_origin.x == Approx(conv::stereo70_fe).margin(origin_tolerance_m));
            CHECK(proj_origin.y == Approx(conv::stereo70_fn).margin(origin_tolerance_m));
        }

        SECTION("Stereo70 False Origin to Krasovsky Origin")
        {
            stereo70::coordinate<T> stereo_origin = {.x = conv::stereo70_fe, .y = conv::stereo70_fn, .z = K::zero};

            // Use internal inverse projection function directly
            geodetic_coord<T> kraso_origin_calc = conv::stereo70_to_krasovsky_geodetic(stereo_origin);

            CHECK(kraso_origin_calc.latitude == Approx(conv::stereo70_lat0_deg).margin(conv::stereo70_origin_coord_tolerance_deg));
            CHECK(kraso_origin_calc.longitude == Approx(conv::stereo70_lon0_deg).margin(conv::stereo70_origin_coord_tolerance_deg));
        }
    }
#endif

#if 1
    TEST_CASE("WGS84 <-> Stereo70 Conversions", "[conversion]")
    {
    #if 0
        // Print Helmert parameters being used for this test case
        const auto& params = conv::pulkovo58_to_wgs84_ellipsoid_params;
        std::cout << "\nINFO: Testing with Helmert Parameters (Pulkovo 1942 -> WGS84):\n";
        std::cout << std::fixed << std::setprecision(4);
        std::cout << " dx=" << params.dx << " m, dy=" << params.dy << " m, dz=" << params.dz << " m\n";
        std::cout << " rx=" << params.rx_sec << "\", ry=" << params.ry_sec << "\", rz=" << params.rz_sec << "\"\n";
        std::cout << " ds=" << params.ds_ppm << " ppm\n";
        std::cout << "Note: Rotation signs assumed Coordinate Frame, flipped internally for Position Vector calc.\n";
        std::cout << std::defaultfloat; // Reset float format
    #endif

        using input_expected_output_pair = std::tuple<wgs84::coordinate<T>, stereo70::coordinate<T>, const char*>;

        static constexpr std::array<input_expected_output_pair, 3u> test_data {
            input_expected_output_pair {{.latitude = 46.76952896129325, .longitude = 23.589875634659435, .height = 346},
                                        {.x = 392437.167, .y = 586510.612, .z = 346},
                                        "Cluj"},
            input_expected_output_pair {{.latitude = 45.79759722839506, .longitude = 24.15208824953577, .height = 428},
                                        {.x = 434216.523, .y = 477889.546, .z = 428},
                                        "Sibiu"},
            input_expected_output_pair {{.latitude = 45.64191842559299, .longitude = 25.588845714078214, .height = 594},
                                        {.x = 546017.285, .y = 460409.378, .z = 594},
                                        "Brasov"}};

        for (auto& tuple: test_data)
        {
            const auto& wgs_in = std::get<0>(tuple);
            const auto& stereo_expected = std::get<1>(tuple);
            INFO("coords from " << std::get<2>(tuple));

            const auto stereo_calc = wgs84_to_stereo70(wgs_in);
            CHECK(std::abs(stereo_calc.x - stereo_expected.x) < stereo_tolerance_m);
            CHECK(std::abs(stereo_calc.y - stereo_expected.y) < stereo_tolerance_m);

            const auto wgs_rt = stereo70_to_wgs84(stereo_calc);
            CHECK(std::abs(wgs_rt.latitude - wgs_in.latitude) < wgs84_lonlat_tolerance_deg);
            CHECK(std::abs(wgs_rt.longitude - wgs_in.longitude) < wgs84_lonlat_tolerance_deg);
            CHECK(std::abs(wgs_rt.height - wgs_in.height) < wgs84_h_tolerance_m);
        }
    }
#endif

    TEST_CASE("Stereo70 -> WGS84 Conversions", "[conversion]")
    {
        using input_expected_output_pair = std::tuple<stereo70::coordinate<T>, wgs84::coordinate<T>, const char*>;

        static constexpr std::array<input_expected_output_pair, 3u> test_data {
            input_expected_output_pair {
                {.x = 392434.50, .y = 586512.00, .z = 346}, {.latitude = 46.769533, .longitude = 23.589875, .height = 346}, "Cluj"},
            input_expected_output_pair {
                {.x = 434218, .y = 477894, .z = 428}, {.latitude = 45.797629, .longitude = 24.152137, .height = 428}, "Sibiu"},
            input_expected_output_pair {
                {.x = 546029.50, .y = 460415.00, .z = 594}, {.latitude = 45.641958, .longitude = 25.589031, .height = 594}, "Brasov"}};

        for (auto& tuple: test_data)
        {
            const auto& stereo_in = std::get<0>(tuple);
            const auto& wgs_expected = std::get<1>(tuple);
            INFO("coords from " << std::get<2>(tuple));

            const auto wgs_rt = stereo70_to_wgs84(stereo_in);
            INFO("lat=" << wgs_rt.latitude << " lon=" << wgs_rt.longitude << " alt=" << wgs_rt.height);
            CHECK(std::abs(wgs_rt.latitude - wgs_expected.latitude) < wgs84_lonlat_tolerance_deg);
            CHECK(std::abs(wgs_rt.longitude - wgs_expected.longitude) < wgs84_lonlat_tolerance_deg);
            CHECK(std::abs(wgs_rt.height - wgs_expected.height) < wgs84_h_tolerance_m);
        }
    }
}
