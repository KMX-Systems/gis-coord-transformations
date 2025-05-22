// Copyright (c) 2025 - present KMX Systems. All rights reserved.
#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <kmx/gis.hpp>
#include <random>

namespace kmx::gis
{
    using T = double;
    using K = constants<T>;
    using conv = coordinate::conversion<T>;
    using stereo70_params = stereo70::projection_params<T>;

    TEST_CASE("Projection Origin Mapping", "[projection][origin]")
    {
        constexpr T origin_tolerance_m = stereo70_params::origin_proj_tol;

        SECTION("Krasovsky Origin to Stereo70 False Origin")
        {
            coordinate::geodetic<T> kraso_origin = {
                .latitude = stereo70_params::lat0_deg, .longitude = stereo70_params::lon0_deg, .altitude = {}};

            // Use internal projection function directly for this test
            const auto proj_origin = conv::krasovsky_geodetic_to_stereo70(kraso_origin);

            CHECK(std::abs(proj_origin.x - stereo70_params::fe) < origin_tolerance_m);
            CHECK(std::abs(proj_origin.y - stereo70_params::fn) < origin_tolerance_m);
        }

        SECTION("Stereo70 False Origin to Krasovsky Origin")
        {
            stereo70::coordinate<T> stereo_origin = stereo70::coordinate<T> {stereo70_params::fe, stereo70_params::fn, 0};

            // Use internal inverse projection function directly
            coordinate::geodetic<T> kraso_origin_calc = conv::stereo70_to_krasovsky_geodetic(stereo_origin);

            CHECK(std::abs(kraso_origin_calc.latitude - stereo70_params::lat0_deg) < stereo70_params::origin_tol_deg);
            CHECK(std::abs(kraso_origin_calc.longitude - stereo70_params::lon0_deg) < stereo70_params::origin_tol_deg);
        }
    }

    void round_trip_test(const helmert_params<T>& params, const T wgs84_lonlat_tolerance_deg, const T altitude_tolerance_m,
                         const char* test_name)
    {
        static constexpr double romania_min_lat = 43.62; // Aproximativ sud (zona Vama Veche)
        static constexpr double romania_max_lat = 48.26; // Aproximativ nord (zona Horodiștea)
        static constexpr double romania_min_lon = 20.26; // Aproximativ vest (zona Beba Veche)
        static constexpr double romania_max_lon = 29.69; // Aproximativ est (zona Sulina)

        std::random_device rd;
        static std::mt19937 generator(rd());
        std::uniform_real_distribution<double> latitude_distribution(romania_min_lat, romania_max_lat);
        std::uniform_real_distribution<double> longitude_distribution(romania_min_lon, romania_max_lon);

        SECTION(test_name)
        {
            for (unsigned i = 1000u; i > 0u; --i)
            {
                const wgs84::coordinate<T> wgs_input {latitude_distribution(generator), longitude_distribution(generator), 0.0};

                const auto stereo_calc = conv::wgs84_to_stereo70(wgs_input, params);
                const auto wgs_output = conv::stereo70_to_wgs84(stereo_calc, params);
                const auto [deltaLat, deltaLon, deltaAlt] = diff(wgs_input, wgs_output);

                CHECK(deltaLat < wgs84_lonlat_tolerance_deg);
                CHECK(deltaLon < wgs84_lonlat_tolerance_deg);
                CHECK(deltaAlt < altitude_tolerance_m);
            }
        }
    }

    TEST_CASE("Round-trip conversion 1 (WGS84 -> Stereo70 -> WGS84)", "[epsg1241]")
    {
        const auto& params = conv::transformation::epsg1241;
        static constexpr T wgs84_lonlat_tolerance_deg = T(1e-11);
        static constexpr T altitude_tolerance_m = T(1e-6);

        round_trip_test(params, wgs84_lonlat_tolerance_deg, altitude_tolerance_m, "[epsg1241]");
    }

    TEST_CASE("Round-trip conversion 2 (WGS84 -> Stereo70 -> WGS84)", "[epsg1838]")
    {
        const auto& params = conv::transformation::epsg1838;
        static constexpr T wgs84_lonlat_tolerance_deg = T(1e-10);
        static constexpr T altitude_tolerance_m = T(1e-5);

        round_trip_test(params, wgs84_lonlat_tolerance_deg, altitude_tolerance_m, "[epsg1838]");
    }

    TEST_CASE("Round-trip conversion 3 (WGS84 -> Stereo70 -> WGS84)", "[ancpi_stereo70_etrs89_approx]")
    {
        const auto& params = conv::transformation::ancpi_stereo70_etrs89_approx;
        static constexpr T wgs84_lonlat_tolerance_deg = T(1e-9);
        static constexpr T altitude_tolerance_m = T(1e-4);

        round_trip_test(params, wgs84_lonlat_tolerance_deg, altitude_tolerance_m, "[ancpi_stereo70_etrs89_approx]");
    }

    TEST_CASE("Round-trip conversion 4 (WGS84 -> Stereo70 -> WGS84)", "[epsg1188]")
    {
        const auto& params = conv::transformation::epsg1188;
        static constexpr T wgs84_lonlat_tolerance_deg = T(1e-9);
        static constexpr T altitude_tolerance_m = T(1e-4);

        round_trip_test(params, wgs84_lonlat_tolerance_deg, altitude_tolerance_m, "[epsg1188]");
    }

    TEST_CASE("Round-trip conversion 5 (WGS84 -> Stereo70 -> WGS84)", "[epsg15861]")
    {
        const auto& params = conv::transformation::epsg15861;
        static constexpr T wgs84_lonlat_tolerance_deg = T(1e-9);
        static constexpr T altitude_tolerance_m = T(1e-4);

        round_trip_test(params, wgs84_lonlat_tolerance_deg, altitude_tolerance_m, "[epsg15861]");
    }

    TEST_CASE("Round-trip conversion 6 (WGS84 -> Stereo70 -> WGS84)", "[pulkovo58_wgs84_ro_approx]")
    {
        const auto& params = conv::transformation::pulkovo58_wgs84_ro_approx;
        static constexpr T wgs84_lonlat_tolerance_deg = T(1e-9);
        static constexpr T altitude_tolerance_m = T(1e-3);

        round_trip_test(params, wgs84_lonlat_tolerance_deg, altitude_tolerance_m, "[pulkovo58_wgs84_ro_approx]");
    }

    TEST_CASE("Edge Cases", "[conversion][edge]")
    {
        SECTION("Invalid Input Handling")
        {
            // Test NaN inputs
            wgs84::coordinate<T> nan_coord {.latitude = std::numeric_limits<T>::quiet_NaN(), .longitude = 0, .altitude = 0};

            REQUIRE_THROWS_AS(conv::wgs84_to_stereo70(nan_coord), std::invalid_argument);

            // Test out of range latitude
            wgs84::coordinate<T> invalid_lat {.latitude = 100, .longitude = 0, .altitude = 0};

            REQUIRE_THROWS_AS(conv::wgs84_to_stereo70(invalid_lat), std::invalid_argument);
        }

        SECTION("Pole Handling")
        {
            // North Pole
            wgs84::coordinate<T> north_pole {.latitude = 90, .longitude = 0, .altitude = 0};

            REQUIRE_NOTHROW(conv::wgs84_to_stereo70(north_pole));

            // South Pole
            wgs84::coordinate<T> south_pole {.latitude = -90, .longitude = 0, .altitude = 0};

            REQUIRE_NOTHROW(conv::wgs84_to_stereo70(south_pole));
        }
    }
}
