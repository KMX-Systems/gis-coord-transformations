// Copyright (c) 2025 - present KMX Systems. All rights reserved.
#ifndef PCH
    #include <algorithm>
    #include <cmath>
    #include <limits>
    #include <numbers>
    #include <stdexcept>
    #include <type_traits>
#endif

namespace kmx::gis
{
    // Core Constants with better organization
    template <typename T = double>
    struct constants
    {
        static_assert(std::is_floating_point_v<T>, "constants requires a floating-point type.");
        constants() = delete; // Prevent instantiation

        // Common Values
        struct value
        {
            static constexpr T zero = T(0);
            static constexpr T one = T(1);
            static constexpr T two = T(2);
            static constexpr T four = T(4);
            static constexpr T six = T(6);
            static constexpr T seven = T(7);
            static constexpr T eight = T(8);
            static constexpr T ten = T(10);
            static constexpr T twelve = T(12.0);
            static constexpr T thirteen = T(13.0);
            static constexpr T fifteen = T(15.0);
            static constexpr T twenty_four = T(24.0);
            static constexpr T twenty_nine = T(29.0);
            static constexpr T forty_eight = T(48.0);
            static constexpr T eighty_one = T(81.0);
            static constexpr T ninety = T(90.0);
            static constexpr T one_twenty = T(120.0);
            static constexpr T one_eighty = T(180.0);
            static constexpr T two_forty = T(240.0);
            static constexpr T three_sixty = T(360.0);

            static constexpr T million = T(1e6);
        };

        // Mathematical Constants
        struct math
        {
            static constexpr T pi = std::numbers::pi_v<T>;
            static constexpr T two_pi = value::two * pi;
            static constexpr T half_pi = pi / value::two;
            static constexpr T quarter_pi = pi / value::four;
            static constexpr T epsilon = std::numeric_limits<T>::epsilon();
        };

        // Conversion Factors
        struct conversions
        {
            static constexpr T sec_per_degree = T(3600);
            static constexpr T deg_to_rad = math::pi / T(180);
            static constexpr T rad_to_deg = T(180) / math::pi;
            static constexpr T arcsec_to_rad = deg_to_rad / sec_per_degree;
            static constexpr T ppm_to_scale = value::one / value::million;
        };

        // Tolerances & Limits
        struct tolerance
        {
            static constexpr T default_geodetic = T(1e-11);
            static constexpr T near_zero = math::epsilon * T(100);
            static constexpr int max_iterations = 10;
        };

        // Series Coefficients
        struct series
        {
            static constexpr T c2_term1_num = value::seven;
            static constexpr T c2_term1_den = T(48);
            static constexpr T c2_term2_num = T(29);
            static constexpr T c2_term2_den = T(240);
            static constexpr T c2_term3_num = T(811);
            static constexpr T c2_term3_den = T(11520);
            static constexpr T c3_term1_num = value::seven;
            static constexpr T c3_term1_den = T(120);
            static constexpr T c3_term2_num = T(81);
            static constexpr T c3_term2_den = T(1120);
            static constexpr T c4_term1_num = T(4279);
            static constexpr T c4_term1_den = T(161280);
        };
    };

    // Core Structures with enhanced documentation
    template <typename T>
    struct ellipsoid
    {
        static_assert(std::is_floating_point_v<T>, "ellipsoid requires a floating-point type.");
        using K = constants<T>;

        T a;     // Semi-major axis
        T inv_f; // Inverse flattening
        T f;     // Flattening
        T e2;    // First eccentricity squared
        T e;     // First eccentricity
        T b;     // Semi-minor axis
        T ep2;   // Second eccentricity squared

        /// @brief Constructs an ellipsoid from semi-major axis and inverse flattening
        /// @throws std::invalid_argument for invalid parameters
        constexpr ellipsoid(const T semi_major, const T inverse_flattening):
            a(semi_major),
            inv_f(inverse_flattening),
            f(inv_f == K::value::zero ? K::value::zero : K::value::one / inv_f),
            e2(f == K::value::zero ? K::value::zero : K::value::two * f - f * f),
            e(std::sqrt(e2)),
            b(a * (K::value::one - f)),
            ep2(e2 == K::value::one ? std::numeric_limits<T>::infinity() : e2 / (K::value::one - e2))
        {
            if (a <= K::value::zero)
            {
                throw std::invalid_argument("Semi-major axis must be positive (got " + std::to_string(a) + ")");
            }

            if (inv_f < K::value::zero)
            {
                throw std::invalid_argument("Inverse flattening cannot be negative (got " + std::to_string(inv_f) + ")");
            }

            if ((K::value::one - e2) <= K::value::zero && e2 != K::value::one)
            {
                throw std::invalid_argument("Invalid ellipsoid parameters (1-e^2 = " + std::to_string(K::value::one - e2) + ")");
            }
        }

        /// @brief Returns whether this ellipsoid is a sphere
        constexpr bool is_sphere() const noexcept { return e2 == K::value::zero; }
    };

    // Enhanced helmert_params with better documentation
    template <typename T>
    struct helmert_params
    {
        static_assert(std::is_floating_point_v<T>, "helmert_params requires a floating-point type.");
        using K = constants<T>;

        T dx {};     // X translation (meters)
        T dy {};     // Y translation (meters)
        T dz {};     // Z translation (meters)
        T rx_sec {}; // X rotation (arcseconds)
        T ry_sec {}; // Y rotation (arcseconds)
        T rz_sec {}; // Z rotation (arcseconds)
        T ds_ppm {}; // Scale factor (parts per million)

        constexpr T rx_rad() const noexcept { return rx_sec * K::conversions::arcsec_to_rad; }
        constexpr T ry_rad() const noexcept { return ry_sec * K::conversions::arcsec_to_rad; }
        constexpr T rz_rad() const noexcept { return rz_sec * K::conversions::arcsec_to_rad; }
        constexpr T scale_factor() const noexcept { return K::value::one + ds_ppm * K::conversions::ppm_to_scale; }

        /// @brief Validates the parameters
        /// @throws std::invalid_argument if parameters are invalid
        constexpr void validate() const
        {
            if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz))
            {
                throw std::invalid_argument("Helmert translation parameters must be finite");
            }
            if (!std::isfinite(rx_sec) || !std::isfinite(ry_sec) || !std::isfinite(rz_sec))
            {
                throw std::invalid_argument("Helmert rotation parameters must be finite");
            }
            if (!std::isfinite(ds_ppm))
            {
                throw std::invalid_argument("Helmert scale parameter must be finite");
            }
        }
    };

    // Coordinate structures with enhanced documentation
    template <typename T>
    struct geodetic_coord
    {
        static_assert(std::is_floating_point_v<T>);
        using K = constants<T>;

        T latitude {};  // Degrees [-90, 90]
        T longitude {}; // Degrees [-180, 180]
        T height {};    // Meters

        /// @brief Normalizes the coordinates to standard ranges
        constexpr void normalize() noexcept
        {
            // Normalize latitude
            latitude = std::clamp(latitude, -K::value::ninety, K::value::ninety);

            // Normalize longitude
            while (longitude > K::value::one_eighty)
            {
                longitude -= K::value::three_sixty;
            }
            while (longitude <= -K::value::one_eighty)
            {
                longitude += K::value::three_sixty;
            }
        }

        /// @brief Validates the coordinates
        /// @throws std::invalid_argument if coordinates are invalid
        constexpr void validate() const
        {
            if (!std::isfinite(latitude) || !std::isfinite(longitude) || !std::isfinite(height))
            {
                throw std::invalid_argument("Geodetic coordinates must be finite");
            }
            if (latitude < -K::value::ninety || latitude > K::value::ninety)
            {
                throw std::invalid_argument("Latitude must be between -90 and 90 degrees (got " + std::to_string(latitude) + ")");
            }
        }
    };

    template <typename T>
    struct xyz_coord
    {
        static_assert(std::is_floating_point_v<T>);

        T x {};
        T y {};
        T z {};

        /// @brief Validates the coordinates
        /// @throws std::invalid_argument if coordinates are invalid
        constexpr void validate() const
        {
            using std::isfinite;
            if (!isfinite(x) || !isfinite(y) || !isfinite(z))
            {
                throw std::invalid_argument("XYZ coordinates must be finite");
            }
        }
    };

    // Specific Coordinate Type Aliases
    template <typename T>
    using geocentric_coord = xyz_coord<T>;

    template <typename T>
    using projected_coord = xyz_coord<T>;

    namespace wgs84
    {
        template <typename T = double>
        using coordinate = geodetic_coord<T>;
    }

    namespace stereo70
    {
        template <typename T = double>
        using coordinate = projected_coord<T>;

        // Projection parameters as a separate struct for better organization
        template <typename T>
        struct projection_params
        {
            static constexpr T lat0_deg = T(46.0); // Latitude of origin
            static constexpr T lon0_deg = T(25.0); // Central meridian
            static constexpr T k0 = T(0.99975);    // Scale factor
            static constexpr T fe = T(500000.0);   // False easting
            static constexpr T fn = T(500000.0);   // False northing

            // Tolerances
            static constexpr T origin_tol_deg = T(1e-9);
            static constexpr T antipodal_tol = T(1e-10);
            static constexpr T origin_proj_tol = T(1e-6);
        };
    }

    // Conversion Logic with improved organization
    template <typename T>
    struct conversion
    {
        using K = constants<T>;
        using stereo70_params = stereo70::projection_params<T>;

        // Ellipsoid definitions
        static constexpr ellipsoid<T> krasovsky1940 {
            T(6378245.0), // a
            T(298.3)      // 1/f
        };

        static constexpr ellipsoid<T> wgs84 {
            T(6378137.0),    // a
            T(298.257223563) // 1/f
        };

        // Default Helmert parameters
        static constexpr helmert_params<T> pulkovo58_to_wgs84 {.dx = T(33.4),
                                                               .dy = T(-146.6),
                                                               .dz = T(-76.3),
                                                               .rx_sec = T(-0.359),
                                                               .ry_sec = T(-0.053),
                                                               .rz_sec = T(0.844),
                                                               .ds_ppm = T(-0.84)};

        // Derived projection parameters with better initialization
        struct derived_params
        {
            T lat0_rad;   // Latitude of origin in radians
            T lon0_rad;   // Central meridian in radians
            T e;          // Eccentricity
            T e2;         // Eccentricity squared
            T e4, e6, e8; // Higher powers of e
            T sin_lat0, cos_lat0;
            T chi0; // Conformal latitude of origin
            T sin_chi0, cos_chi0;
            T n0;             // Intermediate calculation
            T inv_two_n0;     // 1/(2*n0)
            T c1, c2, c3, c4; // Series coefficients

            constexpr derived_params():
                lat0_rad(stereo70_params::lat0_deg * K::conversions::deg_to_rad),
                lon0_rad(stereo70_params::lon0_deg * K::conversions::deg_to_rad),
                e(krasovsky1940.e),
                e2(krasovsky1940.e2),
                e4(e2 * e2),
                e6(e4 * e2),
                e8(e6 * e2),
                sin_lat0(std::sin(lat0_rad)),
                cos_lat0(std::cos(lat0_rad)),
                chi0(calculate_conformal_latitude(lat0_rad, e)),
                sin_chi0(std::sin(chi0)),
                cos_chi0(std::cos(chi0)),
                n0(calculate_n0()),
                inv_two_n0(T(0.5) / n0),
                c1(calculate_c1()),
                c2(calculate_c2()),
                c3(calculate_c3()),
                c4(calculate_c4())
            {
                if (!is_valid())
                {
                    throw std::runtime_error("Invalid derived projection parameters");
                }
            }

        private:
            constexpr T calculate_n0() const
            {
                const T w1_sq = K::value::one - e2 * sin_lat0 * sin_lat0;
                if (w1_sq <= K::value::zero)
                {
                    return std::numeric_limits<T>::quiet_NaN();
                }
                const T denominator = std::sqrt(w1_sq) * cos_chi0;
                if (std::abs(denominator) < K::tolerance::near_zero)
                {
                    return std::numeric_limits<T>::quiet_NaN();
                }
                return stereo70_params::k0 * krasovsky1940.a * cos_lat0 / denominator;
            }

            constexpr T calculate_c1() const { return e2 / T(2) + T(5) * e4 / T(24) + e6 / T(12) + T(13) * e8 / T(360); }

            constexpr T calculate_c2() const
            {
                return (K::series::c2_term1_num / K::series::c2_term1_den) * e4 + (K::series::c2_term2_num / K::series::c2_term2_den) * e6 +
                       (K::series::c2_term3_num / K::series::c2_term3_den) * e8;
            }

            constexpr T calculate_c3() const
            {
                return (K::series::c3_term1_num / K::series::c3_term1_den) * e6 + (K::series::c3_term2_num / K::series::c3_term2_den) * e8;
            }

            constexpr T calculate_c4() const { return (K::series::c4_term1_num / K::series::c4_term1_den) * e8; }

            constexpr bool is_valid() const
            {
                using std::isfinite;
                return isfinite(lat0_rad) && isfinite(lon0_rad) && isfinite(e) && isfinite(e2) && isfinite(sin_lat0) &&
                       isfinite(cos_lat0) && isfinite(chi0) && isfinite(sin_chi0) && isfinite(cos_chi0) && isfinite(n0) &&
                       std::abs(n0) > K::tolerance::near_zero && isfinite(inv_two_n0) && isfinite(c1) && isfinite(c2) && isfinite(c3) &&
                       isfinite(c4);
            }
        };

        static constexpr derived_params proj_params {};

        // Helper functions with improved error handling
        static constexpr T calculate_conformal_latitude(T lat_rad, T e)
        {
            const T sin_lat = std::sin(lat_rad);
            const T esin_lat = e * sin_lat;

            // Handle poles and extreme cases
            if (std::abs(esin_lat) >= K::value::one)
            {
                return std::copysign(K::math::half_pi, lat_rad);
            }

            const T term1 = std::tan(K::math::quarter_pi + lat_rad * T(0.5));
            const T term2_num = K::value::one - esin_lat;
            const T term2_den = K::value::one + esin_lat;

            if (std::abs(term2_den) < K::tolerance::near_zero || term2_num <= K::value::zero)
            {
                return std::numeric_limits<T>::quiet_NaN();
            }

            const T term2 = std::pow(term2_num / term2_den, e * T(0.5));
            const T atan_arg = term1 * term2;

            if (!std::isfinite(atan_arg))
            {
                return std::copysign(K::math::half_pi, atan_arg);
            }

            return T(2) * std::atan(atan_arg) - K::math::half_pi;
        }

        // Core transformation functions with improved robustness
        static geocentric_coord<T> geodetic_to_geocentric(const geodetic_coord<T>& geo, const ellipsoid<T>& ellip)
        {
            geo.validate();

            const T lat_rad = geo.latitude * K::conversions::deg_to_rad;
            const T lon_rad = geo.longitude * K::conversions::deg_to_rad;
            const T sin_lat = std::sin(lat_rad);
            const T cos_lat = std::cos(lat_rad);

            const T n_denominator = K::value::one - ellip.e2 * sin_lat * sin_lat;
            if (n_denominator <= K::value::zero)
            {
                throw std::runtime_error("Invalid geodetic coordinates (N denominator = " + std::to_string(n_denominator) + ")");
            }

            const T N = ellip.a / std::sqrt(n_denominator);
            const T h = geo.height;

            return {.x = (N + h) * cos_lat * std::cos(lon_rad),
                    .y = (N + h) * cos_lat * std::sin(lon_rad),
                    .z = (N * (K::value::one - ellip.e2) + h) * sin_lat};
        }

        static geodetic_coord<T> geocentric_to_geodetic(const geocentric_coord<T>& ecef, const ellipsoid<T>& ellip,
                                                        T tolerance = K::tolerance::default_geodetic,
                                                        int max_iterations = K::tolerance::max_iterations)
        {
            ecef.validate();

            const T p = std::hypot(ecef.x, ecef.y);
            geodetic_coord<T> result;

            // Handle polar cases
            const T pole_threshold = K::tolerance::near_zero * ellip.a;
            if (p < pole_threshold)
            {
                result.longitude = K::value::zero;
                if (std::abs(ecef.z) < pole_threshold)
                {
                    result.latitude = K::value::zero;
                    result.height = -ellip.a;
                }
                else
                {
                    result.latitude = std::copysign(K::value::ninety, ecef.z);
                    result.height = std::abs(ecef.z) - ellip.b;
                }
                return result;
            }

            result.longitude = std::atan2(ecef.y, ecef.x) * K::conversions::rad_to_deg;

            // Initial approximation
            T lat_rad = std::atan2(ecef.z, p * (K::value::one - ellip.e2));
            T N = ellip.a;
            T h = K::value::zero;

            for (int i = 0; i < max_iterations; ++i)
            {
                const T sin_lat = std::sin(lat_rad);
                const T cos_lat = std::cos(lat_rad);

                const T N_denominator = K::value::one - ellip.e2 * sin_lat * sin_lat;
                if (N_denominator <= K::value::zero)
                {
                    throw std::runtime_error("Invalid geocentric coordinates during iteration");
                }

                N = ellip.a / std::sqrt(N_denominator);

                if (std::abs(cos_lat) < K::tolerance::near_zero)
                {
                    // Near pole
                    h = std::abs(ecef.z) - ellip.b;
                    lat_rad = std::copysign(K::math::half_pi, ecef.z);
                    break;
                }
                else
                {
                    h = p / cos_lat - N;
                }

                const T h_plus_N = N + h;
                if (std::abs(h_plus_N) < K::tolerance::near_zero)
                {
                    lat_rad = (std::abs(ecef.z) < K::tolerance::near_zero) ?
                                  K::value::zero :
                                  std::atan2(ecef.z, p * (K::value::one - ellip.e2 * N / h_plus_N));
                }
                else
                {
                    const T lat_denom = p * (K::value::one - ellip.e2 * N / h_plus_N);
                    lat_rad = (std::abs(lat_denom) < K::tolerance::near_zero) ? std::copysign(K::math::half_pi, ecef.z) :
                                                                                std::atan2(ecef.z, lat_denom);
                }

                if (i > 0 && std::abs(std::sin(lat_rad - std::asin(sin_lat))) < tolerance)
                {
                    break;
                }
            }

            result.latitude = lat_rad * K::conversions::rad_to_deg;
            result.height = h;
            result.normalize();
            return result;
        }

        // Projection functions with improved numerical stability
        static projected_coord<T> krasovsky_geodetic_to_stereo70(const geodetic_coord<T>& geo)
        {
            geo.validate();

            // Check for origin point
            if (std::abs(geo.latitude - stereo70_params::lat0_deg) < stereo70_params::origin_tol_deg &&
                std::abs(geo.longitude - stereo70_params::lon0_deg) < stereo70_params::origin_tol_deg)
            {
                return {.x = stereo70_params::fe, .y = stereo70_params::fn, .z = geo.height};
            }

            const T lat_rad = geo.latitude * K::conversions::deg_to_rad;
            const T lon_rad = geo.longitude * K::conversions::deg_to_rad;

            const T chi = calculate_conformal_latitude(lat_rad, proj_params.e);
            if (!std::isfinite(chi))
            {
                throw std::runtime_error("Failed to calculate conformal latitude");
            }

            const T sin_chi = std::sin(chi);
            const T cos_chi = std::cos(chi);

            // Normalize longitude difference
            T lambda_diff = lon_rad - proj_params.lon0_rad;
            while (lambda_diff > K::math::pi)
                lambda_diff -= K::math::two_pi;
            while (lambda_diff <= -K::math::pi)
                lambda_diff += K::math::two_pi;

            const T sin_lambda = std::sin(lambda_diff);
            const T cos_lambda = std::cos(lambda_diff);

            const T denominator = K::value::one + proj_params.sin_chi0 * sin_chi + proj_params.cos_chi0 * cos_chi * cos_lambda;

            if (std::abs(denominator) < stereo70_params::antipodal_tol)
            {
                throw std::runtime_error("Coordinate too close to antipodal point");
            }

            const T k_prime = K::value::two * proj_params.n0 / denominator;

            projected_coord<T> result {.x = stereo70_params::fe + k_prime * cos_chi * sin_lambda,
                                       .y = stereo70_params::fn +
                                            k_prime * (proj_params.cos_chi0 * sin_chi - proj_params.sin_chi0 * cos_chi * cos_lambda),
                                       .z = geo.height};

            if (!std::isfinite(result.x) || !std::isfinite(result.y))
            {
                throw std::runtime_error("Non-finite projection result");
            }

            return result;
        }

        static stereo70::coordinate<T> wgs84_to_stereo70(const wgs84::coordinate<T>& wgs84_coord,
                                                         const helmert_params<T>& params = conversion<T>::pulkovo58_to_wgs84)
        {
            params.validate();

            // 1. WGS84 geodetic -> WGS84 geocentric
            const auto wgs84_ecef = geodetic_to_geocentric(wgs84_coord, wgs84);

            // 2. WGS84 geocentric -> Pulkovo geocentric
            const auto pulkovo_ecef = helmert_transformation_reverse(wgs84_ecef, params);

            // 3. Pulkovo geocentric -> Krasovsky geodetic
            const auto krasovsky_geo = geocentric_to_geodetic(pulkovo_ecef, krasovsky1940);

            // 4. Krasovsky geodetic -> Stereo70 projected
            return krasovsky_geodetic_to_stereo70(krasovsky_geo);
        }

        static wgs84::coordinate<T> stereo70_to_wgs84(const stereo70::coordinate<T>& stereo_coord,
                                                      const helmert_params<T>& params = conversion<T>::pulkovo58_to_wgs84)
        {
            params.validate();
            stereo_coord.validate();

            // 1. Stereo70 projected -> Krasovsky geodetic
            const auto krasovsky_geo = stereo70_to_krasovsky_geodetic(stereo_coord);

            // 2. Krasovsky geodetic -> Pulkovo geocentric
            const auto pulkovo_ecef = geodetic_to_geocentric(krasovsky_geo, krasovsky1940);

            // 3. Pulkovo geocentric -> WGS84 geocentric
            const auto wgs84_ecef = helmert_transformation_forward(pulkovo_ecef, params);

            // 4. WGS84 geocentric -> WGS84 geodetic
            return geocentric_to_geodetic(wgs84_ecef, wgs84);
        }

        static geodetic_coord<T> stereo70_to_krasovsky_geodetic(const projected_coord<T>& proj)
        {
            const T x_rel = proj.x - stereo70_params::fe;
            const T y_rel = proj.y - stereo70_params::fn;
            const T rho = std::hypot(x_rel, y_rel);

            // Handle origin point
            if (rho < stereo70_params::origin_proj_tol)
            {
                return {.latitude = stereo70_params::lat0_deg, .longitude = stereo70_params::lon0_deg, .height = proj.z};
            }

            const T tan_c_half = rho * proj_params.inv_two_n0;
            if (!std::isfinite(tan_c_half))
            {
                throw std::runtime_error("Invalid projection coordinates");
            }

            const T c = T(2) * std::atan(tan_c_half);
            const T sin_c = std::sin(c);
            const T cos_c = std::cos(c);

            const T chi = std::asin(
                std::clamp(cos_c * proj_params.sin_chi0 + (y_rel * sin_c * proj_params.cos_chi0 / rho), -K::value::one, K::value::one));

            const T lambda_diff = std::atan2(x_rel * sin_c, rho * proj_params.cos_chi0 * cos_c - y_rel * proj_params.sin_chi0 * sin_c);

            const T lat_rad = calculate_geodetic_latitude_series(chi, proj_params.c1, proj_params.c2, proj_params.c3, proj_params.c4);

            geodetic_coord<T> result {.latitude = lat_rad * K::conversions::rad_to_deg,
                                      .longitude = (proj_params.lon0_rad + lambda_diff) * K::conversions::rad_to_deg,
                                      .height = proj.z};

            result.normalize();
            return result;
        }

    private:
        // Private helper functions
        static constexpr T calculate_geodetic_latitude_series(T chi, T c1, T c2, T c3, T c4) noexcept
        {
            const T chi2 = T(2) * chi;
            const T chi4 = T(4) * chi;
            const T chi6 = T(6) * chi;
            const T chi8 = T(8) * chi;

            return chi + c1 * std::sin(chi2) + c2 * std::sin(chi4) + c3 * std::sin(chi6) + c4 * std::sin(chi8);
        }

        static geocentric_coord<T> helmert_transformation_forward(const geocentric_coord<T>& source, const helmert_params<T>& params)
        {
            const T rx = params.rx_rad();
            const T ry = params.ry_rad();
            const T rz = params.rz_rad();
            const T scale = params.scale_factor();

            // Apply rotation (small angle approximation)
            const T x_rot = source.x - rz * source.y + ry * source.z;
            const T y_rot = rz * source.x + source.y - rx * source.z;
            const T z_rot = -ry * source.x + rx * source.y + source.z;

            // Apply scale and translation
            return {.x = params.dx + scale * x_rot, .y = params.dy + scale * y_rot, .z = params.dz + scale * z_rot};
        }

        static geocentric_coord<T> helmert_transformation_reverse(const geocentric_coord<T>& target, const helmert_params<T>& params)
        {
            const T rx = params.rx_rad();
            const T ry = params.ry_rad();
            const T rz = params.rz_rad();
            const T scale = params.scale_factor();

            if (std::abs(scale) < K::tolerance::near_zero)
            {
                throw std::runtime_error("Invalid scale factor in reverse Helmert transformation");
            }

            // Remove translation
            const T x_temp = target.x - params.dx;
            const T y_temp = target.y - params.dy;
            const T z_temp = target.z - params.dz;

            // Reverse scale
            const T x_scaled = x_temp / scale;
            const T y_scaled = y_temp / scale;
            const T z_scaled = z_temp / scale;

            // Reverse rotation (small angle approximation)
            return {.x = x_scaled + rz * y_scaled - ry * z_scaled,
                    .y = -rz * x_scaled + y_scaled + rx * z_scaled,
                    .z = ry * x_scaled - rx * y_scaled + z_scaled};
        }
    };
} // namespace kmx::gis
